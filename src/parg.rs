//! The PARG scanner layer — the grammar-artifact loader and the
//! capture-tree utilities every PARG-driven reader shares. The
//! scanners (molfile, the InChI layers, RInChI) run compiled PARG
//! artifacts (grammar/*.json, checksummed in grammar/CHECKSUMS.md)
//! on parsanol's portable engine; the artifacts carry their own
//! sha256 envelope checksum, verified at every load, so a drifted
//! artifact refuses to parse rather than parsing differently.
//!
//! Layout mirrors the asciichem-grammar contract repo: the PARG text
//! contracts (grammar/*.parg) are the source of truth, the compiled
//! artifacts are committed beside them, and re-compiling +
//! re-committing both is part of any grammar change.
//!
//! The walk is zero-copy: the parse result's AST stores the matched
//! text as (offset, length) references into the input, and the
//! segment streams keep those spans — no capture text is
//! materialized. The packrat memo cache is thread-local and reused
//! across parses (a fresh cache per parse dominated the identity
//! pipeline's parse cost by an order of magnitude).

use parsanol::portable::ast::AstNode;
use parsanol::portable::cache::DenseCache;
use parsanol::portable::grammar::Grammar;
use parsanol::{AstArena, PargArtifact, PortableParser};
use std::sync::OnceLock;

/// One positional segment of a parsed line: the input span it covers
/// and the column it starts at. Every PARG line shape tiles its whole
/// line with ordered segments, so the scanners' column reads
/// (`get(a..b)` over the line) slice this stream byte-for-byte.
#[derive(Debug, Clone)]
pub(crate) struct Segments<'input> {
    input: &'input str,
    /// (start column, input offset, byte length), ordered by column.
    items: Vec<(usize, usize, usize)>,
}

impl<'input> Segments<'input> {
    /// The line's verbatim text. Literal granules (the `M  ` / `A  `
    /// prefixes) intern as arena strings and carry no input span, so
    /// the line's start derives from the first spanned segment: its
    /// input offset minus its column. A single trailing carriage
    /// return is dropped, mirroring `str::lines`.
    pub(crate) fn raw(&self) -> &'input str {
        let Some(first) = self.items.iter().find(|(_, _, length)| *length > 0) else {
            return "";
        };
        let start = first.1 - first.0;
        let last = self
            .items
            .iter()
            .map(|(_, offset, length)| offset + length)
            .max()
            .unwrap_or(start);
        if last <= start {
            return "";
        }
        let raw = &self.input[start..last];
        raw.strip_suffix('\r').unwrap_or(raw)
    }

    /// The engine's column read: `line.get(a..b)` over the parsed
    /// line. `None` when the line ends before the column — the
    /// readers' tolerant field defaults apply exactly where they did
    /// when the line was sliced directly.
    pub(crate) fn column(&self, start: usize, end: usize) -> Option<&'input str> {
        for (column, offset, length) in &self.items {
            let finish = column + length;
            if *column <= start && finish >= end {
                return self
                    .input
                    .get(offset + start - column..offset + end - column);
            }
        }
        None
    }

    /// Builds the segment stream from a line-shaped capture node. The
    /// use-site captures wrap the shape rule one or more levels deep
    /// (`line: {body_line: {atom_line: {k0, ...}}}`); single-key
    /// wrappers peel off until the positional granules appear. Fixed
    /// granules carry their start column in the capture name (`k0`,
    /// `k3`, ...); the trailing token run (`units` then `kt`) tiles
    /// the rest of the line contiguously and rides as one segment.
    fn from_node(node: &AstNode, capture: &Capture<'input>) -> Segments<'input> {
        // Collect the shape's named captures: rule wrappers nest the
        // granule hashes one or more levels deep (single-key hashes,
        // and sequences materialize as arrays of those), so the walk
        // flattens named items until the column granules appear.
        fn flatten(node: &AstNode, capture: &Capture<'_>, into: &mut Vec<(String, AstNode)>) {
            match *node {
                AstNode::Hash { pool_index, length } => {
                    let items = capture
                        .arena
                        .get_hash_items(pool_index as usize, length as usize);
                    if items
                        .iter()
                        .any(|(key, _)| is_granule(key) || key == "units" || key == "kt")
                    {
                        into.extend(items);
                    } else {
                        for (_, inner) in items {
                            flatten(&inner, capture, into);
                        }
                    }
                }
                AstNode::Array { pool_index, length } => {
                    for inner in capture
                        .arena
                        .get_array(pool_index as usize, length as usize)
                    {
                        flatten(&inner, capture, into);
                    }
                }
                _ => {}
            }
        }
        let mut named: Vec<(String, AstNode)> = Vec::new();
        flatten(node, capture, &mut named);
        let mut segments: Vec<(usize, usize, usize)> = Vec::new();
        let mut trail: Option<(usize, usize)> = None;
        for (key, inner) in &named {
            if key == "units" || key == "kt" {
                let span = capture.span(inner);
                trail = Some(match trail {
                    Some((start, end)) => (start.min(span.0), end.max(span.1)),
                    None => span,
                });
            } else if let Some(column) = key.strip_prefix('k').and_then(|c| c.parse().ok()) {
                let (offset, end) = capture.span(inner);
                segments.push((column, offset, end - offset));
            }
        }
        segments.sort_by_key(|(column, _, _)| *column);
        if let Some((start, end)) = trail {
            let cursor = segments
                .last()
                .map_or(0, |(column, _, length)| column + length);
            if end > start {
                segments.push((cursor, start, end - start));
            }
        }
        Segments {
            input: capture.input,
            items: segments,
        }
    }
}

/// Whether a capture name is a column granule (`k0`, `k15`, ...).
fn is_granule(key: &str) -> bool {
    key.starts_with('k') && key[1..].parse::<usize>().is_ok()
}

/// Whether a capture name names a layer shape or the raw fallback.
pub(crate) fn is_shape_name(name: &str) -> bool {
    name.ends_with("_seg") || name == "raw"
}

/// A parsed capture tree over one input: the arena-backed AST plus
/// the input itself, handed to a scanner's walk closure.
pub(crate) struct Capture<'input> {
    root: AstNode,
    arena: AstArena,
    input: &'input str,
}

impl<'input> Capture<'input> {
    /// The named children of a capture hash.
    pub(crate) fn hash_items(&self, node: &AstNode) -> Vec<(String, AstNode)> {
        match *node {
            AstNode::Hash { pool_index, length } => self
                .arena
                .get_hash_items(pool_index as usize, length as usize),
            _ => Vec::new(),
        }
    }

    /// A capture's input span: the byte range its matched text covers.
    /// Empty children (lookahead guards and zero-width matches match
    /// nothing) are skipped — a repetition's iteration arrays carry
    /// them as leading elements, and their (0, 0) would poison the
    /// minimum otherwise.
    fn span(&self, node: &AstNode) -> (usize, usize) {
        match *node {
            AstNode::InputRef { offset, length } if length > 0 => {
                (offset as usize, offset as usize + length as usize)
            }
            AstNode::Array { pool_index, length } => {
                let items = self.arena.get_array(pool_index as usize, length as usize);
                self.span_of(items.iter())
            }
            AstNode::Hash { pool_index, length } => {
                let items = self
                    .arena
                    .get_hash_items(pool_index as usize, length as usize);
                self.span_of(items.iter().map(|(_, node)| node))
            }
            _ => (0, 0),
        }
    }

    /// The covering span of the non-empty children.
    fn span_of<'n>(&self, items: impl Iterator<Item = &'n AstNode>) -> (usize, usize) {
        let mut span: Option<(usize, usize)> = None;
        for item in items {
            let (start, end) = self.span(item);
            if end > start {
                span = Some(match span {
                    Some((s, e)) => (s.min(start), e.max(end)),
                    None => (start, end),
                });
            }
        }
        span.unwrap_or((0, 0))
    }

    /// A named child capture.
    pub(crate) fn child(&self, node: &AstNode, key: &str) -> Option<AstNode> {
        self.hash_items(node)
            .into_iter()
            .find(|(name, _)| name == key)
            .map(|(_, inner)| inner)
    }

    /// A named child capture anywhere in a shape node: the AST stores
    /// sequences as arrays whose captured children are single-key
    /// hashes, so the lookup walks both forms.
    pub(crate) fn deep_child(&self, node: &AstNode, key: &str) -> Option<AstNode> {
        for (name, value) in self.hash_items(node) {
            if name == key {
                return Some(value);
            }
        }
        for (_, value) in self.hash_items(node) {
            if let Some(found) = self.deep_child(&value, key) {
                return Some(found);
            }
        }
        for item in self.array_items(node) {
            if self.is_empty_marker(&item) {
                continue;
            }
            if let Some(found) = self.deep_child(&item, key) {
                return Some(found);
            }
        }
        None
    }

    /// A named child capture inside a sequence's array: the AST stores
    /// a rule's captured children as single-key hashes in element
    /// order, nested through optionals' sub-arrays.
    pub(crate) fn array_child(&self, node: &AstNode, key: &str) -> Option<AstNode> {
        for item in self.array_items(node) {
            if let Some(found) = self.descend_child(&item, key, 0) {
                return Some(found);
            }
        }
        None
    }

    fn descend_child(&self, node: &AstNode, key: &str, depth: usize) -> Option<AstNode> {
        if depth > 8 {
            return None;
        }
        let mut items = self.hash_items(node);
        if items.len() == 1 {
            if items[0].0 == key {
                return Some(items.remove(0).1);
            }
            return None;
        }
        for item in self.array_items(node) {
            if let Some(found) = self.descend_child(&item, key, depth + 1) {
                return Some(found);
            }
        }
        None
    }

    /// The ordered elements of an array node (a sequence or
    /// repetition's raw AST form).
    pub(crate) fn array_items(&self, node: &AstNode) -> Vec<AstNode> {
        match *node {
            AstNode::Array { pool_index, length } => {
                self.arena.get_array(pool_index as usize, length as usize)
            }
            _ => Vec::new(),
        }
    }

    /// Whether a repetition-array element is the match's empty seed
    /// marker (not an iteration).
    pub(crate) fn is_empty_marker(&self, node: &AstNode) -> bool {
        !matches!(*node, AstNode::Hash { .. } | AstNode::Array { .. })
    }

    /// A capture's matched text: its input span slice.
    pub(crate) fn text(&self, node: &AstNode) -> &'input str {
        let (start, end) = self.span(node);
        if end > start {
            self.input.get(start..end).unwrap_or("")
        } else {
            ""
        }
    }

    /// The single-key wrap a named rule capture produces
    /// (`{element: {...}}`): the inner value under the sole key.
    pub(crate) fn unwrap_single(&self, node: &AstNode) -> Option<(String, AstNode)> {
        let mut items = self.hash_items(node);
        if items.len() == 1 {
            Some(items.remove(0))
        } else {
            None
        }
    }

    /// The first shape-named single-key hash in a node's subtree
    /// (a choice branch nested through the sequence arrays).
    pub(crate) fn find_shape(&self, node: &AstNode, depth: usize) -> Option<(String, AstNode)> {
        if depth > 8 {
            return None;
        }
        let mut items = self.hash_items(node);
        if items.len() == 1 {
            let (name, value) = items.remove(0);
            if is_shape_name(&name) {
                return Some((name, value));
            }
            return self.find_shape(&value, depth + 1);
        }
        for item in self.array_items(node) {
            if self.is_empty_marker(&item) {
                continue;
            }
            if let Some(found) = self.find_shape(&item, depth + 1) {
                return Some(found);
            }
        }
        None
    }

    /// A repetition's iteration values as captured: multi-key
    /// iteration hashes ride as-is, empty seed markers skipped.
    pub(crate) fn iterations(&self, node: &AstNode) -> Vec<AstNode> {
        self.repetition_items(node)
            .into_iter()
            .filter(|item| !self.is_empty_marker(item))
            .collect()
    }

    /// A repetition capture's iteration values, each resolved to the
    /// named capture inside it (the iteration may wrap it in a
    /// single-key hash or carry it inside a sequence array).
    pub(crate) fn reps(&self, node: &AstNode, key: &str) -> Vec<AstNode> {
        let mut out = Vec::new();
        for item in self.repetition_items(node) {
            if self.is_empty_marker(&item) {
                continue;
            }
            self.collect_key(&item, key, 0, &mut out);
        }
        out
    }

    /// Collects every named-capture descendant in iteration order —
    /// a repetition's iterations nest through sequence arrays and
    /// single-key wraps, and all of them carry the capture.
    #[allow(clippy::only_used_in_recursion)]
    fn collect_key(&self, node: &AstNode, key: &str, depth: usize, out: &mut Vec<AstNode>) {
        if depth > 12 {
            return;
        }
        for (name, value) in self.hash_items(node) {
            if name == key {
                out.push(value);
            } else {
                self.collect_key(&value, key, depth + 1, out);
            }
        }
        for item in self.array_items(node) {
            if self.is_empty_marker(&item) {
                continue;
            }
            self.collect_key(&item, key, depth + 1, out);
        }
    }

    /// A repetition capture's array: single-key capture wrappers peel
    /// off until the iteration array appears.
    fn repetition_items(&self, node: &AstNode) -> Vec<AstNode> {
        let mut current = node.clone();
        loop {
            if matches!(current, AstNode::Array { .. }) {
                return self.array_items(&current);
            }
            let mut items = self.hash_items(&current);
            if items.len() == 1 {
                current = items.remove(0).1;
            } else {
                return Vec::new();
            }
        }
    }

    /// A repetition capture's iteration values, unwrapped from their
    /// single-key hashes and skipping the empty seed marker.
    pub(crate) fn iter_items(&self, node: &AstNode, key: &str) -> Vec<AstNode> {
        let mut out = Vec::new();
        for item in self.array_items(node) {
            if self.is_empty_marker(&item) {
                continue;
            }
            if let Some((name, inner)) = self.unwrap_single(&item) {
                if name == key {
                    out.push(inner);
                }
            }
        }
        out
    }

    /// The line segment stream of a body-unit or line capture.
    pub(crate) fn segments(&self, node: &AstNode) -> Segments<'input> {
        Segments::from_node(node, self)
    }
}

/// Parses `input` against a compiled artifact entry and hands the
/// capture tree to `walk`. The whole input must be consumed — a short
/// capture is a scanner error, never a silently partial parse.
fn parse_entry<'input, T>(
    grammar: &Grammar,
    entry: &'static str,
    input: &'input str,
    walk: impl FnOnce(Capture<'input>) -> crate::Result<T>,
) -> crate::Result<T> {
    // A fresh packrat memo cache per parse: the cache's memo entries
    // reference arena pool state, so reuse across parses (fresh
    // arenas) hands the next parse stale nodes — probed in parsanol
    // 0.18.9, where a hit walks a dead arena. The estimate keeps the
    // table one allocation for the scanners' line-shaped inputs.
    let cache = DenseCache::new(1 << 12);
    let mut arena = AstArena::for_input(input.len());
    let mut parser = PortableParser::new_with_cache(grammar, input, &mut arena, cache);
    let result = parser
        .parse_with_end_pos()
        .map_err(|error| crate::error::Error::new(format!("parg {entry}: {error}")))?;
    if result.end_pos < input.len() {
        return Err(crate::error::Error::new(format!(
            "parg {entry}: input not fully captured: consumed {} of {} bytes",
            result.end_pos,
            input.len()
        )));
    }
    walk(Capture {
        root: result.value,
        arena,
        input,
    })
}

/// Loads one compiled artifact entry, once.
macro_rules! parg_grammar {
    ($function:ident, $bytes:expr, $entry:literal) => {
        pub(crate) fn $function() -> crate::Result<&'static Grammar> {
            static CELL: OnceLock<Result<Grammar, String>> = OnceLock::new();
            static ENTRY: &str = $entry;
            let loaded = CELL.get_or_init(|| {
                PargArtifact::from_json($bytes)
                    .map_err(|error| format!("PARG artifact load failed: {error:?}"))?
                    .grammar(ENTRY)
                    .map_err(|error| format!("PARG entry extraction failed: {error:?}"))
            });
            match loaded {
                Ok(grammar) => Ok(grammar),
                Err(message) => Err(crate::error::Error::new(message.clone())),
            }
        }
    };
}

parg_grammar!(
    molfile_grammar,
    include_str!("../grammar/molfile.json"),
    "mol"
);
parg_grammar!(
    inchi_grammar,
    include_str!("../grammar/inchi.json"),
    "inchi"
);
parg_grammar!(
    rinchi_grammar,
    include_str!("../grammar/rinchi.json"),
    "rinchi"
);

/// Parses the molfile contract's entry and walks the capture tree.
pub(crate) fn with_molfile_capture<T>(
    input: &str,
    walk: impl FnOnce(Capture<'_>) -> crate::Result<T>,
) -> crate::Result<T> {
    parse_entry(molfile_grammar()?, "molfile", input, walk)
}

/// Parses the InChI contract's entry and walks the capture tree.
pub(crate) fn with_inchi_capture<T>(
    input: &str,
    walk: impl FnOnce(Capture<'_>) -> crate::Result<T>,
) -> crate::Result<T> {
    parse_entry(inchi_grammar()?, "inchi", input, walk)
}

/// Parses the RInChI contract's entry and walks the capture tree.
pub(crate) fn with_rinchi_capture<T>(
    input: &str,
    walk: impl FnOnce(Capture<'_>) -> crate::Result<T>,
) -> crate::Result<T> {
    parse_entry(rinchi_grammar()?, "rinchi", input, walk)
}

/// The parse root: a scanner's entry capture.
impl Capture<'_> {
    pub(crate) fn root(&self) -> &AstNode {
        &self.root
    }
}
