//! The day-1 SMILES reader (chemomer#121): the `smiles.parg` capture
//! tree → [`crate::graph::Structure`]. The grammar (day-1 scope:
//! chains, branches, single-digit ring closures, organic + bracket
//! atoms, charges, isotopes, explicit H counts, aromatic markers)
//! captures a flat token stream in source order — atoms (`a`), bond
//! symbols (`b`), ring-closure digits (`r`) — with branch chains
//! nested under `group`; the walk is a state machine over that
//! stream: the branch point is saved around each `group`, ring
//! openings pair by digit value, and the implicit bond reads aromatic
//! (the molfile's bond code 4) between two aromatic atoms, single
//! otherwise.
//!
//! Phase 2: tetrahedral `@`/`@@` markers and the directional `/` `\`
//! bonds ride through — `synthesize_stereo` turns the markers into
//! drawn geometry (a trigonal ligand layout plus a wedge realizing
//! the @/@@ hand per the SMILES neighbor order; the substituents
//! lifted to the side their spelling draws) so the wedge perception
//! and the cis/trans side product read a SMILES structure exactly
//! like a drawn molfile. The alanine enantiomers derive the oracle
//! bytes; the /b signs were verified against the oracle over
//! coordinate drawings. The component dot landed (chemomer#131):
//! dot-disconnected components read as one structure the way the
//! molfile's multi-fragment path does; the identity-multicomponent
//! tranche rides the SMILES gate.

use crate::error::{err, Error, Result};
use crate::graph::{AtomRecord, Structure};

use crate::parg::{with_smiles_capture, Capture};

/// The bond kinds the SMILES reader emits (the V2000 codes).
const BOND_SINGLE: u32 = 1;
const BOND_AROMATIC: u32 = 4;

/// Reads a day-1 SMILES string into a structure.
pub fn read(input: &str) -> Result<Structure> {
    with_smiles_capture(input, |capture| {
        let mut walk = Walk {
            capture: &capture,
            structure: Structure::default(),
            current: None,
            pending_bond: None,
            pending_direction: None,
            stack: Vec::new(),
            rings: std::collections::HashMap::new(),
            aromatic_markers: Vec::new(),
        };
        walk.formula(capture.root())?;
        if let Some(digit) = walk.rings.keys().next() {
            return err(format!("smiles: ring bond {digit} never closed"));
        }
        if walk.structure.atoms.is_empty() {
            return err("smiles: no atoms");
        }
        fold_isotopic_hydrogens(&mut walk.structure);
        synthesize_stereo(&mut walk.structure)?;
        Ok(walk.structure)
    })
}

/// Folds isotopic hydrogens onto their bonded carriers — the molfile
/// reader's parse-time semantics ([`Structure::isotopic_h`], the /h
/// layer is isotope-blind and the /i layer reads this map). An
/// explicit `[2H]`/`[3H]` atom in SMILES bonds to exactly its carrier;
/// the atom leaves the structure and the carrier's D/T counts take it.
fn fold_isotopic_hydrogens(structure: &mut Structure) {
    let carriers: Vec<usize> = structure
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.element == "H" && (atom.isotope == 2 || atom.isotope == 3))
        .map(|(index, _)| index)
        .collect();
    if carriers.is_empty() {
        return;
    }
    let mut removed_bond_positions: Vec<usize> = Vec::new();
    for &hydrogen in &carriers {
        let isotope = structure.atoms[hydrogen].isotope;
        // The carrier: the hydrogen's single bond's other end.
        let Some(bond_position) = structure
            .bonds
            .iter()
            .position(|(a, b, _)| *a == hydrogen || *b == hydrogen)
        else {
            continue;
        };
        removed_bond_positions.push(bond_position);
        let (a, b, _) = structure.bonds[bond_position];
        let carrier = if a == hydrogen { b } else { a };
        let entry = structure.isotopic_h.entry(carrier).or_insert((0u32, 0u32));
        if isotope == 2 {
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
        structure.bonds.remove(bond_position);
    }
    // Remove the hydrogen atoms (descending: the indices stay stable)
    // and reindex the bonds above each removal.
    let removed: std::collections::BTreeSet<usize> = carriers.into_iter().collect();
    let mut kept_atoms = Vec::with_capacity(structure.atoms.len());
    let mut kept_positions = std::mem::take(&mut structure.positions);
    let mut remap = std::collections::HashMap::new();
    for (index, atom) in std::mem::take(&mut structure.atoms).into_iter().enumerate() {
        if removed.contains(&index) {
            continue;
        }
        remap.insert(index, kept_atoms.len());
        kept_positions.push(kept_positions.get(index).copied().unwrap_or((0.0, 0.0)));
        kept_atoms.push(atom);
    }
    structure.atoms = kept_atoms;
    structure.positions = kept_positions;
    structure.bonds = structure
        .bonds
        .drain(..)
        .map(|(a, b, kind)| (remap[&a], remap[&b], kind))
        .collect();
    structure.isotopic_h = structure
        .isotopic_h
        .iter()
        .map(|(carrier, counts)| (remap[carrier], *counts))
        .collect();
    // The SMILES stereo markers ride the same renumbering.
    structure.smiles_parities = structure
        .smiles_parities
        .iter()
        .filter_map(|(old, &mark)| remap.get(old).map(|&new| (new, mark)))
        .collect();
    structure.bond_directions = structure
        .bond_directions
        .iter()
        .filter_map(|(old, &dir)| {
            let removed = removed_bond_positions.iter().filter(|&&p| p < *old).count();
            if removed_bond_positions.contains(old) {
                return None;
            }
            Some((old - removed, dir))
        })
        .collect();
}

struct Walk<'a> {
    capture: &'a Capture<'a>,
    structure: Structure,
    /// The atom the next bond attaches to.
    current: Option<usize>,
    pending_bond: Option<u32>,
    /// The directional spelling (`/` = 1, `\` = 2) riding the next
    /// bond — a stereo marker the synthesis turns into cis/trans
    /// geometry when the bond borders a double bond.
    pending_direction: Option<u8>,
    /// The branch-point stack (one entry per open paren).
    stack: Vec<Option<usize>>,
    /// Open ring bonds by digit value: (atom index, explicit bond
    /// kind, explicit direction).
    rings: std::collections::HashMap<u32, (usize, Option<u32>, Option<u8>)>,
    /// Aromatic markers by atom index (the lowercase organic subset).
    aromatic_markers: Vec<bool>,
}

impl<'a> Walk<'a> {
    /// A capture's matched text: the input slice over its span (the
    /// named-capture nodes wrap their children, so their `text()` is
    /// empty — the span is the read).
    fn text_of(&self, node: &parsanol::portable::ast::AstNode) -> &'a str {
        let (start, end) = self.capture.span(node);
        &self.capture.input()[start..end]
    }

    fn chain(&mut self, node: &parsanol::portable::ast::AstNode) -> crate::Result<()> {
        self.dispatch(node)
    }

    /// The formula walk (grammar v1.1.0, chemomer#131): the first
    /// chain, then each dot-separated component — the component
    /// boundary resets the current atom so the next chain's first
    /// atom bonds to nothing (the molfile multi-fragment semantics;
    /// the identity layer's connected-component analysis derives
    /// per-component).
    fn formula(&mut self, node: &parsanol::portable::ast::AstNode) -> crate::Result<()> {
        // The formula walk (grammar v1.1.0, chemomer#131): descend
        // the unnamed sequence/repetition wrappers to each `chain`
        // capture — exactly one per top-level component (branch
        // interiors nest their chains under the distinct `group`
        // key, so they never surface here). The component boundary
        // resets the current atom so the next chain's first atom
        // bonds to nothing (the molfile multi-fragment semantics;
        // the identity layer's connected-component analysis derives
        // per-component).
        let mut first = true;
        self.walk_formula_element(node, &mut first)
    }

    fn walk_formula_element(
        &mut self,
        node: &parsanol::portable::ast::AstNode,
        first: &mut bool,
    ) -> crate::Result<()> {
        match node {
            parsanol::portable::ast::AstNode::Array { .. } => {
                for item in self.capture.array_items(node) {
                    self.walk_formula_element(&item, first)?;
                }
            }
            parsanol::portable::ast::AstNode::Hash { .. } => {
                let items = self.capture.hash_items(node);
                if items.len() == 1 && items[0].0 == "chain" {
                    // A chain capture itself: one top-level component.
                    if !*first {
                        self.current = None;
                    }
                    *first = false;
                    return self.chain(&items[0].1);
                }
                for (_, value) in items {
                    self.walk_formula_element(&value, first)?;
                }
            }
            _ => {}
        }
        Ok(())
    }

    /// Every named descendant of a capture node, in source order —
    /// the atom capture's value is a Hash or the bracket's unnamed
    /// sequence array, and the inner parts ride at various depths.
    fn flatten_named(
        &self,
        node: &parsanol::portable::ast::AstNode,
        into: &mut Vec<(String, parsanol::portable::ast::AstNode)>,
    ) {
        match node {
            parsanol::portable::ast::AstNode::Hash { .. } => {
                for (key, inner) in self.capture.hash_items(node) {
                    into.push((key.clone(), inner.clone()));
                    self.flatten_named(&inner, into);
                }
            }
            parsanol::portable::ast::AstNode::Array { .. } => {
                for inner in self.capture.array_items(node) {
                    self.flatten_named(&inner, into);
                }
            }
            _ => {}
        }
    }

    /// The generic in-order walk: the capture tree preserves source
    /// order, so the flattened sequence is the SMILES token stream.
    fn dispatch(&mut self, node: &parsanol::portable::ast::AstNode) -> crate::Result<()> {
        match node {
            parsanol::portable::ast::AstNode::Hash { .. } => {
                let items = self.capture.hash_items(node);
                let keys: std::collections::BTreeSet<&str> =
                    items.iter().map(|(k, _)| k.as_str()).collect();
                if keys.contains("a") {
                    // The atom capture's VALUE may be a Hash or the
                    // bracket's unnamed sequence array — either way
                    // the inner shape's named parts (sym/iso/h/chg)
                    // are its named descendants; the chain's remaining
                    // items continue in order.
                    let (_, atom_node) = items.iter().find(|(k, _)| k == "a").expect("a");
                    let mut atom_items = Vec::new();
                    self.flatten_named(atom_node, &mut atom_items);
                    self.atom(&atom_items)?;
                    for (key, inner) in &items {
                        if key != "a" {
                            self.dispatch(inner)?;
                        }
                    }
                    return Ok(());
                }
                if keys.contains("b") {
                    let (_, node) = items.iter().find(|(k, _)| k == "b").expect("b");
                    let symbol = self.text_of(node);
                    let (kind, direction) = bond_kind(symbol)?;
                    self.pending_bond = Some(kind);
                    self.pending_direction = direction;
                    // The remaining items (the atom or ringbond the
                    // bond prefixes) continue in order.
                    for (key, inner) in &items {
                        if key != "b" {
                            self.dispatch(inner)?;
                        }
                    }
                    return Ok(());
                }
                if keys.contains("r") {
                    let (_, node) = items.iter().find(|(k, _)| k == "r").expect("r");
                    let digit = self.text_of(node);
                    self.ring(digit)?;
                    return Ok(());
                }
                if keys.contains("group") {
                    // The branch: the chain inside bonds to the atom
                    // at the branch point; the outer chain resumes
                    // from that same atom afterwards.
                    self.stack.push(self.current);
                    for (key, inner) in &items {
                        if key == "b" {
                            let symbol = self.text_of(inner);
                            let (kind, direction) = bond_kind(symbol)?;
                            self.pending_bond = Some(kind);
                            self.pending_direction = direction;
                        } else if key != "branch" {
                            self.dispatch(inner)?;
                        }
                    }
                    self.current = self.stack.pop().expect("branch stack");
                    return Ok(());
                }
                for (_, inner) in &items {
                    self.dispatch(inner)?;
                }
                Ok(())
            }
            parsanol::portable::ast::AstNode::Array { .. } => {
                for inner in self.capture.array_items(node) {
                    if self.capture.is_empty_marker(&inner) {
                        continue;
                    }
                    self.dispatch(&inner)?;
                }
                Ok(())
            }
            _ => Ok(()),
        }
    }

    fn atom(&mut self, items: &[(String, parsanol::portable::ast::AstNode)]) -> crate::Result<()> {
        let mut symbol = "";
        let mut isotope = 0u16;
        let mut charge = 0i32;
        let mut h_count = 0u32;
        let mut h_pinned = false;
        let mut parity = 0u8;
        for (key, inner) in items {
            let text = self.text_of(inner);
            match key.as_str() {
                "sym" => symbol = text,
                "iso" => {
                    isotope = text
                        .parse()
                        .map_err(|_| Error::new(format!("smiles: bad isotope {text:?}")))?;
                }
                "h" => {
                    let count = text.strip_prefix('H').expect("h prefix");
                    h_count = if count.is_empty() {
                        1
                    } else {
                        count
                            .parse()
                            .map_err(|_| Error::new(format!("smiles: bad H count {text:?}")))?
                    };
                    h_pinned = true;
                }
                "chg" => {
                    charge = parse_charge(text)?;
                }
                "ch" => {
                    parity = parse_chirality(text)?;
                }
                _ => {}
            }
        }
        let (element, aromatic) = element_of(symbol)?;
        let index = self.structure.atoms.len();
        self.structure.atoms.push(AtomRecord {
            element: element.to_string(),
            h_count,
            charge,
            valence_code: 0,
            h_pinned,
            isotope,
        });
        // The pipeline indexes positions per atom (the component
        // split); SMILES carries no 2D coordinates — zeros.
        self.structure.positions.push((0.0, 0.0));
        // Bond to the current atom (implicit: aromatic between two
        // aromatic atoms, single otherwise). The bond carries the
        // pending explicit kind and the pending direction; the
        // direction rides the bond index for the synthesis.
        if let Some(previous) = self.current {
            let kind = self.pending_bond.take().unwrap_or_else(|| {
                if aromatic && self.aromatic_at(previous) {
                    BOND_AROMATIC
                } else {
                    BOND_SINGLE
                }
            });
            let bond_index = self.structure.bonds.len();
            self.structure.bonds.push((previous, index, kind));
            if let Some(direction) = self.pending_direction.take() {
                self.structure.bond_directions.insert(bond_index, direction);
            }
        }
        if parity != 0 {
            self.structure.smiles_parities.insert(index, parity);
        }
        self.aromatic_markers.push(aromatic);
        self.current = Some(index);
        Ok(())
    }

    fn ring(&mut self, digit: &str) -> crate::Result<()> {
        let value: u32 = if let Some(rest) = digit.strip_prefix('%') {
            rest.parse().map_err(|_| bad_ring(digit))?
        } else {
            digit.parse().map_err(|_| bad_ring(digit))?
        };
        let explicit = self.pending_bond.take();
        let direction = self.pending_direction.take();
        let Some(current) = self.current else {
            return err("smiles: ring closure before any atom");
        };
        match self.rings.remove(&value) {
            Some((open_atom, open_kind, open_direction)) => {
                let kind = explicit.or(open_kind).unwrap_or_else(|| {
                    if self.aromatic_at(open_atom) && self.aromatic_at(current) {
                        BOND_AROMATIC
                    } else {
                        BOND_SINGLE
                    }
                });
                let bond_index = self.structure.bonds.len();
                self.structure.bonds.push((open_atom, current, kind));
                if let Some(direction) = direction.or(open_direction) {
                    self.structure.bond_directions.insert(bond_index, direction);
                }
            }
            None => {
                self.rings.insert(value, (current, explicit, direction));
            }
        }
        Ok(())
    }

    fn aromatic_at(&self, index: usize) -> bool {
        self.aromatic_markers.get(index).copied().unwrap_or(false)
    }
}

fn bad_ring(digit: &str) -> Error {
    Error::new(format!("smiles: bad ring id {digit:?}"))
}

/// The bond symbol → the V2000 code plus the directional spelling
/// (`/` = 1, `\` = 2) when it carries one. The direction is the
/// cis/trans marker the synthesis turns into drawn geometry when the
/// bond borders a double bond; away from one it is inert.
fn bond_kind(symbol: &str) -> crate::Result<(u32, Option<u8>)> {
    match symbol {
        "-" | "~" => Ok((BOND_SINGLE, None)),
        ":" => Ok((BOND_AROMATIC, None)),
        "=" => Ok((2, None)),
        "#" | "$" => Ok((3, None)),
        "/" => Ok((BOND_SINGLE, Some(1))),
        "\\" => Ok((BOND_SINGLE, Some(2))),
        other => err(format!("smiles: unknown bond symbol {other:?}")),
    }
}

/// The bracket chirality → the parity mark (1 = `@`/`@TH1`, 2 =
/// `@@`/`@TH2`). The other tetrahedral spellings (@AL/@SP/@TB classes
/// and their offsets) reject with the note — the corpus carries only
/// the two common spellings.
fn parse_chirality(text: &str) -> crate::Result<u8> {
    match text {
        "@" | "@TH1" => Ok(1),
        "@@" | "@TH2" => Ok(2),
        other => err(format!(
            "smiles: chirality spelling {other:?} is not part of the \
             phase-2 contract (@/@TH1/@TH2/@@ only)"
        )),
    }
}

/// Turns the SMILES stereo markers into drawn geometry: the wedge
/// perception reads `stereo_bonds` + `positions`, and the cis/trans
/// side product reads the reference ligands' positions — so the
/// synthesis lays out the parity atoms' ligands trigonally with a
/// wedge realizing the @/@@ hand, and lifts each directional bond's
/// substituent to its side of the adjacent double bond. Everything
/// else stays on the string-order axis; only the marked centers'
/// ligand triangles and the marked bonds' substituents carry
/// geometry. Both marker maps clear — the structure then looks to
/// the pipeline exactly like a drawn molfile.
fn synthesize_stereo(structure: &mut Structure) -> crate::Result<()> {
    if structure.smiles_parities.is_empty() && structure.bond_directions.is_empty() {
        return Ok(());
    }
    // The default layout: the string order as a straight axis.
    for index in 0..structure.positions.len() {
        structure.positions[index] = (index as f64, 0.0);
    }

    // cis/trans: each directional bond lifts its substituent to the
    // side its spelling draws — `/` rises left→right in string order,
    // so the substituent sits below when it precedes the double-bond
    // endpoint and above when it follows; `\` mirrors both.
    let directions = structure.bond_directions.clone();
    for (bond_index, direction) in &directions {
        let Some(&(a, b, _)) = structure.bonds.get(*bond_index) else {
            continue;
        };
        let Some((endpoint, substituent)) = double_bond_endpoint(structure, a, b) else {
            continue;
        };
        let rise: f64 = if *direction == 1 { 1.0 } else { -1.0 };
        let order: f64 = if substituent > endpoint { 1.0 } else { -1.0 };
        let (x, _) = structure.positions[substituent];
        structure.positions[substituent] = (x, rise * order);
    }

    // Tetrahedral: each parity atom's heavy ligands (atom-index
    // order — the SMILES bond order the parity is defined over) sit
    // in a trigonal triangle around the center, and a dash-wedge on
    // the first ligand's bond realizes the hand. The wedge bond is
    // rewritten center-first: the perception reads the code's first
    // atom as the center.
    let parities = structure.smiles_parities.clone();
    for (center, mark) in &parities {
        let ligands: Vec<(usize, usize)> = structure
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(index, &(a, b, _))| {
                if a == *center {
                    Some((index, b))
                } else if b == *center {
                    Some((index, a))
                } else {
                    None
                }
            })
            .collect();
        if ligands.len() < 2 {
            return err(format!(
                "smiles: chirality on atom {center} with fewer than two bonds"
            ));
        }
        let Some(wedge_index) = ligands.iter().find(|(_, n)| n > center).map(|&(i, _)| i) else {
            return err(format!(
                "smiles: chirality on atom {center} whose bonds all close earlier atoms"
            ));
        };
        // The three heavy ligands in neighbor order; a lone heavy
        // pair (CH2-like) is not a stereocenter — the perception's
        // own gates decline it, the synthesis still draws the wedge.
        let neighbors: Vec<usize> = ligands.iter().map(|&(_, n)| n).collect();
        let (cx, _) = structure.positions[*center];
        structure.positions[*center] = (cx, 0.0);
        let angles = [90.0f64, 210.0, 330.0];
        for (slot, &neighbor) in neighbors.iter().take(3).enumerate() {
            let radians = angles[slot].to_radians();
            structure.positions[neighbor] = (cx + radians.cos(), radians.sin());
        }
        let down: u8 = if *mark == 1 { 6 } else { 1 };
        // The wedge code: 1 = up (solid), 6 = dash. The hand flips
        // with the code; @/@@ are enantiomers, so one constant maps
        // the pair — corpus-fitted.
        structure.stereo_bonds.insert(wedge_index, down);
        let (from, to, kind) = structure.bonds[wedge_index];
        if from != *center {
            structure.bonds[wedge_index] = (*center, to, kind);
        }
    }
    structure.smiles_parities.clear();
    structure.bond_directions.clear();
    Ok(())
}

/// The double-bond endpoint this single bond borders, with the bond's
/// substituent (the other end): `Some((axis endpoint, substituent))`
/// when one end of the single bond is a kind-2 bond's atom.
fn double_bond_endpoint(structure: &Structure, a: usize, b: usize) -> Option<(usize, usize)> {
    let is_double_partner = |atom: usize, other: usize| {
        structure.bonds.iter().any(|&(from, to, kind)| {
            kind == 2 && ((from == atom && to != other) || (to == atom && from != other))
        })
    };
    if is_double_partner(a, b) {
        Some((a, b))
    } else if is_double_partner(b, a) {
        Some((b, a))
    } else {
        None
    }
}

fn parse_charge(text: &str) -> crate::Result<i32> {
    let (sign, rest) = match text.chars().next() {
        Some('+') => (1i32, &text[1..]),
        Some('-') => (-1i32, &text[1..]),
        _ => return err(format!("smiles: bad charge {text:?}")),
    };
    let magnitude = if rest.is_empty() {
        1
    } else if rest.chars().all(|c| c == '+' || c == '-') {
        // `++` / `--` spellings: each repeat adds one.
        rest.len() as i32 + 1
    } else {
        rest.parse::<i32>()
            .map_err(|_| Error::new(format!("smiles: bad charge {text:?}")))?
    };
    Ok(sign * magnitude)
}

/// The element symbol: aromatic lowercase markers read as their
/// uppercase elements with the aromatic flag; `*` is the day-1
/// wildcard rejection (the pseudo-element carrier is molfile-side).
fn element_of(symbol: &str) -> crate::Result<(String, bool)> {
    // Aromatic markers read as their uppercase elements with the flag;
    // everything else passes through as spelled (the engine's
    // periodic-table validation rejects non-elements downstream).
    let (element, aromatic) = match symbol {
        "c" => ("C", true),
        "n" => ("N", true),
        "o" => ("O", true),
        "p" => ("P", true),
        "s" => ("S", true),
        "*" => return err("smiles: the wildcard atom (*) is not part of the day-1 contract"),
        other => (other, false),
    };
    Ok((element.to_string(), aromatic))
}

// -- writing (chemomer#130) ----------------------------------------------

/// The organic-subset elements written bare (no bracket) and the
/// aromatic markers written lowercase.
const ORGANIC_UPPER: [&str; 10] = ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"];
const AROMATIC_LOWER: [&str; 8] = ["b", "c", "n", "o", "p", "s", "se", "as"];
const CLOSURE_DIGITS: [&str; 9] = ["1", "2", "3", "4", "5", "6", "7", "8", "9"];

/// Writes a [`Structure`] as a deterministic canonical SMILES: the
/// atoms rank by a Morgan-style refinement (element, charge, isotope,
/// explicit-H count, degree — neighbours' ranks folded in until
/// stable), the DFS spanning tree visits children in rank order, ring
/// bonds carry digits with kind carriers, and bond symbols appear
/// only where the bond departs from its default. Stability is the
/// contract: `write(read(s))` re-reads to itself, and the written
/// form derives the recorded identity.
pub fn write(structure: &Structure) -> Result<String> {
    let expanded = expand_isotopic_h(structure);
    let mut written = String::new();
    for (index, component) in expanded.split_components().iter().enumerate() {
        if index > 0 {
            written.push('.');
        }
        written.push_str(&write_component(component)?);
    }
    Ok(written)
}

/// The inverse of the reader's isotopic-hydrogen fold: carriers'
/// recorded D/T counts re-materialize as explicit `[2H]`/`[3H]`
/// atoms bonded to the carrier, so the written form re-reads to the
/// same folded structure. The plain-H accounting is untouched (the
/// /h layer is isotope-blind; valence carries the gap).
fn expand_isotopic_h(structure: &Structure) -> Structure {
    if structure.isotopic_h.is_empty() {
        return structure.clone();
    }
    let mut expanded = structure.clone();
    let additions: Vec<(usize, u16)> = structure
        .isotopic_h
        .iter()
        .flat_map(|(&carrier, &(d_count, t_count))| {
            std::iter::repeat_n((carrier, 2u16), d_count as usize)
                .chain(std::iter::repeat_n((carrier, 3u16), t_count as usize))
        })
        .collect();
    for &(carrier, isotope) in &additions {
        let index = expanded.atoms.len();
        expanded.atoms.push(AtomRecord {
            element: "H".to_string(),
            h_count: 0,
            charge: 0,
            valence_code: 0,
            h_pinned: false,
            isotope,
        });
        expanded.bonds.push((carrier, index, 1));
        expanded.ride_marks.push(0);
        expanded.positions.push((0.0, 0.0));
    }
    expanded
}

struct Writer<'a> {
    component: &'a Structure,
    adjacency: Vec<Vec<(usize, usize)>>,
    tree_children: Vec<Vec<usize>>,
    digit_for_bond: Vec<Option<&'static str>>,
}

/// The atom-index tie-breaks cannot see which ring bond carries the
/// double (Kekulé-rotation patterns ride the canonicalizer's
/// symmetry ties), so the spanning-tree walk ENUMERATES: at every
/// node the minimum-rank tied children are tried as the first
/// branch, each alternative carrying a cloned state, and the
/// lexicographically smallest complete spelling wins. The minimum is
/// a function of the unordered structure alone — the stability the
/// contract asks for. The budget bounds the alternatives; past it
/// the rank order stands (deterministic, index-based).
const ENUMERATION_CAP: i64 = 32;

#[derive(Clone)]
struct SearchState {
    visited: Vec<bool>,
    parent: Vec<usize>,
    tree_children: Vec<Vec<usize>>,
    /// The ring bonds in discovery order — the digit table's order.
    ring_bonds: Vec<usize>,
    ring_seen: std::collections::BTreeSet<usize>,
}

impl SearchState {
    fn fresh(root: usize, atom_count: usize) -> Self {
        let mut visited = vec![false; atom_count];
        visited[root] = true;
        SearchState {
            visited,
            parent: vec![usize::MAX; atom_count],
            tree_children: vec![Vec::new(); atom_count],
            ring_bonds: Vec::new(),
            ring_seen: std::collections::BTreeSet::new(),
        }
    }

    fn digit_for_bond(&self, bond_count: usize) -> Result<Vec<Option<&'static str>>> {
        let mut table: Vec<Option<&'static str>> = vec![None; bond_count];
        for (slot, &bond) in self.ring_bonds.iter().enumerate() {
            let Some(digit) = CLOSURE_DIGITS.get(slot) else {
                return err("smiles: more than nine ring bonds in one component");
            };
            table[bond] = Some(digit);
        }
        Ok(table)
    }
}

fn write_component(component: &Structure) -> Result<String> {
    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); component.atoms.len()];
    for (bond_index, &(from, to, _)) in component.bonds.iter().enumerate() {
        adjacency[from].push((to, bond_index));
        adjacency[to].push((from, bond_index));
    }

    let ranks = crate::layers::hydrogen::canonical_ranks_for(component);
    if component.atoms.is_empty() {
        return err("smiles: empty component");
    }
    let top_rank = ranks.iter().copied().min().unwrap_or(0);
    let mut best: Option<(String, SearchState)> = None;
    for root in (0..component.atoms.len()).filter(|&i| ranks[i] == top_rank) {
        let state = SearchState::fresh(root, component.atoms.len());
        let candidate = best_subtree(component, &adjacency, &ranks, root, state, ENUMERATION_CAP)?;
        if best
            .as_ref()
            .is_none_or(|(current, _)| candidate.0 < *current)
        {
            best = Some(candidate);
        }
    }
    let Some((_, final_state)) = best else {
        return err("smiles: empty component");
    };
    let root = final_state
        .parent
        .iter()
        .position(|&p| p == usize::MAX)
        .ok_or_else(|| crate::Error::new("smiles: no root".to_string()))?;
    let digit_for_bond = final_state.digit_for_bond(component.bonds.len())?;
    let mut writer = Writer {
        component,
        adjacency,
        tree_children: final_state.tree_children,
        digit_for_bond,
    };
    let mut out = String::new();
    writer.emit(root, &mut out)?;
    Ok(out)
}

/// The enumerating tree search: returns the subtree's spelling (with
/// the subtree's own ring digits) and the completed state. Children
/// record in walk order — tree branches recurse, ring bonds join the
/// discovery order that numbers the digits.
fn best_subtree(
    component: &Structure,
    adjacency: &[Vec<(usize, usize)>],
    ranks: &[u32],
    atom: usize,
    state: SearchState,
    budget: i64,
) -> Result<(String, SearchState)> {
    let mut children: Vec<usize> = adjacency[atom]
        .iter()
        .map(|&(n, _)| n)
        .filter(|&n| n != state.parent[atom])
        .collect();
    children.sort_by_key(|&n| (ranks[n], n));
    children.dedup();

    // The minimum-rank tied children are the alternatives' variable
    // (each tried as the first branch); the rest follow rank order.
    let mut alternatives: Vec<Vec<usize>> = Vec::new();
    if children.len() >= 2 && budget > 0 {
        let min_rank = ranks[children[0]];
        let tied: Vec<usize> = children
            .iter()
            .copied()
            .take_while(|&n| ranks[n] == min_rank)
            .collect();
        // Only ties whose children carry DIFFERENT bond signatures to
        // this atom can change ring closure assignments — the ring
        // entry walking into a kekulized ring. Equal-signature ties
        // (symmetric substituents, aromaticized rings) spell
        // identically whichever way they order, and enumerating them
        // multiplies the search without effect.
        let signature = |child: usize| -> String {
            let bond = adjacency[atom]
                .iter()
                .find(|&&(n, _)| n == child)
                .map(|&(_, b)| b)
                .expect("adjacent bond");
            let (_, _, kind) = component.bonds[bond];
            format!("{kind}")
        };
        if tied.len() >= 2 && tied.len() <= 3 {
            let signatures: std::collections::BTreeSet<String> =
                tied.iter().map(|&child| signature(child)).collect();
            if signatures.len() > 1 {
                for &first in &tied {
                    let mut ordered: Vec<usize> = vec![first];
                    ordered.extend(children.iter().copied().filter(|&n| n != first));
                    alternatives.push(ordered);
                }
            }
        }
    }
    if alternatives.is_empty() {
        alternatives.push(children);
    }

    let mut best: Option<(String, SearchState)> = None;
    for (attempt, order) in alternatives.iter().enumerate() {
        let budget_left = if alternatives.len() > 1 {
            budget - 1
        } else {
            budget
        };
        let mut branch = state.clone();
        let mut pieces: Vec<(usize, String)> = Vec::new();
        for &child in order {
            let symbol = bond_symbol_between(component, adjacency, atom, child);
            let mut child_state = branch.clone();
            match record_child(&mut child_state, adjacency, atom, child) {
                None => {
                    // Ring closure discovered at this child: the digit
                    // joins this atom's token; no subtree piece.
                    branch = child_state;
                }
                Some(tree_child) => {
                    let (piece, done) = best_subtree(
                        component,
                        adjacency,
                        ranks,
                        tree_child,
                        child_state,
                        budget_left,
                    )?;
                    pieces.push((child, format!("{symbol}{piece}")));
                    branch = done;
                }
            }
        }
        let token = atom_token_with(component, adjacency, &branch, atom)?;
        let mut assembled = token;
        let recorded = branch.tree_children[atom].clone();
        let last = recorded.len().saturating_sub(1);
        for (position, &child) in recorded.iter().enumerate() {
            let symbol = bond_symbol_between(component, adjacency, atom, child);
            let piece = child_results_lookup(&pieces, child);
            if position == last {
                assembled.push_str(&symbol);
                assembled.push_str(&piece);
            } else {
                assembled.push('(');
                assembled.push_str(&symbol);
                assembled.push_str(&piece);
                assembled.push(')');
            }
        }
        let _ = attempt;
        let candidate = (assembled, branch);
        if best
            .as_ref()
            .is_none_or(|(current, _)| candidate.0 < *current)
        {
            best = Some(candidate);
        }
    }
    best.ok_or_else(|| crate::Error::new("smiles: no walk".to_string()))
}

fn child_results_lookup(pieces: &[(usize, String)], child: usize) -> String {
    pieces
        .iter()
        .find(|(c, _)| *c == child)
        .map(|(_, p)| p.clone())
        .unwrap_or_default()
}

fn record_child(
    state: &mut SearchState,
    adjacency: &[Vec<(usize, usize)>],
    atom: usize,
    child: usize,
) -> Option<usize> {
    if state.visited[child] {
        let bond = adjacency[atom]
            .iter()
            .find(|&&(n, _)| n == child)
            .map(|&(_, b)| b)
            .expect("adjacent bond");
        if state.ring_seen.insert(bond) {
            state.ring_bonds.push(bond);
        }
        None
    } else {
        state.visited[child] = true;
        state.parent[child] = atom;
        state.tree_children[atom].push(child);
        Some(child)
    }
}

fn bond_symbol_between(
    component: &Structure,
    adjacency: &[Vec<(usize, usize)>],
    from: usize,
    to: usize,
) -> String {
    let bond_index = adjacency[from]
        .iter()
        .find(|&&(n, _)| n == to)
        .map(|&(_, b)| b)
        .expect("adjacent bond");
    let (_, _, kind) = component.bonds[bond_index];
    match kind {
        1 => {
            if is_aromatic_at(component, adjacency, from)
                || is_aromatic_at(component, adjacency, to)
            {
                "-".to_string()
            } else {
                String::new()
            }
        }
        2 => "=".to_string(),
        3 => "#".to_string(),
        4 => {
            if is_aromatic_at(component, adjacency, from)
                && is_aromatic_at(component, adjacency, to)
            {
                String::new()
            } else {
                "-".to_string()
            }
        }
        _ => "-".to_string(),
    }
}

fn is_aromatic_at(component: &Structure, adjacency: &[Vec<(usize, usize)>], atom: usize) -> bool {
    adjacency[atom].iter().any(|&(_, bond_index)| {
        let (_, _, kind) = component.bonds[bond_index];
        kind == 4
    })
}

/// The atom token, reading ring digits off the state's discovery
/// order. A closure discovered HERE (this atom is the second
/// endpoint) and an closure whose opener sits elsewhere both carry
/// the digit at this token; kind carriers ride non-default bonds.
fn atom_token_with(
    component: &Structure,
    adjacency: &[Vec<(usize, usize)>],
    state: &SearchState,
    atom: usize,
) -> Result<String> {
    let record = &component.atoms[atom];
    let aromatic = is_aromatic_at(component, adjacency, atom);
    let organic = ORGANIC_UPPER.contains(&record.element.as_str())
        || (aromatic && AROMATIC_LOWER.contains(&record.element.to_lowercase().as_str()));
    let bracket_needed = !organic
        || record.charge != 0
        || record.isotope != 0
        || record.h_count > 0
        || (record.h_count == 0 && record.h_pinned);

    let mut token = String::new();
    if bracket_needed {
        token.push('[');
        if record.isotope != 0 {
            token.push_str(&record.isotope.to_string());
        }
        token.push_str(&record.element);
        if record.h_count > 0 {
            token.push('H');
            if record.h_count > 1 {
                token.push_str(&record.h_count.to_string());
            }
        } else if record.h_pinned {
            token.push_str("H0");
        }
        if record.charge != 0 {
            let magnitude = record.charge.abs();
            token.push(if record.charge > 0 { '+' } else { '-' });
            if magnitude > 1 {
                token.push_str(&magnitude.to_string());
            }
        }
        token.push(']');
    } else if aromatic {
        token.push_str(&record.element.to_lowercase());
    } else {
        token.push_str(&record.element);
    }

    let digit_table = state.digit_for_bond(component.bonds.len())?;
    let mut carried: Vec<usize> = adjacency[atom]
        .iter()
        .map(|&(_, bond_index)| bond_index)
        .filter(|&bond_index| digit_table[bond_index].is_some())
        .collect();
    carried.sort_by_key(|&bond_index| {
        state
            .ring_bonds
            .iter()
            .position(|&b| b == bond_index)
            .unwrap_or(usize::MAX)
    });
    for bond_index in carried {
        let digit = digit_table[bond_index].unwrap_or("");
        let (from, to, kind) = component.bonds[bond_index];
        let other = if from == atom { to } else { from };
        let both_aromatic = is_aromatic_at(component, adjacency, atom)
            && is_aromatic_at(component, adjacency, other);
        match kind {
            2 => token.push('='),
            3 => token.push('#'),
            1 => {
                if is_aromatic_at(component, adjacency, atom)
                    || is_aromatic_at(component, adjacency, other)
                {
                    token.push('-');
                }
            }
            4 => {
                if !both_aromatic {
                    token.push('-');
                }
            }
            _ => token.push('-'),
        }
        token.push_str(digit);
    }
    Ok(token)
}

impl<'a> Writer<'a> {
    /// Depth-first emission: the atom token, then all children but
    /// the last as branches (rank order), then the last as the
    /// continuation.
    fn emit(&mut self, atom: usize, out: &mut String) -> Result<()> {
        out.push_str(&self.atom_token(atom)?);
        let children = self.tree_children[atom].clone();
        let last_child = children.len().saturating_sub(1);
        for (position, &child) in children.iter().enumerate() {
            let bond_index = self.adjacency[atom]
                .iter()
                .find(|&&(n, _)| n == child)
                .map(|&(_, b)| b)
                .expect("tree bond");
            let symbol = self.bond_symbol(bond_index);
            if position == last_child {
                out.push_str(&symbol);
                self.emit(child, out)?;
            } else {
                out.push('(');
                out.push_str(&symbol);
                self.emit(child, out)?;
                out.push(')');
            }
        }
        Ok(())
    }

    /// The bond symbol between tree parent and child — omitted for
    /// the defaults (single between non-aromatic atoms, aromatic
    /// between two aromatic atoms).
    fn bond_symbol(&self, bond_index: usize) -> String {
        let (from, to, kind) = self.component.bonds[bond_index];
        match kind {
            1 => {
                if self.is_aromatic(from) || self.is_aromatic(to) {
                    "-".to_string()
                } else {
                    String::new()
                }
            }
            2 => "=".to_string(),
            3 => "#".to_string(),
            4 => {
                if self.is_aromatic(from) && self.is_aromatic(to) {
                    String::new()
                } else {
                    "-".to_string()
                }
            }
            _ => "-".to_string(),
        }
    }

    fn is_aromatic(&self, atom: usize) -> bool {
        self.adjacency[atom].iter().any(|&(_, bond_index)| {
            let (_, _, kind) = self.component.bonds[bond_index];
            kind == 4
        })
    }

    fn atom_token(&self, atom: usize) -> Result<String> {
        let record = &self.component.atoms[atom];
        let aromatic = self.is_aromatic(atom);
        let organic = ORGANIC_UPPER.contains(&record.element.as_str())
            || (aromatic && AROMATIC_LOWER.contains(&record.element.to_lowercase().as_str()));
        let bracket_needed = !organic
            || record.charge != 0
            || record.isotope != 0
            || record.h_count > 0
            || (record.h_count == 0 && record.h_pinned);

        let mut token = String::new();
        if bracket_needed {
            token.push('[');
            if record.isotope != 0 {
                token.push_str(&record.isotope.to_string());
            }
            token.push_str(&record.element);
            if record.h_count > 0 {
                token.push('H');
                if record.h_count > 1 {
                    token.push_str(&record.h_count.to_string());
                }
            } else if record.h_pinned {
                token.push_str("H0");
            }
            if record.charge != 0 {
                let magnitude = record.charge.abs();
                token.push(if record.charge > 0 { '+' } else { '-' });
                if magnitude > 1 {
                    token.push_str(&magnitude.to_string());
                }
            }
            token.push(']');
        } else if aromatic {
            token.push_str(&record.element.to_lowercase());
        } else {
            token.push_str(&record.element);
        }

        // Ring-closure digits at this atom, with kind carriers for
        // non-default closure bonds.
        for &(_, bond_index) in &self.adjacency[atom] {
            let Some(digit) = self.digit_for_bond[bond_index] else {
                continue;
            };
            let (from, to, kind) = self.component.bonds[bond_index];
            let other = if from == atom { to } else { from };
            let (atom_aromatic, other_aromatic) = (self.is_aromatic(atom), self.is_aromatic(other));
            let carrier = match kind {
                2 => Some('='),
                3 => Some('#'),
                // Single-or-aromatic closures carry a dash whenever
                // they depart from their default (aromatic only
                // between two aromatic atoms).
                1 => (!other_aromatic || atom_aromatic).then_some('-'),
                4 => (!(atom_aromatic && other_aromatic)).then_some('-'),
                _ => Some('-'),
            };
            if let Some(symbol) = carrier {
                token.push(symbol);
            }
            token.push_str(digit);
        }
        Ok(token)
    }
}
