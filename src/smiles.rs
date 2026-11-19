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

    /// The chain(s) inside a component capture — direct children
    /// only (the dot literal carries no semantics).
    fn component_chains(
        &self,
        node: &parsanol::portable::ast::AstNode,
    ) -> Vec<parsanol::portable::ast::AstNode> {
        let mut chains = Vec::new();
        for (key, value) in self.capture.hash_items(node) {
            if key == "chain" {
                chains.push(value);
            }
        }
        chains
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
