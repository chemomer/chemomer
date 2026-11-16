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
//! Phase 2 (captured, rejected here until the parity carrier lands):
//! tetrahedral `@`/`@@` markers and the directional `/` `\` bonds —
//! the recorded corpus cases carrying them are excluded in the gate
//! with the phase note. Dot-disconnected components reject too (the
//! corpus carries none; the contract phase owns them).

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
            stack: Vec::new(),
            rings: std::collections::HashMap::new(),
            aromatic_markers: Vec::new(),
        };
        walk.chain(capture.root())?;
        if let Some(digit) = walk.rings.keys().next() {
            return err(format!("smiles: ring bond {digit} never closed"));
        }
        if walk.structure.atoms.is_empty() {
            return err("smiles: no atoms");
        }
        fold_isotopic_hydrogens(&mut walk.structure);
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
}

struct Walk<'a> {
    capture: &'a Capture<'a>,
    structure: Structure,
    /// The atom the next bond attaches to.
    current: Option<usize>,
    pending_bond: Option<u32>,
    /// The branch-point stack (one entry per open paren).
    stack: Vec<Option<usize>>,
    /// Open ring bonds by digit value: (atom index, explicit bond kind).
    rings: std::collections::HashMap<u32, (usize, Option<u32>)>,
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
                    self.pending_bond = Some(bond_kind(symbol)?);
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
                            self.pending_bond = Some(bond_kind(symbol)?);
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
                    return err(format!(
                        "smiles: tetrahedral markers ({text}) are the contract's \
                         phase 2 — this reader does not carry parity yet"
                    ));
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
        // aromatic atoms, single otherwise).
        if let Some(previous) = self.current {
            let kind = self.pending_bond.take().unwrap_or_else(|| {
                if aromatic && self.aromatic_at(previous) {
                    BOND_AROMATIC
                } else {
                    BOND_SINGLE
                }
            });
            self.structure.bonds.push((previous, index, kind));
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
        let Some(current) = self.current else {
            return err("smiles: ring closure before any atom");
        };
        match self.rings.remove(&value) {
            Some((open_atom, open_kind)) => {
                let kind = explicit.or(open_kind).unwrap_or_else(|| {
                    if self.aromatic_at(open_atom) && self.aromatic_at(current) {
                        BOND_AROMATIC
                    } else {
                        BOND_SINGLE
                    }
                });
                self.structure.bonds.push((open_atom, current, kind));
            }
            None => {
                self.rings.insert(value, (current, explicit));
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

/// The bond symbol → the V2000 code. Directional `/` `\` and `~` ride
/// phase 2's stereo semantics; day-1 rejects them rather than reading
/// them as single bonds.
fn bond_kind(symbol: &str) -> crate::Result<u32> {
    match symbol {
        "-" | "~" => Ok(BOND_SINGLE),
        ":" => Ok(BOND_AROMATIC),
        "=" => Ok(2),
        "#" => Ok(3),
        "$" => Ok(3),
        "/" | "\\" => err(format!(
            "smiles: directional bond ({symbol}) is the contract's phase 2 — \
             this reader does not carry double-bond stereo yet"
        )),
        other => err(format!("smiles: unknown bond symbol {other:?}")),
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
