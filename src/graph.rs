//! The molecular graph the layers derive from: atoms (element,
//! explicit H count, charge) plus bond orders, with the graph-level
//! queries layers need. Implicit-hydrogen and valence conventions
//! live here because they are graph semantics, not layer semantics.

use crate::error::{err, Result};

#[derive(Debug, Clone, PartialEq)]
pub struct AtomRecord {
    pub element: String,
    /// Explicit H count: MDL code n means n−1 hydrogens (0 = infer
    /// from the valence gap).
    pub h_count: u32,
    pub charge: i32,
    /// The hydrogen count is final as read: the internal-salt
    /// annihilation pins the neutralized anion at its anionic reading
    /// (zero implicit H) — the reference never re-derives H after
    /// moving a charge.
    pub h_pinned: bool,
    /// Recorded mass number (M ISO); zero when the element's standard
    /// isotope. Hydrogen isotopes 2/3 are D/T, never heavy deltas.
    pub isotope: u16,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Structure {
    pub atoms: Vec<AtomRecord>,
    /// V2000 bond types: 1 single, 2 double, 3 triple, 4 aromatic
    /// (aromatic counts as 1.5 in valence sums).
    pub bonds: Vec<(usize, usize, u32)>,
    /// Bond-block stereo codes by bond index (1 = wedge, 6 = dash),
    /// sparse — the tetrahedral layer's input.
    pub stereo_bonds: std::collections::BTreeMap<usize, u8>,
    /// Atom-block (x, y) coordinates — the 2D layout the wedge parity
    /// is read against; empty when the reader carries none.
    pub positions: Vec<(f64, f64)>,
    /// Isotopic hydrogens folded onto their carriers at parse:
    /// carrier index → (D count, T count). The /h layer is
    /// isotope-blind (the D/T add to the carrier's plain count); the
    /// /i layer reads this map.
    pub isotopic_h: std::collections::BTreeMap<usize, (u32, u32)>,
    /// A HARD (DE)PROTONATION phase moved a proton on this structure
    /// (the 5.1b charge-point ride-out or the 5.2b acidic hard
    /// removal) — the reference's FLAG_FORCE_SALT_TAUT: the salt
    /// t-group gates relax to donor-only and the salt merge admits
    /// t-group-carrying candidates without a negative donor.
    pub protons_hard_moved: bool,
    /// Per-bond mobility marks the proton ride left (parallel to
    /// `bonds`, zero = plain): the reference's BNS ride persists the
    /// moved bonds as bond TYPES, not concrete orders — ALT12NS (9,
    /// "1 or 2": either salt role) on single/double-sourced moves,
    /// ALTERN|TRIPLE (7, "2 or 3": neither role) on triple-sourced
    /// ones — while the valence arithmetic reads the moved orders.
    pub ride_marks: Vec<u8>,
    /// The molfile counts-line chiral flag: the structure was marked
    /// chiral at submission (the reference's FLAG_INP_AT_CHIRAL).
    /// The SUCF option reads it and the aux `/rA` letter speaks it.
    pub chiral_flag: bool,
    /// The fold's provenance: for each atom of the folded structure,
    /// its index in the submitted one. Empty when no fold ran (then
    /// the structures align index-for-index) — the aux numbering's
    /// exact original-atom map, hydrogen-bearing originals included.
    pub fold_source: Vec<usize>,
}

/// Allowed valences of the organic subset; the implicit hydrogen
/// count is the smallest allowed valence that absorbs the bond-order
/// sum, so hypervalent S/P chemistry (sulfate's S, phosphate's P)
/// needs no hydrogens where a single fixed valence would invent them.
pub fn allowed_valences(element: &str) -> &'static [f64] {
    match element {
        "C" => &[4.0],
        "N" => &[3.0],
        "O" => &[2.0],
        "S" => &[2.0, 4.0, 6.0],
        "P" => &[3.0, 5.0],
        "F" | "Cl" | "Br" | "I" | "H" => &[1.0],
        _ => &[],
    }
}

/// The valence-gap hydrogen count for one atom.
pub fn valence_gap(element: &str, bond_order_sum: f64) -> u32 {
    for allowed in allowed_valences(element) {
        if *allowed + 1e-9 >= bond_order_sum {
            return (allowed - bond_order_sum).round().max(0.0) as u32;
        }
    }
    0
}

/// The valence-gap hydrogen count for a positively charged atom: the
/// onium table, where the charge itself satisfies the neutral valence
/// (R3S+ sulfonium, R4P+ phosphonium carry no riding hydrogen). Only
/// the group-15/16 onium series realizes its charge as a removable
/// proton; all other cations read gap zero and keep their charge in
/// `/q`.
pub fn cation_valence_gap(element: &str, bond_order_sum: f64) -> u32 {
    let onium: &[f64] = match element {
        "N" | "P" => &[4.0],
        "O" | "S" | "Se" | "Te" => &[3.0],
        _ => return 0,
    };
    for allowed in onium {
        if *allowed + 1e-9 >= bond_order_sum {
            return (allowed - bond_order_sum).round().max(0.0) as u32;
        }
    }
    0
}

pub fn bond_order(kind: u32) -> f64 {
    match kind {
        2 => 2.0,
        3 => 3.0,
        4 => 1.5,
        _ => 1.0,
    }
}

impl Structure {
    /// Folds explicit hydrogen atoms onto their carriers: plain H
    /// and D/T alike add to the carrier's hydrogen count (the /h
    /// layer is isotope-blind), D/T also land in
    /// [`Structure::isotopic_h`] for the /i layer, and the hydrogen
    /// atoms leave the graph — the reference numbers heavy atoms
    /// only.
    pub fn fold_explicit_hydrogens(&mut self) {
        if !self.atoms.iter().any(|atom| atom.element == "H") {
            return;
        }
        let mut isotopic_h = std::mem::take(&mut self.isotopic_h);
        let mut keep = vec![true; self.atoms.len()];
        let mut folds: Vec<(usize, u32, u16)> = Vec::new();
        for (index, atom) in self.atoms.iter().enumerate() {
            if atom.element != "H" {
                continue;
            }
            let carrier = self.bonds.iter().find_map(|&(from, to, _)| {
                if from == index {
                    Some(to)
                } else if to == index {
                    Some(from)
                } else {
                    None
                }
            });
            // A hydrogen survives as a structural atom when it is
            // charged (hydride: InChI=1S/H/q-1), when it has no
            // carrier at all (the valence gap reads it as H2), or
            // when it is an H-H pair's LOWER index — the coalescence
            // TARGET, the higher member folding onto it so the
            // survivor keeps the pair's first original number
            // (InChI=1S/H2/h1H, aux N:1). Otherwise it rides: onto
            // its heavy carrier, or downward within an H-H chain.
            let rides = match carrier {
                None => false,
                Some(carrier) => {
                    atom.charge == 0 && (self.atoms[carrier].element != "H" || carrier < index)
                }
            };
            if !rides {
                continue;
            }
            keep[index] = false;
            folds.push((carrier.unwrap(), 1, atom.isotope));
        }
        for (carrier, count, isotope) in folds {
            self.atoms[carrier].h_count += count;
            match isotope {
                2 => isotopic_h.entry(carrier).or_insert((0, 0)).0 += 1,
                3 => isotopic_h.entry(carrier).or_insert((0, 0)).1 += 1,
                _ => {}
            }
        }
        self.isotopic_h = isotopic_h;
        let mut remap = vec![0usize; self.atoms.len()];
        let mut next = 0usize;
        let mut fold_source = Vec::with_capacity(self.atoms.len());
        for index in 0..self.atoms.len() {
            if keep[index] {
                remap[index] = next;
                next += 1;
                fold_source.push(index);
            }
        }
        self.fold_source = fold_source;
        self.atoms = self
            .atoms
            .iter()
            .enumerate()
            .filter(|(index, _)| keep[*index])
            .map(|(_, atom)| atom.clone())
            .collect();
        self.positions = self
            .positions
            .iter()
            .enumerate()
            .filter(|(index, _)| keep[*index])
            .map(|(_, position)| *position)
            .collect();
        let old_bonds = std::mem::take(&mut self.bonds);
        let old_stereo = std::mem::take(&mut self.stereo_bonds);
        let mut new_stereo = std::collections::BTreeMap::new();
        let mut new_index = 0usize;
        for (bond_index, &(from, to, kind)) in old_bonds.iter().enumerate() {
            if !keep[from] || !keep[to] {
                continue;
            }
            if let Some(&stereo) = old_stereo.get(&bond_index) {
                new_stereo.insert(new_index, stereo);
            }
            self.bonds.push((remap[from], remap[to], kind));
            new_index += 1;
        }
        self.stereo_bonds = new_stereo;
        self.isotopic_h = self
            .isotopic_h
            .iter()
            .map(|(carrier, counts)| (remap[*carrier], *counts))
            .collect();
    }

    pub fn total_charge(&self) -> i32 {
        self.atoms.iter().map(|a| a.charge).sum()
    }

    /// Number of incident bonds per atom (heavy-atom degree).
    pub fn degrees(&self) -> Vec<u32> {
        let mut degrees = vec![0; self.atoms.len()];
        for &(from, to, _) in &self.bonds {
            degrees[from] += 1;
            degrees[to] += 1;
        }
        degrees
    }

    /// Bond-order sums per atom index.
    pub fn bond_order_sums(&self) -> Vec<f64> {
        let mut sums = vec![0.0; self.atoms.len()];
        for &(from, to, kind) in &self.bonds {
            let order = bond_order(kind);
            sums[from] += order;
            sums[to] += order;
        }
        sums
    }

    /// Implicit/explicit hydrogen count per atom: MDL H-count codes
    /// win, otherwise the valence gap (the oracle does the same for
    /// code-0 atoms — a bare `O` line derives H2O). Charged atoms read
    /// their charge-shifted table: anionic heteroatoms take zero
    /// (standardization crosses the `/p` layer), cations read the
    /// onium gap.
    pub fn hydrogens(&self) -> Vec<u32> {
        let sums = self.bond_order_sums();
        self.atoms
            .iter()
            .zip(&sums)
            .map(|(atom, sum)| {
                if atom.h_count > 0 {
                    atom.h_count
                } else if atom.h_pinned || atom.charge < 0 {
                    0
                } else if atom.charge > 0 {
                    cation_valence_gap(&atom.element, *sum)
                } else {
                    valence_gap(&atom.element, *sum)
                }
            })
            .collect()
    }

    /// Element → count, with the hydrogens from [`Self::hydrogens`].
    pub fn composition(&self) -> Vec<(String, u32)> {
        use std::collections::BTreeMap;
        let mut map: BTreeMap<String, u32> = BTreeMap::new();
        for (atom, hydrogens) in self.atoms.iter().zip(self.hydrogens()) {
            *map.entry(atom.element.clone()).or_insert(0) += 1;
            if hydrogens > 0 {
                *map.entry("H".to_string()).or_insert(0) += hydrogens;
            }
        }
        map.into_iter().collect()
    }

    pub fn element(&self, index: usize) -> Result<&str> {
        self.atoms.get(index).map_or_else(
            || err(format!("atom index {index} out of range")),
            |a| Ok(a.element.as_str()),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn structure(atoms: Vec<AtomRecord>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure {
            atoms,
            bonds,
            ..Default::default()
        }
    }

    fn atom(element: &str) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count: 0,
            h_pinned: false,
            isotope: 0,
            charge: 0,
        }
    }

    #[test]
    fn bare_oxygen_composes_water_the_way_the_oracle_does() {
        // A code-0 atom infers H from the valence gap: O with no
        // bonds -> 2 H (recorded identity/breadth water is H2O).
        let s = structure(vec![atom("O")], vec![]);
        assert_eq!(s.composition(), vec![("H".into(), 2), ("O".into(), 1)]);
    }

    #[test]
    fn mdl_hcount_codes_win_over_the_valence_gap() {
        // Ammonium: code 5 -> 4 explicit H (recorded as H3N only
        // after the proton standardization).
        let mut n = atom("N");
        n.h_count = 4;
        n.charge = 1;
        let s = structure(vec![n], vec![]);
        assert_eq!(s.composition(), vec![("H".into(), 4), ("N".into(), 1)]);
    }

    #[test]
    fn anionic_heteroatoms_take_no_implicit_hydrogen() {
        // Acetate: the anionic O carries no H; the mobile H enters via
        // the /p standardization, not here.
        let mut o_anion = atom("O");
        o_anion.charge = -1;
        let s = structure(
            vec![atom("C"), atom("C"), atom("O"), o_anion],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let composition = s.composition();
        assert_eq!(
            composition,
            vec![("C".into(), 2), ("H".into(), 3), ("O".into(), 2)],
            "acetate composition is C2H3O2 before standardization"
        );
    }

    #[test]
    fn aromatic_bonds_count_one_point_five() {
        // Benzene: each C sees 1.5 + 1.5 = 3.0 -> one valence gap -> 1 H.
        let atoms = (0..6).map(|_| atom("C")).collect();
        let bonds = (0..6).map(|i| (i, (i + 1) % 6, 4)).collect();
        let s = structure(atoms, bonds);
        let composition = s.composition();
        assert_eq!(composition, vec![("C".into(), 6), ("H".into(), 6)]);
    }
}
