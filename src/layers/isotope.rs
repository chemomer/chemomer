//! The `/i` isotopic layer: heavy-atom mass deltas and isotopic
//! hydrogens, decoded from the oracle battery (2026-10-25). Heavy
//! deltas read `{rank}{±delta}` against the element's standard mass
//! (`i1+2` for 18O, `i1-2` for 125I). Isotopic H rides its carrier's
//! rank — `{rank}{T<n>}{D<n>}` per carrier (`i1D3,2D`), T before D,
//! counts omitted when 1 — except on a heteroatom carried by no
//! carbon (water, ammonia: the exchangeable/mobile class), where the
//! aggregate form `i/hD2` replaces per-carrier entries.

use crate::error::Result;
use crate::graph::Structure;
use crate::layers::standardization;
use crate::layers::{order, LayerEmitter};

pub struct IsotopeLayer;

impl LayerEmitter for IsotopeLayer {
    fn name(&self) -> &'static str {
        "isotope"
    }

    fn order(&self) -> u8 {
        order::I
    }

    fn segment_kind(&self) -> Option<char> {
        Some('i')
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        standardization::ensure_typable_charges(structure)?;
        let mut entries: Vec<(u32, String)> = Vec::new();
        let mut mobile_d = 0u32;
        let mut mobile_t = 0u32;
        for (index, atom) in structure.atoms.iter().enumerate() {
            if let Some(entry) = heavy_entry(atom.isotope, &atom.element) {
                entries.push((index as u32, entry));
            }
        }
        for (carrier, (d, t)) in &structure.isotopic_h {
            if mobile_carrier(structure, *carrier) {
                mobile_d += d;
                mobile_t += t;
            } else {
                entries.push((*carrier as u32, h_entry(*t, *d)));
            }
        }
        if entries.is_empty() && mobile_d == 0 && mobile_t == 0 {
            return Ok(None);
        }
        let ranks = isotope_informed_ranks(structure);
        entries.sort_by_key(|(index, _)| ranks[*index as usize]);
        let mut segment = String::from("i");
        for (index, entry) in entries {
            if segment.len() > 1 {
                segment.push(',');
            }
            segment.push_str(&ranks[index as usize].to_string());
            segment.push_str(&entry);
        }
        if mobile_d > 0 || mobile_t > 0 {
            segment.push_str("/h");
            if mobile_t > 0 {
                segment.push('T');
                if mobile_t > 1 {
                    segment.push_str(&mobile_t.to_string());
                }
            }
            if mobile_d > 0 {
                segment.push('D');
                if mobile_d > 1 {
                    segment.push_str(&mobile_d.to_string());
                }
            }
        }
        Ok(Some(segment))
    }
}

/// The isotopically-informed canonical numbering: the blind numbering
/// with each symmetry class's marked atoms moved to the class's
/// smallest ranks. A marked atom carries a heavy isotope or a fixed
/// isotopic-H load; its sort key reads (mass delta, T, D) ascending —
/// the reference's `iso_sort_key` packing (ichiisot.c
/// `make_iso_sort_key`, H-counts inside the mass delta), unmarked
/// atoms LAST (`EMPTY_ISO_SORT_KEY`). The reference realizes this as
/// a second canonicalization seeded from the base classes with the
/// iso keys as the auxiliary rank (ichican2.c CanonGraph06, the aux
/// ranking built from `nSymmRankBase` + `iso_sort_keyBase`); within a
/// class the base connection table ties, so the iso rows decide, and
/// they minimize with the marked atoms at the class's smallest
/// canonical numbers — identical to reassigning the class's rank
/// slots, marked atoms first. Unmarked members keep their blind
/// relative order.
pub fn isotope_informed_ranks(structure: &Structure) -> Vec<u32> {
    let ranks = crate::layers::hydrogen::canonical_ranks_for(structure);
    let marks: Vec<Option<(i32, u32, u32)>> = (0..structure.atoms.len())
        .map(|index| mark_of(structure, index))
        .collect();
    if marks.iter().all(Option::is_none) {
        return ranks;
    }
    let (adjacency, invariants) = crate::canonical::graph_of(structure);
    let counts = standardization::adjusted_hydrogens(structure);
    let mask = crate::canonical::tautomer::endpoint_mask(structure, &counts);
    let groups = crate::canonical::group_membership(structure, &counts);
    let (_, symmetry) = crate::canonical::canonical_ranks_and_symmetry(
        &adjacency,
        &invariants,
        &crate::canonical::tautomer::without_mobile_h(&counts, &mask),
        &groups,
    );
    let mut informed = ranks.clone();
    // One class at a time: its rank slots (sorted), its members
    // (marked ascending by key, then unmarked by blind rank).
    let mut classes: Vec<Vec<usize>> = Vec::new();
    for atom in 0..structure.atoms.len() {
        let position = classes
            .iter()
            .position(|members| symmetry[members[0]] == symmetry[atom]);
        match position {
            Some(index) => classes[index].push(atom),
            None => classes.push(vec![atom]),
        }
    }
    for members in &classes {
        let mut slots: Vec<u32> = members.iter().map(|&atom| ranks[atom]).collect();
        slots.sort_unstable();
        let mut ordered: Vec<usize> = members.clone();
        // Marked atoms first (ascending iso key — `EMPTY_ISO_SORT_KEY`
        // holds the blind order's tail, never the head).
        ordered.sort_by_key(|&atom| match marks[atom] {
            Some(key) => (0u8, key, ranks[atom]),
            None => (1u8, (0, 0, 0), ranks[atom]),
        });
        for (&atom, &slot) in ordered.iter().zip(&slots) {
            informed[atom] = slot;
        }
    }
    informed
}

/// The atom's isotopic mark: `Some((mass delta, T, D))` for a heavy
/// isotope or a fixed isotopic-H carrier, `None` for a plain atom.
fn mark_of(structure: &Structure, index: usize) -> Option<(i32, u32, u32)> {
    let atom = &structure.atoms[index];
    let delta = if atom.isotope != 0 && atom.element != "H" {
        atom.isotope as i32 - standard_mass(&atom.element)?
    } else {
        0
    };
    let (d, t) = structure.isotopic_h.get(&index).copied().unwrap_or((0, 0));
    if delta == 0 && d == 0 && t == 0 {
        return None;
    }
    Some((delta, t, d))
}

/// The recorded mass minus the element's standard mass, signed —
/// `None` when the atom is the standard isotope or an H (H isotopes
/// ride the carrier).
fn heavy_entry(isotope: u16, element: &str) -> Option<String> {
    if isotope == 0 || element == "H" {
        return None;
    }
    let delta = isotope as i32 - standard_mass(element)?;
    Some(format!("{delta:+}"))
}

fn h_entry(t: u32, d: u32) -> String {
    let mut entry = String::new();
    if t > 0 {
        entry.push('T');
        if t > 1 {
            entry.push_str(&t.to_string());
        }
    }
    if d > 0 {
        entry.push('D');
        if d > 1 {
            entry.push_str(&d.to_string());
        }
    }
    entry
}

fn mobile_carrier(structure: &Structure, carrier: usize) -> bool {
    let atom = &structure.atoms[carrier];
    matches!(atom.element.as_str(), "O" | "S" | "Se" | "Te" | "N")
        && !structure.bonds.iter().any(|&(from, to, _)| {
            let other = if from == carrier {
                to
            } else if to == carrier {
                from
            } else {
                return false;
            };
            structure.atoms[other].element == "C"
        })
}

/// The element's standard (lowest-mass stable) isotope.
fn standard_mass(element: &str) -> Option<i32> {
    Some(match element {
        "H" => 1,
        "B" => 11,
        "C" => 12,
        "N" => 14,
        "O" => 16,
        "F" => 19,
        "P" => 31,
        "S" => 32,
        "Cl" => 35,
        "Br" => 79,
        "I" => 127,
        "Se" => 80,
        "Te" => 128,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;

    fn atom(element: &str, h_count: u32, isotope: u16) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count,
            h_pinned: false,
            charge: 0,
            isotope,
        }
    }

    fn structure(atoms: Vec<AtomRecord>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure {
            atoms,
            bonds,
            ..Default::default()
        }
    }

    #[test]
    fn heavy_deltas_are_signed_against_the_standard_mass() {
        // Oracle battery: [18OH2] -> i1+2; [125IH] -> i1-2;
        // [13CH4] -> i1+1 (fixtures breadth2-1116/1073/c13-methane).
        let s = structure(vec![atom("O", 2, 18)], vec![]);
        assert_eq!(IsotopeLayer.emit(&s).unwrap().as_deref(), Some("i1+2"));

        let s = structure(vec![atom("I", 1, 125)], vec![]);
        assert_eq!(IsotopeLayer.emit(&s).unwrap().as_deref(), Some("i1-2"));

        let s = structure(vec![atom("C", 4, 13)], vec![]);
        assert_eq!(IsotopeLayer.emit(&s).unwrap().as_deref(), Some("i1+1"));
    }

    #[test]
    fn folded_deuterium_reads_per_carrier_and_mobile() {
        // CD3OD (oracle battery): fixed D3 on the carbon rank plus one
        // on the oxygen rank — i1D3,2D. Water's D2 rides the mobile
        // aggregate — i/hD2 (fixture identity/d2o).
        let mut s = structure(vec![atom("C", 0, 0), atom("O", 1, 0)], vec![(0, 1, 1)]);
        s.isotopic_h.insert(0, (3, 0));
        s.isotopic_h.insert(1, (1, 0));
        assert_eq!(IsotopeLayer.emit(&s).unwrap().as_deref(), Some("i1D3,2D"));

        let mut s = structure(vec![atom("O", 2, 0)], vec![]);
        s.isotopic_h.insert(0, (2, 0));
        assert_eq!(IsotopeLayer.emit(&s).unwrap().as_deref(), Some("i/hD2"));
    }

    #[test]
    fn entries_sort_by_canonical_rank_not_atom_order() {
        // fixture identity/breadth2-1984: i1D,2D,3D,4D — the recorded
        // carriers sit at ranks 1..4 in rank order, not input order.
        let mut s = structure(
            vec![
                atom("C", 1, 0),
                atom("C", 1, 0),
                atom("C", 1, 0),
                atom("C", 1, 0),
            ],
            vec![(0, 1, 2), (1, 2, 1), (2, 3, 1), (3, 0, 1)],
        );
        s.isotopic_h.insert(2, (1, 0));
        s.isotopic_h.insert(0, (1, 0));
        s.isotopic_h.insert(3, (1, 0));
        s.isotopic_h.insert(1, (1, 0));
        assert_eq!(
            IsotopeLayer.emit(&s).unwrap().as_deref(),
            Some("i1D,2D,3D,4D")
        );
    }

    #[test]
    fn tritium_precedes_deuterium_and_counts_append() {
        // Oracle battery CD2HOT: i1TD2,2D; T2O: i/hT2.
        let mut s = structure(vec![atom("C", 0, 0), atom("O", 1, 0)], vec![(0, 1, 1)]);
        s.isotopic_h.insert(0, (2, 1));
        s.isotopic_h.insert(1, (1, 0));
        assert_eq!(IsotopeLayer.emit(&s).unwrap().as_deref(), Some("i1TD2,2D"));

        let mut s = structure(vec![atom("O", 2, 0)], vec![]);
        s.isotopic_h.insert(0, (0, 2));
        assert_eq!(IsotopeLayer.emit(&s).unwrap().as_deref(), Some("i/hT2"));
    }

    #[test]
    fn the_hydrogen_fold_makes_explicit_h_invisible_to_the_skeleton() {
        // D2O parses to ONE atom: the D's fold onto the oxygen —
        // /h h1H2 (isotope-blind) and the /c layer goes silent.
        let molfile = "d2o\n  gem\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\nM  ISO  2   2   2   3   2\nM  END\n";
        let mut s = crate::molfile::parse(molfile).unwrap();
        s.fold_explicit_hydrogens();
        assert_eq!(s.atoms.len(), 1);
        assert_eq!(s.atoms[0].h_count, 2);
        assert_eq!(s.isotopic_h.get(&0), Some(&(2, 0)));
        assert_eq!(s.bonds.len(), 0);
    }
}
