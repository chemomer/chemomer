//! Double-bond stereo (`/b`) — M4 slice 2, decoded from oracle
//! probes. A stereogenic double bond: both endpoints carry two
//! DISTINCT ligands (heavy neighbor ranks plus at most one implicit
//! H — a terminal =CH2 has two identical H and is out), and the bond
//! is not inside a ring smaller than 8. The reference ligands are the
//! highest canonical ranks at each end; the sign is + when they lie
//! on OPPOSITE sides of the drawn double bond, − when on the same
//! side (the probes' E/Z mapping). The segment lists each pair with
//! the endpoint ranks descending; /b comes between /p and /t.
//!
//! Non-defined members: a crossed/either code (3 or 4) on the double
//! bond marks UNKNOWN parity — labelled `u` under the SLUUD knob and
//! `?` otherwise; a stereogerm whose drawn sides are degenerate is
//! UNDEFINED (`?`). Both ride along a group that holds at least one
//! well-defined bond, and an all-non-defined group drops unless the
//! SUU knob keeps it (the probes: `-SUU` on the collinear/crossed
//! butene records `b4-3?`, `-SUU -SLUUD` records `b4-3u`).

use super::stereo::Emission;
use crate::graph::Structure;

/// One ligand: its canonical rank (0 for the implicit hydrogen) and,
/// for heavy ligands, its drawn position.
type Ligand = (u32, Option<(f64, f64)>);

/// The derived `/b` segment; empty when no double bond carries stereo.
pub fn segments(
    structure: &Structure,
    counts: &[u32],
    ranks: &[u32],
    emission: Emission,
) -> Vec<(char, String)> {
    let mut defined: Vec<(u32, u32, char)> = Vec::new();
    let mut unknown: Vec<(u32, u32, char)> = Vec::new();
    let mut undefined: Vec<(u32, u32)> = Vec::new();
    for (index, &(from, to, kind)) in structure.bonds.iter().enumerate() {
        if kind != 2 {
            continue;
        }
        if matches!(structure.stereo_bonds.get(&index).copied(), Some(3 | 4)) {
            // The crossed/either code: UNKNOWN parity, never geometry.
            unknown.push((
                ranks[from].max(ranks[to]),
                ranks[from].min(ranks[to]),
                if emission.label_unknown { 'u' } else { '?' },
            ));
            continue;
        }
        let Some(bonds_geometry) = stereogerm(structure, counts, ranks, from, to) else {
            continue;
        };
        match bond_sign(ranks, from, to, bonds_geometry) {
            Some((sign, hi, lo)) => defined.push((hi, lo, sign)),
            None => undefined.push((ranks[from].max(ranks[to]), ranks[from].min(ranks[to]))),
        }
    }
    if defined.is_empty() && !emission.keep_undefined {
        // No well-defined bond parity: the group drops unless SUU
        // keeps the unknown/undefined members.
        return Vec::new();
    }
    let mut bonds: Vec<(u32, u32, char)> = defined;
    bonds.extend(unknown);
    bonds.extend(undefined.into_iter().map(|(hi, lo)| (hi, lo, '?')));
    if bonds.is_empty() {
        return Vec::new();
    }
    bonds.sort_unstable();
    let text: String = bonds
        .iter()
        .map(|(hi, lo, sign)| format!("{hi}-{lo}{sign}"))
        .collect::<Vec<_>>()
        .join(",");
    vec![('b', format!("b{text}"))]
}

/// The stereogerm verdict for one double bond: distinct ligands at
/// both ends, not small-ring, and the reference ligands' drawn
/// positions present — with the coordinates the side product reads
/// (bond direction, left and right reference ligands).
type Geometry = ((f64, f64), (f64, f64), (f64, f64), (f64, f64));

fn stereogerm(
    structure: &Structure,
    counts: &[u32],
    ranks: &[u32],
    a: usize,
    b: usize,
) -> Option<Geometry> {
    if in_small_ring(structure, a, b) {
        return None;
    }
    let ligands_of = |atom: usize, other: usize| -> Option<Vec<Ligand>> {
        let heavy: Vec<(u32, Option<(f64, f64)>)> = structure
            .bonds
            .iter()
            .filter_map(|&(from, to, _)| {
                let neighbor = if from == atom && to != other {
                    Some(to)
                } else if to == atom && from != other {
                    Some(from)
                } else {
                    None
                }?;
                Some((ranks[neighbor], structure.positions.get(neighbor).copied()))
            })
            .collect();
        let mut ligands = heavy;
        if counts[atom] > 0 {
            if counts[atom] > 1 {
                return None; // =CH2 or =NH2: two identical hydrogens
            }
            ligands.push((0, None)); // the implicit H ranks last
        }
        Some(ligands)
    };
    let left = ligands_of(a, b)?;
    let right = ligands_of(b, a)?;
    for ligands in [&left, &right] {
        if ligands.len() != 2 {
            return None;
        }
        let mut ranked: Vec<u32> = ligands.iter().map(|(rank, _)| *rank).collect();
        ranked.sort_unstable();
        ranked.dedup();
        if ranked.len() != 2 {
            return None;
        }
    }
    // Reference ligands: the highest ranks (the implicit H, rank 0,
    // never wins — every stereogenic ligand set here has a heavy
    // member).
    let reference = |ligands: &[Ligand]| -> Ligand {
        ligands
            .iter()
            .max_by_key(|(rank, _)| *rank)
            .copied()
            .unwrap_or((0, None))
    };
    let (_, left_at) = reference(&left);
    let (_, right_at) = reference(&right);
    let (ax, ay) = structure.positions.get(a).copied()?;
    let (bx, by) = structure.positions.get(b).copied()?;
    let (lx, ly) = left_at?;
    let (rx, ry) = right_at?;
    Some(((ax, ay), (bx, by), (lx, ly), (rx, ry)))
}

/// The sign for one stereogerm double bond, with its endpoint ranks
/// descending; `None` when the drawn sides carry no stereo
/// information (the UNDEFINED parity the SUU knob prints as `?`).
fn bond_sign(ranks: &[u32], a: usize, b: usize, geometry: Geometry) -> Option<(char, u32, u32)> {
    let ((ax, ay), (bx, by), (lx, ly), (rx, ry)) = geometry;
    // A degenerate layout cannot carry stereo — the oracle ignores
    // it exactly the same way (coordinate-free fixtures record no
    // /b at all).
    let degenerate = |dx: f64, dy: f64| dx * dx + dy * dy < 1e-9;
    if degenerate(bx - ax, by - ay) || degenerate(lx - ax, ly - ay) || degenerate(rx - ax, ry - ay)
    {
        return None;
    }
    // Sides measured against the same double-bond direction; a zero
    // side (collinear substituent) carries no stereo information —
    // the oracle records nothing for those layouts.
    let side = |px: f64, py: f64| (bx - ax) * (py - ay) - (by - ay) * (px - ax);
    let product = side(lx, ly) * side(rx, ry);
    if product.abs() < 1e-9 {
        return None;
    }
    Some((
        if product < 0.0 { '+' } else { '-' },
        ranks[a].max(ranks[b]),
        ranks[a].min(ranks[b]),
    ))
}

/// The shortest cycle through the bond; `true` when it is smaller
/// than 8 — the reference leaves small-ring double bonds unstereod.
fn in_small_ring(structure: &Structure, a: usize, b: usize) -> bool {
    // BFS over the graph without b: the cycle is min over b's
    // neighbors of (distance from a) + 2.
    use std::collections::VecDeque;
    let adjacency = |atom: usize| -> Vec<usize> {
        structure
            .bonds
            .iter()
            .filter_map(|&(from, to, _)| {
                if from == atom {
                    Some(to)
                } else if to == atom {
                    Some(from)
                } else {
                    None
                }
            })
            .filter(|&next| next != b)
            .collect()
    };
    let mut distance = vec![usize::MAX; structure.atoms.len()];
    let mut queue = VecDeque::new();
    distance[a] = 0;
    queue.push_back(a);
    while let Some(atom) = queue.pop_front() {
        for next in adjacency(atom) {
            if distance[next] == usize::MAX {
                distance[next] = distance[atom] + 1;
                queue.push_back(next);
            }
        }
    }
    structure
        .bonds
        .iter()
        .filter(|&&(from, to, _)| from == b || to == b)
        .map(|&(from, to, _)| if from == b { to } else { from })
        .filter(|&neighbor| neighbor != a) // the bond itself is not a cycle
        .map(|neighbor| distance[neighbor].saturating_add(2))
        .min()
        .map(|cycle| cycle < 8)
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;

    fn atom(element: &str, h_count: u32) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count,
            h_pinned: false,
            isotope: 0,
            charge: 0,
        }
    }

    fn butene(coords: Vec<(f64, f64)>) -> Structure {
        Structure {
            atoms: vec![atom("C", 3), atom("C", 1), atom("C", 1), atom("C", 3)],
            bonds: vec![(0, 1, 1), (1, 2, 2), (2, 3, 1)],
            stereo_bonds: Default::default(),
            positions: coords,
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            chiral_flag: false,
        }
    }

    #[test]
    fn trans_and_cis_2_butene_match_the_oracle_signs() {
        // Oracle: trans -> b4-3+, cis -> b4-3- over ranks [1,3,4,2].
        let ranks = [1u32, 3, 4, 2];
        let trans = butene(vec![
            (12.0, 85.0),
            (60.0, 50.0),
            (120.0, 50.0),
            (168.0, 15.0),
        ]);
        assert_eq!(
            segments(
                &trans,
                &[3, 1, 1, 3],
                &ranks,
                crate::canonical::stereo::Emission::default()
            ),
            vec![('b', "b4-3+".to_string())]
        );
        let cis = butene(vec![
            (12.0, 15.0),
            (60.0, 50.0),
            (120.0, 50.0),
            (168.0, 15.0),
        ]);
        assert_eq!(
            segments(
                &cis,
                &[3, 1, 1, 3],
                &ranks,
                crate::canonical::stereo::Emission::default()
            ),
            vec![('b', "b4-3-".to_string())]
        );
    }

    #[test]
    fn terminal_methylene_is_not_stereogenic() {
        // Propene: the =CH2 end carries two identical H.
        let structure = Structure {
            atoms: vec![atom("C", 3), atom("C", 1), atom("C", 2)],
            bonds: vec![(0, 1, 1), (1, 2, 2)],
            stereo_bonds: Default::default(),
            positions: vec![(12.0, 22.0), (72.0, 50.0), (132.0, 50.0)],
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            chiral_flag: false,
        };
        assert!(segments(
            &structure,
            &[3, 1, 2],
            &[1, 2, 3],
            crate::canonical::stereo::Emission::default()
        )
        .is_empty());
    }

    #[test]
    fn small_ring_double_bonds_are_not_stereogenic() {
        // Cyclohexene skeleton: the double bond sits in a 6-ring.
        let coords: Vec<(f64, f64)> = (0..6)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 6.0;
                (60.0 + 30.0 * angle.cos(), 50.0 + 30.0 * angle.sin())
            })
            .collect();
        let bonds: Vec<(usize, usize, u32)> = (0..6)
            .map(|i| (i, (i + 1) % 6, if i == 0 { 2 } else { 1 }))
            .collect();
        let structure = Structure {
            atoms: (0..6)
                .map(|i| atom("C", if i % 2 == 0 { 1 } else { 2 }))
                .collect(),
            bonds,
            stereo_bonds: Default::default(),
            positions: coords,
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            chiral_flag: false,
        };
        assert!(segments(
            &structure,
            &[1, 2, 1, 2, 1, 2],
            &[1, 2, 3, 4, 5, 6],
            crate::canonical::stereo::Emission::default()
        )
        .is_empty());
    }
}
