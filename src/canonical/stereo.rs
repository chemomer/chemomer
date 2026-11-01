//! Tetrahedral stereo (`/t`, `/m`, `/s`) — M4 slice 1, ported from
//! the reference (ichister.c `set_stereo_atom_parity`; dev-time
//! source read, never compiled or called). Input contract: a wedge
//! (1) or dash (6) bond-stereo code whose FIRST atom is the
//! stereocenter — the gem writer's shape; the atom-parity field is
//! ignored for 2D input, and degenerate coordinates void wedges. A
//! wiggly bond (code 4) marks its center UNKNOWN — the parity the
//! standard reads `?` and the SLUUD option reads `u`.
//!
//! Per center: the heavy ligands in neighbor order are lifted to 3D
//! (the wedge carries z = ± its drawn length), the implicit hydrogen
//! closes the tetrahedron at −(sum), and the parity is the sign of
//! the triple product of ligands 1, 2, H viewed from ligand 0,
//! inverted once for the appended H, composed with the
//! input→canonical ligand permutation and one global convention
//! constant. Emission: the first defined center's `/t` char is
//! normalized `-`, other defined centers show `p_i ⊕ p_first`
//! (`-` same, `+` differ), undefined candidate centers `?`; `/m` is
//! the first center's parity and `/s1` rides it; two same-orbit
//! centers with opposite parity are the meso case — bare `/t`, no
//! `/m` `/s`. A group with no defined parity at all drops unless the
//! SUU option keeps it (then the bare `/t` prints, no `/m` `/s`).

use crate::graph::Structure;

/// The emission knobs the generation options map onto: whether a
/// stereo group with no well-defined parity still prints (SUU) and
/// whether unknown (wiggly) stereo carries its own label `u` (SLUUD).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Emission {
    pub keep_undefined: bool,
    pub label_unknown: bool,
}

/// One center's derived state.
struct Center {
    atom: usize,
    /// Heavy-neighbor atom indices in input order (atom-index order).
    neighbors: Vec<usize>,
    parity: Option<u32>,
    /// A wiggly bond (code 4) points at this center: UNKNOWN parity,
    /// labelled `u` under SLUUD and `?` otherwise.
    unknown: bool,
}

impl Center {
    /// The heavy ligands in neighbor order — the reference pushes
    /// neighbors in molfile bond order (mol2atom.c), which the parity
    /// is defined over; atom-number order diverges from it exactly on
    /// centers whose first-listed bonds differ from their numbering.
    fn heavy_by_neighbor_order(&self) -> Vec<usize> {
        self.neighbors.clone()
    }

    /// The permutation sequence the parity composes with: the heavy
    /// ligands' input positions reordered rank-ascending (the implicit
    /// H is rank-last and adds no transposition).
    fn permutation_sequence(&self, ranks: &[u32]) -> Vec<usize> {
        let mut remaining: Vec<(usize, u32)> = self
            .neighbors
            .iter()
            .enumerate()
            .map(|(position, &atom)| (position, ranks[atom]))
            .collect();
        remaining.sort_by_key(|&(_, rank)| rank);
        remaining
            .into_iter()
            .map(|(position, _)| position)
            .collect()
    }
}

/// The composition constant between the reference's atom parity and
/// our canonical parity bit, one global convention flip fitted on the
/// alanine pair.
const PARITY_CONSTANT: u32 = 1;

/// The absolute parity convention the AUX stereo rows speak (the
/// reference's t_parity): witnessed against the instrumented oracle —
/// l-alanine's absolute pass reads parity 2, its inverted pass 1.
const AUX_PARITY_CONSTANT: u32 = 0;

/// The derived segments in `/t`, `/m`, `/s` order; empty when the
/// structure carries no wedge parity at all and the SUU knob does not
/// keep the undefined group.
pub fn segments(
    structure: &Structure,
    hydrogens: &[u32],
    ranks: &[u32],
    symmetry: &[u32],
    emission: Emission,
) -> Vec<(char, String)> {
    let wedged = wedged_centers(structure);
    let defined_exists = wedged.iter().any(|&(_, _, down)| down.is_some());
    if !defined_exists && !emission.keep_undefined {
        return Vec::new();
    }
    let mut centers = candidates(structure, hydrogens, ranks);
    for &(atom, bond_index, down) in &wedged {
        if let Some(center) = centers.iter_mut().find(|c| c.atom == atom) {
            let Some(down) = down else {
                // The wiggly bond: UNKNOWN parity, no geometry read.
                center.unknown = true;
                continue;
            };
            // Heavy ligands in the reference's neighbor order, the
            // implicit H appended last — the list the atom parity is
            // defined over.
            let ligands = center.heavy_by_neighbor_order();
            let Some(bit) = atom_parity(structure, &ligands, bond_index, down) else {
                continue;
            };
            let perm = permutation_parity(center.permutation_sequence(ranks));
            center.parity = Some(bit ^ perm ^ PARITY_CONSTANT);
        }
    }
    let Some(first) = centers.iter().find(|c| c.parity.is_some()) else {
        // No defined parity anywhere — wedges all degenerate, wiggly
        // or absent: the group drops unless the SUU knob keeps it
        // (then the bare `/t` prints over every candidate — unknown
        // centres labelled, undefined ones `?` — with no `/m` `/s`:
        // all-non-well-def rows leave nCompInv2Abs at 0, SP3_ONLY).
        // No perceivable center at all is SP3_NONE: nothing prints.
        if centers.is_empty() || !emission.keep_undefined {
            return Vec::new();
        }
        return vec![('t', render_t(&centers, ranks, None, emission.label_unknown))];
    };
    let first_parity = first.parity.expect("checked");

    // Meso: same-orbit defined centers with opposite parity — the
    // relative case takes neither /m nor /s.
    let mut meso = false;
    for (i, a) in centers.iter().enumerate() {
        for b in centers.iter().skip(i + 1) {
            if let (Some(pa), Some(pb)) = (a.parity, b.parity) {
                if symmetry[a.atom] == symmetry[b.atom] && pa != pb {
                    meso = true;
                }
            }
        }
    }

    let mut out = vec![(
        't',
        render_t(&centers, ranks, Some(first_parity), emission.label_unknown),
    )];
    if !meso {
        out.push(('m', format!("m{first_parity}")));
        out.push(('s', "s1".to_string()));
    }
    out
}

/// The `/t` text: per center its rank and mark — defined centers
/// relative to `first` (None renders every defined center `-`, the
/// all-undefined bare group), unknown centers `u`/`?` by the SLUUD
/// knob, undefined candidates `?`.
fn render_t(centers: &[Center], ranks: &[u32], first: Option<u32>, label_unknown: bool) -> String {
    let mut t = String::from("t");
    for (index, center) in centers.iter().enumerate() {
        let mark = match center.parity {
            None if center.unknown => {
                if label_unknown {
                    'u'
                } else {
                    '?'
                }
            }
            None => '?',
            Some(p) => match first {
                Some(first_parity) if index > 0 && p != first_parity => '+',
                _ => '-',
            },
        };
        if index > 0 {
            t.push(',');
        }
        t.push_str(&format!("{}{}", ranks[center.atom], mark));
    }
    t
}

/// Permutation parity (1 = odd) of the sequence that orders the heavy
/// ligands rank-ascending; the implicit H is rank-last and adds no
/// transposition.
fn permutation_parity(sequence: Vec<usize>) -> u32 {
    u32::try_from(inversions(&sequence) % 2).expect("parity fits u32")
}

fn inversions(sequence: &[usize]) -> usize {
    let mut count = 0;
    for i in 0..sequence.len() {
        for j in (i + 1)..sequence.len() {
            if sequence[i] > sequence[j] {
                count += 1;
            }
        }
    }
    count
}
/// The reference's atom parity (ichister.c `set_stereo_atom_parity`):
/// the heavy ligands in neighbor[] order are lifted to 3D — the wedge
/// bond carries z = ± its drawn length, in-plane neighbors z = 0 — and
/// each is normalized; the implicit hydrogen closes the tetrahedron at
/// −(sum of the three); with the origin moved to ligand 0, the parity
/// is the sign of the triple product of ligands 1, 2, and the H
/// (negative = odd), inverted once because the H was appended last.
/// `None` for degenerate geometry.
fn atom_parity(
    structure: &Structure,
    ligands: &[usize],
    wedge_bond: usize,
    down: u32,
) -> Option<u32> {
    let position = |atom: usize| structure.positions.get(atom).copied();
    let (cx, cy) = position(structure.bonds[wedge_bond].0)?;
    let sign = if down == 1 { -1.0 } else { 1.0 };
    let mut coord = [[0.0f64; 3]; 3];
    for (slot, &atom) in ligands.iter().take(3).enumerate() {
        let (x, y) = position(atom)?;
        let (dx, dy) = (x - cx, y - cy);
        let len = (dx * dx + dy * dy).sqrt();
        if len < 1e-9 {
            return None;
        }
        let z = if structure.bonds[wedge_bond].1 == atom {
            sign * len
        } else {
            0.0
        };
        coord[slot] = [dx / len, dy / len, z / len];
    }
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    // The implicit hydrogen closes the tetrahedron.
    let hydrogen = [
        -(coord[0][0] + coord[1][0] + coord[2][0]),
        -(coord[0][1] + coord[1][1] + coord[2][1]),
        -(coord[0][2] + coord[1][2] + coord[2][2]),
    ];
    // Origin at ligand 0; orientation of ligands 1, 2, H viewed from it.
    let v1 = sub(coord[1], coord[0]);
    let v2 = sub(coord[2], coord[0]);
    let v3 = sub(hydrogen, coord[0]);
    let triple = dot(cross(v1, v2), v3);
    if triple.abs() < 1e-12 {
        return None;
    }
    // Negative triple = odd; then the appended-H inversion (3 − parity).
    Some(u32::from(triple < 0.0) ^ 1)
}
/// Atoms the wedge codes point FROM (the code's first atom), with the
/// bond index and the direction (`Some(1)` = dash, `Some(0)` = up,
/// `None` = the wiggly code 4 — the UNKNOWN parity mark).
fn wedged_centers(structure: &Structure) -> Vec<(usize, usize, Option<u32>)> {
    let mut out = Vec::new();
    for (&bond_index, &code) in &structure.stereo_bonds {
        let Some(&(from, _, _)) = structure.bonds.get(bond_index) else {
            continue;
        };
        let down = match code {
            6 => Some(1),
            1 => Some(0),
            4 => None,
            _ => continue,
        };
        out.push((from, bond_index, down));
    }
    out.sort();
    out
}

/// The wedge-defined centers with their geometry parity bit and heavy
/// ligands (neighbor order) — the per-terminal stereo-row inputs the
/// aux numbering's informed search compares. The geometry bit is
/// terminal-independent; only the ligand permutation composes with
/// the canonical ranks. Wiggly-marked centers carry no geometry bit
/// and stay out (the rows speak defined parity only).
pub fn center_inputs(
    structure: &Structure,
    hydrogens: &[u32],
    ranks: &[u32],
) -> Vec<(usize, Vec<usize>, u32)> {
    let wedged = wedged_centers(structure);
    if wedged.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::new();
    for &(atom, bond_index, down) in &wedged {
        let Some(down) = down else {
            continue;
        };
        let Some(center) = candidates(structure, hydrogens, ranks)
            .into_iter()
            .find(|c| c.atom == atom)
        else {
            continue;
        };
        let ligands = center.heavy_by_neighbor_order();
        let Some(bit) = atom_parity(structure, &ligands, bond_index, down) else {
            continue;
        };
        out.push((atom, ligands, bit));
    }
    out.sort_by_key(|&(atom, _, _)| ranks[atom]);
    out
}

/// The (center rank, parity) pairs of [`terminal_rows`] — the compare
/// the reference's CopyLinearCTStereoToINChIStereo diff runs.
pub fn terminal_pairs(
    centers: &[(usize, Vec<usize>, u32)],
    ranks: &[u32],
    flip: bool,
) -> Vec<(u32, u32)> {
    let mut ordered: Vec<&(usize, Vec<usize>, u32)> = centers
        .iter()
        .filter(|(atom, _, _)| ranks[*atom] > 0)
        .collect();
    ordered.sort_by_key(|&&(atom, _, _)| ranks[atom]);
    let mut out = Vec::new();
    for &(atom, ref ligands, bit) in ordered {
        let mut sequence: Vec<(usize, u32)> = ligands
            .iter()
            .enumerate()
            .map(|(position, &ligand)| (position, ranks[ligand]))
            .collect();
        sequence.sort_by_key(|&(_, rank)| rank);
        let perm = permutation_parity(sequence.into_iter().map(|(position, _)| position).collect());
        let parity = 2 - ((bit ^ perm ^ AUX_PARITY_CONSTANT) ^ u32::from(flip));
        out.push((ranks[atom], parity));
    }
    out
}

/// The stereo-informed ranking verdict — the machinery the AuxInfo
/// speaks (the reference's nCanonOrdStereo + nCompInv2Abs): the
/// winning pass's ranks and symmetry, both passes' terminal pairs,
/// and the inversion sign (`-1` the inverted pass won and its
/// numbering was adopted, `0` the rows tie, `1` the absolute pass won
/// on rows alone). `None` when the structure carries no wedge-defined
/// center. Callers assemble the exact ranking inputs (tautomeric
/// structures strip the mobile H and drop the group rows) so the aux
/// and the identifier speak one verdict.
pub(crate) struct StereoVerdict {
    pub ranks: Vec<u32>,
    pub symmetry: Vec<u32>,
    pub main_pairs: Vec<(u32, u32)>,
    pub inv_pairs: Vec<(u32, u32)>,
    pub sign: i32,
}

/// Computes [`StereoVerdict`] over `counts` (the adjusted hydrogen
/// counts the centers perceive against) and `ranks` (the base
/// canonical ranks [`crate::layers::hydrogen::canonical_ranks_for`]
/// supplies). `strip_mobile_h` mirrors the aux's tautomeric branch:
/// the ranking reads the mobile-H-stripped counts and no group rows.
pub(crate) fn informed_verdict(
    structure: &Structure,
    counts: &[u32],
    ranks: &[u32],
    strip_mobile_h: bool,
) -> Option<StereoVerdict> {
    let centers = center_inputs(structure, counts, ranks);
    if centers.is_empty() {
        return None;
    }
    let (adjacency, invariants) = crate::canonical::graph_of(structure);
    let (hydrogens, groups): (Vec<u32>, Vec<(u32, Vec<usize>)>) = if strip_mobile_h {
        let mask = crate::canonical::tautomer::endpoint_mask(structure, counts);
        (
            crate::canonical::tautomer::without_mobile_h(counts, &mask),
            Vec::new(),
        )
    } else {
        (
            counts.to_vec(),
            crate::canonical::group_membership(structure, counts),
        )
    };
    let (ranks, symmetry, (main_pairs, inv_pairs, sign)) =
        crate::canonical::canonical_ranks_and_symmetry_stereo(
            &adjacency,
            &invariants,
            &hydrogens,
            &groups,
            &centers,
        );
    Some(StereoVerdict {
        ranks,
        symmetry,
        main_pairs,
        inv_pairs,
        sign,
    })
}

/// One terminal's stereo rows over `ranks`: per center (ascending
/// canonical number) the center's rank, the parity number (1 = odd,
/// 2 = even — the geometry bit composed with the ligand permutation
/// and the convention constant, flipped once for the inverted pass),
/// then the ligand ranks ascending. The rows extend the terminal
/// compare the way the reference's LinearCTStereoCarb section does.
pub fn terminal_rows(centers: &[(usize, Vec<usize>, u32)], ranks: &[u32], flip: bool) -> Vec<u32> {
    let mut out = Vec::new();
    let mut ordered: Vec<&(usize, Vec<usize>, u32)> = centers
        .iter()
        .filter(|(atom, _, _)| ranks[*atom] > 0)
        .collect();
    ordered.sort_by_key(|&&(atom, _, _)| ranks[atom]);
    for &(atom, ref ligands, bit) in ordered {
        let mut sequence: Vec<(usize, u32)> = ligands
            .iter()
            .enumerate()
            .map(|(position, &ligand)| (position, ranks[ligand]))
            .collect();
        sequence.sort_by_key(|&(_, rank)| rank);
        let perm = permutation_parity(sequence.into_iter().map(|(position, _)| position).collect());
        let parity = 2 - ((bit ^ perm ^ AUX_PARITY_CONSTANT) ^ u32::from(flip));
        out.push(ranks[atom]);
        out.push(parity);
        let mut ligand_ranks: Vec<u32> = ligands.iter().map(|&l| ranks[l]).collect();
        ligand_ranks.sort_unstable();
        out.extend(ligand_ranks);
    }
    out
}

/// Every candidate center: three heavy neighbors plus one implicit H,
/// or four heavy neighbors, with all ligand ranks distinct.
fn candidates(structure: &Structure, hydrogens: &[u32], ranks: &[u32]) -> Vec<Center> {
    let mut out = Vec::new();
    for (atom, &hydrogen_count) in hydrogens.iter().enumerate() {
        let neighbors: Vec<usize> = structure
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
            .collect();
        // Three heavy neighbors plus one implicit hydrogen, or four
        // heavy neighbors.
        let tetrahedral = (hydrogen_count == 1 && neighbors.len() == 3)
            || (hydrogen_count == 0 && neighbors.len() == 4);
        if !tetrahedral {
            continue;
        }
        let mut ligands: Vec<u32> = neighbors.iter().map(|&a| ranks[a]).collect();
        ligands.sort_unstable();
        ligands.dedup();
        if ligands.len() != neighbors.len() {
            continue;
        }
        out.push(Center {
            atom,
            neighbors,
            parity: None,
            unknown: false,
        });
    }
    out.sort_by_key(|center| ranks[center.atom]);
    out
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

    /// 2-butanol: atom1 is the center; the wedge rides the bond to the
    /// OH (atom2) unless the layout moves it. Recorded oracle answers:
    /// the canonical numbering is c1-3-4(2)5 (ranks [5,4,2,3,1]).
    fn butanol(wedge: u8, positions: Vec<(f64, f64)>, wedge_bond: usize) -> Structure {
        Structure {
            atoms: vec![
                atom("C", 3),
                atom("C", 1),
                atom("O", 2),
                atom("C", 2),
                atom("C", 3),
            ],
            bonds: vec![(0, 1, 1), (1, 2, 1), (1, 3, 1), (3, 4, 1)],
            stereo_bonds: if wedge == 0 {
                Default::default()
            } else {
                [(wedge_bond, wedge)].into_iter().collect()
            },
            positions,
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            chiral_flag: false,
        }
    }

    fn check(name: &str, structure: &Structure, expected: &[(&str, &str)]) {
        let hydrogens = crate::graph::Structure::hydrogens(structure);
        let (adjacency, invariants) = crate::canonical::graph_of(structure);
        let symmetry = crate::canonical::base_symmetry(&adjacency, &invariants);
        let ranks = crate::canonical::canonical_ranks(&adjacency, &invariants, &hydrogens, &[]);
        let got = segments(
            structure,
            &hydrogens,
            &ranks,
            &symmetry,
            Emission::default(),
        );
        let want: Vec<(char, String)> = expected
            .iter()
            .map(|(k, v)| (k.chars().next().unwrap(), v.to_string()))
            .collect();
        assert_eq!(got, want, "{name}");
    }

    #[test]
    fn the_probe_grid_reproduces_the_oracle_recordings() {
        // The identity-probe layout: wedge centre->OH up gives /m1,
        // down /m0 (oracle-recorded).
        let layout = vec![
            (12.0, 22.0),
            (72.0, 67.0),
            (132.0, 22.0),
            (72.0, 12.0),
            (132.0, 62.0),
        ];
        check(
            "original up",
            &butanol(1, layout.clone(), 1),
            &[("t", "t4-"), ("m", "m1"), ("s", "s1")],
        );
        check(
            "original down",
            &butanol(6, layout.clone(), 1),
            &[("t", "t4-"), ("m", "m0"), ("s", "s1")],
        );
    }

    #[test]
    fn mirrored_layouts_flip_the_parity_like_the_oracle() {
        // Probe grid A (CH3 west, CH2 north) and its mirror A-m
        // (CH2 south): the same wedge code flips /m between them.
        let a = vec![
            (27.0, 67.0),
            (72.0, 67.0),
            (117.0, 67.0),
            (72.0, 112.0),
            (87.0, 129.0),
        ];
        let a_mirrored = vec![
            (27.0, 67.0),
            (72.0, 67.0),
            (117.0, 67.0),
            (72.0, 22.0),
            (87.0, 5.0),
        ];
        check(
            "A up",
            &butanol(1, a.clone(), 1),
            &[("t", "t4-"), ("m", "m0"), ("s", "s1")],
        );
        check(
            "A down",
            &butanol(6, a.clone(), 1),
            &[("t", "t4-"), ("m", "m1"), ("s", "s1")],
        );
        check(
            "A-m up",
            &butanol(1, a_mirrored.clone(), 1),
            &[("t", "t4-"), ("m", "m1"), ("s", "s1")],
        );
        check(
            "A-m down",
            &butanol(6, a_mirrored, 1),
            &[("t", "t4-"), ("m", "m0"), ("s", "s1")],
        );
    }

    #[test]
    fn no_wedge_carries_no_stereo() {
        let layout = vec![
            (12.0, 22.0),
            (72.0, 67.0),
            (132.0, 22.0),
            (72.0, 12.0),
            (132.0, 62.0),
        ];
        // Default emission drops the all-undefined group...
        assert!(segments(
            &butanol(0, layout.clone(), 1),
            &[3, 1, 2, 2, 3],
            &[5, 4, 2, 3, 1],
            &[1, 2, 3, 4, 5],
            Emission::default()
        )
        .is_empty());
        // ...while the SUU knob keeps it as the bare /t of undefined
        // candidates (the oracle's -SUU recording on the flat probe:
        // t4?).
        assert_eq!(
            segments(
                &butanol(0, layout, 1),
                &[3, 1, 2, 2, 3],
                &[5, 4, 2, 3, 1],
                &[1, 2, 3, 4, 5],
                Emission {
                    keep_undefined: true,
                    label_unknown: false,
                }
            ),
            vec![('t', "t4?".to_string())]
        );
    }
}
