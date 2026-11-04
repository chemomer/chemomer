//! The DoNotAddH aux numbering's stereo pass — the reference's
//! `map_stereo_bonds4` walk over undefined-geometry double-bond marks
//! (ichimap4.c; dev-time source read, never compiled or called).
//!
//! The aux `/N:` reads `nCanonOrdStereo` whenever the stereo
//! canonicalization ran (ichimak2.c takes it in both the with-stereo
//! and the no-stereo branch), and the pass runs on every structure —
//! center-free included — because `bCanonStereo` only asks for the CT
//! buffers to exist (ichicano.c:1463). Under DoNotAddH the zero
//! coordinates leave every eligible double bond at AB_PARITY_UNDF, so
//! the pass is a pure numbering walk: both mapping sides seed from the
//! base winner's equivalence classes (the `FixCanonEquivalenceInfo`
//! rebuild of nRank/nAtomNumber), the marked bonds enumerate in
//! canonical order (`Next_SB_At_CanonRanks2` — the partner rank keeps
//! the SMALLEST adjacent candidate, not the marked bond's own), each
//! entry maps onto the to-side's candidate bonds with tie-break pins
//! (`map_an_atom2`: both sides' atoms pin to the cell's bottom rank,
//! the `DifferentiateRanks2` refinement cascading after each pin), and
//! the finish is `BreakAllTies` on both sides plus the position-wise
//! copy of the base ranks onto the to-side order — the numbering the
//! aux speaks.
//!
//! Witnessed mechanism (breadth-27, instrumented oracle): the ene bond
//! enumerates as (cr1, cr2) = (3, 1), so the from-side maps its rank-1
//! atom onto the to-side's rank-2 bond end; the pin flips the {1,3}
//! class, the refinement cascade reorders the rest, and the aux reads
//! the flipped order against the base pass's own.
//!
//! Scope guards: the port covers the UNDF-only shape (all coordinates
//! zero — the whole corpus); a drawn coordinate keeps the plain
//! corpus numbering. The CurTree/nSymmStereo candidate skip is
//! omitted: its unions fire only on repeat complete mappings, which
//! cannot overwrite the numbering the FIRST complete mapping fixes,
//! and the walk stops there.

use super::partition::Partition;
use crate::graph::Structure;

/// The strict stereo-bond candidacy gate (the reference's
/// `set_stereo_bonds_parity` counting core): every isolated double
/// bond ([`super::stereo::isolated_double_bonds`]) both of whose ends
/// carry valence + drawn hydrogen in the stereo-bond window [2, 3]
/// (MIN/MAX_NUM_STEREO_BOND_NEIGH), an eligible element
/// (`bCanAtomHaveAStereoBond`: C/Si/Ge neutral, N neutral or +1), no
/// middle-cumulene shape at the registering end, exactly one double
/// bond at the registering end, a legal opposite-end neighborhood
/// (single/double bonds only, one double-bond type, two doubles only
/// for an allene-middle end), and no ring smaller than
/// MIN_SB_RING_SIZE = 8. Zero coordinates leave every surviving mark
/// at AB_PARITY_UNDF. Returns the marks as unordered atom pairs.
pub(crate) fn undf_bond_marks(structure: &Structure, counts: &[u32]) -> Vec<(usize, usize)> {
    // The UNDF-only shape: the reference computes the parity from the
    // drawn geometry, and zero coordinates leave every mark undefined.
    // Any drawn coordinate means the geometry path — out of scope.
    if structure
        .positions
        .iter()
        .any(|&(x, y)| x != 0.0 || y != 0.0)
    {
        return Vec::new();
    }
    let isolated = super::stereo::isolated_double_bonds(structure);
    let neighbors = adjacency_of(structure);
    let bond_kind = |atom: usize, other: usize| -> u32 {
        structure
            .bonds
            .iter()
            .find_map(|&(from, to, kind)| {
                if (from == atom && to == other) || (to == atom && from == other) {
                    Some(kind)
                } else {
                    None
                }
            })
            .unwrap_or(0)
    };
    let eligible_element = |atom: usize| -> bool {
        let element = structure.atoms[atom].element.as_str();
        let charge = structure.atoms[atom].charge;
        matches!(
            (element, charge),
            ("C", 0) | ("Si", 0) | ("Ge", 0) | ("N", 0) | ("N", 1)
        )
    };
    let in_window = |atom: usize| -> bool {
        let valence = neighbors[atom].len() as u32 + counts[atom];
        (2..=3).contains(&valence)
    };
    // The middle-cumulene shape: two double bonds, no hydrogen — the
    // registering end skips (`set_stereo_bonds_parity`'s early out).
    let middle_cumulene = |atom: usize| -> bool {
        counts[atom] == 0
            && neighbors[atom].len() == 2
            && neighbors[atom]
                .iter()
                .all(|&other| bond_kind(atom, other) == 2)
    };
    // The bond indices for the aromatic exclusion.
    let bond_index = |atom: usize, other: usize| -> Option<usize> {
        structure.bonds.iter().position(|&(from, to, _)| {
            (from == atom && to == other) || (to == atom && from == other)
        })
    };
    let mut out = Vec::new();
    for high in 0..structure.atoms.len() {
        // The registration end: the reference collects a bond only at
        // its higher-numbered end (`at_1 > at_2`).
        for &low in &neighbors[high] {
            if low >= high {
                continue;
            }
            let Some(bond) = bond_index(high, low) else {
                continue;
            };
            if structure.bonds[bond].2 != 2 || !isolated.contains(&bond) {
                continue;
            }
            // The C=C slice: the reference's normalized state types the
            // C=N/N=N-adjacent bonds of amidine/azine systems with the
            // mark_alt alternation marks (BOND_ALT_13 et al. — get
            // _allowed_stereo_bond_type maps them to zero, the ONE_BAD
            // wrong-type), and nothing in the local structure predicts
            // where the flow puts them (the instrumented amidine
            // trio: the same N=C(=S) shape typed bt1 on one fixture
            // and bt0 on another). The carbon-carbon marks — where the
            // corpus witnesses the walk — stay; the C=N/N=N marks are
            // the alternation state's debt, dated in the ledger.
            if structure.atoms[high].element != "C" || structure.atoms[low].element != "C" {
                continue;
            }
            // The registering end: element, window, shape, one double.
            if !eligible_element(high) || !in_window(high) || middle_cumulene(high) {
                continue;
            }
            let doubles_here = neighbors[high]
                .iter()
                .filter(|&&next| bond_kind(high, next) == 2)
                .count();
            if doubles_here != 1 {
                continue;
            }
            // The opposite end: element, window, legal neighborhood.
            if !eligible_element(low) || !in_window(low) {
                continue;
            }
            let mut doubles_opposite = 0;
            let mut illegal = false;
            for &next in &neighbors[low] {
                match bond_kind(low, next) {
                    2 => doubles_opposite += 1,
                    1 => {}
                    _ => illegal = true,
                }
            }
            if illegal || doubles_opposite == 0 || doubles_opposite > 2 {
                continue;
            }
            if doubles_opposite == 2 && neighbors[low].len() != 2 {
                continue;
            }
            // No small ring through the bond (min_sb_ring_size = 8).
            // No small ring through the bond (min_sb_ring_size = 8).
            if super::stereo::in_small_ring_stereo(structure, high, low) {
                continue;
            }
            out.push((high, low));
        }
    }
    out
}

/// The stereo-pass numbering over `marks` — the two-sided mapping walk
/// seeded from the base winner's classes. `base_ranks` is the plain
/// corpus canonicalization's discrete ranks; `symmetry` its orbit
/// classes. Returns the final ranks per atom, or `None` when the walk
/// finds no complete mapping (the plain numbering stands).
pub(crate) fn stereo_pass_ranks(
    adjacency: &[Vec<usize>],
    base_ranks: &[u32],
    symmetry: &[u32],
    marks: &[(usize, usize)],
) -> Option<Vec<u32>> {
    let n = adjacency.len();
    if marks.is_empty() {
        return None;
    }
    // Both directions registered (save_a_stereo_bond saves both ends).
    let mut bond_neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(a, b) in marks {
        if !bond_neighbors[a].contains(&b) {
            bond_neighbors[a].push(b);
            bond_neighbors[b].push(a);
        }
    }
    // The pre-filled CT: entries (hi rank, lo rank) of the marks in
    // the base canonical order (FillSingleStereoDescriptors: atoms
    // walked in canonical order, each bond at its larger-rank end,
    // an atom's bonds ordered by the partner's rank).
    let mut ct: Vec<(u32, u32)> = Vec::new();
    for position in 0..n {
        let atom = canonical_atom(base_ranks, position);
        let mut partners: Vec<u32> = bond_neighbors[atom]
            .iter()
            .map(|&b| base_ranks[b])
            .filter(|&rank| rank < base_ranks[atom])
            .collect();
        partners.sort_unstable();
        for rank in partners {
            ct.push((base_ranks[atom], rank));
        }
    }
    // The seeds: both stacks carry the class partition (the reference's
    // FixCanonEquivalenceInfo rebuild of nRank/nAtomNumber).
    let seed: Vec<(u32, u32)> = (0..n).map(|atom| (symmetry[atom], 0)).collect();
    let mut state = Walk {
        from: Partition::initial_pair(&seed),
        to: Partition::initial_pair(&seed),
        used: bond_neighbors.iter().map(Vec::len).collect(),
    };
    let mut numbering = None;
    if map_bonds(
        adjacency,
        base_ranks,
        &ct,
        &bond_neighbors,
        &mut state,
        0,
        &mut numbering,
    ) {
        numbering
    } else {
        None
    }
}

/// The atom holding canonical number `position + 1`.
fn canonical_atom(base_ranks: &[u32], position: usize) -> usize {
    (0..base_ranks.len())
        .find(|&atom| base_ranks[atom] == position as u32 + 1)
        .expect("base ranks are a permutation")
}

fn adjacency_of(structure: &Structure) -> Vec<Vec<usize>> {
    let n = structure.atoms.len();
    let mut adjacency = vec![Vec::new(); n];
    for &(from, to, _) in &structure.bonds {
        adjacency[from].push(to);
        adjacency[to].push(from);
    }
    adjacency
}

/// The walk state: the two mapping sides and the to-side bond
/// counters (`bAtomUsedForStereo` — the marks' per-atom counts).
struct Walk {
    from: Partition,
    to: Partition,
    used: Vec<usize>,
}

impl Clone for Walk {
    fn clone(&self) -> Walk {
        Walk {
            from: self.from.clone(),
            to: self.to.clone(),
            used: self.used.clone(),
        }
    }
}

/// The bond-level mapping (map_stereo_bonds4's body): enumerate the
/// next CT entry, walk the to-side candidates, pin-and-refine both
/// sides per mapped endpoint, recurse. `numbering` receives the FIRST
/// complete mapping's result (`bFirstCT`).
#[allow(clippy::too_many_arguments)]
fn map_bonds(
    adjacency: &[Vec<usize>],
    base_ranks: &[u32],
    ct: &[(u32, u32)],
    bond_neighbors: &[Vec<usize>],
    state: &mut Walk,
    entry: usize,
    numbering: &mut Option<Vec<u32>>,
) -> bool {
    let n = adjacency.len();
    if entry >= ct.len() {
        // All bonds mapped: BreakAllTies on both sides, then the
        // position-wise copy — the to-side order carrying the base
        // ranks of the from-side atoms at the same positions
        // (nCanonRankTo[n2] = nCanonRankFrom[n1]).
        let mut from = state.from.clone();
        let mut to = state.to.clone();
        break_all_ties(&mut from, adjacency);
        break_all_ties(&mut to, adjacency);
        let mut ranks = vec![0u32; n];
        for position in 0..n {
            ranks[to.at_number()[position]] = base_ranks[from.at_number()[position]];
        }
        *numbering = Some(ranks);
        return true;
    }
    // The enumeration (Next_SB_At_CanonRanks2): seeded from the
    // previous CT entry, cr1 walks upward and cr2 keeps the smallest
    // adjacent partner rank; canon_rank2_inp resets whenever cr1
    // advances.
    let mut canon_rank1_inp = if entry > 0 { ct[entry - 1].0 } else { 0 };
    let mut canon_rank2_inp = if entry > 0 { ct[entry - 1].1 } else { 0 };
    if canon_rank1_inp < 2 {
        canon_rank1_inp = 2;
        canon_rank2_inp = 0;
    }
    let mut found: Option<(u32, u32)> = None;
    let mut cr1 = canon_rank1_inp;
    while cr1 <= n as u32 {
        let from_atom = canonical_atom(base_ranks, cr1 as usize - 1);
        let r1 = state.from.rank_of(from_atom);
        let mut cr2_best = n as u32;
        // The to-side cell of rank r1, positions top-down.
        let mut pos = r1 as usize;
        while pos > 0 {
            pos -= 1;
            if state.to.rank_of(state.to.at_number()[pos]) != r1 {
                break;
            }
            let s1 = state.to.at_number()[pos];
            if state.used[s1] == 0 {
                continue;
            }
            for &s2 in &bond_neighbors[s1] {
                if state.used[s2] == 0 {
                    continue;
                }
                let r2 = state.to.rank_of(s2);
                // The from-side cell of rank r2: find the member
                // adjacent to from_atom with the smallest base rank.
                let mut pos2 = r2 as usize;
                while pos2 > 0 {
                    pos2 -= 1;
                    if state.from.rank_of(state.from.at_number()[pos2]) != r2 {
                        break;
                    }
                    let n2 = state.from.at_number()[pos2];
                    if adjacency[from_atom].contains(&n2)
                        && base_ranks[n2] < cr2_best
                        && base_ranks[n2] > canon_rank2_inp
                    {
                        cr2_best = base_ranks[n2];
                    }
                }
            }
        }
        if cr2_best >= cr1 {
            cr1 += 1;
            canon_rank2_inp = 0;
        } else {
            found = Some((cr1, cr2_best));
            break;
        }
    }
    let Some((cr1, cr2)) = found else {
        return false;
    };
    // The CT compare: reject the level when the enumerated entry is
    // larger than the pre-filled one (c < 0: stored < candidate).
    let (stored1, stored2) = ct[entry];
    if stored1 < cr1 || (stored1 == cr1 && stored2 < cr2) {
        return false;
    }
    let from_atom1 = canonical_atom(base_ranks, cr1 as usize - 1);
    let from_atom2 = canonical_atom(base_ranks, cr2 as usize - 1);
    let at_rank1 = state.from.rank_of(from_atom1);
    // The to-side candidates: the cell of the from-atom's class rank,
    // positions top-down (iMax = at_rank1 - 1 in the reference).
    let mut pos = at_rank1 as usize;
    while pos > 0 {
        pos -= 1;
        if state.to.rank_of(state.to.at_number()[pos]) != at_rank1 {
            break;
        }
        let to_atom1 = state.to.at_number()[pos];
        if bond_neighbors[to_atom1].is_empty() {
            continue;
        }
        for &to_atom2 in &bond_neighbors[to_atom1] {
            if state.from.rank_of(from_atom2) != state.to.rank_of(to_atom2) {
                continue;
            }
            // All_SB_Same: when every equivalent candidate bond around
            // the cell carries the same marks, the mapping (and its
            // refinements) is skipped entirely.
            let all_same = all_sb_same(state, base_ranks, cr1, cr2, bond_neighbors, adjacency);
            let saved = state.clone();
            if !all_same {
                map_an_atom(state, from_atom1, to_atom1, adjacency);
                // The endpoint-2 rank consistency after endpoint-1's
                // refinement (the charge/bond-isomerism guard).
                if state.from.rank_of(from_atom2) != state.to.rank_of(to_atom2) {
                    *state = saved;
                    continue;
                }
                map_an_atom(state, from_atom2, to_atom2, adjacency);
                // The bond's to-side counters drop while the sub-mapping
                // runs (an atom with two marks cannot map twice).
                state.used[to_atom1] -= 1;
                state.used[to_atom2] -= 1;
            }
            let mapped = map_bonds(
                adjacency,
                base_ranks,
                ct,
                bond_neighbors,
                state,
                entry + 1,
                numbering,
            );
            if !all_same {
                state.used[to_atom1] += 1;
                state.used[to_atom2] += 1;
            }
            *state = saved;
            if mapped {
                return true;
            }
        }
    }
    false
}

/// `map_an_atom2`: the shared-rank pin. When the atom's class holds
/// more than one member, both sides' respective atoms pin to the
/// cell's bottom rank (r - i1 + 1, i1 the cell size walked from
/// position r - 2 down) and the equitable refinement cascades.
fn map_an_atom(state: &mut Walk, from_atom: usize, to_atom: usize, adjacency: &[Vec<usize>]) {
    let r = state.from.rank_of(from_atom);
    if r != state.to.rank_of(to_atom) {
        return;
    }
    let imax = r as usize - 1;
    let mut i1 = 1usize;
    while i1 <= imax && state.from.rank_of(state.from.at_number()[imax - i1]) == r {
        i1 += 1;
    }
    if i1 > 1 {
        let new_rank = r - i1 as u32 + 1;
        state.from.pin_and_refine(from_atom, new_rank, adjacency);
        state.to.pin_and_refine(to_atom, new_rank, adjacency);
    }
}

/// `All_SB_Same`: true when every to-cell member's same-rank neighbor
/// pairs carry the marks (so the mapping would repeat the same parity
/// and the reference skips it); false when some pair is unmarked (a
/// mapping onto a non-stereo bond exists). The reference's
/// program-error path resolves to `false` — the mapping proceeds.
fn all_sb_same(
    state: &Walk,
    base_ranks: &[u32],
    cr1: u32,
    cr2: u32,
    bond_neighbors: &[Vec<usize>],
    adjacency: &[Vec<usize>],
) -> bool {
    let n1 = canonical_atom(base_ranks, cr1 as usize - 1);
    let r1 = state.from.rank_of(n1);
    let r2 = state
        .from
        .rank_of(canonical_atom(base_ranks, cr2 as usize - 1));
    // The first to-cell member carrying a marked partner of rank r2.
    let mut seed = None;
    let mut pos = r1 as usize;
    while pos > 0 {
        pos -= 1;
        if state.to.rank_of(state.to.at_number()[pos]) != r1 {
            break;
        }
        let s1 = state.to.at_number()[pos];
        if bond_neighbors[s1]
            .iter()
            .any(|&s2| state.to.rank_of(s2) == r2)
        {
            seed = Some(s1);
            break;
        }
    }
    if seed.is_none() {
        return false;
    }
    // Every to-cell member: each rank-r2 neighbor pair must be a
    // marked bond (both directions), else the equivalence fails.
    let mut num_equal = 0;
    let mut pos = r1 as usize;
    while pos > 0 {
        pos -= 1;
        if state.to.rank_of(state.to.at_number()[pos]) != r1 {
            break;
        }
        let atom = state.to.at_number()[pos];
        for &neighbor in &adjacency[atom] {
            if state.to.rank_of(neighbor) != r2 {
                continue;
            }
            if !bond_neighbors[atom].contains(&neighbor)
                || !bond_neighbors[neighbor].contains(&atom)
            {
                return false;
            }
            num_equal += 1;
        }
    }
    num_equal > 0
}
/// `BreakAllTies`: walk the canonical order; each adjacent tied pair
/// pins its first member to the current position rank and re-refines.
fn break_all_ties(partition: &mut Partition, adjacency: &[Vec<usize>]) {
    let n = partition.at_number().len();
    for i in 1..n {
        let a = partition.at_number()[i - 1];
        let b = partition.at_number()[i];
        if partition.rank_of(a) == partition.rank_of(b) {
            partition.pin_and_refine(a, i as u32, adjacency);
        }
    }
}
