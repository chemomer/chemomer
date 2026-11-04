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

/// The strict stereo-bond candidacy gate — the reference's
/// `set_stereo_bonds_parity` counting core over the mark_alt
/// alternation state. The bond types the normalization leaves
/// ([`super::tautomer::mark_alt_bond_types`]) drive every reading:
/// through `get_allowed_stereo_bond_type` (1.07.5 flags) a DOUBLE (2)
/// or an ALTERN (4) bond is a stereo-bond candidate, SINGLE (1) and
/// TAUTOM (8) pass unmarked, and TRIPLE/ALT_123/ALT_13/ALT_23/ALT12NS
/// (3/5/6/7/9) read ZERO — the ONE_BAD wrong-type whose neighbors a
/// valence-2 atom may carry once, never twice. The gates, in the
/// reference's order: the atom windows ([2,3] heavy + drawn H), the
/// element gate (C/Si/Ge neutral, N neutral or +1), the c-point
/// rejection, the middle-cumulene early out; per candidate the
/// opposite end's neighborhood (wrong-typed bonds counted, exactly one
/// of alt/double present, at most two doubles with the N(V) exception),
/// the type_N oxide pair, the allene chain walk, and the parity loop's
/// ring gate — `is_bond_in_Nmax_memb_ring` only when both ends share
/// the MarkRingSystemsInp ring system, rejecting cycles under
/// MIN_SB_RING_SIZE = 8 (a spiro-split or terminal bond evades it
/// whole). Zero coordinates leave every surviving mark at
/// AB_PARITY_UNDF. Returns the marks as unordered atom pairs.
pub(crate) fn undf_bond_marks(
    structure: &Structure,
    counts: &[u32],
    ranks: &[u32],
    classes: &[u32],
) -> Vec<(usize, usize)> {
    use super::tautomer::endpoint_valence;
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
    let types = super::tautomer::mark_alt_bond_types(structure, counts);
    let c_point = super::tautomer::settled_c_points(structure, counts);
    let ring = super::tautomer::ring_systems(structure);
    let neighbors = adjacency_of(structure);
    let bond_index = |atom: usize, other: usize| -> Option<usize> {
        structure.bonds.iter().position(|&(from, to, _)| {
            (from == atom && to == other) || (to == atom && from == other)
        })
    };
    // get_allowed_stereo_bond_type (1.07.5: ADD_6MEMB_AROM_BOND_PARITY
    // on, TAUTOM attachments allowed).
    let allowed = |bond: usize| -> u32 {
        match types[bond] {
            1 | 2 | 4 | 8 => u32::from(types[bond]),
            _ => 0,
        }
    };
    // chem_bonds_valence over the normalized types, doubled to stay
    // integral: 1/2/3 plain, the alternating blends 1.5 (bt4), the
    // range mids (bt5/bt6: 2, bt7: 2.5), TAUTOM/ALT12NS 1.5 — the
    // MABT-deck readings (a pyrophosphate P reads cbv 5 over 1+8+8+1).
    let blend2 = |bond: usize| -> i32 {
        match types[bond] {
            1 => 2,
            2 => 4,
            3 => 6,
            4 => 3,
            5 | 6 => 4,
            7 => 5,
            8 | 9 => 3,
            _ => 0,
        }
    };
    let cbv2 = |atom: usize| -> i32 {
        neighbors[atom]
            .iter()
            .filter_map(|&next| bond_index(atom, next))
            .map(blend2)
            .sum()
    };
    let degree = |atom: usize| -> u32 { neighbors[atom].len() as u32 };
    let in_window = |atom: usize| -> bool { (2..=3).contains(&(degree(atom) + counts[atom])) };
    let can_have_sb = |atom: usize| -> bool {
        let element = structure.atoms[atom].element.as_str();
        let charge = structure.atoms[atom].charge;
        matches!(
            (element, charge),
            ("C", 0) | ("Si", 0) | ("Ge", 0) | ("N", 0) | ("N", 1)
        )
    };
    // bIsOxide: a DOUBLE bond to a terminal O-like atom (valence-2
    // endpoint, uncharged, H-free), or a TAUTOM/ALT12NS bond to a
    // terminal valence-2 endpoint.
    let is_oxide = |atom: usize| -> bool {
        neighbors[atom].iter().any(|&next| {
            let Some(bond) = bond_index(atom, next) else {
                return false;
            };
            let other_degree = neighbors[next].len();
            match types[bond] {
                2 => {
                    other_degree == 1
                        && structure.atoms[next].charge == 0
                        && counts[next] == 0
                        && endpoint_valence(&structure.atoms[next].element) == Some(2)
                }
                8 | 9 => {
                    other_degree == 1 && endpoint_valence(&structure.atoms[next].element) == Some(2)
                }
                _ => false,
            }
        })
    };
    // bIsSuitableHeteroInpAtom: 0 the =N-/=NH double neighbor, 1 the
    // =O, -1 nothing.
    let suitable_hetero = |atom: usize| -> i32 {
        let record = &structure.atoms[atom];
        if record.charge != 0 {
            return -1;
        }
        let Some(valence) = endpoint_valence(&record.element) else {
            return -1;
        };
        if cbv2(atom) / 2 + counts[atom] as i32 != valence as i32 {
            return -1;
        }
        match valence {
            2 if counts[atom] == 0 && neighbors[atom].len() == 1 => 0,
            3 if (neighbors[atom].len() == 1 && counts[atom] == 1)
                || (neighbors[atom].len() == 2 && counts[atom] == 0) =>
            {
                1
            }
            _ => -1,
        }
    };
    // The N(V) shape type: -N= (valence 2, cbv 3) or >N(V) (valence 3,
    // cbv 5) on a neutral hydrogen-free endpoint-valence-3 atom.
    let type_n = |atom: usize| -> i32 {
        let record = &structure.atoms[atom];
        if counts[atom] != 0 || record.charge != 0 {
            return 0;
        }
        if endpoint_valence(&record.element) != Some(3) {
            return 0;
        }
        if degree(atom) == 2 && cbv2(atom) == 6 {
            1
        } else if degree(atom) == 3 && cbv2(atom) == 10 {
            2
        } else {
            0
        }
    };
    let middle_cumulene = |atom: usize| -> bool {
        counts[atom] == 0
            && neighbors[atom].len() == 2
            && neighbors[atom]
                .iter()
                .all(|&next| bond_index(atom, next).is_some_and(|bond| allowed(bond) == 2))
    };
    let allene_end = |atom: usize| -> bool {
        matches!(
            (
                structure.atoms[atom].element.as_str(),
                structure.atoms[atom].charge
            ),
            ("C", 0) | ("Si", 0) | ("Ge", 0)
        )
    };
    let mut out = Vec::new();
    for at_1 in 0..structure.atoms.len() {
        if !in_window(at_1) || !can_have_sb(at_1) || c_point[at_1].is_some() {
            continue;
        }
        // The middle-cumulene early out: an =C= atom is no bond-parity
        // source.
        if middle_cumulene(at_1) {
            continue;
        }
        let mut num_2s_1 = 0u32;
        let mut num_alt_1 = 0u32;
        let mut num_wrong_1 = 0u32;
        let mut num_chains = 0u32;
        let mut chain_length = 0usize;
        let mut hetero_n = 0u32;
        let mut hetero_o = 0u32;
        let mut candidates: Vec<usize> = Vec::new();
        let mut abort = false;
        for &start in &neighbors[at_1] {
            let mut at_2 = start;
            let Some(bond) = bond_index(at_1, at_2) else {
                continue;
            };
            match allowed(bond) {
                4 => {
                    num_alt_1 += 1;
                    // An alternation bond joins candidates only within
                    // its ring system.
                    if ring.system[at_1] != ring.system[at_2] || !can_have_sb(at_2) {
                        continue;
                    }
                }
                2 => {
                    num_2s_1 += 1;
                    // The cumulene chain walk: an =C= partner extends
                    // the chain while the far atom stays =Y=; a chain
                    // must end at a terminal-allene atom or the bond
                    // drops out. The chain end REPLACES the candidate
                    // partner.
                    if allene_end(at_1) {
                        let mut cursor = at_2;
                        let mut back = at_1;
                        let mut length = 0usize;
                        let mut b_allene;
                        loop {
                            b_allene = counts[cursor] == 0
                                && neighbors[cursor].len() == 2
                                && neighbors[cursor].iter().all(|&next| {
                                    bond_index(cursor, next).is_some_and(|bond| allowed(bond) == 2)
                                });
                            if !(b_allene && allene_end(cursor)) {
                                break;
                            }
                            let k = neighbors[cursor][0] == back;
                            back = cursor;
                            cursor = neighbors[cursor][if k { 1 } else { 0 }];
                            length += 1;
                        }
                        if length > 0 {
                            num_chains += 1;
                            if b_allene || !allene_end(cursor) {
                                continue;
                            }
                            chain_length = length;
                            at_2 = cursor;
                        }
                    }
                    if chain_length == 0 {
                        let n2sh = suitable_hetero(at_2);
                        if n2sh == 0 {
                            hetero_n += 1;
                        } else if n2sh == 1 {
                            hetero_o += 1;
                        }
                        if !can_have_sb(at_2) {
                            continue;
                        }
                    }
                }
                1 | 8 => continue,
                _ => {
                    num_wrong_1 += 1;
                    continue;
                }
            }
            // The registration end: the higher-numbered atom only.
            if at_1 <= at_2 {
                continue;
            }
            if !in_window(at_2) {
                continue;
            }
            // The opposite end's neighborhood.
            let mut num_2s_2 = 0u32;
            let mut num_alt_2 = 0u32;
            let mut num_wrong_2 = 0u32;
            let mut opposite_hetero_n = 0u32;
            let mut opposite_hetero_o = 0u32;
            let mut found_back = false;
            for &next in &neighbors[at_2] {
                let Some(next_bond) = bond_index(at_2, next) else {
                    continue;
                };
                let next_type = allowed(next_bond);
                if next_type == 0 {
                    num_wrong_2 += 1;
                    continue;
                }
                if next_type == 2 {
                    num_2s_2 += 1;
                    let n2sh = suitable_hetero(next);
                    if n2sh == 0 {
                        opposite_hetero_n += 1;
                    } else if n2sh == 1 {
                        opposite_hetero_o += 1;
                    }
                } else if next_type == 4 {
                    num_alt_2 += 1;
                }
                if next == at_1 {
                    found_back = true;
                }
            }
            if num_wrong_2 > 1 || (num_wrong_2 > 0 && degree(at_2) <= 2) {
                continue;
            }
            if i32::from(num_alt_2 > 0) + i32::from(num_2s_2 > 0) != 1 || !found_back {
                continue;
            }
            if c_point[at_2].is_some() {
                continue;
            }
            if num_2s_2 > 2 {
                continue;
            }
            let (tn1, tn2) = (type_n(at_1), type_n(at_2));
            if tn1 | tn2 == 3 && ((tn1 == 2 && !is_oxide(at_1)) || (tn2 == 2 && !is_oxide(at_2))) {
                continue;
            }
            if num_2s_2 == 2 {
                // The N(V) twin-double exception (the O=C=N / N=C=N /
                // O=N=N shapes) — otherwise two doubles at the opposite
                // end kill the bond (`set_stereo_bonds_parity` has no
                // valence-2 pass-through here).
                let record = &structure.atoms[at_2];
                let nv = opposite_hetero_n | opposite_hetero_o == 1
                    && degree(at_2) + counts[at_2] == 3
                    && cbv2(at_2) / 2 + counts[at_2] as i32 == 5
                    && endpoint_valence(&record.element) == Some(3)
                    && (tn2 == 0 || is_oxide(at_2));
                if !nv {
                    continue;
                }
            }
            if chain_length > 0 && num_alt_2 > 0 {
                abort = true; // no alternation bonds inside cumulenes
                break;
            }
            candidates.push(at_2);
        }
        if abort
            || num_chains > 1
            || num_wrong_1 > 1
            || (num_wrong_1 > 0 && degree(at_1) <= 2)
            || i32::from(num_alt_1 > 0) + i32::from(num_2s_1 > 0) != 1
            || candidates.is_empty()
        {
            continue;
        }
        if chain_length > 0 && (candidates.len() != 1 || num_alt_1 > 0 || chain_length > 2) {
            continue;
        }
        if num_2s_1 > 1 {
            // The at_1 side's N(V) twin-double exception.
            let record = &structure.atoms[at_1];
            let nv = hetero_n | hetero_o == 1
                && type_n(at_1) == 2
                && degree(at_1) + counts[at_1] == 3
                && cbv2(at_1) / 2 + counts[at_1] as i32 == 5
                && endpoint_valence(&record.element) == Some(3);
            if !nv {
                continue;
            }
        }
        // The parity loop's ring gate: is_bond_in_Nmax_memb_ring runs
        // only when both ends share the ring system, and rejects the
        // bonds in cycles under MIN_SB_RING_SIZE (8).
        for &at_2 in &candidates {
            if ring.system[at_1] == ring.system[at_2]
                && super::stereo::in_small_ring_stereo(structure, at_1, at_2)
            {
                continue;
            }
            if std::env::var_os("INCHI_TRACE_BNS").is_some() {
                eprintln!("UNDJPUSH at1={} at2={}", at_1 + 1, at_2 + 1);
            }
            out.push((at_1, at_2));
        }
    }
    // The reference's RemoveKnownNonStereoBondParities fixpoint
    // (ichicans.c, post-canonicalization): a mark whose registering
    // end carries THREE heavy neighbors and whose OTHER two neighbors
    // are constitutionally equivalent (in a different ring system than
    // the atom, over parity-free equivalent branches) cannot speak
    // stereo — the C=N+(Me)2 iminium mark dies, the valence-2 imine
    // and azo marks survive.
    prune_known_non_stereo(structure, &mut out, classes, &ring.system, ranks);
    out
}

/// `PathsHaveIdenticalKnownParities` over parity-free branches: a
/// parallel DFS pairing both branches' next neighbors in ascending
/// canon-rank order; every step's pair must share the class rank and
/// visited state, and an edge that crosses a mark fails (the marks'
/// parities are UNDF — never well-defined — so the reference's
/// well-defined-parity requirement rejects the crossing).
#[allow(clippy::too_many_arguments)]
fn paths_identical(
    neighbors: &[Vec<usize>],
    classes: &[u32],
    marks: &std::collections::BTreeSet<(usize, usize)>,
    ranks: &[u32],
    prev1: usize,
    cur1: usize,
    prev2: usize,
    cur2: usize,
    visited1: &mut [u32],
    visited2: &mut [u32],
    length: u32,
) -> u32 {
    let length = length + 1;
    visited1[cur1] = length;
    visited2[cur2] = length;
    if neighbors[cur1].len() != neighbors[cur2].len() {
        return 0;
    }
    if neighbors[cur1].len() == 1 {
        return length;
    }
    // The paired-neighbor iterator: the smallest canon rank above the
    // last, the parent excluded.
    let next_of = |cur: usize, prev: usize, last: &mut u32| -> Option<usize> {
        neighbors[cur]
            .iter()
            .copied()
            .filter(|&next| next != prev)
            .map(|next| (ranks[next], next))
            .filter(|&(rank, _)| rank > *last)
            .min_by_key(|&(rank, _)| rank)
            .map(|(rank, next)| {
                *last = rank;
                next
            })
    };
    let edge_marked =
        |cur: usize, next: usize| -> bool { marks.contains(&(cur.min(next), cur.max(next))) };
    let mut last1 = 0u32;
    let mut last2 = 0u32;
    let steps = neighbors[cur1].len() - 1;
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for _ in 0..steps {
        let (Some(n1), Some(n2)) = (
            next_of(cur1, prev1, &mut last1),
            next_of(cur2, prev2, &mut last2),
        ) else {
            return 0;
        };
        if classes[n1] != classes[n2] || visited1[n1] != visited2[n2] {
            return 0;
        }
        if edge_marked(cur1, n1) || edge_marked(cur2, n2) {
            // Both edges would have to be marks with well-defined
            // parities — the UNDF marks never are.
            return 0;
        }
        pairs.push((n1, n2));
    }
    for (n1, n2) in pairs {
        if visited1[n1] == 0
            && paths_identical(
                neighbors, classes, marks, ranks, cur1, n1, cur2, n2, visited1, visited2, length,
            ) == 0
        {
            return 0;
        }
    }
    length
}

/// The removal fixpoint.
fn prune_known_non_stereo(
    structure: &Structure,
    out: &mut Vec<(usize, usize)>,
    classes: &[u32],
    ring_system: &[usize],
    ranks: &[u32],
) {
    let neighbors = adjacency_of(structure);
    let n = structure.atoms.len();
    loop {
        let marks: std::collections::BTreeSet<(usize, usize)> = out.iter().copied().collect();
        let mut removed = false;
        for index in (0..out.len()).rev() {
            let (a, b) = out[index];
            for &(atom, other) in &[(a, b), (b, a)] {
                if neighbors[atom].len() != 3 {
                    continue;
                }
                let rest: Vec<usize> = neighbors[atom]
                    .iter()
                    .copied()
                    .filter(|&next| next != other)
                    .collect();
                if rest.len() != 2 || classes[rest[0]] != classes[rest[1]] {
                    continue;
                }
                if ring_system[atom] == ring_system[rest[0]] {
                    continue;
                }
                let mut visited1 = vec![0u32; n];
                let mut visited2 = vec![0u32; n];
                if paths_identical(
                    &neighbors,
                    classes,
                    &marks,
                    ranks,
                    atom,
                    rest[0],
                    atom,
                    rest[1],
                    &mut visited1,
                    &mut visited2,
                    0,
                ) > 0
                {
                    out.remove(index);
                    removed = true;
                    break;
                }
            }
            if removed {
                break;
            }
        }
        if !removed {
            break;
        }
    }
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
    let base = Walk {
        from: Partition::initial_pair(&seed),
        to: Partition::initial_pair(&seed),
    };
    // The stack of mapping states: slot 0 is the seed, each accepted
    // pin pushes a refined copy of both sides (`map_an_atom2`'s
    // pRankStack push, `bAddStack`).
    let mut stacks: Vec<Walk> = vec![base];
    let mut used: Vec<usize> = bond_neighbors.iter().map(Vec::len).collect();
    let mut numbering = None;
    let mut better = false;
    if map_bonds(
        adjacency,
        base_ranks,
        &mut ct,
        &bond_neighbors,
        &mut stacks,
        &mut used,
        0,
        0,
        &mut better,
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

/// One stack slot: the two mapping sides (`pRankStack1`/`pRankStack2`
/// rank+atomNumber pairs). Slots above the seed are pushed copies
/// carrying accepted pins; the walk reads the slot its level's stack
/// pointer names.
#[derive(Clone)]
struct Walk {
    from: Partition,
    to: Partition,
}

/// The bond-level mapping (map_stereo_bonds4's body): enumerate the
/// next CT entry, walk the to-side candidates, run the acceptance
/// state machine per candidate — pin-and-refine through stack-slot
/// pushes (the to-side bond counters drop while a sub-mapping runs),
/// recurse. `numbering` receives the FIRST complete mapping's result
/// (`bFirstCT`); `better` threads the reference's `bStereoIsBetter`
/// (an earlier entry improved, so later entries compare unbounded).
///
/// The acceptance mechanics (witnessed on the biguanide fixture's
/// single-pin trace): a candidate whose cell bonds all carry the same
/// marks (`All_SB_Same`, over the level's slot) accepts WITHOUT any
/// pin — the recursion rides the slot unchanged. Otherwise the first
/// endpoint pins once per candidate first-atom (the reference's
/// `LastMappedTo1` keeps the pushed slot across the second-endpoint
/// candidates of the same first atom), a second `All_SB_Same` over the
/// PINNED slot skips only the second endpoint's pin, and the known
/// parity (UNDF never reaches the calc branch) accepts with the
/// endpoint pins the state machine left behind. The CT compare gates
/// each accept: a candidate past the stored entry rejects unless an
/// earlier entry improved, and an accepted smaller entry raises
/// `better` for everything above it — the lexicographic discipline
/// that let the reference accept its second bond past a worse stored
/// entry while the first stayed better.
#[allow(clippy::too_many_arguments)]
fn map_bonds(
    adjacency: &[Vec<usize>],
    base_ranks: &[u32],
    ct: &mut [(u32, u32)],
    bond_neighbors: &[Vec<usize>],
    stacks: &mut Vec<Walk>,
    used: &mut [usize],
    base: usize,
    entry: usize,
    better: &mut bool,
    numbering: &mut Option<Vec<u32>>,
) -> bool {
    let n = adjacency.len();
    if entry >= ct.len() {
        // All bonds mapped: BreakAllTies on both sides, then the
        // position-wise copy — the to-side order carrying the base
        // ranks of the from-side atoms at the same positions
        // (nCanonRankTo[n2] = nCanonRankFrom[n1]).
        let mut from = stacks[base].from.clone();
        let mut to = stacks[base].to.clone();
        break_all_ties(&mut from, adjacency);
        break_all_ties(&mut to, adjacency);
        let mut ranks = vec![0u32; n];
        for position in 0..n {
            ranks[to.at_number()[position]] = base_ranks[from.at_number()[position]];
        }
        *numbering = Some(ranks);
        return true;
    }
    let prev_bond = ct[entry];
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
        let r1 = stacks[base].from.rank_of(from_atom);
        let mut cr2_best = n as u32;
        // The to-side cell of rank r1, positions top-down.
        let mut pos = r1 as usize;
        while pos > 0 {
            pos -= 1;
            if stacks[base].to.rank_of(stacks[base].to.at_number()[pos]) != r1 {
                break;
            }
            let s1 = stacks[base].to.at_number()[pos];
            if used[s1] == 0 {
                continue;
            }
            for &s2 in &bond_neighbors[s1] {
                if used[s2] == 0 {
                    continue;
                }
                let r2 = stacks[base].to.rank_of(s2);
                // The from-side cell of rank r2: find the member
                // adjacent to from_atom with the smallest base rank.
                let mut pos2 = r2 as usize;
                while pos2 > 0 {
                    pos2 -= 1;
                    if stacks[base].from.rank_of(stacks[base].from.at_number()[pos2]) != r2 {
                        break;
                    }
                    let n2 = stacks[base].from.at_number()[pos2];
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
    // The level bound: an enumerated pair past the stored entry rejects
    // the level — unless an earlier entry improved (bStereoIsBetter
    // lifts the bound for the rest of the CT).
    if !*better && (ct[entry].0 < cr1 || (ct[entry].0 == cr1 && ct[entry].1 < cr2)) {
        ct[entry] = prev_bond;
        return false;
    }
    let from_atom1 = canonical_atom(base_ranks, cr1 as usize - 1);
    let from_atom2 = canonical_atom(base_ranks, cr2 as usize - 1);
    let at_rank1 = stacks[base].from.rank_of(from_atom1);
    // The candidate scan: the to-side cell of the from-atom's class
    // rank, positions top-down (iMax = at_rank1 - 1). `last_mapped_to1`
    // keeps the first endpoint's pushed slot across the second-endpoint
    // candidates of the same first atom (the reference's
    // LastMappedTo1); `sp1` is that slot pointer (nStackPtr[1]).
    let mut last_mapped_to1: Option<usize> = None;
    let mut sp1 = base;
    let mut pos = at_rank1 as usize;
    while pos > 0 {
        pos -= 1;
        if stacks[base].to.rank_of(stacks[base].to.at_number()[pos]) != at_rank1 {
            break;
        }
        let to_atom1 = stacks[base].to.at_number()[pos];
        if bond_neighbors[to_atom1].is_empty() {
            continue; // at_to1 does not belong to a stereo bond
        }
        for &to_atom2 in &bond_neighbors[to_atom1] {
            if stacks[base].from.rank_of(from_atom2) != stacks[base].to.rank_of(to_atom2) {
                continue; // at_from2 cannot be mapped on at_to2
            }
            // All_SB_Same over the level's slot: when every equivalent
            // candidate bond around the cell carries the same marks the
            // mapping runs with NO pins at all (the reference's outer
            // same-parity branch — ClearPreviousMappings and recurse on
            // the slot unchanged).
            let all_same =
                all_sb_same(&stacks[base], base_ranks, cr1, cr2, bond_neighbors, adjacency);
            let mut istk = base;
            let mut same2 = false;
            if !all_same {
                if last_mapped_to1 != Some(to_atom1) {
                    // map_an_atom2 for the first endpoint: a tied cell
                    // pushes a refined copy of both sides; a singleton
                    // pins nothing (bAddStack = 0).
                    sp1 = match number_of_ties(&stacks[istk], from_atom1, to_atom1) {
                        Err(()) => return false,
                        Ok(None) => istk,
                        Ok(Some(new_rank)) => {
                            let pushed = pin_slot(stacks, istk, from_atom1, to_atom1, new_rank, adjacency);
                            pushed
                        }
                    };
                    last_mapped_to1 = Some(to_atom1);
                }
                istk = sp1;
                // The endpoint-2 rank consistency after endpoint-1's
                // refinement (the charge/bond-isomerism guard).
                if stacks[istk].from.rank_of(from_atom2) != stacks[istk].to.rank_of(to_atom2) {
                    continue;
                }
                // The second All_SB_Same over the PINNED slot: skip
                // only the second endpoint's pin.
                if istk > base {
                    same2 = all_sb_same(&stacks[istk], base_ranks, cr1, cr2, bond_neighbors, adjacency);
                }
                if !same2 {
                    istk = match number_of_ties(&stacks[istk], from_atom2, to_atom2) {
                        Err(()) => return false,
                        Ok(None) => istk,
                        Ok(Some(new_rank)) => {
                            pin_slot(stacks, istk, from_atom2, to_atom2, new_rank, adjacency)
                        }
                    };
                }
            }
            // The known-parity accept (UNDF never reaches the calc
            // branch): compare the candidate against the stored entry —
            // reject when past it unless an earlier entry improved;
            // write the entry; recurse on the state the pins left.
            let c = compare_ct(ct[entry], cr1, cr2);
            if c < 0 && !*better {
                continue; // reject: not a minimal CT
            }
            let mut better_was_set_here = false;
            let mut prev_bond2 = prev_bond;
            if c > 0 && !*better {
                *better = true;
                better_was_set_here = true;
                prev_bond2 = ct[entry];
            }
            ct[entry] = (cr1, cr2);
            if !all_same {
                // the to-side counters drop only around a pinning
                // sub-mapping (the reference guards with
                // !bAllParitiesIdentical)
                used[to_atom1] -= 1;
                used[to_atom2] -= 1;
            }
            let mapped = map_bonds(
                adjacency,
                base_ranks,
                ct,
                bond_neighbors,
                stacks,
                used,
                istk,
                entry + 1,
                better,
                numbering,
            );
            if !all_same {
                used[to_atom1] += 1;
                used[to_atom2] += 1;
            }
            if mapped {
                return true;
            }
            if better_was_set_here {
                // rollback: this level's improvement died with the
                // failed sub-mapping
                *better = false;
                ct[entry] = prev_bond2;
            }
            if all_same || same2 {
                break; // j2 cycle: at_to2 needs no repeat
            }
        }
    }
    // The level failed: restore the stored entry (prevBond).
    ct[entry] = prev_bond;
    false
}

/// `map_an_atom2`'s pin: both sides refine onto fresh stack slots (the
/// pushed copy; a repeat pin at the same depth overwrites the scratch
/// slot, the way the reference's stack reuse does).
fn pin_slot(
    stacks: &mut Vec<Walk>,
    slot: usize,
    from_atom: usize,
    to_atom: usize,
    new_rank: u32,
    adjacency: &[Vec<usize>],
) -> usize {
    let mut pushed = stacks[slot].clone();
    pushed.from.pin_and_refine(from_atom, new_rank, adjacency);
    pushed.to.pin_and_refine(to_atom, new_rank, adjacency);
    let next = slot + 1;
    if next < stacks.len() {
        stacks[next] = pushed;
    } else {
        stacks.push(pushed);
    }
    next
}

/// `NumberOfTies`: the shared mapping rank's cell size, and the pin
/// rank (r - i1 + 1) when the cell holds a tie. Both sides' cells must
/// agree on the size (the reference's CT_MAPCOUNT_ERR guards); a
/// singleton cell pins nothing.
fn number_of_ties(
    state: &Walk,
    from_atom: usize,
    to_atom: usize,
) -> Result<Option<u32>, ()> {
    let r = state.from.rank_of(from_atom);
    if r != state.to.rank_of(to_atom) {
        return Err(()); // atoms cannot be mapped onto each other
    }
    let imax = r as usize - 1;
    let mut i1 = 1usize;
    while i1 <= imax && state.from.rank_of(state.from.at_number()[imax - i1]) == r {
        i1 += 1;
    }
    let mut i2 = 1usize;
    while i2 <= imax && state.to.rank_of(state.to.at_number()[imax - i2]) == r {
        i2 += 1;
    }
    if i2 != i1 {
        return Err(()); // must be identical number of equal ranks
    }
    if i1 > 1 {
        Ok(Some(r - i1 as u32 + 1))
    } else {
        Ok(None)
    }
}

/// The stored-entry compare (`CompareLinCtStereoDoubleToValues` over
/// the at_num1/at_num2 pair; the parity byte is UNDF on both sides of
/// every corpus comparison): 1 when the stored entry is greater (the
/// candidate improves it), -1 when past it, 0 when equal.
fn compare_ct(stored: (u32, u32), cr1: u32, cr2: u32) -> i32 {
    if stored.0 > cr1 {
        return 1;
    }
    if stored.0 < cr1 {
        return -1;
    }
    if stored.1 > cr2 {
        return 1;
    }
    if stored.1 < cr2 {
        return -1;
    }
    0
}

/// `All_SB_Same` (the UNDF simplification): true when every to-cell
/// member's same-rank neighbor pairs carry the marks (so a pinning
/// mapping would repeat the same parity); false when some pair is
/// unmarked (a mapping onto a non-stereo bond exists). The reference's
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
