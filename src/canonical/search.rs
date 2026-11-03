//! The canonicalization search, mirroring the reference's base flow
//! (ichican2.c `GetBaseCanonRanking` + `CanonGraph`).
//!
//! Stage 1 is hydrogen-blind: invariants, equitable refinement,
//! individualization, and the smallest connection table wins
//! (`CtFullCompare` under `INCHI_CANON_MIN`). Every table-equal
//! terminal is an automorphic image of the first, and the position-wise
//! mappings merge into orbits; each orbit's minimal canonical rank is
//! its symmetry class. That partition is then rebuilt from the classes
//! and searched again when the rebuild changed anything — the
//! reference's `FixCanonEquivalenceInfo` + second `CanonGraph` call,
//! which is what makes the result independent of the input atom order.
//!
//! Stage 2 (the reference's `CanonGraph03`) is the hydrogen-aware pass,
//! grown from the partition `FixCanonEquivalenceInfo` left behind: the
//! auxiliary `(symmetry, encoded H)` refinement orders candidates
//! within cells (its `CellGetMinNode` with `nAuxRank`), and the
//! terminal compare adds the encoded hydrogen layer after the
//! connection table. It runs only when some symmetry class holds
//! different hydrogen counts, or when a mobile group needs the t-group
//! rows; otherwise the H-free numbering already is canonical and the
//! reference copies it verbatim.

use super::partition::Partition;
use super::table::ConnectionTable;

/// The reference's `numHNoTautH` encoding (`BASE_H_NUMBER` /
/// `EMPTY_H_NUMBER`, ichican2.c): a real count is
/// `count + BASE_H_NUMBER`, and an atom with no attached hydrogen is
/// `EMPTY_H_NUMBER`, which is LARGER than every real count. The
/// terminal compare minimizes the encoded vector, so a zero-H atom
/// sorts after every hydrogen-bearing atom — the opposite of comparing
/// raw counts, and the reason the canonical numbering drives the
/// unsaturation last.
pub(crate) const BASE_H_NUMBER: u32 = 16_383;
pub(crate) const EMPTY_H_NUMBER: u32 = 32_766;

pub(crate) fn encoded_hydrogen(hydrogens: u32) -> u32 {
    if hydrogens == 0 {
        EMPTY_H_NUMBER
    } else {
        hydrogens + BASE_H_NUMBER
    }
}

/// Development tracing, enabled with `INCHI_TRACE_CANON=1`. The output
/// mirrors the reference's `PartitionColorVertex` / L-loop dumps so the
/// two search trees can be diffed level by level.
fn tracing() -> bool {
    static TRACING: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *TRACING.get_or_init(|| std::env::var_os("INCHI_TRACE_CANON").is_some())
}

fn trace_partition(tag: &str, partition: &Partition) {
    if !tracing() {
        return;
    }
    let n = partition.ranks().len();
    let mut line = String::new();
    line.push_str(tag);
    line.push_str(" rank=");
    for atom in 0..n {
        line.push_str(&format!("{},", partition.rank_of(atom)));
    }
    line.push_str(" atnum=");
    for &atom in partition.at_number() {
        line.push_str(&format!("{},", atom + 1));
    }
    eprintln!("{line}");
}

/// The H-free canonical search plus the reference's fix-and-rerun: the
/// winning numbering's orbit minima become symmetry classes; those
/// classes build a fresh seed partition; the search runs again from it
/// whenever the rebuild changed the partition (`bChangeCurrRank`) or
/// the stored class minima were not the contiguous class ranks
/// (`bChangeSymmRank`). The reference does not iterate this fix; one
/// rerun is all `FixCanonEquivalenceInfo` triggers.
///
/// `h_seed` is the class partition as `FixCanonEquivalenceInfo` leaves
/// it — built from the FIRST pass's orbit minima. The reference hands
/// exactly that partition (never the rerun's own numbering) to the
/// hydrogen-aware `CanonGraph03`, while the rerun's numbering supplies
/// the symmetry minima the aux ranks are built from.
struct HfreeCanon {
    final_search: StageOne,
    h_seed: Partition,
    /// The partition the final H-free walk started from — the rerun's
    /// `h_seed` when the fix triggered, the first refined partition
    /// otherwise. The stereo-informed H-free walk starts here too.
    root: Partition,
}

fn hfree_canonical(adjacency: &[Vec<usize>], invariants: &[(u16, u16)]) -> HfreeCanon {
    hfree_canonical_bounded(adjacency, invariants, adjacency.len())
}

/// [`hfree_canonical`] with the search's individualization bounded to
/// the first `atom_bound` nodes — the t-group-tail variant never
/// branches on the tail cells.
fn hfree_canonical_bounded(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    atom_bound: usize,
) -> HfreeCanon {
    let mut initial = Partition::initial(invariants);
    initial.refine(adjacency);
    let first = search_stage_one(adjacency, &initial, "S1 enter pi0", atom_bound);
    let first_symmetry = orbit_min_symmetry(&first);
    let h_seed = class_partition(&first_symmetry);
    // `bChangeCurrRank` compares the rebuilt class ranks against the
    // ORIGINAL seed's ranks (the reference never re-runs the fix, so
    // the comparison baseline stays the first equitable partition).
    let ranks_changed = h_seed.ranks() != initial.ranks();
    let final_search = if ranks_changed || !class_ranks_are_contiguous(&first_symmetry) {
        search_stage_one(adjacency, &h_seed, "S1 rerun pi0", atom_bound)
    } else {
        first
    };
    let root = if ranks_changed || !class_ranks_are_contiguous(&first_symmetry) {
        h_seed.clone()
    } else {
        initial.clone()
    };
    HfreeCanon {
        final_search,
        h_seed,
        root,
    }
}

/// Orbit-minimal canonical rank per atom — the reference's
/// `nSymmRank`: for every orbit, the smallest winning rank among its
/// members, shared by the whole orbit.
fn orbit_min_symmetry(base: &StageOne) -> Vec<u32> {
    let n = base.partition.ranks().len();
    let mut min_by_root = vec![u32::MAX; n];
    for atom in 0..n {
        let root = base.orbits.representative(atom);
        min_by_root[root] = min_by_root[root].min(base.partition.rank_of(atom));
    }
    (0..n)
        .map(|atom| min_by_root[base.orbits.representative(atom)])
        .collect()
}

/// The seed partition rebuilt from symmetry classes
/// (`SortedEquInfoToRanks`): classes ordered by their minimal rank,
/// atoms within a class by atom number, ranks the cell-end positions.
fn class_partition(symmetry: &[u32]) -> Partition {
    let invariants: Vec<(u32, u32)> = symmetry.iter().map(|&value| (value, 0)).collect();
    Partition::initial_pair(&invariants)
}

/// `FixCanonEquivalenceInfo`'s `bChangeSymmRank`: true when a class's
/// stored minimum is not the rank its class has in the rebuilt
/// partition (class starting at position `i` must carry minimum `i+1`).
fn class_ranks_are_contiguous(symmetry: &[u32]) -> bool {
    let mut classes: Vec<(u32, usize)> = symmetry
        .iter()
        .enumerate()
        .map(|(atom, &value)| (value, atom))
        .collect();
    classes.sort_unstable();
    let mut start = 0usize;
    let mut index = 0usize;
    while index < classes.len() {
        let value = classes[index].0;
        if value as usize != start + 1 {
            return false;
        }
        while index < classes.len() && classes[index].0 == value {
            index += 1;
        }
        start = index;
    }
    true
}

/// True when every symmetry class carries one hydrogen count. The
/// reference skips the hydrogen-aware canonicalization in that case:
/// no class can split on hydrogen, so the H-free numbering stands.
fn hydrogens_uniform(symmetry: &[u32], hydrogens: &[u32]) -> bool {
    let mut classes: Vec<(u32, u32)> = symmetry
        .iter()
        .zip(hydrogens)
        .map(|(&class, &count)| (class, count))
        .collect();
    classes.sort_unstable();
    classes
        .windows(2)
        .all(|pair| pair[0].0 != pair[1].0 || pair[0].1 == pair[1].1)
}

/// Symmetry classes per atom: the orbit-minimal base ranks the stereo
/// layer reads orbits from, after the reference's rebuild-and-rerun.
pub fn base_symmetry(adjacency: &[Vec<usize>], invariants: &[(u16, u16)]) -> Vec<u32> {
    let hfree = hfree_canonical(adjacency, invariants);
    orbit_min_symmetry(&hfree.final_search)
}

/// Final canonical rank (1-based) per atom, plus the orbit-minimal
/// symmetry class per atom of the same search — the pair the
/// identifier layers and the aux equivalence classes speak. The
/// symmetry is H-aware: it comes from the stage-two orbits when the
/// hydrogen-aware pass runs, from the H-free orbits otherwise.
pub fn canonical_ranks_and_symmetry(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
) -> (Vec<u32>, Vec<u32>) {
    let n = adjacency.len();
    let hfree = hfree_canonical(adjacency, invariants);
    let symmetry = orbit_min_symmetry(&hfree.final_search);
    // No class splits on hydrogen and no mobile group needs the t-group
    // rows: the H-free numbering is the final answer, exactly the
    // reference's copy of nCanonRankNoTautH.
    if groups.is_empty() && hydrogens_uniform(&symmetry, hydrogens) {
        let ranks = (0..n)
            .map(|atom| hfree.final_search.partition.rank_of(atom))
            .collect();
        return (ranks, symmetry);
    }
    // Stage-2 auxiliary partition: (symmetry rank, encoded hydrogen
    // count), then equitable refinement — the reference's nRankAux.
    let aux_invariants: Vec<(u32, u32)> = (0..n)
        .map(|atom| (symmetry[atom], encoded_hydrogen(hydrogens[atom])))
        .collect();
    let mut aux = Partition::initial_pair(&aux_invariants);
    aux.refine(adjacency);
    let aux_rank: Vec<u32> = (0..n).map(|atom| aux.rank_of(atom)).collect();

    // THE SEARCH TREE GROWS FROM THE FIX PARTITION — the reference's
    // pi[0] handed to CanonGraph03, i.e. the classes built from the
    // FIRST H-free pass (the rerun, when it happens, never replaces
    // them). The (symmetry, H) aux partition only orders the
    // candidates; the classes themselves bound the tree.
    let base = &hfree.h_seed;

    if tracing() {
        let mut line = String::from("S2 symmetry rank=");
        for &value in &symmetry {
            line.push_str(&format!("{value},"));
        }
        line.push_str(" aux=");
        for &value in &aux_rank {
            line.push_str(&format!("{value},"));
        }
        eprintln!("{line}");
        trace_partition("S2 enter pi0", base);
        let mut line = String::from("S2 hydrogens=");
        for &value in hydrogens {
            line.push_str(&format!("{value},"));
        }
        eprintln!("{line}");
    }

    let (informed, informed_orbits) =
        search_stage_two(adjacency, base, &aux_rank, hydrogens, groups, n);
    let ranks: Vec<u32> = (0..n).map(|atom| informed.rank_of(atom)).collect();
    // The stage-two orbits: the walk's theta partition (ichican2.c
    // CanonGraph's UnorderedPartition) — union-find over the cycles of
    // every transposition the walk witnesses between table+hydrogen
    // equal terminals, against the first leaf (L7) and against the
    // best leaf (L8), with the walk trying only candidates that are
    // their own orbit minimum (L14's mcr skip). The symmetry the aux
    // equivalence classes speak is the orbit MINIMUM over theta,
    // mirroring `orbit_min_symmetry`.
    let informed_symmetry = orbit_minima(&ranks, &informed_orbits);
    (ranks, informed_symmetry)
}

/// Orbit-minimal canonical rank per atom over an explicit orbit
/// partition: for every orbit, the smallest winning rank among its
/// members, shared by the whole orbit.
fn orbit_minima(ranks: &[u32], orbits: &Orbits) -> Vec<u32> {
    let n = ranks.len();
    let mut min_by_root = vec![u32::MAX; n];
    for (atom, &rank) in ranks.iter().enumerate() {
        let root = orbits.representative(atom);
        min_by_root[root] = min_by_root[root].min(rank);
    }
    (0..n)
        .map(|atom| min_by_root[orbits.representative(atom)])
        .collect()
}

/// The stereo-informed ranking verdict: ranks and symmetry of the
/// winning pass, plus the inversion comparison `(main_pairs,
/// inv_pairs, sign)` — sign is -1 when the inverted pass won (its
/// numbering was adopted), 0 when the stereo rows came out equal, 1
/// when the absolute pass won on rows alone.
pub type StereoRanking = (Vec<u32>, Vec<u32>, (Vec<(u32, u32)>, Vec<(u32, u32)>, i32));

/// Final canonical rank (1-based) per atom, plus the orbit-minimal
/// symmetry class per atom — with the terminal compare extended by
/// the stereo rows (the reference's stereo canonicalization, whose
/// numbering the aux speaks). Returns the winning pass's ranks and
/// symmetry plus the inversion verdict (see [`StereoRanking`]).
pub fn canonical_ranks_and_symmetry_stereo(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
    centers: &[(usize, Vec<usize>, u32)],
) -> StereoRanking {
    let run = |flip: bool| {
        let (ranks, symmetry) =
            search_with_stereo(adjacency, invariants, hydrogens, groups, centers, flip);
        let pairs = super::stereo::terminal_pairs(centers, &ranks, flip);
        (ranks, symmetry, pairs)
    };
    let (abs_ranks, abs_symmetry, abs_pairs) = run(false);
    let (inv_ranks, inv_symmetry, inv_pairs) = run(true);
    // The reference's diff: centers' canonical numbers first, then
    // parities — the first difference carries the sign; equal rows
    // mean a trivial inversion.
    let mut sign = 0;
    for index in 0..abs_pairs.len().min(inv_pairs.len()) {
        if abs_pairs[index].0 != inv_pairs[index].0 {
            sign = if inv_pairs[index].0 > abs_pairs[index].0 {
                1
            } else {
                -1
            };
            break;
        }
        if abs_pairs[index].1 != inv_pairs[index].1 {
            sign = if inv_pairs[index].1 > abs_pairs[index].1 {
                1
            } else {
                -1
            };
            break;
        }
    }
    if sign < 0 {
        (inv_ranks, inv_symmetry, (inv_pairs, abs_pairs, -1))
    } else {
        (abs_ranks, abs_symmetry, (abs_pairs, inv_pairs, sign))
    }
}

/// The H-aware search with the stereo rows appended to the terminal
/// compare (after the connection table, the encoded hydrogens and the
/// t-group rows). The base flow is mirrored exactly — stage 2 runs
/// when it would for the stereo-free structure, and otherwise the
/// H-FREE tree is walked (atom-number candidate order) with the
/// stereo compare added — so the stereo rows only ever decide ties
/// the base leaves standing.
fn search_with_stereo(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
    centers: &[(usize, Vec<usize>, u32)],
    flip: bool,
) -> (Vec<u32>, Vec<u32>) {
    let n = adjacency.len();
    let hfree = hfree_canonical(adjacency, invariants);
    let symmetry = orbit_min_symmetry(&hfree.final_search);
    let use_stage2 = !groups.is_empty() || !hydrogens_uniform(&symmetry, hydrogens);
    if !use_stage2 {
        let mut state = StereoStageOneState {
            first: None,
            best: None,
            orbits: Orbits::new(n),
        };
        descend_hfree_stereo(&hfree.root, adjacency, centers, flip, &mut state);
        let winner = state.best.expect("a discrete partition is always reached");
        let ranks = (0..n).map(|atom| winner.rank_of(atom)).collect();
        let mut min_by_root = vec![u32::MAX; n];
        for atom in 0..n {
            let root = state.orbits.representative(atom);
            min_by_root[root] = min_by_root[root].min(winner.rank_of(atom));
        }
        let informed_symmetry = (0..n)
            .map(|atom| min_by_root[state.orbits.representative(atom)])
            .collect();
        return (ranks, informed_symmetry);
    }
    let aux_invariants: Vec<(u32, u32)> = (0..n)
        .map(|atom| (symmetry[atom], encoded_hydrogen(hydrogens[atom])))
        .collect();
    let mut aux = Partition::initial_pair(&aux_invariants);
    aux.refine(adjacency);
    let aux_rank: Vec<u32> = (0..n).map(|atom| aux.rank_of(atom)).collect();
    let base = &hfree.h_seed;
    let (informed, informed_orbits) =
        search_stage_two_stereo(adjacency, base, &aux_rank, hydrogens, groups, centers, flip);
    let ranks: Vec<u32> = (0..n).map(|atom| informed.rank_of(atom)).collect();
    // The walk's theta partition, as in the stereo-free path.
    let informed_symmetry = orbit_minima(&ranks, &informed_orbits);
    (ranks, informed_symmetry)
}

/// The H-free walk's state (the stereo-informed variant).
struct StereoStageOneState {
    first: Option<Partition>,
    best: Option<Partition>,
    orbits: Orbits,
}

fn descend_hfree_stereo(
    partition: &Partition,
    adjacency: &[Vec<usize>],
    centers: &[(usize, Vec<usize>, u32)],
    flip: bool,
    state: &mut StereoStageOneState,
) {
    match partition.first_nonsingleton_cell() {
        None => {
            let table = ConnectionTable::build(partition.ranks(), adjacency);
            if let Some(first) = &state.first {
                if table.entries() == ConnectionTable::build(first.ranks(), adjacency).entries() {
                    for position in 0..adjacency.len() {
                        state
                            .orbits
                            .union(first.at_number()[position], partition.at_number()[position]);
                    }
                }
            } else {
                state.first = Some(partition.clone());
            }
            let better = match &state.best {
                None => true,
                Some(best_partition) => {
                    let best_table = ConnectionTable::build(best_partition.ranks(), adjacency);
                    match table.compare(&best_table) {
                        std::cmp::Ordering::Greater => true,
                        std::cmp::Ordering::Less => false,
                        std::cmp::Ordering::Equal => {
                            // The stereo rows keep the reference's own
                            // (minimizing) ordering: witnessed on the
                            // dichlorobutane pair, whose table-equal
                            // terminals the stereo rows split.
                            super::stereo::terminal_rows(centers, partition.ranks(), flip)
                                < super::stereo::terminal_rows(
                                    centers,
                                    best_partition.ranks(),
                                    flip,
                                )
                        }
                    }
                }
            };
            if better {
                state.best = Some(partition.clone());
            }
        }
        Some((start, end)) => {
            let mut candidates: Vec<usize> = partition.at_number()[start..end].to_vec();
            candidates.sort_unstable();
            for v in candidates {
                let mut branch = partition.clone();
                branch.individualize(adjacency, v);
                descend_hfree_stereo(&branch, adjacency, centers, flip, state);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn search_stage_two_stereo(
    adjacency: &[Vec<usize>],
    aux: &Partition,
    aux_rank: &[u32],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
    centers: &[(usize, Vec<usize>, u32)],
    flip: bool,
) -> (Partition, Orbits) {
    let n = adjacency.len();
    let mut best: Option<Partition> = None;
    let mut first: Option<Partition> = None;
    let mut orbits = Orbits::new(n);
    descend_stereo(
        aux,
        adjacency,
        aux_rank,
        hydrogens,
        groups,
        centers,
        flip,
        &mut first,
        &mut best,
        &mut orbits,
    );
    (
        best.expect("a discrete partition is always reached"),
        orbits,
    )
}

#[allow(clippy::too_many_arguments)]
fn descend_stereo(
    partition: &Partition,
    adjacency: &[Vec<usize>],
    aux_rank: &[u32],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
    centers: &[(usize, Vec<usize>, u32)],
    flip: bool,
    first: &mut Option<Partition>,
    best: &mut Option<Partition>,
    orbits: &mut Orbits,
) {
    let cell = match partition.first_nonsingleton_cell() {
        Some(span) => span,
        None => {
            let ranks = partition.ranks();
            let table = ConnectionTable::build(ranks, adjacency);
            let matches = |other: &Partition| {
                table.entries() == ConnectionTable::build(other.ranks(), adjacency).entries()
                    && (0..adjacency.len()).all(|position| {
                        hydrogens[other.at_number()[position]]
                            == hydrogens[partition.at_number()[position]]
                    })
            };
            if let Some(first_partition) = first {
                if matches(first_partition) {
                    for position in 0..adjacency.len() {
                        orbits.union(
                            first_partition.at_number()[position],
                            partition.at_number()[position],
                        );
                    }
                }
            } else {
                *first = Some(partition.clone());
            }
            let counts = |p: &Partition| -> Vec<u32> {
                (0..adjacency.len())
                    .map(|atom| encoded_hydrogen(hydrogens[p.at_number()[atom]]))
                    .collect()
            };
            let group_rows = |p: &Partition| -> Vec<(u32, Vec<u32>)> {
                let ranks = p.ranks();
                groups
                    .iter()
                    .map(|(count, members)| {
                        let mut member_ranks: Vec<u32> =
                            members.iter().map(|&m| ranks[m]).collect();
                        member_ranks.sort_unstable();
                        (*count, member_ranks)
                    })
                    .collect()
            };
            let stereo_rows =
                |p: &Partition| super::stereo::terminal_rows(centers, p.ranks(), flip);
            let better = match &best {
                None => true,
                Some(best_partition) => {
                    let best_table = ConnectionTable::build(best_partition.ranks(), adjacency);
                    match table.compare(&best_table) {
                        std::cmp::Ordering::Greater => true,
                        std::cmp::Ordering::Less => false,
                        std::cmp::Ordering::Equal => {
                            match counts(partition).cmp(&counts(best_partition)) {
                                std::cmp::Ordering::Less => true,
                                std::cmp::Ordering::Greater => false,
                                std::cmp::Ordering::Equal => {
                                    match group_rows(partition).cmp(&group_rows(best_partition)) {
                                        std::cmp::Ordering::Greater => true,
                                        std::cmp::Ordering::Less => false,
                                        std::cmp::Ordering::Equal => {
                                            // The stereo rows keep the
                                            // reference's own (minimizing)
                                            // ordering: witnessed on the
                                            // dichlorobutane pair, whose
                                            // table-equal terminals the
                                            // stereo rows split.
                                            stereo_rows(partition) < stereo_rows(best_partition)
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            };
            if better {
                *best = Some(partition.clone());
            }
            return;
        }
    };
    let mut candidates: Vec<usize> = partition.at_number()[cell.0..cell.1].to_vec();
    candidates.sort_by_key(|&atom| (aux_rank[atom], atom));
    // The L14 theta skip covers only backtracking-iterated candidates;
    // a freshly entered cell's first vertex is the L2 descent's
    // unconditional `CellGetMinNode(..., 0, ...)` (see the stereo-free
    // walk).
    for (index, &v) in candidates.iter().enumerate() {
        if index > 0 && orbits.representative(v) != v {
            continue;
        }
        let mut branch = partition.clone();
        branch.individualize(adjacency, v);
        descend_stereo(
            &branch, adjacency, aux_rank, hydrogens, groups, centers, flip, first, best, orbits,
        );
    }
}

/// Final canonical rank (1-based) per atom. The mobile groups (rank
/// free: count + member atoms) ride along as the terminal compare's
/// final rows — the reference's t-group CT section, which settles
/// the ties the connection table and the hydrogen layer leave equal.
pub fn canonical_ranks(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
) -> Vec<u32> {
    canonical_ranks_and_symmetry(adjacency, invariants, hydrogens, groups).0
}

/// The connection table over the final numbering.
pub fn canonical_table(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
) -> ConnectionTable {
    let ranks = canonical_ranks(adjacency, invariants, hydrogens, groups);
    ConnectionTable::build(&ranks, adjacency)
}

/// The t-group-aware canonicalization — the reference's TAUT pass
/// (ichican2.c `GetBaseCanonRanking`, the `TAUT_YES` run): the atoms
/// enter pre-ranked by the H-free base pass's equitable partition
/// (`pAtomInvariant[i].val[0] = pRankStack[0][i]` under
/// `bTautGroupsOnly`), the t-groups join as appended nodes with their
/// own invariants (`val[0] = num_at_tg`, `NUM_CONNECTIONS = MAXVAL+1`,
/// the endpoint count and the `num[]` pair at the TAUT_ISO slots), the
/// whole set refines over the augmented graph, and the search
/// individualizes ATOMS only (`PartitionGetFirstCell` over `n`) while
/// the t-group tail rides. The connection table grows the t-group rows
/// (`tgroup_rank, endpoint ranks ascending` — the reference's CT
/// format comment), which is what pulls the group's endpoints against
/// a non-endpoint twin the base numbering left tied (the DoNotAddH
/// bicarbonate witness: the acid pair ranks ahead of the H-less O).
///
/// Returns the atom numbering and the atom equivalence classes of the
/// winning pass — the pair the identifier's `/h` digits and the aux
/// `/N:`/`/E:` speak for a tautomeric structure.
pub fn taut_ranks_and_symmetry(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    hydrogens: &[u32],
    groups: &[(u32, u32, Vec<usize>)],
) -> (Vec<u32>, Vec<u32>) {
    let n = adjacency.len();
    let g = groups.len();
    let total = n + g;

    // The base pass's equitable partition supplies the atom seed ranks.
    let mut base_initial = Partition::initial(invariants);
    base_initial.refine(adjacency);

    // The augmented graph: t-group node n+k links its members.
    let mut augmented: Vec<Vec<usize>> = adjacency.to_vec();
    augmented.resize(total, Vec::new());
    let mut member_of = vec![None; n];
    for (k, (_, _, members)) in groups.iter().enumerate() {
        let tg = n + k;
        for &m in members {
            augmented[m].push(tg);
            augmented[tg].push(m);
            member_of[m] = Some(k);
        }
    }

    // The seed invariants, the reference's two-deck layout packed.
    // atoms: (base rank, 0); t-group nodes: (num_at_tg, MAXVAL+1) and
    // (0, 0, endpoints, num_H, num_minus) behind them.
    let max_hill = invariants.iter().map(|(hill, _)| *hill).max().unwrap_or(0);
    let tg_first: u32 = (total.max(max_hill as usize + 1) as u32) << 8 | 21;
    let mut seed: Vec<(u32, u32)> = (0..n)
        .map(|i| {
            let endpoint = member_of[i];
            match endpoint {
                Some(k) => {
                    let (count, charge, members) = &groups[k];
                    let fields = (*members as u32).min(4095) << 24
                        | (*count as u32).min(4095) << 12
                        | (*charge as u32).min(4095);
                    (base_initial.rank_of(i) as u32, fields)
                }
                None => (base_initial.rank_of(i) as u32, 0),
            }
        })
        .collect();
    for (k, (_, _, members)) in groups.iter().enumerate() {
        let (count, charge) = (groups[k].0, groups[k].1);
        let fields = (*members as u32).min(4095) << 24
            | (count as u32).min(4095) << 12
            | (charge as u32).min(4095);
        seed.push((tg_first, fields));
    }
    let _ = max_hill;

    // The taut pass: bounded stage one over the augmented graph, then
    // the hydrogen-aware stage when a class splits on the non-taut
    // counts.
    let hfree = hfree_canonical_bounded(&augmented, &seed, n);
    let symmetry = orbit_min_symmetry(&hfree.final_search);
    let ranks: Vec<u32>;
    if hydrogens_uniform(&symmetry, hydrogens) {
        ranks = (0..total)
            .map(|node| hfree.final_search.partition.rank_of(node))
            .collect();
    } else {
        let aux_invariants: Vec<(u32, u32)> = (0..total)
            .map(|node| (symmetry[node], encoded_hydrogen(hydrogens.get(node).copied().unwrap_or(0))))
            .collect();
        let mut aux = Partition::initial_pair(&aux_invariants);
        aux.refine(&augmented);
        let aux_rank: Vec<u32> = (0..total).map(|node| aux.rank_of(node)).collect();
        let base = &hfree.h_seed;
        let (informed, _informed_orbits) =
            search_stage_two(&augmented, base, &aux_rank, hydrogens, &[], n);
        ranks = (0..total)
            .map(|node| informed.rank_of(node))
            .collect();
    }

    // The atom prefix is the identifier numbering; the symmetry reads
    // the stage-one orbits restricted to atoms.
    let atom_ranks: Vec<u32> = ranks[..n].to_vec();
    let atom_symmetry: Vec<u32> = symmetry[..n].to_vec();
    (atom_ranks, atom_symmetry)
}

/// Disjoint sets over atoms, merged by position-wise mappings of
/// automorphic terminal labelings.
pub(crate) struct Orbits {
    parent: Vec<usize>,
}

impl Orbits {
    pub(crate) fn new(n: usize) -> Orbits {
        Orbits {
            parent: (0..n).collect(),
        }
    }

    /// The orbit's canonical representative — the SMALLEST member, the
    /// reference's mcr (ichican2.c `nGetMcr2`, "minimal class
    /// representative"). The union keeps the minimum rooted, so the
    /// find is the mcr directly.
    pub(crate) fn representative(&self, atom: usize) -> usize {
        let mut root = atom;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        root
    }

    pub(crate) fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.representative(a), self.representative(b));
        if ra == rb {
            return;
        }
        // Root the union at the smaller member: the representative
        // stays the class minimum, which is what the candidate skip
        // tests (`try v only when v is its own theta mcr`).
        if ra < rb {
            self.parent[rb] = ra;
        } else {
            self.parent[ra] = rb;
        }
    }
}

struct StageOne {
    partition: Partition,
    orbits: Orbits,
}

fn search_stage_one(
    adjacency: &[Vec<usize>],
    root: &Partition,
    trace_tag: &str,
    atom_bound: usize,
) -> StageOne {
    let n = adjacency.len();
    trace_partition(trace_tag, root);
    let mut state = StageOneState {
        first: None,
        best: None,
        orbits: Orbits::new(n),
    };
    descend_hfree(root, adjacency, &mut state, atom_bound);
    StageOne {
        partition: state.best.expect("a discrete partition is always reached"),
        orbits: state.orbits,
    }
}

struct StageOneState {
    first: Option<Partition>,
    best: Option<Partition>,
    orbits: Orbits,
}

/// A fully discrete copy of a leaf partition: any remaining tied cells
/// (the t-group tail the search never branches on) split by node
/// order, so the connection table sees distinct ranks. All leaves
/// resolve identically, so the comparisons keep their meaning.
fn resolved_leaf(
    partition: &Partition,
    adjacency: &[Vec<usize>],
    atom_bound: usize,
) -> Partition {
    if partition.first_nonsingleton_cell_within(atom_bound).is_none()
        && partition.first_nonsingleton_cell().is_none()
    {
        return partition.clone();
    }
    let mut resolved = partition.clone();
    while let Some((start, end)) = resolved.first_nonsingleton_cell() {
        let mut candidates: Vec<usize> = resolved.at_number()[start..end].to_vec();
        candidates.sort_unstable();
        for v in candidates.into_iter().skip(1) {
            resolved.individualize(adjacency, v);
        }
    }
    let _ = atom_bound;
    resolved
}

fn descend_hfree(
    partition: &Partition,
    adjacency: &[Vec<usize>],
    state: &mut StageOneState,
    atom_bound: usize,
) {
    match partition.first_nonsingleton_cell_within(atom_bound) {
        None => {
            let leaf = resolved_leaf(partition, adjacency, atom_bound);
            let table = ConnectionTable::build(leaf.ranks(), adjacency);
            if let Some(first) = &state.first {
                if table.entries() == ConnectionTable::build(first.ranks(), adjacency).entries() {
                    // Automorphic image of the first terminal: merge
                    // the position-wise mapping into the orbits.
                    for position in 0..adjacency.len() {
                        state
                            .orbits
                            .union(first.at_number()[position], leaf.at_number()[position]);
                    }
                }
                if tracing() {
                    eprintln!(
                        "S1 terminal ct={:?} equal_first={}",
                        table.entries(),
                        table.entries()
                            == ConnectionTable::build(first.ranks(), adjacency).entries()
                    );
                }
            } else {
                state.first = Some(leaf.clone());
                if tracing() {
                    eprintln!("S1 terminal ct={:?} FIRST", table.entries());
                }
            }
            let better = match &state.best {
                None => true,
                Some(best_partition) => {
                    table.compare(&ConnectionTable::build(best_partition.ranks(), adjacency))
                        == std::cmp::Ordering::Greater
                }
            };
            if better {
                state.best = Some(leaf);
            }
        }
        Some((start, end)) => {
            if tracing() {
                eprintln!("S1 cell=[{start},{end})");
            }
            // The H-free pass has no auxiliary rank, and the
            // reference's CellGetMinNode then walks the cell by
            // increasing ATOM NUMBER, not by the partition's
            // within-cell order. The order decides which of several
            // table-equal maxima is the winner, and with it the rank
            // values of the symmetry classes.
            let mut candidates: Vec<usize> = partition.at_number()[start..end].to_vec();
            candidates.sort_unstable();
            for v in candidates {
                let mut branch = partition.clone();
                branch.individualize(adjacency, v);
                if tracing() {
                    trace_partition(&format!("S1 indiv v={}", v + 1), &branch);
                }
                descend_hfree(&branch, adjacency, state, atom_bound);
            }
        }
    }
}

fn search_stage_two(
    adjacency: &[Vec<usize>],
    aux: &Partition,
    aux_rank: &[u32],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
    atom_bound: usize,
) -> (Partition, Orbits) {
    let n = adjacency.len();
    let mut best: Option<Partition> = None;
    let mut first: Option<Partition> = None;
    let mut orbits = Orbits::new(n);
    descend_hinformed(
        aux,
        adjacency,
        aux_rank,
        hydrogens,
        groups,
        &mut first,
        &mut best,
        &mut orbits,
        atom_bound,
    );
    (
        best.expect("a discrete partition is always reached"),
        orbits,
    )
}

#[allow(clippy::too_many_arguments)]
fn descend_hinformed(
    partition: &Partition,
    adjacency: &[Vec<usize>],
    aux_rank: &[u32],
    hydrogens: &[u32],
    groups: &[(u32, Vec<usize>)],
    first: &mut Option<Partition>,
    best: &mut Option<Partition>,
    orbits: &mut Orbits,
    atom_bound: usize,
) {
    let cell = match partition.first_nonsingleton_cell_within(atom_bound) {
        Some(span) => span,
        None => {
            let leaf = resolved_leaf(partition, adjacency, atom_bound);
            let ranks = leaf.ranks();
            let table = ConnectionTable::build(ranks, adjacency);
            // Table-equal terminals are automorphic images: fold the
            // position-wise mapping into the orbit partition.
            let matches = |other: &Partition| {
                table.entries() == ConnectionTable::build(other.ranks(), adjacency).entries()
                    && (0..adjacency.len()).all(|position| {
                        hydrogens[other.at_number()[position]]
                            == hydrogens[leaf.at_number()[position]]
                    })
            };
            // The winner compare's full equality — table, hydrogen
            // layer, t-group rows all level — the reference's layered
            // compare reading zero (ichican2.c L8's `nu^gamma == rho`).
            let equals = |other: &Partition| -> bool {
                let other_table = ConnectionTable::build(other.ranks(), adjacency);
                if table.entries() != other_table.entries() {
                    return false;
                }
                let h_rows = |p: &Partition| -> Vec<u32> {
                    (0..adjacency.len())
                        .map(|position| encoded_hydrogen(hydrogens[p.at_number()[position]]))
                        .collect()
                };
                if h_rows(&leaf) != h_rows(other) {
                    return false;
                }
                let rows = |p: &Partition| -> Vec<(u32, Vec<u32>)> {
                    let ranks = p.ranks();
                    groups
                        .iter()
                        .map(|(count, members)| {
                            let mut member_ranks: Vec<u32> =
                                members.iter().map(|&m| ranks[m]).collect();
                            member_ranks.sort_unstable();
                            (*count, member_ranks)
                        })
                        .collect()
                };
                rows(&leaf) == rows(other)
            };
            if let Some(first_partition) = first {
                if matches(first_partition) {
                    if tracing() {
                        eprintln!("S2 merge-FIRST nu_atnum={:?}", leaf.at_number());
                    }
                    for position in 0..adjacency.len() {
                        orbits.union(
                            first_partition.at_number()[position],
                            leaf.at_number()[position],
                        );
                    }
                } else if let Some(best_partition) = best {
                    // The L8 branch: the leaf missed zeta but EQUALS
                    // the best terminal — its transposition against
                    // rho joins theta too (ichican2.c: nu^gamma ==
                    // rho → PartitionGetTransposition(&pi[k-1], &rho)
                    // → L10). Generators the first-terminal star
                    // never witnesses compose here — the reflections
                    // the first-equal unions alone stop short of
                    // (the (1,2,5,6)(3,4) cyclohexadiene witness).
                    if equals(best_partition) {
                        if tracing() {
                            eprintln!("S2 merge-RHO nu_atnum={:?}", leaf.at_number());
                        }
                        for position in 0..adjacency.len() {
                            orbits.union(
                                best_partition.at_number()[position],
                                leaf.at_number()[position],
                            );
                        }
                    }
                }
            } else {
                *first = Some(leaf.clone());
            }
            let counts = |p: &Partition| -> Vec<u32> {
                (0..adjacency.len())
                    .map(|atom| encoded_hydrogen(hydrogens[p.at_number()[atom]]))
                    .collect()
            };
            // The t-group CT rows: (count, ascending member ranks)
            // per group — the compare's final section.
            let rows = |p: &Partition| -> Vec<(u32, Vec<u32>)> {
                let ranks = p.ranks();
                groups
                    .iter()
                    .map(|(count, members)| {
                        let mut member_ranks: Vec<u32> =
                            members.iter().map(|&m| ranks[m]).collect();
                        member_ranks.sort_unstable();
                        (*count, member_ranks)
                    })
                    .collect()
            };
            let better = match &best {
                None => true,
                Some(best_partition) => {
                    let best_table = ConnectionTable::build(best_partition.ranks(), adjacency);
                    match table.compare(&best_table) {
                        std::cmp::Ordering::Greater => true,
                        std::cmp::Ordering::Less => false,
                        std::cmp::Ordering::Equal => {
                            match counts(&leaf).cmp(&counts(best_partition)) {
                                // The M3 corpus decode: at table-equal
                                // terminals the SMALLER encoded
                                // hydrogen vector wins. The encoding
                                // is what places the zero-H atoms
                                // last; on raw counts the winner
                                // would invert (see `encoded_hydrogen`).
                                std::cmp::Ordering::Less => true,
                                std::cmp::Ordering::Greater => false,
                                std::cmp::Ordering::Equal => rows(&leaf) > rows(best_partition),
                            }
                        }
                    }
                }
            };
            if better {
                *best = Some(leaf);
                if tracing() {
                    eprintln!(
                        "S2 terminal BETTER nu_atnum={:?} ct={:?} hydrogens={:?}",
                        partition.at_number(),
                        table.entries(),
                        counts(partition)
                    );
                }
            } else if tracing() {
                eprintln!(
                    "S2 terminal keep nu_atnum={:?} ct={:?} hydrogens={:?}",
                    partition.at_number(),
                    table.entries(),
                    counts(partition)
                );
            }
            return;
        }
    };
    if tracing() {
        eprintln!("S2 cell=[{},{}] arms", cell.0, cell.1);
    }
    // Candidates within the cell ordered by auxiliary rank, then atom
    // number — the reference's aux-informed CellGetMinNode. The theta
    // skip (ichican2.c L14) applies ONLY to the candidates backtracking
    // iterates after the branch's first one: a freshly entered cell's
    // first vertex is `CellGetMinNode(..., 0, ...)` inside the L2
    // descent, chosen before any theta test — the branch it opens is
    // explored even when theta already relates it to a tried one, and
    // the equal-terminal pair it witnesses is a generator the skip
    // would otherwise hide (the breadth-439 witness: the (5,6) branch
    // opened after theta held {5,6} still yields the second gamma that
    // completes the orbit partition).
    let mut candidates: Vec<usize> = partition.at_number()[cell.0..cell.1].to_vec();
    candidates.sort_by_key(|&atom| (aux_rank[atom], atom));
    for (index, &v) in candidates.iter().enumerate() {
        if index > 0 && orbits.representative(v) != v {
            if tracing() {
                eprintln!("S2 skip v={} (orbit)", v + 1);
            }
            continue;
        }
        let mut branch = partition.clone();
        branch.individualize(adjacency, v);
        if tracing() {
            trace_partition(&format!("S2 indiv v={}", v + 1), &branch);
        }
        descend_hinformed(
            &branch,
            adjacency,
            aux_rank,
            hydrogens,
            groups,
            first,
            best,
            orbits,
            atom_bound,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::canonical::invariants::atom_invariants;

    fn run(elements: &[&str], bonds: &[(usize, usize)], hydrogens: &[u32]) -> Vec<u32> {
        let n = elements.len();
        let mut adjacency = vec![Vec::new(); n];
        for &(a, b) in bonds {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        let elements: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
        let invariants = atom_invariants(&elements, &adjacency);
        canonical_ranks(&adjacency, &invariants, hydrogens, &[])
    }

    fn table_of(elements: &[&str], bonds: &[(usize, usize)], hydrogens: &[u32]) -> Vec<u32> {
        let n = elements.len();
        let mut adjacency = vec![Vec::new(); n];
        for &(a, b) in bonds {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        let elements: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
        let invariants = atom_invariants(&elements, &adjacency);
        canonical_table(&adjacency, &invariants, hydrogens, &[])
            .entries()
            .to_vec()
    }

    #[test]
    fn acetaldehyde_numbers_methyl_then_carbonyl_chain() {
        // CH3-CH=O: the oracle records /c 1-2-3 with the methyl H3.
        assert_eq!(
            table_of(&["C", "C", "O"], &[(0, 1), (1, 2)], &[3, 1, 0]),
            vec![1, 2, 1, 3, 2]
        );
    }

    #[test]
    fn benzene_walks_two_then_four() {
        // The oracle records /c 1-2-4-6-5-3-1 around the ring.
        let bonds: Vec<(usize, usize)> = (0..6).map(|i| (i, (i + 1) % 6)).collect();
        assert_eq!(
            table_of(&["C"; 6], &bonds, &[1; 6]),
            vec![1, 2, 1, 3, 1, 4, 2, 5, 3, 6, 4, 5]
        );
    }

    #[test]
    fn isoprene_numbers_the_ch2s_before_the_methyl() {
        // The H counts order the hub's otherwise-equivalent
        // neighbors through the auxiliary ranking: the CH3 at rank
        // 2, the CH at rank 5, the CH2s at ranks 3-4 — the
        // descending-hydrogen layout the corpus records for this
        // skeleton family.
        let ranks = run(
            &["C", "C", "C", "C", "C"],
            &[(0, 1), (1, 2), (3, 4), (1, 3)],
            &[3, 0, 2, 1, 2],
        );
        // The H-auxiliary ranking orders the hub's neighbors: the
        // CH (atom 3, H1) at rank 4, the hub (H0) at rank 5, the
        // CH3/CH2/CH2 at ranks 1-3.
        assert_eq!(ranks, vec![3, 5, 2, 4, 1]);
    }

    #[test]
    fn ranks_are_a_permutation() {
        let bonds: Vec<(usize, usize)> = (0..6).map(|i| (i, (i + 1) % 6)).collect();
        let mut ranks = run(&["C"; 6], &bonds, &[1; 6]);
        ranks.sort_unstable();
        assert_eq!(ranks, vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn zero_hydrogen_atoms_number_after_hydrogenated_ones() {
        // Cyclopropanone carrying an OH on each remaining carbon. The
        // heavy skeleton is symmetric (every ring carbon bears one
        // oxygen), so the H-free search cannot separate the ketone
        // carbon and its oxygen from the CHOH/OH atoms; only the
        // encoded hydrogen layer can. `EMPTY_H_NUMBER` sorts above
        // every real count and the compare minimizes, so the two
        // zero-H atoms land after the four hydrogen-bearing ones. On
        // raw counts the polarity inverts and the ketone atoms take
        // the lowest numbers — the breadth-39 failure shape.
        let ranks = run(
            &["C", "C", "C", "O", "O", "O"],
            &[(0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 5)],
            &[0, 1, 1, 0, 1, 1],
        );
        // The invariant cells separate the carbons from the oxygens
        // (the C block is numbered first), so the encoded-hydrogen
        // preference shows within each block: the ketone carbon after
        // both CHOH carbons, the ketone oxygen after both OH oxygens.
        assert!(ranks[0] > ranks[1] && ranks[0] > ranks[2], "{ranks:?}");
        assert!(ranks[3] > ranks[4] && ranks[3] > ranks[5], "{ranks:?}");
        assert_eq!(ranks, vec![3, 1, 2, 6, 4, 5]);
    }

    #[test]
    fn encoded_hydrogen_puts_the_empty_count_last() {
        assert!(encoded_hydrogen(0) > encoded_hydrogen(3));
        assert!(encoded_hydrogen(3) > encoded_hydrogen(1));
    }
}
