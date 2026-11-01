//! The Fixed-H layer's canonical machinery — the `/f` block of the
//! full (non-standard) InChI, ported from the reference's
//! GetBaseCanonRankings fixed-H section (ichican2.c 6246-6400) and its
//! writers (ichiprt1.c `OutputINCHI_FixedHLayerWithSublayers`,
//! ichimake.c `CompINChILayers`). Dev-time source read, never
//! compiled or linked; the corpus arbitrates.
//!
//! The fixed-H structure is the PREPROCESSED structure copied before
//! the taut pass (ichimake.c:3904: the TAUT_NON memcpy sits before the
//! TAUT_YES `mark_alt_bonds_and_taut_groups`) — the (DE)PROTONATION's
//! moved protons never touch it. That snapshot then runs the SAME
//! perception crippled (t_group_info = NULL, ichimake.c:4090):
//! MarkChargeGroups, BnsTestAndMarkAltBonds and the radical
//! adjustment still run (ichi_bns.c:5423-5487), while the
//! (DE)PROTONATION and MarkTautomerGroups return early
//! (ichitaut.c:4413) — a normalization whose observable surface (atom
//! H counts, the connection table, charge totals) is invariant under
//! positive-charge migration and alternation marking, so the /f
//! outputs read straight off the snapshot.
//!
//! Fixed-H H counts are read off during canonicalization
//! (ichican2.c:6256-6272): a taut endpoint takes its SNAPSHOT count
//! (BASE_H_NUMBER/EMPTY_H_NUMBER encoding), a non-endpoint that
//! differs between the two structures takes the signed delta,
//! everything else is EMPTY (zero after CleanNumH). The fixed-H
//! numbering (CanonGraph07) grows the same equitable tree as the
//! notaut pass; the auxiliary partition refined from
//! (nSymmRankBase[i], NumHfixed[i]) (ichican2.c:6310-6325) only
//! ORDERS the candidates within cells, and the terminal compare
//! extends the notaut rows with the encoded deltas. The reference's
//! consistency check (ichican2.c:6625-6636) requires the resulting
//! connection table to EQUAL the notaut one and asserts it here
//! rather than skipping.

use crate::error::{err, Result};
use crate::graph::Structure;

use super::partition::Partition;
use super::search::{encoded_hydrogen, Orbits, BASE_H_NUMBER, EMPTY_H_NUMBER};
use super::table::ConnectionTable;

/// The fixed-H analysis for one structure: the canonical numbering and
/// symmetry of the fixed-H representation, its per-atom H deltas
/// (CleanNumH'd `NumHfixed`), its Hill formula and total charge.
#[derive(Clone)]
pub struct FixedH {
    /// The fixed-H canonical rank (1-based) per atom — the numbering
    /// the /f /h deltas and the aux /F: block speak.
    pub ranks: Vec<u32>,
    /// The fixed-H orbit-minimal symmetry class per atom (the aux
    /// /F: equivalence input).
    pub symmetry: Vec<u32>,
    /// Per-atom fixed-H delta in atom order: a taut endpoint's
    /// absolute snapshot count, a differing non-endpoint's signed
    /// delta, zero elsewhere (CleanNumH of NumHfixed).
    pub deltas: Vec<i32>,
    /// The fixed-H representation's Hill formula (the snapshot's own
    /// composition; the /f formula segment).
    pub formula: String,
    /// The fixed-H representation's total charge (the snapshot's; the
    /// /f /q segment's value).
    pub total_charge: i32,
    /// Whether the aux /F: block carries the `it:m` marker — the
    /// corpus speaks it wherever the main block's verdict adopted its
    /// inverted pass (`it:im`).
    pub inverted: bool,
}

/// The fixed-H deltas (ichican2.c:6256-6272) between the working
/// structure (`base`, the mobile-H representation the standard layers
/// speak — pinned counts when the (DE)PROTONATION fired) and the
/// pre-taut snapshot (`other`).
pub fn fixed_h_deltas(base: &[u32], other: &[u32], endpoints: &[bool]) -> Vec<i32> {
    base.iter()
        .zip(other)
        .zip(endpoints)
        .map(|((&b, &o), &endpoint)| {
            if endpoint {
                o as i32
            } else if o != b {
                o as i32 - b as i32
            } else {
                0
            }
        })
        .collect()
}

/// The fixed-H analysis over the mobile-H representation (`working`)
/// and its pre-taut `snapshot`. `centers` are the wedge-defined stereo
/// centers (the main pass's inputs); `main_inverted` says the main aux
/// verdict adopted its inverted pass — the /F: `it:` marker mirrors it.
pub fn analyze(
    snapshot: &Structure,
    working: &Structure,
    centers: &[(usize, Vec<usize>, u32)],
    main_inverted: bool,
) -> Result<FixedH> {
    let (adjacency, invariants) = crate::canonical::graph_of(working);
    let counts = crate::layers::standardization::adjusted_hydrogens(working);
    let mask = crate::canonical::tautomer::endpoint_mask(working, &counts);
    let deltas = fixed_h_deltas(&working.hydrogens(), &snapshot.hydrogens(), &mask);
    let total_charge = snapshot.total_charge();
    let formula = hill_formula(snapshot);
    let inverted = !centers.is_empty() && main_inverted;

    // The fixed-H canonicalization proper runs when the mobile-H
    // representation carries REGISTERED t-groups (nLenLinearCTTautomer
    // > 0); otherwise the TAUT_NON results ARE the notaut ones
    // (ichican2.c:6731-6754 rearrangement). A pinned charge-only
    // group — the (DE)PROTONATION-declined zwitterion's carboxylate —
    // registers too (its mobile unit is the charge), so its delta rows
    // still reorder the numbering (breadth2-429).
    let groups_rows = crate::canonical::group_membership(working, &counts);
    let tautomeric = !groups_rows.is_empty();

    let (ranks, symmetry) = if !tautomeric {
        // Non-tautomeric: the fixH results are the notaut results —
        // the numbering the aux speaks on this path.
        let (ranks, symmetry) = crate::canonical::canonical_ranks_and_symmetry(
            &adjacency,
            &invariants,
            &counts,
            &groups_rows,
        );
        (ranks, symmetry)
    } else {
        let notaut = crate::canonical::tautomer::without_mobile_h(&counts, &mask);
        // nSymmRankBase: the base (mobile-H) canonicalization's orbit
        // minima — the corpus-validated /c search supplies them.
        let (_base_ranks, base_symmetry) = crate::canonical::canonical_ranks_and_symmetry(
            &adjacency,
            &invariants,
            &notaut,
            &groups_rows,
        );
        // Ct_NoTautH: the notaut numbering the consistency check
        // compares against (CanonGraph03, no t-group rows).
        let notaut_ranks =
            crate::canonical::canonical_ranks_and_symmetry(&adjacency, &invariants, &notaut, &[]).0;
        if centers.is_empty() {
            let state = fixed_h_search(
                &adjacency,
                &invariants,
                &base_symmetry,
                &notaut,
                &deltas,
                &[],
                false,
            );
            let ranks: Vec<u32> = state
                .best
                .expect("a discrete partition is always reached")
                .ranks()
                .to_vec();
            assert_fixed_h_consistency(&adjacency, &ranks, &notaut_ranks)?;
            let symmetry = orbit_min_symmetry(&ranks, &state.orbits);
            (ranks, symmetry)
        } else {
            // The stereo-informed fixed-H numbering (the TAUT_NON
            // pass's nCanonOrdStereo): the seeded search with the
            // stereo rows appended, both parities, winner adopted.
            let (ranks, symmetry, inverted_pass) = fixed_h_search_stereo(
                &adjacency,
                &invariants,
                &base_symmetry,
                &notaut,
                &deltas,
                centers,
            );
            assert_fixed_h_consistency(&adjacency, &ranks, &notaut_ranks)?;
            let _ = inverted_pass;
            (ranks, symmetry)
        }
    };

    Ok(FixedH {
        ranks,
        symmetry,
        deltas,
        formula,
        total_charge,
        inverted,
    })
}

/// The Hill formula of a structure's own composition — the fixed-H
/// representation's formula reads the SNAPSHOT's atoms and counts
/// (the C's szHillFormula of the TAUT_NON INChI).
fn hill_formula(structure: &Structure) -> String {
    let composition = structure.composition();
    let has_carbon = composition.iter().any(|(e, _)| e == "C");
    let rank = |element: &str| match (has_carbon, element) {
        (true, "C") => 0,
        (true, "H") => 1,
        (true, _) => 2,
        (false, _) => 0,
    };
    let mut entries: Vec<&(String, u32)> = composition.iter().collect();
    entries.sort_by(|a, b| rank(&a.0).cmp(&rank(&b.0)).then(a.0.cmp(&b.0)));
    entries
        .iter()
        .filter(|(_, count)| *count > 0)
        .map(|(element, count)| {
            if *count > 1 {
                format!("{element}{count}")
            } else {
                element.clone()
            }
        })
        .collect()
}

/// The encoded delta for the compare rows and the seed: zero (EMPTY)
/// sorts after every real count, a negative delta below every
/// positive one — mirroring the reference's `NumHfixed` encoding
/// (`value + BASE_H_NUMBER`, `EMPTY_H_NUMBER` for nothing).
fn encoded_delta(delta: i32) -> u32 {
    if delta == 0 {
        EMPTY_H_NUMBER
    } else {
        (i64::from(BASE_H_NUMBER) + i64::from(delta)) as u32
    }
}

/// Orbit-minimal symmetry over a search's terminal ranks and orbits.
fn orbit_min_symmetry(ranks: &[u32], orbits: &Orbits) -> Vec<u32> {
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

/// The seeded fixed-H search state.
struct FixedHSearchState {
    first: Option<Partition>,
    best: Option<Partition>,
    orbits: Orbits,
}

/// The seeded fixed-H search: the tree grows from the EQUITABLE
/// skeleton partition (the reference's pi[0] — pRankStack[0] is never
/// written back by an earlier pass, ichican2.c:5396-5417), candidates
/// order by the (nSymmRankBase, encoded NumHfixed) auxiliary rank
/// refined equitably (ichican2.c:6310-6325 — the aux partition only
/// orders candidates, exactly like CanonGraph03/05's nRankAux), and
/// the terminal compare is (connection table, encoded notaut H,
/// encoded deltas, stereo rows) with no t-group rows — the
/// reference's CanonGraph07 over pCD[iOther].
#[allow(clippy::too_many_arguments)]
fn fixed_h_search(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    base_symmetry: &[u32],
    notaut: &[u32],
    deltas: &[i32],
    centers: &[(usize, Vec<usize>, u32)],
    flip: bool,
) -> FixedHSearchState {
    let n = adjacency.len();
    let mut root = Partition::initial(invariants);
    root.refine(adjacency);
    let aux_invariants: Vec<(u32, u32)> = (0..n)
        .map(|atom| (base_symmetry[atom], encoded_delta(deltas[atom])))
        .collect();
    let mut aux = Partition::initial_pair(&aux_invariants);
    aux.refine(adjacency);
    let aux_rank: Vec<u32> = aux.ranks().to_vec();

    let mut state = FixedHSearchState {
        first: None,
        best: None,
        orbits: Orbits::new(n),
    };
    descend(
        &root, adjacency, &aux_rank, notaut, deltas, centers, flip, &mut state,
    );
    state
}

#[allow(clippy::too_many_arguments)]
fn descend(
    partition: &Partition,
    adjacency: &[Vec<usize>],
    aux_rank: &[u32],
    notaut: &[u32],
    deltas: &[i32],
    centers: &[(usize, Vec<usize>, u32)],
    flip: bool,
    state: &mut FixedHSearchState,
) {
    let cell = match partition.first_nonsingleton_cell() {
        Some(span) => span,
        None => {
            let ranks = partition.ranks();
            let table = ConnectionTable::build(ranks, adjacency);
            // Table+H+delta-equal terminals are automorphic images:
            // fold the position-wise mapping into the orbits.
            let matches = |other: &Partition| {
                table.entries() == ConnectionTable::build(other.ranks(), adjacency).entries()
                    && (0..adjacency.len()).all(|position| {
                        notaut[other.at_number()[position]]
                            == notaut[partition.at_number()[position]]
                            && deltas[other.at_number()[position]]
                                == deltas[partition.at_number()[position]]
                    })
            };
            if let Some(first) = &state.first {
                if matches(first) {
                    for position in 0..adjacency.len() {
                        state
                            .orbits
                            .union(first.at_number()[position], partition.at_number()[position]);
                    }
                }
            } else {
                state.first = Some(partition.clone());
            }
            let h_rows = |p: &Partition| -> Vec<u32> {
                (0..adjacency.len())
                    .map(|position| encoded_hydrogen(notaut[p.at_number()[position]]))
                    .collect()
            };
            let delta_rows = |p: &Partition| -> Vec<u32> {
                (0..adjacency.len())
                    .map(|position| encoded_delta(deltas[p.at_number()[position]]))
                    .collect()
            };
            let stereo_rows = |p: &Partition| -> Vec<u32> {
                super::stereo::terminal_rows(centers, p.ranks(), flip)
            };
            let better = match &state.best {
                None => true,
                Some(best_partition) => {
                    let best_table = ConnectionTable::build(best_partition.ranks(), adjacency);
                    match table.compare(&best_table) {
                        std::cmp::Ordering::Greater => true,
                        std::cmp::Ordering::Less => false,
                        std::cmp::Ordering::Equal => {
                            match h_rows(partition).cmp(&h_rows(best_partition)) {
                                std::cmp::Ordering::Less => true,
                                std::cmp::Ordering::Greater => false,
                                std::cmp::Ordering::Equal => {
                                    match delta_rows(partition).cmp(&delta_rows(best_partition)) {
                                        // The smaller encoded delta vector
                                        // wins — the H-holding endpoint
                                        // numbers first (the reference's
                                        // INCHI_CANON_MIN over NumHfixed).
                                        std::cmp::Ordering::Less => true,
                                        std::cmp::Ordering::Greater => false,
                                        std::cmp::Ordering::Equal => {
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
                state.best = Some(partition.clone());
            }
            return;
        }
    };
    let mut candidates: Vec<usize> = partition.at_number()[cell.0..cell.1].to_vec();
    candidates.sort_by_key(|&atom| (aux_rank[atom], atom));
    // Every candidate explores: the reference prunes with its
    // level-scoped stabilizer chain (Omega/Phi intersected per depth,
    // ichican2.c L14/L17), which a flat orbit partition cannot
    // reproduce — a merge discovered under one individualization
    // prefix does not license skipping under another. The orbit
    // merges stay for the symmetry classes; only the pruning goes.
    for &v in &candidates {
        let mut branch = partition.clone();
        branch.individualize(adjacency, v);
        descend(
            &branch, adjacency, aux_rank, notaut, deltas, centers, flip, state,
        );
    }
}

/// The stereo-informed fixed-H numbering: the seeded search run over
/// both parities, the winner adopted by the reference's diff (the
/// `canonical_ranks_and_symmetry_stereo` pattern over the fixed-H
/// compare rows). Returns the winner's ranks and symmetry plus
/// whether the inverted pass won.
fn fixed_h_search_stereo(
    adjacency: &[Vec<usize>],
    invariants: &[(u16, u16)],
    base_symmetry: &[u32],
    notaut: &[u32],
    deltas: &[i32],
    centers: &[(usize, Vec<usize>, u32)],
) -> (Vec<u32>, Vec<u32>, bool) {
    let run = |flip: bool| {
        let state = fixed_h_search(
            adjacency,
            invariants,
            base_symmetry,
            notaut,
            deltas,
            centers,
            flip,
        );
        let ranks: Vec<u32> = state
            .best
            .expect("a discrete partition is always reached")
            .ranks()
            .to_vec();
        let symmetry = orbit_min_symmetry(&ranks, &state.orbits);
        let pairs = super::stereo::terminal_pairs(centers, &ranks, flip);
        (ranks, symmetry, pairs)
    };
    let (abs_ranks, abs_symmetry, abs_pairs) = run(false);
    let (inv_ranks, inv_symmetry, inv_pairs) = run(true);
    // The reference's diff: centers' canonical numbers first, then
    // parities — the first difference carries the sign.
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
        (inv_ranks, inv_symmetry, true)
    } else {
        (abs_ranks, abs_symmetry, false)
    }
}

/// The reference's consistency check (ichican2.c:6625-6636): the
/// fixed-H connection table must EQUAL the notaut one. A violation is
/// the reference's CT_CANON_ERR — the identifier cannot be produced.
fn assert_fixed_h_consistency(
    adjacency: &[Vec<usize>],
    fixh_ranks: &[u32],
    notaut_ranks: &[u32],
) -> Result<()> {
    let fixh = ConnectionTable::build(fixh_ranks, adjacency);
    let notaut = ConnectionTable::build(notaut_ranks, adjacency);
    if fixh.entries() != notaut.entries() {
        return err(format!(
            "fixedh: Ct_FixH != Ct_NoTautH (CT_CANON_ERR): {:?} vs {:?}",
            fixh.entries(),
            notaut.entries()
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;

    fn atom(element: &str, h_count: u32, charge: i32) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count,
            charge,
            h_pinned: false,
            isotope: 0,
        }
    }

    fn structure(atoms: Vec<AtomRecord>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure {
            atoms,
            bonds,
            ..Default::default()
        }
    }

    fn pin_all(structure: &mut Structure) {
        for record in structure.atoms.iter_mut() {
            record.h_pinned = true;
        }
    }

    fn acetate_snapshot() -> Structure {
        structure(
            vec![
                atom("C", 3, 0),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", 0, -1),
            ],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        )
    }

    fn acetate_working() -> Structure {
        // The (DE)PROTONATION-pinned mobile rep: the acid H rides the
        // group, the charge left with the proton.
        let mut s = acetate_snapshot();
        pin_all(&mut s);
        s.atoms[3].h_count = 1;
        s.atoms[3].charge = 0;
        s
    }

    #[test]
    fn acetate_deltas_read_endpoints_absolute() {
        let snapshot = acetate_snapshot();
        let working = acetate_working();
        let counts = crate::layers::standardization::adjusted_hydrogens(&working);
        let mask = crate::canonical::tautomer::endpoint_mask(&working, &counts);
        let deltas = fixed_h_deltas(&working.hydrogens(), &snapshot.hydrogens(), &mask);
        // Both oxygens are taut endpoints of the (H,3,4) group; the
        // snapshot holds no H on either -> both EMPTY -> zero.
        assert_eq!(deltas, vec![0, 0, 0, 0]);
    }

    #[test]
    fn acid_holding_oxygen_takes_its_absolute_count() {
        // CH3-C(=O)-OH: the drawn OH oxygen (atom 4) is an endpoint
        // and keeps its snapshot H -> delta +1; the carbonyl endpoint
        // reads EMPTY -> zero.
        let snapshot = structure(
            vec![
                atom("C", 3, 0),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", 1, 0),
            ],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let working = snapshot.clone();
        let counts = crate::layers::standardization::adjusted_hydrogens(&working);
        let mask = crate::canonical::tautomer::endpoint_mask(&working, &counts);
        let deltas = fixed_h_deltas(&working.hydrogens(), &snapshot.hydrogens(), &mask);
        assert_eq!(deltas, vec![0, 0, 0, 1]);
    }

    #[test]
    fn ammonium_delta_is_the_proton_the_standardization_removed() {
        // [NH4+]: at_base (pinned) 3 H, the snapshot 4 -> non-endpoint
        // delta +1.
        let snapshot = structure(vec![atom("N", 4, 1)], vec![]);
        let mut working = structure(vec![atom("N", 3, 0)], vec![]);
        pin_all(&mut working);
        let counts = crate::layers::standardization::adjusted_hydrogens(&working);
        let mask = crate::canonical::tautomer::endpoint_mask(&working, &counts);
        let deltas = fixed_h_deltas(&working.hydrogens(), &snapshot.hydrogens(), &mask);
        assert_eq!(deltas, vec![1]);
    }

    #[test]
    fn hydroxide_delta_is_negative() {
        // [OH-] -> the parent-path mobile rep H2O: the O keeps its H
        // fixed (not an endpoint) and loses one -> delta -1.
        let snapshot = structure(vec![atom("O", 1, -1)], vec![]);
        let mut working = structure(vec![atom("O", 2, 0)], vec![]);
        pin_all(&mut working);
        let counts = crate::layers::standardization::adjusted_hydrogens(&working);
        let mask = crate::canonical::tautomer::endpoint_mask(&working, &counts);
        let deltas = fixed_h_deltas(&working.hydrogens(), &snapshot.hydrogens(), &mask);
        assert_eq!(deltas, vec![-1]);
    }

    #[test]
    fn encoded_delta_orders_negative_then_positive_then_empty() {
        assert!(encoded_delta(-1) < encoded_delta(1));
        assert!(encoded_delta(1) < encoded_delta(2));
        assert!(encoded_delta(2) < encoded_delta(0));
    }

    #[test]
    fn acetate_analysis_ranks_and_charges() {
        let snapshot = acetate_snapshot();
        let working = acetate_working();
        let fixed = analyze(&snapshot, &working, &[], false).expect("analyzes");
        // The /f formula reads the snapshot: CH3-C(=O)-O(-) with no
        // riding H -> C2H3O2 (one H short of the mobile rep).
        assert_eq!(fixed.formula, "C2H3O2");
        assert_eq!(fixed.total_charge, -1);
        // No nonzero deltas -> the /f /h layer stays empty.
        assert!(fixed.deltas.iter().all(|&d| d == 0));
    }

    #[test]
    fn acid_analysis_puts_the_holding_oxygen_first() {
        // The drawn carboxylic acid: the fixH numbering ranks the
        // H-holding oxygen before the carbonyl one (aspirin's F
        // numbering swaps 12/13; glycine's swaps 4/5).
        let snapshot = structure(
            vec![
                atom("C", 3, 0),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", 1, 0),
            ],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let working = snapshot.clone();
        let fixed = analyze(&snapshot, &working, &[], false).expect("analyzes");
        // The H-holder (atom 3, the OH) ranks before the carbonyl
        // (atom 2).
        assert!(fixed.ranks[2] > fixed.ranks[3], "{:?}", fixed.ranks);
        assert_eq!(fixed.deltas, vec![0, 0, 0, 1]);
    }
}
