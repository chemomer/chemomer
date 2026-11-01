//! The ±1 proton standardization rules shared by the layers the
//! rules couple: the formula layer adjusts its hydrogens, the
//! proton layer emits the `/p` segment, the charge layer emits
//! `/q`. Decoded empirically from the oracle battery + corpus
//! (2026-10-25): a cation is proton-related iff it carries a
//! removable hydrogen or is a valence-deficient N/P center (an
//! all-heavy onium keeps its charge in `/q`); an anionic O/S/Se/Te/N
//! is proton-related only when it joins a mobile-H group (the parent
//! stores the H in the group — a lone alkoxide keeps `/q-1`); only a
//! BONDED opposite-charge pair annihilates (the terminal ion-pair
//! preprocessing — nitro, N-oxide); a distant zwitterion pair keeps
//! its charges and hides them in the identifier on its own.
//!
//! The gate (2026-10-31) bounds PER-ATOM charges at ±1 — the
//! reference's GetAtomChargeType guard — not the structure total:
//! the multi-charge families (multi-onium `/q+2`/`/q+4`, polyanions
//! `/q-2`/`/q-3`/`/q-4`) ride the declined-(DE)PROTONATION pinned
//! path, whose `ATTOT_TOT_CHARGE` quota arithmetic never assumed a
//! ±1 total.

use crate::canonical::tautomer;
use crate::error::{err, Result};
use crate::graph::{AtomRecord, Structure};

/// The shared gate message: a per-atom charge beyond ±1 is untypable
/// for the M1 layers (the reference's GetAtomChargeType refuses to
/// class the atom — ATT_NONE) and the structure is beyond their
/// authority. The TOTAL charge is deliberately unbounded.
pub const BEYOND_M1: &str = "beyond the ±1 corpus gate";

/// The charge-typability gate — the ONE place the policy lives;
/// formula, proton, charge, isotope and hydrogen layers all call it
/// so none can silently half-decide an untypable structure.
///
/// The guard is the reference's GetAtomChargeType guard
/// (ichi_bns.c): a PER-ATOM charge beyond ±1 is untypable — the C
/// function exits without a class and every charge-exchange quota
/// skips the atom — so no M1 layer may decide such a structure. The
/// TOTAL is deliberately unbounded: the reference's quota arithmetic
/// sums ATTOT_TOT_CHARGE as a plain signed total, and the corpus
/// records the multi-charge families assembling byte-identically
/// through the declined-(DE)PROTONATION pinned path — the double
/// quaternary ammonium `/q+2` (breadth-433), the tetrakis-onium
/// `/q+4` (breadth3-27) and the polyanions `/q-2` `/q-3` `/q-4`
/// (breadth-1072, breadth2-1593, breadth2-1076). The former
/// `|total| > 1` rejection was the M1-era corpus scope, not a rule
/// of the reference: the 2026-10-31 census of the 29 recorded-but-
/// uncounted fixtures found per-atom `|charge|` = 1 on every one.
pub fn ensure_typable_charges(structure: &Structure) -> Result<()> {
    if let Some((index, charge)) = structure
        .atoms
        .iter()
        .enumerate()
        .find_map(|(index, atom)| (atom.charge.abs() > 1).then_some((index, atom.charge)))
    {
        return err(format!(
            "standardization: atom {index} charge {charge:+} {BEYOND_M1}"
        ));
    }
    Ok(())
}

/// The decoded charge plan: how the structure's charges split
/// between the proton balance (`p`) and the fixed-charge identity
/// (`q`), with which side drove the proton shift.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChargePlan {
    pub q: i32,
    pub p: i32,
    pub p_cation: bool,
    pub p_anion: bool,
}

pub(crate) fn cation_is_proton_related(atom: &AtomRecord, degree: u32, bond_sum: f64) -> bool {
    atom.h_count > 0 || proton_deficient_pnicogen(atom, degree, bond_sum)
}

pub fn charge_plan(structure: &Structure) -> ChargePlan {
    let sums = structure.bond_order_sums();
    let degrees = structure.degrees();
    let adjusted = adjusted_hydrogens(structure);
    let mask = tautomer::endpoint_mask(structure, &adjusted);
    let parent_anion = first_charged_heteroatom(structure, false);
    let mut plan = ChargePlan {
        q: 0,
        p: 0,
        p_cation: false,
        p_anion: false,
    };
    for (index, atom) in structure.atoms.iter().enumerate() {
        if atom.charge > 0 {
            if cation_is_proton_related(atom, degrees[index], sums[index]) {
                plan.p += 1;
                plan.p_cation = true;
            } else {
                plan.q += 1;
            }
        } else if atom.charge < 0 {
            let hetero = matches!(atom.element.as_str(), "O" | "S" | "Se" | "Te" | "N");
            let parent_path = parent_anion == Some(index);
            if hetero && (parent_path || mask[index]) {
                plan.p -= 1;
                plan.p_anion = true;
            } else {
                plan.q -= 1;
            }
        }
    }
    plan
}

/// The reference's terminal ion-pair handling (remove_ion_pairs /
/// fix_odd_things preprocessing): an anionic O/S/Se/Te cancels with a
/// FIXED (proton-unrelated) cation DIRECTLY BONDED to it — nitro
/// N⁺–O⁻ and N-oxide (>N⁺–O⁻ ⇒ >N=O class). Distant zwitterion pairs
/// (betaine's quaternary-N / carboxylate, 3+ bonds apart) SURVIVE the
/// whole reference flow — BNS-ENTRY still reads both charges — and the
/// identifier hides them on its own: the net charge is zero (/q
/// omitted) and the charge-only t-group filters out of /h. Canceling
/// a distant pair would erase the −1 the mobile-group machinery must
/// fold into a super-group. A proton-related cation never
/// annihilates: glycine's NH3⁺/COO⁻ pair stays, its charges netting
/// in the proton balance. Runs before any layer, so an annihilated
/// anion registers no riding hydrogen and no mobile group.
pub fn normalize_internal_salts(structure: &mut Structure) {
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); structure.atoms.len()];
    for &(from, to, _) in &structure.bonds {
        neighbors[from].push(to);
        neighbors[to].push(from);
    }
    let sums = structure.bond_order_sums();
    let degrees = structure.degrees();
    let mut consumed = vec![false; structure.atoms.len()];
    let anions: Vec<usize> = structure
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| {
            atom.charge < 0 && matches!(atom.element.as_str(), "O" | "S" | "Se" | "Te" | "N")
        })
        .map(|(index, _)| index)
        .collect();
    for anion in anions {
        if consumed[anion] {
            continue;
        }
        // The directly-bonded fixed cation — the terminal ion pair.
        let found = neighbors[anion].iter().copied().find(|&site| {
            structure.atoms[site].charge > 0
                && !consumed[site]
                && !cation_is_proton_related(&structure.atoms[site], degrees[site], sums[site])
        });
        if let Some(site) = found {
            // The reference never re-derives H after moving a charge:
            // the neutralized anion stays at its anionic reading
            // (zero implicit H), the cation loses its + to the pair,
            // and the pair bond reads the neutral double (the
            // >N⁺–O⁻ ⇒ >N=O promotion — the nitro/N-oxide standard
            // skeleton the tautomer arithmetic reads).
            structure.atoms[anion].charge = 0;
            structure.atoms[anion].h_pinned = true;
            structure.atoms[site].charge = 0;
            for bond in structure.bonds.iter_mut() {
                let (from, to) = (bond.0, bond.1);
                if (from == anion && to == site) || (from == site && to == anion) {
                    bond.2 = 2;
                }
            }
            consumed[anion] = true;
            consumed[site] = true;
        }
    }
}

fn proton_deficient_pnicogen(atom: &AtomRecord, degree: u32, bond_sum: f64) -> bool {
    matches!(atom.element.as_str(), "N" | "P") && degree <= 3 && bond_sum <= 3.0
}

/// The proton-layer value for this structure: `Some("p±n")` when the
/// standardization moved protons, `None` when the charge is part of
/// the identity (quaternary N / lone alkoxide / carbanion). A fully
/// pinned structure never consults the ±1 heuristics: its verdict
/// rode in through the (DE)PROTONATION shift — a pinned structure
/// with a zero shift is the phase's DECLINED decision, and no /p
/// segment may appear.
pub fn proton_layer(structure: &Structure) -> Option<String> {
    if structure.atoms.iter().all(|atom| atom.h_pinned) {
        return None;
    }
    let p = charge_plan(structure).p;
    (p != 0).then(|| format!("p{p:+}"))
}

/// The charge-layer value: the fixed charges the normalization did
/// not move, in total.
pub fn charge_layer(structure: &Structure) -> Option<String> {
    // A fully pinned structure — the (DE)PROTONATION-standardized
    // state — reads its fixed charge directly: the skeleton's total.
    // The /p layer carries the proton shift and net = q + p returns
    // the submitted structure's charge.
    if structure.atoms.iter().all(|atom| atom.h_pinned) {
        let q = structure.total_charge();
        return (q != 0).then(|| format!("q{q:+}"));
    }
    let q = charge_plan(structure).q;
    (q != 0).then(|| format!("q{q:+}"))
}

/// The hydrogen adjustment the formula layer applies: an anion gains
/// its mobile H, a proton-carrying cation loses one.
pub fn hydrogen_adjustment(structure: &Structure) -> i32 {
    // Pinned structures carry their final hydrogen counts; the
    // legacy ±1 heuristics would double-apply.
    if structure.atoms.iter().all(|atom| atom.h_pinned) {
        return 0;
    }
    let plan = charge_plan(structure);
    match (plan.p_cation, plan.p_anion) {
        (true, false) => -1,
        (false, true) => 1,
        _ => 0,
    }
}

/// The hydrogen count per atom with the ±1 standardization applied at
/// its site. Anionic O/N split two ways: carried by carbon they take
/// the salt path — a riding hydrogen that survives only when the
/// mobile-H machinery finds the anion a group (the parent stores the
/// H in the group; a lone alkoxide keeps none) — while carried by
/// anything else (hydroxide, nitrate) they take the parent path: the
/// riding hydrogen stays FIXED, never mobile. A proton-related cation
/// loses one. The `/h` layer reads these; the formula layer's total
/// is their sum.
pub fn adjusted_hydrogens(structure: &Structure) -> Vec<u32> {
    // A fully pinned structure already carries the standardized
    // (DE)PROTONATION state — no further adjustment.
    if structure.atoms.iter().all(|atom| atom.h_pinned) {
        return structure.hydrogens();
    }
    let mut counts = structure.hydrogens();
    let parent_anion = first_charged_heteroatom(structure, false);
    let salt_anion = first_charged_heteroatom(structure, true);
    if let Some(index) = parent_anion.or(salt_anion) {
        counts[index] += 1;
    }
    if let Some(index) = first_cation_site(structure) {
        counts[index] = counts[index].saturating_sub(1);
    }
    if let Some(index) = salt_anion {
        let mask = tautomer::endpoint_mask(structure, &counts);
        if !mask[index] {
            counts[index] -= 1;
        }
    }
    counts
}

/// The first anionic O/N, optionally restricted to those carried by
/// carbon (the salt candidates; the reference's GetSaltChargeType
/// anchor condition).
fn first_charged_heteroatom(structure: &Structure, carbon_carried: bool) -> Option<usize> {
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); structure.atoms.len()];
    for &(from, to, _) in &structure.bonds {
        neighbors[from].push(to);
        neighbors[to].push(from);
    }
    structure
        .atoms
        .iter()
        .enumerate()
        .position(|(index, atom)| {
            if atom.charge >= 0 || !(atom.element == "O" || atom.element == "N") {
                return false;
            }
            let carried = neighbors[index]
                .iter()
                .any(|&n| structure.atoms[n].element == "C");
            carried == carbon_carried
        })
}

fn first_cation_site(structure: &Structure) -> Option<usize> {
    let sums = structure.bond_order_sums();
    let degrees = structure.degrees();
    structure
        .atoms
        .iter()
        .enumerate()
        .position(|(index, atom)| {
            atom.charge > 0
                && (atom.h_count > 0
                    || proton_deficient_pnicogen(atom, degrees[index], sums[index]))
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;

    fn atom(element: &str, charge: i32, h_count: u32) -> AtomRecord {
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
            stereo_bonds: Default::default(),
            positions: Default::default(),
            isotopic_h: Default::default(),
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            fold_source: Vec::new(),
            chiral_flag: false,
        }
    }

    #[test]
    fn protonated_amines_are_proton_related() {
        // Oracle battery: [NH4+] and C[NH+](C)C carry /p+1; the /h
        // parent loses the riding H.
        let s = structure(vec![atom("N", 1, 4)], vec![]);
        assert_eq!(proton_layer(&s), Some("p+1".to_string()));
        assert_eq!(charge_layer(&s), None);
        assert_eq!(adjusted_hydrogens(&s), vec![3]);
    }

    #[test]
    fn quaternary_nitrogen_is_fixed() {
        // C[N+](C)(C)C: /q+1, no /p, no hydrogen moved.
        let s = structure(
            vec![
                atom("N", 1, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
            ],
            vec![(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
        );
        assert_eq!(proton_layer(&s), None);
        assert_eq!(charge_layer(&s), Some("q+1".to_string()));
    }

    #[test]
    fn sp3_pnicogen_cations_without_hydrogen_are_proton_related() {
        // Quinuclidinium-class N+ (three single bonds, no H): the
        // charge is the proton deficit (fixture identity/breadth-233).
        let s = structure(
            vec![
                atom("N", 1, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
            ],
            vec![(0, 1, 1), (0, 2, 1), (0, 3, 1)],
        );
        assert_eq!(proton_layer(&s), Some("p+1".to_string()));
        assert_eq!(charge_layer(&s), None);
    }

    #[test]
    fn lone_alkoxide_is_fixed_and_keeps_no_riding_hydrogen() {
        // [O-]C: /q-1; the formula stays CH3O — the riding H exists
        // only inside a mobile group.
        let s = structure(vec![atom("O", -1, 0), atom("C", 0, 3)], vec![(0, 1, 1)]);
        assert_eq!(proton_layer(&s), None);
        assert_eq!(charge_layer(&s), Some("q-1".to_string()));
        assert_eq!(hydrogen_adjustment(&s), 0);
        assert_eq!(adjusted_hydrogens(&s), vec![0, 3]);
    }

    #[test]
    fn carboxylate_anion_joins_the_group_and_shifts_the_proton() {
        // CC(=O)[O-]: /p-1, formula gains the parent H, the anionic O
        // keeps its riding H (the mobile group holds it).
        let s = structure(
            vec![
                atom("C", 0, 3),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", -1, 0),
            ],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        assert_eq!(proton_layer(&s), Some("p-1".to_string()));
        assert_eq!(charge_layer(&s), None);
        assert_eq!(hydrogen_adjustment(&s), 1);
    }

    #[test]
    fn zwitterionic_glycine_nets_to_zero() {
        // [NH3+]CC(=O)[O-]: p = +1 (N) - 1 (carboxylate) = 0 — no
        // charge layers, identical to neutral glycine.
        let s = structure(
            vec![
                atom("N", 1, 3),
                atom("C", 0, 2),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", -1, 0),
            ],
            vec![(0, 1, 1), (1, 2, 1), (2, 3, 2), (2, 4, 1)],
        );
        assert_eq!(proton_layer(&s), None);
        assert_eq!(charge_layer(&s), None);
    }

    #[test]
    fn distant_zwitterion_charges_survive() {
        // C[N+](C)(C)CC(=O)[O-]: the reference keeps the pair's
        // charges (BNS-ENTRY reads both — the oracle trace); only a
        // BONDED ion pair cancels. The identifier hides them on its
        // own: net zero (no /q) and the charge-only t-group filtered
        // from /h once the structure pins.
        let mut s = structure(
            vec![
                atom("N", 1, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 2),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", -1, 0),
            ],
            vec![
                (0, 1, 1),
                (0, 2, 1),
                (0, 3, 1),
                (0, 4, 1),
                (4, 5, 1),
                (5, 6, 2),
                (5, 7, 1),
            ],
        );
        normalize_internal_salts(&mut s);
        assert_eq!(s.atoms[0].charge, 1, "the distant N+ survives");
        assert_eq!(s.atoms[7].charge, -1, "the distant carboxylate survives");
    }

    #[test]
    fn protonated_cation_never_annihilates() {
        // Glycine's N is proton-related: the pair stays and nets in
        // the proton balance instead.
        let mut s = structure(
            vec![
                atom("N", 1, 3),
                atom("C", 0, 2),
                atom("C", 0, 0),
                atom("O", 0, 0),
                atom("O", -1, 0),
            ],
            vec![(0, 1, 1), (1, 2, 1), (2, 3, 2), (2, 4, 1)],
        );
        normalize_internal_salts(&mut s);
        assert_eq!(s.atoms[0].charge, 1);
        assert_eq!(s.atoms[4].charge, -1);
    }

    #[test]
    fn multi_onium_totals_assemble_from_the_pinned_state() {
        // The (DE)PROTONATION declined on the double quaternary N+
        // (no removable proton anywhere): the pinned total emits
        // /q+2 with no /p (the identity/breadth-433 class — the
        // 2026-10-31 census). The total is unbounded; the gate is
        // per-atom.
        let mut s = structure(
            vec![
                atom("N", 1, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("N", 1, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
            ],
            vec![
                (0, 1, 1),
                (0, 2, 1),
                (0, 3, 1),
                (4, 5, 1),
                (4, 6, 1),
                (4, 7, 1),
            ],
        );
        for record in s.atoms.iter_mut() {
            record.h_pinned = true;
        }
        assert_eq!(charge_layer(&s), Some("q+2".to_string()));
        assert_eq!(proton_layer(&s), None);
        assert_eq!(hydrogen_adjustment(&s), 0);
    }

    #[test]
    fn per_atom_charges_beyond_plus_minus_one_stay_beyond() {
        // GetAtomChargeType's guard: a per-atom |charge| > 1 is
        // untypable (ATT_NONE); the polyanion totals — every atom
        // within ±1 — stay inside the gate (breadth-1072 emits
        // /q-2 through the pinned path).
        let pair = structure(vec![atom("O", -1, 0), atom("O", -1, 0)], vec![(0, 1, 1)]);
        assert!(ensure_typable_charges(&pair).is_ok());
        let dianion = structure(vec![atom("O", -2, 0), atom("C", 0, 3)], vec![(0, 1, 1)]);
        assert!(ensure_typable_charges(&dianion).is_err());
    }
}
