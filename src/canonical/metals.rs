//! The metal element data — the reference's ElData metal columns
//! (util.c:103-277): the METAL/METAL2 `nType` bits, the
//! `bSkipAddingH` read-time hydrogen refusal, and the per-charge
//! valence rows the salt/coordination classification reads. The
//! valence rows themselves live in [`super::protons::el_valences`]
//! (the one cValence table); this module holds the flags and the
//! classification helpers built on them, and — unit #107/#113 — the
//! standard-path disconnection transforms.
//!
//! Bit semantics (incomdef.h:115-117): `METAL` (1) marks an element
//! whose LOWEST registered valence is the normal one; `METAL2` (3)
//! marks one whose lowest AND next valences are both normal. The
//! salt check compares a metal's bond count against the valence row
//! selected by those bits.

use super::protons::el_valences;
use crate::graph::Structure;

/// The reference's `METAL` bit (incomdef.h:115).
pub(crate) const METAL: u32 = 1;
/// The reference's `METAL2` bit — METAL plus "the next valence is
/// normal too" (incomdef.h:116); `IS_METAL` is the 3-bit mask.
pub(crate) const METAL2: u32 = 3;

/// The ElData `nType` column for the element's periodic number: 0 for
/// non-metals, [`METAL`] or [`METAL2`] for metals (util.c ElData).
/// The transuranics reproduce the reference's placeholder rows (all
/// plain `METAL` with a single valence of 1).
pub(crate) fn el_type(el: i32) -> u32 {
    match el {
        3 | 4 | 11 | 12 | 13 | 19 | 20 | 21 | 22 | 23 | 24 | 29 | 30 | 31 | 37 | 38 | 39 | 40
        | 41 | 42 | 43 | 44 | 45 | 46 | 47 | 48 | 49 | 51 | 55 | 56 | 57 | 60 | 61 | 64 | 66
        | 67 | 68 | 71 | 72 | 73 | 79 | 83 | 87 | 88 | 89 | 96 | 97 | 98 | 99 | 100 | 101 | 102
        | 103 | 104 | 105 | 106 | 107 | 108 | 109 | 111 | 112 | 113 | 114 | 115 | 116 | 117
        | 118 => METAL,
        25 | 26 | 27 | 28 | 50 | 58 | 59 | 62 | 63 | 65 | 69 | 70 | 74 | 75 | 76 | 77 | 78 | 80
        | 81 | 82 | 84 | 90 | 91 | 92 | 93 | 94 | 95 | 110 => METAL2,
        _ => 0,
    }
}

/// The `IS_METAL` test (util.c:679 `get_el_type`).
pub(crate) fn is_metal_el(el: i32) -> bool {
    el_type(el) & METAL2 != 0
}

/// The ElData `bSkipAddingH` column: elements the reader never adds
/// implicit hydrogen to (util.c:912 — the read-time inference is
/// skipped when this is set, independently of DoNotAddH). The
/// transuranic rows carry the reference's placeholder 1.
pub(crate) fn skip_adding_h(el: i32) -> bool {
    matches!(
        el,
        21 | 22
            | 23
            | 24
            | 25
            | 26
            | 27
            | 28
            | 29
            | 30
            | 39
            | 40
            | 41
            | 42
            | 43
            | 44
            | 45
            | 46
            | 47
            | 48
            | 57
            | 58
            | 59
            | 60
            | 61
            | 62
            | 63
            | 64
            | 65
            | 66
            | 67
            | 68
            | 69
            | 70
            | 71
            | 72
            | 73
            | 74
            | 75
            | 76
            | 77
            | 78
            | 79
            | 80
            | 89
            | 90
            | 91
            | 92
            | 93
            | 94
            | 95
            | 96
            | 97
            | 98
            | 99
            | 100
            | 101
            | 102
            | 103
            | 104
            | 105
            | 106
            | 107
            | 108
            | 109
            | 110
            | 111
            | 112
            | 113
            | 114
            | 115
            | 116
            | 117
            | 118
    )
}

/// The metal valence the salt/coordination checks compare a bond
/// count against: `get_el_valence(el, charge, index)` (util.c:439) —
/// the `index`-th registered valence at the charge state, zero when
/// unregistered.
pub(crate) fn metal_valence(el: i32, charge: i32, index: usize) -> i32 {
    let row = el_valences(el, charge);
    row.get(index).copied().unwrap_or(0) as i32
}

/// Per-atom classification inputs the transforms read: element
/// numbers, the aromatic-aware chemical bond valence (mol2atom.c's
/// aromatic blend: 1 aromatic bond reads 1, 2 read 1.5 each (+3), 3
/// read +4 total — mol2atom.c:862-886), heavy-atom degrees and
/// derived hydrogen counts. The structure's hydrogens are
/// MATERIALIZED first (pinned `h_count`), so the transforms and the
/// layers read one frozen H state — the reference's `num_H` was
/// materialized at read time and the disconnect never re-derives it.
pub(crate) struct MetalContext {
    pub(crate) element: Vec<i32>,
    /// The aromatic-blended chemical bond valence, doubled to keep
    /// the arithmetic integral (aromatic bond = 3 here: one aromatic
    /// reads 2, two read 3 each (+6), three read 8 total).
    pub(crate) cbv2: Vec<i32>,
    pub(crate) degrees: Vec<u32>,
    pub(crate) h: Vec<u32>,
}

impl MetalContext {
    pub(crate) fn new(structure: &Structure) -> Self {
        let n = structure.atoms.len();
        let mut aromatic_count = vec![0u32; n];
        let mut cbv2 = vec![0i32; n];
        let mut degrees = vec![0u32; n];
        for &(from, to, order) in &structure.bonds {
            degrees[from] += 1;
            degrees[to] += 1;
            let weight = if order == 4 { 0 } else { 2 * order as i32 };
            cbv2[from] += weight;
            cbv2[to] += weight;
            if order == 4 {
                aromatic_count[from] += 1;
                aromatic_count[to] += 1;
            }
        }
        for index in 0..n {
            cbv2[index] += match aromatic_count[index] {
                0 => 0,
                1 => 2,                                // replaced with single bonds (+1 each)
                2 => 6,                                // =A- blend: +3 each
                3 => 8,                                // =A< blend: +4 total
                _ => 2 * aromatic_count[index] as i32, // 1 or >3: singles
            };
        }
        MetalContext {
            element: structure
                .atoms
                .iter()
                .map(|atom| super::protons::el_number(&atom.element).unwrap_or(-1))
                .collect(),
            cbv2,
            degrees,
            h: structure.hydrogens(),
        }
    }
}

/// The per-atom read-time implicit-hydrogen count (util.c:862
/// `get_num_H`): the smallest registered valence at the atom's
/// CHARGE state that absorbs the chemical valence (anion rows
/// included — a bare O⁻ reads one H), with the reference's N(IV) and
/// S(III) special cases. The `bSkipAddingH` elements never derive
/// (`util.c:912`), charges outside the cValence table keep
/// inp_num_H — zero (util.c:66-68) — and the H-count code field
/// never feeds this (`MOL_FMT_QUERY` is absent). Radicals are not
/// read by the corpus reader, so the radical branch does not arise.
/// `cbv` is the chemical bond valence INCLUDING one per folded
/// explicit H (a bonded H atom in the reference).
pub(crate) fn get_num_h(el: i32, charge: i32, cbv: i32) -> u32 {
    if skip_adding_h(el) {
        return 0;
    }
    if !(-2..=2).contains(&charge) {
        return 0;
    }
    let mut val = 0i32;
    for &v in &el_valences(el, charge) {
        let v = i32::from(v);
        if v == 0 {
            break; // row exhausted
        }
        if v < cbv {
            continue; // keep walking to a valence that absorbs cbv
        }
        val = v;
        break;
    }
    // do not add H to N(IV), S(III) (util.c 946-950)
    if el == 7 && charge == 0 && val == 5 {
        val = 3;
    }
    if el == 16 && charge == 0 && val == 4 && cbv == 3 {
        val = 3;
    }
    (val - cbv).max(0) as u32
}

/// The read-time implicit-hydrogen derivation over a whole structure
/// (util.c:862 `get_num_H`, fed by mol2atom.c `calculate_valences`).
pub(crate) fn read_time_hydrogens(structure: &Structure) -> Vec<u32> {
    let context = MetalContext::new(structure);
    (0..structure.atoms.len())
        .map(|i| {
            let drawn = structure.atoms[i].h_count;
            let el = context.element[i];
            let charge = structure.atoms[i].charge;
            // chem_bonds_valence: the aromatic-blended bond valence
            // plus one per folded explicit H (a bonded H atom in the
            // reference).
            let cbv = context.cbv2[i] / 2 + drawn as i32;
            drawn + get_num_h(el, charge, cbv)
        })
        .collect()
}

/// The rA echo's unusual-valence digit (ichiprt1.c:2494
/// `needed_unusual_el_valence`): the chemical valence the read-time
/// derivation cannot reconstruct unambiguously. Zero when exactly one
/// registered valence fits the bonds without implicit H; −1 (no
/// digit) when the chemical valence itself reads zero.
pub(crate) fn needed_unusual_valence(
    el: i32,
    charge: i32,
    cbv: i32,
    num_h: u32,
    _positions: u32,
) -> i32 {
    let chem_valence = cbv + num_h as i32;
    let num_h_expected = get_num_h(el, charge, cbv) as i32;
    // The `positions` argument mirrors the reference's num_bonds
    // parameter, consulted only in the bonds_valence !=
    // actual_bonds_valence gate that never fires here: both read the
    // same aromatic-blended chemical valence (ichi_bns.c
    // nBondsValenceInpAt).
    if !(-2..=2).contains(&charge)
        || metal_valence(el, charge, 0) == 0
        || skip_adding_h(el)
        || num_h_expected != num_h as i32
    {
        if num_h as i32 == 0 && num_h_expected == 0 {
            return 0; // no H
        }
        return chem_valence; // needs H atoms added
    }
    let mut found_known = 0;
    let mut found = 0;
    let mut exact = false;
    for i in 0..5 {
        let known = metal_valence(el, charge, i);
        if known > 0 && cbv <= known {
            found_known += 1;
            if known <= chem_valence {
                found += 1;
            }
            if known == chem_valence {
                exact = true;
                break;
            }
        }
    }
    if exact && found == 1 && found_known == 1 {
        0
    } else if chem_valence != 0 {
        chem_valence
    } else {
        -1
    }
}

/// Freezes the read-time implicit hydrogens onto every UNPINNED atom
/// (`h_count = derived`, `h_pinned = true`): the reference's num_H
/// was fixed at read time and the disconnect never re-derives it —
/// the ammine ligand the reader gave its 2 hydrogens keeps them as a
/// bare NH2⁻ component, and the bare =O the reader gave none stays a
/// bare O. DoNotAddH pins every atom drawn-only first, so nothing
/// derives there. The `bSkipAddingH` elements (util.c:912) never took
/// read-time implicit H — their drawn count stands, pinned.
pub fn materialize_hydrogens(structure: &mut Structure) {
    let derived = read_time_hydrogens(structure);
    for (atom, h) in structure.atoms.iter_mut().zip(&derived) {
        if !atom.h_pinned {
            atom.h_count = *h;
            atom.h_pinned = true;
        }
    }
}

/// Whether the (N)H4-O-C ammonium-salt pattern holds for atom `i`
/// (strutil.c:2300 `bIsAmmoniumSalt`): a nitrogen whose heavy-atom
/// degree plus hydrogen count reads 5, every heavy neighbor either
/// the ONE bridging O (2 heavy neighbors, neutral) or the ONE
/// bridging halogen (single-bonded, neutral, hydrogen-free), any
/// other neighbor rejecting. Returns the bridging atom's index.
fn ammonium_salt_site(context: &MetalContext, structure: &Structure, i: usize) -> Option<usize> {
    if context.element[i] != 7 {
        return None; // not nitrogen
    }
    let hydrogens = context.h[i];
    if context.degrees[i] + hydrogens != 5 {
        return None;
    }
    let mut site = None;
    let mut found = 0u32;
    for neighbor in neighbors_of(structure, i) {
        let (degree, charge, el) = (
            context.degrees[neighbor],
            structure.atoms[neighbor].charge,
            context.element[neighbor],
        );
        // The common per-neighbor reject (strutil.c 2287-2292): any
        // neighbor carrying implicit H rejects; a charged neighbor
        // rejects unless the O bridge whose charge cancels the N's.
        if context.h[neighbor] != 0
            || (charge != 0 && (el != 8 || charge + structure.atoms[i].charge != 0))
        {
            return None;
        }
        if matches!(el, 9 | 17 | 35 | 53)
            && degree == 1
            && context.cbv2[neighbor] == 2
            && charge == 0
            && found == 0
        {
            found += 1; // the N-X halogen bridge (first only)
            site = Some(neighbor);
        } else if el == 8
            && degree == 2
            && site.is_none()
            && carbon_side(structure, context, neighbor, i)
        {
            found += 1; // the N-O-C oxygen bridge
            site = Some(neighbor);
        } else {
            return None; // any other neighbor rejects
        }
    }
    // the final gate (strutil.c 2320-2323): exactly one O/halide
    // bridge and num_H == 4 — the quaternary ammonium.
    if found == 1 && context.h[i] == 4 {
        site
    } else {
        None
    }
}

/// The O bridge's far side must be an UNCHARGED carbon (strutil.c
/// 2440-2450: `at[iC].el_number == EL_NUMBER_C && !charge`, the H and
/// valence tests commented out upstream).
fn carbon_side(structure: &Structure, context: &MetalContext, o: usize, n: usize) -> bool {
    neighbors_of(structure, o)
        .iter()
        .copied()
        .find(|&other| other != n)
        .is_some_and(|carbon| context.element[carbon] == 6 && structure.atoms[carbon].charge == 0)
}

/// Whether the M-X / M-O-C(=) metal-salt pattern holds for metal
/// atom `i` (strutil.c:2511 `bIsMetalSalt`): a bonded, hydrogen-free
/// metal at a REGISTERED valence for its charge state, every
/// neighbor a neutral single-bonded halide or the M-O-C oxygen (2
/// heavy neighbors, neutral, hydrogen-free) whose far side is an
/// uncharged hydrogen-free carbon at chemical valence 4 with a
/// multiple bond in play.
fn is_metal_salt(context: &MetalContext, structure: &Structure, i: usize) -> bool {
    let el = context.element[i];
    let (type_bits, valence, charge) = (el_type(el), context.degrees[i], structure.atoms[i].charge);
    if valence == 0 || type_bits == 0 || context.h[i] != 0 {
        return false;
    }
    // The valence gate (strutil.c 2533-2546): neutral metals read
    // their lowest (METAL) or next-lowest (METAL2) valence; cationic
    // metals read the lowest valence at their charge.
    let normal_valence = match charge {
        0 => {
            let lowest = metal_valence(el, 0, 0) as u32;
            (type_bits & 1 != 0 && valence == lowest)
                || (type_bits & 2 != 0 && valence == metal_valence(el, 0, 1) as u32 && valence != 0)
        }
        charge if charge > 0 => {
            type_bits & 1 != 0 && valence == metal_valence(el, charge, 0) as u32
        }
        _ => false,
    };
    if !normal_valence || valence != context.cbv2[i] as u32 / 2 {
        return false; // multi bonds or unregistered valence reject
    }
    for neighbor in neighbors_of(structure, i) {
        let (degree, charge, el) = (
            context.degrees[neighbor],
            structure.atoms[neighbor].charge,
            context.element[neighbor],
        );
        let halide = matches!(el, 9 | 17 | 35 | 53)
            && degree == 1
            && context.cbv2[neighbor] == 2
            && charge == 0
            && context.h[neighbor] == 0;
        if halide {
            continue;
        }
        let oxide = el == 8
            && degree == 2
            && charge == 0
            && context.h[neighbor] == 0
            && context.cbv2[neighbor] == 4 // two single bonds
            && carbon_side_multiple(structure, context, neighbor, i);
        if !oxide {
            return false;
        }
    }
    true
}

/// The M-O-C carbon: hydrogen-free, uncharged, chemical valence 4
/// with its degree below it (`chem_bonds_valence == 4` and
/// `valence != chem_bonds_valence` — strutil.c 2560-2570). Explicit
/// hydrogens cannot exist post-fold.
fn carbon_side_multiple(
    structure: &Structure,
    context: &MetalContext,
    o: usize,
    metal: usize,
) -> bool {
    neighbors_of(structure, o)
        .iter()
        .copied()
        .find(|&other| other != metal)
        .is_some_and(|carbon| {
            context.element[carbon] == 6
                && structure.atoms[carbon].charge == 0
                && context.h[carbon] == 0
                && context.cbv2[carbon] == 8
                && context.degrees[carbon] != 4
        })
}

/// The shared neighbor list (bond order irrelevant).
fn neighbors_of(structure: &Structure, atom: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for &(from, to, _) in &structure.bonds {
        if from == atom {
            out.push(to);
        } else if to == atom {
            out.push(from);
        }
    }
    out
}

/// Removes the first bond between `a` and `b`.
fn remove_bond(structure: &mut Structure, a: usize, b: usize) -> bool {
    match structure
        .bonds
        .iter()
        .position(|&(from, to, _)| (from == a && to == b) || (from == b && to == a))
    {
        Some(position) => {
            structure.bonds.remove(position);
            true
        }
        None => false,
    }
}

/// Runs `f` over a fresh [`MetalContext`] (every salt test reads the
/// structure's CURRENT state — the salt pass mutates as it walks,
/// strutil.c's live `at` array).
fn context_fresh<T>(structure: &Structure, f: impl FnOnce(&MetalContext) -> T) -> T {
    f(&MetalContext::new(structure))
}

/// The salt disconnection (strutil.c:2586 `DisconnectSalts`): walks
/// the atoms in index order, each atom taking the ammonium pattern
/// first, then the metal-salt one; the walk mutates as it goes.
/// Returns whether anything disconnected. The structure's hydrogens
/// must be materialized ([`materialize_hydrogens`]).
pub fn disconnect_salts(structure: &mut Structure) -> bool {
    let mut changed = false;
    for i in 0..structure.atoms.len() {
        // The walk's guard (strutil.c 2597-2601): a bonded atom with
        // only single bonds and no radical reaches the tests.
        if context_fresh(structure, |context| {
            context.degrees[i] == 0 || context.degrees[i] * 2 != context.cbv2[i] as u32
        }) {
            continue;
        }
        // bIsAmmoniumSalt + DisconnectAmmoniumSalt (strutil.c 2300,
        // 2405): break the N-X bridge, annihilate cancelling charges,
        // move one N hydrogen onto the bridge atom.
        let site = context_fresh(structure, |context| {
            ammonium_salt_site(context, structure, i)
        });
        if let Some(site) = site {
            let n_charge = structure.atoms[i].charge;
            let site_charge = structure.atoms[site].charge;
            if n_charge != 0 && n_charge + site_charge == 0 {
                structure.atoms[i].charge = 0;
                structure.atoms[site].charge = 0;
            }
            if remove_bond(structure, i, site) && structure.atoms[i].h_count > 0 {
                structure.atoms[i].h_count -= 1;
                structure.atoms[site].h_count += 1;
            }
            changed = true;
            continue;
        }
        // bIsMetalSalt + DisconnectMetalSalt (strutil.c 2511, 2612):
        // every bond of the salt metal breaks; each partner takes -1,
        // the metal +1 per bond. NEVER reconnectable (TECHMAN 6.1) —
        // the /r layer never covers a salt-only disconnection.
        if context_fresh(structure, |context| is_metal_salt(context, structure, i)) {
            for partner in neighbors_of(structure, i) {
                if remove_bond(structure, i, partner) {
                    structure.atoms[partner].charge = -1;
                    structure.atoms[i].charge += 1;
                }
            }
            changed = true;
        }
    }
    changed
}

/// `bIsMetalToDisconnect` with the reference's default
/// `bCheckMetalValence = 0` (mode.h:588): a metal with any bond or
/// hydrogen marks — unless an ammonium or metal SALT site claims it
/// (those take the salt pass, strutil.c 2933-2941).
fn metal_marks(context: &MetalContext, structure: &Structure, i: usize) -> bool {
    if !is_metal_el(context.element[i]) {
        return false;
    }
    if context.degrees[i] + context.h[i] == 0 {
        return false; // nothing to disconnect
    }
    let salt =
        ammonium_salt_site(context, structure, i).is_some() || is_metal_salt(context, structure, i);
    !salt
}

/// The coordination-disconnection PROBE (strutil.c:2762
/// `bMayDisconnectMetals`): whether any metal atom would mark. The
/// RecMet layer needs the answer BEFORE the destructive pass — the
/// reconnected copy is taken at exactly this point in the reference
/// (runichi3.c:672-683).
pub fn may_disconnect_metals(structure: &Structure) -> bool {
    let context = MetalContext::new(structure);
    (0..structure.atoms.len()).any(|i| metal_marks(&context, structure, i))
}

/// The coordination disconnection (strutil.c:2865 `DisconnectMetals`):
/// marks the metals the salt pass did not claim, breaks every
/// metal-ligand bond (non-metal ligands first, then metal-metal),
/// and renormalizes each heteroatom ligand's charge to its natural
/// valence (strutil.c:3112 `DisconnectOneLigand`) — the charge delta
/// riding to the metal. Returns whether anything disconnected.
pub fn disconnect_metals(structure: &mut Structure) -> bool {
    let context = MetalContext::new(structure);
    let marked: Vec<bool> = (0..structure.atoms.len())
        .map(|i| metal_marks(&context, structure, i))
        .collect();
    if !marked.iter().any(|&m| m) {
        return false;
    }
    // Non-metal ligands first (strutil.c 3035-3050), then M-M (3060-
    // 3075): each ligand disconnects ONCE, taking ALL its bonds to
    // marked metals with it.
    let mut changed = false;
    for pass in [false, true] {
        for i in 0..structure.atoms.len() {
            if !marked[i] {
                continue;
            }
            for ligand in neighbors_of(structure, i) {
                if marked[ligand] == pass && disconnect_one_ligand(structure, &marked, i, ligand) {
                    changed = true;
                }
            }
        }
    }
    changed
}

/// `DisconnectOneLigand` (strutil.c:3112): break EVERY bond between
/// the ligand and marked metals, then renormalize the ligand charge
/// to a registered valence, the delta riding to THE one metal it was
/// bonded to. The reference's bridging-ligand defect reproduces as
/// written: a ligand bonded to TWO marked metals renormalizes never
/// (the `1 == num_metal_neigh` guard, strutil.c 3269-3272), so its
/// charge survives the split unchanged.
fn disconnect_one_ligand(
    structure: &mut Structure,
    marked: &[bool],
    metal: usize,
    ligand: usize,
) -> bool {
    let metal_bonds: Vec<usize> = structure
        .bonds
        .iter()
        .enumerate()
        .filter_map(|(index, &(from, to, _))| {
            ((from == ligand && marked[to]) || (to == ligand && marked[from])).then_some(index)
        })
        .collect();
    if metal_bonds.is_empty() {
        return false;
    }
    let num_metal_neigh = metal_bonds.len();
    let ligand_aromatic_bonds = structure
        .bonds
        .iter()
        .filter(|&&(from, to, order)| (from == ligand || to == ligand) && order > 3)
        .count();
    let removed_aromatic = metal_bonds
        .iter()
        .filter(|&&index| structure.bonds[index].2 > 3)
        .count();
    for index in metal_bonds.iter().rev() {
        structure.bonds.remove(*index);
    }
    // The renormalization gate (strutil.c 3214-3220): a heteroatom
    // ligand whose REMAINING aromatic bonds read 0, 2 or 3 — any
    // other count, or a non-heteroatom (carbon, a metal), keeps its
    // charge.
    let ligand_el = super::protons::el_number(&structure.atoms[ligand].element).unwrap_or(-1);
    let heteroat = matches!(
        ligand_el,
        9 | 17 | 35 | 53 | 85 | 7 | 15 | 33 | 8 | 16 | 34 | 52 | 5
    );
    let remaining_aromatic = ligand_aromatic_bonds - removed_aromatic;
    if (remaining_aromatic != 0 && remaining_aromatic != 2 && remaining_aromatic != 3) || !heteroat
    {
        return true;
    }
    // The natural-valence search (strutil.c 3221-3250): the ligand's
    // REMAINING chemical valence (chem_bonds_valence: single 1,
    // double 2, triple 3; the aromatic blend of the remaining bonds —
    // 2 aromatics read 1.5 each, 3 read 4 total, other counts gated
    // above) plus the read-time hydrogen count, matched against the
    // registered valences at charges -1, 0, +1 in that order, the
    // FIRST match naming the new charge. A valence of 0 renormalizes
    // a HALOGEN to -1 (any other element stays).
    let mut aromatic = 0u32;
    let mut cbv: i32 = 0;
    for &(from, to, order) in &structure.bonds {
        if from == ligand || to == ligand {
            if order > 3 {
                aromatic += 1;
            } else {
                cbv += order as i32;
            }
        }
    }
    cbv += match aromatic {
        2 => 3, // two remaining aromatic bonds read 1.5 each
        3 => 4, // three read 4 total
        _ => 0, // single bonds only (1 and >3 gated above)
    };
    let val = cbv + structure.atoms[ligand].h_count as i32;
    let new_charge = if val == 0 {
        matches!(ligand_el, 9 | 17 | 35 | 53 | 85).then_some(-1)
    } else {
        (-1..=1).find(|&charge| metal_valence(ligand_el, charge, 0) == val)
    };
    let Some(new_charge) = new_charge else {
        return true;
    };
    // The transfer (strutil.c 3260-3293): only when the charge
    // actually changes and the ligand had exactly one marked metal.
    if new_charge == structure.atoms[ligand].charge || num_metal_neigh != 1 {
        return true;
    }
    // The =N= carve-out (strutil.c 3271-3275): a cumulene nitrogen
    // whose two surviving bonds read the same type stays uncharged.
    let neighbors = neighbors_of(structure, ligand);
    if new_charge == 1 && val == 4 && neighbors.len() == 2 && cbv == 4 {
        let order_to = |neighbor: usize| -> Option<u32> {
            structure.bonds.iter().find_map(|&(from, to, order)| {
                ((from == ligand && to == neighbor) || (to == ligand && from == neighbor))
                    .then_some(order)
            })
        };
        if order_to(neighbors[0]) == order_to(neighbors[1]) {
            return true; // do not charge the disconnected =N= (2004-10-27)
        }
    }
    structure.atoms[metal].charge -= new_charge - structure.atoms[ligand].charge;
    structure.atoms[ligand].charge = new_charge;
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_metal_flags_reproduce_the_eldata_columns() {
        // Spot rows read straight off util.c's ElData dump.
        assert_eq!(el_type(3), METAL); // Li
        assert_eq!(el_type(78), METAL2); // Pt
        assert_eq!(el_type(50), METAL2); // Sn
        assert_eq!(el_type(51), METAL); // Sb — the reference marks it METAL
        assert_eq!(el_type(6), 0); // C
        assert_eq!(el_type(85), 0); // At
        assert!(is_metal_el(78));
        assert!(is_metal_el(49)); // In — missing from the proton module's old list
        assert!(!is_metal_el(85));
        assert!(skip_adding_h(26)); // Fe
        assert!(!skip_adding_h(3)); // Li
        assert!(!skip_adding_h(50)); // Sn
    }

    #[test]
    fn the_metal_valence_rows_match_the_reference() {
        // get_el_valence spot rows: the valence the salt check keys on.
        assert_eq!(metal_valence(3, 0, 0), 1); // Li neutral 1
        assert_eq!(metal_valence(29, 0, 0), 1); // Cu neutral {1,2}
        assert_eq!(metal_valence(29, 0, 1), 2);
        assert_eq!(metal_valence(78, 0, 0), 2); // Pt neutral {2,4}
        assert_eq!(metal_valence(78, 0, 1), 4);
        assert_eq!(metal_valence(7, -1, 0), 2); // N(-) {2} — the ammine renorm
        assert_eq!(metal_valence(8, -1, 0), 1); // O(-) {1} — the oxy renorm
        assert_eq!(metal_valence(26, 3, 0), 0); // charge 3 out of range
    }

    #[test]
    fn the_aromatic_blend_reproduces_mol2atom() {
        // The recorded phenoxide pattern (identity/metal-sodium-
        // phenoxide-bonded): the O carries one Na single + one
        // aromatic bond — the 1-aromatic case replaces the aromatic
        // with a single, cbv 2; the ipso C carries THREE aromatic
        // bonds — the 3-aromatic blend reads cbv 4, the value the
        // M-O-C salt test keys on.
        let s = crate::parse_molfile(
            "pheno\n  chemomer\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    1.0    1.0    0.0 Na  0  0  0  0  0  0  0  0  0  0  0  0  0\n    2.0    1.0    0.0 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n    3.0    1.0    0.0 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n    4.0    1.0    0.0 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n    5.0    1.0    0.0 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0  0  0  0\n  2  3  4  0  0  0  0  0  0  0  0\n  3  4  4  0  0  0  0  0  0  0  0\n  3  5  4  0  0  0  0  0  0  0  0\nM  END\n",
        )
        .expect("parses");
        let context = MetalContext::new(&s);
        assert_eq!(context.cbv2[0], 2); // Na: one single
        assert_eq!(context.cbv2[1], 4); // O: single + aromatic-as-single = 2
        assert_eq!(context.cbv2[2], 8); // ipso C: three aromatics = 4
    }
}
