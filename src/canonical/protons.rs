//! The (DE)PROTONATION normalization (the reference's
//! RemoveNPProtonsAndAcidCharges, TECHMAN 5.1/5.2): charged
//! structures have acidic protons removed from neutral heteroatoms
//! (net charge > 0) or protons added to anionic heteroatoms (net
//! charge < 0) so the standardized skeleton reads as the conjugate
//! form; the shifted protons surface later as the /p layer.
//!
//! Scope of this slice: the SIMPLE variants (5.1a SimpleRemoveHplusNPO,
//! 5.2a SimpleRemoveAcidicProtons, SimpleAddAcidicProtons) with the
//! reference's priority masks (AR_SIMPLE/AA_SIMPLE step 1; steps 2-3
//! are compiled out under AR_ANY_OH=0/AA_ANY_O_Minus=0). The HARD
//! variants (BNS-checked removal, TECHMAN 5.1b/5.2b, on under
//! bTgFlagHardAddRenProtons=1) are pending; they extend the same
//! classification with alt-path tests.

use crate::graph::Structure;

/// AT_* type bits (inpdef.h ATT_*).
mod att {
    pub const NONE: u16 = 0x0000;
    pub const ACIDIC_CO: u16 = 0x0001;
    pub const ACIDIC_S: u16 = 0x0002;
    pub const OO: u16 = 0x0004;
    pub const ZOO: u16 = 0x0008;
    pub const NO: u16 = 0x0010;
    pub const N_O: u16 = 0x0020;
    pub const ATOM_N: u16 = 0x0040;
    pub const ATOM_P: u16 = 0x0080;
    pub const OTHER_NEG_O: u16 = 0x0100;
    pub const OTHER_ZO: u16 = 0x0200;
    pub const OH_MINUS: u16 = 0x0400;
    pub const O_PLUS: u16 = 0x0800;
    pub const PROTON: u16 = 0x1000;
    pub const HAL_ANION: u16 = 0x2000;
    pub const HAL_ACID: u16 = 0x4000;
}

/// ATBIT_* mask bits: bit i = the atom counts in ATTOT total i.
mod atbit {
    use super::attot;
    pub const NP_PLUS: u32 = 1 << attot::NP_PLUS;
    pub const NP_PROTON: u32 = 1 << attot::NP_PROTON;
    pub const NP_H: u32 = 1 << attot::NP_H;
    pub const N_MINUS: u32 = 1 << attot::N_MINUS;
    pub const NP: u32 = 1 << attot::NP;
    pub const ON: u32 = 1 << attot::ON;
    pub const COH: u32 = 1 << attot::COH;
    pub const CSH: u32 = 1 << attot::CSH;
    pub const ZOH: u32 = 1 << attot::ZOH;
    pub const OOH: u32 = 1 << attot::OOH;
    pub const ZOOH: u32 = 1 << attot::ZOOH;
    pub const NOH: u32 = 1 << attot::NOH;
    pub const N_OH: u32 = 1 << attot::N_OH;
    pub const CO: u32 = 1 << attot::CO;
    pub const ZO: u32 = 1 << attot::ZO;
    pub const NO2: u32 = 1 << attot::NO2;
    pub const N_O: u32 = 1 << attot::N_O;
    pub const CO_MINUS: u32 = 1 << attot::CO_MINUS;
    pub const CS_MINUS: u32 = 1 << attot::CS_MINUS;
    pub const ZO_MINUS: u32 = 1 << attot::ZO_MINUS;
    pub const OO_MINUS: u32 = 1 << attot::OO_MINUS;
    pub const ZOO_MINUS: u32 = 1 << attot::ZOO_MINUS;
    pub const NO_MINUS: u32 = 1 << attot::NO_MINUS;
    pub const N_O_MINUS: u32 = 1 << attot::N_O_MINUS;
    pub const O_MINUS: u32 = 1 << attot::O_MINUS;
    pub const OH_PLUS: u32 = 1 << attot::OH_PLUS;
    pub const O_PLUS: u32 = 1 << attot::O_PLUS;
    pub const PROTON: u32 = 1 << attot::PROTON;
    pub const HAL_ANION: u32 = 1 << attot::HAL_ANION;
    pub const HAL_ACID: u32 = 1 << attot::HAL_ACID;
    pub const ERRORS: u32 = 1 << attot::ERRORS;
}

/// ATTOT_* total slots (only the ones the driver reads).
mod attot {
    pub const NP_PLUS: usize = 0;
    pub const NP_PROTON: usize = 1;
    #[allow(dead_code)]
    pub const NP_H: usize = 2;
    pub const N_MINUS: usize = 3;
    #[allow(dead_code)]
    pub const NP: usize = 4;
    #[allow(dead_code)]
    pub const ON: usize = 5;
    #[allow(dead_code)]
    pub const COH: usize = 6;
    #[allow(dead_code)]
    pub const CSH: usize = 7;
    #[allow(dead_code)]
    pub const ZOH: usize = 8;
    #[allow(dead_code)]
    pub const OOH: usize = 9;
    #[allow(dead_code)]
    pub const ZOOH: usize = 10;
    #[allow(dead_code)]
    pub const NOH: usize = 11;
    #[allow(dead_code)]
    pub const N_OH: usize = 12;
    #[allow(dead_code)]
    pub const CO: usize = 13;
    #[allow(dead_code)]
    pub const ZO: usize = 14;
    #[allow(dead_code)]
    pub const NO2: usize = 15;
    #[allow(dead_code)]
    pub const N_O: usize = 16;
    #[allow(dead_code)]
    pub const CO_MINUS: usize = 17;
    #[allow(dead_code)]
    pub const CS_MINUS: usize = 18;
    #[allow(dead_code)]
    pub const ZO_MINUS: usize = 19;
    #[allow(dead_code)]
    pub const OO_MINUS: usize = 20;
    #[allow(dead_code)]
    pub const ZOO_MINUS: usize = 21;
    #[allow(dead_code)]
    pub const NO_MINUS: usize = 22;
    #[allow(dead_code)]
    pub const N_O_MINUS: usize = 23;
    #[allow(dead_code)]
    pub const O_MINUS: usize = 24;
    #[allow(dead_code)]
    pub const OH_PLUS: usize = 25;
    #[allow(dead_code)]
    pub const O_PLUS: usize = 26;
    #[allow(dead_code)]
    pub const PROTON: usize = 27;
    #[allow(dead_code)]
    pub const HAL_ANION: usize = 28;
    #[allow(dead_code)]
    pub const HAL_ACID: usize = 29;
    pub const ERRORS: usize = 30;
    pub const TOT_CHARGE: usize = 31;
    pub const NUM_CHARGES: usize = 32;
    pub const ARRAY_LEN: usize = 33;
}

/// PR_SIMPLE: protons off charged N,P and OH(+).
const PR_SIMPLE_MSK: u32 = atbit::NP_PROTON | atbit::OH_PLUS;
const PR_SIMPLE_TYP: u16 = att::ATOM_N | att::ATOM_P | att::O_PLUS;

/// AR_SIMPLE step 1 (steps 2-3 compiled out: AR_ANY_OH=0).
const AR_SIMPLE_MSK1: u32 =
    atbit::COH | atbit::CSH | atbit::OOH | atbit::ZOOH | atbit::NOH | atbit::HAL_ACID;
const AR_SIMPLE_TYP1: u16 =
    att::ACIDIC_CO | att::ACIDIC_S | att::OO | att::ZOO | att::NO | att::HAL_ACID;

/// AA_SIMPLE step 1 (steps 2-3 compiled out: AA_ANY_O_Minus=0).
const AA_SIMPLE_MSK1: u32 = atbit::CO_MINUS
    | atbit::CS_MINUS
    | atbit::OO_MINUS
    | atbit::ZOO_MINUS
    | atbit::NO_MINUS
    | atbit::O_MINUS
    | atbit::HAL_ANION;
const AA_SIMPLE_TYP1: u16 =
    att::ACIDIC_CO | att::ACIDIC_S | att::OO | att::ZOO | att::NO | att::OH_MINUS | att::HAL_ANION;

const AR_STEPS: [(u16, u32); 1] = [(AR_SIMPLE_TYP1, AR_SIMPLE_MSK1)];
const AA_STEPS: [(u16, u32); 1] = [(AA_SIMPLE_TYP1, AA_SIMPLE_MSK1)];

/// Element number for the atoms this phase classifies; None for
/// anything outside the table (treated as a non-Z, non-metal other).
pub(crate) fn el_number(element: &str) -> Option<i32> {
    let n = match element {
        "H" | "D" | "T" => 1,
        "He" => 2,
        "Li" => 3,
        "Be" => 4,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "Ne" => 10,
        "Na" => 11,
        "Mg" => 12,
        "Al" => 13,
        "Si" => 14,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "Ar" => 18,
        "K" => 19,
        "Ca" => 20,
        "Sc" => 21,
        "Ti" => 22,
        "V" => 23,
        "Cr" => 24,
        "Mn" => 25,
        "Fe" => 26,
        "Co" => 27,
        "Ni" => 28,
        "Cu" => 29,
        "Zn" => 30,
        "Ga" => 31,
        "Ge" => 32,
        "As" => 33,
        "Se" => 34,
        "Br" => 35,
        "Kr" => 36,
        "Rb" => 37,
        "Sr" => 38,
        "Y" => 39,
        "Zr" => 40,
        "Nb" => 41,
        "Mo" => 42,
        "Tc" => 43,
        "Ru" => 44,
        "Rh" => 45,
        "Pd" => 46,
        "Ag" => 47,
        "Cd" => 48,
        "In" => 49,
        "Sn" => 50,
        "Sb" => 51,
        "Te" => 52,
        "I" => 53,
        "Xe" => 54,
        "Cs" => 55,
        "Ba" => 56,
        "La" => 57,
        "Lu" => 71,
        "Hf" => 72,
        "Ta" => 73,
        "W" => 74,
        "Re" => 75,
        "Os" => 76,
        "Ir" => 77,
        "Pt" => 78,
        "Au" => 79,
        "Hg" => 80,
        "Tl" => 81,
        "Pb" => 82,
        "Bi" => 83,
        "Po" => 84,
        "At" => 85,
        "Rn" => 86,
        "Fr" => 87,
        "Ra" => 88,
        "Ac" => 89,
        "Th" => 90,
        _ => return None,
    };
    Some(n)
}

/// METAL/METAL2 nType flags, from the reference's ElData table.
fn is_metal(el: i32) -> bool {
    matches!(
        el,
        3 | 4
            | 11
            | 12
            | 13
            | 19
            | 20
            | 21
            | 22
            | 23
            | 24
            | 25
            | 26
            | 27
            | 28
            | 29
            | 30
            | 31
            | 37
            | 38
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
            | 55
            | 56
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
            | 81
            | 82
            | 83
            | 84
            | 85
            | 86
            | 87
            | 88
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
    )
}

/// Z atoms: the centers the O/S/Se/Te terminal classification keys on.
/// The reference's `get_unusual_el_valence`: the chemical valence
/// (`bonds_valence + num_h`) unless exactly one registered valence for
/// the (element, charge) state fits between the bond count and it —
/// the aux `/CRV:` suffix that flags bonds the reader could not
/// otherwise reconstruct.
pub(crate) fn unusual_el_valence(
    el: i32,
    charge: i32,
    bonds_valence: i32,
    num_h: u32,
    num_bonds: u32,
) -> i32 {
    if num_bonds == 0 && num_h == 0 {
        return 0;
    }
    if !(-2..=2).contains(&charge) {
        return if bonds_valence == num_bonds as i32 {
            0
        } else {
            bonds_valence
        };
    }
    if el_valences(el, charge)[0] == 0 && bonds_valence == num_bonds as i32 {
        return 0;
    }
    let chem_valence = bonds_valence + num_h as i32;
    let mut num_found = 0;
    let mut exact_found = false;
    for known in el_valences(el, charge) {
        if known == 0 {
            continue;
        }
        let known = i32::from(known);
        if num_bonds as i32 <= known && known <= chem_valence {
            num_found += 1;
            if known == chem_valence {
                exact_found = true;
                break;
            }
        }
    }
    if exact_found && num_found == 1 {
        0
    } else {
        chem_valence
    }
}

fn is_z_atom(el: i32) -> bool {
    matches!(
        el,
        5 | 6 | 7 | 14 | 15 | 33 | 51 | 16 | 34 | 52 | 17 | 35 | 53
    )
}

/// Valences by (element, charge state -2..+2), from the reference's
/// cValence table. Zero = the state is not allowed.
pub(crate) fn el_valences(el: i32, charge: i32) -> [u8; 5] {
    let idx = (charge.clamp(-2, 2) + 2) as usize;
    let row: [[u8; 5]; 5] = match el {
        // H/D/T: one valence of 1 at the neutral state only (ElData
        // util.c:107) — the per-charge table never folds into a
        // single charge's valence list.
        1 => [[0; 5], [0; 5], [1, 0, 0, 0, 0], [0; 5], [0; 5]],
        2 | 10 | 18 | 36 | 54 | 86 => [[0; 5]; 5], // noble gases
        // Alkali metals: one neutral valence only (ElData rows).
        3 | 11 | 19 | 37 | 55 | 87 => [[0; 5], [0; 5], [1, 0, 0, 0, 0], [0; 5], [0; 5]],
        4 => [[0; 5], [0; 5], [2, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0; 5]], // Be
        5 => [
            [3, 0, 0, 0, 0],
            [4, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [2, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
        ], // B
        6 => [
            [2, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [4, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [2, 0, 0, 0, 0],
        ], // C
        7 => [
            [1, 0, 0, 0, 0],
            [2, 0, 0, 0, 0],
            [3, 5, 0, 0, 0],
            [4, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
        ], // N
        8 => [
            [0; 5],
            [1, 0, 0, 0, 0],
            [2, 0, 0, 0, 0],
            [3, 5, 0, 0, 0],
            [4, 0, 0, 0, 0],
        ], // O
        9 => [
            [0; 5],
            [0; 5],
            [1, 0, 0, 0, 0],
            [2, 0, 0, 0, 0],
            [3, 5, 0, 0, 0],
        ], // F
        12 => [[0; 5], [0; 5], [2, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0; 5]], // Mg
        13 => [
            [3, 5, 0, 0, 0],
            [4, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [2, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
        ], // Al
        14 => [
            [2, 0, 0, 0, 0],
            [3, 5, 0, 0, 0],
            [4, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [2, 0, 0, 0, 0],
        ], // Si
        15 => [
            [1, 3, 5, 7, 0],
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
            [4, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
        ], // P
        16 => [
            [0; 5],
            [1, 3, 5, 7, 0],
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
            [4, 0, 0, 0, 0],
        ], // S
        17 | 35 | 53 => [
            [0; 5],
            [0; 5],
            [1, 3, 5, 7, 0],
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
        ], // Cl/Br/I
        20 => [[0; 5], [0; 5], [2, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0; 5]], // Ca
        31 => [
            [3, 5, 0, 0, 0],
            [4, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [0; 5],
            [1, 0, 0, 0, 0],
        ], // Ga
        32 => [
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
            [4, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [0; 5],
        ], // Ge
        33 => [
            [1, 3, 5, 7, 0],
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
            [4, 0, 0, 0, 0],
            [3, 0, 0, 0, 0],
        ], // As
        34 => [
            [0; 5],
            [1, 3, 5, 7, 0],
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
            [4, 0, 0, 0, 0],
        ], // Se
        // The metal rows (ElData, util.c:103-277). Only the neutral
        // column carries registered valences for most metals; the
        // main-group metals carry the charged-state rows their salt
        // chemistry reads.
        21 | 39 | 43 | 57 | 60 | 61 | 64 | 66 | 67 | 68 | 71 | 89 | 96 | 98 | 99 | 100 | 101 => {
            [[0; 5], [0; 5], [3, 0, 0, 0, 0], [0; 5], [0; 5]]
        } // Sc Y Tc La Nd Pm Gd Dy Ho Er Lu Ac Cm Cf Es Fm Md
        22 => [[0; 5], [0; 5], [3, 4, 0, 0, 0], [0; 5], [0; 5]], // Ti
        23 => [[0; 5], [0; 5], [2, 3, 4, 5, 0], [0; 5], [0; 5]], // V
        24 => [[0; 5], [0; 5], [2, 3, 6, 0, 0], [0; 5], [0; 5]], // Cr
        25 => [[0; 5], [0; 5], [2, 3, 4, 6, 0], [0; 5], [0; 5]], // Mn
        26 => [[0; 5], [0; 5], [2, 3, 4, 6, 0], [0; 5], [0; 5]], // Fe
        27 | 28 => [[0; 5], [0; 5], [2, 3, 0, 0, 0], [0; 5], [0; 5]], // Co Ni
        29 => [[0; 5], [0; 5], [1, 2, 0, 0, 0], [0; 5], [0; 5]], // Cu
        30 => [[0; 5], [0; 5], [2, 0, 0, 0, 0], [0; 5], [0; 5]], // Zn
        38 => [[0; 5], [0; 5], [2, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0; 5]], // Sr
        40 => [[0; 5], [0; 5], [4, 0, 0, 0, 0], [0; 5], [0; 5]], // Zr
        41 => [[0; 5], [0; 5], [3, 5, 0, 0, 0], [0; 5], [0; 5]], // Nb
        42 => [[0; 5], [0; 5], [3, 4, 5, 6, 0], [0; 5], [0; 5]], // Mo
        44 => [[0; 5], [0; 5], [2, 3, 4, 6, 0], [0; 5], [0; 5]], // Ru
        45 => [[0; 5], [0; 5], [2, 3, 4, 0, 0], [0; 5], [0; 5]], // Rh
        46 => [[0; 5], [0; 5], [2, 4, 0, 0, 0], [0; 5], [0; 5]], // Pd
        47 => [[0; 5], [0; 5], [1, 0, 0, 0, 0], [0; 5], [0; 5]], // Ag
        48 => [[0; 5], [0; 5], [2, 0, 0, 0, 0], [0; 5], [0; 5]], // Cd
        49 => [
            [3, 5, 0, 0, 0],
            [2, 4, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [0; 5],
            [1, 0, 0, 0, 0],
        ], // In
        50 => [
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
            [2, 4, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [0; 5],
        ], // Sn
        51 => [
            [1, 3, 5, 7, 0],
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
            [2, 4, 0, 0, 0],
            [3, 0, 0, 0, 0],
        ], // Sb
        56 => [[0; 5], [0; 5], [2, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0; 5]], // Ba
        58 => [[0; 5], [0; 5], [3, 4, 0, 0, 0], [0; 5], [0; 5]], // Ce
        59 => [[0; 5], [0; 5], [3, 4, 0, 0, 0], [0; 5], [0; 5]], // Pr
        62 | 63 | 69 | 70 => [[0; 5], [0; 5], [2, 3, 0, 0, 0], [0; 5], [0; 5]], // Sm Eu Tm Yb
        65 => [[0; 5], [0; 5], [3, 4, 0, 0, 0], [0; 5], [0; 5]], // Tb
        72 => [[0; 5], [0; 5], [4, 0, 0, 0, 0], [0; 5], [0; 5]], // Hf
        73 => [[0; 5], [0; 5], [5, 0, 0, 0, 0], [0; 5], [0; 5]], // Ta
        74 => [[0; 5], [0; 5], [3, 4, 5, 6, 0], [0; 5], [0; 5]], // W
        75 => [[0; 5], [0; 5], [2, 4, 6, 7, 0], [0; 5], [0; 5]], // Re
        76 | 77 => [[0; 5], [0; 5], [2, 3, 4, 6, 0], [0; 5], [0; 5]], // Os Ir
        78 => [[0; 5], [0; 5], [2, 4, 0, 0, 0], [0; 5], [0; 5]], // Pt
        79 => [[0; 5], [0; 5], [1, 3, 0, 0, 0], [0; 5], [0; 5]], // Au
        80 => [[0; 5], [0; 5], [1, 2, 0, 0, 0], [0; 5], [0; 5]], // Hg
        81 => [
            [3, 5, 0, 0, 0],
            [2, 4, 0, 0, 0],
            [1, 3, 0, 0, 0],
            [0; 5],
            [0; 5],
        ], // Tl
        82 => [
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
            [2, 4, 0, 0, 0],
            [3, 0, 0, 0, 0],
            [0; 5],
        ], // Pb
        83 => [
            [1, 3, 5, 7, 0],
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
            [2, 4, 0, 0, 0],
            [3, 0, 0, 0, 0],
        ], // Bi
        84 => [
            [0; 5],
            [1, 3, 5, 7, 0],
            [2, 4, 6, 0, 0],
            [3, 5, 0, 0, 0],
            [2, 4, 0, 0, 0],
        ], // Po
        88 => [[0; 5], [0; 5], [2, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0; 5]], // Ra
        90 => [[0; 5], [0; 5], [3, 4, 0, 0, 0], [0; 5], [0; 5]], // Th
        91 => [[0; 5], [0; 5], [3, 4, 5, 0, 0], [0; 5], [0; 5]], // Pa
        92..=95 => [[0; 5], [0; 5], [3, 4, 5, 6, 0], [0; 5], [0; 5]], // U Np Pu Am
        97 => [[0; 5], [0; 5], [3, 4, 0, 0, 0], [0; 5], [0; 5]], // Bk
        102..=118 => [[0; 5], [0; 5], [1, 0, 0, 0, 0], [0; 5], [0; 5]], // No..Og placeholders
        _ => return [0; 5],
    };
    row[idx]
}

/// The per-atom view the classifier reads, owning the standardized
/// H/charge state. The doubled bond-order sum (aromatic = 3) keeps
/// the arithmetic integral.
pub(crate) struct ProtonContext {
    pub(crate) element: Vec<i32>,
    pub(crate) neighbors: Vec<Vec<usize>>,
    pub(crate) cbv2: Vec<i32>,
    /// The bond list (from, to, doubled order) the pairing network
    /// reads.
    pub(crate) bonds: Vec<(usize, usize, i32)>,
    pub(crate) h: Vec<u32>,
    pub(crate) charge: Vec<i32>,
}

impl ProtonContext {
    pub(crate) fn new(structure: &Structure) -> Self {
        let n = structure.atoms.len();
        let mut neighbors = vec![Vec::new(); n];
        let mut cbv2 = vec![0i32; n];
        for &(from, to, order) in &structure.bonds {
            neighbors[from].push(to);
            neighbors[to].push(from);
            let o = if order == 4 { 3 } else { 2 * order as i32 };
            cbv2[from] += o;
            cbv2[to] += o;
        }
        ProtonContext {
            element: structure
                .atoms
                .iter()
                .map(|atom| el_number(&atom.element).unwrap_or(-1))
                .collect(),
            neighbors,
            cbv2,
            bonds: structure
                .bonds
                .iter()
                .map(|&(from, to, order)| (from, to, if order == 4 { 3 } else { 2 * order as i32 }))
                .collect(),
            h: structure.hydrogens(),
            charge: structure.atoms.iter().map(|atom| atom.charge).collect(),
        }
    }

    fn valence(&self, atom: usize) -> usize {
        self.neighbors[atom].len()
    }

    /// detect_unusual_el_valence: chem valence (cbv + H) must equal
    /// some allowed valence for the charge state.
    fn unusual_valence(&self, atom: usize) -> bool {
        let h = self.h[atom] as i32;
        if self.valence(atom) == 0 && h == 0 {
            return false;
        }
        let charge = self.charge[atom];
        if !(-2..=2).contains(&charge) {
            return self.cbv2[atom] != 2 * self.valence(atom) as i32;
        }
        if el_valences(self.element[atom], charge) == [0; 5]
            && self.cbv2[atom] == 2 * self.valence(atom) as i32
        {
            return false;
        }
        let chem_valence2 = self.cbv2[atom] + 2 * h;
        !el_valences(self.element[atom], charge)
            .iter()
            .any(|&v| v > 0 && 2 * v as i32 == chem_valence2)
    }
}

/// GetAtomChargeType: (type, mask) for the atom under the current
/// standardized H/charge state.
pub(crate) fn atom_charge_type(context: &ProtonContext, atom: usize) -> (u16, u32) {
    let (el, charge, h) = (context.element[atom], context.charge[atom], context.h[atom]);
    let valence = context.valence(atom);
    let mut typ = att::NONE;
    let mut mask = 0u32;

    // The reference's first two exits (ichi_bns.c:2892-2897): metals
    // never class — their charges never enter the proton-exchange —
    // and per-atom charges beyond ±1 are untypable.
    if super::metals::is_metal_el(el) {
        return (typ, mask);
    }
    if !(-1..=1).contains(&charge) {
        return (typ, mask);
    }
    if valence == 0 && charge == 1 && h == 0 && el == 1 {
        return (att::PROTON, atbit::PROTON);
    }
    if valence == 0 && charge == -1 && h == 0 && matches!(el, 9 | 17 | 35 | 53) {
        return (att::HAL_ANION, atbit::HAL_ANION);
    }
    if (valence == 0 && charge == 0 && h == 1 && matches!(el, 9 | 17 | 35 | 53))
        || (valence == 0 && charge == 0 && h == 2 && matches!(el, 8 | 16 | 34 | 52))
    {
        return (att::HAL_ACID, atbit::HAL_ACID);
    }
    if context.unusual_valence(atom) {
        return (typ, mask);
    }

    let mut num_z = 0usize;
    let mut num_o = 0usize;
    let mut unsat_n_has_terminal_o = 0usize;
    for &neigh in &context.neighbors[atom] {
        if context.charge[neigh] != 0 && charge != 0 {
            // adjacent-ion rule: the atom is neither donor nor acceptor
            return (att::NONE, 0);
        }
        if context.unusual_valence(neigh) {
            return (typ, mask);
        }
        if is_z_atom(context.element[neigh]) {
            num_z += 1;
        }
        num_o += (context.element[neigh] == 8) as usize;
        if el == 7 && valence == 2 && charge == 0 {
            let nv = context.valence(neigh);
            if nv == 1
                && context.cbv2[neigh] == 4
                && matches!(context.element[neigh], 8 | 16 | 34 | 52)
            {
                unsat_n_has_terminal_o += 1;
            }
        }
    }

    if matches!(el, 8 | 16 | 34 | 52) {
        // O, S, Se, Te terminal chemistry
        if charge == 1 {
            if h > 0 {
                return (att::O_PLUS, atbit::OH_PLUS);
            }
            return (att::O_PLUS, atbit::O_PLUS);
        }
        if valence > 1 {
            return (typ, mask);
        }
        if valence == 1 && num_z == 0 && num_o == 0 {
            return (typ, mask); // #C3: charges only, no donor/acceptor
        }
        if valence == 0 {
            if charge == -1 && h == 1 {
                return (att::OH_MINUS, atbit::O_MINUS);
            }
            return (typ, mask);
        }
        let neigh = context.neighbors[atom][0];
        let (nel, ncharge) = (context.element[neigh], context.charge[neigh]);
        let nvalence = context.valence(neigh);
        let ncbv2 = context.cbv2[neigh];
        if ncharge == 0 && nel == 6 && ncbv2 > 2 * nvalence as i32 {
            // =C-OH, =C-O(-), -C=O
            typ = att::ACIDIC_CO;
            if h == 1 {
                mask |= atbit::COH;
            } else if charge == -1 {
                mask |= atbit::CO_MINUS;
            } else if h == 0 && charge == 0 {
                mask |= atbit::CO;
            } else {
                mask |= atbit::ERRORS;
            }
        } else if ncharge == 0 && matches!(nel, 8 | 16 | 34 | 52) && ncbv2 == 2 * nvalence as i32 {
            // -O-OH, -O-O(-)
            typ = att::OO;
            if h == 1 {
                mask |= atbit::OOH;
            } else if charge == -1 {
                mask |= atbit::OO_MINUS;
            } else {
                mask |= atbit::ERRORS;
            }
        } else if ncharge == 0 && nel == 6 && ncbv2 == 2 * nvalence as i32 && el != 8 {
            // >C-SH, >C-S(-)
            typ = att::ACIDIC_S;
            if h == 1 {
                mask |= atbit::CSH;
            } else if charge == -1 {
                mask |= atbit::CS_MINUS;
            } else {
                mask |= atbit::ERRORS;
            }
        } else if nel == 7
            && nvalence == 2
            && (context.h[neigh] == 0 || (context.h[neigh] == 1 && ncharge == 1))
        {
            typ = att::NO;
            if h == 1 {
                mask |= atbit::NOH;
            } else if charge == -1 {
                mask |= atbit::NO_MINUS;
            } else if ncharge == 1 || ncharge == 0 {
                mask |= atbit::NO2;
            } else {
                mask |= atbit::ERRORS;
            }
        } else if nel == 7 {
            typ = att::N_O;
            if h == 1 {
                mask |= atbit::N_OH;
            } else if charge == -1 {
                mask |= atbit::N_O_MINUS;
            } else if ncharge == 1 {
                mask |= atbit::N_O;
            } else {
                mask |= atbit::ERRORS;
            }
        } else if nel != 6 && nel != 8 && !is_metal(nel) && ncbv2 > 2 * nvalence as i32 {
            // Z \ {N, C}: P, As, Sb, S, Se, Te, Cl, Br, I
            if context.cbv2[atom] == 2 && is_zox(context, atom, 0) > 0 {
                typ = att::ZOO;
                if h == 1 {
                    mask |= atbit::ZOOH;
                } else if charge == -1 {
                    mask |= atbit::ZOO_MINUS;
                } else {
                    mask |= atbit::ERRORS;
                }
            } else {
                typ = att::OTHER_ZO;
                if h == 1 {
                    mask |= atbit::ZOH;
                } else if charge == -1 {
                    mask |= atbit::ZO_MINUS;
                } else if h == 0 {
                    mask |= atbit::ZO;
                } else {
                    mask |= atbit::ERRORS;
                }
            }
        } else if charge == -1 && !is_metal(nel) {
            typ = att::OTHER_NEG_O;
            mask |= atbit::O_MINUS;
        }
        return (typ, mask);
    }

    if matches!(el, 7 | 15)
        && (0..=3).contains(&valence)
        && context.cbv2[atom] + 2 * h as i32 == 2 * (3 + charge)
    {
        // N, P (neutral valence 3); atoms whose neighbors are all
        // non-Z are ignored here (the C's goto count_mask_bits with
        // no type)
        if valence > 0 && num_z == 0 {
            return (att::NONE, 0);
        }
        typ = if el == 7 { att::ATOM_N } else { att::ATOM_P };
        match charge {
            -1 => {
                if el == 7 {
                    mask |= atbit::N_MINUS;
                    if h > 0 {
                        mask |= atbit::NP_H;
                    }
                }
            }
            0 => {
                if h > 0 {
                    mask |= atbit::NP_H;
                } else if unsat_n_has_terminal_o == 1 {
                    mask |= atbit::ON;
                } else {
                    mask |= atbit::NP;
                }
            }
            1 => {
                if h > 0 {
                    mask |= atbit::NP_PROTON;
                } else if context.cbv2[atom] > 2 * valence as i32 {
                    mask |= atbit::NP_PLUS;
                } else {
                    typ = att::NONE; // onium cations >N(+)< ignored
                }
            }
            _ => mask |= atbit::ERRORS,
        }
    }
    (typ, mask)
}

/// IsZOX: terminal =O/S/Se/Te count on the (single) Z neighbor.
fn is_zox(context: &ProtonContext, atom: usize, ord: usize) -> i32 {
    let z = context.neighbors[atom][ord];
    let mut count = 0;
    for &neigh in &context.neighbors[z] {
        if neigh == atom {
            continue;
        }
        if context.valence(neigh) == 1
            && context.cbv2[neigh] == 4
            && context.charge[neigh] == 0
            && matches!(context.element[neigh], 8 | 16 | 34 | 52)
        {
            count += 1;
        }
    }
    count
}

/// Recompute the totals the driver reads. The C maintains them
/// incrementally through the classifier's delta parameter; a fresh
/// scan is equivalent and keeps the mutation sites dumb.
fn totals(context: &ProtonContext) -> [i64; attot::ARRAY_LEN] {
    let mut totals = [0i64; attot::ARRAY_LEN];
    for atom in 0..context.element.len() {
        let (_, mask) = atom_charge_type(context, atom);
        if mask != 0 && mask & atbit::ERRORS == 0 {
            for (slot, total) in totals.iter_mut().enumerate().take(32) {
                if mask & (1 << slot) != 0 {
                    *total += 1;
                }
            }
        }
        if context.charge[atom] != 0 {
            totals[attot::TOT_CHARGE] += context.charge[atom] as i64;
            totals[attot::NUM_CHARGES] += 1;
        }
    }
    totals
}

/// bHasAcidicHydrogen: a neutral atom carrying H that matches the
/// acidic-removal masks (the merge-phase donor class).
pub(crate) fn has_acidic_hydrogen(context: &ProtonContext, atom: usize) -> bool {
    if context.charge[atom] != 0 || context.h[atom] == 0 {
        return false;
    }
    let (typ, mask) = atom_charge_type(context, atom);
    typ & AR_SIMPLE_TYP1 != 0 && mask & AR_SIMPLE_MSK1 != 0
}

/// bHasAcidicMinus: an anion matching the acidic-addition masks (the
/// merge-phase acceptor class).
pub(crate) fn has_acidic_minus(context: &ProtonContext, atom: usize) -> bool {
    if context.charge[atom] != -1 {
        return false;
    }
    let (typ, mask) = atom_charge_type(context, atom);
    typ & AA_SIMPLE_TYP1 != 0 && mask & AA_SIMPLE_MSK1 != 0
}

/// The standardization result: adjusted H counts and charges per atom
/// plus the net proton shift for the /p accounting. `bonds` carries
/// the bond orders the HARD removal's π alternation moved (the
/// reference's BNS ride shifts orders along the route — the acid's
/// =O reads single and the donor-side π lands one bond over); empty
/// when no hard pass fired.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct StandardizedProtons {
    pub h: Vec<u32>,
    pub charge: Vec<i32>,
    pub removed_protons: u32,
    pub added_protons: u32,
    /// (a, b, new order, ride mark) for the bonds whose order moved;
    /// `a < b`. The mark is the bond TYPE the reference's ride leaves
    /// (extr_ct.h): ALT12NS = 9 ("1 or 2", either salt role) on
    /// single/double-sourced moves, ALTERN|TRIPLE = 7 ("2 or 3",
    /// neither role) on triple-sourced ones.
    pub bonds: Vec<(usize, usize, u32, u8)>,
    /// A HARD (DE)PROTONATION phase moved a proton (the 5.1b
    /// charge-point ride-out or the 5.2b acidic hard removal) — the
    /// reference's FLAG_PROTON_NP_HARD_REMOVED / AC_HARD_* bits,
    /// whose disjunction is FLAG_FORCE_SALT_TAUT: the salt t-group
    /// machinery's later gates relax to donor-only.
    pub hard_moved: bool,
}

/// RemoveNPProtonsAndAcidCharges (simple variants). Returns None when
/// the structure carries no charges (the no-op gate) or nothing
/// changed.
pub fn standardize(structure: &Structure) -> Option<StandardizedProtons> {
    let mut context = ProtonContext::new(structure);
    let initial = totals(&context);
    if initial[attot::NUM_CHARGES] == 0 {
        return None;
    }

    let mut removed = 0u32;
    let mut added = 0u32;
    let mut hard_moved = false;

    // 5.1a: protons off charged N,P and OH(+)
    if initial[attot::NP_PROTON] + initial[attot::OH_PLUS] > 0 {
        removed += simple_remove_hplus_npo(&mut context);
    }

    // 5.1b: the hard removal from cationic N/P (the reference's
    // HardRemoveHplusNP under bCancelChargesAlways=1): a (+) charge
    // group (every N-classed MSK_NP atom, charged or not) rides out
    // against a PR_HARD H-star over the whole structure — the REM_
    // PROTON route moves the proton (and the pi alternation) from the
    // star's donors and cancels the charge at its destination, where
    // AddOrRemoveExplOrImplH(+1) may pin an extra fixed hydrogen (the
    // diazonium ride). Gated on a nonzero NP_Plus total, as in the C.
    // It runs BEFORE the 5.2 acid phase — the reference never re-runs
    // it afterwards.
    let after_51a = totals(&context);
    let mut hard_bonds: Vec<(usize, usize, u32, u8)> = Vec::new();
    if after_51a[attot::NP_PLUS] > 0 {
        // The /p accounting is the reference's net proton counter
        // (tni.nNumRemovedProtons): every H that leaves the structure
        // counts up, every H the ride pins at its destination counts
        // back down — the cancelled (+) itself never does.
        let hydrogens_before: i64 = context.h.iter().map(|&h| h as i64).sum();
        let (rode_out, bonds) = hard_remove_hplus_np(&mut context);
        hard_bonds = bonds;
        let hydrogens_after: i64 = context.h.iter().map(|&h| h as i64).sum();
        hard_moved |= rode_out > 0;
        if hydrogens_before > hydrogens_after {
            removed += (hydrogens_before - hydrogens_after) as u32;
        } else if hydrogens_after > hydrogens_before {
            added += (hydrogens_after - hydrogens_before) as u32;
        }
    }

    // 5.2a: acidic protons off neutral atoms (net +), protons onto
    // anions (net -)
    let after = totals(&context);
    if after[attot::TOT_CHARGE] > 0 {
        removed += simple_remove_acidic_protons(&mut context, after[attot::TOT_CHARGE]);
    } else if after[attot::TOT_CHARGE] < 0 {
        added += simple_add_acidic_protons(&mut context, -after[attot::TOT_CHARGE]);
    }

    // 5.2b: the HARD acidic removal (the reference's
    // HardRemoveAcidicProtons under bTgFlagHardAddRenProtons) — when
    // the simple masks cannot drain the positive total, a proton still
    // leaves: it rides the skeleton from a donor-star member onto the
    // acidic =O/N-O star member (the REM_PROTON route over the two
    // virtual stars), where the driver's follow-up simple pass strips
    // it as a charge. Net effect per move: the donor stays neutral,
    // the acid carries the -1.
    let after = totals(&context);
    let mut moved_bonds = hard_bonds;
    if after[attot::TOT_CHARGE] > 0 {
        let (moved, bonds) = hard_remove_acidic_protons(&mut context, after[attot::TOT_CHARGE]);
        removed += moved;
        moved_bonds.extend(bonds);
        hard_moved |= moved > 0;
    }

    if removed == 0 && added == 0 {
        // The (DE)PROTONATION ran on a charged structure and declined:
        // the charged state IS the standardized state (the reference
        // leaves e.g. an N-carried O(-) at its fixed charge — /q-1
        // with no riding hydrogen). Returning Some pins the structure
        // so the legacy ±1 heuristics never re-open the decision.
        return Some(StandardizedProtons {
            h: context.h,
            charge: context.charge,
            removed_protons: 0,
            added_protons: 0,
            bonds: moved_bonds,
            hard_moved: false,
        });
    }
    Some(StandardizedProtons {
        h: context.h,
        charge: context.charge,
        removed_protons: removed,
        added_protons: added,
        bonds: moved_bonds,
        hard_moved,
    })
}

/// PR_HARD masks (ichi_bns.c): the (+) charge group reads MSK_NP, the
/// H-group and the (-) group read MSK_NP | MSK_OS — every charge-type
/// bit except the onium OH(+)/O(+), the bare proton and the halogens.
const PR_HARD_MSK_NP: u32 =
    atbit::NP_PLUS | atbit::NP_PROTON | atbit::NP_H | atbit::N_MINUS | atbit::NP;
const PR_HARD_MSK_OS: u32 = atbit::COH
    | atbit::CSH
    | atbit::ZOH
    | atbit::OOH
    | atbit::ZOOH
    | atbit::NOH
    | atbit::N_OH
    | atbit::CO
    | atbit::ZO
    | atbit::NO2
    | atbit::N_O
    | atbit::CO_MINUS
    | atbit::CS_MINUS
    | atbit::ZO_MINUS
    | atbit::OO_MINUS
    | atbit::ZOO_MINUS
    | atbit::NO_MINUS
    | atbit::N_O_MINUS;
/// PR_HARD_TYP_H: ATTYP_N | ATTYP_OS — the N class plus the acidic
/// CO/S, OO, ZOO, NO and other-ZO classes (NOT the plain O-minus).
const PR_HARD_TYP_H: u16 =
    att::ATOM_N | att::ACIDIC_CO | att::ACIDIC_S | att::OO | att::ZOO | att::NO | att::OTHER_ZO;

/// HardRemoveHplusNP (TECHMAN 5.1b, bCancelChargesAlways=1): while a
/// REM_PROTON route runs from the PR_HARD H-star to the (+) charge
/// group, each success applies the found alternating paths — the H
/// deltas land on the t-group star members (AddOrRemoveExplOrImplH
/// on both ends: donors shed their proton, the destination may pin a
/// fixed hydrogen), the charge deltas land on the charge-group
/// members (the (+) leaves the cation), and the moved atom-atom bonds
/// keep their new orders. The moved orders are RETURNED (the
/// reference's BNS_EF_SAVE_ALL persists them into the structure): the
/// downstream salt-tautomer machinery reclassifies the H-less
/// endpoint by its raised bond order — an oxime O whose proton rode
/// out reads order 2 and re-enters the salt pairing as an acceptor.
/// Every pass re-reads both stars from the updated H/charge/bond
/// state.
fn hard_remove_hplus_np(context: &mut ProtonContext) -> (u32, Vec<(usize, usize, u32, u8)>) {
    let n = context.element.len();
    let mut rounds = 0;
    let mut all_moved: Vec<(usize, usize, u32, u8)> = Vec::new();
    loop {
        // The (+) charge group: N-classed MSK_NP atoms at charge 0 or
        // +1 (CreateCGroupInBnStruct's OPPOSITE_CHARGE_IN_CGROUP=0
        // gate); the charged member's star edge flows 0, the neutral
        // members' flow 1.
        let plus: Vec<(usize, bool)> = (0..n)
            .filter_map(|atom| {
                let (typ, mask) = atom_charge_type(context, atom);
                if typ & att::ATOM_N == 0 || mask & PR_HARD_MSK_NP == 0 {
                    return None;
                }
                match context.charge[atom] {
                    0 => Some((atom, false)),
                    1 => Some((atom, true)),
                    _ => None,
                }
            })
            .collect();
        // The H-star: the PR_HARD classes under the
        // CreateTGroupInBnStruct gates — the neutral-valence break
        // ends the membership scan, the per-member edge reads
        // cap = neutral valence - bonds (raised for -NH2-) and
        // flow = the capped hydrogen count.
        let mut h_star: Vec<(usize, i32, i32)> = Vec::new();
        for atom in 0..n {
            let (typ, mask) = atom_charge_type(context, atom);
            if typ & PR_HARD_TYP_H == 0 || mask & (PR_HARD_MSK_NP | PR_HARD_MSK_OS) == 0 {
                continue;
            }
            let neutral_valence =
                context.cbv2[atom] / 2 + context.h[atom] as i32 - context.charge[atom];
            if !matches!(neutral_valence, 2 | 3) {
                break;
            }
            let valence = context.neighbors[atom].len() as i32;
            let mut cap = neutral_valence - valence;
            if neutral_valence == 3 && valence > 1 {
                cap += 1;
            }
            let flow = cap.min(context.h[atom] as i32);
            h_star.push((atom, cap, flow));
        }
        if plus.is_empty() || h_star.is_empty() {
            break;
        }
        // The ride network: the skeleton at the standard flow
        // convention, then the three stars in the reference's
        // creation order — (+) group, (-) group, H-group.
        let bonds: Vec<(usize, usize, i32)> = context
            .bonds
            .iter()
            .map(|&(from, to, order2)| {
                let flow = match order2 {
                    4 => 1,
                    6 => 2,
                    _ => 0,
                };
                (from, to, flow)
            })
            .collect();
        let max_flow: Vec<i32> = (0..n)
            .map(|atom| (context.cbv2[atom] / 2 - context.neighbors[atom].len() as i32).max(0))
            .collect();
        // The (-) group: the same PR_HARD classes at charge 0 or -1,
        // without the H-star's valence gates.
        let minus: Vec<(usize, bool)> = (0..n)
            .filter_map(|atom| {
                let (typ, mask) = atom_charge_type(context, atom);
                if typ & PR_HARD_TYP_H == 0 || mask & (PR_HARD_MSK_NP | PR_HARD_MSK_OS) == 0 {
                    return None;
                }
                match context.charge[atom] {
                    0 => Some((atom, false)),
                    -1 => Some((atom, true)),
                    _ => None,
                }
            })
            .collect();
        let (net, donor_star, acceptor_star) =
            crate::canonical::bns::build_hnp(&crate::canonical::bns::HnpSpec {
                bonds,
                max_flow,
                plus,
                minus,
                h_star,
            });
        let ride = crate::canonical::bns::star_ride(&net, donor_star, acceptor_star);
        if std::env::var_os("INCHI_TRACE_PROTON").is_some() {
            eprintln!(
                "HNP round={rounds} route={} tg={:?} cg={:?} bonds={:?}",
                ride.route, ride.tg, ride.cg, ride.bonds
            );
        }
        if !ride.route {
            break;
        }
        // Apply: the H ride onto the t-group members, the charge ride
        // off the charge-group members, and the pi alternation the
        // route moved keeps its new orders.
        for (atom, delta) in &ride.tg {
            let updated = context.h[*atom] as i32 + delta;
            context.h[*atom] = updated.max(0) as u32;
        }
        for (atom, delta) in &ride.cg {
            context.charge[*atom] -= delta;
        }
        for (a, b, delta) in &ride.bonds {
            let Some(entry) = context.bonds.iter_mut().find(|entry| {
                let (from, to) = (entry.0, entry.1);
                (from == *a && to == *b) || (from == *b && to == *a)
            }) else {
                continue;
            };
            // The mark rides on the SOURCE order: a triple demotes to
            // "2 or 3" (the reference's ALTERN|TRIPLE — a bond type
            // the salt role sets never match), anything else to
            // ALT12NS "1 or 2".
            let mark = if entry.2 == 6 { 7u8 } else { 9u8 };
            entry.2 = (entry.2 + 2 * delta).max(2);
            // The ride's pi alternation persists (the reference's
            // BNS_EF_SAVE_ALL): the new orders leave the 5.1b pass
            // with the proton shift, the downstream tautomer machinery
            // reads them.
            let order = (entry.2 / 2) as u32;
            all_moved.push(((*a).min(*b), (*a).max(*b), order, mark));
            context.cbv2[*a] += 2 * delta;
            context.cbv2[*b] += 2 * delta;
        }
        rounds += 1;
    }
    (rounds, all_moved)
}

/// The hard-removal star memberships (HardRemoveAcidicProtons'
/// CreateTGroupInBnStruct scans over the AR_HARD classes): the donor
/// star holds the NP/OS-classed heteroatoms minus the acidic =O/N-O
/// mask bits, the acid star exactly those =C=O / nitro-type oxygens;
/// both under the reference's neutral-valence gate (chem valence + H
/// - charge reads 2 for O/S/Se/Te, 3 for N/P — and a class-passing
///   atom outside the gate ends that scan, as in the C).
pub(crate) fn hard_removal_stars(context: &ProtonContext) -> (Vec<bool>, Vec<bool>) {
    /// AR_HARD_TYP_HN: the N/P and O/S/Se/Te charge-type atoms.
    const TYP_HN: u16 = att::ATOM_N
        | att::ATOM_P
        | att::ACIDIC_CO
        | att::ACIDIC_S
        | att::OO
        | att::ZOO
        | att::NO
        | att::OTHER_ZO;
    /// AR_HARD_MSK_HN: the NP/OS masks minus the acidic =O/N-O bits.
    const MSK_HN: u32 = (atbit::NP_PLUS
        | atbit::NP_PROTON
        | atbit::NP_H
        | atbit::N_MINUS
        | atbit::NP
        | atbit::COH
        | atbit::CSH
        | atbit::ZOH
        | atbit::OOH
        | atbit::ZOOH
        | atbit::NOH
        | atbit::N_OH
        | atbit::ZO
        | atbit::N_O
        | atbit::CO_MINUS
        | atbit::CS_MINUS
        | atbit::ZO_MINUS
        | atbit::OO_MINUS
        | atbit::ZOO_MINUS
        | atbit::NO_MINUS
        | atbit::N_O_MINUS)
        & !(atbit::CO | atbit::NO2);
    /// AR_HARD_TYP_HA / AR_HARD_MSK_HA: the acidic carbonyl / nitro
    /// oxygen the donor proton rides onto (the port's NO2 bit is the
    /// reference's ATBIT_NO, the −N=O slot).
    const TYP_HA: u16 = att::ACIDIC_CO | att::NO;
    const MSK_HA: u32 = atbit::CO | atbit::NO2;

    let n = context.element.len();
    let mut donors = vec![false; n];
    let mut acids = vec![false; n];
    // The neutral-valence gate and the reference's scan-break quirk:
    // per star, the first class-passing atom without a readable
    // neutral valence stops the membership scan.
    let neutral_valence =
        |atom: usize| context.cbv2[atom] / 2 + context.h[atom] as i32 - context.charge[atom];
    for atom in (0..n).map(|atom| (atom, atom_charge_type(context, atom))) {
        let (typ, mask) = atom.1;
        if typ & TYP_HA != 0 && mask & MSK_HA != 0 {
            if !matches!(neutral_valence(atom.0), 2 | 3) {
                break;
            }
            acids[atom.0] = true;
        }
    }
    for atom in (0..n).map(|atom| (atom, atom_charge_type(context, atom))) {
        let (typ, mask) = atom.1;
        if typ & TYP_HN != 0 && mask & MSK_HN != 0 {
            if !matches!(neutral_valence(atom.0), 2 | 3) {
                break;
            }
            donors[atom.0] = true;
        }
    }
    (donors, acids)
}

/// The REM_PROTON route test over a fresh network for this pass: one
/// t-group star over the H-donor heteroatoms, one over the acidic
/// =O/N-O oxygens, both as the reference builds them for the hard
/// removal (AR_HARD classes, neutral-valence gate, member-edge
/// capacity raises). The donor proton must ride the skeleton from its
/// star onto the acidic star for the removal to be viable.
pub(crate) fn rem_proton_route(
    context: &ProtonContext,
    donors: &[bool],
    acids: &[bool],
    holder: usize,
    acceptor: usize,
) -> (bool, Vec<(usize, usize, i32)>) {
    let n = context.element.len();
    let mut bond_adjacency: Vec<Vec<(usize, i32, bool)>> = vec![Vec::new(); n];
    let mut single_like: Vec<Vec<bool>> = vec![Vec::new(); n];
    for &(from, to, order2) in &context.bonds {
        // The tautomer build's flow convention: order minus one
        // (single 0, double 1, triple 2; aromatic reads single).
        let flow = match order2 {
            4 => 1,
            6 => 2,
            _ => 0,
        };
        bond_adjacency[from].push((to, flow, false));
        bond_adjacency[to].push((from, flow, false));
        single_like[from].push(true);
        single_like[to].push(true);
    }
    let max_flow: Vec<i32> = (0..n)
        .map(|atom| (context.cbv2[atom] / 2 - context.neighbors[atom].len() as i32).max(0))
        .collect();
    let group_of: Vec<Option<usize>> = (0..n)
        .map(|atom| {
            if donors[atom] {
                Some(0)
            } else if acids[atom] {
                Some(1)
            } else {
                None
            }
        })
        .collect();
    // A star member pushes only when it actually holds hydrogen (the
    // reference's member-edge flow is the capped H count).
    let donor: Vec<bool> = (0..n)
        .map(|atom| donors[atom] && context.h[atom] > 0)
        .collect();
    let net = crate::canonical::bns::Network::build(&crate::canonical::bns::Build {
        bond_adjacency: &bond_adjacency,
        max_flow: &max_flow,
        group_of: &group_of,
        donor: &donor,
        single_like: &single_like,
        c_point_of: &[None; 0],
        c_charged: &[false; 0],
    });
    crate::canonical::bns::rem_proton_route_moved(&net, holder, acceptor)
}

/// The hard acidic-proton removal (TECHMAN 5.2b): up to `quota`
/// protons ride from donor-star members onto acidic-star members —
/// each move leaves the donor neutral and the acid carrying the -1
/// (the driver's ride + follow-up simple strip, collapsed to its net
/// state). The route's π alternation moves WITH the proton (the
/// reference's BNS_EF_SAVE_ALL | BNS_EF_ALTR_NS ride): the reported
/// atom-atom bonds change order by the signed delta, which both
/// updates this context's valence arithmetic and returns to the
/// caller so the structure can re-read the new orders. Every pass
/// re-reads the stars: the previous move changed the H/charge/bond
/// state the AR_HARD classification reads. Returns (moves, bonds).
fn hard_remove_acidic_protons(
    context: &mut ProtonContext,
    quota: i64,
) -> (u32, Vec<(usize, usize, u32, u8)>) {
    let mut moved = 0;
    let mut all_moved: Vec<(usize, usize, u32, u8)> = Vec::new();
    while (moved as i64) < quota {
        let (donors, acids) = hard_removal_stars(context);
        let acid_sites: Vec<usize> = (0..acids.len()).filter(|&atom| acids[atom]).collect();
        let mut found = None;
        'pair: for holder in 0..context.h.len() {
            if !donors[holder] || context.h[holder] == 0 {
                continue;
            }
            for &acceptor in &acid_sites {
                let (route, delta) = rem_proton_route(context, &donors, &acids, holder, acceptor);
                if route {
                    found = Some((holder, acceptor, delta));
                    break 'pair;
                }
            }
        }
        let Some((holder, acceptor, route_bonds)) = found else {
            break;
        };
        if std::env::var_os("INCHI_TRACE_PROTON").is_some() {
            eprintln!(
                "HARD-REM holder={holder} acceptor={acceptor} moved={moved} bonds={route_bonds:?}"
            );
        }
        context.h[holder] -= 1;
        context.charge[acceptor] -= 1;
        // The π alternation rides with the proton: shift the orders of
        // every atom-atom bond the search moved.
        for (a, b, delta) in route_bonds {
            let Some(entry) = context.bonds.iter_mut().find(|entry| {
                let (from, to) = (entry.0, entry.1);
                (from == a && to == b) || (from == b && to == a)
            }) else {
                continue;
            };
            let mark = if entry.2 == 6 { 7u8 } else { 9u8 };
            let order2 = (entry.2 + 2 * delta).max(2);
            entry.2 = order2;
            let order = (order2 / 2) as u32;
            all_moved.push((a.min(b), a.max(b), order, mark));
            // Keep the doubled valence sums the classifier reads.
            context.cbv2[a] += 2 * delta;
            context.cbv2[b] += 2 * delta;
        }
        moved += 1;
    }
    (moved, all_moved)
}

/// SimpleRemoveHplusNPO: every N/P-H(+) and OH(+) loses its proton.
fn simple_remove_hplus_npo(context: &mut ProtonContext) -> u32 {
    let mut removed = 0;
    for atom in 0..context.h.len() {
        let (typ, mask) = atom_charge_type(context, atom);
        if typ & PR_SIMPLE_TYP != 0 && mask & PR_SIMPLE_MSK != 0 {
            context.h[atom] -= 1;
            context.charge[atom] -= 1;
            removed += 1;
        }
    }
    removed
}

/// Steps shared by the remove/add scans: count eligible atoms per
/// mask step, cut the scan at the step where the quota is reached,
/// then apply in priority order.
fn step_quota(num_per_step: &[i64], quota: i64) -> Option<(i64, usize)> {
    let total: i64 = num_per_step.iter().sum();
    if total == 0 {
        return None;
    }
    let max_step = num_per_step
        .iter()
        .scan(0i64, |acc, &n| {
            *acc += n;
            Some(*acc)
        })
        .position(|cum| cum >= quota)
        .unwrap_or(num_per_step.len() - 1);
    Some((total, max_step))
}

/// SimpleRemoveAcidicProtons: up to `quota` acidic protons off
/// neutral atoms, mask-priority order.
fn simple_remove_acidic_protons(context: &mut ProtonContext, quota: i64) -> u32 {
    let mut num_per_step = vec![0i64; AR_STEPS.len()];
    for atom in 0..context.h.len() {
        if context.charge[atom] != 0 || context.h[atom] == 0 {
            continue;
        }
        let (typ, mask) = atom_charge_type(context, atom);
        for (step, &(step_typ, step_msk)) in AR_STEPS.iter().enumerate() {
            if typ & step_typ != 0 && mask & step_msk != 0 {
                num_per_step[step] += 1;
                break;
            }
        }
    }
    let Some((_, max_step)) = step_quota(&num_per_step, quota) else {
        return 0;
    };

    let mut removed: u32 = 0;
    for atom in 0..context.h.len() {
        if removed as i64 >= quota {
            break;
        }
        if context.charge[atom] != 0 || context.h[atom] == 0 {
            continue;
        }
        let (typ, mask) = atom_charge_type(context, atom);
        for (step, &(step_typ, step_msk)) in AR_STEPS.iter().enumerate() {
            if step > max_step || num_per_step[step] == 0 {
                continue;
            }
            if typ & step_typ != 0 && mask & step_msk != 0 {
                context.h[atom] -= 1;
                context.charge[atom] -= 1;
                num_per_step[step] -= 1;
                removed += 1;
                break;
            }
        }
    }
    removed
}

/// SimpleAddAcidicProtons: up to `quota` protons onto anionic
/// heteroatoms, mask-priority order.
fn simple_add_acidic_protons(context: &mut ProtonContext, quota: i64) -> u32 {
    let mut num_per_step = vec![0i64; AA_STEPS.len()];
    for atom in 0..context.h.len() {
        if context.charge[atom] != -1 {
            continue;
        }
        let (typ, mask) = atom_charge_type(context, atom);
        for (step, &(step_typ, step_msk)) in AA_STEPS.iter().enumerate() {
            if typ & step_typ != 0 && mask & step_msk != 0 {
                num_per_step[step] += 1;
                break;
            }
        }
    }
    let Some((_, max_step)) = step_quota(&num_per_step, quota) else {
        return 0;
    };

    let mut added: u32 = 0;
    for atom in 0..context.h.len() {
        if added as i64 >= quota {
            break;
        }
        if context.charge[atom] != -1 {
            continue;
        }
        let (typ, mask) = atom_charge_type(context, atom);
        for (step, &(step_typ, step_msk)) in AA_STEPS.iter().enumerate() {
            if step > max_step || num_per_step[step] == 0 {
                continue;
            }
            if typ & step_typ != 0 && mask & step_msk != 0 {
                context.h[atom] += 1;
                context.charge[atom] += 1;
                num_per_step[step] -= 1;
                added += 1;
                break;
            }
        }
    }
    added
}

#[cfg(test)]
mod tests {
    use super::{hard_removal_stars, standardize, ProtonContext};
    use crate::graph::{AtomRecord, Structure};

    fn structure(atoms: Vec<(&str, i32, u32)>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure {
            atoms: atoms
                .into_iter()
                .map(|(element, charge, h_count)| AtomRecord {
                    element: element.to_string(),
                    h_count,
                    charge,
                    h_pinned: false,
                    isotope: 0,
                })
                .collect(),
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
    fn the_neutral_structure_is_a_no_op() {
        let s = structure(vec![("C", 0, 0), ("O", 0, 0)], vec![(0, 1, 1)]);
        assert_eq!(standardize(&s), None);
    }

    #[test]
    fn the_anion_gains_its_proton() {
        // Acetate: the standardized skeleton is the acid; the added
        // proton rides /p as -1.
        let s = structure(
            vec![("C", 0, 3), ("C", 0, 0), ("O", 0, 0), ("O", -1, 0)],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let std = standardize(&s).expect("acetate standardizes");
        assert_eq!(std.added_protons, 1);
        assert_eq!(std.removed_protons, 0);
        assert_eq!(std.charge[3], 0);
        assert_eq!(std.h[3], 1);
        assert!(std.h.iter().all(|_| true));
    }

    #[test]
    fn the_cationic_acid_loses_its_acidic_proton() {
        // Ammonium acetate (covalent reading): the +1 on nitrogen
        // drives an acidic proton off the carboxyl — the oracle's b53
        // mechanism on a minimal body.
        let s = structure(
            vec![
                ("N", 1, 3),
                ("C", 0, 0),
                ("C", 0, 0),
                ("O", 0, 0),
                ("O", 0, 0),
            ],
            vec![(0, 1, 1), (1, 2, 1), (2, 3, 2), (2, 4, 1)],
        );
        let std = standardize(&s).expect("cationic acid standardizes");
        assert_eq!(std.removed_protons, 1);
        // 5.1a removed the N-H(+) proton; the acid proton stays.
        assert_eq!(std.charge[0], 0);
        assert_eq!(std.h[0], 2);
    }

    #[test]
    fn ammonium_loses_the_proton_on_the_charged_nitrogen() {
        let s = structure(vec![("N", 1, 4)], vec![]);
        let std = standardize(&s).expect("ammonium standardizes");
        assert_eq!(std.removed_protons, 1);
        assert_eq!(std.h[0], 3);
        assert_eq!(std.charge[0], 0);
    }

    #[test]
    fn the_dianion_gains_two_protons() {
        // Sulfate: the two single-bonded oxygens are the anionic
        // sites; one proton each.
        let s = structure(
            vec![
                ("S", 0, 0),
                ("O", 0, 0),
                ("O", 0, 0),
                ("O", -1, 0),
                ("O", -1, 0),
            ],
            vec![(0, 1, 2), (0, 2, 2), (0, 3, 1), (0, 4, 1)],
        );
        let std = standardize(&s).expect("sulfate standardizes");
        assert_eq!(std.added_protons, 2);
        assert!(std.charge.iter().all(|&charge| charge == 0));
    }

    #[test]
    fn the_hard_removal_stars_split_donor_from_acidic_oxygen() {
        // The amide/cation shape (fixture identity/breadth-580): the
        // donor star holds the N-H amide nitrogen, the acid star the
        // carbonyl oxygen — the two virtual endpoints of the
        // REM_PROTON route.
        let s = structure(
            vec![
                ("C", 0, 0),
                ("N", 0, 0),
                ("C", 0, 0),
                ("O", 0, 0),
                ("C", 0, 0),
                ("N", 1, 0),
            ],
            vec![(0, 1, 1), (1, 2, 1), (2, 3, 2), (2, 4, 1), (4, 5, 1)],
        );
        let context = ProtonContext::new(&s);
        let (donors, acids) = hard_removal_stars(&context);
        assert!(donors[1], "the amide N-H joins the donor star");
        assert!(acids[3], "the =O joins the acid star");
        assert!(!donors[3] && !acids[1], "the pair sits in disjoint stars");
        assert!(donors.iter().zip(&acids).all(|(&d, &a)| !(d && a)));
    }
}
