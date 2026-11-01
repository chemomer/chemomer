//! The V2000 molfile reader — the identity pipeline's input edge.
//! Encapsulated on purpose: callers see a `Structure`, never the
//! column arithmetic; only the subset the corpus identity molfiles
//! use (the gem writer's output) is supported, and unsupported
//! constructs are reported rather than guessed.

use crate::error::{err, Result};
use crate::graph::{AtomRecord, Structure};

const CONTEXT: &str = "molfile";

/// Reads the V2000 subset: counts line, atom block, bond block, and
/// `M  CHG` / `M  ISO` property lines. Isotope mass numbers land on the
/// atoms (the /i layer and the hydrogen fold read them).
pub fn parse(input: &str) -> Result<Structure> {
    let lines: Vec<&str> = input.lines().collect();
    if lines.len() < 5 {
        return err(format!("{CONTEXT}: molfile too short"));
    }

    let counts = lines[3];
    let atom_count: usize = counts
        .get(0..3)
        .unwrap_or("")
        .trim()
        .parse()
        .map_err(|_| bad_counts(counts))?;
    let bond_count: usize = counts
        .get(3..6)
        .unwrap_or("")
        .trim()
        .parse()
        .map_err(|_| bad_counts(counts))?;
    // The chiral flag (counts column 12–15): nonzero marks the
    // structure chiral (the reference's FLAG_INP_AT_CHIRAL, read at
    // mol2atom.c:483) — the SUCF option and the aux /rA letter read
    // it.
    let chiral_flag: u32 = counts.get(12..15).unwrap_or("").trim().parse().unwrap_or(0);

    let mut atoms = Vec::with_capacity(atom_count);
    let mut positions = Vec::with_capacity(atom_count);
    for index in 0..atom_count {
        let line = line_or(&lines, 4 + index, "atom block")?;
        let element = line.get(31..34).unwrap_or("").trim().to_string();
        if element.is_empty() {
            return err(format!("{CONTEXT}: atom {index} has no element symbol"));
        }
        let code: u32 = line.get(42..45).unwrap_or("").trim().parse().unwrap_or(0);
        // The MDL charge code (cols 37-40): the oracle SDfile flavor
        // carries charges here (3 = +1, 5 = −1, …); the gem's
        // molfiles keep it zero and use M CHG instead, which wins
        // when both appear.
        let charge_code: u32 = line.get(37..40).unwrap_or("").trim().parse().unwrap_or(0);
        let charge = match charge_code {
            0 => 0,
            1 => 3,
            2 => 2,
            3 => 1,
            5 => -1,
            6 => -2,
            7 => -3,
            4 => {
                return err(format!(
                    "{CONTEXT}: atom {index} radical charge code unsupported"
                ));
            }
            _ => {
                return err(format!(
                    "{CONTEXT}: atom {index} unknown charge code {charge_code}"
                ));
            }
        };
        let x: f64 = line.get(0..10).unwrap_or("").trim().parse().unwrap_or(0.0);
        let y: f64 = line.get(10..20).unwrap_or("").trim().parse().unwrap_or(0.0);
        positions.push((x, y));
        atoms.push(AtomRecord {
            element,
            h_count: code.saturating_sub(1),
            charge,
            h_pinned: false,
            isotope: 0,
        });
    }

    let mut bonds = Vec::with_capacity(bond_count);
    let mut stereo_bonds = std::collections::BTreeMap::new();
    for index in 0..bond_count {
        let line = line_or(&lines, 4 + atom_count + index, "bond block")?;
        let from: usize = line
            .get(0..3)
            .unwrap_or("")
            .trim()
            .parse()
            .map_err(|_| bad(line))?;
        let to: usize = line
            .get(3..6)
            .unwrap_or("")
            .trim()
            .parse()
            .map_err(|_| bad(line))?;
        let kind: u32 = line
            .get(6..9)
            .unwrap_or("")
            .trim()
            .parse()
            .map_err(|_| bad(line))?;
        // The stereo column may be absent (9-char bond lines — the
        // corpus's merged molfiles): absent reads as none.
        let stereo: u8 = line.get(9..12).unwrap_or("0").trim().parse().unwrap_or(0);
        if stereo != 0 {
            stereo_bonds.insert(index, stereo);
        }
        bonds.push((from - 1, to - 1, kind));
    }

    // The alias block: `A  n` lines replace atom n's symbol with an
    // enriched one — the oracle SDfile's isotope flavor (OD2 = oxygen
    // carrying two deuteriums). The element keeps its parsed symbol's
    // leading letters; trailing D/T counts fold into the isotope map
    // (the /h layer stays isotope-blind, so the carrier's count is
    // the valence gap).
    let mut isotopic_h: std::collections::BTreeMap<usize, (u32, u32)> = Default::default();
    let property_lines = &lines[4 + atom_count + bond_count..];
    for (position, line) in property_lines.iter().enumerate() {
        if !line.starts_with("A  ") {
            continue;
        }
        let index: usize = line[3..].trim().parse().map_err(|_| bad(line))?;
        let symbol = property_lines
            .get(position + 1)
            .ok_or_else(|| bad(line))?
            .trim();
        let atom = atoms.get_mut(index - 1).ok_or_else(|| bad(line))?;
        let bytes = symbol.as_bytes();
        let mut end = 1;
        if bytes.get(1).is_some_and(|b| b.is_ascii_lowercase()) {
            end = 2;
        }
        if bytes.len() < end {
            return err(format!("{CONTEXT}: alias {symbol:?} has no element"));
        }
        atom.element = symbol[..end].to_string();
        let mut deuteriums = 0u32;
        let mut tritiums = 0u32;
        let mut cursor = end;
        while cursor < bytes.len() {
            match bytes[cursor] {
                b'D' | b'T' => {
                    let is_deuterium = bytes[cursor] == b'D';
                    cursor += 1;
                    let start = cursor;
                    while cursor < bytes.len() && bytes[cursor].is_ascii_digit() {
                        cursor += 1;
                    }
                    let count: u32 = if cursor > start {
                        symbol[start..cursor].parse().map_err(|_| bad(line))?
                    } else {
                        1
                    };
                    if is_deuterium {
                        deuteriums += count;
                    } else {
                        tritiums += count;
                    }
                }
                _ => {
                    return err(format!("{CONTEXT}: alias {symbol:?} marker unsupported"));
                }
            }
        }
        if deuteriums > 0 || tritiums > 0 {
            isotopic_h.insert(index - 1, (deuteriums, tritiums));
        }
    }

    for line in property_lines {
        let is_charge = line.starts_with("M  CHG");
        if !(is_charge || line.starts_with("M  ISO")) {
            continue;
        }
        let values: Vec<&str> = line[6..].split_whitespace().collect();
        for pair in values[1..].chunks(2) {
            if pair.len() != 2 {
                continue;
            }
            let index: usize = pair[0].parse().map_err(|_| bad(line))?;
            let value: i32 = pair[1].parse().map_err(|_| bad(line))?;
            let atom = atoms.get_mut(index - 1).ok_or_else(|| bad(line))?;
            if is_charge {
                atom.charge = value;
            } else {
                atom.isotope = value.max(0) as u16;
            }
        }
    }

    Ok(Structure {
        atoms,
        bonds,
        stereo_bonds,
        positions,
        isotopic_h,
        protons_hard_moved: false,
        ride_marks: Vec::new(),
        fold_source: Vec::new(),
        chiral_flag: chiral_flag != 0,
    })
}

fn line_or<'a>(lines: &'a [&'a str], index: usize, block: &str) -> Result<&'a str> {
    lines
        .get(index)
        .copied()
        .ok_or_else(|| Error::new(format!("{CONTEXT}: truncated {block}")))
}

use crate::error::Error;

fn bad_counts(line: &str) -> Error {
    Error::new(format!("{CONTEXT}: malformed counts line: {line:?}"))
}

fn bad(line: &str) -> Error {
    Error::new(format!("{CONTEXT}: malformed line: {line:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    #[test]
    fn reads_elements_hcounts_and_charges() {
        let molfile = "ammonium\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 N   0  0  0  5  0  0  0  0  0  0  0\nM  CHG  1   1   1\nM  END\n";
        let structure = parse(molfile).expect("parses");
        assert_eq!(structure.atoms[0].element, "N");
        // MDL code 5 means 4 explicit hydrogens.
        assert_eq!(structure.atoms[0].h_count, 4);
        assert_eq!(structure.atoms[0].charge, 1);
    }

    #[test]
    fn reads_bonds_as_zero_based_types() {
        let molfile = "ethane\n  AsciiChem\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0  0  0  0\nM  END\n";
        let structure = parse(molfile).expect("parses");
        assert_eq!(structure.bonds, vec![(0, 1, 1)]);
    }

    #[test]
    fn truncated_molfile_is_an_error_not_a_guess() {
        assert!(parse("only a header\n").is_err());
    }

    #[test]
    fn water_fixture_constant_parses() {
        let structure = parse(WATER).expect("parses");
        assert_eq!(structure.atoms[0].element, "O");
    }
}

#[cfg(test)]
mod oracle_flavor {
    use super::*;

    const HYDROXIDE: &str = "Structure #1. \n  InChIV10                                      \n\n  1  0  0  0  0  0  0  0  0  2 V2000\n    0.0000    0.0000    0.0000 O   0  5  0     0  0  0  0  0  0\nM  CHG  1   1  -1\nM  END\n$$$$\n";

    const D2O: &str = "Structure #1. \n  InChIV10                                      \n\n  1  0  0  0  0  0  0  0  0  3 V2000\n    0.0000    0.0000    0.0000 O   0  0  0     0  0  0  0  0  0\nA  1\nOD2\nM  END\n$$$$\n";

    const GEM_WATER: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    #[test]
    fn the_oracle_charge_code_column_reads_as_a_charge() {
        // Hydroxide's charge rides the atom-line code (5 = −1); the
        // M CHG line repeats it and wins.
        let structure = parse(HYDROXIDE).expect("parses");
        assert_eq!(structure.atoms[0].charge, -1);
    }

    #[test]
    fn the_oracle_alias_line_folds_isotopic_hydrogens() {
        let structure = parse(D2O).expect("parses");
        assert_eq!(structure.atoms[0].element, "O");
        assert_eq!(structure.isotopic_h.get(&0), Some(&(2, 0)));
    }

    #[test]
    fn the_gem_flavor_still_parses_without_the_new_columns() {
        let structure = parse(GEM_WATER).expect("parses");
        assert_eq!(structure.atoms[0].charge, 0);
        assert!(structure.isotopic_h.is_empty());
    }
}
