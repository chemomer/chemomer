//! The V2000/V3000 molfile reader — the identity pipeline's input
//! edge. Encapsulated on purpose: callers see a `Structure`, never
//! the column arithmetic; only the subset the corpus identity
//! molfiles use (the gem writer's output) is supported, and
//! unsupported constructs are reported rather than guessed.

use crate::error::{err, Result};
use crate::graph::{AtomRecord, Structure};

const CONTEXT: &str = "molfile";

/// Reads the V2000 subset (counts line, atom block, bond block, and
/// `M  CHG` / `M  ISO` property lines) or the V3000 CTAB block —
/// dispatched on the counts line's version tag, the way the reference
/// reads the version string after the eleven fixed fields and arms
/// its V3000 reader on "V3000" (mol_fmt1.c:630-645). Isotope mass
/// numbers land on the atoms (the /i layer and the hydrogen fold read
/// them).
///
/// This is the DEFAULT read (the reference's `treat_polymers` off,
/// mol_fmt1.c:1365-1372): the whole V2000 Sgroup property family
/// (`M  STY/SST/SLB/SCN/SAL/SBL/SDI/SMT/SBT` — superatom
/// abbreviations included) is skipped wholesale, the structure reads
/// as plain atoms + bonds, and the only trace is the reference's
/// "Ignore polymer data" warning — no polymer data reaches the
/// structure, so no sgroup line can fail it either.
pub fn parse(input: &str) -> Result<Structure> {
    read(input, false)
}

/// [`parse`] under the `-Polymers` option (the reference's
/// `treat_polymers` on): the Sgroup property family parses into the
/// polymer units, and the superatom/abbreviation types (SUP — the
/// drawn-expanded group whose label rides `M  SMT`; TEMPLATE in
/// modern writers) fail the read the way the reference's reader does,
/// with every failing line's raw text accumulated into the message
/// (mol_fmt1.c MolfileReadSgroupOfPolymer). The V3000 flavor carries
/// the equivalent as SGROUP blocks — those stay discarded in EVERY
/// mode (the reference's V3000 reader never converts them into
/// polymer units; probe-pinned: a V3000 SUP reads plain under
/// `-Polymers` too).
pub fn parse_with_polymers(input: &str) -> Result<Structure> {
    read(input, true)
}

fn read(input: &str, treat_polymers: bool) -> Result<Structure> {
    let lines: Vec<&str> = input.lines().collect();
    if lines.len() < 5 {
        return err(format!("{CONTEXT}: molfile too short"));
    }

    let counts = lines[3];
    if counts_version(counts) == Some("V3000") {
        return parse_v3000(&lines[4..]);
    }
    parse_v2000(input, treat_polymers)
}

/// The counts line's CTFile version tag — the whitespace token at
/// column 33, after the eleven 3-character fields
/// (mol_fmt1.c:630-636 reads exactly that field and dispatches).
fn counts_version(counts: &str) -> Option<&str> {
    counts.get(33..)?.split_whitespace().last()
}

/// Splits an SDfile into its records' verbatim texts, on the `$$$$`
/// terminator lines — the batch reader's whole face. Each record keeps
/// everything up to its terminator (header, connection table, any SDF
/// data fields), so a record parses exactly as it would standing
/// alone; empty gaps between terminators produce no records.
pub fn parse_sdf(input: &str) -> Vec<String> {
    let mut records = Vec::new();
    let mut current = String::new();
    for line in input.lines() {
        if line.trim_end() == "$$$$" {
            if !current.trim().is_empty() {
                records.push(std::mem::take(&mut current));
            }
        } else {
            current.push_str(line);
            current.push('\n');
        }
    }
    if !current.trim().is_empty() {
        records.push(current);
    }
    records
}

/// Reads the V2000 positional fields from the PARG capture tree: the
/// counts line, atom block, bond block, and the property-line tail
/// are the grammar's ordered line captures; the charge-code table,
/// the tolerant field defaults, and the polymer semantics stay here.
/// `treat_polymers` decides the Sgroup family's fate: skipped
/// wholesale (the default read) or parsed into the polymer units
/// (the `-Polymers` read).
fn parse_v2000(input: &str, treat_polymers: bool) -> Result<Structure> {
    crate::parg::with_molfile_capture(input, |capture| {
        let mol = capture
            .child(capture.root(), "molfile")
            .ok_or_else(|| Error::new(format!("{CONTEXT}: molfile capture missing")))?;
        // The counts ride the column-named granules; every read goes
        // through the same get(a..b) slices the reader always used,
        // so a short line's tolerant defaults land identically.
        let counts = capture
            .array_child(&mol, "counts")
            .ok_or_else(|| Error::new(format!("{CONTEXT}: counts capture missing")))?;
        let counts_segments = capture.segments(&counts);
        let atom_count: usize = counts_segments
            .column(0, 3)
            .unwrap_or("")
            .trim()
            .parse()
            .map_err(|_| bad_counts(counts_segments.raw()))?;
        let bond_count: usize = counts_segments
            .column(3, 6)
            .unwrap_or("")
            .trim()
            .parse()
            .map_err(|_| bad_counts(counts_segments.raw()))?;
        // The chiral flag (counts column 12–15): nonzero marks the
        // structure chiral (the reference's FLAG_INP_AT_CHIRAL, read
        // at mol2atom.c:483) — the SUCF option and the aux /rA letter
        // read it.
        let chiral_flag: u32 = counts_segments
            .column(12, 15)
            .unwrap_or("")
            .trim()
            .parse()
            .unwrap_or(0);

        // Each body unit descends its capture wrappers (the shape
        // rule nests one or more levels deep) to its positional
        // granules. A unit without granules is a blank line: an empty
        // segment stream (its column reads return None, exactly what
        // slicing an empty tail line did). Only MISSING units — the
        // counts naming more lines than the file carries — stay None,
        // the truncated-block error the direct slice raised.
        let body = capture
            .array_child(&mol, "body")
            .map(|body| capture.array_items(&body))
            .unwrap_or_default();
        // The repetition's array carries an empty-marker element
        // (the match's zero-length seed) before the iterations; the
        // iterations are arrays themselves, blank lines included.
        let mut lines: Vec<Option<crate::parg::Segments>> = Vec::with_capacity(body.len());
        for unit in &body {
            if capture.is_empty_marker(unit) {
                continue;
            }
            lines.push(Some(capture.segments(unit)));
        }
        if lines.len() < atom_count + bond_count {
            lines.resize(atom_count + bond_count, None);
        }

        let mut atoms = Vec::with_capacity(atom_count);
        let mut positions = Vec::with_capacity(atom_count);
        for (index, line) in lines.iter_mut().enumerate().take(atom_count) {
            let segments = line
                .as_ref()
                .ok_or_else(|| Error::new(format!("{CONTEXT}: truncated atom block")))?;
            let element = segments.column(31, 34).unwrap_or("").trim().to_string();
            if element.is_empty() {
                return err(format!("{CONTEXT}: atom {index} has no element symbol"));
            }
            // The MDL charge code (cols 37-40): the oracle SDfile flavor
            // carries charges here (3 = +1, 5 = −1, …); the gem's
            // molfiles keep it zero and use M CHG instead, which wins
            // when both appear.
            let charge_code: u32 = segments
                .column(37, 40)
                .unwrap_or("")
                .trim()
                .parse()
                .unwrap_or(0);
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
            let x: f64 = segments
                .column(0, 10)
                .unwrap_or("")
                .trim()
                .parse()
                .unwrap_or(0.0);
            let y: f64 = segments
                .column(10, 20)
                .unwrap_or("")
                .trim()
                .parse()
                .unwrap_or(0.0);
            positions.push((x, y));
            atoms.push(AtomRecord {
                // The MDL star atom reads as the pseudo-element Zz at
                // parse (mol_fmt1.c MolfileTreatPseudoElementAtoms:
                // `*` ===> "Zz" when pseudo-elements are allowed); the
                // engine's option gate decides permission.
                element: if element == "*" {
                    "Zz".to_string()
                } else {
                    element
                },
                // The H-count code field is NEVER read by the release
                // engine: mol_fmt1.c walks it only under MOL_FMT_QUERY,
                // and mol_fmt.h:79 sets MOL_FMT_QUERY = MOL_FMT_ABSENT —
                // inp_num_H reaches get_num_H (util.c:862) as zero
                // (probe-proven: a bare N carrying H-count 2 still reads
                // H3N, oracle recording valence-codes/hcount-ignored).
                h_count: 0,
                // The valence code column, by contrast, IS read in the
                // release path: the same MolfileReadField walk reads it
                // unconditionally (mol_fmt1.c:728), calculate_valences
                // hands it to get_num_H (mol2atom.c:963) and codes 1-15
                // re-derive the implicit H (bonded N with code 1 reads
                // CH3N, bare N with 15 reads bare N, a singly-bonded C
                // with 15 reads H14 — oracle-recorded probes; the parked
                // sodium-citrate quaternary C with code 2 reads zero H
                // either way, max(0, 2−4) = the gap, the one case that
                // masked the read when this unit was parked). Column
                // arithmetic is the walk's cumulative offsets: element
                // 31..34, mass-difference 34..36, charge 36..39, parity
                // 39..42, the skipped query pair 42..48, valence 48..51.
                // Malformed or out-of-S_CHAR-range bytes store 0 there
                // (util.c MolfileReadField) — the same tolerant default.
                valence_code: segments
                    .column(48, 51)
                    .unwrap_or("")
                    .trim()
                    .parse::<i32>()
                    .map_or(0, |v| i16::try_from(v).unwrap_or(0)),
                charge,
                h_pinned: false,
                isotope: 0,
            });
        }

        let mut bonds = Vec::with_capacity(bond_count);
        let mut stereo_bonds = std::collections::BTreeMap::new();
        for index in 0..bond_count {
            let segments = lines[atom_count + index]
                .as_ref()
                .ok_or_else(|| Error::new(format!("{CONTEXT}: truncated bond block")))?;
            let raw = segments.raw();
            let from: usize = segments
                .column(0, 3)
                .unwrap_or("")
                .trim()
                .parse()
                .map_err(|_| bad(raw))?;
            let to: usize = segments
                .column(3, 6)
                .unwrap_or("")
                .trim()
                .parse()
                .map_err(|_| bad(raw))?;
            let kind: u32 = segments
                .column(6, 9)
                .unwrap_or("")
                .trim()
                .parse()
                .map_err(|_| bad(raw))?;
            // The stereo column may be absent (9-char bond lines — the
            // corpus's merged molfiles): absent reads as none.
            let stereo: u8 = segments
                .column(9, 12)
                .unwrap_or("")
                .trim()
                .parse()
                .unwrap_or(0);
            if stereo != 0 {
                stereo_bonds.insert(index, stereo);
            }
            bonds.push((from - 1, to - 1, kind));
        }

        // The property tail: every body line past the blocks, as raw
        // text. The alias / property / sgroup scans walk it the way
        // they walked the tail's line slice.
        let tail: Vec<&str> = lines[atom_count + bond_count..]
            .iter()
            .map(|segments| segments.as_ref().map_or("", crate::parg::Segments::raw))
            .collect();
        let isotopic_h = apply_alias_lines(&tail, &mut atoms)?;
        apply_property_lines(&tail, &mut atoms)?;
        let polymer = if treat_polymers {
            parse_sgroup_lines(&tail, atom_count, bond_count)?
        } else {
            None
        };

        Ok(Structure {
            atoms,
            bonds,
            stereo_bonds,
            positions,
            isotopic_h,
            protons_hard_moved: false,
            ride_marks: Vec::new(),
            taut_ext: Default::default(),
            fold_source: Vec::new(),
            chiral_flag: chiral_flag != 0,
            polymer,
        })
    })
}

/// Reads the V3000 CTAB block: `M  V30`-prefixed logical lines after
/// the quasi-counts line — `BEGIN CTAB`, the `COUNTS` line, the atom
/// and bond blocks, the optional tail blocks (SGROUP, OBJ3D, LINKNODE,
/// COLLECTION — read and skipped, the reference stores enhanced stereo
/// only to ignore it), and `END CTAB` (mol_fmt3.c). Graph semantics —
/// element case, bond orders 1-3 plus aromatic 4, the V2000 stereo
/// codes, charges and isotope masses — are the reference's: free-form
/// whitespace fields, `KEY=value` attributes in the V30 form of the
/// V2000 `M  CHG` / `M  ISO` property lines (mol_fmt3.c:1071-1122).
fn parse_v3000(lines: &[&str]) -> Result<Structure> {
    let mut v30 = V30Lines { lines, pos: 0 };

    expect_v30_line(&mut v30, "BEGIN CTAB")?;

    let counts_line = next_v30_line(&mut v30)?
        .ok_or_else(|| Error::new(format!("{CONTEXT} V3000: truncated before the COUNTS line")))?;
    let fields: Vec<&str> = counts_line.split_whitespace().collect();
    if fields.first() != Some(&"COUNTS") {
        return err(format!(
            "{CONTEXT} V3000: cannot interpret the COUNTS line: {counts_line:?}"
        ));
    }
    let atom_count: usize = fields
        .get(1)
        .and_then(|field| field.parse().ok())
        .ok_or_else(|| Error::new(format!("{CONTEXT} V3000: COUNTS carries no atom count")))?;
    let bond_count: usize = fields
        .get(2)
        .and_then(|field| field.parse().ok())
        .ok_or_else(|| Error::new(format!("{CONTEXT} V3000: COUNTS carries no bond count")))?;
    // The chiral flag is the COUNTS line's fifth number (mol_fmt3.c:530).
    let chiral_flag: u32 = fields
        .get(5)
        .and_then(|field| field.parse().ok())
        .unwrap_or(0);

    expect_v30_line(&mut v30, "BEGIN ATOM")?;
    let mut atoms = Vec::with_capacity(atom_count);
    let mut positions = Vec::with_capacity(atom_count);
    // Original atom index -> structure position: V3000 bonds may name
    // the atoms in any order (mol_fmt3.c get_actual_atom_number).
    let mut index_map = std::collections::HashMap::with_capacity(atom_count);
    for position in 0..atom_count {
        let line = next_v30_line(&mut v30)?
            .ok_or_else(|| Error::new(format!("{CONTEXT} V3000: truncated atom block")))?;
        let (index, atom, position_xy) = parse_v3000_atom(&line)?;
        index_map.insert(index, position);
        positions.push(position_xy);
        atoms.push(atom);
    }
    expect_v30_line(&mut v30, "END ATOM")?;

    let mut bonds = Vec::with_capacity(bond_count);
    let mut stereo_bonds = std::collections::BTreeMap::new();
    // The reference skips the whole bond block when the counts say
    // zero bonds — no BEGIN BOND line need exist (mol_fmt3.c:1274-1277).
    if bond_count > 0 {
        expect_v30_line(&mut v30, "BEGIN BOND")?;
        for bond_index in 0..bond_count {
            let line = next_v30_line(&mut v30)?
                .ok_or_else(|| Error::new(format!("{CONTEXT} V3000: truncated bond block")))?;
            let (bond, stereo) = parse_v3000_bond(&line, &index_map)?;
            if let Some(stereo) = stereo {
                stereo_bonds.insert(bond_index, stereo);
            }
            bonds.push(bond);
        }
        expect_v30_line(&mut v30, "END BOND")?;
    }

    // The tail: optional SGROUP, OBJ3D, LINKNODE and COLLECTION
    // blocks — read and skipped (the reference reads the enhanced
    // stereo collections only to ignore them, mol_fmt3.c:1602-1727).
    loop {
        let line = next_v30_line(&mut v30)?
            .ok_or_else(|| Error::new(format!("{CONTEXT} V3000: no END CTAB marker")))?;
        match line.as_str() {
            "BEGIN SGROUP" => skip_v30_block(&mut v30, "END SGROUP")?,
            "BEGIN OBJ3D" => skip_v30_block(&mut v30, "END OBJ3D")?,
            "LINKNODE" => {}
            "BEGIN COLLECTION" => skip_v30_block(&mut v30, "END COLLECTION")?,
            "END CTAB" => break,
            other => {
                return err(format!(
                    "{CONTEXT} V3000: cannot interpret the CTAB line {other:?}"
                ));
            }
        }
    }

    // After the CTAB the physical remainder is the property-line tail:
    // the reference still walks `M  CHG` / `M  ISO` there (the property
    // loop runs while the version string is non-empty, mol_fmt1.c:995).
    apply_property_lines(&v30.lines[v30.pos..], &mut atoms)?;

    Ok(Structure {
        atoms,
        bonds,
        stereo_bonds,
        positions,
        isotopic_h: Default::default(),
        protons_hard_moved: false,
        ride_marks: Vec::new(),
        taut_ext: Default::default(),
        fold_source: Vec::new(),
        chiral_flag: chiral_flag != 0,
        polymer: None,
    })
}

/// One `M  V30 <index> <element> x y z <aamap> [KEY=value ...]` atom
/// line: positional fields then attributes (mol_fmt3.c:956-1222).
/// `CHG=` and `MASS=` carry the V30 form of the V2000 `M  CHG` /
/// `M  ISO` property lines — charge and absolute isotope mass number.
fn parse_v3000_atom(line: &str) -> Result<(usize, AtomRecord, (f64, f64))> {
    let mut tokens = line.split_whitespace();
    let index: usize = v30_number(&mut tokens, line)?;
    let symbol = tokens.next().ok_or_else(|| bad(line))?;
    let x: f64 = v30_number(&mut tokens, line)?;
    let y: f64 = v30_number(&mut tokens, line)?;
    let _z: f64 = v30_number(&mut tokens, line)?;
    let _aamap: usize = v30_number(&mut tokens, line)?;

    let mut atom = AtomRecord {
        // Two-letter symbols arrive case-normalized (mol_fmt3.c:1054-1057
        // lowercases a trailing uppercase letter: "CL" reads "Cl");
        // the star atom reads Zz here too.
        element: if symbol == "*" {
            "Zz".to_string()
        } else {
            normalize_v3000_symbol(symbol)
        },
        h_count: 0,
        charge: 0,
        // The VAL= attribute loop below fills this (V3000 spells zero
        // valence -1, normalized to the V2000 code 15 there).
        valence_code: 0,
        h_pinned: false,
        isotope: 0,
    };
    for attribute in tokens {
        let Some((key, value)) = attribute.split_once('=') else {
            return err(format!(
                "{CONTEXT} V3000: atom attribute without a value: {attribute:?}"
            ));
        };
        match key {
            "CHG" => atom.charge = v30_value(value, line)?,
            "MASS" => atom.isotope = v30_value(value, line)?,
            // The V3000 valence attribute feeds the same
            // get_num_H branch as the V2000 column (mol_fmt3.c:1134-
            // 1145 reads VAL= into ctab.atoms[].valence): -1 is the
            // V3000 spelling of zero valence and normalizes to the
            // V2000 code 15, every other value takes the code's seat
            // through the C `char` cast (VAL=-2 probe reads zero H —
            // the truthy-negative gap clamps to zero, exactly what
            // the i16 field carries).
            "VAL" => {
                let raw: i32 = v30_value(value, line)?;
                atom.valence_code = if raw == -1 { 15 } else { (raw as i8) as i16 };
            }
            // Read-and-ignored attributes: the reference consumes the
            // value and moves on (query and reaction hints never
            // reach the identity pipeline — the V2000 H-count column
            // probe (hcount ignored) covers HCOUNT's seat too).
            "HCOUNT" | "RAD" | "CFG" | "STBOX" | "INVRET" | "EXACHG" | "SUBST" | "UNSAT"
            | "RBCNT" | "ATTCHPT" | "RGROUPS" | "ATTCHORD" | "CLASS" | "SEQID" => {}
            other => {
                return err(format!("{CONTEXT} V3000: unknown atom attribute {other:?}"));
            }
        }
    }
    Ok((index, atom, (x, y)))
}

/// One `M  V30 <index> <type> <atom1> <atom2> [KEY=value ...]` bond
/// line. Types are the V2000 orders (1-3 plain, 4 aromatic);
/// `CFG=` carries the V30 form of the V2000 stereo column — 1 up,
/// 2 either, 3 down, remapped to the V2000 codes 1/4/6
/// (mol_fmt3.c:1407-1428).
fn parse_v3000_bond(
    line: &str,
    index_map: &std::collections::HashMap<usize, usize>,
) -> Result<((usize, usize, u32), Option<u8>)> {
    let mut tokens = line.split_whitespace();
    let _index: usize = v30_number(&mut tokens, line)?;
    let bond_type: u32 = v30_number(&mut tokens, line)?;
    if !(1..=4).contains(&bond_type) {
        return err(format!(
            "{CONTEXT} V3000: unsupported bond type {bond_type}"
        ));
    }
    let from: usize = v30_number(&mut tokens, line)?;
    let to: usize = v30_number(&mut tokens, line)?;
    let mut stereo = None;
    for attribute in tokens {
        let Some((key, value)) = attribute.split_once('=') else {
            return err(format!(
                "{CONTEXT} V3000: bond attribute without a value: {attribute:?}"
            ));
        };
        match key {
            "CFG" => {
                let code: u8 = v30_value(value, line)?;
                stereo = Some(match code {
                    1 => 1,
                    2 => 4,
                    3 => 6,
                    other => {
                        return err(format!(
                            "{CONTEXT} V3000: unsupported bond CFG code {other}"
                        ));
                    }
                });
            }
            // Read-and-ignored attributes (query, display, reaction).
            "TOPO" | "RXCTR" | "STBOX" | "DISP" | "ATTACH" => {}
            // Haptic bonds (ENDPTS) are stored-and-ignored by the
            // reference — a construct the identity pipeline never
            // sees; reported rather than guessed here.
            "ENDPTS" => {
                return err(format!("{CONTEXT} V3000: haptic bonds are not supported"));
            }
            other => {
                return err(format!("{CONTEXT} V3000: unknown bond attribute {other:?}"));
            }
        }
    }
    let from = *index_map.get(&from).ok_or_else(|| bad(line))?;
    let to = *index_map.get(&to).ok_or_else(|| bad(line))?;
    Ok(((from, to, bond_type), stereo))
}

/// Two-letter symbols arrive case-normalized: a trailing uppercase
/// letter lowercases ("CL" reads "Cl") — mol_fmt3.c:1054-1057.
fn normalize_v3000_symbol(symbol: &str) -> String {
    let mut chars: Vec<char> = symbol.chars().collect();
    if chars.len() == 2 && chars[1].is_ascii_uppercase() {
        chars[1] = chars[1].to_ascii_lowercase();
    }
    chars.into_iter().collect()
}

/// The `M  V30 ` logical-line reader: every physical line carries the
/// 7-character prefix (stripped here), and a trailing `-` continues
/// onto the next physical line (mol_fmt3.c
/// get_V3000_input_line_to_strbuf:1889-1927).
struct V30Lines<'a> {
    lines: &'a [&'a str],
    pos: usize,
}

fn next_v30_line(v30: &mut V30Lines<'_>) -> Result<Option<String>> {
    let Some(line) = v30.lines.get(v30.pos).copied() else {
        return Ok(None);
    };
    v30.pos += 1;
    let Some(stripped) = line.strip_prefix("M  V30 ") else {
        return err(format!(
            "{CONTEXT} V3000: V30 line lacks the `M  V30 ` prefix: {line:?}"
        ));
    };
    let mut logical = stripped.trim_end().to_string();
    while logical.ends_with('-') {
        logical.pop();
        let Some(continuation) = v30.lines.get(v30.pos).copied() else {
            return err(format!(
                "{CONTEXT} V3000: continued line has no continuation: {line:?}"
            ));
        };
        v30.pos += 1;
        let Some(stripped) = continuation.strip_prefix("M  V30 ") else {
            return err(format!(
                "{CONTEXT} V3000: V30 line lacks the `M  V30 ` prefix: {continuation:?}"
            ));
        };
        logical.push_str(stripped.trim_end());
    }
    Ok(Some(logical))
}

fn expect_v30_line(v30: &mut V30Lines<'_>, marker: &str) -> Result<()> {
    match next_v30_line(v30)? {
        Some(line) if line == marker => Ok(()),
        Some(line) => err(format!(
            "{CONTEXT} V3000: expected {marker:?}, found {line:?}"
        )),
        None => err(format!("{CONTEXT} V3000: truncated before {marker:?}")),
    }
}

fn skip_v30_block(v30: &mut V30Lines<'_>, marker: &str) -> Result<()> {
    loop {
        match next_v30_line(v30)? {
            Some(line) if line == marker => return Ok(()),
            Some(_) => {}
            None => {
                return err(format!("{CONTEXT} V3000: no {marker:?} marker"));
            }
        }
    }
}

fn v30_number<'a, T: std::str::FromStr>(
    tokens: &mut impl Iterator<Item = &'a str>,
    line: &str,
) -> Result<T> {
    let field = tokens.next().ok_or_else(|| bad(line))?;
    field.parse().map_err(|_| bad(line))
}

fn v30_value<T: std::str::FromStr>(value: &str, line: &str) -> Result<T> {
    value.parse().map_err(|_| bad(line))
}

/// The alias block: `A  n` lines replace atom n's symbol with an
/// enriched one — the oracle SDfile's isotope flavor (OD2 = oxygen
/// carrying two deuteriums). The element keeps its parsed symbol's
/// leading letters; trailing D/T counts fold into the isotope map
/// (the /h layer stays isotope-blind, so the carrier's count is
/// the valence gap).
fn apply_alias_lines(
    lines: &[&str],
    atoms: &mut [AtomRecord],
) -> Result<std::collections::BTreeMap<usize, (u32, u32)>> {
    let mut isotopic_h: std::collections::BTreeMap<usize, (u32, u32)> = Default::default();
    for (position, line) in lines.iter().enumerate() {
        if !line.starts_with("A  ") {
            continue;
        }
        let index: usize = line[3..].trim().parse().map_err(|_| bad(line))?;
        let symbol = lines.get(position + 1).ok_or_else(|| bad(line))?.trim();
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
    Ok(isotopic_h)
}

/// The `M  CHG` / `M  ISO` property-line tail shared by both flavors.
fn apply_property_lines(lines: &[&str], atoms: &mut [AtomRecord]) -> Result<()> {
    for line in lines {
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
    Ok(())
}

/// The Sgroup property lines under the `-Polymers` read
/// (mol_fmt1.c MolfileReadSgroupOfPolymer, 1451–1900): `M  STY`
/// declares units (SRU/MON/COP/MOD/CRO/MER), `M  SST` subtypes
/// (ALT/RAN/BLK), `M  SCN` connection patterns (HT/HH/EU), `M  SLB`
/// unique-id labels (read; the /Z: echo never prints them), `M  SAL`
/// atom lists, `M  SBL` crossing-bond lists, `M  SMT` subscripts,
/// `M  SDI` bracket coordinates (first pair of lines per unit). The
/// lines may arrive in any order — SST/SCN/SLB/SAL/SBL/SDI/SMT
/// reference units by Sgroup number, and a reference to an undeclared
/// unit fails the line exactly like the reference's
/// MolFmtSgroups_GetIndexBySgroupId miss.
///
/// The superatom/abbreviation types (SUP — TEMPLATE in modern
/// writers — and every other CTFile type word: GRA, MUL, DAT, COM,
/// MIX, FOR, ANY, GEN) are NOT polymer types: the reference's STY
/// reader maps only the six and fails the line otherwise
/// (mol_fmt1.c:1498-1503). A failure frees the WHOLE sgroups table
/// (`MolFmtSgroups_Free` memsets it zero, mol_fmt4.c:751-765) — and
/// the freed table never comes back (`MolFmtSgroups_Append`'s ReAlloc
/// refuses a zeroed struct), so every later Sgroup line fails too and
/// its raw text joins the accumulated message; the read rejects at
/// the end (the executable's "no InChI" error over `Could not
/// interpret Molfile polymer data:` + the failing lines, joined "; ",
/// duplicate line texts dropped — ichierr.c AddErrorMessage).
///
/// A file with no Sgroup lines reads `Ok(None)`.
fn parse_sgroup_lines(
    lines: &[&str],
    atom_count: usize,
    bond_count: usize,
) -> Result<Option<crate::graph::Polymer>> {
    let mut units: Vec<crate::graph::PolymerUnit> = Vec::new();
    // SDI lines arrive pairwise per unit: the first fills xbr1, the
    // second xbr2 (mol_fmt1.c SDI handler's brkt1-not-yet-here test).
    let mut sdi_seen: Vec<usize> = Vec::new();
    // The failure state: the table died (this and every later Sgroup
    // line fails) and the failing lines so far.
    let mut table_dead = false;
    let mut failed_lines: Vec<String> = Vec::new();
    let find_unit = |units: &[crate::graph::PolymerUnit], id: usize| {
        units.iter().position(|unit| unit.id == id)
    };
    for line in lines {
        if !line.starts_with("M  S") {
            continue;
        }
        let kind = &line[3..6.min(line.len())];
        let fields: Vec<&str> = line[6..].split_whitespace().collect();
        // A failed line: its raw text joins the message and the table
        // dies — units already collected from EARLIER lines included
        // (the reference frees them all at err_exit).
        macro_rules! fail_line {
            () => {{
                table_dead = true;
                units.clear();
                sdi_seen.clear();
                failed_lines.push((*line).to_string());
            }};
        }
        // The dead-table walk. The reference still evaluates each
        // tag's own shape on the way: STY/SST/SCN/SLB read their count
        // in the dispatch chain — a line whose count cannot read is
        // silently ignored (the chain condition fails and nothing else
        // matches); SAL/SBL/SDI/SMT run handler bodies that fail on
        // ANY dead-table reference, unreadable fields included. Our
        // whitespace-split parse draws the same line for every
        // well-formed spelling, which is what the corpus carries.
        match kind {
            "STY" => {
                // `M  STY nn8 (sg, type)...` — the count field then
                // (number, 3-letter-type) pairs, at most 8 per line.
                if table_dead {
                    if let Some(count) = dead_count(&fields) {
                        if fields.len() > 2 * count {
                            fail_line!();
                        }
                    }
                    continue;
                }
                // The count rides the dispatch chain's width-3 field:
                // a line whose count cannot read (or whose pairs are
                // missing) never enters the branch — silently ignored
                // (mol_fmt1.c:1460-1463), the shape behind the oracle's
                // treatment of the malformed spellings.
                let Some(count) = fields.first().and_then(|f| f.parse::<usize>().ok()) else {
                    continue;
                };
                if count == 0 || count > 8 || fields.len() < 1 + 2 * count {
                    continue;
                }
                let mut rejected = false;
                for pair in fields[1..1 + 2 * count].chunks(2) {
                    let id = id_of(pair[0], line)?;
                    let unit_type = match pair[1] {
                        "SRU" => 1u8,
                        "MON" => 2,
                        "COP" => 3,
                        "MOD" => 4,
                        "CRO" => 5,
                        "MER" => 6,
                        // SUP and every other abbreviation type word:
                        // the failure that kills the whole table — the
                        // pairs the line declared BEFORE the bad one
                        // die with it.
                        _ => {
                            rejected = true;
                            break;
                        }
                    };
                    match find_unit(&units, id) {
                        Some(index) => units[index].unit_type = unit_type,
                        None => units.push(crate::graph::PolymerUnit {
                            id,
                            unit_type,
                            ..Default::default()
                        }),
                    }
                }
                if rejected {
                    fail_line!();
                }
            }
            "SST" => {
                if table_dead {
                    if let Some(count) = dead_count(&fields) {
                        if fields.len() > 2 * count {
                            fail_line!();
                        }
                    }
                    continue;
                }
                let Some(count) = fields.first().and_then(|f| f.parse::<usize>().ok()) else {
                    continue;
                };
                if count == 0 || count > 8 || fields.len() < 1 + 2 * count {
                    continue;
                }
                let mut rejected = false;
                for pair in fields[1..1 + 2 * count].chunks(2) {
                    let subtype = match pair[1] {
                        "ALT" => 1u8,
                        "RAN" => 2,
                        "BLK" | "BLO" => 3,
                        _ => {
                            rejected = true;
                            break;
                        }
                    };
                    match find_unit(&units, id_of(pair[0], line)?) {
                        Some(index) => units[index].subtype = subtype,
                        None => {
                            rejected = true;
                            break;
                        }
                    }
                }
                if rejected {
                    fail_line!();
                }
            }
            "SCN" => {
                if table_dead {
                    if let Some(count) = dead_count(&fields) {
                        if fields.len() > 2 * count {
                            fail_line!();
                        }
                    }
                    continue;
                }
                let Some(count) = fields.first().and_then(|f| f.parse::<usize>().ok()) else {
                    continue;
                };
                if count == 0 || count > 8 || fields.len() < 1 + 2 * count {
                    continue;
                }
                let mut rejected = false;
                for pair in fields[1..1 + 2 * count].chunks(2) {
                    let conn = match pair[1] {
                        "HT" => 1u8,
                        "HH" => 2,
                        "EU" => 3,
                        _ => {
                            rejected = true;
                            break;
                        }
                    };
                    match find_unit(&units, id_of(pair[0], line)?) {
                        Some(index) => units[index].conn = conn,
                        None => {
                            rejected = true;
                            break;
                        }
                    }
                }
                if rejected {
                    fail_line!();
                }
            }
            "SLB" => {
                // `M  SLB nn8 (sg, label)...` — the label lands in the
                // reference's unit struct but the /Z: echo never prints
                // it; only the dead/undeclared reference matters here.
                if table_dead {
                    if let Some(count) = dead_count(&fields) {
                        if fields.len() > 2 * count {
                            fail_line!();
                        }
                    }
                    continue;
                }
                let Some(count) = fields.first().and_then(|f| f.parse::<usize>().ok()) else {
                    continue;
                };
                if count == 0 || count > 8 || fields.len() < 1 + 2 * count {
                    continue;
                }
                let mut rejected = false;
                for pair in fields[1..1 + 2 * count].chunks(2) {
                    let declared = find_unit(&units, id_of(pair[0], line)?).is_some();
                    let label = pair.get(1).and_then(|value| value.parse::<i32>().ok());
                    if !declared || label.is_none() {
                        rejected = true;
                        break;
                    }
                }
                if rejected {
                    fail_line!();
                }
            }
            "SAL" => {
                // The handler-body shape: an unreadable id/count/pair
                // fails the line like any dead reference (the reference's
                // `else { fail = 1; }` paths, mol_fmt1.c:1705-1713).
                if table_dead {
                    fail_line!();
                    continue;
                }
                let Some(id) = fields.first().and_then(|f| f.parse::<usize>().ok()) else {
                    fail_line!();
                    continue;
                };
                let Some(count) = fields.get(1).and_then(|f| f.parse::<usize>().ok()) else {
                    fail_line!();
                    continue;
                };
                if fields.len() < 2 + count {
                    fail_line!();
                    continue;
                }
                let index = match find_unit(&units, id) {
                    Some(index) => index,
                    None => {
                        fail_line!();
                        continue;
                    }
                };
                let mut rejected = false;
                for field in &fields[2..2 + count] {
                    let Ok(atom) = field.parse::<usize>() else {
                        rejected = true;
                        break;
                    };
                    if atom == 0 || atom > atom_count {
                        rejected = true;
                        break;
                    }
                    units[index].alist.push(atom - 1);
                }
                if rejected {
                    fail_line!();
                }
            }
            "SBL" => {
                if table_dead {
                    fail_line!();
                    continue;
                }
                let Some(id) = fields.first().and_then(|f| f.parse::<usize>().ok()) else {
                    fail_line!();
                    continue;
                };
                let Some(count) = fields.get(1).and_then(|f| f.parse::<usize>().ok()) else {
                    fail_line!();
                    continue;
                };
                if fields.len() < 2 + count {
                    fail_line!();
                    continue;
                }
                let index = match find_unit(&units, id) {
                    Some(index) => index,
                    None => {
                        fail_line!();
                        continue;
                    }
                };
                let mut rejected = false;
                for field in &fields[2..2 + count] {
                    let Ok(bond) = field.parse::<usize>() else {
                        rejected = true;
                        break;
                    };
                    if bond == 0 || bond > bond_count {
                        rejected = true;
                        break;
                    }
                    units[index].blist.push(bond - 1);
                }
                if rejected {
                    fail_line!();
                }
            }
            "SMT" => {
                if table_dead {
                    fail_line!();
                    continue;
                }
                let Some(id) = fields.first().and_then(|f| f.parse::<usize>().ok()) else {
                    fail_line!();
                    continue;
                };
                if fields.len() < 2 {
                    fail_line!();
                    continue;
                }
                let label = fields[1..].join(" ");
                match find_unit(&units, id) {
                    Some(index) => units[index].smt = label,
                    None => fail_line!(),
                }
            }
            "SDI" => {
                if table_dead {
                    fail_line!();
                    continue;
                }
                let Some(id) = fields.first().and_then(|f| f.parse::<usize>().ok()) else {
                    fail_line!();
                    continue;
                };
                if fields.len() < 6 || fields.get(1) != Some(&"4") {
                    fail_line!();
                    continue;
                }
                let mut coords = [0.0f64; 4];
                let mut rejected = false;
                for (slot, field) in coords.iter_mut().zip(fields[2..6].iter()) {
                    match field.parse::<f64>() {
                        Ok(value) => *slot = value,
                        Err(_) => {
                            rejected = true;
                            break;
                        }
                    }
                }
                if rejected {
                    fail_line!();
                    continue;
                }
                let index = match find_unit(&units, id) {
                    Some(index) => index,
                    None => {
                        fail_line!();
                        continue;
                    }
                };
                let seen = sdi_seen.iter().filter(|&&s| s == id).count();
                if seen == 0 {
                    units[index].xbr1 = Some(coords);
                } else {
                    units[index].xbr2 = Some(coords);
                }
                sdi_seen.push(id);
            }
            _ => {}
        }
    }
    if !failed_lines.is_empty() {
        return err(cascade_message(&failed_lines));
    }
    // Units with no STY line never materialize (SAL/SBL-only files are
    // rejected by the reference's GetIndexBySgroupId miss, surfaced
    // above); a declaration pass with zero units means no polymer.
    Ok((!units.is_empty()).then_some(crate::graph::Polymer { units }))
}

/// The STY/SST/SCN/SLB count read over a dead table: the dispatch
/// chain's width-3 window; `None` reads as the silent-ignore shape.
fn dead_count(fields: &[&str]) -> Option<usize> {
    fields
        .first()
        .and_then(|field| field.parse::<usize>().ok())
        .filter(|count| (1..=8).contains(count))
}

fn id_of(field: &str, line: &str) -> Result<usize> {
    field.parse().map_err(|_| bad(line))
}

/// The reference's accumulated read-error message (ichierr.c
/// AddErrorMessage): the prefix once, then every failing line — the
/// line after the ':'-terminated prefix joins with a single space,
/// the rest with "; ", and a line already present (the
/// already_have_this_message boundary check) is dropped.
fn cascade_message(lines: &[String]) -> String {
    let mut message = String::from("Could not interpret Molfile polymer data:");
    for line in lines {
        if already_have(&message, line) {
            continue;
        }
        if message.ends_with(':') {
            message.push(' ');
        } else {
            message.push_str("; ");
        }
        message.push_str(line);
    }
    message
}

/// ichierr.c already_have_this_message: the message contains `line`
/// as a whole separator-bounded member.
fn already_have(message: &str, line: &str) -> bool {
    let mut search = 0;
    while let Some(found) = message[search..].find(line) {
        let start = search + found;
        let end = start + line.len();
        let bounded_before = start == 0
            || (start >= 2
                && (&message[start - 2..start] == "; " || &message[start - 2..start] == ": "));
        let bounded_after = end == message.len() || message[end..].starts_with("; ");
        if bounded_before && bounded_after {
            return true;
        }
        search = start + 1;
    }
    false
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

    /// Propane drawn with its CH3 end marked as a superatom
    /// abbreviation (STY SUP + SAL + SMT) — the shape the RInChI
    /// distribution's R005a_with_agents.rxn agents carry.
    const PROPANE_SUP: &str = "propane ch3 sup\n  chemomer\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   24.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   36.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  STY  1   1 SUP\nM  SAL   1  2   1   2\nM  SMT   1 CH3\nM  END\n";

    /// The superatom is metadata, not a construct: the default read
    /// (polymers off) skips the whole Sgroup family — the structure
    /// reads as plain atoms + bonds, the reference's "Ignore polymer
    /// data" behavior (mol_fmt1.c:1365-1372; oracle probe over the
    /// recorded sup/propane-abbrev tranche fixture).
    #[test]
    fn debug_slb() {
        let molfile = "sru with slb\n  chemomer\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   24.0000   24.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   36.0000   36.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  STY  1   1 SRU\nM  SAL   1  2   1   2\nM  SBL   1  1   1\nM  SLB  1   9   7\nM  END\n";
        let _ = crate::parg::with_molfile_capture(molfile, |capture| {
            let mol = capture.child(capture.root(), "molfile").expect("mol");
            for (i, item) in capture.array_items(&mol).iter().enumerate() {
                println!("[{i}] {item:?}");
                for (key, inner) in capture.hash_items(item) {
                    println!("    {key} -> {inner:?}");
                }
            }
            Ok(())
        });
    }

    #[test]
    fn the_default_read_skips_superatom_lines() {
        let structure = parse(PROPANE_SUP).expect("parses");
        assert_eq!(structure.atoms.len(), 3);
        assert_eq!(structure.bonds, vec![(0, 1, 1), (1, 2, 1)]);
        assert!(structure.polymer.is_none());
    }

    /// The -Polymers read fails the superatom the way the reference's
    /// reader does: the STY line fails (mol_fmt1.c:1498), the table
    /// frees, and the SAL/SMT lines referencing the freed unit fail
    /// after it — the accumulated "Could not interpret Molfile polymer
    /// data:" message carries every failing line verbatim
    /// (ichierr.c AddErrorMessage; oracle-recorded message).
    #[test]
    fn the_polymers_read_rejects_superatoms_with_the_accumulated_lines() {
        let error = parse_with_polymers(PROPANE_SUP).expect_err("rejects");
        assert_eq!(
            error.message(),
            "Could not interpret Molfile polymer data: M  STY  1   1 SUP; \
             M  SAL   1  2   1   2; M  SMT   1 CH3"
        );
    }

    /// The freed table never comes back (MolFmtSgroups_Free memsets it
    /// zero and Append's ReAlloc refuses the zeroed struct): a later
    /// VALID SRU declaration fails too, and the read still rejects —
    /// the sup/sup-before-sru tranche fixture's recording.
    #[test]
    fn a_valid_polymer_after_a_superatom_still_rejects() {
        let molfile = "sup then sru\n  chemomer\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   24.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   36.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   48.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\nM  STY  1   1 SUP\nM  SAL   1  2   1   2\nM  SMT   1 MeO\nM  STY  1   2 SRU\nM  SAL   2  2   3   4\nM  SBL   2  1   3\nM  SCN   2 HT\nM  END\n";
        let error = parse_with_polymers(molfile).expect_err("rejects");
        assert_eq!(
            error.message(),
            "Could not interpret Molfile polymer data: M  STY  1   1 SUP; \
             M  SAL   1  2   1   2; M  SMT   1 MeO; M  STY  1   2 SRU; \
             M  SAL   2  2   3   4; M  SBL   2  1   3"
        );
    }

    /// The polymer family itself still reads under -Polymers — the
    /// superatom rejection must not have disturbed it.
    #[test]
    fn the_polymers_read_still_parses_polymer_units() {
        let molfile = "sru chain\n  chemomer\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  STY  1   1 SRU\nM  SAL   1  2   1   2\nM  END\n";
        let structure = parse_with_polymers(molfile).expect("parses");
        let polymer = structure.polymer.expect("units");
        assert_eq!(polymer.units.len(), 1);
        assert_eq!(polymer.units[0].unit_type, 1);
        assert_eq!(polymer.units[0].alist, vec![0, 1]);
    }

    /// An undeclared SLB reference fails under -Polymers (the
    /// sup/slb-dead-reference tranche fixture) and stays ignored under
    /// the default read.
    #[test]
    fn an_undeclared_slb_reference_fails_under_polymers() {
        let molfile = "sru with slb\n  chemomer\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   24.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   36.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  STY  1   1 SRU\nM  SAL   1  2   1   2\nM  SBL   1  1   1\nM  SLB  1   9   7\nM  END\n";
        assert!(parse(molfile).is_ok());
        let error = parse_with_polymers(molfile).expect_err("rejects");
        assert_eq!(
            error.message(),
            "Could not interpret Molfile polymer data: M  SLB  1   9   7"
        );
    }

    #[test]
    fn reads_elements_hcounts_and_charges() {
        let molfile = "ammonium\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 N   0  0  0  5  0  0  0  0  0  0  0\nM  CHG  1   1   1\nM  END\n";
        let structure = parse(molfile).expect("parses");
        assert_eq!(structure.atoms[0].element, "N");
        // The H-count code field is NEVER read (mol_fmt.h:79
        // MOL_FMT_QUERY = MOL_FMT_ABSENT — the reader jumps the
        // field): the implicit H derive at materialization from the
        // charge-shifted valence rows instead. Code 5 on a bare N+
        // derives the same 4 hydrogens there.
        assert_eq!(structure.atoms[0].h_count, 0);
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

    // The V3000 flavor: the same water, one CTAB block.
    const WATER_V3000: &str = "water\n  AsciiChem\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 O 12 12 0 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n";

    #[test]
    fn reads_the_v3000_block_into_the_same_structure() {
        let v3000 = parse(WATER_V3000).expect("parses");
        let v2000 = parse(WATER).expect("parses");
        assert_eq!(v3000, v2000);
    }

    #[test]
    fn v3000_carries_charges_isotopes_and_case_normalizes_symbols() {
        let molfile = "dichloro\n  AsciiChem\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 CL 12 12 0 0 MASS=37\nM  V30 2 N 20 12 0 0 CHG=1\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
        let structure = parse(molfile).expect("parses");
        assert_eq!(structure.atoms[0].element, "Cl");
        assert_eq!(structure.atoms[0].isotope, 37);
        assert_eq!(structure.atoms[1].element, "N");
        assert_eq!(structure.atoms[1].charge, 1);
        assert_eq!(structure.bonds, vec![(0, 1, 1)]);
    }

    #[test]
    fn v3000_aromatic_bonds_and_wedge_cfg_map_to_the_v2000_codes() {
        let molfile = "wedge\n  AsciiChem\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 3 0 0 1\nM  V30 BEGIN ATOM\nM  V30 1 C 12 12 0 0\nM  V30 2 C 20 14 0 0\nM  V30 3 O 4 14 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 CFG=1\nM  V30 2 4 2 3\nM  V30 3 1 1 3 CFG=3\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
        let structure = parse(molfile).expect("parses");
        assert_eq!(structure.bonds, vec![(0, 1, 1), (1, 2, 4), (0, 2, 1)]);
        assert_eq!(structure.stereo_bonds.get(&0), Some(&1));
        assert_eq!(structure.stereo_bonds.get(&2), Some(&6));
        assert!(structure.chiral_flag);
    }

    #[test]
    fn v3000_continuation_lines_join_into_one_logical_line() {
        let molfile = "wrapped\n  AsciiChem\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 12 12 0 0\nM  V30 2 O 20 12 0 -\nM  V30 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
        let structure = parse(molfile).expect("parses");
        assert_eq!(structure.atoms[1].element, "O");
        assert_eq!(structure.atoms[1].isotope, 0);
        assert_eq!(structure.bonds, vec![(0, 1, 1)]);
    }

    #[test]
    fn v3000_tail_blocks_are_skipped_to_end_ctab() {
        let molfile = "collected\n  AsciiChem\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 12 12 0 0\nM  V30 END ATOM\nM  V30 BEGIN COLLECTION\nM  V30 MDLV30/STEABS ATOMS=(1 1)\nM  V30 END COLLECTION\nM  V30 END CTAB\nM  END\n";
        let structure = parse(molfile).expect("parses");
        assert_eq!(structure.atoms[0].element, "C");
        assert!(!structure.chiral_flag);
    }

    #[test]
    fn v3000_unsupported_constructs_report_rather_than_guess() {
        let unknown = "odd\n  AsciiChem\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 12 12 0 0 FOO=3\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n";
        assert!(parse(unknown).is_err(), "unknown atom attribute");
        let haptic = "haptic\n  AsciiChem\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 Fe 12 12 0 0\nM  V30 2 C 20 12 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 ENDPTS=(5 2 3 4 5 6) ATTACH=ALL\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
        assert!(parse(haptic).is_err(), "haptic bond");
        let unprefixed =
            "bad\n  AsciiChem\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nBEGIN CTAB\nM  END\n";
        assert!(parse(unprefixed).is_err(), "missing M  V30 prefix");
    }

    #[test]
    fn v3000_zero_bond_counts_need_no_bond_block() {
        let molfile = "lone\n  AsciiChem\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 Na 12 12 0 0 CHG=1\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n";
        let structure = parse(molfile).expect("parses");
        assert!(structure.bonds.is_empty());
        assert_eq!(structure.atoms[0].charge, 1);
    }

    #[test]
    fn v3000_elements_read_through_the_same_engine() {
        let engine = crate::engine::Engine::from_molfile(WATER_V3000).expect("engine");
        let inchi = engine.inchi_string().expect("inchi");
        assert_eq!(inchi, "InChI=1S/H2O/h1H2");
    }

    #[test]
    fn parse_sdf_splits_on_the_terminator_lines() {
        let sdf = format!("{WATER}$$$$\n{WATER}$$$$\n");
        let records = parse_sdf(&sdf);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0], WATER.to_string());
        assert_eq!(records[1], WATER.to_string());
    }

    #[test]
    fn sdf_data_fields_ride_their_record() {
        let with_fields = format!("{WATER}> <ID>\nsdf/water\n\n$$$$\n");
        let records = parse_sdf(&with_fields);
        assert_eq!(records.len(), 1);
        assert!(records[0].contains("> <ID>"));
        // The record still parses: the data fields trail the
        // connection table, where the property-line loop ignores them.
        let engine = crate::engine::Engine::from_molfile(&records[0]).expect("engine");
        assert_eq!(engine.inchi_string().expect("inchi"), "InChI=1S/H2O/h1H2");
    }

    #[test]
    fn sdf_empty_gaps_and_a_missing_final_terminator_read() {
        let sdf = format!("{WATER}$$$$\n\n\n$$$$\n{WATER}");
        let records = parse_sdf(&sdf);
        assert_eq!(records.len(), 2, "the empty gap produces no record");
        assert!(parse_sdf("").is_empty());
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
