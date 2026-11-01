//! InChI-string → structure reconstruction: the reverse direction of
//! the identity pipeline.
//!
//! The wire grammars mirror the reference reader (INCHI_BASE/src/
//! ichiread.c, dev-time reading only — never vendored):
//! [`ParseSegmentFormula`] (ichiread.c:9810) for the Hill-order atom
//! list, [`ParseSegmentConnections`] (ichiread.c:9084) for the /c
//! branch-stack walk, [`ParseSegmentMobileH`] (ichiread.c:8230) for
//! the /h fixed counts and `(Hn,atoms)` t-groups,
//! [`ParseSegmentCharge`] (ichiread.c:8032), [`ParseSegmentProtons`]
//! (ichiread.c:7426), [`ParseSegmentIsoAtoms`] (ichiread.c:5806) and
//! [`ParseSegmentIsoExchgH`] (ichiread.c:5541) for /q /p /i.
//!
//! What the reference leaves to its BNS structure restore
//! (ichirvr1..7.c, the deliberately lossy mobile-H/charge placement)
//! this module reconstructs by observable rule — explicit H counts,
//! charges and a valence-feasible Kekulé assignment — and the gate
//! (tests/gate.rs) pins those choices against the oracle's own
//! recorded reconstruction: our forward of our reconstruction must
//! equal the oracle's re-forward of the oracle's reconstruction.

use crate::error::{err, Result};
use crate::graph::{AtomRecord, Structure};

const CONTEXT: &str = "reverse";

/// The standard isotope masses the /i layer's ±deltas shift from —
/// the same table the isotope layer emits against.
fn standard_mass(element: &str) -> Option<u16> {
    Some(match element {
        "H" => 1,
        "B" => 11,
        "C" => 12,
        "N" => 14,
        "O" => 16,
        "F" => 19,
        "P" => 31,
        "S" => 32,
        "Cl" => 35,
        "Br" => 79,
        "I" => 127,
        "Se" => 80,
        "Te" => 128,
        "Si" => 28,
        _ => return None,
    })
}

/// The allowed valences per element and charge — the reference's
/// charged-valence rows (ichirvr1.c:118-469: anions read V−1,
/// group-13/14 cations read V−1, onium group-15/16/17 cations read
/// V+1), with the hypervalent chalcogen/pnictogen sets the oracle's
/// reconstructions realize (nitrate's N(V), sulfate's S(VI)).
fn charged_valences(element: &str, charge: i32) -> Option<Vec<f64>> {
    let mut base: Vec<f64> = match element {
        "H" => vec![1.0],
        "B" => vec![3.0],
        "C" => vec![4.0],
        "N" => vec![3.0, 5.0],
        "O" => vec![2.0],
        "F" | "Cl" | "Br" | "I" => vec![1.0],
        "Si" => vec![4.0],
        "P" => vec![3.0, 5.0],
        "S" | "Se" | "Te" => vec![2.0, 4.0, 6.0],
        _ => return None,
    };
    if charge == 0 {
        return Some(base);
    }
    // Anions shift every valence down; cations shift by the onium
    // direction: group-15/16/17 gain a bond site, group-13/14 (and H)
    // lose one.
    let onium = matches!(
        element,
        "N" | "P" | "O" | "S" | "Se" | "Te" | "F" | "Cl" | "Br" | "I"
    );
    let delta = if charge < 0 || onium { charge } else { -charge };
    for value in &mut base {
        *value += delta as f64;
    }
    base.retain(|value| *value >= 0.0);
    if base.is_empty() {
        None
    } else {
        Some(base)
    }
}

/// A parsed `(Hn,atoms)` mobile-H group (ichiread.c's t-group record).
#[derive(Debug, Clone)]
struct MobileGroup {
    hydrogens: u32,
    endpoints: Vec<usize>,
}

/// A parsed /i entry: per-atom mass delta and/or isotopic-H counts.
#[derive(Debug, Clone, Default)]
struct IsoEntry {
    atom: usize,
    delta: Option<i32>,
    deuteriums: u32,
    tritiums: u32,
}

/// The parsed layers of one single-component standard InChI.
#[derive(Debug, Clone, Default)]
struct Layers {
    /// Heavy elements in Hill order, expanded per atom (canonical
    /// numbering 1..n; H is implicit and never numbered).
    elements: Vec<String>,
    formula_hydrogens: u32,
    bonds: Vec<(usize, usize)>,
    /// Fixed /h counts per canonical atom (0-based index here).
    fixed_h: Vec<u32>,
    groups: Vec<MobileGroup>,
    charge: i32,
    protons: i32,
    isotopes: Vec<IsoEntry>,
    /// The /i segment's exchangeable-H counts (the `hD2` form after
    /// an empty /i): (deuteriums, tritiums).
    exchanged: (u32, u32),
}

/// Reconstructs a molecular [`Structure`] from a standard InChI
/// string — the reverse direction of the identity pipeline. Scope:
/// single-component standard InChIs with formula / /c / /h (fixed and
/// mobile groups) / /q / /p / /i layers; compressed (abc) layers,
/// Fixed-H (`/f`) layers, multi-component strings and the stereo
/// layers (`/b` `/t` `/m` `/s`) are refused with a named error.
pub fn structure_from_inchi(inchi: &str) -> Result<Structure> {
    let layers = parse_layers(inchi)?;
    reconstruct(layers)
}

/// Parses the identifier into its layers (ichiread.c's segment
/// dispatch over the '/'-separated body).
fn parse_layers(inchi: &str) -> Result<Layers> {
    let rest = inchi.strip_prefix("InChI=1S/").ok_or_else(|| {
        crate::error::Error::new(format!(
            "{CONTEXT}: not a standard InChI (InChI=1S/) string: {inchi:?}"
        ))
    })?;
    let segments: Vec<&str> = rest.split('/').collect();
    if segments.is_empty() {
        return err(format!("{CONTEXT}: empty identifier body"));
    }

    let mut layers = Layers::default();
    parse_formula(segments[0], &mut layers)?;

    let mut seen_isotope = false;
    for segment in &segments[1..] {
        let head = segment.as_bytes().first().copied().unwrap_or(b'\0');
        match head {
            b'c' => parse_connections(segment, &mut layers)?,
            b'h' if !seen_isotope => parse_hydrogen(segment, &mut layers)?,
            b'q' => parse_charge(segment, &mut layers)?,
            b'p' => parse_protons(segment, &mut layers)?,
            b'i' => {
                seen_isotope = true;
                parse_isotope(segment, &mut layers)?;
            }
            // The isotopic-layer exchangeable-H record: `hD2` after an
            // (empty) /i segment (ParseSegmentIsoExchgH).
            b'h' => parse_exchanged_hydrogen(segment, &mut layers)?,
            b'b' | b't' | b'm' | b's' => {
                return err(format!(
                    "{CONTEXT}: stereo layers (/b /t /m /s) are a follow-up slice: {inchi:?}"
                ));
            }
            b'f' => {
                return err(format!(
                    "{CONTEXT}: Fixed-H layer (/f) is out of the reverse slice: {inchi:?}"
                ));
            }
            _ => {
                return err(format!(
                    "{CONTEXT}: unrecognized layer segment {segment:?} in {inchi:?}"
                ));
            }
        }
    }
    Ok(layers)
}

/// Formula layer → heavy element list (canonical numbering) plus the
/// formula H count. The reference numbers heavy atoms only, in the
/// formula's Hill order, ignoring H entirely (ichiread.c:10129-10169);
/// an H-only formula is the single-atom bridge case (ichiread.c:10171).
fn parse_formula(segment: &str, layers: &mut Layers) -> Result<()> {
    if segment.contains('.') {
        return err(format!(
            "{CONTEXT}: multi-component formula {segment:?} is out of the reverse slice"
        ));
    }
    let bytes = segment.as_bytes();
    let mut index = 0;
    let mut hydrogens = 0u32;
    while index < bytes.len() {
        let upper = bytes[index];
        if !upper.is_ascii_uppercase() {
            return err(format!("{CONTEXT}: malformed formula layer {segment:?}"));
        }
        let mut symbol = String::new();
        symbol.push(upper as char);
        index += 1;
        if index < bytes.len() && bytes[index].is_ascii_lowercase() {
            symbol.push(bytes[index] as char);
            index += 1;
        }
        let mut count = 1u32;
        let start = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        if index > start {
            count = segment[start..index].parse().map_err(|_| {
                crate::error::Error::new(format!(
                    "{CONTEXT}: malformed count in formula {segment:?}"
                ))
            })?;
        }
        if count == 0 {
            return err(format!("{CONTEXT}: zero count in formula {segment:?}"));
        }
        if symbol == "H" {
            hydrogens += count;
        } else {
            if standard_mass(&symbol).is_none() {
                return err(format!(
                    "{CONTEXT}: element {symbol:?} is outside the reverse slice's table"
                ));
            }
            for _ in 0..count {
                layers.elements.push(symbol.clone());
            }
        }
    }
    layers.formula_hydrogens = hydrogens;
    if layers.elements.is_empty() && hydrogens > 0 {
        // H-only formula: the hydrogen becomes atom 1 (the reference's
        // EL_NUMBER_H single-atom fill).
        layers.elements.push("H".to_string());
    }
    Ok(())
}

/// /c layer → undirected bonds. The reference's stack walk
/// (ichiread.c:9502-9617): the first number seeds the chain; `-`
/// chains from the current atom; `(` opens a branch off it whose first
/// number bonds back to the opener; `,` bonds another branch member to
/// the same opener; `)` pops one level, bonds the next number to the
/// popped base and continues the chain from that number.
fn parse_connections(segment: &str, layers: &mut Layers) -> Result<()> {
    let body = &segment[1..];
    if body.contains(';') {
        return err(format!(
            "{CONTEXT}: multi-component /c layer is out of the reverse slice"
        ));
    }
    if body.contains('*') {
        return err(format!(
            "{CONTEXT}: component-multiplied /c layer is out of the reverse slice"
        ));
    }
    if body.bytes().any(|b| b.is_ascii_uppercase()) {
        return err(format!(
            "{CONTEXT}: compressed (abc) /c layer is out of the reverse slice"
        ));
    }
    let bytes = body.as_bytes();
    let mut index = 0usize;
    // Depth stack of "current base" atoms; capacity matches the
    // reference's LAST_AT_LEN bracket without its arbitrary bound.
    let mut stack: Vec<usize> = Vec::new();
    // State: None (expect first number), or the last separator seen.
    #[derive(PartialEq, Clone, Copy)]
    enum State {
        Start,
        Number,
        Open,
        Close,
        Comma,
        Dash,
    }
    let mut state = State::Start;
    let mut highest = 0usize;
    while index < bytes.len() {
        let byte = bytes[index];
        match byte {
            b'(' | b')' | b',' | b'-' => {
                if state != State::Number {
                    return err(format!("{CONTEXT}: misplaced separator in /c {segment:?}"));
                }
                state = match byte {
                    b'(' => State::Open,
                    b')' => State::Close,
                    b',' => State::Comma,
                    _ => State::Dash,
                };
                index += 1;
            }
            b'0'..=b'9' => {
                let start = index;
                while index < bytes.len() && bytes[index].is_ascii_digit() {
                    index += 1;
                }
                let atom: usize = body[start..index].parse().map_err(|_| {
                    crate::error::Error::new(format!(
                        "{CONTEXT}: bad atom number in /c {segment:?}"
                    ))
                })?;
                if atom == 0 || atom > layers.elements.len() + 1 {
                    return err(format!(
                        "{CONTEXT}: /c atom {atom} out of range 1..={} in {segment:?}",
                        layers.elements.len()
                    ));
                }
                highest = highest.max(atom);
                match state {
                    State::Start => {
                        stack.push(atom);
                        state = State::Number;
                    }
                    State::Open => {
                        let base = *stack.last().ok_or_else(|| {
                            crate::error::Error::new(format!(
                                "{CONTEXT}: branch without base in /c {segment:?}"
                            ))
                        })?;
                        layers.bonds.push((base - 1, atom - 1));
                        stack.push(atom);
                        state = State::Number;
                    }
                    State::Close => {
                        if stack.len() < 2 {
                            return err(format!("{CONTEXT}: unbalanced ')' in /c {segment:?}"));
                        }
                        stack.pop();
                        let base = *stack.last().ok_or_else(|| {
                            crate::error::Error::new(format!(
                                "{CONTEXT}: pop without base in /c {segment:?}"
                            ))
                        })?;
                        layers.bonds.push((base - 1, atom - 1));
                        *stack.last_mut().expect("depth") = atom;
                        state = State::Number;
                    }
                    State::Comma => {
                        if stack.len() < 2 {
                            return err(format!(
                                "{CONTEXT}: ',' outside a branch in /c {segment:?}"
                            ));
                        }
                        let base = stack[stack.len() - 2];
                        layers.bonds.push((base - 1, atom - 1));
                        *stack.last_mut().expect("depth") = atom;
                        state = State::Number;
                    }
                    State::Dash => {
                        let base = *stack.last().ok_or_else(|| {
                            crate::error::Error::new(format!(
                                "{CONTEXT}: chain without base in /c {segment:?}"
                            ))
                        })?;
                        layers.bonds.push((base - 1, atom - 1));
                        *stack.last_mut().expect("depth") = atom;
                        state = State::Number;
                    }
                    State::Number => {
                        return err(format!(
                            "{CONTEXT}: adjacent atom numbers in /c {segment:?}"
                        ));
                    }
                }
            }
            _ => {
                return err(format!("{CONTEXT}: unexpected character in /c {segment:?}"));
            }
        }
    }
    if highest > layers.elements.len() {
        // Bridging hydrogens: atoms the /c layer numbers beyond the
        // formula's heavy count are hydrogen (ichiread.c:9412-9433).
        while layers.elements.len() < highest {
            layers.elements.push("H".to_string());
        }
    }
    // Duplicate bonds are the reference's syntax error
    // (FIX_GAF_2020_25726: "same nbr twice").
    let mut seen = std::collections::BTreeSet::new();
    for &(from, to) in &layers.bonds {
        let key = (from.min(to), from.max(to));
        if !seen.insert(key) {
            return err(format!(
                "{CONTEXT}: duplicate bond {from}-{to} in /c {segment:?}"
            ));
        }
    }
    Ok(())
}

/// /h layer → fixed per-atom H counts and mobile groups. The immobile
/// prefix lists `<atoms>H<count>` tokens (ranges allowed); everything
/// from the first `(` is the `(H<count>,atoms...)` group sequence
/// (ichiread.c:8568-8641 immobile, 8848-8990 groups).
fn parse_hydrogen(segment: &str, layers: &mut Layers) -> Result<()> {
    let body = &segment[1..];
    if body.contains(';') || body.contains('*') {
        return err(format!(
            "{CONTEXT}: multi-component /h layer is out of the reverse slice"
        ));
    }
    let n = layers.elements.len();
    layers.fixed_h = vec![0; n];
    let group_start = body.find('(').unwrap_or(body.len());
    let (immobile, groups) = body.split_at(group_start);

    // Immobile tokens: "1-2,3H2,4,5H3" — atoms before the 'H', count
    // after it. Lowercase 'h' reads negative (the FixedH subtractive
    // form; refused here — slice 1 sees 'H' only).
    let mut index = 0usize;
    let bytes = immobile.as_bytes();
    while index < bytes.len() {
        // Atom list: "a", "a-b", or comma-joined lists.
        let mut covered: Vec<usize> = Vec::new();
        loop {
            let start = index;
            while index < bytes.len() && bytes[index].is_ascii_digit() {
                index += 1;
            }
            if index == start {
                return err(format!("{CONTEXT}: expected atom number in /h {segment:?}"));
            }
            let first: usize = immobile[start..index].parse().map_err(|_| {
                crate::error::Error::new(format!("{CONTEXT}: bad atom in /h {segment:?}"))
            })?;
            let mut last = first;
            if index < bytes.len() && bytes[index] == b'-' {
                index += 1;
                let rstart = index;
                while index < bytes.len() && bytes[index].is_ascii_digit() {
                    index += 1;
                }
                last = immobile[rstart..index].parse().map_err(|_| {
                    crate::error::Error::new(format!("{CONTEXT}: bad range in /h {segment:?}"))
                })?;
            }
            if first == 0 || last < first || last > n {
                return err(format!(
                    "{CONTEXT}: /h atom range {first}-{last} out of 1..{n} in {segment:?}"
                ));
            }
            covered.extend(first..=last);
            if index < bytes.len() && bytes[index] == b',' {
                index += 1;
                continue;
            }
            break;
        }
        // The H marker with optional count.
        if index >= bytes.len() || (bytes[index] != b'H' && bytes[index] != b'h') {
            return err(format!(
                "{CONTEXT}: expected H count after atom list in /h {segment:?}"
            ));
        }
        let negative = bytes[index] == b'h';
        index += 1;
        let cstart = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        let count: u32 = if index > cstart {
            immobile[cstart..index].parse().map_err(|_| {
                crate::error::Error::new(format!("{CONTEXT}: bad H count in /h {segment:?}"))
            })?
        } else {
            1
        };
        if negative {
            return err(format!(
                "{CONTEXT}: negative /h counts are FixedH-layer territory: {segment:?}"
            ));
        }
        if count > i32::MAX as u32 {
            return err(format!("{CONTEXT}: absurd H count in /h {segment:?}"));
        }
        for atom in covered {
            layers.fixed_h[atom - 1] = count;
        }
        if index < bytes.len() && bytes[index] == b',' {
            index += 1;
        }
    }

    // Mobile groups: "(H<count>,a,b,...)(H<count2>,...)".
    let mut rest = groups;
    while !rest.is_empty() {
        let bytes = rest.as_bytes();
        if bytes.first() != Some(&b'(') {
            return err(format!("{CONTEXT}: expected (H group in /h {segment:?}"));
        }
        let mut index = 1usize;
        // H with optional count.
        if bytes.get(index) != Some(&b'H') {
            return err(format!("{CONTEXT}: group must start 'H' in /h {segment:?}"));
        }
        index += 1;
        let cstart = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        let hydrogens: u32 = if index > cstart {
            rest[cstart..index].parse().map_err(|_| {
                crate::error::Error::new(format!("{CONTEXT}: bad group H count in /h {segment:?}"))
            })?
        } else {
            1
        };
        // Optional minus charge inside the group — FixedH territory,
        // refused in this slice.
        let mut minus = 0u32;
        if bytes.get(index) == Some(&b'-') {
            index += 1;
            let mstart = index;
            while index < bytes.len() && bytes[index].is_ascii_digit() {
                index += 1;
            }
            minus = if index > mstart {
                rest[mstart..index].parse().map_err(|_| {
                    crate::error::Error::new(format!(
                        "{CONTEXT}: bad group minus count in /h {segment:?}"
                    ))
                })?
            } else {
                1
            };
        }
        if minus > 0 {
            return err(format!(
                "{CONTEXT}: minus-in-group is FixedH territory: {segment:?}"
            ));
        }
        if bytes.get(index) != Some(&b',') {
            return err(format!(
                "{CONTEXT}: expected ',' before endpoints in /h {segment:?}"
            ));
        }
        index += 1;
        let mut endpoints = Vec::new();
        loop {
            let start = index;
            while index < bytes.len() && bytes[index].is_ascii_digit() {
                index += 1;
            }
            if index == start {
                return err(format!(
                    "{CONTEXT}: expected endpoint atom in /h {segment:?}"
                ));
            }
            let atom: usize = rest[start..index].parse().map_err(|_| {
                crate::error::Error::new(format!("{CONTEXT}: bad endpoint in /h {segment:?}"))
            })?;
            if atom == 0 || atom > n {
                return err(format!(
                    "{CONTEXT}: endpoint {atom} out of 1..{n} in /h {segment:?}"
                ));
            }
            endpoints.push(atom - 1);
            if bytes.get(index) == Some(&b',') {
                index += 1;
                continue;
            }
            break;
        }
        if bytes.get(index) != Some(&b')') {
            return err(format!("{CONTEXT}: unterminated group in /h {segment:?}"));
        }
        index += 1;
        layers.groups.push(MobileGroup {
            hydrogens,
            endpoints,
        });
        rest = &rest[index..];
    }
    Ok(())
}

/// /q layer → component total charge (ichiread.c:8096-8130).
fn parse_charge(segment: &str, layers: &mut Layers) -> Result<()> {
    let body = &segment[1..];
    if body.contains(';') || body.contains('*') {
        return err(format!(
            "{CONTEXT}: multi-component /q layer is out of the reverse slice"
        ));
    }
    let value = signed_value(body, segment)?;
    layers.charge = value;
    Ok(())
}

/// /p layer → proton shift (ichiread.c:7450-7465: a single ±N).
fn parse_protons(segment: &str, layers: &mut Layers) -> Result<()> {
    let body = &segment[1..];
    if body.contains(';') {
        return err(format!(
            "{CONTEXT}: multi-component /p layer is out of the reverse slice"
        ));
    }
    let value = signed_value(body, segment)?;
    if value == 0 {
        return err(format!("{CONTEXT}: /p 0 is not a valid proton shift"));
    }
    layers.protons = value;
    Ok(())
}

fn signed_value(body: &str, segment: &str) -> Result<i32> {
    let bytes = body.as_bytes();
    let (sign, start) = match bytes.first() {
        Some(b'+') => (1i32, 1usize),
        Some(b'-') => (-1i32, 1usize),
        _ => {
            return err(format!("{CONTEXT}: expected +/- value in {segment:?}"));
        }
    };
    if !bytes.get(start).is_some_and(|b| b.is_ascii_digit()) {
        return err(format!("{CONTEXT}: expected digits in {segment:?}"));
    }
    let magnitude: i32 = body[start..].parse().map_err(|_| {
        crate::error::Error::new(format!("{CONTEXT}: value out of range in {segment:?}"))
    })?;
    Ok(sign * magnitude)
}

/// /i layer → per-atom isotope entries (ichiread.c:6073-6177: entries
/// `atom[±delta][T[n]][D[n]][H[n]]`, comma separated). An empty body
/// records "no atom isotopes"; a following `h`-segment then carries
/// the exchangeable-H counts.
fn parse_isotope(segment: &str, layers: &mut Layers) -> Result<()> {
    let body = &segment[1..];
    if body.is_empty() || body == "_" {
        return Ok(());
    }
    if body.contains(';') || body.contains('*') {
        return err(format!(
            "{CONTEXT}: multi-component /i layer is out of the reverse slice"
        ));
    }
    for token in body.split(',') {
        let bytes = token.as_bytes();
        let mut index = 0usize;
        let start = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        if index == start {
            return err(format!("{CONTEXT}: expected atom number in /i {segment:?}"));
        }
        let atom: usize = token[start..index].parse().map_err(|_| {
            crate::error::Error::new(format!("{CONTEXT}: bad atom in /i {segment:?}"))
        })?;
        if atom == 0 || atom > layers.elements.len() {
            return err(format!(
                "{CONTEXT}: /i atom {atom} out of 1..{} in {segment:?}",
                layers.elements.len()
            ));
        }
        let mut entry = IsoEntry {
            atom: atom - 1,
            ..Default::default()
        };
        if bytes.get(index) == Some(&b'+') || bytes.get(index) == Some(&b'-') {
            let sign = if bytes[index] == b'+' { 1 } else { -1 };
            index += 1;
            let dstart = index;
            while index < bytes.len() && bytes[index].is_ascii_digit() {
                index += 1;
            }
            if index == dstart {
                return err(format!(
                    "{CONTEXT}: expected delta digits in /i {segment:?}"
                ));
            }
            let magnitude: i32 = token[dstart..index].parse().map_err(|_| {
                crate::error::Error::new(format!("{CONTEXT}: bad delta in /i {segment:?}"))
            })?;
            entry.delta = Some(sign * magnitude);
        }
        while index < bytes.len() {
            let count_of = |_marker: u8, position: &mut usize| -> Result<u32> {
                *position += 1;
                let begin = *position;
                while *position < bytes.len() && bytes[*position].is_ascii_digit() {
                    *position += 1;
                }
                let count: u32 = if *position > begin {
                    token[begin..*position].parse().map_err(|_| {
                        crate::error::Error::new(format!(
                            "{CONTEXT}: bad isotope count in /i {segment:?}"
                        ))
                    })?
                } else {
                    1
                };
                Ok(count)
            };
            match bytes[index] {
                b'T' => {
                    let count = count_of(b'T', &mut index)?;
                    entry.tritiums += count;
                }
                b'D' => {
                    let count = count_of(b'D', &mut index)?;
                    entry.deuteriums += count;
                }
                _ => {
                    return err(format!(
                        "{CONTEXT}: unexpected isotope marker in /i {segment:?}"
                    ));
                }
            }
        }
        if entry.delta.is_none() && entry.deuteriums == 0 && entry.tritiums == 0 {
            return err(format!("{CONTEXT}: empty /i entry in {segment:?}"));
        }
        layers.isotopes.push(entry);
    }
    Ok(())
}

/// The exchangeable-H record after an (empty) /i: `hT2D1`-style
/// counts of exchanged hydrogens (ichiread.c:5589-5609).
fn parse_exchanged_hydrogen(segment: &str, layers: &mut Layers) -> Result<()> {
    let body = &segment[1..];
    let bytes = body.as_bytes();
    let mut index = 0usize;
    let mut deuteriums = 0u32;
    let mut tritiums = 0u32;
    while index < bytes.len() {
        match bytes[index] {
            b'T' => {
                index += 1;
                let start = index;
                while index < bytes.len() && bytes[index].is_ascii_digit() {
                    index += 1;
                }
                let count: u32 = if index > start {
                    body[start..index].parse().map_err(|_| {
                        crate::error::Error::new(format!("{CONTEXT}: bad count in {segment:?}"))
                    })?
                } else {
                    1
                };
                tritiums += count;
            }
            b'D' => {
                index += 1;
                let start = index;
                while index < bytes.len() && bytes[index].is_ascii_digit() {
                    index += 1;
                }
                let count: u32 = if index > start {
                    body[start..index].parse().map_err(|_| {
                        crate::error::Error::new(format!("{CONTEXT}: bad count in {segment:?}"))
                    })?
                } else {
                    1
                };
                deuteriums += count;
            }
            _ => {
                return err(format!(
                    "{CONTEXT}: unexpected exchangeable-H marker in {segment:?}"
                ));
            }
        }
    }
    layers.exchanged = (deuteriums, tritiums);
    Ok(())
}

/// The smallest allowed valence that can still absorb `used` bond
/// orders plus hydrogens (the reference's valence-gap read).
fn valence_ceiling(range: &[f64], used: f64) -> f64 {
    for allowed in range {
        if *allowed + 1e-9 >= used {
            return *allowed;
        }
    }
    *range.last().expect("non-empty")
}

/// The valence floor: the smallest allowed valence at or above the
/// atom's mandatory hydrogens.
fn valence_floor(range: &[f64], hydrogens: u32) -> f64 {
    let mandatory = hydrogens as f64;
    for allowed in range {
        if *allowed + 1e-9 >= mandatory {
            return *allowed;
        }
    }
    *range.last().expect("non-empty")
}

/// A single atom's reconstruction state during assembly.
#[derive(Debug, Clone)]
struct AtomState {
    element: String,
    charge: i32,
    hydrogens: u32,
}

impl AtomState {
    fn allowed(&self) -> Result<Vec<f64>> {
        charged_valences(&self.element, self.charge).ok_or_else(|| {
            crate::error::Error::new(format!(
                "{CONTEXT}: no valence set for {self:?} (element outside the reverse slice)"
            ))
        })
    }

    /// Whether the atom can absorb one more bond order given its
    /// charge and hydrogens.
    fn can_absorb(&self, sum: f64) -> bool {
        match self.allowed() {
            Ok(range) => {
                let used = sum + self.hydrogens as f64;
                let ceiling = valence_ceiling(&range, used);
                used <= ceiling + 1e-9 || used + 1.0 <= *range.last().expect("non-empty") + 1e-9
            }
            Err(_) => false,
        }
    }
}

/// Assembles the structure from the parsed layers: fixed H, mobile-H
/// placement, the /p proton shift, /q charge placement, isotopes, and
/// the Kekulé assignment over the /c skeleton.
fn reconstruct(layers: Layers) -> Result<Structure> {
    let n = layers.elements.len();
    if n == 0 {
        return err(format!("{CONTEXT}: identifier has no atoms"));
    }
    let mut atoms: Vec<AtomState> = layers
        .elements
        .iter()
        .map(|element| AtomState {
            element: element.clone(),
            charge: 0,
            hydrogens: 0,
        })
        .collect();
    for (index, count) in layers.fixed_h.iter().enumerate() {
        atoms[index].hydrogens = *count;
    }

    // The /p proton shift: protons leave (p<0) from the mobile pools
    // first, then from fixed-H heteroatoms; they arrive (p>0) on the
    // first atoms with valence room, spreading across carriers. Each
    // moved proton charges its site ±1 (the reference's
    // AddProtonAndIsoHBalanceToMobHStruct observable behavior).
    let mut groups = layers.groups;
    let mut depleted = vec![0usize; groups.len()];
    if layers.protons < 0 {
        for _ in 0..(-layers.protons) {
            // Removals rotate across the group's last endpoints (the
            // oracle's carbonate: -1 on each of the last two O).
            let pool = groups
                .iter_mut()
                .enumerate()
                .find(|(_, group)| group.hydrogens > 0);
            if let Some((group_index, group)) = pool {
                group.hydrogens -= 1;
                let taken = depleted[group_index];
                let site_index = group
                    .endpoints
                    .len()
                    .saturating_sub(1 + taken.min(group.endpoints.len() - 1));
                let site = group.endpoints[site_index];
                atoms[site].charge -= 1;
                depleted[group_index] = taken + 1;
            } else {
                let site = atoms
                    .iter()
                    .position(|atom| atom.hydrogens > 0)
                    .ok_or_else(|| {
                        crate::error::Error::new(format!(
                            "{CONTEXT}: /p{} has no proton to remove",
                            layers.protons
                        ))
                    })?;
                atoms[site].hydrogens -= 1;
                atoms[site].charge -= 1;
            }
        }
    } else if layers.protons > 0 {
        for _ in 0..layers.protons {
            let bond_sums = bond_sums_of(&layers.bonds, n);
            let mut candidates: Vec<usize> = Vec::new();
            for group in &groups {
                for &endpoint in &group.endpoints {
                    if !candidates.contains(&endpoint) {
                        candidates.push(endpoint);
                    }
                }
            }
            for index in 0..n {
                if !candidates.contains(&index) {
                    candidates.push(index);
                }
            }
            let site = candidates
                .into_iter()
                .filter(|&index| {
                    let mut probe = atoms[index].clone();
                    probe.hydrogens += 1;
                    probe.charge += 1;
                    probe.can_absorb(bond_sums[index])
                })
                .min_by_key(|&index| (atoms[index].hydrogens, index));
            let site = site.ok_or_else(|| {
                crate::error::Error::new(format!(
                    "{CONTEXT}: /p+{} has no atom with valence room",
                    layers.protons
                ))
            })?;
            atoms[site].hydrogens += 1;
            atoms[site].charge += 1;
        }
    }

    // The /q charge rides the first atom whose charged assignment is
    // valence-feasible (the reference's charge-point choice: nitrate's
    // first O, the quaternary N) — charges spread across uncharged
    // atoms first, one unit per step, feasibility proven by an actual
    // Kekulé attempt.
    let mut charge_left = layers.charge;
    if charge_left != 0 {
        let step = if charge_left > 0 { 1 } else { -1 };
        while charge_left != 0 {
            let mut placed = false;
            for prefer_uncharged in [true, false] {
                for index in 0..n {
                    if prefer_uncharged && atoms[index].charge != 0 {
                        continue;
                    }
                    let mut probe = atoms[index].clone();
                    probe.charge += step;
                    if charged_valences(&probe.element, probe.charge).is_none() {
                        continue;
                    }
                    let mut trial = atoms.clone();
                    trial[index].charge += step;
                    if kekulize(&layers.bonds, &trial).is_ok() {
                        atoms[index].charge += step;
                        charge_left -= step;
                        placed = true;
                        break;
                    }
                }
                if placed {
                    break;
                }
            }
            if !placed {
                return err(format!(
                    "{CONTEXT}: /q{} has no atom whose charge yields a valence-feasible assignment",
                    layers.charge
                ));
            }
        }
    }

    // Mobile-H placement: the pool spreads across the group's
    // endpoints, one each per round, while the site's charged valence
    // has room (the BNS's observable hydroxylation pattern).
    for group in &groups {
        let mut pool = group.hydrogens;
        while pool > 0 {
            let mut placed = false;
            let bond_sums = bond_sums_of(&layers.bonds, n);
            for &endpoint in &group.endpoints {
                if pool == 0 {
                    break;
                }
                let mut probe = atoms[endpoint].clone();
                probe.hydrogens += 1;
                if probe.can_absorb(bond_sums[endpoint]) {
                    atoms[endpoint].hydrogens += 1;
                    pool -= 1;
                    placed = true;
                }
            }
            if !placed {
                return err(format!(
                    "{CONTEXT}: mobile group has {} H left with no endpoint room",
                    pool
                ));
            }
        }
    }

    // Isotopes: atom-numbered mass deltas become absolute mass
    // numbers; entry D/T counts and the exchangeable-H record fold
    // onto H-bearing carriers in canonical order (the /h layer is
    // isotope-blind, so the counts ride the plain hydrogens).
    let mut isotopic_h: std::collections::BTreeMap<usize, (u32, u32)> = Default::default();
    let exchange = layers.exchanged;
    for entry in &layers.isotopes {
        let atom = &atoms[entry.atom];
        if let Some(delta) = entry.delta {
            let base = standard_mass(&atom.element).ok_or_else(|| {
                crate::error::Error::new(format!(
                    "{CONTEXT}: no standard mass for element {:?}",
                    atom.element
                ))
            })?;
            let mass = base as i32 + delta;
            if mass <= 0 {
                return err(format!(
                    "{CONTEXT}: isotope delta {} leaves a non-positive mass",
                    delta
                ));
            }
        }
        if entry.deuteriums > 0 || entry.tritiums > 0 {
            fold_isotopic_h(
                &mut atoms,
                &mut isotopic_h,
                entry.atom,
                entry.deuteriums,
                entry.tritiums,
            )?;
        }
    }
    if exchange.0 > 0 || exchange.1 > 0 {
        let mut remaining_d = exchange.0;
        let mut remaining_t = exchange.1;
        for index in 0..n {
            if remaining_d > 0 && atoms[index].hydrogens > 0 {
                let d = remaining_d.min(atoms[index].hydrogens);
                fold_isotopic_h(&mut atoms, &mut isotopic_h, index, d, 0)?;
                remaining_d -= d;
            }
            if remaining_t > 0 && atoms[index].hydrogens > 0 {
                let t = remaining_t.min(atoms[index].hydrogens);
                fold_isotopic_h(&mut atoms, &mut isotopic_h, index, 0, t)?;
                remaining_t -= t;
            }
        }
        if remaining_d > 0 || remaining_t > 0 {
            return err(format!(
                "{CONTEXT}: exchangeable H left over ({} D, {} T) with no carriers",
                remaining_d, remaining_t
            ));
        }
    }

    // Bond orders: single everywhere, then raise orders while any
    // atom sits under its charged valence floor — matching the
    // deficient atoms pairwise over existing bonds (the observable
    // Kekulé rule), one raise at a time.
    let orders = kekulize(&layers.bonds, &atoms)?;

    // Formula consistency: the placed hydrogens (plus D/T) must
    // account for the formula's H count. The reference tolerates the
    // bridging-H slack, so only a deficit errors.
    // The formula layer counts the pre-/p formula: placed = formula + p.
    // The placed hydrogens already include the isotopic ones (the /h
    // layer is isotope-blind; the fold keeps h_count the full count).
    let placed_h: u32 = atoms.iter().map(|a| a.hydrogens).sum::<u32>();
    let expected = layers.formula_hydrogens as i64 + layers.protons as i64;
    if placed_h as i64 != expected {
        return err(format!(
            "{CONTEXT}: placed {} H but formula layer claims {} (+/p {})",
            placed_h, layers.formula_hydrogens, layers.protons
        ));
    }

    let records = atoms
        .into_iter()
        .enumerate()
        .map(|(index, atom)| {
            let isotope = layers
                .isotopes
                .iter()
                .find(|entry| entry.atom == index && entry.delta.is_some())
                .and_then(|entry| {
                    standard_mass(&layers.elements[index])
                        .map(|base| (base as i32 + entry.delta.unwrap()) as u16)
                })
                .unwrap_or(0);
            AtomRecord {
                element: atom.element,
                h_count: atom.hydrogens,
                charge: atom.charge,
                h_pinned: false,
                isotope,
            }
        })
        .collect();
    let bonds = layers
        .bonds
        .iter()
        .cloned()
        .zip(orders)
        .map(|((from, to), order)| (from, to, order))
        .collect();
    Ok(Structure {
        atoms: records,
        bonds,
        stereo_bonds: Default::default(),
        positions: Vec::new(),
        isotopic_h,
        protons_hard_moved: false,
        ride_marks: Vec::new(),
        fold_source: Vec::new(),
        chiral_flag: false,
    })
}

/// Bond-order sums over the current orders.
fn ordered_sums(bonds: &[(usize, usize)], orders: &[u32], atoms: usize) -> Vec<f64> {
    let mut sums = vec![0.0; atoms];
    for (&(from, to), &order) in bonds.iter().zip(orders) {
        let order = f64::from(order);
        sums[from] += order;
        sums[to] += order;
    }
    sums
}

/// Plain single-bond sums (before order assignment).
fn bond_sums_of(bonds: &[(usize, usize)], atoms: usize) -> Vec<f64> {
    ordered_sums(bonds, &vec![1; bonds.len()], atoms)
}

/// Folds isotopic hydrogens onto a carrier's plain count.
fn fold_isotopic_h(
    atoms: &mut [AtomState],
    isotopic_h: &mut std::collections::BTreeMap<usize, (u32, u32)>,
    carrier: usize,
    deuteriums: u32,
    tritiums: u32,
) -> Result<()> {
    if atoms[carrier].hydrogens < deuteriums + tritiums {
        return err(format!(
            "{CONTEXT}: {} isotopic H exceed the carrier's plain count",
            deuteriums + tritiums
        ));
    }
    let slot = isotopic_h.entry(carrier).or_insert((0, 0));
    slot.0 += deuteriums;
    slot.1 += tritiums;
    Ok(())
}

/// The Kekulé assignment: raise bond orders while any atom is under
/// its charged valence floor, pairing deficient atoms over existing
/// bonds. A maximum matching pass covers rings (benzene needs a
/// perfect matching over its six deficient carbons); paired raises
/// then carry atoms needing two more orders (nitriles' C≡N).
fn kekulize(bonds: &[(usize, usize)], atoms: &[AtomState]) -> Result<Vec<u32>> {
    let mut orders = vec![1u32; bonds.len()];
    assign_orders(bonds, &mut orders, atoms)?;
    Ok(orders)
}

fn assign_orders(bonds: &[(usize, usize)], orders: &mut [u32], atoms: &[AtomState]) -> Result<()> {
    let mut sums = ordered_sums(bonds, orders, atoms.len());
    loop {
        // Deficiency: how much order each atom still needs to reach
        // its valence floor.
        let mut deficient: Vec<bool> = vec![false; atoms.len()];
        let mut any = false;
        for (index, atom) in atoms.iter().enumerate() {
            let range = atom.allowed()?;
            let floor = valence_floor(&range, atom.hydrogens);
            if sums[index] + atom.hydrogens as f64 + 1e-9 < floor {
                deficient[index] = true;
                any = true;
            }
        }
        if !any {
            // Every atom must also sit at or under its charged
            // ceiling — a carbocation charged onto a tetravalent C
            // passes the floor check but overflows here.
            for (index, atom) in atoms.iter().enumerate() {
                if !valid_after(atom, sums[index]) {
                    return err(format!(
                        "{CONTEXT}: atom {} ({}) count {} + {} H lands outside its charged valence range",
                        index + 1, atom.element, sums[index], atom.hydrogens
                    ));
                }
            }
            return Ok(());
        }

        // Maximum matching over the deficient-atom subgraph: greedy
        // pairs first, then augmenting paths (rings need them — a
        // greedy raise on a six-ring strands the far side).
        let deficient_atoms: Vec<usize> = (0..atoms.len()).filter(|&i| deficient[i]).collect();
        let mut matched: std::collections::BTreeMap<usize, usize> = Default::default();
        for &atom in &deficient_atoms {
            if matched.contains_key(&atom) {
                continue;
            }
            for &(from, to) in bonds {
                let other = if from == atom && deficient[to] {
                    Some(to)
                } else if to == atom && deficient[from] {
                    Some(from)
                } else {
                    None
                };
                if let Some(other) = other {
                    if !matched.contains_key(&other) && other != atom {
                        matched.insert(atom, other);
                        matched.insert(other, atom);
                        break;
                    }
                }
            }
        }
        // Augmenting paths over alternating matched/unmatched bonds.
        for &root in &deficient_atoms {
            if matched.contains_key(&root) {
                continue;
            }
            let mut visited: Vec<usize> = vec![root];
            augment(root, bonds, &deficient, &mut matched, &mut visited);
        }

        // Raise the matched pairs and any remaining same-bond pair
        // (the triple-bond case: both endpoints still short after the
        // first raise re-evaluates).
        let mut raised = false;
        for (atom, other) in matched.iter() {
            if *atom < *other {
                if let Some(position) = bonds.iter().position(|&(from, to)| {
                    (from == *atom && to == *other) || (from == *other && to == *atom)
                }) {
                    if orders[position] < 3 {
                        orders[position] += 1;
                        let order = f64::from(orders[position]);
                        sums[*atom] += 1.0;
                        sums[*other] += 1.0;
                        let _ = order;
                        raised = true;
                    }
                }
            }
        }
        // Paired raises on a shared bond for leftover deficiencies
        // (both endpoints deficient on the same edge, e.g. C≡N's
        // second raise).
        for (position, &(from, to)) in bonds.iter().enumerate() {
            if orders[position] >= 3 {
                continue;
            }
            let range_from = atoms[from].allowed()?;
            let range_to = atoms[to].allowed()?;
            let floor_from = valence_floor(&range_from, atoms[from].hydrogens);
            let floor_to = valence_floor(&range_to, atoms[to].hydrogens);
            let short_from = sums[from] + atoms[from].hydrogens as f64 + 1e-9 < floor_from;
            let short_to = sums[to] + atoms[to].hydrogens as f64 + 1e-9 < floor_to;
            if short_from && short_to {
                orders[position] += 1;
                sums[from] += 1.0;
                sums[to] += 1.0;
                raised = true;
            }
        }
        // One-sided raises: a deficient atom raises its bond to a
        // partner that can absorb the order (sulfate's S=O from the
        // unhydroxylated O; nitrate's N=O2/N=O3), only where BOTH
        // endpoints stay inside their charged ranges.
        if !raised {
            for (position, &(from, to)) in bonds.iter().enumerate() {
                if orders[position] >= 3 {
                    continue;
                }
                let raise_helps = |short_side: usize, absorber: usize| -> bool {
                    let Ok(range) = atoms[short_side].allowed() else {
                        return false;
                    };
                    let floor = valence_floor(&range, atoms[short_side].hydrogens);
                    let short =
                        sums[short_side] + atoms[short_side].hydrogens as f64 + 1e-9 < floor;
                    short
                        && atoms[absorber].can_absorb(sums[absorber])
                        && valid_after(&atoms[short_side], sums[short_side] + 1.0)
                        && valid_after(&atoms[absorber], sums[absorber] + 1.0)
                };
                if raise_helps(from, to) || raise_helps(to, from) {
                    orders[position] += 1;
                    sums[from] += 1.0;
                    sums[to] += 1.0;
                    raised = true;
                }
            }
        }
        if !raised {
            return err(format!(
                "{CONTEXT}: no valence-feasible bond assignment exists for the /c skeleton"
            ));
        }
    }
}

/// Whether `used` bond orders plus the atom's hydrogens stay inside
/// its charged valence range.
fn valid_after(atom: &AtomState, used: f64) -> bool {
    match atom.allowed() {
        Ok(range) => {
            let total = used + atom.hydrogens as f64;
            let floor = valence_floor(&range, atom.hydrogens);
            let ceiling = *range.last().expect("non-empty");
            total >= floor - 1e-9 && total <= ceiling + 1e-9
        }
        Err(_) => false,
    }
}

/// One augmenting-path attempt: walk alternating bonds from `root`
/// toward an unmatched deficient atom, flipping matches. Returns the
/// flip target when the path length allows (small molecules keep this
/// depth-first walk trivial).
fn augment(
    root: usize,
    bonds: &[(usize, usize)],
    deficient: &[bool],
    matched: &mut std::collections::BTreeMap<usize, usize>,
    visited: &mut Vec<usize>,
) -> bool {
    for &(from, to) in bonds {
        let next = if from == root {
            Some(to)
        } else if to == root {
            Some(from)
        } else {
            None
        };
        let Some(next) = next else { continue };
        if !deficient[next] || visited.contains(&next) {
            continue;
        }
        visited.push(next);
        let exposed = match matched.get(&next) {
            None => true,
            Some(&partner) => {
                partner != root && augment(partner, bonds, deficient, matched, visited)
            }
        };
        if exposed {
            // Flip along the walk: the reassignments overwrite the
            // previous partner mappings.
            matched.insert(root, next);
            matched.insert(next, root);
            return true;
        }
    }
    false
}

/// Renders a reconstructed structure as the V2000 molfile the engine
/// path reads (the same subset the corpus identity molfiles use).
/// Isotopic hydrogens render as explicit D/T atoms bonded to their
/// carrier (the corpus fixture convention) so the reader refolds them
/// into the carrier's count and isotope map identically.
pub fn structure_to_molfile(structure: &Structure) -> String {
    // Isotopic hydrogens leave the carrier's rendered count and come
    // back as explicit D/T atoms (isotope mass 2/3).
    let mut plain_h: Vec<u32> = structure.atoms.iter().map(|atom| atom.h_count).collect();
    let mut deuterium_atoms: Vec<(usize, u16)> = Vec::new();
    for (&carrier, &(d, t)) in &structure.isotopic_h {
        plain_h[carrier] = plain_h[carrier].saturating_sub(d + t);
        for _ in 0..d {
            deuterium_atoms.push((carrier, 2));
        }
        for _ in 0..t {
            deuterium_atoms.push((carrier, 3));
        }
    }
    let na = structure.atoms.len() + deuterium_atoms.len();
    let nb = structure.bonds.len() + deuterium_atoms.len();
    let mut out = String::from("reverse\n  chemomer-reverse\n\n");
    out.push_str(&format!("{:3}{:3}  0  0  0  0  0  0  0999 V2000\n", na, nb));
    for (index, atom) in structure.atoms.iter().enumerate() {
        let plain = plain_h[index];
        let code = if plain > 0 { plain + 1 } else { 0 };
        out.push_str(&format!(
            "   12.0000   12.0000    0.0000 {:<3} 0  0  {:>2}  0  0  0  0  0  0  0  0\n",
            atom.element, code
        ));
    }
    for _ in &deuterium_atoms {
        out.push_str(&format!(
            "   12.0000   12.0000    0.0000 {:<3} 0  0  {:>2}  0  0  0  0  0  0  0  0\n",
            "H", 0
        ));
    }
    for &(from, to, kind) in &structure.bonds {
        out.push_str(&format!(
            "{:3}{:3}{:3}  0  0  0  0  0  0  0\n",
            from + 1,
            to + 1,
            kind
        ));
    }
    for (position, (carrier, _mass)) in deuterium_atoms.iter().enumerate() {
        let atom_number = structure.atoms.len() + position + 1;
        out.push_str(&format!(
            "{:3}{:3}  1  0  0  0  0  0  0  0\n",
            carrier + 1,
            atom_number
        ));
    }
    let charged: Vec<(usize, i32)> = structure
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.charge != 0)
        .map(|(index, atom)| (index + 1, atom.charge))
        .collect();
    if !charged.is_empty() {
        out.push_str(&format!("M  CHG{:3}", charged.len()));
        for (index, charge) in &charged {
            out.push_str(&format!(" {:3}{:4}", index, charge));
        }
        out.push('\n');
    }
    let isotopic: Vec<(usize, u16)> = structure
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.isotope != 0)
        .map(|(index, atom)| (index + 1, atom.isotope))
        .collect();
    for (index, mass) in &isotopic {
        out.push_str(&format!("M  ISO  1{:4}{:4}\n", index, mass));
    }
    for (position, (_carrier, mass)) in deuterium_atoms.iter().enumerate() {
        let atom_number = structure.atoms.len() + position + 1;
        out.push_str(&format!("M  ISO  1{:4}{:4}\n", atom_number, mass));
    }
    out.push_str("M  END\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Round-trips a reconstruction through the molfile render into
    /// the engine and returns the assembled identifier.
    fn forward(inchi: &str) -> String {
        let structure = structure_from_inchi(inchi).expect("reconstructs");
        let molfile = structure_to_molfile(&structure);
        let reparsed = crate::molfile::parse(&molfile).expect("re-reads");
        crate::Engine::from_molfile(&structure_to_molfile(&reparsed))
            .and_then(|engine| engine.inchi_string())
            .map_err(|error| format!("{error:?}"))
            .expect("forwards")
    }

    #[test]
    fn water_reconstructs_and_forwards_itself() {
        assert_eq!(forward("InChI=1S/H2O/h1H2"), "InChI=1S/H2O/h1H2");
    }

    #[test]
    fn methane_reconstructs_and_forwards_itself() {
        assert_eq!(forward("InChI=1S/CH4/h1H4"), "InChI=1S/CH4/h1H4");
    }

    #[test]
    fn ethanol_reconstructs_the_chain() {
        let structure = structure_from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3").expect("parses");
        assert_eq!(structure.atoms.len(), 3);
        assert_eq!(structure.atoms[2].element, "O");
        assert_eq!(structure.atoms[2].h_count, 1);
        assert_eq!(structure.bonds, vec![(0, 1, 1), (1, 2, 1)]);
        assert_eq!(
            forward("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"),
            "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        );
    }

    #[test]
    fn pyrrole_kekule_keeps_nitrogen_hydrogen() {
        let structure = structure_from_inchi("InChI=1S/C4H5N/c1-2-4-5-3-1/h1-5H").expect("parses");
        assert_eq!(structure.atoms[4].element, "N");
        assert_eq!(structure.atoms[4].h_count, 1);
        // Every carbon carries exactly one H; two ring bonds doubled.
        let doubles = structure.bonds.iter().filter(|&&(_, _, k)| k == 2).count();
        assert_eq!(doubles, 2);
        assert_eq!(
            forward("InChI=1S/C4H5N/c1-2-4-5-3-1/h1-5H"),
            "InChI=1S/C4H5N/c1-2-4-5-3-1/h1-5H"
        );
    }

    #[test]
    fn naphthalene_paren_branches_and_fusion() {
        let inchi = "InChI=1S/C10H8/c1-2-6-10-8-4-3-7-9(10)5-1/h1-8H";
        let structure = structure_from_inchi(inchi).expect("parses");
        assert_eq!(structure.atoms.len(), 10);
        assert_eq!(structure.bonds.len(), 11);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn acetic_acid_mobile_h_lands_and_forwards() {
        let inchi = "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)";
        let structure = structure_from_inchi(inchi).expect("parses");
        // Mobile pool: one H over the two O; charge nowhere.
        assert_eq!(structure.total_charge(), 0);
        let h_total: u32 = structure.atoms.iter().map(|a| a.h_count).sum();
        assert_eq!(h_total, 4);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn acetate_proton_removal_charges_the_site() {
        let inchi = "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)/p-1";
        let structure = structure_from_inchi(inchi).expect("parses");
        assert_eq!(structure.total_charge(), -1);
        let h_total: u32 = structure.atoms.iter().map(|a| a.h_count).sum();
        assert_eq!(h_total, 3);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn ammonium_proton_addition_gives_onium() {
        let inchi = "InChI=1S/H3N/h1H3/p+1";
        let structure = structure_from_inchi(inchi).expect("parses");
        assert_eq!(structure.atoms[0].h_count, 4);
        assert_eq!(structure.atoms[0].charge, 1);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn hydroxide_removes_a_fixed_hydrogen() {
        let inchi = "InChI=1S/H2O/h1H2/p-1";
        let structure = structure_from_inchi(inchi).expect("parses");
        assert_eq!(structure.atoms[0].h_count, 1);
        assert_eq!(structure.atoms[0].charge, -1);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn carbonate_dianion_spreads_the_charge() {
        let inchi = "InChI=1S/CH2O3/c2-1(3)4/h(H2,2,3,4)/p-2";
        let structure = structure_from_inchi(inchi).expect("parses");
        assert_eq!(structure.total_charge(), -2);
        let h_total: u32 = structure.atoms.iter().map(|a| a.h_count).sum();
        assert_eq!(h_total, 0);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn sulfuric_acid_mobile_h3_over_four_oxygens() {
        let inchi = "InChI=1S/H4O4S/c1-5(2,3)4/h5H,(H3,1,2,3,4)";
        let structure = structure_from_inchi(inchi).expect("parses");
        assert_eq!(structure.atoms[4].element, "S");
        assert_eq!(structure.atoms[4].h_count, 1);
        let h_total: u32 = structure.atoms.iter().map(|a| a.h_count).sum();
        assert_eq!(h_total, 4);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn nitrate_q_minus_one_lands_on_oxygen() {
        let inchi = "InChI=1S/NO3/c2-1(3)4/q-1";
        let structure = structure_from_inchi(inchi).expect("parses");
        assert_eq!(structure.total_charge(), -1);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn c13_methane_isotope_delta() {
        let inchi = "InChI=1S/CH4/h1H4/i1+1";
        let structure = structure_from_inchi(inchi).expect("parses");
        assert_eq!(structure.atoms[0].isotope, 13);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn d2o_exchanged_hydrogens_fold() {
        let inchi = "InChI=1S/H2O/h1H2/i/hD2";
        let structure = structure_from_inchi(inchi).expect("parses");
        assert_eq!(structure.isotopic_h.get(&0), Some(&(2, 0)));
        // The /h layer is isotope-blind: the carrier keeps the full
        // count (the engine's fold convention).
        assert_eq!(structure.atoms[0].h_count, 2);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn tetrahexylammonium_q_plus_one_rides_nitrogen() {
        let inchi = "InChI=1S/C16H36N/c1-5-9-13-17(14-10-6-2,15-11-7-3)16-12-8-4/h5-16H2,1-4H3/q+1";
        let structure = structure_from_inchi(inchi).expect("parses");
        let charged: Vec<_> = structure
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.charge != 0)
            .collect();
        assert_eq!(charged.len(), 1);
        assert_eq!(charged[0].1.element, "N");
        assert_eq!(charged[0].1.charge, 1);
        assert_eq!(forward(inchi), inchi);
    }

    #[test]
    fn stereo_layers_are_a_follow_up_slice() {
        let error =
            structure_from_inchi("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1")
                .expect_err("stereo refused");
        assert!(error.message().contains("stereo layers"));
    }

    #[test]
    fn multi_component_formula_refused() {
        let error =
            structure_from_inchi("InChI=1S/CH4.H2O/c;1-2/h;1H2").expect_err("multi refused");
        assert!(error.message().contains("multi-component"));
    }

    #[test]
    fn non_standard_prefix_refused() {
        let error = structure_from_inchi("InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3")
            .expect_err("non-standard refused");
        assert!(error.message().contains("not a standard InChI"));
    }

    #[test]
    fn glycine_amino_acid_reconstruction() {
        let inchi = "InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)";
        let structure = structure_from_inchi(inchi).expect("parses");
        // N3 with two H on the aminomethyl carbon C1.
        assert_eq!(structure.atoms[2].element, "N");
        assert_eq!(structure.atoms[2].h_count, 2);
        assert_eq!(forward(inchi), inchi);
    }
}
