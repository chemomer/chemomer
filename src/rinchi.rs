//! The RInChI (Reaction InChI) reaction-identity layer — a reaction
//! (reactants + products + agents, each a molecule) gets its own
//! identifier assembled from the per-component standard InChIs.
//!
//! Ported from the IUPAC RInChI reference (rinchi_reaction.cpp +
//! rinchi_consts.cpp + rinchi_hashing.cpp; IUPAC/InChI-Trust Licence
//! No. 1.0, source read at development time, never compiled, linked,
//! or vendored). The corpus ships only recorded oracle strings.
//!
//! Grammar of the identifier (version 1.00, standard InChI
//! constituents only):
//!
//! `RInChI=1.00.1S/<group1><><group2><><group3>/d<flag>`
//!
//! - a *group* holds one role's components, each component's
//!   identifier minus the `InChI=1S/` prefix, joined by `!`
//!   (DELIM_COMP); groups join with `<>` (DELIM_GROUP).
//! - within a role the components sort lexically by their full
//!   InChI string; roles output in lexical order too — when the
//!   products group sorts before the reactants group the reaction
//!   outputs reversed and the direction flag reads `-` (backward)
//!   instead of `+` (forward); equilibrium reads `=` regardless.
//! - trailing roles carry no content: with no agents only the first
//!   two groups output; a role whose every component is a
//!   No-Structure contributes nothing to the groups.
//! - `/u` appears when any role holds No-Structures (a component
//!   with no drawn structure, or a placeholder `A`/`X`/`R` atom):
//!   `/u<r>-<p>-<a>` counts them per role in output order.
//!
//! The RAuxInfo line mirrors the assembly over the per-component
//! `AuxInfo=1/` bodies (`RAuxInfo=1.00.1/<group1><><group2>`); an
//! empty aux body contributes a bare `/`.
//!
//! The RInChIKey comes in three shapes; the short key is the hashed
//! one: `Short-RInChIKey=SA-<F|B|E>UHFF-<major r>-<major p>-<major
//! a>-<minor r>-<minor p>-<minor a>-<nostruct chars>` where each
//! major hash is a 10-char truncated SHA-256 (the InChIKey base-26
//! tables) over the role's major layers — formula plus `c`/`h`/`q`,
//! one per component, joined by `!` — and each minor hash is the
//! role's aggregated protonation letter (`N` at zero, ±1..±12 as
//! `O`..`Z`/`M`..`B`, `A` beyond) over a 4-char hash of the
//! remaining layers. The long key concatenates whole per-component
//! InChIKeys (No-Structures hash to `MOSFIJXAXDLOML-UHFFFAOYSA-N`),
//! and the web key hashes the role-independent set of distinct
//! component identifiers (17-char major, proton letter + 12-char
//! minor, `SA` tail).

use crate::error::{err, Result};
use crate::key;

const CONTEXT: &str = "rxn";

// String constants — rinchi_consts.cpp.
const RINCHI_STD_HEADER: &str = "RInChI=1.00.1S/";
const RINCHI_AUXINFO_HEADER: &str = "RAuxInfo=1.00.1/";
const DELIM_LAYER: char = '/';
const DELIM_COMP: char = '!';
const DELIM_GROUP: &str = "<>";
const DIRECTION_TAG: &str = "/d";
const DIRECTION_FORWARD: char = '+';
const DIRECTION_REVERSE: char = '-';
const DIRECTION_EQUILIBRIUM: char = '=';
const NOSTRUCT_TAG: &str = "/u";
const NOSTRUCT_DELIM: char = '-';
const INCHI_STD_HEADER: &str = "InChI=1S/";
const INCHI_AUXINFO_HEADER: &str = "AuxInfo=1/";
const NOSTRUCT_INCHI: &str = "InChI=1S/";
const NOSTRUCT_AUXINFO: &str = "AuxInfo=1/";
/// The No-Structure InChIKey as used in the long key (rinchi_consts.cpp).
const NOSTRUCT_RINCHI_LONGKEY: &str = "MOSFIJXAXDLOML-UHFFFAOYSA-N";

// RInChI key constants.
const RINCHI_LONG_KEY_HEADER: &str = "Long-RInChIKey=";
const RINCHI_SHORT_KEY_HEADER: &str = "Short-RInChIKey=";
const RINCHI_WEB_KEY_HEADER: &str = "Web-RInChIKey=";
const RINCHI_KEY_VERSION_ID: &str = "SA";
const KEY_DELIM_BLOCK: &str = "-";
const KEY_DELIM_COMP: &str = "-";
const KEY_DELIM_GROUP: &str = "--";

// Empty-string hash constants — rinchi_hashing.cpp; the unit tests
// re-derive them from the native hasher.
const HASH_04_EMPTY: &str = "UHFF";
const HASH_12_EMPTY: &str = "UHFFFADPSCTJ";

// MDL RXNfile tags — mdl_rxnfile.h.
const MDL_TAG_RXN_BEGIN: &str = "$RXN";
const MDL_TAG_RXN_COMPONENT_START: &str = "$MOL";
const MDL_TAG_MOLFILE_END: &str = "M  END";

/// The atom symbols that mark a placeholder structure — a component
/// made of one of these is a No-Structure (mdl_molfile_reader.cpp's
/// `is_special_atom`).
const SPECIAL_ATOMS: [&str; 5] = ["A", "X", "R#", "R", "*"];

/// Does the reaction run from reactants to products, or is it an
/// equilibrium? (rinchi_reaction.h's ReactionDirectionality; the
/// RXNfile format carries no direction field, so parsing always
/// yields directional.)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Directionality {
    Directional,
    Equilibrium,
}

/// One reaction component: a molecule as drawn by its V2000 molfile
/// block (the `$MOL` body, header lines through `M  END`), or a
/// No-Structure placeholder.
#[derive(Debug, Clone)]
pub struct RinchiComponent {
    molfile: String,
}

impl RinchiComponent {
    /// The component's V2000 molfile block exactly as parsed.
    pub fn molfile(&self) -> &str {
        &self.molfile
    }

    /// A No-Structure: the block holds no atoms, or a single
    /// placeholder atom (`A`, `X`, `R#`, `R`, `*`) — both are the
    /// reference's "unknown structure" (mdl_molfile_reader.cpp's
    /// special-atom rule, which bypasses the InChI library entirely).
    fn is_no_structure(&self) -> bool {
        let lines: Vec<&str> = self.molfile.lines().collect();
        let Some(counts) = lines.get(3) else {
            return false;
        };
        let atom_count: usize = counts.get(0..3).unwrap_or("").trim().parse().unwrap_or(0);
        if atom_count == 0 {
            return true;
        }
        if atom_count > 1 {
            return false;
        }
        let symbol = lines
            .get(4)
            .map(|line| line.get(31..34).unwrap_or("").trim());
        matches!(symbol, Some(symbol) if SPECIAL_ATOMS.contains(&symbol))
    }

    /// The component's identity triple: standard InChI, AuxInfo, and
    /// InChIKey — computed through the existing engine for molecules,
    /// fixed constants for No-Structures (the reference never feeds a
    /// No-Structure to the InChI library; inchi_generator.cpp).
    fn identity(&self) -> Result<Identity> {
        if self.is_no_structure() {
            return Ok(Identity {
                inchi: NOSTRUCT_INCHI.to_string(),
                aux: NOSTRUCT_AUXINFO.to_string(),
                key: NOSTRUCT_RINCHI_LONGKEY.to_string(),
            });
        }
        let engine = crate::Engine::from_molfile(&self.molfile)?;
        Ok(Identity {
            inchi: engine.inchi_string()?,
            aux: engine.aux_string()?,
            key: engine.inchikey()?,
        })
    }
}

/// A computed component identity.
#[derive(Debug, Clone)]
struct Identity {
    inchi: String,
    aux: String,
    key: String,
}

/// A reaction: reactants, products, and agents, each list holding
/// whole molecules (multi-fragment molfile blocks included).
#[derive(Debug, Clone)]
pub struct RinchiReaction {
    reactants: Vec<RinchiComponent>,
    products: Vec<RinchiComponent>,
    agents: Vec<RinchiComponent>,
    directionality: Directionality,
}

impl RinchiReaction {
    pub fn reactants(&self) -> &[RinchiComponent] {
        &self.reactants
    }

    pub fn products(&self) -> &[RinchiComponent] {
        &self.products
    }

    pub fn agents(&self) -> &[RinchiComponent] {
        &self.agents
    }

    pub fn directionality(&self) -> Directionality {
        self.directionality
    }

    /// Forces the equilibrium reading (the reference CLI's
    /// `/force-equilibrium`).
    pub fn set_directionality(&mut self, directionality: Directionality) {
        self.directionality = directionality;
    }

    /// The assembled reaction identifier.
    pub fn rinchi_string(&self) -> Result<String> {
        Ok(rinchi_string_of(self.directionality, &self.groups()?))
    }

    /// The RAuxInfo line beside the identifier.
    pub fn rinchi_auxinfo(&self) -> Result<String> {
        Ok(rinchi_auxinfo_of(&self.groups()?))
    }

    /// The long key: whole per-component InChIKeys per role
    /// (No-Structures contribute the fixed No-Structure key).
    pub fn rinchi_long_key(&self) -> Result<String> {
        rinchi_long_key_of(self.directionality, &self.groups()?)
    }

    /// The short key: truncated hashes over the per-role layer
    /// strings — the hashed, database-friendly shape.
    pub fn rinchi_short_key(&self) -> Result<String> {
        rinchi_short_key_of(self.directionality, &self.groups()?)
    }

    /// The web key: one hash over the distinct component identifiers
    /// of the whole reaction, role-independent.
    pub fn rinchi_web_key(&self) -> Result<String> {
        rinchi_web_key_of(&self.groups()?)
    }

    /// The per-role assembly cache — rinchi_reaction.cpp's
    /// `update_cache`: per role, the No-Structures counted and split
    /// off, the rest sorted by InChI string, the group strings
    /// built, and the output order decided by lexical comparison of
    /// the products group against the reactants group.
    fn groups(&self) -> Result<Groups> {
        let mut groups = Groups {
            rinchi_groups: Default::default(),
            rauxinfo_groups: Default::default(),
            identities: Default::default(),
            nostruct_counts: [0; 3],
            output_order: [0, 1, 2],
            reverse_output: false,
        };
        let roles = [&self.reactants, &self.products, &self.agents];
        for (index, role) in roles.iter().enumerate() {
            let mut sorted = Vec::new();
            for component in *role {
                if component.is_no_structure() {
                    groups.nostruct_counts[index] += 1;
                } else {
                    sorted.push(component.identity()?);
                }
            }
            sorted.sort_by(|a, b| a.inchi.cmp(&b.inchi));
            let mut rinchi_group = String::new();
            let mut rauxinfo_group = String::new();
            for (position, identity) in sorted.iter().enumerate() {
                if position > 0 {
                    rinchi_group.push(DELIM_COMP);
                    rauxinfo_group.push(DELIM_COMP);
                }
                // InChI minus the leading "InChI=1S/".
                rinchi_group.push_str(identity.inchi.strip_prefix(INCHI_STD_HEADER).ok_or_else(
                    || {
                        bad(&format!(
                            "component is not a standard InChI: {:?}",
                            identity.inchi
                        ))
                    },
                )?);
                // AuxInfo minus the leading "AuxInfo=1/"; an empty
                // body contributes a bare layer delimiter.
                if identity.aux.is_empty() {
                    rauxinfo_group.push(DELIM_LAYER);
                } else {
                    rauxinfo_group.push_str(
                        identity
                            .aux
                            .strip_prefix(INCHI_AUXINFO_HEADER)
                            .ok_or_else(|| {
                                bad(&format!(
                                    "component has a non-AuxInfo line: {:?}",
                                    identity.aux
                                ))
                            })?,
                    );
                }
            }
            groups.rinchi_groups[index] = rinchi_group;
            groups.rauxinfo_groups[index] = rauxinfo_group;
            groups.identities[index] = sorted;
        }
        // Is the output reversed? Products group sorts before the
        // reactants group.
        groups.reverse_output = groups.rinchi_groups[1] < groups.rinchi_groups[0];
        if groups.reverse_output {
            groups.output_order.swap(0, 1);
        }
        Ok(groups)
    }
}

#[derive(Debug)]
struct Groups {
    rinchi_groups: [String; 3],
    rauxinfo_groups: [String; 3],
    identities: [Vec<Identity>; 3],
    nostruct_counts: [usize; 3],
    output_order: [usize; 3],
    reverse_output: bool,
}

impl Groups {
    /// How many groups output — the last output position holding any
    /// content (rinchi_reaction.cpp's num_output_groups); the key
    /// shapes also count No-Structure-bearing roles.
    fn num_output_groups(&self, consider_nostruct_counts: bool) -> usize {
        let mut result = 1;
        for position in 0..3 {
            let index = self.output_order[position];
            if !self.rinchi_groups[index].is_empty() {
                result = position + 1;
            }
            if consider_nostruct_counts && self.nostruct_counts[index] != 0 {
                result = position + 1;
            }
        }
        result
    }
}

/// The identifier body — rinchi_reaction.cpp's rinchi_string.
fn rinchi_string_of(directionality: Directionality, groups: &Groups) -> String {
    let mut result = String::from(RINCHI_STD_HEADER);
    let group_count = groups.num_output_groups(false);
    for position in 0..group_count {
        result.push_str(&groups.rinchi_groups[groups.output_order[position]]);
        if position < group_count - 1 {
            result.push_str(DELIM_GROUP);
        }
    }
    result.push_str(&direction_flag(directionality, groups.reverse_output));
    if groups.nostruct_counts.iter().any(|&count| count > 0) {
        result.push_str(NOSTRUCT_TAG);
        for position in 0..3 {
            result.push_str(&groups.nostruct_counts[groups.output_order[position]].to_string());
            if position < 2 {
                result.push(NOSTRUCT_DELIM);
            }
        }
    }
    result
}

/// The RAuxInfo body — rinchi_reaction.cpp's rinchi_auxinfo.
fn rinchi_auxinfo_of(groups: &Groups) -> String {
    let mut result = String::from(RINCHI_AUXINFO_HEADER);
    let group_count = groups.num_output_groups(false);
    for position in 0..group_count {
        result.push_str(&groups.rauxinfo_groups[groups.output_order[position]]);
        if position < group_count - 1 {
            result.push_str(DELIM_GROUP);
        }
    }
    result
}

/// The long key body — rinchi_reaction.cpp's rinchi_long_key.
fn rinchi_long_key_of(directionality: Directionality, groups: &Groups) -> Result<String> {
    let mut result = format!(
        "{RINCHI_LONG_KEY_HEADER}{RINCHI_KEY_VERSION_ID}{KEY_DELIM_BLOCK}{}{}",
        direction_code(directionality, groups.reverse_output),
        // The reserved "additional reaction layers" block hashes to
        // the empty string's first four characters, run straight
        // against the direction letter (rinchi_reaction.cpp).
        &HASH_12_EMPTY[..4],
    );
    result.push_str(KEY_DELIM_BLOCK);
    let empty_len = result.len();
    let group_count = groups.num_output_groups(true);
    for position in 0..group_count {
        let index = groups.output_order[position];
        let mut items: Vec<&str> = groups.identities[index]
            .iter()
            .map(|identity| identity.key.as_str())
            .collect();
        items.extend(std::iter::repeat_n(
            NOSTRUCT_RINCHI_LONGKEY,
            groups.nostruct_counts[index],
        ));
        result.push_str(&items.join(KEY_DELIM_COMP));
        if position < group_count - 1 {
            result.push_str(KEY_DELIM_GROUP);
        }
    }
    // Strip the trailing reserved block when no groups output — the
    // reference's empty-reaction shape.
    if result.len() == empty_len {
        result.truncate(empty_len - KEY_DELIM_BLOCK.len());
    }
    Ok(result)
}

/// The short key body — rinchi_reaction.cpp's rinchi_short_key.
fn rinchi_short_key_of(directionality: Directionality, groups: &Groups) -> Result<String> {
    let mut result = format!(
        "{RINCHI_SHORT_KEY_HEADER}{RINCHI_KEY_VERSION_ID}{KEY_DELIM_BLOCK}{}{HASH_04_EMPTY}",
        direction_code(directionality, groups.reverse_output),
    );
    // All three major hashes first, then all three minors — the
    // reference accumulates two blocks (rinchi_short_key).
    let mut all_major = String::new();
    let mut all_minor = String::new();
    for position in 0..3 {
        let index = groups.output_order[position];
        let layers = layers_of(&groups.identities[index])?;
        all_major.push(KEY_DELIM_BLOCK.as_bytes()[0] as char);
        all_major.push_str(&layers.major_hash());
        all_minor.push(KEY_DELIM_BLOCK.as_bytes()[0] as char);
        all_minor.push_str(&layers.minor_hash());
    }
    result.push_str(&all_major);
    result.push_str(&all_minor);
    result.push(KEY_DELIM_BLOCK.as_bytes()[0] as char);
    for position in 0..3 {
        result.push(nostruct_count_char(
            groups.nostruct_counts[groups.output_order[position]],
        ));
    }
    Ok(result)
}

/// The web key body — rinchi_reaction.cpp's rinchi_web_key.
fn rinchi_web_key_of(groups: &Groups) -> Result<String> {
    // The unique component InChIs across all roles — No-Structures
    // included (the reference walks every component list) — sorted.
    let mut inchis: Vec<String> = Vec::new();
    for index in 0..3 {
        for identity in &groups.identities[index] {
            inchis.push(identity.inchi.clone());
        }
        for _ in 0..groups.nostruct_counts[index] {
            inchis.push(NOSTRUCT_INCHI.to_string());
        }
    }
    inchis.sort();
    inchis.dedup();
    let mut layers = Layers::default();
    for inchi in &inchis {
        layers.append(inchi)?;
    }
    Ok(format!(
        "{RINCHI_WEB_KEY_HEADER}{}{KEY_DELIM_BLOCK}{}{RINCHI_KEY_VERSION_ID}",
        layers.major_hash_ext(),
        layers.minor_hash_ext(),
    ))
}

fn direction_flag(directionality: Directionality, reverse_output: bool) -> String {
    let flag = match directionality {
        Directionality::Directional if reverse_output => DIRECTION_REVERSE,
        Directionality::Directional => DIRECTION_FORWARD,
        Directionality::Equilibrium => DIRECTION_EQUILIBRIUM,
    };
    format!("{DIRECTION_TAG}{flag}")
}

fn direction_code(directionality: Directionality, reverse_output: bool) -> char {
    match directionality {
        Directionality::Directional if reverse_output => 'B',
        Directionality::Directional => 'F',
        Directionality::Equilibrium => 'E',
    }
}

/// The No-Structure count letter for the short key: 0 reads `Z`,
/// counts 1..24 read `A`..`X`, anything beyond caps at `Y`.
fn nostruct_count_char(count: usize) -> char {
    match count {
        0 => 'Z',
        1..=24 => char::from_u32(b'A' as u32 + count as u32 - 1).expect("count letter"),
        _ => 'Y',
    }
}

/// The aggregated proton-count letter (the InChIKey proton flag
/// arithmetic over a reaction sum): `N` at zero, ±1..±12 as
/// `O`..`Z` / `M`..`B`, `A` beyond.
fn proton_count_char(count: i32) -> char {
    if !(-12..=12).contains(&count) {
        return 'A';
    }
    char::from_u32((b'N' as i32 + count) as u32).unwrap_or('A')
}

/// The per-role layer split — rinchi_reaction.cpp's InChILayers:
/// each component's identifier splits into major layers (formula,
/// then `c`/`h`/`q` while they last), a proton sum (`/p` values
/// aggregate), and everything from the first other layer on is
/// minor; components join with `!` inside each string.
#[derive(Default)]
struct Layers {
    majors: String,
    minors: String,
    proton_count: i32,
}

impl Layers {
    fn append(&mut self, inchi: &str) -> Result<()> {
        if inchi.is_empty() {
            return Ok(());
        }
        if !inchi.is_ascii() {
            return err(format!("{CONTEXT}: non-ASCII identifier {inchi:?}"));
        }
        let bytes = inchi.as_bytes();
        let Some(delim_pos) = inchi.find(DELIM_LAYER) else {
            return err(format!(
                "{CONTEXT}: invalid InChI string - no layers: {inchi:?}"
            ));
        };
        if delim_pos != 8 {
            return err(format!(
                "{CONTEXT}: invalid InChI string - no layers: {inchi:?}"
            ));
        }
        if bytes[delim_pos - 1] != b'S' {
            return err(format!(
                "{CONTEXT}: only standard InChIs are supported: {inchi:?}"
            ));
        }
        if bytes[delim_pos - 2] != b'1' {
            return err(format!(
                "{CONTEXT}: only InChI version 1 supported: {inchi:?}"
            ));
        }

        let mut token_start = delim_pos + 1;
        let mut major_layers = String::new();
        let mut minor_layers = String::new();
        let mut first_layer = true;
        let mut in_major = true;
        loop {
            let end = inchi[token_start..]
                .find(DELIM_LAYER)
                .map(|position| position + token_start);
            // Each layer token carries its own leading slash.
            let layer = match end {
                Some(position) => {
                    let slice = &inchi[token_start - 1..position];
                    token_start = position + 1;
                    slice
                }
                None => &inchi[token_start - 1..],
            };
            if first_layer {
                // The formula leads — always present.
                major_layers.push_str(layer);
            } else if layer.len() >= 2 {
                // Shorter tokens are the bare-slash empty layers of
                // e.g. "InChI=1S//"; they are ignored.
                if in_major {
                    match layer.as_bytes()[1] {
                        b'c' | b'h' | b'q' => major_layers.push_str(layer),
                        b'p' => {
                            self.proton_count += layer[2..]
                                .parse::<i32>()
                                .map_err(|_| bad(&format!("unparseable proton layer {layer:?}")))?;
                        }
                        _ => {
                            minor_layers.push_str(layer);
                            in_major = false;
                        }
                    }
                } else {
                    minor_layers.push_str(layer);
                }
            }
            first_layer = false;
            if end.is_none() {
                break;
            }
        }

        // Drop the leading slash the formula token carried; an empty
        // major side reads as a bare delimiter.
        if !major_layers.is_empty() {
            major_layers.remove(0);
        }
        if major_layers.is_empty() {
            major_layers.push(DELIM_LAYER);
        }
        if !minor_layers.is_empty() {
            minor_layers.remove(0);
        }

        if !self.majors.is_empty() {
            self.majors.push(DELIM_COMP);
        }
        self.majors.push_str(&major_layers);
        if !self.minors.is_empty() {
            self.minors.push(DELIM_COMP);
        }
        self.minors.push_str(&minor_layers);
        Ok(())
    }

    /// The major hash: 10-char truncated SHA-256 over the major
    /// layers.
    fn major_hash(&self) -> String {
        hash10char(&self.majors)
    }

    /// The major hash extended: 17-char truncated SHA-256 (the web
    /// key's shape).
    fn major_hash_ext(&self) -> String {
        hash17char(&self.majors)
    }

    /// The minor hash: the proton letter over a 4-char truncated
    /// SHA-256 of the minor layers.
    fn minor_hash(&self) -> String {
        format!(
            "{}{}",
            proton_count_char(self.proton_count),
            hash04char(&self.minors)
        )
    }

    /// The minor hash extended: the proton letter over a 12-char
    /// truncated SHA-256.
    fn minor_hash_ext(&self) -> String {
        format!(
            "{}{}",
            proton_count_char(self.proton_count),
            hash12char(&self.minors)
        )
    }
}

fn layers_of(identities: &[Identity]) -> Result<Layers> {
    let mut layers = Layers::default();
    for identity in identities {
        layers.append(&identity.inchi)?;
    }
    Ok(layers)
}

// The truncated-SHA-256 shapes — rinchi_hashing.cpp over the
// InChIKey base-26 tables (see crate::key): n-char prefixes of the
// 14-bit triplet / 9-bit dublet bit stream.

fn hash04char(input: &str) -> String {
    let digest = key::sha256(input.as_bytes());
    let mut result = key::triplet(&digest, 0);
    result.push_str(&key::triplet(&digest, 14));
    result.truncate(4);
    result
}

fn hash10char(input: &str) -> String {
    let mut result = hash12char(input);
    result.truncate(10);
    result
}

fn hash12char(input: &str) -> String {
    let digest = key::sha256(input.as_bytes());
    let mut result = key::triplet(&digest, 0);
    for bit in [14, 28, 42] {
        result.push_str(&key::triplet(&digest, bit));
    }
    result
}

fn hash17char(input: &str) -> String {
    let digest = key::sha256(input.as_bytes());
    let mut result = hash12char(input);
    result.push_str(&key::dublet(&digest, 56));
    // Bits 64.. — the reference reads a fresh triplet at byte 8.
    result.push_str(&key::triplet(&digest, 64));
    result
}

/// Parses an MDL V2000 RXNfile: the `$RXN` header, the three
/// title lines, the 6- or 9-character component count line, then
/// one `$MOL`-delimited molfile block per component — reactants
/// first, then products, then agents. Each block is read as a whole
/// molecule; a block holding no atoms, or a single placeholder atom,
/// parses as a No-Structure.
pub fn parse_rxn(input: &str) -> Result<RinchiReaction> {
    let lines: Vec<&str> = input.lines().collect();
    if lines.first().copied().unwrap_or("") != MDL_TAG_RXN_BEGIN {
        return err(format!(
            "{CONTEXT}: RXN files must begin with a '{MDL_TAG_RXN_BEGIN}' line"
        ));
    }
    // Skip the name, source-and-timestamp, and comment lines.
    let counts_line = lines.get(4).copied().unwrap_or("");
    let trimmed = counts_line.trim_end_matches(' ');
    if trimmed.len() != 6 && trimmed.len() != 9 {
        return err(format!(
            "{CONTEXT}: invalid component count line - must be 6 or 9 characters long: {counts_line:?}"
        ));
    }
    let parse_count = |field: &str| -> Result<usize> {
        field
            .trim()
            .parse::<usize>()
            .map_err(|_| bad(&format!("malformed component count line: {counts_line:?}")))
    };
    let reactant_count = parse_count(&trimmed[0..3])?;
    let product_count = parse_count(&trimmed[3..6])?;
    let agent_count = if trimmed.len() == 9 {
        parse_count(&trimmed[6..9])?
    } else {
        0
    };

    let mut cursor = 5;
    let mut reactants = Vec::with_capacity(reactant_count);
    for _ in 0..reactant_count {
        reactants.push(read_component(&lines, &mut cursor)?);
    }
    let mut products = Vec::with_capacity(product_count);
    for _ in 0..product_count {
        products.push(read_component(&lines, &mut cursor)?);
    }
    let mut agents = Vec::with_capacity(agent_count);
    for _ in 0..agent_count {
        agents.push(read_component(&lines, &mut cursor)?);
    }
    Ok(RinchiReaction {
        reactants,
        products,
        agents,
        directionality: Directionality::Directional,
    })
}

/// Reads one `$MOL`-delimited block: the tag line, then everything
/// through the block's `M  END` — the component's molfile text.
fn read_component(lines: &[&str], cursor: &mut usize) -> Result<RinchiComponent> {
    let tag = lines.get(*cursor).copied().unwrap_or("");
    if tag != MDL_TAG_RXN_COMPONENT_START {
        return err(format!(
            "{CONTEXT}: reaction components must be delimited by a '{MDL_TAG_RXN_COMPONENT_START}' line"
        ));
    }
    *cursor += 1;
    let mut block: Vec<&str> = Vec::new();
    loop {
        let line = lines
            .get(*cursor)
            .ok_or_else(|| bad("truncated $MOL block - no M  END"))?;
        *cursor += 1;
        let ends = *line == MDL_TAG_MOLFILE_END;
        block.push(line);
        if ends {
            break;
        }
    }
    Ok(RinchiComponent {
        molfile: block.join("\n"),
    })
}

fn bad(message: &str) -> crate::error::Error {
    crate::error::Error::new(format!("{CONTEXT}: {message}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    fn rxn(reactants: &[&str], products: &[&str], agents: &[&str]) -> RinchiReaction {
        let components = |blocks: &[&str]| -> Vec<RinchiComponent> {
            blocks
                .iter()
                .map(|molfile| RinchiComponent {
                    molfile: (*molfile).to_string(),
                })
                .collect()
        };
        RinchiReaction {
            reactants: components(reactants),
            products: components(products),
            agents: components(agents),
            directionality: Directionality::Directional,
        }
    }

    /// A groups cache built straight from recorded identities — the
    /// oracle's per-component strings, bypassing the molfile edge.
    fn groups_of(
        reactants: Vec<Identity>,
        products: Vec<Identity>,
        agents: Vec<Identity>,
    ) -> Groups {
        let mut groups = Groups {
            rinchi_groups: Default::default(),
            rauxinfo_groups: Default::default(),
            identities: [reactants, products, agents],
            nostruct_counts: [0; 3],
            output_order: [0, 1, 2],
            reverse_output: false,
        };
        for index in 0..3 {
            let mut sorted = groups.identities[index].clone();
            sorted.sort_by(|a, b| a.inchi.cmp(&b.inchi));
            let mut rinchi_group = String::new();
            let mut rauxinfo_group = String::new();
            for (position, identity) in sorted.iter().enumerate() {
                if position > 0 {
                    rinchi_group.push(DELIM_COMP);
                    rauxinfo_group.push(DELIM_COMP);
                }
                rinchi_group.push_str(identity.inchi.strip_prefix(INCHI_STD_HEADER).expect("std"));
                rauxinfo_group.push_str(
                    identity
                        .aux
                        .strip_prefix(INCHI_AUXINFO_HEADER)
                        .expect("aux"),
                );
            }
            groups.rinchi_groups[index] = rinchi_group;
            groups.rauxinfo_groups[index] = rauxinfo_group;
            groups.identities[index] = sorted;
        }
        groups.reverse_output = groups.rinchi_groups[1] < groups.rinchi_groups[0];
        if groups.reverse_output {
            groups.output_order.swap(0, 1);
        }
        groups
    }

    #[test]
    fn the_empty_hash_constants_come_from_the_native_hasher() {
        // rinchi_hashing.cpp's empty-string constants.
        assert_eq!(hash04char(""), HASH_04_EMPTY);
        assert_eq!(hash10char(""), "UHFFFADPSC");
        assert_eq!(hash12char(""), HASH_12_EMPTY);
        assert_eq!(hash17char(""), "UHFFFADPSCTJAUYIS");
    }

    #[test]
    fn the_proton_letters_follow_the_key_arithmetic() {
        assert_eq!(proton_count_char(0), 'N');
        assert_eq!(proton_count_char(1), 'O');
        assert_eq!(proton_count_char(12), 'Z');
        assert_eq!(proton_count_char(-1), 'M');
        assert_eq!(proton_count_char(-12), 'B');
        assert_eq!(proton_count_char(13), 'A');
        assert_eq!(proton_count_char(-13), 'A');
    }

    #[test]
    fn the_nostruct_count_letters_run_z_a_to_x() {
        assert_eq!(nostruct_count_char(0), 'Z');
        assert_eq!(nostruct_count_char(1), 'A');
        assert_eq!(nostruct_count_char(24), 'X');
        assert_eq!(nostruct_count_char(25), 'Y');
        assert_eq!(nostruct_count_char(100), 'Y');
    }

    #[test]
    fn layers_split_into_major_minor_and_proton_sum() {
        let mut layers = Layers::default();
        layers.append("InChI=1S/H2O/h1H2").expect("layers");
        layers
            .append("InChI=1S/C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3/t3-,4?/m0/s1")
            .expect("layers");
        layers
            .append("InChI=1S/Na.H2O/h;1H2/q+1;/p-1")
            .expect("layers");
        assert_eq!(
            layers.majors,
            "H2O/h1H2!C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3!Na.H2O/h;1H2/q+1;"
        );
        assert_eq!(layers.minors, "t3-,4?/m0/s1!");
        assert_eq!(layers.proton_count, -1);
    }

    #[test]
    fn water_to_water_runs_forward() {
        let reaction = rxn(&[WATER], &[WATER], &[]);
        assert_eq!(
            reaction.rinchi_string().expect("rinchi"),
            "RInChI=1.00.1S/H2O/h1H2<>H2O/h1H2/d+"
        );
    }

    /// The reference distribution's own sample reaction (R005a.rxn)
    /// recorded from the built oracle (`rinchi_cmdline`): an epoxide
    /// opening to tetrahydrofuran with NaOH. The component identities
    /// and all five output lines are the oracle's recording.
    #[test]
    fn the_assembly_reproduces_the_r005a_oracle_lines() {
        let identity = |inchi: &str, aux: &str| -> Identity {
            Identity {
                inchi: inchi.to_string(),
                aux: aux.to_string(),
                key: key::inchikey(inchi).expect("key"),
            }
        };
        let reactants = vec![
            identity(
                "InChI=1S/C4H9BrO/c1-3(5)4(2)6/h3-4,6H,1-2H3/t3-,4+/m1/s1",
                "AuxInfo=1/0/N:4,1,3,2,6,5/it:im/rA:6nCCCCOBr/rB:s1;s2;s3;N2;P3;/rC:-.825,-.7557,0;-.4125,-.0412,0;.4125,-.0412,0;.825,.6733,0;-.626,.7557,0;.825,-.7557,0;",
            ),
            identity(
                "InChI=1S/Na.H2O/h;1H2/q+1;/p-1",
                "AuxInfo=1/1/N:1;2/rA:2nNaO/rB:s1;/rC:-.4125,0,0;.4125,0,0;",
            ),
        ];
        let products = vec![identity(
            "InChI=1S/C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3/t3-,4?/m0/s1",
            "AuxInfo=1/0/N:4,1,3,2,5/E:(1,2)(3,4)/it:im/rA:5nCCCCO/rB:N1;s2;P3;s2s3;/rC:-1.127,-.5635,0;-.4125,-.151,0;.4125,-.151,0;1.127,-.5635,0;0,.5635,0;",
        )];
        let groups = groups_of(reactants, products, Vec::new());
        assert!(groups.reverse_output, "products group sorts first");
        assert_eq!(
            rinchi_string_of(Directionality::Directional, &groups),
            "RInChI=1.00.1S/C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3/t3-,4?/m0/s1<>C4H9BrO/c1-3(5)4(2)6/h3-4,6H,1-2H3/t3-,4+/m1/s1!Na.H2O/h;1H2/q+1;/p-1/d-"
        );
        assert_eq!(
            rinchi_auxinfo_of(&groups),
            "RAuxInfo=1.00.1/0/N:4,1,3,2,5/E:(1,2)(3,4)/it:im/rA:5nCCCCO/rB:N1;s2;P3;s2s3;/rC:-1.127,-.5635,0;-.4125,-.151,0;.4125,-.151,0;1.127,-.5635,0;0,.5635,0;<>0/N:4,1,3,2,6,5/it:im/rA:6nCCCCOBr/rB:s1;s2;s3;N2;P3;/rC:-.825,-.7557,0;-.4125,-.0412,0;.4125,-.0412,0;.825,.6733,0;-.626,.7557,0;.825,-.7557,0;!1/N:1;2/rA:2nNaO/rB:s1;/rC:-.4125,0,0;.4125,0,0;"
        );
        assert_eq!(
            rinchi_long_key_of(Directionality::Directional, &groups).expect("long key"),
            "Long-RInChIKey=SA-BUHFF-PQXKWPLDPFFDJP-WUCPZUCCSA-N--JCYSVJNMXBWPHS-DMTCNVIQSA-N-HEMHJVSKTPXQMS-UHFFFAOYSA-M"
        );
        assert_eq!(
            rinchi_short_key_of(Directionality::Directional, &groups).expect("short key"),
            "Short-RInChIKey=SA-BUHFF-PQXKWPLDPF-LSZRVNUIRJ-UHFFFADPSC-NIBFB-MCQMH-NUHFF-ZZZ"
        );
        assert_eq!(
            rinchi_web_key_of(&groups).expect("web key"),
            "Web-RInChIKey=UFEYRRLMNYPCDUWUI-MSZCBGVVHVRIXSA"
        );
    }

    #[test]
    fn a_placeholder_component_reads_as_a_no_structure() {
        let placeholder = "unknown\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
        let reaction = rxn(&[WATER], &[placeholder], &[]);
        // Recorded from the built oracle: the no-structure product
        // group contributes no InChI body, the /u counts carry it,
        // and the empty products group sorts before the reactants
        // group — the reaction reads reversed.
        assert_eq!(
            reaction.rinchi_string().expect("rinchi"),
            "RInChI=1.00.1S/<>H2O/h1H2/d-/u1-0-0"
        );
        assert_eq!(
            reaction.rinchi_auxinfo().expect("rauxinfo"),
            "RAuxInfo=1.00.1/<>0/N:1/rA:1nO/rB:/rC:12,12,0;"
        );
        assert_eq!(
            reaction.rinchi_long_key().expect("long key"),
            "Long-RInChIKey=SA-BUHFF-MOSFIJXAXDLOML-UHFFFAOYSA-N--XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        );
        assert_eq!(
            reaction.rinchi_short_key().expect("short key"),
            "Short-RInChIKey=SA-BUHFF-UHFFFADPSC-XLYOFNOQVP-UHFFFADPSC-NUHFF-NUHFF-NUHFF-AZZ"
        );
        assert_eq!(
            reaction.rinchi_web_key().expect("web key"),
            "Web-RInChIKey=STNSBOCECDZATIDBV-NUHFFFADPSCTJSA"
        );
    }

    #[test]
    fn a_rxnfile_parses_into_roles() {
        let block = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END";
        let input = format!("$RXN\n\n  chemomer\n\n  1  1\n$MOL\n{block}\n$MOL\n{block}\n");
        let reaction = parse_rxn(&input).expect("parses");
        assert_eq!(reaction.reactants().len(), 1);
        assert_eq!(reaction.products().len(), 1);
        assert_eq!(reaction.agents().len(), 0);
        assert_eq!(reaction.directionality(), Directionality::Directional);
        assert_eq!(
            reaction.rinchi_string().expect("rinchi"),
            "RInChI=1.00.1S/H2O/h1H2<>H2O/h1H2/d+"
        );
    }

    #[test]
    fn a_missing_rxn_tag_is_rejected_not_guessed() {
        assert!(parse_rxn("not a reaction\n").is_err());
    }
}
