//! The polymer `/z` layer: the repeat-unit records the reference
//! writes after the charge/proton layers (ichiprt1.c
//! OutputINCHI_PolymerLayer 3880-4333) — per unit, the header digits,
//! the canonical-numbered atom list, and either the crossing bonds or
//! the frame-shift pattern, units sorted by their canonical atom
//! lists (longer first). The module also carries the aux `/Z:` echo
//! (OutputAUXINFO_PolymerInfo 4923) and the legacy Zz-hide text edit
//! (`EditINCHI_HidePolymerZz` 5039) the v1.05 modes apply.

use crate::canonical::polymer::{self, AtomProps, Closing, Treated, TreatedUnit};
use crate::error::{err, Result};
use crate::graph::{Polymer, PolymerUnit, Structure};
use crate::layers::hydrogen::canonical_ranks_for;
use crate::layers::order;
use crate::options::{Options, PolymersMode};

/// The layer over a submitted structure and its treated polymer. The
/// units speak the submitted numbering; the emit translates through
/// the final canonical numbering of the structure the layers read.
pub struct PolymerLayer {
    original: Structure,
    original_polymer: Polymer,
    treated: Treated,
    options: Options,
    /// Whether this layer speaks for an interim (105+) pass — the
    /// senior-first backbone reorder applies there too.
    pub interim_mode: bool,
}

impl PolymerLayer {
    pub fn new(
        original: Structure,
        original_polymer: Polymer,
        treated: Treated,
        options: Options,
    ) -> Self {
        PolymerLayer {
            original,
            original_polymer,
            treated,
            options,
            interim_mode: false,
        }
    }

    /// The submitted-numbering -> 1-based canonical-number map the
    /// layers' structure speaks (the aux `/N:` recipe: atoms in
    /// canonical order, numbered 1.., original indices through the
    /// fold source).
    pub fn canonical_numbering(working: &Structure) -> std::collections::BTreeMap<usize, usize> {
        let ranks = canonical_ranks_for(working);
        let mut order: Vec<usize> = (0..working.atoms.len()).collect();
        order.sort_by_key(|&atom| ranks[atom]);
        let mut map = std::collections::BTreeMap::new();
        for (cano, atom) in order.into_iter().enumerate() {
            let original = if working.fold_source.is_empty() {
                atom
            } else {
                working.fold_source[atom]
            };
            map.insert(original, cano + 1);
        }
        map
    }

    /// The whole-structure numbering over sorted sibling components:
    /// each sibling's canonical order continues the global counter
    /// (the reference's InternallyGetCanoNumsAndComponentNums walks
    /// the concatenated `/N:` list).
    pub fn combined_numbering(siblings: &[&Structure]) -> std::collections::BTreeMap<usize, usize> {
        let mut map = std::collections::BTreeMap::new();
        let mut cano = 0;
        for sibling in siblings {
            let ranks = canonical_ranks_for(sibling);
            let mut order: Vec<usize> = (0..sibling.atoms.len()).collect();
            order.sort_by_key(|&atom| ranks[atom]);
            for atom in order {
                let original = if sibling.fold_source.is_empty() {
                    atom
                } else {
                    sibling.fold_source[atom]
                };
                map.insert(original, cano + 1);
                cano += 1;
            }
        }
        map
    }

    /// The aux `/Z:` echo: per unit in declaration order, the header
    /// digits, the SMT subscript (or `n`), the ORIGINAL-numbered atom
    /// list, and the crossing bonds' endpoint atoms.
    pub fn aux_z_block(&self) -> String {
        let mut out = String::new();
        let any_structure = self
            .treated
            .units
            .iter()
            .any(|unit| unit.end1 != usize::MAX || !unit.bkbonds.is_empty());
        let units = &self.original_polymer.units;
        for (index, unit) in units.iter().enumerate() {
            let crossing = self
                .treated
                .units
                .get(index)
                .map(|state| state.end1 != usize::MAX)
                .unwrap_or(false);
            let (ty, sub, conn) = polymer::header_digits(unit, crossing, any_structure);
            out.push_str(&format!("{ty}{sub}{conn}-"));
            out.push_str(if unit.smt.is_empty() { "n" } else { &unit.smt });
            out.push('-');
            push_compressed_range(
                &mut out,
                &unit.alist.iter().map(|a| a + 1).collect::<Vec<_>>(),
            );
            let echo_pairs = self
                .treated
                .units
                .get(index)
                .map(|state| state.echo_pairs.clone())
                .unwrap_or_default();
            if !echo_pairs.is_empty() {
                out.push('(');
                for (position, &(a, b)) in echo_pairs.iter().enumerate() {
                    if position > 0 {
                        out.push(',');
                    }
                    out.push_str(&format!("{a},{b}"));
                }
                out.push(')');
            }
            if index + 1 < units.len() {
                out.push(';');
            }
        }
        out
    }

    /// The `/z` segment over a numbering map (submitted 0-based index
    /// -> 1-based canonical number).
    pub fn emit_segment(
        &self,
        numbering: &std::collections::BTreeMap<usize, usize>,
    ) -> Result<Option<String>> {
        if !self.treated.has_units() {
            return Ok(None);
        }
        let props = self.canonical_props(numbering);
        let working = self.working_set(numbering)?;
        if working.is_empty() {
            return Ok(None);
        }
        // Units in modified-alist lexicographic order (longer first).
        let mut order: Vec<usize> = (0..working.len()).collect();
        order.sort_by(|&a, &b| compare_atom_lists_mod(&working[a].alist, &working[b].alist));

        let star_total = self.treated.pzz.len();
        let mut used_stars = 0usize;
        let n_expl_h = self
            .original
            .atoms
            .iter()
            .filter(|atom| atom.element == "H")
            .count();
        let num_atoms = self.original.atoms.len();
        let legacy = self.options.polymers == PolymersMode::Legacy;
        let mut out = String::new();
        for (position, index) in order.iter().enumerate() {
            let unit = &working[*index];
            out.push_str(&format!(
                "{}{}{}-",
                unit.header.0, unit.header.1, unit.header.2
            ));
            push_compressed_range(&mut out, &unit.alist);
            let crossing =
                unit.nb == 2 && (unit.cyclizable == Closing::NotApplicable || !unit.cyclized);
            if crossing {
                // "(cap-end,cap-end)", the first bond pointing to the
                // more senior end under the modern rule, the higher
                // cap number under the legacy one.
                let pairs = unit.blist_pairs();
                let (a1, a2, a3, a4) = (pairs[0].0, pairs[0].1, pairs[1].0, pairs[1].1);
                let swap = if legacy {
                    a3 < a1
                } else {
                    polymer::is_first_rank_lower(a2, a4, &props) == 1
                };
                if swap {
                    out.push_str(&format!("({}-{},{}-{})", a3, a4, a1, a2));
                } else {
                    out.push_str(&format!("({}-{},{}-{})", a1, a2, a3, a4));
                }
            } else if unit.nb <= 2
                && (unit.cyclizable != Closing::NotApplicable || !unit.bkbonds.is_empty())
            {
                // The frame-shift pattern "(cap1,cap2-...)". The caps
                // speak the star pool: the last star positions of the
                // canonical numbering.
                let mut a1 = 0usize;
                let mut a2 = 0usize;
                if unit.cap1 > 0 || unit.cap2 > 0 {
                    if unit.cap1 > 0 {
                        a1 = num_atoms - n_expl_h - star_total + used_stars + 1;
                        used_stars += 1;
                    }
                    if unit.cap2 > 0 {
                        a2 = num_atoms - n_expl_h - star_total + used_stars + 1;
                        used_stars += 1;
                    }
                }
                out.push_str(&format!("({a1},{a2}-"));
                match unit.cyclizable {
                    Closing::Diradical => {
                        out.push_str(&format!("{})", unit.end1));
                    }
                    Closing::HigherOrder => {
                        let (lo, hi) = (unit.end1.min(unit.end2), unit.end1.max(unit.end2));
                        out.push_str(&format!("{lo}.{hi})"));
                    }
                    Closing::Ring => {
                        if unit.bkbonds.is_empty() {
                            let (lo, hi) = (unit.end1.min(unit.end2), unit.end1.max(unit.end2));
                            out.push_str(&format!("{lo},{hi})"));
                        } else {
                            // Min-atom-order sort, then — modern and
                            // 105+ treatments under an active frame
                            // shift — the senior bond to the front.
                            let mut bonds = unit.bkbonds.clone();
                            bonds.sort_by(|x, y| {
                                let xw = (x.0.min(x.1), x.0.max(x.1));
                                let yw = (y.0.min(y.1), y.0.max(y.1));
                                xw.cmp(&yw)
                            });
                            let reorder = matches!(self.options.polymers, PolymersMode::Modern)
                                || self.interim_mode;
                            if reorder
                                && self.options.frame_shift_edits()
                                && unit.cap1 >= 1
                                && unit.cap2 >= 1
                            {
                                bonds = polymer::sort_backbone_bonds(&bonds, &props);
                            }
                            for (k, bond) in bonds.iter().enumerate() {
                                out.push_str(&format!(
                                    "{},{}{}",
                                    bond.0.min(bond.1),
                                    bond.0.max(bond.1),
                                    if k + 1 == bonds.len() { ")" } else { "," }
                                ));
                            }
                        }
                    }
                    Closing::NotApplicable => {}
                }
            }
            if position + 1 < order.len() {
                out.push(';');
            }
        }
        Ok(Some(format!("z{out}")))
    }

    /// The seniority props keyed by 1-based canonical number (the
    /// ring perception runs over the submitted structure; the eranks
    /// over its elements).
    fn canonical_props(
        &self,
        numbering: &std::collections::BTreeMap<usize, usize>,
    ) -> Vec<AtomProps> {
        let raw = polymer::atom_props(&self.original, &self.treated);
        let mut keyed = vec![AtomProps::default(); numbering.len() + 1];
        for (original, &cano) in numbering {
            keyed[cano] = raw[*original];
        }
        keyed
    }

    /// The working copy of the units in canonical numbers
    /// (`OAD_Polymer_PrepareWorkingSet`): alists translated (removed
    /// atoms dropped) and sorted, blists cap-first and bond-sorted,
    /// caps and ends translated.
    fn working_set(
        &self,
        numbering: &std::collections::BTreeMap<usize, usize>,
    ) -> Result<Vec<WorkingUnit>> {
        let mut out = Vec::new();
        let any_structure_marker = self
            .treated
            .units
            .iter()
            .any(|unit| unit.end1 != usize::MAX || !unit.bkbonds.is_empty());
        for (index, unit) in self.original_polymer.units.iter().enumerate() {
            let state = &self.treated.units[index];
            let cano = |atom: usize| -> Option<usize> { numbering.get(&atom).copied() };
            let mut working = WorkingUnit {
                header: polymer::header_digits(
                    unit,
                    state.end1 != usize::MAX,
                    any_structure_marker,
                ),
                alist: unit.alist.iter().filter_map(|&a| cano(a)).collect(),
                nb: unit.blist.len(),
                blist: Vec::new(),
                cyclizable: state.cyclizable,
                cyclized: state.cyclized,
                cap1: state.cap1,
                cap2: state.cap2,
                end1: state.end1,
                end2: state.end2,
                bkbonds: state.bkbonds.clone(),
            };
            working.alist.sort_unstable();
            if working.nb == 2 {
                // The crossing pairs speak cap-first, ordered by the
                // caps (the reference's OrderBondAtomsAndBonds... on
                // the treated ends).
                let (Some(ccap1), Some(cend1)) = (cano(state.cap1), cano(state.end1)) else {
                    return err(
                        "Can not proceed as one of PU crossing bond ends leads to explicit H",
                    );
                };
                let (Some(ccap2), Some(cend2)) = (cano(state.cap2), cano(state.end2)) else {
                    return err(
                        "Can not proceed as one of PU crossing bond ends leads to explicit H",
                    );
                };
                let mut pairs = vec![(ccap1, cend1), (ccap2, cend2)];
                if pairs[0].0 > pairs[1].0 {
                    pairs.swap(0, 1);
                }
                working.blist = pairs.into_iter().flat_map(|(a, b)| [a, b]).collect();
                working.cap1 = ccap1;
                working.cap2 = ccap2;
                working.end1 = cend1;
                working.end2 = cend2;
                working.bkbonds = state
                    .bkbonds
                    .iter()
                    .filter_map(|&(a, b)| match (cano(a), cano(b)) {
                        (Some(ca), Some(cb)) => Some((ca.min(cb), ca.max(cb))),
                        _ => None,
                    })
                    .collect();
            }
            out.push(working);
        }
        Ok(out)
    }
}

/// One translated unit the emitter prints.
struct WorkingUnit {
    header: (u8, u8, u8),
    alist: Vec<usize>,
    nb: usize,
    cyclizable: Closing,
    cyclized: bool,
    cap1: usize,
    cap2: usize,
    end1: usize,
    end2: usize,
    bkbonds: Vec<(usize, usize)>,
    blist: Vec<usize>,
}

impl WorkingUnit {
    fn blist_pairs(&self) -> Vec<(usize, usize)> {
        self.blist.chunks(2).map(|p| (p[0], p[1])).collect()
    }
}

/// `OAD_PolymerUnit_CompareAtomListsMod`: longer lists first, then
/// element-wise lexicographic.
fn compare_atom_lists_mod(a: &[usize], b: &[usize]) -> std::cmp::Ordering {
    a.len().cmp(&b.len()).reverse().then_with(|| a.cmp(b))
}

/// `print_sequence_of_nums_compressing_ranges`: runs of consecutive
/// numbers read `first-last`, the rest comma-joined.
fn push_compressed_range(out: &mut String, nums: &[usize]) {
    if nums.is_empty() {
        return;
    }
    let mut range = 0usize;
    for k in 0..nums.len() - 1 {
        if nums[k + 1] == nums[k] + 1 {
            if range > 0 {
                range += 1;
            } else {
                out.push_str(&format!("{}-", nums[k]));
                range = 1;
            }
        } else {
            range = 0;
            out.push_str(&format!("{},", nums[k]));
        }
    }
    out.push_str(&nums[nums.len() - 1].to_string());
}

impl crate::layers::LayerEmitter for PolymerLayer {
    fn name(&self) -> &'static str {
        "polymer"
    }

    fn order(&self) -> u8 {
        order::Z
    }

    fn segment_kind(&self) -> Option<char> {
        Some('z')
    }

    fn emit(&self, working: &Structure) -> Result<Option<String>> {
        if !self.treated.has_units() {
            return Ok(None);
        }
        let numbering = Self::canonical_numbering(working);
        self.emit_segment(&numbering)
    }
}

/// The v1.05 Zz-hide post-edit (`EditINCHI_HidePolymerZz`), applied
/// semantically to the assembled identifier and aux lines: the
/// pure-Zz formula terms drop, the detached star components' trailing
/// empty fields drop (never after `p`/`s`/`f`/`z` layers), and the
/// frame-shift patterns lose their `(cap,cap-` prefixes.
pub fn hide_polymer_zz(identifier: &str, aux: &str, n_pzz: usize) -> (String, String) {
    if !identifier.contains("/z") {
        return (identifier.to_string(), aux.to_string());
    }
    (hide_line(identifier, n_pzz), hide_line(aux, n_pzz))
}

/// The aux line's share of the hide post-edit — the caller (the aux
/// assembly) already knows the /z layer stands.
pub fn hide_aux_line(aux: &str, n_pzz: usize) -> String {
    hide_line(aux, n_pzz)
}

fn hide_line(line: &str, n_pzz: usize) -> String {
    let (prefix, body) = match line.split_once('/') {
        Some((prefix, body)) => (prefix, body),
        None => return line.to_string(),
    };
    let mut segments: Vec<String> = body.split('/').map(str::to_string).collect();
    for (index, segment) in segments.iter_mut().enumerate() {
        if index == 0 {
            // The formula: drop the pure-Zz component terms.
            *segment = eat_zz_formula_term(segment);
            continue;
        }
        match segment.chars().next() {
            Some('z') => {
                *segment = hide_z_pattern(segment);
            }
            Some('m') => {
                *segment = eat_trailing(segment, '.', n_pzz, 1);
            }
            Some('p') | Some('s') | Some('f') => {}
            Some(_) | None => {
                *segment = eat_trailing(segment, ';', n_pzz, 0);
            }
        }
    }
    format!("{prefix}/{}", segments.join("/"))
}

/// The /z layer's pattern cut: "(cap,cap-rest" reads "(rest" — the
/// leading number pair before the first '-' vanishes (the reference's
/// skip-until-'-' state machine).
fn hide_z_pattern(segment: &str) -> String {
    let mut out = String::new();
    let mut rest = segment;
    while let Some(open) = rest.find('(') {
        let after = &rest[open + 1..];
        let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
        let next = after.chars().nth(digits.len());
        if digits.is_empty() || next == Some('-') {
            // "(cap-end,...)": kept verbatim.
            out.push_str(&rest[..open + 1 + digits.len() + 1]);
            rest = &rest[open + 1 + digits.len() + 1..];
            continue;
        }
        // "(cap,cap-...": skip to the '-'.
        out.push_str(&rest[..=open]);
        match after.find('-') {
            Some(dash) => {
                rest = &after[dash + 1..];
            }
            None => {
                rest = "";
                break;
            }
        }
    }
    out.push_str(rest);
    out
}

/// Eats up to `count` trailing delimiter chars (empty fields), never
/// reaching into a non-empty field.
fn eat_trailing(segment: &str, delimiter: char, count: usize, keep: usize) -> String {
    let chars: Vec<char> = segment.chars().collect();
    let mut end = chars.len();
    let mut eaten = 0;
    while end > keep && eaten < count && chars[end - 1] == delimiter {
        end -= 1;
        eaten += 1;
    }
    chars[..end].iter().collect()
}

/// Drops the trailing pure-Zz component term (optionally multiplied).
fn eat_zz_formula_term(formula: &str) -> String {
    let mut terms: Vec<&str> = formula.split('.').collect();
    if let Some(last) = terms.last() {
        let stripped = last.trim_start_matches(|c: char| c.is_ascii_digit());
        let base = stripped.trim_end_matches(|c: char| c.is_ascii_digit());
        if base == "Zz" && stripped.starts_with("Zz") {
            terms.pop();
        }
    }
    terms.join(".")
}

/// The modern frame shift over a polymer structure: for every
/// eligible unit (bistar, two or more backbone bonds), the interim
/// senior backbone bond (most-senior-first, its senior atom first)
/// names the two ends the caps re-attach to, and the senior bond is
/// cut while the original ends join. The seniority runs in the
/// interim canonical numbering (the reference reads the interim `/z`
/// pattern and maps through its AuxInfo `/N:`); `interim_ranks` is
/// the interim structure's canonical rank vector, `interim` shares
/// the units' atom indices with the submitted structure.
pub fn apply_frame_shift(
    structure: &mut Structure,
    treated: &mut Treated,
    interim: &Structure,
    interim_treated: &Treated,
    interim_ranks: &[u32],
) -> bool {
    // The stereo-bond codes ride their BONDS, not their positions:
    // the edits below remove and append bond entries, shifting every
    // later index, so the wedge map re-keys by endpoints once the
    // edits settle (a code whose bond the edit removed goes with it —
    // the caps carry no wedges in anything witnessed).
    let stereo_by_ends: std::collections::BTreeMap<(usize, usize), u8> = structure
        .stereo_bonds
        .iter()
        .filter_map(|(&index, &code)| {
            structure.bonds.get(index).map(|&(a, b, _)| ((a.min(b), a.max(b)), code))
        })
        .collect();
    // Submitted index -> 1-based interim canonical number, and the
    // reverse (ties keep atom order, the aux numbering's stable sort).
    let mut order: Vec<usize> = (0..interim.atoms.len()).collect();
    order.sort_by_key(|&atom| interim_ranks[atom]);
    let mut cano_of = vec![0usize; interim.atoms.len()];
    let mut orig_of_cano = vec![0usize; interim.atoms.len() + 1];
    for (cano, &atom) in order.iter().enumerate() {
        cano_of[atom] = cano + 1;
        let original = if interim.fold_source.is_empty() {
            atom
        } else {
            interim.fold_source[atom]
        };
        orig_of_cano[cano + 1] = original;
    }
    // Seniority props keyed by canonical number.
    let raw = polymer::atom_props(interim, interim_treated);
    let mut props = vec![AtomProps::default(); interim.atoms.len() + 1];
    for (atom, &cano) in cano_of.iter().enumerate() {
        props[cano] = raw[atom];
    }
    let mut applied = false;
    for index in 0..treated.units.len() {
        let state = &interim_treated.units[index];
        // The reference's senior-bond read parses the interim 105+
        // identifier's /z patterns (GetFrameShiftInfoFrom105PlusInChI):
        // only the cyclized ring-comma list "(a,b,..." speaks a senior
        // bond. The crossing "(a-b,...", the higher-order "(a.b)" and
        // the diradical "(a)" forms reject — a ring-closing CRU (its
        // ends adjacent on a pre-existing ring) never settles.
        if !state.cyclized || state.cyclizable != polymer::Closing::Ring || state.bkbonds.len() < 2
        {
            continue;
        }
        let bk_cano: Vec<(usize, usize)> = state
            .bkbonds
            .iter()
            .map(|&(a, b)| {
                let (ca, cb) = (cano_of[a], cano_of[b]);
                (ca.min(cb), ca.max(cb))
            })
            .collect();
        let sorted = polymer::sort_backbone_bonds(&bk_cano, &props);
        let (end1, end2) = (orig_of_cano[sorted[0].0], orig_of_cano[sorted[0].1]);
        let edits = polymer::frame_shift_edits(&treated.units[index], end1, end2);
        if edits.is_empty() {
            continue;
        }
        for edit in &edits {
            polymer::apply_edit(structure, edit);
        }
        let state_now = &mut treated.units[index];
        state_now.end1 = end1;
        state_now.end2 = end2;
        for edit in &edits {
            let (rm_a, rm_b, add_a, add_b) = *edit;
            for pair in &mut state_now.echo_pairs {
                let hit = (*pair == (rm_a + 1, rm_b + 1)) || (*pair == (rm_b + 1, rm_a + 1));
                if hit {
                    *pair = (add_a + 1, add_b + 1);
                }
            }
        }
        applied = true;
    }
    if applied && !stereo_by_ends.is_empty() {
        structure.stereo_bonds = stereo_by_ends
            .into_iter()
            .filter_map(|(ends, code)| {
                structure
                    .bonds
                    .iter()
                    .position(|&(a, b, _)| (a.min(b), a.max(b)) == ends)
                    .map(|index| (index, code))
            })
            .collect();
    }
    applied
}

/// Rebuilds a unit's crossing-bond pairs after the frame shift: the
/// bonds now run cap -> senior end.
pub fn rewire_unit_blist(unit: &mut PolymerUnit, state: &TreatedUnit) {
    if unit.blist.len() < 4 {
        return;
    }
    let targets = [(state.cap1, state.end1), (state.cap2, state.end2)];
    for (cap, end) in targets {
        // Move each bond entry holding this cap to the (possibly new)
        // end partner.
        for chunk in unit.blist.chunks_mut(2) {
            let (a, b) = (chunk[0], chunk[1]);
            if a == cap {
                chunk[1] = end;
            } else if b == cap {
                chunk[0] = end;
            }
        }
    }
}

impl PolymerLayer {
    /// Whether any unit data exists (the aux /Z: echo rides only with
    /// units).
    pub fn aux_has_units(&self) -> bool {
        self.treated.has_units()
    }

    /// The polymer-related star count (the hide's n_pzz).
    pub fn star_count(&self) -> usize {
        self.treated.pzz.len()
    }
}
