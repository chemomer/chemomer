//! The `/f` Fixed-H layer of the full (non-standard) InChI — the
//! segment-difference writer over the mobile-H and fixed-H
//! representations (ichiprt1.c `OutputINCHI_FixedHLayerWithSublayers`,
//! fed by ichimake.c `CompINChILayers`' DIFL_F column). Dev-time
//! source read, never compiled or linked; the corpus arbitrates.
//!
//! Layout after the main layers: `/f` + [formula] + `/h` + `/q` +
//! `/b` `/t` `/m` `/s` + the FI `/i` block. A segment equal to its
//! namesake in the mobile-H representation is omitted (EQL2PRECED →
//! INCHI_SEGM_OMIT, ichimake.h DIFV bits) — EXCEPT the formula, the
//! F layer's first segment, which an active layer promotes to the
//! empty marker (MarkUnusedAndEmptyLayers: an EQL/BOTH_EMPTY first
//! segment becomes IS_EMPTY), rendering the bare `/f`. No `/p`
//! re-enters inside `/f` (the protons segment is mobile-H only,
//! strutil.h DIFS_p comment). The whole block vanishes when no
//! segment differs — MarkUnusedAndEmptyLayers omits the layer when
//! nothing carries a NEQ/IS_EMPTY bit and the FI layer is empty —
//! and non-iso content drives the aux /F: repetition flags.
//!
//! Witnessed shapes: acetate `/fC2H3O2/q-1` (formula one H short, no
//! /h — endpoint deltas all EMPTY), ammonium `/fH4N/h1H/q+1`,
//! hydroxide `/fHO/h1h/q-1` (negative deltas print lowercase `h`,
//! MakeHString), alanine `/f/h5H` (equal formula retracts to bare
//! `/f`), d2o `/f/i1D2` (FI only — no aux /F: block follows).

use crate::canonical::fixedh::FixedH;
use crate::error::Result;
use crate::graph::Structure;
use crate::layers::formula::FormulaLayer;
use crate::layers::isotope::IsotopeLayer;
use crate::layers::standardization;
use crate::layers::stereo::{BLayer, MLayer, SLayer, TLayer};
use crate::layers::{order, LayerEmitter};
use crate::options::Options;

pub struct FixedHLayer {
    fixed: FixedH,
    options: Options,
}

impl FixedHLayer {
    #[must_use]
    pub fn new(fixed: FixedH, options: Options) -> Self {
        Self { fixed, options }
    }

    /// Whether the identifier's /f block carries NON-isotopic content
    /// (formula/h/q/stereo) — the flag that lets the AuxInfo speak the
    /// /F: repetition block (ichiprt1.c's bNonTautNonIsoIdentifierNotEmpty
    /// gating the aux second pass; an FI-only /f block leaves it clear).
    pub fn non_iso_content(&self, structure: &Structure) -> Result<bool> {
        Ok(self.segments(structure)?.non_iso_content)
    }

    /// The computed segment-difference pieces.
    fn segments(&self, structure: &Structure) -> Result<Segments> {
        standardization::ensure_typable_charges(structure)?;
        let fixed = &self.fixed;
        // The mobile formula (the namesake the fixed formula compares
        // against — CompINChILayers' CompareHillFormulasNoH plus the
        // equal-total-H condition, i.e. full string equality).
        let mobile_formula = FormulaLayer.emit(structure)?.unwrap_or_default();
        let formula_differs = fixed.formula != mobile_formula;

        // /f /h — the nonzero deltas spoken at fixed-H canonical
        // numbers (MakeHString over nNum_H_fixed).
        let mut by_rank: Vec<(u32, i32)> = fixed
            .deltas
            .iter()
            .enumerate()
            .filter(|&(_, &delta)| delta != 0)
            .map(|(atom, &delta)| (fixed.ranks[atom], delta))
            .collect();
        by_rank.sort_unstable();
        let h_text = fixed_h_string(&by_rank);

        // /f /q — the fixed-H representation's total charge against
        // the mobile one (CompINChILayers DIFS_q). The equal-charge
        // case is EQL2PRECED, whose segment action is OMIT (bit 1
        // carries neither OMIT_F nor FILL bits): no /q prints at all —
        // the `qm` repetition mark only ever renders in the
        // multi-component list writers. The empty fixed charge over a
        // charged mobile one is IS_EMPTY: the bare slash (breadth-372).
        let mobile_charge = structure.total_charge();
        let fixed_charge = fixed.total_charge;
        let q_segment = match (mobile_charge, fixed_charge) {
            (0, 0) => None,
            (m, f) if m != 0 && f == m => None, // EQL2PRECED: omitted
            (_, 0) => Some(String::new()),      // IS_EMPTY: the bare slash
            (_, f) => Some(format!("{f:+}")),
        };

        // /f stereo — the fixed representation's rows against the
        // main block's; the recorded corpus never re-enters them.
        let stereo_main = self.main_stereo(structure)?;
        let stereo_fixed = self.fixed_stereo(structure)?;
        let stereo_diff = stereo_main != stereo_fixed;

        // The FI block — the fixed representation's isotopic entries
        // (every isotopic H fixed on its carrier) against the mobile
        // /i text.
        let mi_text = IsotopeLayer.emit(structure)?;
        let fi_text = self.fixed_isotopes(structure)?;
        let fi_segment = match (&mi_text, &fi_text) {
            (None, None) => None,
            (Some(_), None) => Some(String::new()), // IS_EMPTY
            (a, b) if a == b => None,               // EQL2PRECED
            (_, Some(text)) => Some(text.clone()),
        };

        // The segment-difference bits that can activate the F layer
        // (MarkUnusedAndEmptyLayers): a NEQ or IS_EMPTY value. A
        // /q alone can never summon the block: the equal-charge case
        // omits (breadth-11's /q+1 keeps no /f) and only the IS_EMPTY
        // bare slash or a differing value carries a bit.
        let non_iso_content =
            formula_differs || h_text.is_some() || q_segment.is_some() || stereo_diff;

        Ok(Segments {
            formula_differs,
            h_text,
            q_segment,
            stereo_main,
            stereo_fixed,
            stereo_diff,
            fi_segment,
            non_iso_content,
        })
    }

    /// The main block's stereo segments (the options carry SNon and
    /// the stereo-type modes).
    fn main_stereo(&self, structure: &Structure) -> Result<Vec<(char, String)>> {
        let mut out = Vec::new();
        for layer in [
            &BLayer::with_options(self.options.clone()) as &dyn LayerEmitter,
            &TLayer::with_options(self.options.clone()),
            &MLayer::with_options(self.options.clone()),
            &SLayer::with_options(self.options.clone()),
        ] {
            if let Some(segment) = layer.emit(structure)? {
                out.push((layer.segment_kind().expect("stereo kind"), segment));
            }
        }
        Ok(out)
    }

    /// The fixed representation's stereo rows: the same perceptions
    /// (bond stereo, then the tetrahedral set) over the fixed-H
    /// numbering.
    fn fixed_stereo(&self, structure: &Structure) -> Result<Vec<(char, String)>> {
        if self.options.no_stereo {
            return Ok(Vec::new());
        }
        let counts = standardization::adjusted_hydrogens(structure);
        let emission = crate::canonical::stereo::Emission {
            keep_undefined: self.options.include_unknown_stereo,
            label_unknown: self.options.distinct_unknown_undefined,
        };
        let mut segments = crate::canonical::bond_stereo::segments(
            structure,
            &counts,
            &self.fixed.ranks,
            emission,
        );
        segments.extend(crate::canonical::stereo::segments(
            structure,
            &counts,
            &self.fixed.ranks,
            &self.fixed.symmetry,
            emission,
        ));
        Ok(segments)
    }

    /// The fixed representation's isotopic text: every isotopic H
    /// fixed on its carrier (no exchangeable aggregate), entries at
    /// the fixed-H numbering.
    fn fixed_isotopes(&self, structure: &Structure) -> Result<Option<String>> {
        let mut entries: Vec<(u32, String)> = Vec::new();
        for (index, atom) in structure.atoms.iter().enumerate() {
            if let Some(entry) = heavy_entry(atom.isotope, &atom.element) {
                entries.push((self.fixed.ranks[index], entry));
            }
        }
        for (carrier, (d, t)) in &structure.isotopic_h {
            entries.push((self.fixed.ranks[*carrier], h_entry(*t, *d)));
        }
        entries.retain(|(_, entry)| !entry.is_empty());
        if entries.is_empty() {
            return Ok(None);
        }
        entries.sort_by_key(|(rank, _)| *rank);
        let mut segment = String::from("i");
        for (position, (rank, entry)) in entries.iter().enumerate() {
            if position > 0 {
                segment.push(',');
            }
            segment.push_str(&rank.to_string());
            segment.push_str(entry);
        }
        Ok(Some(segment))
    }
}

fn heavy_entry(isotope: u16, element: &str) -> Option<String> {
    if isotope == 0 || element == "H" {
        return None;
    }
    let delta = isotope as i32 - standard_mass(element)?;
    Some(format!("{delta:+}"))
}

fn h_entry(t: u32, d: u32) -> String {
    let mut entry = String::new();
    if t > 0 {
        entry.push('T');
        if t > 1 {
            entry.push_str(&t.to_string());
        }
    }
    if d > 0 {
        entry.push('D');
        if d > 1 {
            entry.push_str(&d.to_string());
        }
    }
    entry
}

fn standard_mass(element: &str) -> Option<i32> {
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
        _ => return None,
    })
}

/// The computed segment-difference pieces (CompINChILayers' DIFL_F
/// column for one component).
struct Segments {
    formula_differs: bool,
    h_text: Option<String>,
    q_segment: Option<String>,
    stereo_main: Vec<(char, String)>,
    stereo_fixed: Vec<(char, String)>,
    stereo_diff: bool,
    fi_segment: Option<String>,
    non_iso_content: bool,
}

impl LayerEmitter for FixedHLayer {
    fn name(&self) -> &'static str {
        "FixedHLayer"
    }

    fn order(&self) -> u8 {
        order::O
    }

    fn segment_kind(&self) -> Option<char> {
        Some('f')
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        let segments = self.segments(structure)?;

        // MarkUnusedAndEmptyLayers: the whole F layer omits when no
        // segment carries a NEQ/IS_EMPTY bit and the FI layer is
        // empty. The q_segment values already carry only bits (the
        // equal-charge case resolved to None above).
        let active = segments.formula_differs
            || segments.h_text.is_some()
            || segments.q_segment.is_some()
            || segments.stereo_diff
            || segments.fi_segment.is_some();
        if !active {
            return Ok(None);
        }

        // The formula segment (IL_FMLF | IL_FIXH): the layer's first
        // segment — an equal formula promotes to the empty marker and
        // renders the bare `/f`.
        let mut out = String::from("f");
        if segments.formula_differs {
            out.push_str(&self.fixed.formula);
        }
        // /f /h (IL_HFIX).
        if let Some(text) = &segments.h_text {
            out.push('/');
            out.push('h');
            out.push_str(text);
        }
        // /f /q (IL_CHRG): the differing value, or the bare slash when
        // the fixed representation lost the charge (IS_EMPTY).
        if let Some(text) = &segments.q_segment {
            out.push('/');
            out.push('q');
            out.push_str(text);
        }
        // /f /b /t /m /s — re-entered only when the fixed rows differ
        // from the main ones (the recorded corpus keeps them equal).
        if segments.stereo_diff {
            for (kind, text) in &segments.stereo_fixed {
                let main = segments
                    .stereo_main
                    .iter()
                    .find(|(main_kind, _)| main_kind == kind);
                match main {
                    Some((_, main_text)) if main_text == text => {}
                    Some(_) => {
                        out.push('/');
                        out.push(*kind);
                        out.push_str(text);
                    }
                    None => {
                        out.push('/');
                        out.push(*kind);
                    }
                }
            }
            for (kind, _text) in &segments.stereo_main {
                if !segments.stereo_fixed.iter().any(|(k, _)| k == kind) {
                    // IS_EMPTY: the fixed representation lost the row.
                    out.push('/');
                    out.push(*kind);
                }
            }
        }
        // The FI block (IL_ATMS under the fixed-H tag).
        if let Some(text) = &segments.fi_segment {
            out.push('/');
            out.push_str(text);
        }
        Ok(Some(out))
    }
}

/// MakeHString over the fixed deltas: counts ascending (negative
/// counts first, lowercase `h`), ranks ascending within a count,
/// consecutive runs collapsed (`32,34,36H,22H2`, `1h`).
fn fixed_h_string(by_rank: &[(u32, i32)]) -> Option<String> {
    if by_rank.is_empty() {
        return None;
    }
    let counts: Vec<i32> = {
        let mut counts: Vec<i32> = by_rank.iter().map(|&(_, count)| count).collect();
        counts.sort_unstable();
        counts.dedup();
        counts
    };
    let mut segment = String::new();
    for count in counts {
        if count == 0 {
            continue; // MakeHString skips the zero-H group.
        }
        let mut ranks: Vec<u32> = by_rank
            .iter()
            .filter(|&(_, c)| *c == count)
            .map(|&(rank, _)| rank)
            .collect();
        ranks.sort_unstable();
        if !segment.is_empty() {
            segment.push(',');
        }
        push_runs(&ranks, &mut segment);
        let letter = if count > 0 { 'H' } else { 'h' };
        let magnitude = count.unsigned_abs();
        segment.push(letter);
        if magnitude > 1 {
            segment.push_str(&magnitude.to_string());
        }
    }
    if segment.is_empty() {
        return None;
    }
    Some(segment)
}

fn push_runs(ranks: &[u32], out: &mut String) {
    let mut index = 0;
    while index < ranks.len() {
        let start = ranks[index];
        let mut end = start;
        while index + 1 < ranks.len() && ranks[index + 1] == end + 1 {
            index += 1;
            end = ranks[index];
        }
        if start == end {
            out.push_str(&start.to_string());
        } else {
            out.push_str(&format!("{start}-{end}"));
        }
        index += 1;
        if index < ranks.len() {
            out.push(',');
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delta_strings_group_negative_first() {
        // hydroxide: the -1 delta prints lowercase (fHO/h1h).
        let by_rank = vec![(1, -1)];
        assert_eq!(fixed_h_string(&by_rank).as_deref(), Some("1h"));
        // breadth-53's shape: three single deltas then a +2.
        let by_rank = vec![(22, 2), (32, 1), (34, 1), (36, 1)];
        assert_eq!(fixed_h_string(&by_rank).as_deref(), Some("32,34,36H,22H2"));
        // tartaric's shape: two holders in one run.
        let by_rank = vec![(7, 1), (9, 1)];
        assert_eq!(fixed_h_string(&by_rank).as_deref(), Some("7,9H"));
        // aspirin: a single holder.
        let by_rank = vec![(11, 1)];
        assert_eq!(fixed_h_string(&by_rank).as_deref(), Some("11H"));
    }

    #[test]
    fn empty_deltas_render_no_h_segment() {
        assert!(fixed_h_string(&[]).is_none());
        assert!(fixed_h_string(&[(3, 0)]).is_none());
    }
}
