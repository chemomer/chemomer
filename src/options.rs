//! The generation-option surface: the reference's InChI 1.07.5 option
//! flags (ichiparm.c), carried by the [`crate::Engine`] and consulted
//! where each behavior lives — the stereo layers at their emit step,
//! the hydrogen perception at the engine's normalization entry. The
//! default value is the standard identifier: the default pipeline is
//! bit-identical to the pre-options engine because every flag reads
//! false and the layers produce today's bytes under
//! [`Options::default`].
//!
//! Flag semantics witnessed against the clean 1.07.5 oracle binary:
//!
//! | flag      | version     | stereo layers          | notes |
//! |-----------|-------------|------------------------|-------|
//! | `SNon`    | `InChI=1S/` | `/b /t /m /s` withheld | standard-legal; perception still runs (the canonical numbering keeps the stereo-informed verdict) |
//! | `SRel`    | `InChI=1/`  | `/t` + `/s2`, no `/m`  | rows are the lexicographically smaller of the absolute/inverted passes |
//! | `SRac`    | `InChI=1/`  | `/t` + `/s3`, no `/m`  | same row rule |
//! | `SUCF`    | `InChI=1/`  | per molfile chiral flag | `/s1` when set, `/s2` when clear |
//! | `SUU`     | `InChI=1/`  | undefined centers/bonds print (`?`) | default drops the whole group when nothing is well-defined |
//! | `SLUUD`   | `InChI=1/`  | unknown prints `u` | distinguishes wiggly (unknown) from undefined |
//! | `DoNotAddH` | `InChI=1S/` | — | only drawn hydrogens exist; H-count codes are ignored too |

/// The option flags one generation runs under. Builder methods mirror
/// the reference's parse-order toggles exactly, so a chain of
/// builders reads like the option line it encodes.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Options {
    /// `SNon` — the stereo layers never emit; the aux
    /// inverted-stereo block is cleared. Standard-legal: the version
    /// stays `InChI=1S/` (ichiparm.c:248–253).
    pub no_stereo: bool,
    /// `SRel` — relative stereo (ichiparm.c:272–286).
    pub relative_stereo: bool,
    /// `SRac` — racemic stereo (ichiparm.c:287–301).
    pub racemic_stereo: bool,
    /// `SUCF` — the molfile chiral flag decides between absolute
    /// (`/s1`) and relative (`/s2`) (ichiparm.c:302–316,
    /// runichi2.c:538–552). Needs the chiral flag the molfile reader
    /// carries in the header counts line.
    pub chiral_flag_stereo: bool,
    /// `SUU` — keep the omitted unknown/undefined stereo: candidate
    /// centers and bonds without well-defined parity print instead of
    /// the whole group dropping (ichiparm.c:336–342,
    /// ichimak2.c `UnmarkAllUndefinedUnknownStereo`).
    pub include_unknown_stereo: bool,
    /// `SLUUD` — unknown (wiggly) and undefined stereo carry distinct
    /// labels: `u` where the standard prints `?` (ichiparm.c:343–349).
    pub distinct_unknown_undefined: bool,
    /// `DoNotAddH` — no implicit-hydrogen inference: only drawn
    /// hydrogens exist, and the molfile H-count codes are ignored
    /// (ichiparm.c:259–263, mol2atom.c `calculate_valences`).
    pub do_not_add_h: bool,
}

impl Options {
    /// `SNon` — also turns `SUU`/`SLUUD` off: no stereo mode leaves no
    /// unknown/undefined stereo to speak of (the reference's
    /// postprocess clears the two bits whenever `REQ_MODE_STEREO` is
    /// absent, ichiparm.c:2038–2043).
    #[must_use]
    pub fn with_no_stereo(mut self) -> Self {
        self.no_stereo = true;
        self.include_unknown_stereo = false;
        self.distinct_unknown_undefined = false;
        self
    }

    /// `SRel` — XORs the sibling stereo-type modes off, exactly like
    /// the reference's parse toggles: each of `SRel`/`SRac`/`SUCF`
    /// clears the other two, so the last builder call wins.
    #[must_use]
    pub fn with_relative_stereo(mut self) -> Self {
        self.racemic_stereo = false;
        self.chiral_flag_stereo = false;
        self.relative_stereo = true;
        self
    }

    /// `SRac` — see [`Options::with_relative_stereo`] for the XOR rule.
    #[must_use]
    pub fn with_racemic_stereo(mut self) -> Self {
        self.relative_stereo = false;
        self.chiral_flag_stereo = false;
        self.racemic_stereo = true;
        self
    }

    /// `SUCF` — see [`Options::with_relative_stereo`] for the XOR rule.
    #[must_use]
    pub fn with_chiral_flag_stereo(mut self) -> Self {
        self.relative_stereo = false;
        self.racemic_stereo = false;
        self.chiral_flag_stereo = true;
        self
    }

    /// `SUU` — keep the omitted unknown/undefined stereo.
    #[must_use]
    pub fn with_include_unknown_stereo(mut self) -> Self {
        self.include_unknown_stereo = true;
        self
    }

    /// `SLUUD` — distinct labels for unknown vs undefined stereo.
    #[must_use]
    pub fn with_distinct_unknown_undefined(mut self) -> Self {
        self.distinct_unknown_undefined = true;
        self
    }

    /// `DoNotAddH` — only drawn hydrogens exist.
    #[must_use]
    pub fn with_do_not_add_h(mut self) -> Self {
        self.do_not_add_h = true;
        self
    }

    /// Whether the assembled identifier is the standard one
    /// (`InChI=1S/`): `SRel`/`SRac`/`SUCF`/`SUU`/`SLUUD` each turn the
    /// std flag off; `SNon` and `DoNotAddH` leave it standing (the
    /// reference's `bStdFormat`: the perception options "DO NOT TURN
    /// OFF Std flag", the format options do).
    #[must_use]
    pub fn is_standard(&self) -> bool {
        !(self.relative_stereo
            || self.racemic_stereo
            || self.chiral_flag_stereo
            || self.include_unknown_stereo
            || self.distinct_unknown_undefined)
    }

    /// Whether the identifier speaks a stereo TYPE segment (`/s`) with
    /// the non-absolute value: `Some('2')` under `SRel`, `Some('3')`
    /// under `SRac`, `None` otherwise (the absolute `/s1` and the
    /// withheld `/s` of `SNon` are not decided here — the `/s` layer
    /// and the per-structure verdict decide those).
    #[must_use]
    pub fn requested_stereo_type(&self) -> Option<char> {
        if self.relative_stereo {
            Some('2')
        } else if self.racemic_stereo {
            Some('3')
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_assemble_the_standard_identifier() {
        let options = Options::default();
        assert!(options.is_standard());
        assert_eq!(options.requested_stereo_type(), None);
    }

    #[test]
    fn perception_flags_keep_the_std_flag() {
        // SNon and DoNotAddH are standard-legal.
        assert!(Options::default().with_no_stereo().is_standard());
        assert!(Options::default().with_do_not_add_h().is_standard());
    }

    #[test]
    fn format_flags_turn_the_std_flag_off() {
        for options in [
            Options::default().with_relative_stereo(),
            Options::default().with_racemic_stereo(),
            Options::default().with_chiral_flag_stereo(),
            Options::default().with_include_unknown_stereo(),
            Options::default().with_distinct_unknown_undefined(),
        ] {
            assert!(!options.is_standard(), "{options:?}");
        }
    }

    #[test]
    fn the_stereo_type_modes_xor_each_other_off() {
        // The last flag on the option line wins, mirroring the
        // reference's parse-order toggles.
        let options = Options::default()
            .with_relative_stereo()
            .with_racemic_stereo();
        assert!(options.racemic_stereo);
        assert!(!options.relative_stereo);
        assert!(!options.chiral_flag_stereo);
        let options = options.with_chiral_flag_stereo();
        assert!(options.chiral_flag_stereo);
        assert!(!options.racemic_stereo);
    }

    #[test]
    fn snon_clears_the_unknown_undefined_bits() {
        let options = Options::default()
            .with_include_unknown_stereo()
            .with_distinct_unknown_undefined()
            .with_no_stereo();
        assert!(options.no_stereo);
        assert!(!options.include_unknown_stereo);
        assert!(!options.distinct_unknown_undefined);
        // SNon is standard-legal even after SUU/SLUUD turned the flag
        // off in the reference's sense: the version reads InChI=1/ only
        // when a FORMAT option survives, and here none does.
        assert!(options.is_standard());
    }

    #[test]
    fn the_requested_stereo_type_reads_two_and_three() {
        assert_eq!(
            Options::default()
                .with_relative_stereo()
                .requested_stereo_type(),
            Some('2')
        );
        assert_eq!(
            Options::default()
                .with_racemic_stereo()
                .requested_stereo_type(),
            Some('3')
        );
    }
}
