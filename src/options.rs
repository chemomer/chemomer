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
//!
//! The key-output flags (`XHash1`/`XHash2`/`MergeHash`) sit beside the
//! table: they never touch the assembled identifier — they extend the
//! [`crate::key`] composition the engine prints (see
//! [`crate::Engine::inchikey_output`]).
//!
//! The experimental tautomerism flags (`KET`, `15T`, `PT_*`,
//! ichiparm.c 1988–2005) ride the same surface: each enables one
//! extra perception pass (the keto-enol C/O anchors, the 1,5
//! shift-walk, the six PT endpoint families), and each turns the
//! std flag off — the identifier reads `InChI=1/`.

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
    /// `XHash1` — extend the key's FIRST block hash to 256 bits: the
    /// major digest's un-consumed bits as 48 lowercase hex characters
    /// (ichiparm.c:236–240, ikey_base26.c `get_xtra_hash_major_hex`).
    /// Key-output only: it never touches the identifier. The reference
    /// warns the flag is inert without `Key` ("Hash extension(s) not
    /// generated: InChIKey not requested", ichiparm.c:2061–2063); the
    /// engine's key output reads these flags the same way.
    pub xhash1: bool,
    /// `XHash2` — extend the key's SECOND block hash to 256 bits: the
    /// minor digest's un-consumed bits as 56 lowercase hex characters
    /// (ichiparm.c:241–245, ikey_base26.c `get_xtra_hash_minor_hex`).
    pub xhash2: bool,
    /// `MergeHash` — combine the key with the requested extensions
    /// onto one `InChIHash=` line, dash-joined, instead of the separate
    /// `InChIKey=`/`XHash1=`/`XHash2=` lines (ichiparm.c:455–458,
    /// ichimain.c:1376–1412). An absent extension leaves a bare dash.
    pub merge_hash: bool,
    /// `FixedH` — the /f layer of the full (non-standard) InChI: the
    /// pre-normalization structure speaks its own fixed-H numbering,
    /// H deltas and charges after the standard layers (ichiparm.c
    /// "FixedH — Include Fixed H layer"; bStdFormat off, ichiparm.c
    /// 248-253).
    pub fixed_h: bool,
    /// `RecMet` — append the reconnected identifier (`/r` + the
    /// identifier over the metal-BONDED structure) after the main
    /// (disconnected) one (ichiparm.c 972-990, mode.h:586
    /// RECONNECT_METALS). Non-standard: the version prefix reads
    /// `InChI=1/`. The /r part exists ONLY when the standard pass
    /// actually disconnected a coordination metal (the reference's
    /// TG_FLAG_DISCONNECT_COORD_DONE + RECONNECT_COORD pair,
    /// runichi.c:980) — a salt-only disconnection never reconnects.
    pub rec_met: bool,
    /// `KET` — keto-enol tautomerism: the extra 1,3 registration pass
    /// over C/O endpoint pairs (TG_FLAG_KETO_ENOL_TAUT, ichiparm.c
    /// 367–371). Non-standard: the version prefix reads `InChI=1/`.
    pub keto_enol: bool,
    /// `15T` — 1,5 tautomerism across alternating paths (one endpoint
    /// off-ring): the shift-walk pass (TG_FLAG_1_5_TAUT, ichiparm.c
    /// 373–378). Non-standard.
    pub t15: bool,
    /// `PT_22_00` — the imine/enamine endpoint pass
    /// `[#1:1][CX4:2][NX2:3]=[CX3:4]` (ichiparm.c 383–388).
    pub pt_22_00: bool,
    /// `PT_16_00` — the N-oxide/oxime endpoint pass
    /// `[#1:1][O;!R:2][N:3]=[CX3:4]` (ichiparm.c 390–395).
    pub pt_16_00: bool,
    /// `PT_06_00` — the general hetero-ene endpoint pass with Q = N,C,P
    /// (ichiparm.c 397–402).
    pub pt_06_00: bool,
    /// `PT_39_00` — the nitronic-acid endpoint pass
    /// `[CX3,NX2:1]=[NX3+:2]([O-:3])[CX4:4]` (ichiparm.c 404–409).
    pub pt_39_00: bool,
    /// `PT_13_00` — the thio-cumulene/alkyne endpoint pass
    /// `[O,S,Se,Te;X1:1]=[C:2]=[C:3][#1:4]` (ichiparm.c 411–416).
    pub pt_13_00: bool,
    /// `PT_18_00` — the cyanic-acid endpoint pass
    /// `[#1:1][O:2][C:3]#[N:4]` (ichiparm.c 418–423).
    pub pt_18_00: bool,
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

    /// `XHash1` — the key's first block hash extends to 256 bits.
    #[must_use]
    pub fn with_xhash1(mut self) -> Self {
        self.xhash1 = true;
        self
    }

    /// `XHash2` — the key's second block hash extends to 256 bits.
    #[must_use]
    pub fn with_xhash2(mut self) -> Self {
        self.xhash2 = true;
        self
    }

    /// `MergeHash` — key and requested extensions join onto one
    /// `InChIHash=` line.
    #[must_use]
    pub fn with_merge_hash(mut self) -> Self {
        self.merge_hash = true;
        self
    }

    /// `FixedH` — the full-InChI /f layer. Non-standard: the version
    /// prefix reads `InChI=1/`.
    #[must_use]
    pub fn with_fixed_h(mut self) -> Self {
        self.fixed_h = true;
        self
    }

    /// `RecMet` — the reconnected identifier rides after the main
    /// one. Non-standard: the version prefix reads `InChI=1/`.
    #[must_use]
    pub fn with_rec_met(mut self) -> Self {
        self.rec_met = true;
        self
    }

    /// `KET` — keto-enol tautomerism. Non-standard.
    #[must_use]
    pub fn with_ket(mut self) -> Self {
        self.keto_enol = true;
        self
    }

    /// `15T` — 1,5 tautomerism across alternating paths. Non-standard.
    #[must_use]
    pub fn with_15t(mut self) -> Self {
        self.t15 = true;
        self
    }

    /// `PT_22_00` — the imine/enamine pass. Non-standard.
    #[must_use]
    pub fn with_pt_22_00(mut self) -> Self {
        self.pt_22_00 = true;
        self
    }

    /// `PT_16_00` — the N-oxide/oxime pass. Non-standard.
    #[must_use]
    pub fn with_pt_16_00(mut self) -> Self {
        self.pt_16_00 = true;
        self
    }

    /// `PT_06_00` — the general hetero-ene pass. Non-standard.
    #[must_use]
    pub fn with_pt_06_00(mut self) -> Self {
        self.pt_06_00 = true;
        self
    }

    /// `PT_39_00` — the nitronic-acid pass. Non-standard.
    #[must_use]
    pub fn with_pt_39_00(mut self) -> Self {
        self.pt_39_00 = true;
        self
    }

    /// `PT_13_00` — the thio-cumulene/alkyne pass. Non-standard.
    #[must_use]
    pub fn with_pt_13_00(mut self) -> Self {
        self.pt_13_00 = true;
        self
    }

    /// `PT_18_00` — the cyanic-acid pass. Non-standard.
    #[must_use]
    pub fn with_pt_18_00(mut self) -> Self {
        self.pt_18_00 = true;
        self
    }

    /// Whether the assembled identifier is the standard one
    /// (`InChI=1S/`): `SRel`/`SRac`/`SUCF`/`SUU`/`SLUUD`, `FixedH`,
    /// `RecMet` and the experimental tautomerism flags each turn the
    /// std flag off; `SNon` and `DoNotAddH` leave it standing (the
    /// reference's `bStdFormat`: the perception options "DO NOT TURN
    /// OFF Std flag", the format options do).
    #[must_use]
    pub fn is_standard(&self) -> bool {
        !(self.relative_stereo
            || self.racemic_stereo
            || self.chiral_flag_stereo
            || self.include_unknown_stereo
            || self.distinct_unknown_undefined
            || self.fixed_h
            || self.rec_met
            || self.keto_enol
            || self.t15
            || self.pt_22_00
            || self.pt_16_00
            || self.pt_06_00
            || self.pt_39_00
            || self.pt_13_00
            || self.pt_18_00)
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
            Options::default().with_fixed_h(),
            Options::default().with_rec_met(),
            Options::default().with_ket(),
            Options::default().with_15t(),
            Options::default().with_pt_22_00(),
            Options::default().with_pt_16_00(),
            Options::default().with_pt_06_00(),
            Options::default().with_pt_39_00(),
            Options::default().with_pt_13_00(),
            Options::default().with_pt_18_00(),
        ] {
            assert!(!options.is_standard(), "{options:?}");
        }
    }

    #[test]
    fn the_tautomer_flags_stay_independent() {
        // Each builder sets exactly its own bit; combinations carry
        // both (-KET -15T is a legal option line).
        let options = Options::default().with_ket().with_15t();
        assert!(options.keto_enol && options.t15);
        assert!(!options.pt_22_00);
        let options = options.with_pt_18_00();
        assert!(options.keto_enol && options.t15 && options.pt_18_00);
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
