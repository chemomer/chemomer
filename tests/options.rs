//! The option spot checks: oracle recordings probed by hand against
//! the clean 1.07.5 binary (2026-11-01, /tmp probes) for the modes the
//! shared corpus does not carry a tranche for yet. Byte equality per
//! assembled identifier; the fixtures' molfiles stand in as constants
//! so the checks run without the corpus checkout. A recorded tranche
//! promotes these to corpus gates.

use chemomer::Options;

const L_ALANINE: &str = "\n  AsciiChem\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n   12.0000   22.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   22.0000    0.0000 C   0  0  2  2  0  0  0  0  0  0  0  0  0\n  132.0000   67.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  6  0  0  0  0  0  0  0\n  2  3  1  0  0  0  0  0  0  0  0\n  4  5  2  0  0  0  0  0  0  0  0\n  2  4  1  0  0  0  0  0  0  0  0\n  4  6  1  0  0  0  0  0  0  0  0\nM  END\n";

const S_2_BUTANOL: &str = "\n  AsciiChem\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n   12.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   17.0000    0.0000 C   0  0  1  2  0  0  0  0  0  0  0  0  0\n  192.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0  0  0  0\n  3  2  1  1  0  0  0  0  0  0  0\n  3  4  1  0  0  0  0  0  0  0  0\n  3  5  1  0  0  0  0  0  0  0  0\nM  END\n";

const MESO_TARTARIC_ACID: &str = "\n  AsciiChem\n\n 10  9  0  0  0  0  0  0  0  0999 V2000\n   12.0000   22.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   22.0000    0.0000 C   0  0  1  2  0  0  0  0  0  0  0  0  0\n  132.0000  117.0000    0.0000 C   0  0  1  2  0  0  0  0  0  0  0  0  0\n  192.0000  162.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000  112.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  252.0000  107.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  252.0000  157.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1  0  0  0  0  0  0  0\n  3  2  1  1  0  0  0  0  0  0  0\n  3  4  1  0  0  0  0  0  0  0  0\n  5  6  2  0  0  0  0  0  0  0  0\n  8  9  2  0  0  0  0  0  0  0  0\n  3  5  1  0  0  0  0  0  0  0  0\n  5  7  1  0  0  0  0  0  0  0  0\n  2  8  1  0  0  0  0  0  0  0  0\n  8 10  1  0  0  0  0  0  0  0  0\nM  END\n";

const RR_2_3_DICHLOROBUTANE: &str = "\n  AsciiChem\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n   12.0000   22.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   22.0000    0.0000 C   0  0  1  2  0  0  0  0  0  0  0  0  0\n  132.0000   67.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   17.0000    0.0000 C   0  0  1  2  0  0  0  0  0  0  0  0  0\n  192.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1  0  0  0  0  0  0  0\n  2  3  1  0  0  0  0  0  0  0  0\n  4  5  1  1  0  0  0  0  0  0  0\n  2  4  1  0  0  0  0  0  0  0  0\n  4  6  1  0  0  0  0  0  0  0  0\nM  END\n";

const E_2_BUTENE: &str = "e-2-butene\n  AsciiChem\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   12.0000   85.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   60.0000   50.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  120.0000   50.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  168.0000   15.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0\n  2  3  2  0  0  0  0  0\n  3  4  1  0  0  0  0  0\nM  END\n";

fn identifier(molfile: &str, options: Options) -> String {
    chemomer::Engine::with_options(molfile, options)
        .expect("parses")
        .inchi_string()
        .expect("assembles")
}

#[test]
fn relative_and_racemic_reshape_the_stereo_tail() {
    // Oracle: InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1
    let body = "C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)";
    assert_eq!(
        identifier(L_ALANINE, Options::default()),
        format!("InChI=1S/{body}/t2-/m1/s1")
    );
    assert_eq!(
        identifier(L_ALANINE, Options::default().with_relative_stereo()),
        format!("InChI=1/{body}/t2-/s2")
    );
    assert_eq!(
        identifier(L_ALANINE, Options::default().with_racemic_stereo()),
        format!("InChI=1/{body}/t2-/s3")
    );
    // Oracle: InChI=1S/C4H10O/c1-3-4(2)5/h4-5H,3H2,1-2H3/t4-/m0/s1
    let body = "C4H10O/c1-3-4(2)5/h4-5H,3H2,1-2H3";
    assert_eq!(
        identifier(S_2_BUTANOL, Options::default().with_relative_stereo()),
        format!("InChI=1/{body}/t4-/s2")
    );
    assert_eq!(
        identifier(S_2_BUTANOL, Options::default().with_racemic_stereo()),
        format!("InChI=1/{body}/t4-/s3")
    );
}

#[test]
fn the_tied_inversion_verdict_emits_no_stereo_type() {
    // The meso (rows-tie) case: /t stands bare, no /m, no /s — the
    // only change the modes make is the version prefix. Oracle:
    // InChI=1S/C4H6O6/.../t1-,2+ under -SRel stays /t1-,2+ with
    // InChI=1/.
    let body = "C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)";
    assert_eq!(
        identifier(MESO_TARTARIC_ACID, Options::default()),
        format!("InChI=1S/{body}/t1-,2+")
    );
    assert_eq!(
        identifier(
            MESO_TARTARIC_ACID,
            Options::default().with_relative_stereo()
        ),
        format!("InChI=1/{body}/t1-,2+")
    );
    assert_eq!(
        identifier(MESO_TARTARIC_ACID, Options::default().with_racemic_stereo()),
        format!("InChI=1/{body}/t1-,2+")
    );
    let body = "C4H8Cl2/c1-3(5)4(2)6/h3-4H,1-2H3";
    assert_eq!(
        identifier(
            RR_2_3_DICHLOROBUTANE,
            Options::default().with_relative_stereo()
        ),
        format!("InChI=1/{body}/t3-,4+")
    );
    assert_eq!(
        identifier(
            RR_2_3_DICHLOROBUTANE,
            Options::default().with_racemic_stereo()
        ),
        format!("InChI=1/{body}/t3-,4+")
    );
}

#[test]
fn bond_stereo_survives_the_relative_racemic_modes() {
    // Oracle: -SRel on the E fixture keeps /b4-3+ and records no /s
    // (bond-only stereo carries no /t set — SP3_NONE) — the /b layer
    // is outside the SRel/SRac reshape, and SNon withholds it.
    let body = "C4H8/c1-3-4-2/h3-4H,1-2H3";
    assert_eq!(
        identifier(E_2_BUTENE, Options::default()),
        format!("InChI=1S/{body}/b4-3+")
    );
    assert_eq!(
        identifier(E_2_BUTENE, Options::default().with_relative_stereo()),
        format!("InChI=1/{body}/b4-3+")
    );
    assert_eq!(
        identifier(E_2_BUTENE, Options::default().with_racemic_stereo()),
        format!("InChI=1/{body}/b4-3+")
    );
    assert_eq!(
        identifier(E_2_BUTENE, Options::default().with_no_stereo()),
        format!("InChI=1S/{body}")
    );
}

const WIGGLY_ALANINE: &str = "\n  AsciiChem\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n   12.0000   22.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   22.0000    0.0000 C   0  0  2  2  0  0  0  0  0  0  0  0  0\n  132.0000   67.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  4  0  0  0  0  0  0  0\n  2  3  1  0  0  0  0  0  0  0  0\n  4  5  2  0  0  0  0  0  0  0  0\n  2  4  1  0  0  0  0  0  0  0  0\n  4  6  1  0  0  0  0  0  0  0  0\nM  END\n";

const TRULY_FLAT_BUTANOL: &str = "\n  AsciiChem\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n   12.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   17.0000    0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0  0\n  192.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0  0  0  0\n  3  2  1  0  0  0  0  0  0  0  0\n  3  4  1  0  0  0  0  0  0  0  0\n  3  5  1  0  0  0  0  0  0  0  0\nM  END\n";

const CROSSED_2_BUTENE: &str = "crossed-2-butene\n  probe\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   12.0000   22.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   22.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   32.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  3  0  0  0  0\n  1  3  1  0  0  0  0  0\n  2  4  1  0  0  0  0  0\nM  END\n";

const MIXED_BOND: &str = "mixed-bond\n  probe\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n   12.0000   22.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   22.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   62.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  252.0000  102.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  312.0000  102.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0  0\n  2  3  2  0  0  0  0  0\n  3  4  1  0  0  0  0  0\n  4  5  2  3  0  0  0  0\n  5  6  1  0  0  0  0  0\nM  END\n";

const MIXED_DCB: &str = "mixed-dcb\n  probe\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n   12.0000   22.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   22.0000    0.0000 C   0  0  1  2  0  0  0  0  0  0  0  0  0\n  132.0000   67.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   17.0000    0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0  0\n  192.0000   12.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1  0  0  0  0\n  2  3  1  0  0  0  0  0\n  4  5  1  0  0  0  0  0\n  2  4  1  0  0  0  0  0\n  4  6  1  0  0  0  0  0\nM  END\n";

#[test]
fn suu_keeps_the_omitted_unknown_undefined_stereo() {
    // Oracle (-SUU on the wiggly-alanine probe): the wiggly center is
    // UNKNOWN parity — the default drops the group, SUU prints t2?,
    // and SUU+SLUUD relabels t2u. SLUUD alone changes nothing (the
    // drop decision belongs to SUU).
    let body = "C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)";
    assert_eq!(
        identifier(WIGGLY_ALANINE, Options::default()),
        format!("InChI=1S/{body}")
    );
    assert_eq!(
        identifier(
            WIGGLY_ALANINE,
            Options::default().with_include_unknown_stereo()
        ),
        format!("InChI=1/{body}/t2?")
    );
    assert_eq!(
        identifier(
            WIGGLY_ALANINE,
            Options::default()
                .with_include_unknown_stereo()
                .with_distinct_unknown_undefined()
        ),
        format!("InChI=1/{body}/t2u")
    );
    assert_eq!(
        identifier(
            WIGGLY_ALANINE,
            Options::default().with_distinct_unknown_undefined()
        ),
        format!("InChI=1/{body}")
    );
    // Oracle (-SUU on the flat probe): an unwedged candidate center
    // prints as UNDEFINED.
    let body = "C4H10O/c1-3-4(2)5/h4-5H,3H2,1-2H3";
    assert_eq!(
        identifier(TRULY_FLAT_BUTANOL, Options::default()),
        format!("InChI=1S/{body}")
    );
    assert_eq!(
        identifier(
            TRULY_FLAT_BUTANOL,
            Options::default().with_include_unknown_stereo()
        ),
        format!("InChI=1/{body}/t4?")
    );
    // Oracle (-SUU on the crossed probe): the crossed (code 3) bond is
    // UNKNOWN parity — b4-3? kept, b4-3u under SUU+SLUUD.
    let body = "C4H8/c1-3-4-2/h3-4H,1-2H3";
    assert_eq!(
        identifier(CROSSED_2_BUTENE, Options::default()),
        format!("InChI=1S/{body}")
    );
    assert_eq!(
        identifier(
            CROSSED_2_BUTENE,
            Options::default().with_include_unknown_stereo()
        ),
        format!("InChI=1/{body}/b4-3?")
    );
    assert_eq!(
        identifier(
            CROSSED_2_BUTENE,
            Options::default()
                .with_include_unknown_stereo()
                .with_distinct_unknown_undefined()
        ),
        format!("InChI=1/{body}/b4-3u")
    );
}

#[test]
fn sluud_relabels_printed_unknown_stereo_without_keeping_it() {
    // Oracle (the mixed E/crossed probe): the crossed bond rides the
    // group the defined bond keeps — `b5-3+,6-4?` in every default
    // run, `6-4u` once SLUUD distinguishes the unknown parity.
    let body = "C6H10/c1-3-5-6-4-2/h3-6H,1-2H3";
    assert_eq!(
        identifier(MIXED_BOND, Options::default()),
        format!("InChI=1S/{body}/b5-3+,6-4?")
    );
    assert_eq!(
        identifier(
            MIXED_BOND,
            Options::default().with_distinct_unknown_undefined()
        ),
        format!("InChI=1/{body}/b5-3+,6-4u")
    );
    assert_eq!(
        identifier(MIXED_BOND, Options::default().with_include_unknown_stereo()),
        format!("InChI=1/{body}/b5-3+,6-4?")
    );
    // Oracle (the one-wedge dichlorobutane probe): the UNDEFINED
    // candidate rides as `?` in every mode — SLUUD never relabels it.
    let body = "C4H8Cl2/c1-3(5)4(2)6/h3-4H,1-2H3";
    for options in [
        Options::default(),
        Options::default().with_include_unknown_stereo(),
        Options::default().with_distinct_unknown_undefined(),
    ] {
        // The group reads identically in every mode (`?` never
        // relabels); only the version prefix moves — SUU and SLUUD
        // each turn the std flag off.
        let prefix = if options.is_standard() {
            "InChI=1S"
        } else {
            "InChI=1"
        };
        let expected = format!("{prefix}/{body}/t3-,4?/m1/s1");
        assert_eq!(
            identifier(MIXED_DCB, options.clone()),
            expected,
            "{options:?}"
        );
    }
}
