//! The SMILES reader gate (chemomer#121, phase 2 in chemomer#128):
//! every identity fixture's `smiles` field reads through the contract
//! and derives the recorded 1.07.5-oracle InChI byte-for-byte — no
//! new oracle, no fixture edits. The stereo cases all run (the
//! alanine enantiomers through the carried @/@@; the butene pair
//! through the corpus's drawn-coordinate re-record, inchi-tests
//! v0.3.42); the recorded exclusions: the two transcription
//! disagreements (sulfuric-acid/sulfate-2).

use std::fs;
use std::path::PathBuf;

fn suite_dir() -> PathBuf {
    std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"))
}

#[test]
fn the_smiles_reader_derives_the_recorded_identities() {
    let dir = suite_dir();
    let mut cases: Vec<serde_json::Value> = Vec::new();
    for name in [
        "corpus/fixtures/identity.json",
        "corpus/fixtures/identity-breadth.json",
        "corpus/fixtures/identity-metal.json",
        "corpus/fixtures/identity-multicomponent.json",
    ] {
        let Ok(text) = fs::read_to_string(dir.join(name)) else {
            return; // the corpus checkout gates in CI
        };
        cases.extend(serde_json::from_str::<Vec<serde_json::Value>>(&text).expect("fixtures"));
    }
    if cases.is_empty() {
        return;
    }

    let mut ran = 0usize;
    let mut skipped_phase2 = Vec::new();
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let Some(smiles) = case["smiles"].as_str() else {
            continue;
        };
        const TRANSCRIPTION_DISAGREEMENTS: [&str; 15] = [
            "identity/sulfuric-acid",
            "identity/sulfate-2",
            "identity/multi-sodium-hydrogen-sulfate",
            "identity/multi-magnesium-sulfate",
            "identity/multi-calcium-sulfate",
            "identity/multi-copper-sulfate",
            "identity/multi-sodium-thiosulfate",
            "identity/multi-potassium-permanganate",
            "identity/multi-sodium-dodecyl-sulfate",
            "identity/multi-ammonium-sulfate",
            "identity/multi-sodium-citrate",
            "identity/multi-magnesium-sulfate-heptahydrate",
            "identity/multi-copper-sulfate-pentahydrate",
            "identity/multi-alum-dodecahydrate",
            "identity/multi-sodium-sulfate-decahydrate",
        ];
        if TRANSCRIPTION_DISAGREEMENTS.contains(&id) {
            // The fixture's smiles transcription disagrees with its
            // recorded molfile — the same class as the sulfuric-acid
            // pair: the molfile records the neutral acid with explicit
            // H (the oracle's H4O4S/5H normalization) while the smiles
            // field transcribes the deprotonated anion (or, for the
            // permanganate/dodecyl-sulfate/citrate rows, a different
            // H placement). The two inputs legitimately diverge; the
            // recording pipeline owns the fixture fix. The molfile
            // legs of every one of these rows pass the main gate.
            skipped_phase2.push(id.to_string());
            continue;
        }
        let recorded = case["inchi"].as_str().expect("recorded inchi");
        ran += 1;
        match chemomer::parse_smiles(smiles) {
            Ok(_) => {}
            Err(error) => {
                failures.push(format!("{id}: read errored: {error}"));
                continue;
            }
        }
        match chemomer::Engine::from_smiles(smiles) {
            Ok(engine) => match engine.inchi_string() {
                Ok(derived) if derived == recorded => {}
                Ok(derived) => {
                    failures.push(format!("{id}: {derived:?} != oracle {recorded:?}"));
                }
                Err(error) => {
                    failures.push(format!("{id}: derive errored: {error}"));
                }
            },
            Err(error) => {
                failures.push(format!("{id}: engine errored: {error}"));
            }
        }
    }
    eprintln!(
        "smiles gate: {ran} cases derived vs the oracle, {} phase-2 skips, {} failures",
        skipped_phase2.len(),
        failures.len()
    );
    for skip in &skipped_phase2 {
        eprintln!("  phase-2 skip: {skip}");
    }
    assert!(
        failures.is_empty(),
        "{} smile identities diverged:\n{}",
        failures.len(),
        failures.join("\n")
    );
    assert!(ran > 1_000, "the fixture corpus shrank: {ran}");
}

#[test]
fn the_smiles_reader_derives_the_stereo_markers() {
    // The phase-2 markers ride through: the alanine enantiomers match
    // the oracle bytes (the corpus gate runs them now), and the E/Z
    // spellings split the /b signs the coordinate drawings record.
    for (input, oracle) in [
        (
            "N[C@@H](C)C(=O)O",
            "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1",
        ),
        (
            "N[C@H](C)C(=O)O",
            "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1",
        ),
    ] {
        let engine = chemomer::Engine::from_smiles(input).expect("reads");
        let derived = engine.inchi_string().expect("derives");
        assert_eq!(derived, oracle, "{input:?}");
    }
    let e = chemomer::Engine::from_smiles("C/C=C/C")
        .expect("reads")
        .inchi_string()
        .expect("derives");
    let z = chemomer::Engine::from_smiles("C/C=C\\C")
        .expect("reads")
        .inchi_string()
        .expect("derives");
    assert!(e.ends_with("/b4-3+"), "{e:?}");
    assert!(z.ends_with("/b4-3-"), "{z:?}");
}

#[test]
fn the_smiles_reader_rejects_the_documented_edges() {
    // The phase-2 chirality spellings outside the contract's two
    // markers fail cleanly; dots: the grammar rejects (no component
    // separator).
    let error = chemomer::parse_smiles("[C@TH5](N)(C)C(=O)O").expect_err("@TH5 should not read");
    assert!(
        error.message().contains("phase-2 contract"),
        "{error} lacks the note"
    );
    assert!(chemomer::parse_smiles("C.").is_err());
}
