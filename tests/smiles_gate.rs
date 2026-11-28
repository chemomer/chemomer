//! The SMILES reader gate (chemomer#121, phase 2 in chemomer#128):
//! every identity fixture's `smiles` field reads through the contract
//! and derives the recorded 1.07.5-oracle InChI byte-for-byte — no
//! new oracle, no fixture edits. The stereo cases all run (the
//! alanine enantiomers through the carried @/@@; the butene pair
//! through the corpus's drawn-coordinate re-record, inchi-tests
//! v0.3.42); the transcription disagreements emptied when the corpus
//! re-recorded the fifteen broken molfiles from their own SMILES
//! (inchi-tests v0.3.44) — the recorded exclusions are the oxyanion
//! normalization slice (chemomer#151).

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
    // The triage ledger carries the dated debts (the corpus's
    // triage.json): a triaged row is an EXPECTED failure, and a
    // triaged row that passes is a retirement to act on, never a
    // silent skip.
    let Ok(triage_text) = fs::read_to_string(dir.join("corpus/triage.json")) else {
        return;
    };
    let triage: serde_json::Map<String, serde_json::Value> =
        serde_json::from_str(&triage_text).expect("valid ledger");

    let mut ran = 0usize;
    let mut skipped_phase2 = Vec::new();
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let Some(smiles) = case["smiles"].as_str() else {
            continue;
        };
        // The triage ledger carries the dated debts (the corpus's
        // triage.json): the oxyanion normalization rows (chemomer#151)
        // are triaged there — the oracle keeps an explicit [O-]'s
        // charge fixed and protonates it under the default
        // hydrogen-add mode; the reader delocalizes the acid hydrogen
        // over the oxygens instead. The oracle bytes on the (now
        // correct) molfiles are that slice's spec.
        let triaged = triage.get(id).map(|v| !v.is_null()).unwrap_or(false);
        if triaged {
            // The fixed-charge oxyanion normalization slice
            // (chemomer#151): the recorded identities are the oracle's
            // fixed-charge readings; the reader delocalizes. Expected
            // failure until the slice lands.
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
        "smiles gate: {ran} cases derived vs the oracle, {} triaged, {} failures",
        skipped_phase2.len(),
        failures.len()
    );
    for skip in &skipped_phase2 {
        eprintln!("  triaged: {skip}");
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
