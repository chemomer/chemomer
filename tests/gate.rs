//! The gate: the same inchi-tests corpus the reference engine passes.
//! Contract: FULL byte equality — the engine's assembled identifier
//! equals the oracle's recording exactly, every layer present in the
//! oracle's order (formula, /c, /h, /q, /p, stereo, /i), plus the
//! native InChIKey. Counted skips are the triaged fixtures (dated,
//! decoded-reason debts) and a PER-ATOM |charge| > 1 (the
//! GetAtomChargeType guard; multi-charge totals — the multi-onium
//! /q+2..+4 and polyanion /q-2..-4 families — are compared since the
//! 2026-10-31 census). No silent skips; a triaged fixture that
//! passes is a retirement to act on.

use std::fs;
use std::path::PathBuf;

fn suite_dir() -> PathBuf {
    std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"))
}

/// Every identity fixture record, files in name order.
fn load_cases(dir: &std::path::Path) -> Vec<serde_json::Value> {
    let fixtures = dir.join("corpus/fixtures");
    let mut files: Vec<PathBuf> = fs::read_dir(&fixtures)
        .expect("corpus fixtures dir")
        .filter_map(|entry| {
            let path = entry.expect("dir entry").path();
            let name = path
                .file_name()
                .expect("name")
                .to_string_lossy()
                .to_string();
            (name.starts_with("identity") && name.ends_with(".json")).then_some(path)
        })
        .collect();
    files.sort();
    assert!(!files.is_empty(), "no identity fixture files");
    let mut cases = Vec::new();
    for file in files {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read fixtures"))
                .expect("valid fixtures");
        cases.extend(parsed);
    }
    cases
}

/// A dated ledger (triage.json, aux-triage.json); a missing file
/// reads as an empty one.
fn load_ledger(dir: &std::path::Path, name: &str) -> serde_json::Value {
    match fs::read_to_string(dir.join("corpus").join(name)) {
        Ok(text) => serde_json::from_str(&text).expect("valid ledger"),
        Err(_) => serde_json::Value::Null,
    }
}

#[test]
fn the_identifier_equals_the_oracle_recording() {
    let dir = suite_dir();
    let cases = load_cases(&dir);
    let triage = load_ledger(&dir, "triage.json");

    let mut pass = 0;
    let mut skipped = 0;
    let mut retired = 0;
    let mut beyond = 0;
    let mut failures = Vec::new();
    let mut retirements = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        // Every fixture is compared; the triage file marks EXPECTED
        // failures (dated, decoded debts). A triaged fixture that
        // passes is a retirement to act on, never a silent skip.
        let expected_failure = triage.get(id).map(|_| ()).is_some();
        let recorded = case["inchi"].as_str().expect("inchi");

        let molfile = case["molfile"].as_str().expect("molfile");

        match chemomer::Engine::from_molfile(molfile).and_then(|engine| engine.inchi_string()) {
            Ok(assembled) => {
                // Full byte equality: every layer present, in the
                // oracle's order, byte-identical.
                let mut ok = assembled == recorded;
                let mut failures_now = Vec::new();
                if !ok {
                    failures_now.push(format!("assembled {assembled:?} != oracle {recorded:?}"));
                }
                // The native hasher must reproduce every recorded
                // InChIKey from the recorded identifier.
                let key = chemomer::key::inchikey(recorded).expect("keys");
                if key != case["inchikey"].as_str().unwrap_or("") {
                    ok = false;
                    failures_now.push(format!(
                        "key {key} != oracle {}",
                        case["inchikey"].as_str().unwrap_or("")
                    ));
                }
                if ok {
                    if expected_failure {
                        retired += 1;
                        retirements.push(id.to_string());
                    } else {
                        pass += 1;
                    }
                } else if !expected_failure {
                    for line in &failures_now {
                        failures.push(format!("{id}: {line}"));
                    }
                } else {
                    skipped += 1;
                }
            }
            Err(error) => {
                if error.message().contains("beyond the ±1 corpus gate") {
                    beyond += 1;
                } else if !expected_failure {
                    failures.push(format!("{id}: {error}"));
                } else {
                    skipped += 1;
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} mismatch(es):\n{}\n(pass: {pass}, expected-failure skips: {skipped},          RETIRED (triaged, now passing — update the triage file): {retired}, \
         beyond-±1: {beyond}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "gate: pass {pass}, expected-failure skips {skipped}, retired {retired}, \
         beyond-±1 {beyond}, total {}",
        cases.len()
    );
    assert!(
        retirements.is_empty(),
        "{} triaged fixture(s) now PASS — retire their triage entries:\n{}",
        retirements.len(),
        retirements.join("\n")
    );
}

/// The SNon gate: the corpus records every fixture's `-SNon` oracle
/// run (`snon.inchi` — stereo layers withheld, version stays
/// `InChI=1S/`). The engine under [`chemomer::Options::with_no_stereo`]
/// must equal the recording byte-for-byte. Expected failures are the
/// union of the two dated ledgers (the identifier triage and the aux
/// ledger, whose debts share machinery the stereo-off path also
/// touches); a triaged identifier-ledger fixture that passes is a
/// retirement to act on, never a silent skip. Fixtures without the
/// field (a corpus carrying no SNon tranche yet) are skipped, the way
/// the aux gate skips missing aux recordings.
#[test]
fn the_snon_output_equals_the_oracle_recording() {
    let dir = suite_dir();
    let cases = load_cases(&dir);
    let triage = load_ledger(&dir, "triage.json");
    let aux_triage = load_ledger(&dir, "aux-triage.json");

    let mut pass = 0;
    let mut skipped = 0;
    let mut without_recording = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let Some(recorded) = case["snon"]
            .as_object()
            .and_then(|s| s.get("inchi"))
            .and_then(|v| v.as_str())
        else {
            without_recording += 1;
            continue;
        };
        // Expected failure: either dated ledger marks the id. The
        // RETIREMENT check reads the identifier ledger only — an
        // aux-ledger fixture passing the SNon identifier says nothing
        // about its aux debt.
        let identifier_debt = triage.get(id).map(|_| ()).is_some();
        let aux_debt = aux_triage.get(id).map(|_| ()).is_some();
        let expected_failure = identifier_debt || aux_debt;
        let molfile = case["molfile"].as_str().expect("molfile");

        let assembled =
            chemomer::Engine::with_options(molfile, chemomer::Options::default().with_no_stereo())
                .and_then(|engine| engine.inchi_string());
        match assembled {
            Ok(assembled) => {
                if assembled == recorded {
                    if identifier_debt {
                        // Its identifier is still indebted: the
                        // identifier gate owns retirement; the mode
                        // passing is noted, not acted on.
                        skipped += 1;
                    } else {
                        pass += 1;
                    }
                } else if !expected_failure {
                    failures.push(format!(
                        "{id}: assembled {assembled:?} != oracle {recorded:?}"
                    ));
                } else {
                    skipped += 1;
                }
            }
            Err(error) => {
                if error.message().contains("beyond the ±1 corpus gate") {
                    continue;
                }
                if !expected_failure {
                    failures.push(format!("{id}: {error}"));
                } else {
                    skipped += 1;
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} SNon mismatch(es):\n{}\n(SNon pass: {pass}, expected-failure skips: {skipped}, \
         no snon recording: {without_recording}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "SNon gate: pass {pass}, expected-failure skips {skipped}, no recording {without_recording}"
    );
}

/// The FixedH gate: the corpus records every fixture's `-FixedH`
/// oracle run (`fixedh.inchi` — the full identifier with its /f block,
/// plain `InChI=1/` prefix) and its auxiliary line (`fixedh.aux` — the
/// `/F:` pass marker plus repetition markers). The engine under
/// [`chemomer::Options::with_fixed_h`] must equal both recordings
/// byte-for-byte. Expected failures are the union of the two dated
/// ledgers; the retirement check reads the identifier ledger only.
/// Fixtures without the field are skipped, the way the other mode
/// gates skip missing recordings.
#[test]
fn the_fixedh_output_equals_the_oracle_recording() {
    let dir = suite_dir();
    let cases = load_cases(&dir);
    let triage = load_ledger(&dir, "triage.json");
    let aux_triage = load_ledger(&dir, "aux-triage.json");

    let mut pass = 0;
    let mut skipped = 0;
    let mut without_recording = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let Some(recorded) = case["fixedh"]
            .as_object()
            .and_then(|s| s.get("inchi"))
            .and_then(|v| v.as_str())
        else {
            without_recording += 1;
            continue;
        };
        let identifier_debt = triage.get(id).map(|_| ()).is_some();
        let aux_debt = aux_triage.get(id).map(|_| ()).is_some();
        let expected_failure = identifier_debt || aux_debt;
        let molfile = case["molfile"].as_str().expect("molfile");

        let assembled =
            chemomer::Engine::with_options(molfile, chemomer::Options::default().with_fixed_h())
                .and_then(|engine| engine.inchi_string());
        match assembled {
            Ok(assembled) => {
                if assembled == recorded {
                    if identifier_debt {
                        // Its identifier is still indebted: the
                        // identifier gate owns retirement; the mode
                        // passing is noted, not acted on.
                        skipped += 1;
                    } else {
                        pass += 1;
                    }
                } else if !expected_failure {
                    failures.push(format!(
                        "{id}: assembled {assembled:?} != oracle {recorded:?}"
                    ));
                } else {
                    skipped += 1;
                }
            }
            Err(error) => {
                if error.message().contains("beyond the ±1 corpus gate") {
                    continue;
                }
                if !expected_failure {
                    failures.push(format!("{id}: {error}"));
                } else {
                    skipped += 1;
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} FixedH mismatch(es):\n{}\n(FixedH pass: {pass}, expected-failure skips: {skipped}, \
         no fixedh recording: {without_recording}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "FixedH gate: pass {pass}, expected-failure skips {skipped}, \
         no recording {without_recording}"
    );
}

/// The DoNotAddH gate: the corpus records every fixture's
/// `-DoNotAddH` oracle run (`donotaddh.inchi` — only drawn hydrogens
/// exist, H-count codes ignored, the valence-gap inference never
/// runs; version stays `InChI=1S/`). Expected failures and the
/// retirement check follow the SNon gate's union-of-ledgers pattern.
/// The tranche recording lives in the inchi-tests PR branch — until
/// it merges, fixtures without the field skip (the way the aux gate
/// skips missing aux recordings).
#[test]
fn the_do_not_addh_output_equals_the_oracle_recording() {
    let dir = suite_dir();
    let cases = load_cases(&dir);
    let triage = load_ledger(&dir, "triage.json");
    let aux_triage = load_ledger(&dir, "aux-triage.json");
    let dna_triage = load_ledger(&dir, "donotaddh-triage.json");

    let mut pass = 0;
    let mut skipped = 0;
    let mut without_recording = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let Some(recorded) = case["donotaddh"]
            .as_object()
            .and_then(|s| s.get("inchi"))
            .and_then(|v| v.as_str())
        else {
            without_recording += 1;
            continue;
        };
        let identifier_debt = triage.get(id).map(|_| ()).is_some();
        let aux_debt = aux_triage.get(id).map(|_| ()).is_some();
        let dna_debt = dna_triage.get(id).map(|_| ()).is_some();
        let expected_failure = identifier_debt || aux_debt || dna_debt;
        let molfile = case["molfile"].as_str().expect("molfile");

        let assembled = chemomer::Engine::with_options(
            molfile,
            chemomer::Options::default().with_do_not_add_h(),
        )
        .and_then(|engine| engine.inchi_string());
        match assembled {
            Ok(assembled) => {
                if assembled == recorded {
                    if identifier_debt {
                        // Its identifier is still indebted: the
                        // identifier gate owns retirement; the mode
                        // passing is noted, not acted on.
                        skipped += 1;
                    } else {
                        pass += 1;
                    }
                } else if !expected_failure {
                    failures.push(format!(
                        "{id}: assembled {assembled:?} != oracle {recorded:?}"
                    ));
                } else {
                    skipped += 1;
                }
            }
            Err(error) => {
                if error.message().contains("beyond the ±1 corpus gate") {
                    continue;
                }
                if !expected_failure {
                    failures.push(format!("{id}: {error}"));
                } else {
                    skipped += 1;
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} DoNotAddH mismatch(es):\n{}\n(DoNotAddH pass: {pass}, expected-failure skips: {skipped}, \
         no recording: {without_recording}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "DoNotAddH gate: pass {pass}, expected-failure skips {skipped}, no recording {without_recording}"
    );
}

/// The AuxInfo gate: the engine's auxiliary line — canonical
/// numbering, equivalence classes, charges, reversibility echo — must
/// equal the oracle's recording byte-for-byte on the same corpus.
#[test]
fn the_aux_info_equals_the_oracle_recording() {
    let dir = suite_dir();
    let cases = load_cases(&dir);
    // The aux debts live in their own dated ledger: an entry here
    // marks an EXPECTED aux failure only — the identifier gate's
    // triage.json must not learn about them, or their passing
    // identifiers read as retirements.
    let triage = load_ledger(&dir, "aux-triage.json");

    let mut pass = 0;
    let mut skipped = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let expected_failure = triage.get(id).map(|_| ()).is_some();
        let Some(recorded) = case["aux"].as_str() else {
            continue;
        };
        let molfile = case["molfile"].as_str().expect("molfile");
        match chemomer::Engine::from_molfile(molfile).and_then(|engine| engine.aux_string()) {
            Ok(assembled) => {
                if assembled == recorded {
                    pass += 1;
                } else if !expected_failure {
                    failures.push(format!(
                        "{id}: assembled {assembled:?} != oracle {recorded:?}"
                    ));
                } else {
                    skipped += 1;
                }
            }
            Err(error) => {
                if error.message().contains("beyond the ±1 corpus gate") {
                    continue;
                }
                if !expected_failure {
                    failures.push(format!("{id}: {error}"));
                } else {
                    skipped += 1;
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} aux mismatch(es):\n{}\n(aux pass: {pass}, expected-failure skips: {skipped}, total with aux: {})",
        failures.len(),
        failures.join("\n"),
        pass + skipped + failures.len()
    );
    println!("aux gate: pass {pass}, expected-failure skips {skipped}");
}

#[test]
fn the_rinchi_equals_the_oracle_recording() {
    let dir = suite_dir();
    let Ok(text) = fs::read_to_string(dir.join("corpus/fixtures/rinchi-reactions.json")) else {
        return;
    };
    let cases: Vec<serde_json::Value> = serde_json::from_str(&text).expect("valid rinchi fixtures");
    assert!(!cases.is_empty(), "rinchi fixture file is empty");

    let mut pass = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let rxnfile = case["rxnfile"].as_str().expect("rxnfile");
        let reaction = match chemomer::parse_rxnfile(rxnfile) {
            Ok(reaction) => reaction,
            Err(error) => {
                failures.push(format!("{id}: {error}"));
                continue;
            }
        };
        let mut ok = true;
        let check = |what: &str,
                     value: chemomer::Result<String>,
                     recorded: Option<&str>,
                     ok: &mut bool,
                     failures: &mut Vec<String>| {
            let Some(recorded) = recorded else {
                return;
            };
            match value {
                Ok(assembled) if assembled == recorded => {}
                Ok(assembled) => {
                    *ok = false;
                    failures.push(format!("{id}: {what} {assembled:?} != oracle {recorded:?}"));
                }
                Err(error) => {
                    *ok = false;
                    failures.push(format!("{id}: {what} errored: {error}"));
                }
            }
        };
        check(
            "rinchi",
            reaction.rinchi_string(),
            case["rinchi"].as_str(),
            &mut ok,
            &mut failures,
        );
        check(
            "rauxinfo",
            reaction.rinchi_auxinfo(),
            case["rauxinfo"].as_str(),
            &mut ok,
            &mut failures,
        );
        // Keys record without their headers (the corpus convention);
        // the comparison re-attaches them.
        check(
            "rinchikey",
            reaction.rinchi_short_key().map(|key| {
                key.strip_prefix("Short-RInChIKey=")
                    .map(str::to_owned)
                    .unwrap_or(key)
            }),
            case["rinchikey"].as_str(),
            &mut ok,
            &mut failures,
        );
        check(
            "longKey",
            reaction.rinchi_long_key().map(|key| {
                key.strip_prefix("Long-RInChIKey=")
                    .map(str::to_owned)
                    .unwrap_or(key)
            }),
            case["longKey"].as_str(),
            &mut ok,
            &mut failures,
        );
        check(
            "webKey",
            reaction.rinchi_web_key().map(|key| {
                key.strip_prefix("Web-RInChIKey=")
                    .map(str::to_owned)
                    .unwrap_or(key)
            }),
            case["webKey"].as_str(),
            &mut ok,
            &mut failures,
        );
        // The per-component recording (RXN role order): every
        // molecule's standard InChI must reproduce; No-Structure
        // components never reach the engine.
        let roles: Vec<&chemomer::RinchiComponent> = reaction
            .reactants()
            .iter()
            .chain(reaction.products().iter())
            .chain(reaction.agents().iter())
            .collect();
        let recorded_components = case["components"].as_array().cloned().unwrap_or_default();
        if roles.len() != recorded_components.len() {
            ok = false;
            failures.push(format!(
                "{id}: {} component(s) parsed, {} recorded",
                roles.len(),
                recorded_components.len()
            ));
        } else {
            for (component, recorded) in roles.iter().zip(&recorded_components) {
                let recorded_inchi = recorded["inchi"].as_str().expect("component inchi");
                if recorded_inchi == "InChI=1S/" {
                    continue;
                }
                match chemomer::Engine::from_molfile(component.molfile()) {
                    Ok(engine) => match engine.inchi_string() {
                        Ok(inchi) if inchi == recorded_inchi => {}
                        Ok(inchi) => {
                            ok = false;
                            failures.push(format!(
                                "{id}: component inchi {inchi:?} != oracle {recorded_inchi:?}"
                            ));
                        }
                        Err(error) => {
                            ok = false;
                            failures.push(format!("{id}: component errored: {error}"));
                        }
                    },
                    Err(error) => {
                        ok = false;
                        failures.push(format!("{id}: component molfile errored: {error}"));
                    }
                }
            }
        }
        if ok {
            pass += 1;
        }
    }

    assert!(
        failures.is_empty(),
        "{} rinchi mismatch(es):\n{}\n(rinchi pass: {pass}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!("rinchi gate: pass {pass}, total {}", cases.len());
}

/// The decode gate: every recorded RInChI string reads back into a
/// reaction whose re-assembly reproduces the recording byte-for-byte —
/// the identifier line, the RAuxInfo companion, and all three key
/// shapes (the RInChI string is canonical, so decode∘encode is the
/// identity). The decoded components' carried identifiers must equal
/// the recorded per-component InChIs (order-free: the reader sorts
/// within a role, and the byte-exact direction flag proves the role
/// recovery itself — a swapped decode re-encodes with the opposite
/// flag).
#[test]
fn the_rinchi_decode_round_trips() {
    let dir = suite_dir();
    let Ok(text) = fs::read_to_string(dir.join("corpus/fixtures/rinchi-reactions.json")) else {
        println!("decode gate: no rinchi-reactions tranche recorded; skipped");
        return;
    };
    let cases: Vec<serde_json::Value> = serde_json::from_str(&text).expect("valid rinchi fixtures");
    assert!(!cases.is_empty(), "rinchi fixture file is empty");

    let mut pass = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let recorded = case["rinchi"].as_str().expect("rinchi");
        let recorded_aux = case["rauxinfo"].as_str().expect("rauxinfo");
        let mut ok = true;

        // The two-line decode: every recorded output line re-assembles
        // byte-for-byte from the decoded reaction alone.
        let check = |what: &str,
                     assembled: chemomer::Result<String>,
                     expected: &str,
                     ok: &mut bool,
                     failures: &mut Vec<String>| {
            match assembled {
                Ok(value) if value == expected => {}
                Ok(value) => {
                    *ok = false;
                    failures.push(format!(
                        "{id}: decode re-encodes {what} {value:?} != oracle {expected:?}"
                    ));
                }
                Err(error) => {
                    *ok = false;
                    failures.push(format!("{id}: decode {what} errored: {error}"));
                }
            }
        };
        match chemomer::rinchi::reaction_from_rinchi_with_auxinfo(recorded, recorded_aux) {
            Err(error) => {
                ok = false;
                failures.push(format!("{id}: decode errored: {error}"));
            }
            Ok(reaction) => {
                check(
                    "rinchi",
                    reaction.rinchi_string(),
                    recorded,
                    &mut ok,
                    &mut failures,
                );
                check(
                    "rauxinfo",
                    reaction.rinchi_auxinfo(),
                    recorded_aux,
                    &mut ok,
                    &mut failures,
                );
                check(
                    "rinchikey",
                    reaction.rinchi_short_key().map(|key| {
                        key.strip_prefix("Short-RInChIKey=")
                            .map(str::to_owned)
                            .unwrap_or(key)
                    }),
                    case["rinchikey"].as_str().expect("rinchikey"),
                    &mut ok,
                    &mut failures,
                );
                check(
                    "longKey",
                    reaction.rinchi_long_key().map(|key| {
                        key.strip_prefix("Long-RInChIKey=")
                            .map(str::to_owned)
                            .unwrap_or(key)
                    }),
                    case["longKey"].as_str().expect("longKey"),
                    &mut ok,
                    &mut failures,
                );
                check(
                    "webKey",
                    reaction.rinchi_web_key().map(|key| {
                        key.strip_prefix("Web-RInChIKey=")
                            .map(str::to_owned)
                            .unwrap_or(key)
                    }),
                    case["webKey"].as_str().expect("webKey"),
                    &mut ok,
                    &mut failures,
                );

                // The decoded components carry the recorded per-role
                // identities (No-Structures included, their bare
                // identifier recorded like the corpus does).
                let mut decoded: Vec<String> = reaction
                    .reactants()
                    .iter()
                    .chain(reaction.products().iter())
                    .chain(reaction.agents().iter())
                    .map(|component| component.inchi_string().expect("decoded inchi"))
                    .collect();
                let mut recorded_components: Vec<String> = case["components"]
                    .as_array()
                    .expect("components")
                    .iter()
                    .map(|component| {
                        component["inchi"]
                            .as_str()
                            .expect("component inchi")
                            .to_owned()
                    })
                    .collect();
                decoded.sort();
                recorded_components.sort();
                if decoded != recorded_components {
                    ok = false;
                    failures.push(format!(
                        "{id}: decoded components {decoded:?} != oracle {recorded_components:?}"
                    ));
                }
            }
        }

        // The identifier-only decode (no RAuxInfo companion at hand)
        // reads the same reaction.
        match chemomer::reaction_from_rinchi(recorded) {
            Err(error) => {
                ok = false;
                failures.push(format!("{id}: bare decode errored: {error}"));
            }
            Ok(reaction) => check(
                "bare rinchi",
                reaction.rinchi_string(),
                recorded,
                &mut ok,
                &mut failures,
            ),
        }

        if ok {
            pass += 1;
        }
    }

    assert!(
        failures.is_empty(),
        "{} decode mismatch(es):\n{}\n(decode pass: {pass}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!("decode gate: pass {pass}, total {}", cases.len());
}

/// The key-extension gate: the corpus tranche records the oracle's
/// 256-bit hash extensions (`inchikey_xhash1`, `inchikey_xhash2`) and
/// the merged `InChIHash=` line (`inchikey_merged`) for a family-
/// spanning subset of fixtures. The native computation reads the SAME
/// digests the truncated key reads — every recorded field must equal
/// our extension byte-for-byte, and the recorded standard key must
/// still match the native key over the recorded identifier (the
/// per-record invariant the recording script checks). Fixtures without
/// the fields skip — no triage interplay: the gate is hash-level over
/// the recorded identifiers, so the engine's milestone state never
/// enters it.
#[test]
fn the_key_extensions_match_the_oracle() {
    let dir = suite_dir();
    let cases = load_cases(&dir);

    let mut pass = 0;
    let mut without_recording = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let Some(recorded_xhash1) = case["inchikey_xhash1"].as_str() else {
            without_recording += 1;
            continue;
        };
        let id = case["id"].as_str().expect("id");
        let recorded_inchi = case["inchi"].as_str().expect("inchi");
        let recorded_xhash2 = case["inchikey_xhash2"].as_str().expect("xhash2 field");
        let recorded_merged = case["inchikey_merged"].as_str().expect("merged field");
        let recorded_key = case["inchikey"].as_str().expect("standard key");

        let (key, xhash1, xhash2) = match chemomer::key::key_extensions(recorded_inchi) {
            Ok(parts) => parts,
            Err(error) => {
                failures.push(format!("{id}: {error}"));
                continue;
            }
        };
        // The invariant: the extensions ride the same digests, so the
        // standard key must come out identical to the recorded one.
        if key != recorded_key {
            failures.push(format!("{id}: key {key} != recorded {recorded_key}"));
            continue;
        }
        if xhash1 != recorded_xhash1 {
            failures.push(format!("{id}: xhash1 {xhash1} != oracle {recorded_xhash1}"));
        }
        if xhash2 != recorded_xhash2 {
            failures.push(format!("{id}: xhash2 {xhash2} != oracle {recorded_xhash2}"));
        }
        let merged = chemomer::key::merged_line(&key, Some(&xhash1), Some(&xhash2));
        if merged != recorded_merged {
            failures.push(format!(
                "{id}: merged {merged:?} != oracle {recorded_merged:?}"
            ));
        }
        if xhash1 == recorded_xhash1 && xhash2 == recorded_xhash2 && merged == recorded_merged {
            pass += 1;
        }
    }

    assert!(
        failures.is_empty(),
        "{} key-extension mismatch(es):\n{}\n(key-extension pass: {pass}, \
         no recording: {without_recording}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "key-extension gate: pass {pass}, no recording {without_recording}, total {}",
        cases.len()
    );
}

/// The FixedH aux gate: the auxiliary line under `-FixedH` — the `/F:`
/// pass marker and the repetition markers (`m` = same as main) — must
/// equal the oracle's `fixedh.aux` recording byte-for-byte.
#[test]
fn the_fixedh_aux_equals_the_oracle_recording() {
    let dir = suite_dir();
    let cases = load_cases(&dir);
    let triage = load_ledger(&dir, "triage.json");
    let aux_triage = load_ledger(&dir, "aux-triage.json");

    let mut pass = 0;
    let mut skipped = 0;
    let mut without_recording = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let Some(recorded) = case["fixedh"]
            .as_object()
            .and_then(|s| s.get("aux"))
            .and_then(|v| v.as_str())
        else {
            without_recording += 1;
            continue;
        };
        let expected_failure =
            triage.get(id).map(|_| ()).is_some() || aux_triage.get(id).map(|_| ()).is_some();
        let molfile = case["molfile"].as_str().expect("molfile");
        match chemomer::Engine::with_options(molfile, chemomer::Options::default().with_fixed_h())
            .and_then(|engine| engine.aux_string())
        {
            Ok(assembled) => {
                if assembled == recorded {
                    pass += 1;
                } else if !expected_failure {
                    failures.push(format!(
                        "{id}: assembled {assembled:?} != oracle {recorded:?}"
                    ));
                } else {
                    skipped += 1;
                }
            }
            Err(error) => {
                if error.message().contains("beyond the ±1 corpus gate") {
                    continue;
                }
                if !expected_failure {
                    failures.push(format!("{id}: {error}"));
                } else {
                    skipped += 1;
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} FixedH aux mismatch(es):\n{}\n(FixedH aux pass: {pass}, expected-failure skips: {skipped}, \
         no fixedh recording: {without_recording}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "FixedH aux gate: pass {pass}, expected-failure skips {skipped}, \
         no recording {without_recording}"
    );
}

/// The RecMet gate: the corpus records every metal fixture's
/// `-RecMet` oracle run — the main identifier under the non-standard
/// prefix with the reconnected `/r` layers appended (only when a
/// coordination disconnect actually happened), plus the aux line
/// carrying the embedded `/R:` short-aux block of the reconnected
/// pass. Expected failures are the union of the two dated ledgers.
/// Fixtures without the field (no metal anywhere) are skipped.
#[test]
fn the_recmet_output_equals_the_oracle_recording() {
    let dir = suite_dir();
    let cases = load_cases(&dir);
    let triage = load_ledger(&dir, "triage.json");
    let aux_triage = load_ledger(&dir, "aux-triage.json");

    let mut pass = 0;
    let mut skipped = 0;
    let mut without_recording = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let Some(recmet) = case["recmet"].as_object() else {
            without_recording += 1;
            continue;
        };
        let expected_failure =
            triage.get(id).map(|_| ()).is_some() || aux_triage.get(id).map(|_| ()).is_some();
        let molfile = case["molfile"].as_str().expect("molfile");
        let assembled =
            chemomer::Engine::with_options(molfile, chemomer::Options::default().with_rec_met());
        for (kind, recorded) in [
            ("inchi", recmet.get("inchi").and_then(|v| v.as_str())),
            ("aux", recmet.get("aux").and_then(|v| v.as_str())),
        ] {
            let Some(recorded) = recorded else {
                continue;
            };
            let line = match kind {
                "inchi" => assembled.as_ref().map(|engine| engine.recmet_string()),
                _ => assembled.as_ref().map(|engine| engine.aux_string()),
            };
            let line = match line {
                Ok(line) => line,
                Err(error) => Err(error.clone()),
            };
            match line {
                Ok(line) => {
                    if line == recorded {
                        pass += 1;
                    } else if !expected_failure {
                        failures.push(format!(
                            "{id}: recmet {kind} assembled {line:?} != oracle {recorded:?}"
                        ));
                    } else {
                        skipped += 1;
                    }
                }
                Err(error) => {
                    if error.message().contains("beyond the ±1 corpus gate") {
                        continue;
                    }
                    if !expected_failure {
                        failures.push(format!("{id}: recmet {kind} {error}"));
                    } else {
                        skipped += 1;
                    }
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} RecMet mismatch(es):\n{}\n(RecMet pass: {pass}, expected-failure skips: {skipped}, \
         no recmet recording: {without_recording}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "RecMet gate: pass {pass}, expected-failure skips {skipped}, \
         no recording {without_recording}"
    );
}

/// The reverse-reconstruction gate: for every recorded fixture, the
/// oracle's InChI2Struct reconstruction — its SDfile, verbatim —
/// parsed and re-forwarded by this engine must reproduce the oracle's
/// re-forward InChI byte-for-byte, and this engine's own
/// reconstruction of the source InChI must forward to the same
/// identifier. The reference's mobile-H/charge restore is
/// deliberately lossy (ichirvr1..7.c); the ground truth is the
/// oracle's own recorded reconstruction, never a chemistry argument.
///
/// Stereo recordings relax only the sdf leg, and only for the
/// reference's own serialization loss: its InChI2Struct parities live
/// in memory (ichirvr2.c sets 0D parity descriptors) and its SDfile
/// writer emits neither wedges nor coordinates (every atom at the
/// origin), so no reader can recover the stereo layers from that
/// file. For a provably 0D recording the sdf leg witnesses the
/// non-stereo layers; the stereo layers are witnessed by the reverse
/// leg, which forwards OUR reconstruction — wedges and drawn
/// coordinates included — to the FULL recorded identifier.
/// Skipped entirely when the tranche file is absent (no triage
/// interplay — the reverse slice records its own corpus).
#[test]
fn the_reverse_reconstruction_matches_the_oracle() {
    let dir = suite_dir();
    let Ok(text) = fs::read_to_string(dir.join("corpus/fixtures/reverse-reconstruction.json"))
    else {
        println!("reverse gate: no reverse-reconstruction tranche recorded; skipped");
        return;
    };
    let cases: Vec<serde_json::Value> = serde_json::from_str(&text).expect("valid reverse tranche");

    /// Forwards a structure through the engine path and returns the
    /// assembled identifier.
    fn forward(structure: chemomer::Structure) -> Result<String, String> {
        let molfile = chemomer::reverse::structure_to_molfile(&structure);
        chemomer::Engine::from_molfile(&molfile)
            .and_then(|engine| engine.inchi_string())
            .map_err(|error| error.to_string())
    }

    /// Whether the identifier carries stereo segments (/b /t /m /s).
    fn has_stereo_layers(inchi: &str) -> bool {
        inchi
            .split('/')
            .skip(1)
            .any(|segment| matches!(segment.as_bytes().first(), Some(b'b' | b't' | b'm' | b's')))
    }

    /// The identifier without its stereo segments — the layers a 0D
    /// serialization can still witness.
    fn strip_stereo_layers(inchi: &str) -> String {
        let mut out = String::new();
        for (index, segment) in inchi.split('/').enumerate() {
            if index > 0 && matches!(segment.as_bytes().first(), Some(b'b' | b't' | b'm' | b's')) {
                continue;
            }
            if index > 0 {
                out.push('/');
            }
            out.push_str(segment);
        }
        out
    }

    /// Whether the SDfile serializes a 0D reconstruction: every atom
    /// at the origin and no bond stereo code — the shape in which the
    /// oracle's own writer loses its in-memory parities.
    fn serializes_zero_dimensional(sdf: &str) -> bool {
        let lines: Vec<&str> = sdf.lines().collect();
        if lines.len() < 5 {
            return false;
        }
        let na: usize = lines[3].get(0..3).unwrap_or("").trim().parse().unwrap_or(0);
        let nb: usize = lines[3].get(3..6).unwrap_or("").trim().parse().unwrap_or(0);
        let mut atoms = 0;
        for line in lines.iter().take(4 + na).skip(4) {
            let (Some(x), Some(y), Some(z)) = (line.get(0..10), line.get(10..20), line.get(20..30))
            else {
                return false;
            };
            let zero = |field: &str| {
                field
                    .trim()
                    .parse::<f64>()
                    .map(|value| value == 0.0)
                    .unwrap_or(false)
            };
            if !(zero(x) && zero(y) && zero(z)) {
                return false;
            }
            atoms += 1;
        }
        if atoms != na {
            return false;
        }
        lines.iter().skip(4 + na).take(nb).all(|line| {
            line.get(9..12)
                .map(|field| field.trim() == "0")
                .unwrap_or(true)
        })
    }

    let mut pass = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let inchi = case["inchi"].as_str().expect("inchi");
        let sdf = case["sdf"].as_str().expect("sdf");
        let recorded_forward = case["forward_inchi"].as_str().expect("forward_inchi");

        // The oracle's own reconstruction re-forwards identically —
        // modulo the reference's own 0D stereo serialization loss
        // (documented above; the assembled string then carries no
        // stereo layers at all, and the stripped comparison still
        // pins every other layer).
        if let Err(line) = chemomer::Engine::from_molfile(sdf)
            .and_then(|engine| engine.inchi_string())
            .map_err(|error| error.to_string())
            .and_then(|assembled| {
                if assembled == recorded_forward {
                    return Ok(());
                }
                if has_stereo_layers(recorded_forward)
                    && !has_stereo_layers(&assembled)
                    && serializes_zero_dimensional(sdf)
                    && strip_stereo_layers(&assembled) == strip_stereo_layers(recorded_forward)
                {
                    return Ok(());
                }
                Err(format!(
                    "sdf forward {assembled:?} != oracle {recorded_forward:?}"
                ))
            })
        {
            failures.push(format!("{id}: {line}"));
            continue;
        }

        // Our reconstruction of the source InChI forwards to the same
        // identifier — the reverse direction's contract.
        match chemomer::structure_from_inchi(inchi) {
            Err(error) => failures.push(format!("{id}: reverse error: {error}")),
            Ok(structure) => match forward(structure) {
                Err(line) => failures.push(format!("{id}: reverse forward error: {line}")),
                Ok(assembled) => {
                    if assembled != recorded_forward {
                        failures.push(format!(
                            "{id}: reverse forward {assembled:?} != oracle {recorded_forward:?}"
                        ));
                        continue;
                    }
                    pass += 1;
                }
            },
        }
    }

    assert!(
        failures.is_empty(),
        "{} reverse-reconstruction mismatch(es):\n{}\n(pass: {pass}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!("reverse gate: pass {pass}, total {}", cases.len());
}

/// The experimental tautomerism gate: the corpus's tautomer-extremes
/// tranche records every hand-built tautomer-rich fixture under each
/// flag combination (ket, t15, ket15t, pt22, pt16, pt06, pt39, pt13,
/// pt18); the engine under the matching options must equal each
/// recording byte-for-byte. The tranche file is optional — a corpus
/// without it skips (no triage interplay: these recordings are all
/// flag-path outputs the standard pipeline never produces).
#[test]
fn the_tautomer_extremes_match_the_oracle() {
    let dir = suite_dir();
    let Ok(text) = fs::read_to_string(dir.join("corpus/fixtures/tautomer-extremes.json")) else {
        return;
    };
    let cases: Vec<serde_json::Value> =
        serde_json::from_str(&text).expect("valid tautomer-extremes fixtures");
    assert!(
        !cases.is_empty(),
        "the tautomer-extremes file carries no fixtures"
    );

    /// One flag mode: the recording key and the option builder.
    type FlagMode = (&'static str, fn(chemomer::Options) -> chemomer::Options);

    const MODES: [FlagMode; 9] = [
        ("ket", |options| options.with_ket()),
        ("t15", |options| options.with_15t()),
        ("ket15t", |options| options.with_ket().with_15t()),
        ("pt22", |options| options.with_pt_22_00()),
        ("pt16", |options| options.with_pt_16_00()),
        ("pt06", |options| options.with_pt_06_00()),
        ("pt39", |options| options.with_pt_39_00()),
        ("pt13", |options| options.with_pt_13_00()),
        ("pt18", |options| options.with_pt_18_00()),
    ];

    let mut checked = 0usize;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let molfile = case["molfile"].as_str().expect("molfile");
        for (mode, flag) in MODES {
            let recorded = case[mode]["inchi"].as_str().expect("mode inchi recording");
            let options = flag(chemomer::Options::default());
            match chemomer::Engine::with_options(molfile, options)
                .and_then(|engine| engine.inchi_string())
            {
                Ok(assembled) => {
                    checked += 1;
                    if assembled != recorded {
                        failures.push(format!(
                            "{id} [{mode}]: assembled {assembled:?} != oracle {recorded:?}"
                        ));
                    }
                }
                Err(error) => failures.push(format!("{id} [{mode}]: {error}")),
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} tautomer-extremes mismatch(es):\n{}\n(mode checks: {checked}, fixtures: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "tautomer-extremes gate: {checked} mode checks over {} fixtures",
        cases.len()
    );
}

/// One identity tranche walk, shared by the input-breadth gates: every
/// `{id, molfile, inchi, inchikey}` fixture must reproduce the oracle's
/// recorded identifier byte-for-byte through `Engine::from_molfile`,
/// with the native key over the recorded identifier (the identifier
/// gate's per-record invariant). A missing tranche file reads as an
/// empty one — the gate skips, the way the rinchi gate skips a missing
/// tranche. Beyond-±1 fixtures (the per-atom GetAtomChargeType bound)
/// count separately, never as failures.
fn check_identity_tranche(dir: &std::path::Path, file: &str, label: &str) {
    let Ok(text) = fs::read_to_string(dir.join("corpus/fixtures").join(file)) else {
        println!("{label} gate: no tranche file, skipped");
        return;
    };
    let cases: Vec<serde_json::Value> = serde_json::from_str(&text)
        .unwrap_or_else(|error| panic!("valid {file} fixtures: {error}"));
    assert!(!cases.is_empty(), "the {file} file carries no fixtures");

    let mut pass = 0;
    let mut beyond = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let molfile = case["molfile"].as_str().expect("molfile");
        let recorded = case["inchi"].as_str().expect("inchi");
        let recorded_key = case["inchikey"].as_str().expect("inchikey");
        match chemomer::Engine::from_molfile(molfile).and_then(|engine| engine.inchi_string()) {
            Ok(assembled) => {
                if assembled == recorded {
                    let key = chemomer::key::inchikey(recorded).expect("keys");
                    if key == recorded_key {
                        pass += 1;
                    } else {
                        failures.push(format!("{id}: key {key} != oracle {recorded_key}"));
                    }
                } else {
                    failures.push(format!(
                        "{id}: assembled {assembled:?} != oracle {recorded:?}"
                    ));
                }
            }
            Err(error) => {
                if error.message().contains("beyond the ±1 corpus gate") {
                    beyond += 1;
                } else {
                    failures.push(format!("{id}: {error}"));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} {label} mismatch(es):\n{}\n({label} pass: {pass}, beyond-±1: {beyond}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "{label} gate: pass {pass}, beyond-±1 {beyond}, total {}",
        cases.len()
    );
}

/// The V3000 input gate: the corpus's v3000-molfiles tranche records a
/// sample of existing fixtures re-authored as V3000 molfiles (the
/// CHG=/MASS= attributes, CFG= wedges, free-form V30 lines), each with
/// the oracle's identifier over that V3000 input — cross-checked at
/// recording time against the source fixture's V2000 identity. The
/// engine's reader dispatches on the counts-line version tag, so the
/// same `Engine::from_molfile` path must reproduce both recordings
/// byte-for-byte.
#[test]
fn the_v3000_input_matches_the_oracle() {
    let dir = suite_dir();
    check_identity_tranche(&dir, "v3000-molfiles.json", "v3000");
}

/// The LargeMolecules gate: the corpus's large-molecules tranche
/// records alkane chains past the reference's normally-allowed input
/// size (ichisize.h:51, the 1024-atom cap its reader enforces without
/// -LargeMolecules) with the oracle's `-LargeMolecules` identifiers —
/// the beta `InChI=1B` flavor the reference writes for every structure
/// past that bound. Our reader carries no atom cap and the assembler
/// applies the same atom-count rule, so plain `Engine::from_molfile`
/// reproduces both recordings.
#[test]
fn the_large_molecules_match_the_oracle() {
    let dir = suite_dir();
    check_identity_tranche(&dir, "large-molecules.json", "large-molecules");
}

/// The SDfile batch gate: the corpus's sdf-batch tranche records SDFs
/// verbatim with the oracle's per-record identifiers, in record order.
/// `parse_sdf` splits the records on the $$$$ terminator lines and each
/// record pipes through `Engine::from_molfile` standing alone — the
/// engine-side piece of a batch read.
#[test]
fn the_sdf_batch_matches_the_oracle() {
    let dir = suite_dir();
    let Ok(text) = fs::read_to_string(dir.join("corpus/fixtures/sdf-batch.json")) else {
        return;
    };
    let cases: Vec<serde_json::Value> =
        serde_json::from_str(&text).expect("valid sdf-batch fixtures");
    assert!(!cases.is_empty(), "the sdf-batch file carries no fixtures");

    let mut pass = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let sdf = case["sdf"].as_str().expect("sdf");
        let recorded: Vec<&str> = case["inchis"]
            .as_array()
            .expect("inchis")
            .iter()
            .map(|value| value.as_str().expect("inchi"))
            .collect();
        let records = chemomer::parse_sdf(sdf);
        if records.len() != recorded.len() {
            failures.push(format!(
                "{id}: parse_sdf split {} record(s), oracle recorded {}",
                records.len(),
                recorded.len()
            ));
            continue;
        }
        let mut ok = true;
        for (index, (record, recorded_inchi)) in records.iter().zip(&recorded).enumerate() {
            match chemomer::Engine::from_molfile(record).and_then(|engine| engine.inchi_string()) {
                Ok(assembled) if assembled == *recorded_inchi => {}
                Ok(assembled) => {
                    ok = false;
                    failures.push(format!(
                        "{id}: record {index} assembled {assembled:?} != oracle {recorded_inchi:?}"
                    ));
                }
                Err(error) => {
                    ok = false;
                    failures.push(format!("{id}: record {index} errored: {error}"));
                }
            }
        }
        if ok {
            pass += 1;
        }
    }
    assert!(
        failures.is_empty(),
        "{} sdf-batch mismatch(es):\n{}\n(sdf-batch pass: {pass}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!("sdf-batch gate: pass {pass}, total {}", cases.len());
}

/// The polymer gate: the corpus's polymers tranche records CRU
/// molfiles with the oracle's per-mode identifier and aux line —
/// `polymers` (the modern treatment: the /z layer, the InChI=1B
/// flavor, the /Z: aux echo and the /U/ unedited-aux appendix),
/// `polymers105` (the legacy cyclized flavor with the Zz hidden),
/// the frame-shift variants, `npzz` — plus the recorded rejections
/// (an `error` message in place of the identifier). Each recorded
/// mode pipes through `Engine::with_options`; a missing tranche file
/// reads as an empty one — the gate skips.
#[test]
fn the_polymer_layers_match_the_oracle() {
    let dir = suite_dir();
    let Ok(text) = fs::read_to_string(dir.join("corpus/fixtures/polymers.json")) else {
        println!("polymer gate: no polymers.json tranche recorded; skipped");
        return;
    };
    let cases: Vec<serde_json::Value> = serde_json::from_str(&text)
        .unwrap_or_else(|error| panic!("valid polymers.json fixtures: {error}"));
    assert!(
        !cases.is_empty(),
        "the polymers.json file carries no fixtures"
    );

    fn mode_options(mode: &str) -> chemomer::Options {
        use chemomer::options::{FrameShift, PolymersMode};
        let mut options = chemomer::Options::default();
        for part in mode.split('+') {
            match part {
                "polymers" => options.polymers = PolymersMode::Modern,
                "polymers105" => options.polymers = PolymersMode::Legacy,
                "npzz" => options.npzz = true,
                "fs-none" => options.frame_shift = FrameShift::None,
                "fs-stars" => options.frame_shift = FrameShift::Opened,
                "foldcru" => options.fold_cru = true,
                _ => {}
            }
        }
        options
    }

    let mut pass = 0;
    let mut rejects = 0;
    let mut failures = Vec::new();
    let mut total = 0;
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let molfile = case["molfile"].as_str().expect("molfile");
        let modes = case["modes"].as_object().expect("modes");
        for (mode, recording) in modes {
            total += 1;
            let recorded_error = recording.get("error").and_then(|v| v.as_str());
            let Some(recorded_inchi) = recording.get("inchi").and_then(|v| v.as_str()) else {
                // The recorded rejection: the engine must refuse the
                // structure with the oracle's message.
                let expected = recorded_error.expect("error or inchi recording");
                match chemomer::Engine::with_options(molfile, mode_options(mode)) {
                    Ok(_) => failures.push(format!(
                        "{id} [{mode}]: assembled an identifier; the oracle rejects with {expected:?}"
                    )),
                    Err(error) if error.message() == expected => rejects += 1,
                    Err(error) => failures.push(format!(
                        "{id} [{mode}]: error {:?} != oracle {expected:?}",
                        error.message()
                    )),
                }
                continue;
            };
            let recorded_aux = recording.get("aux").and_then(|v| v.as_str());
            match chemomer::Engine::with_options(molfile, mode_options(mode)) {
                Ok(engine) => match (engine.inchi_string(), engine.aux_string()) {
                    (Ok(assembled), Ok(aux)) => {
                        if assembled != recorded_inchi {
                            failures.push(format!(
                                "{id} [{mode}]: assembled {assembled:?} != oracle {recorded_inchi:?}"
                            ));
                        } else if let Some(recorded) = recorded_aux {
                            if aux != recorded {
                                failures.push(format!(
                                    "{id} [{mode}]: aux {aux:?} != oracle {recorded:?}"
                                ));
                            } else {
                                pass += 1;
                            }
                        } else {
                            pass += 1;
                        }
                    }
                    (Ok(_), Err(error)) => {
                        failures.push(format!("{id} [{mode}]: aux errored: {error}"))
                    }
                    (Err(error), _) => failures.push(format!("{id} [{mode}]: {error}")),
                },
                Err(error) => failures.push(format!("{id} [{mode}]: {error}")),
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} polymer mismatch(es):\n{}\n(polymer pass: {pass}, oracle rejections reproduced: {rejects}, total: {total})",
        failures.len(),
        failures.join("\n"),
    );
    println!("polymer gate: pass {pass}, rejections {rejects}, total {total}");
}

/// The DoNotAddH aux gate: the auxiliary line under `-DoNotAddH` —
/// the canonical numbering, equivalence and CRV fields over the
/// un-augmented structure, plus the `/rA` echo whose unusual-valence
/// digit reads the INPUT's hydrogen state (the reference's
/// `WriteOrigAtoms` over orig_at_data, whose num_H the bDoNotAddH
/// read left at the declared count — zero corpus-wide) — must equal
/// the oracle's `donotaddh.aux` recording byte-for-byte. Expected
/// failures are the union of the four dated ledgers: the identifier
/// triage, the aux ledger, the DoNotAddH identifier ledger (whose
/// debts span both sides of the mode), and the DoNotAddH AUX ledger —
/// the CRV/numbering shadow of the reference's zero-H BNS
/// π-restoration, a debt of its own. Fixtures without the field are
/// skipped, the way the other mode gates skip missing recordings.
#[test]
fn the_do_not_addh_aux_equals_the_oracle_recording() {
    let dir = suite_dir();
    let cases = load_cases(&dir);
    let triage = load_ledger(&dir, "triage.json");
    let aux_triage = load_ledger(&dir, "aux-triage.json");
    let dna_triage = load_ledger(&dir, "donotaddh-triage.json");
    let dna_aux_triage = load_ledger(&dir, "donotaddh-aux-triage.json");

    let mut pass = 0;
    let mut skipped = 0;
    let mut without_recording = 0;
    let mut failures = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id");
        let Some(recorded) = case["donotaddh"]
            .as_object()
            .and_then(|s| s.get("aux"))
            .and_then(|v| v.as_str())
        else {
            without_recording += 1;
            continue;
        };
        let expected_failure = triage.get(id).map(|_| ()).is_some()
            || aux_triage.get(id).map(|_| ()).is_some()
            || dna_triage.get(id).map(|_| ()).is_some()
            || dna_aux_triage.get(id).map(|_| ()).is_some();
        let molfile = case["molfile"].as_str().expect("molfile");
        match chemomer::Engine::with_options(
            molfile,
            chemomer::Options::default().with_do_not_add_h(),
        )
        .and_then(|engine| engine.aux_string())
        {
            Ok(assembled) => {
                if assembled == recorded {
                    pass += 1;
                } else if !expected_failure {
                    failures.push(format!(
                        "{id}: assembled {assembled:?} != oracle {recorded:?}"
                    ));
                } else {
                    skipped += 1;
                }
            }
            Err(error) => {
                if error.message().contains("beyond the ±1 corpus gate") {
                    continue;
                }
                if !expected_failure {
                    failures.push(format!("{id}: {error}"));
                } else {
                    skipped += 1;
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} DoNotAddH aux mismatch(es):\n{}\n(DoNotAddH aux pass: {pass}, expected-failure skips: {skipped}, \
         no recording: {without_recording}, total: {})",
        failures.len(),
        failures.join("\n"),
        cases.len()
    );
    println!(
        "DoNotAddH aux gate: pass {pass}, expected-failure skips {skipped}, \
         no recording {without_recording}"
    );
}
