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

#[test]
fn the_identifier_equals_the_oracle_recording() {
    let dir = suite_dir();
    let fixtures = dir.join("corpus/fixtures");
    let mut cases = Vec::new();
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
    for file in files {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read fixtures"))
                .expect("valid fixtures");
        cases.extend(parsed);
    }
    let triage: serde_json::Value = match fs::read_to_string(dir.join("corpus/triage.json")) {
        Ok(text) => serde_json::from_str(&text).expect("valid triage"),
        Err(_) => serde_json::Value::Null,
    };

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

/// The RInChI gate: the reaction tranche (corpus/fixtures/
/// rinchi-reactions.json, recorded by the RInChI oracle) must
/// reproduce byte-for-byte — the assembled identifier, the RAuxInfo
/// line, all three RInChIKey shapes, and each component's standard
/// InChI. Skips quietly when the tranche has not landed in the
/// corpus checkout (it ships through its own PR); no triage
/// interplay — a reaction fixture either matches or fails.
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

/// The AuxInfo gate: the engine's auxiliary line — canonical
/// numbering, equivalence classes, charges, reversibility echo — must
/// equal the oracle's recording byte-for-byte on the same corpus.
#[test]
fn the_aux_info_equals_the_oracle_recording() {
    let dir = suite_dir();
    let fixtures = dir.join("corpus/fixtures");
    let mut cases = Vec::new();
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
    for file in files {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read fixtures"))
                .expect("valid fixtures");
        cases.extend(parsed);
    }
    // The aux debts live in their own dated ledger: an entry here
    // marks an EXPECTED aux failure only — the identifier gate's
    // triage.json must not learn about them, or their passing
    // identifiers read as retirements.
    let triage: serde_json::Value = match fs::read_to_string(dir.join("corpus/aux-triage.json")) {
        Ok(text) => serde_json::from_str(&text).expect("valid aux triage"),
        Err(_) => serde_json::Value::Null,
    };

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
