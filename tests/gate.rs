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
