// Probe: the retirement candidates — every multi/metal ledgered
// fixture against EVERY recording it carries (the ledger interplay:
// removing an aux entry exposes the fixture to all the mode gates).
use std::fs;

fn suite_dir() -> std::path::PathBuf {
    std::env::var("INCHI_TESTS")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| std::path::PathBuf::from("../inchi-tests"))
}

fn main() {
    let dir = suite_dir();
    let fixtures = dir.join("corpus/fixtures");
    let mut files: Vec<std::path::PathBuf> = fs::read_dir(&fixtures)
        .expect("corpus fixtures dir")
        .filter_map(|e| {
            let p = e.expect("dir entry").path();
            let n = p.file_name().expect("name").to_string_lossy().to_string();
            (n.starts_with("identity") && n.ends_with(".json")).then_some(p)
        })
        .collect();
    files.sort();
    let ledger: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(dir.join("corpus/aux-triage.json")).expect("aux ledger"),
    )
    .expect("valid ledger");

    let mut checked = 0usize;
    let mut clean: Vec<String> = Vec::new();
    let mut dirty = 0usize;
    for file in &files {
        let cases: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(file).expect("read")).expect("valid");
        for case in &cases {
            let id = case["id"].as_str().expect("id");
            let Some(entry) = ledger.get(id) else {
                continue;
            };
            if !id.contains("/multi-") && !id.contains("/metal-") {
                continue;
            }
            let _ = entry;
            let molfile = case["molfile"].as_str().expect("molfile");
            let mut problems: Vec<String> = Vec::new();

            // 1. the standard identifier + key (always compared).
            let recorded = case["inchi"].as_str().expect("inchi");
            match chemomer::Engine::from_molfile(molfile).and_then(|e| e.inchi_string()) {
                Ok(a) if a == recorded => {}
                Ok(a) => problems.push(format!("std inchi {a:?} != {recorded:?}")),
                Err(e) => problems.push(format!("std inchi err {e}")),
            }
            let key = chemomer::key::inchikey(recorded).expect("keys");
            if key != case["inchikey"].as_str().unwrap_or("") {
                problems.push(format!(
                    "key {key} != oracle {}",
                    case["inchikey"].as_str().unwrap_or("")
                ));
            }

            // 2. the standard aux.
            if let Some(recorded) = case["aux"].as_str() {
                match chemomer::Engine::from_molfile(molfile).and_then(|e| e.aux_string()) {
                    Ok(a) if a == recorded => {}
                    Ok(a) => problems.push(format!("aux {a:?} != {recorded:?}")),
                    Err(e) => problems.push(format!("aux err {e}")),
                }
            }

            // 3. SNon.
            if let Some(recorded) = case["snon"]
                .as_object()
                .and_then(|s| s.get("inchi"))
                .and_then(|v| v.as_str())
            {
                match chemomer::Engine::with_options(
                    molfile,
                    chemomer::Options::default().with_no_stereo(),
                )
                .and_then(|e| e.inchi_string())
                {
                    Ok(a) if a == recorded => {}
                    Ok(a) => problems.push(format!("snon {a:?} != {recorded:?}")),
                    Err(e) => problems.push(format!("snon err {e}")),
                }
            }

            // 4. FixedH identifier + aux.
            if let Some(recorded) = case["fixedh"]
                .as_object()
                .and_then(|s| s.get("inchi"))
                .and_then(|v| v.as_str())
            {
                let engine = chemomer::Engine::with_options(
                    molfile,
                    chemomer::Options::default().with_fixed_h(),
                );
                match engine.as_ref().ok().and_then(|e| e.inchi_string().ok()) {
                    Some(a) if a == recorded => {}
                    Some(a) => problems.push(format!("fixedh inchi {a:?} != {recorded:?}")),
                    None => problems.push("fixedh inchi engine err".to_string()),
                }
                if let Some(recorded_aux) = case["fixedh"]
                    .as_object()
                    .and_then(|s| s.get("aux"))
                    .and_then(|v| v.as_str())
                {
                    match engine.as_ref().ok().and_then(|e| e.aux_string().ok()) {
                        Some(a) if a == recorded_aux => {}
                        Some(a) => problems.push(format!("fixedh aux {a:?} != {recorded_aux:?}")),
                        None => problems.push("fixedh aux engine err".to_string()),
                    }
                }
            }

            // 5. RecMet identifier + aux.
            if let Some(recmet) = case["recmet"].as_object() {
                let assembled = chemomer::Engine::with_options(
                    molfile,
                    chemomer::Options::default().with_rec_met(),
                );
                if let Some(recorded) = recmet.get("inchi").and_then(|v| v.as_str()) {
                    match assembled.as_ref().ok().map(|e| e.inchi_string()) {
                        Some(Ok(a)) if a == recorded => {}
                        Some(Ok(a)) => problems.push(format!("recmet inchi {a:?} != {recorded:?}")),
                        Some(Err(e)) => problems.push(format!("recmet inchi err {e}")),
                        None => problems.push("recmet engine absent".to_string()),
                    }
                }
                if let Some(recorded) = recmet.get("aux").and_then(|v| v.as_str()) {
                    match assembled.as_ref().ok().map(|e| e.aux_string()) {
                        Some(Ok(a)) if a == recorded => {}
                        Some(Ok(a)) => problems.push(format!("recmet aux {a:?} != {recorded:?}")),
                        Some(Err(e)) => problems.push(format!("recmet aux err {e}")),
                        None => problems.push("recmet engine absent".to_string()),
                    }
                }
            }

            // 6. DoNotAddH.
            if let Some(recorded) = case["donotaddh"]
                .as_object()
                .and_then(|s| s.get("inchi"))
                .and_then(|v| v.as_str())
            {
                match chemomer::Engine::with_options(
                    molfile,
                    chemomer::Options::default().with_do_not_add_h(),
                )
                .and_then(|e| e.inchi_string())
                {
                    Ok(a) if a == recorded => {}
                    Ok(a) => problems.push(format!("donotaddh {a:?} != {recorded:?}")),
                    Err(e) => problems.push(format!("donotaddh err {e}")),
                }
            }

            checked += 1;
            if problems.is_empty() {
                clean.push(id.to_string());
            } else {
                dirty += 1;
                println!("{id}: DIRTY ({} recording problem(s))", problems.len());
                for problem in &problems {
                    println!("  - {problem}");
                }
            }
        }
    }
    println!(
        "\ncandidates checked: {checked}, clean: {}, dirty: {dirty}",
        clean.len()
    );
    for id in &clean {
        println!("  CLEAN {id}");
    }
}
