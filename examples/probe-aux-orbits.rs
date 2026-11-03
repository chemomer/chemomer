//! Aux /E: census: for every aux-triaged fixture, compare the
//! assembled auxiliary line with the oracle recording — the /E:
//! retirement radar.
use std::fs;
use std::path::PathBuf;

use chemomer::Engine;

fn main() {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    let aux_triage: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(dir.join("corpus/aux-triage.json")).expect("aux triage"),
    )
    .expect("valid");
    let mut files: Vec<PathBuf> = fs::read_dir(dir.join("corpus/fixtures"))
        .expect("dir")
        .filter_map(|e| {
            let p = e.expect("entry").path();
            let n = p.file_name()?.to_string_lossy().to_string();
            (n.starts_with("identity") && n.ends_with(".json")).then_some(p)
        })
        .collect();
    files.sort();
    let mut passing = Vec::new();
    let mut failing: Vec<(String, String, String)> = Vec::new();
    let mut total = 0usize;
    for file in files {
        let parsed: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read")).expect("json");
        for case in parsed {
            let id = case["id"].as_str().expect("id").to_string();
            if aux_triage.get(&id).is_none() {
                continue;
            }
            let Some(recorded) = case["aux"].as_str() else {
                continue;
            };
            total += 1;
            let Ok(engine) = Engine::from_molfile(case["molfile"].as_str().unwrap()) else {
                failing.push((id, "ENGINE ERROR".into(), recorded.into()));
                continue;
            };
            let Ok(assembled) = engine.aux_string() else {
                failing.push((id, "AUX ERROR".into(), recorded.into()));
                continue;
            };
            if assembled == recorded {
                passing.push(id);
            } else {
                failing.push((id, assembled, recorded.into()));
            }
        }
    }
    println!("aux-triaged {total}, now passing {}:", passing.len());
    for id in &passing {
        println!("  {id}");
    }
    println!("aux-triaged, still failing {}:", failing.len());
    for (id, assembled, oracle) in &failing {
        println!("  {id}\n    ours   {assembled}\n    oracle {oracle}");
        // The /E: fields side by side when both carry one.
        let ours_e = assembled.split("/E:").nth(1).map(|s| {
            s.split('/').next().unwrap_or("").to_string()
        });
        let oracle_e = oracle.split("/E:").nth(1).map(|s| {
            s.split('/').next().unwrap_or("").to_string()
        });
        if let (Some(o), Some(r)) = (ours_e, oracle_e) {
            println!("    /E: ours {o}\n    /E: oracle   {r}");
        }
    }
}
