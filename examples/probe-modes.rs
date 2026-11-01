// Probe: every recorded mode for the ex-triaged fixtures.
use chemomer::{Engine, Options};

fn cmp(label: &str, id: &str, ours: chemomer::Result<String>, recorded: Option<&str>) -> bool {
    match (ours, recorded) {
        (Ok(a), Some(r)) => {
            if a == r {
                true
            } else {
                println!("  FAIL {label} {id}:");
                println!("    ours:   {a}");
                println!("    oracle: {r}");
                false
            }
        }
        (Ok(a), None) => {
            println!("  FAIL {label} {id}: no recording, ours {a}");
            false
        }
        (Err(e), Some(r)) => {
            println!("  FAIL {label} {id}: ERR {} oracle {r}", e.message());
            false
        }
        (Err(e), None) => {
            println!("  FAIL {label} {id}: ERR {}", e.message());
            false
        }
    }
}

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let ids: Vec<String> = std::fs::read_to_string("/tmp/triaged-ids.txt")
        .unwrap()
        .lines()
        .map(String::from)
        .collect();
    let mut ok = 0;
    let mut bad = 0;
    for name in [
        "identity-metal.json",
        "identity-multicomponent.json",
        "identity.json",
    ] {
        let text = std::fs::read_to_string(format!("{dir}/corpus/fixtures/{name}")).unwrap();
        let cases: Vec<serde_json::Value> = serde_json::from_str(&text).unwrap();
        for case in cases {
            let id = case["id"].as_str().unwrap_or("").to_string();
            if !ids.contains(&id) {
                continue;
            }
            let molfile = case["molfile"].as_str().unwrap_or("");
            let mut all = true;
            all &= cmp(
                "std",
                &id,
                Engine::from_molfile(molfile).and_then(|e| e.inchi_string()),
                case.get("inchi").and_then(|v| v.as_str()),
            );
            all &= cmp(
                "std-aux",
                &id,
                Engine::from_molfile(molfile).and_then(|e| e.aux_string()),
                case.get("aux").and_then(|v| v.as_str()),
            );
            if let Some(rec) = case.get("donotaddh") {
                all &= cmp(
                    "dna",
                    &id,
                    Engine::with_options(molfile, Options::default().with_do_not_add_h())
                        .and_then(|e| e.inchi_string()),
                    rec.get("inchi").and_then(|v| v.as_str()),
                );
                all &= cmp(
                    "dna-aux",
                    &id,
                    Engine::with_options(molfile, Options::default().with_do_not_add_h())
                        .and_then(|e| e.aux_string()),
                    rec.get("aux").and_then(|v| v.as_str()),
                );
            }
            if let Some(rec) = case.get("fixedh") {
                all &= cmp(
                    "fixedh",
                    &id,
                    Engine::with_options(molfile, Options::default().with_fixed_h())
                        .and_then(|e| e.inchi_string()),
                    rec.get("inchi").and_then(|v| v.as_str()),
                );
                all &= cmp(
                    "fixedh-aux",
                    &id,
                    Engine::with_options(molfile, Options::default().with_fixed_h())
                        .and_then(|e| e.aux_string()),
                    rec.get("aux").and_then(|v| v.as_str()),
                );
            }
            if let Some(rec) = case.get("snon") {
                all &= cmp(
                    "snon",
                    &id,
                    Engine::with_options(molfile, Options::default().with_no_stereo())
                        .and_then(|e| e.inchi_string()),
                    rec.get("inchi").and_then(|v| v.as_str()),
                );
                all &= cmp(
                    "snon-aux",
                    &id,
                    Engine::with_options(molfile, Options::default().with_no_stereo())
                        .and_then(|e| e.aux_string()),
                    rec.get("aux").and_then(|v| v.as_str()),
                );
            }
            if let Some(rec) = case.get("recmet") {
                all &= cmp(
                    "recmet",
                    &id,
                    Engine::with_options(molfile, Options::default().with_rec_met())
                        .and_then(|e| e.recmet_string()),
                    rec.get("inchi").and_then(|v| v.as_str()),
                );
                all &= cmp(
                    "recmet-aux",
                    &id,
                    Engine::with_options(molfile, Options::default().with_rec_met())
                        .and_then(|e| e.aux_string()),
                    rec.get("aux").and_then(|v| v.as_str()),
                );
            }
            if all {
                ok += 1;
            } else {
                bad += 1;
            }
        }
    }
    println!("fixtures all-modes-ok {ok}, with-failures {bad}");
}
