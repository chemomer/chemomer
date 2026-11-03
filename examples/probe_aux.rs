//! Development probe (NOT committed): dump the aux diff for fixtures
//! named on the command line, or the full pass/skip census when none.

use std::fs;

fn suite_dir() -> std::path::PathBuf {
    std::env::var("INCHI_TESTS")
        .expect("set INCHI_TESTS to the inchi-tests worktree")
        .into()
}

fn load_cases(dir: &std::path::Path) -> Vec<serde_json::Value> {
    let mut cases = Vec::new();
    let entries = [
        "identity.json",
        "identity-breadth.json",
        "identity-breadth2.json",
        "identity-breadth3.json",
        "identity-ez.json",
        "identity-metal.json",
        "identity-multicomponent.json",
        "identity-stereo.json",
        "large-molecules.json",
        "tautomer-extremes.json",
        "v3000-molfiles.json",
        "sdf-batch.json",
    ];
    for name in entries {
        let path = dir.join("corpus/fixtures").join(name);
        let Ok(text) = fs::read_to_string(&path) else {
            continue;
        };
        let value: serde_json::Value = serde_json::from_str(&text).expect("fixture json");
        let Some(array) = value.as_array() else {
            continue;
        };
        cases.extend(array.iter().cloned());
    }
    cases
}

fn load_ledger(dir: &std::path::Path, name: &str) -> std::collections::HashSet<String> {
    let path = dir.join("corpus").join(name);
    let Ok(text) = fs::read_to_string(&path) else {
        return Default::default();
    };
    let value: serde_json::Value = serde_json::from_str(&text).expect("ledger json");
    value
        .as_object()
        .map(|map| map.keys().cloned().collect())
        .unwrap_or_default()
}

fn main() {
    let dir = suite_dir();
    let cases = load_cases(&dir);
    // The four ledgers the dna-aux gate unions (triage, aux, dna, dna-aux).
    let ledgers: std::collections::HashSet<String> = [
        "triage.json",
        "aux-triage.json",
        "donotaddh-triage.json",
        "donotaddh-aux-triage.json",
    ]
    .iter()
    .flat_map(|name| load_ledger(&dir, name))
    .collect();
    // `dna` selects the DoNotAddH mode (the 211-family's gate).
    let dna = std::env::args()
        .nth(1)
        .map(|arg| arg == "dna")
        .unwrap_or(false);
    let targets: Vec<String> = if dna {
        std::env::args().skip(2).collect()
    } else {
        std::env::args().skip(1).collect()
    };
    let record = |case: &serde_json::Value| -> Option<String> {
        if dna {
            case["donotaddh"]
                .as_object()
                .and_then(|s| s.get("aux"))
                .and_then(|v| v.as_str())
                .map(|s| s.to_string())
        } else {
            case["aux"].as_str().map(|s| s.to_string())
        }
    };
    let assemble = |molfile: &str| -> chemomer::Result<String> {
        let options = if dna {
            chemomer::Options::default().with_do_not_add_h()
        } else {
            chemomer::Options::default()
        };
        chemomer::Engine::with_options(molfile, options).and_then(|engine| engine.aux_string())
    };

    let mut pass = 0usize;
    let mut skip = 0usize;
    let mut unexpected = 0usize;
    let mut skip_ids: Vec<String> = Vec::new();
    for case in &cases {
        let id = case["id"].as_str().expect("id").to_string();
        let Some(recorded) = record(case) else {
            continue;
        };
        let ledgered = ledgers.contains(&id);
        let molfile = case["molfile"].as_str().expect("molfile");
        let assembled = assemble(molfile);
        let ok = matches!(&assembled, Ok(text) if *text == recorded);
        if !targets.is_empty() {
            if targets.iter().any(|t| id.contains(t.as_str())) {
                println!("=== {id}  ledgered={ledgered}");
                println!("oracle:   {recorded}");
                match &assembled {
                    Ok(text) => println!("engine:   {text}"),
                    Err(error) => println!("engine:   ERR {error}"),
                }
                // Field-by-field diff.
                if let Ok(text) = &assembled {
                    let a: Vec<&str> = recorded.split('/').collect();
                    let b: Vec<&str> = text.split('/').collect();
                    for i in 0..a.len().max(b.len()) {
                        let (x, y) = (a.get(i), b.get(i));
                        if x != y {
                            println!("  field {i}: oracle {x:?} != engine {y:?}");
                        }
                    }
                }
            }
            continue;
        }
        if ok {
            pass += 1;
        } else if ledgered {
            skip += 1;
            println!("LEDGERED-FAIL {id}");
            skip_ids.push(id);
            if let Ok(text) = &assembled {
                let a: Vec<&str> = recorded.split('/').collect();
                let b: Vec<&str> = text.split('/').collect();
                for i in 0..a.len().max(b.len()) {
                    let (x, y) = (a.get(i), b.get(i));
                    if x != y {
                        let kind = match i {
                            0 => "header",
                            1 => "N",
                            2 => "E",
                            3 => "gE",
                            4 => "CRV",
                            5 => "rA",
                            _ => "tail",
                        };
                        println!("  {kind}: {x:?} != {y:?}");
                    }
                }
            }
        } else {
            unexpected += 1;
            println!("UNEXPECTED-FAIL {id}");
            match &assembled {
                Ok(text) => println!("  oracle: {recorded}\n  engine: {text}"),
                Err(error) => println!("  error:  {error}"),
            }
        }
    }
    println!("pass {pass}, ledgered-fail {skip}, unexpected-fail {unexpected}");
    println!("skip set: {}", skip_ids.join(" "));
}
