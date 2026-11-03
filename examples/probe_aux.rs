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

fn main() {
    let dir = suite_dir();
    let cases = load_cases(&dir);
    let ledger_text =
        fs::read_to_string(dir.join("corpus/aux-triage.json")).expect("aux ledger");
    let ledger: serde_json::Value = serde_json::from_str(&ledger_text).expect("ledger json");
    let targets: Vec<String> = std::env::args().skip(1).collect();

    let mut pass = 0usize;
    let mut skip = 0usize;
    let mut unexpected = 0usize;
    for case in &cases {
        let id = case["id"].as_str().expect("id").to_string();
        let Some(recorded) = case["aux"].as_str() else {
            continue;
        };
        let ledgered = ledger.get(&id).is_some();
        let molfile = case["molfile"].as_str().expect("molfile");
        let assembled = chemomer::Engine::from_molfile(molfile)
            .and_then(|engine| engine.aux_string());
        let ok = matches!(&assembled, Ok(text) if text == recorded);
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
}
