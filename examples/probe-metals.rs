// Probe: assembled vs recorded for the still-triaged metal fixtures.
use chemomer::{Engine, Options};

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let want: Vec<String> = std::env::args().skip(1).collect();
    for name in [
        "identity-metal.json",
        "identity-multicomponent.json",
        "identity.json",
    ] {
        let text = std::fs::read_to_string(format!("{dir}/corpus/fixtures/{name}")).unwrap();
        let cases: Vec<serde_json::Value> = serde_json::from_str(&text).unwrap();
        for case in cases {
            let id = case["id"].as_str().unwrap_or("");
            if !want.iter().any(|w| id.contains(w)) {
                continue;
            }
            let molfile = case["molfile"].as_str().unwrap_or("");
            println!("== {id}");
            match Engine::from_molfile(molfile).and_then(|e| e.inchi_string()) {
                Ok(a) => {
                    let recorded = case["inchi"].as_str().unwrap_or("");
                    if a == recorded {
                        println!("  PASS main");
                    } else {
                        println!("  FAIL main:");
                        println!("    ours:   {a}");
                        println!("    oracle: {recorded}");
                    }
                }
                Err(e) => println!("  ERR main: {}", e.message()),
            }
            if let Some(recmet) = case.get("recmet") {
                let recorded = recmet.get("inchi").and_then(|v| v.as_str()).unwrap_or("");
                match Engine::with_options(molfile, Options::default().with_rec_met())
                    .and_then(|e| e.recmet_string())
                {
                    Ok(a) => {
                        if a == recorded {
                            println!("  PASS recmet");
                        } else {
                            println!("  FAIL recmet:");
                            println!("    ours:   {a}");
                            println!("    oracle: {recorded}");
                        }
                    }
                    Err(e) => println!("  ERR recmet: {}", e.message()),
                }
            }
        }
    }
}
