// Probe: the citrate deprotonated-hydroxyl witness fixture
// (identity/multi-sodium-citrate) across EVERY recorded layer,
// including the ungated donotaddh.aux / snon.aux recordings.
use chemomer::{Engine, Options};

fn cmp(label: &str, ours: chemomer::Result<String>, recorded: Option<&str>) -> bool {
    match (ours, recorded) {
        (Ok(a), Some(r)) => {
            if a == r {
                println!("  ok   {label}");
                true
            } else {
                println!("  FAIL {label}:");
                println!("    ours:   {a}");
                println!("    oracle: {r}");
                false
            }
        }
        (Ok(a), None) => {
            println!("  FAIL {label}: no recording, ours {a}");
            false
        }
        (Err(e), Some(r)) => {
            println!("  FAIL {label}: ERR {} oracle {r}", e.message());
            false
        }
        (Err(e), None) => {
            println!("  FAIL {label}: ERR {}", e.message());
            false
        }
    }
}

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let text = std::fs::read_to_string(format!(
        "{dir}/corpus/fixtures/identity-multicomponent.json"
    ))
    .unwrap();
    let cases: Vec<serde_json::Value> = serde_json::from_str(&text).unwrap();
    let mut all_ok = true;
    for case in &cases {
        let id = case["id"].as_str().unwrap_or("");
        if !id.contains("citrate") && !id.contains("citric") {
            continue;
        }
        println!("fixture {id}");
        let molfile = case["molfile"].as_str().unwrap_or("");
        all_ok &= cmp(
            "std",
            Engine::from_molfile(molfile).and_then(|e| e.inchi_string()),
            case.get("inchi").and_then(|v| v.as_str()),
        );
        all_ok &= cmp(
            "std-aux",
            Engine::from_molfile(molfile).and_then(|e| e.aux_string()),
            case.get("aux").and_then(|v| v.as_str()),
        );
        if let Some(rec) = case.get("snon") {
            all_ok &= cmp(
                "snon",
                Engine::with_options(molfile, Options::default().with_no_stereo())
                    .and_then(|e| e.inchi_string()),
                rec.get("inchi").and_then(|v| v.as_str()),
            );
            all_ok &= cmp(
                "snon-aux (ungated)",
                Engine::with_options(molfile, Options::default().with_no_stereo())
                    .and_then(|e| e.aux_string()),
                rec.get("aux").and_then(|v| v.as_str()),
            );
        }
        if let Some(rec) = case.get("donotaddh") {
            all_ok &= cmp(
                "dna",
                Engine::with_options(molfile, Options::default().with_do_not_add_h())
                    .and_then(|e| e.inchi_string()),
                rec.get("inchi").and_then(|v| v.as_str()),
            );
            all_ok &= cmp(
                "dna-aux (ungated)",
                Engine::with_options(molfile, Options::default().with_do_not_add_h())
                    .and_then(|e| e.aux_string()),
                rec.get("aux").and_then(|v| v.as_str()),
            );
        }
        if let Some(rec) = case.get("fixedh") {
            all_ok &= cmp(
                "fixedh",
                Engine::with_options(molfile, Options::default().with_fixed_h())
                    .and_then(|e| e.inchi_string()),
                rec.get("inchi").and_then(|v| v.as_str()),
            );
            all_ok &= cmp(
                "fixedh-aux",
                Engine::with_options(molfile, Options::default().with_fixed_h())
                    .and_then(|e| e.aux_string()),
                rec.get("aux").and_then(|v| v.as_str()),
            );
        }
        if let Some(rec) = case.get("recmet") {
            all_ok &= cmp(
                "recmet",
                Engine::with_options(molfile, Options::default().with_rec_met())
                    .and_then(|e| e.recmet_string()),
                rec.get("inchi").and_then(|v| v.as_str()),
            );
            all_ok &= cmp(
                "recmet-aux",
                Engine::with_options(molfile, Options::default().with_rec_met())
                    .and_then(|e| e.aux_string()),
                rec.get("aux").and_then(|v| v.as_str()),
            );
        }
    }
    println!("ALL OK: {all_ok}");
}
