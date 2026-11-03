// Full ours-vs-oracle dump for every donotaddh-aux-ledger fixture.
use chemomer::{Engine, Options};

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let filter = std::env::var("ONLY").unwrap_or_default();
    let mut fixtures = std::collections::BTreeMap::new();
    for entry in std::fs::read_dir(format!("{dir}/corpus/fixtures")).unwrap() {
        let p = entry.unwrap().path();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        if !(name.starts_with("identity") && name.ends_with(".json")) {
            continue;
        }
        let cases: Vec<serde_json::Value> =
            serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
        for c in cases {
            fixtures.insert(c["id"].as_str().unwrap().to_string(), c);
        }
    }
    let ledger: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(format!("{dir}/corpus/donotaddh-aux-triage.json")).unwrap(),
    )
    .unwrap();
    for id in ledger.as_object().unwrap().keys() {
        if !filter.is_empty() && !id.contains(&filter) {
            continue;
        }
        let c = &fixtures[id];
        let molfile = c["molfile"].as_str().unwrap();
        let recorded = c["donotaddh"]["aux"].as_str().unwrap_or("");
        match Engine::with_options(molfile, Options::default().with_do_not_add_h())
            .and_then(|e| e.aux_string())
        {
            Ok(ours) if ours != recorded => {
                println!("=== {id}");
                println!("  ours:   {ours}");
                println!("  oracle: {recorded}");
            }
            _ => {}
        }
    }
}
