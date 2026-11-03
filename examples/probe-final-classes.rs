// Probe: ours-vs-oracle for the Class 3 (dna-aux) and Class 4 (dna
// identifier) ledger fixtures.
use chemomer::{Engine, Options};

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let class = std::env::var("CLASS").unwrap_or_default();
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
    let (ledger_path, want_aux, want_id) = if class == "aux" {
        (format!("{dir}/corpus/donotaddh-aux-triage.json"), true, false)
    } else if class == "id" {
        (format!("{dir}/corpus/donotaddh-triage.json"), true, true)
    } else {
        panic!("set CLASS=aux or CLASS=id");
    };
    let ledger: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&ledger_path).unwrap()).unwrap();
    for id in ledger.as_object().unwrap().keys() {
        if !filter.is_empty() && !id.contains(&filter) {
            continue;
        }
        let c = &fixtures[id];
        let molfile = c["molfile"].as_str().unwrap();
        let engine =
            Engine::with_options(molfile, Options::default().with_do_not_add_h()).unwrap();
        if want_id {
            let recorded = c["donotaddh"]["inchi"].as_str().unwrap_or("");
            let ours = engine.inchi_string().unwrap();
            println!("=== {id} identifier");
            println!("  ours:   {ours}");
            println!("  oracle: {recorded}");
            println!("  match:  {}", ours == recorded);
        }
        if want_aux {
            let recorded = c["donotaddh"]["aux"].as_str().unwrap_or("");
            let ours = engine.aux_string().unwrap();
            println!("=== {id} aux");
            println!("  ours:   {ours}");
            println!("  oracle: {recorded}");
            println!("  match:  {}", ours == recorded);
        }
    }
}
