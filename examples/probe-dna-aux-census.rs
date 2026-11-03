// Census: every fixture's DoNotAddH aux recording vs the engine,
// classified by the first diverging aux field (N: E: CRV: rA rB rC).
use chemomer::{Engine, Options};

fn first_field_divergence(ours: &str, oracle: &str) -> String {
    let split = |text: &str| -> Vec<(String, String)> {
        let body = text.strip_prefix("AuxInfo=1").unwrap_or(text);
        let mut fields = Vec::new();
        for raw in body.split('/') {
            if raw.is_empty() {
                continue;
            }
            let (tag, rest) = match raw.split_once(':') {
                Some((tag, rest)) => (tag.to_string(), rest.to_string()),
                None => (String::new(), raw.to_string()),
            };
            fields.push((tag, rest));
        }
        fields
    };
    let (a, b) = (split(ours), split(oracle));
    for (index, (ta, va)) in a.iter().enumerate() {
        match b.get(index) {
            Some((tb, vb)) if tb == ta && vb == va => {}
            Some((tb, _)) if tb != ta => return format!("tag-mismatch {ta} vs {tb}"),
            _ => return format!("/{ta}"),
        }
    }
    if a.len() != b.len() {
        return format!("field-count {} vs {}", a.len(), b.len());
    }
    "none".to_string()
}

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let fixtures = std::fs::read_dir(format!("{dir}/corpus/fixtures")).unwrap();
    let mut files: Vec<_> = fixtures
        .map(|e| e.unwrap().path())
        .filter(|p| {
            let name = p.file_name().unwrap().to_string_lossy().to_string();
            name.starts_with("identity") && name.ends_with(".json")
        })
        .collect();
    files.sort();
    let mut ledger = std::collections::HashSet::new();
    let text = std::fs::read_to_string(format!("{dir}/corpus/donotaddh-aux-triage.json")).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&text).unwrap();
    if let Some(map) = parsed.as_object() {
        for key in map.keys() {
            ledger.insert(key.clone());
        }
    }
    let mut pass = 0u32;
    let mut red = 0u32;
    let mut classes: std::collections::BTreeMap<String, Vec<String>> = Default::default();
    let mut ledgered_red = 0u32;
    let mut ledgered_pass: Vec<String> = Vec::new();
    let mut unledgered_red: Vec<String> = Vec::new();
    for file in &files {
        let text = std::fs::read_to_string(file).unwrap();
        let cases: Vec<serde_json::Value> = serde_json::from_str(&text).unwrap();
        for case in &cases {
            let id = case["id"].as_str().unwrap_or("").to_string();
            let Some(recorded) = case
                .get("donotaddh")
                .and_then(|d| d.get("aux"))
                .and_then(|v| v.as_str())
            else {
                continue;
            };
            let molfile = case["molfile"].as_str().unwrap_or("");
            let ours = Engine::with_options(molfile, Options::default().with_do_not_add_h())
                .and_then(|e| e.aux_string());
            match ours {
                Ok(ours) if ours == recorded => {
                    pass += 1;
                    if ledger.contains(&id) {
                        ledgered_pass.push(id);
                    }
                }
                Ok(ours) => {
                    red += 1;
                    let class = first_field_divergence(&ours, recorded);
                    classes.entry(class).or_default().push(id.clone());
                    if ledger.contains(&id) {
                        ledgered_red += 1;
                    } else {
                        unledgered_red.push(id);
                    }
                }
                Err(error) => {
                    red += 1;
                    classes
                        .entry(format!("ERR {}", error.message()))
                        .or_default()
                        .push(id.clone());
                    if ledger.contains(&id) {
                        ledgered_red += 1;
                    } else {
                        unledgered_red.push(id);
                    }
                }
            }
        }
    }
    println!("dna-aux census: pass {pass}, red {red} (ledgered red {ledgered_red}, ledgered-now-pass {})", ledgered_pass.len());
    if !ledgered_pass.is_empty() {
        println!("RETIRED (ledgered, now passing):");
        for id in &ledgered_pass {
            println!("  {id}");
        }
    }
    println!("divergence classes:");
    for (class, ids) in &classes {
        println!(
            "  {} [{}]: {}...",
            class,
            ids.len(),
            ids.first().map(String::as_str).unwrap_or("")
        );
    }
    println!("unledgered red ({}):", unledgered_red.len());
    for id in unledgered_red.iter().take(30) {
        println!("  {id}");
    }
    if std::env::var("DUMP").is_ok() {
        let want: std::collections::HashSet<&str> =
            unledgered_red.iter().map(String::as_str).collect();
        for file in &files {
            let text = std::fs::read_to_string(file).unwrap();
            let cases: Vec<serde_json::Value> = serde_json::from_str(&text).unwrap();
            for case in &cases {
                let id = case["id"].as_str().unwrap_or("");
                if !want.contains(id) {
                    continue;
                }
                let Some(recorded) = case
                    .get("donotaddh")
                    .and_then(|d| d.get("aux"))
                    .and_then(|v| v.as_str())
                else {
                    continue;
                };
                let molfile = case["molfile"].as_str().unwrap_or("");
                if let Ok(ours) =
                    Engine::with_options(molfile, Options::default().with_do_not_add_h())
                        .and_then(|e| e.aux_string())
                {
                    if ours != recorded {
                        println!("=== {id}\n  ours:   {ours}\n  oracle: {recorded}");
                    }
                }
            }
        }
    }
}
