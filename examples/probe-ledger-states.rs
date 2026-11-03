// Per-ledger-fixture gate states: for every ledgered fixture, which
// gates pass and which fail (aux, dna identifier, dna aux).
use chemomer::{Engine, Options};
use std::collections::BTreeMap;

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let mut fixtures = std::collections::BTreeMap::new();
    let dpath = format!("{dir}/corpus/fixtures");
    for entry in std::fs::read_dir(&dpath).unwrap() {
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
    let ledgers = [
        ("aux", format!("{dir}/corpus/aux-triage.json")),
        ("dna-id", format!("{dir}/corpus/donotaddh-triage.json")),
        ("dna-aux", format!("{dir}/corpus/donotaddh-aux-triage.json")),
    ];
    let mut members: BTreeMap<&str, Vec<String>> = Default::default();
    for (name, path) in &ledgers {
        let parsed: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        members.insert(name, parsed.as_object().unwrap().keys().cloned().collect());
    }
    let all: std::collections::BTreeSet<String> =
        members.values().flat_map(|v| v.iter().cloned()).collect();
    println!("union: {}", all.len());
    for id in &all {
        let c = &fixtures[id];
        let molfile = c["molfile"].as_str().unwrap();
        let tags: Vec<String> = ledgers
            .iter()
            .filter(|(n, _)| members[n].iter().any(|m| m == id))
            .map(|(n, _)| n.to_string())
            .collect();
        let mut states = Vec::new();
        // aux gate
        if let Some(rec) = c.get("aux").and_then(|v| v.as_str()) {
            let ours = Engine::from_molfile(molfile).and_then(|e| e.aux_string());
            states.push(match ours {
                Ok(o) if o == rec => "aux:PASS".into(),
                Ok(_) => "aux:FAIL".into(),
                Err(e) => format!("aux:ERR({e})"),
            });
        }
        // dna identifier gate
        if let Some(rec) = c
            .get("donotaddh")
            .and_then(|d| d.get("inchi"))
            .and_then(|v| v.as_str())
        {
            let ours = Engine::with_options(molfile, Options::default().with_do_not_add_h())
                .and_then(|e| e.inchi_string());
            states.push(match ours {
                Ok(o) if o == rec => "dna-id:PASS".into(),
                Ok(o) => format!("dna-id:FAIL({o} vs {rec})"),
                Err(e) => format!("dna-id:ERR({e})"),
            });
        }
        // dna aux gate
        if let Some(rec) = c
            .get("donotaddh")
            .and_then(|d| d.get("aux"))
            .and_then(|v| v.as_str())
        {
            let ours = Engine::with_options(molfile, Options::default().with_do_not_add_h())
                .and_then(|e| e.aux_string());
            states.push(match ours {
                Ok(o) if o == rec => "dna-aux:PASS".into(),
                Ok(o) => format!("dna-aux:FAIL({}..)", &o[..o.len().min(60)]),
                Err(e) => format!("dna-aux:ERR({e})"),
            });
        }
        println!("{id} [{}] {}", tags.join(","), states.join(" | "));
    }
}
