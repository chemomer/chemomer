//! The four rows that looped write∘read before the canonical
//! renumbering — stability now rides the fixed-point writer.
use std::fs;
use std::path::PathBuf;

fn main() {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    let mut cases: Vec<serde_json::Value> = Vec::new();
    for name in [
        "corpus/fixtures/identity.json",
        "corpus/fixtures/identity-breadth.json",
        "corpus/fixtures/identity-metal.json",
        "corpus/fixtures/identity-multicomponent.json",
    ] {
        let text = fs::read_to_string(dir.join(name)).expect(name);
        cases.extend(serde_json::from_str::<Vec<serde_json::Value>>(&text).expect("fixtures"));
    }
    const LOOPS: [&str; 4] = [
        "identity/breadth-180",
        "identity/breadth-251",
        "identity/breadth-1327",
        "identity/breadth-1351",
    ];
    for case in &cases {
        let id = case["id"].as_str().unwrap_or("");
        if !LOOPS.contains(&id) {
            continue;
        }
        let Some(smiles) = case["smiles"].as_str() else {
            continue;
        };
        let structure = chemomer::parse_smiles(smiles).expect(id);
        let Ok(w1) = chemomer::write_smiles(&structure) else {
            println!("{id}: WRITE-DECLINED");
            continue;
        };
        let stable = chemomer::parse_smiles(&w1)
            .ok()
            .and_then(|re| chemomer::write_smiles(&re).ok())
            .is_some_and(|again| again == w1);
        let identity = chemomer::Engine::from_smiles(&w1)
            .and_then(|e| e.inchi_string())
            .is_ok_and(|derived| derived == case["inchi"].as_str().unwrap_or(""));
        println!("{id}: stable={stable} identity={identity} -> {w1}");
    }
}
