//! The rows whose aromatic-N spelling moved identity on the oracle
//! gate: with the donor-substituent exclusion they must spell
//! kekulized again and match byte-for-byte.
use std::fs;
use std::path::PathBuf;

fn main() {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    let mut cases: Vec<serde_json::Value> = Vec::new();
    let text = fs::read_to_string(dir.join("corpus/fixtures/identity-breadth.json"))
        .expect("identity-breadth");
    cases.extend(serde_json::from_str::<Vec<serde_json::Value>>(&text).expect("fixtures"));
    const DRIFT: [&str; 11] = [
        "identity/breadth-31",
        "identity/breadth-151",
        "identity/breadth-189",
        "identity/breadth-203",
        "identity/breadth-316",
        "identity/breadth-406",
        "identity/breadth-452",
        "identity/breadth-727",
        "identity/breadth-744",
        "identity/breadth-773",
        "identity/breadth-952",
    ];
    let mut exact = 0;
    let mut total = 0;
    for case in &cases {
        let id = case["id"].as_str().unwrap_or("");
        if !DRIFT.contains(&id) {
            continue;
        }
        total += 1;
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
        let aromatic = w1.contains(['c', 'n']) && w1.chars().any(|ch| ch.is_lowercase());
        let exact_now = chemomer::Engine::from_smiles(&w1)
            .and_then(|e| e.inchi_string())
            .is_ok_and(|d| d == case["inchi"].as_str().unwrap_or(""));
        if exact_now {
            exact += 1;
        }
        println!(
            "{id}: stable={stable} identity={} aromatic={aromatic} w={w1}",
            if exact_now { "EXACT" } else { "DRIFT" }
        );
    }
    println!("exact: {exact}/{total}");
}
