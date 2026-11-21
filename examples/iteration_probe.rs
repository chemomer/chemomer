//! Print the write∘read iteration sequence for the two remaining
//! unstable rows (breadth-534, breadth-1152).
use std::fs;
use std::path::PathBuf;

fn main() {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    let mut cases: Vec<serde_json::Value> = Vec::new();
    let text =
        fs::read_to_string(dir.join("corpus/fixtures/identity-breadth.json"))
            .expect("identity-breadth");
    cases.extend(serde_json::from_str::<Vec<serde_json::Value>>(&text).expect("fixtures"));
    for case in &cases {
        let id = case["id"].as_str().unwrap_or("");
        if id != "identity/breadth-534" && id != "identity/breadth-1152" {
            continue;
        }
        let smiles = case["smiles"].as_str().unwrap_or("");
        let structure = chemomer::parse_smiles(smiles).expect(id);
        println!("{id} input: {smiles}");
        let mut current = chemomer::write_smiles(&structure).expect("write");
        println!("  w0: {current}");
        for round in 1..25 {
            let re = chemomer::parse_smiles(&current).expect("re-read");
            let next = chemomer::write_smiles(&re).expect("re-write");
            if next == current {
                println!("  round {round}: FIXED POINT");
                break;
            }
            println!("  round {round}: {next}");
            current = next;
        }
    }
}
