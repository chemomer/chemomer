//! Probe the 13 KEKULE_RESIDUAL rows against the renumbered engine:
//! stability (write∘read∘write) and identity through the recorded oracle.
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
    const RESIDUAL: [&str; 13] = [
        "identity/breadth-237",
        "identity/breadth-273",
        "identity/breadth-364",
        "identity/breadth-371",
        "identity/breadth-403",
        "identity/breadth-453",
        "identity/breadth-534",
        "identity/breadth-661",
        "identity/breadth-685",
        "identity/breadth-862",
        "identity/breadth-885",
        "identity/breadth-1152",
        "identity/breadth-1429",
    ];
    let mut closed = 0;
    for case in &cases {
        let id = case["id"].as_str().unwrap_or("");
        if !RESIDUAL.contains(&id) {
            continue;
        }
        let Some(smiles) = case["smiles"].as_str() else {
            continue;
        };
        let recorded = case["inchi"].as_str().unwrap_or("");
        let mut line = format!("{id}: ");
        let Ok(structure) = chemomer::parse_smiles(smiles) else {
            line.push_str("PARSE-FAIL");
            println!("{line}");
            continue;
        };
        let Ok(w1) = chemomer::write_smiles(&structure) else {
            line.push_str("WRITE-DECLINED");
            println!("{line}");
            continue;
        };
        match chemomer::parse_smiles(&w1) {
            Err(e) => line.push_str(&format!("REPARSE-FAIL ({e})")),
            Ok(re) => {
                let w2 = chemomer::write_smiles(&re).unwrap_or_default();
                if w2 == w1 {
                    line.push_str("stable ");
                } else {
                    line.push_str("UNSTABLE ");
                }
                match chemomer::Engine::from_smiles(&w1).and_then(|e| e.inchi_string()) {
                    Ok(derived) if derived == recorded => {
                        line.push_str("identity-EXACT");
                        closed += 1;
                    }
                    Ok(derived) => line.push_str(&format!("identity-DRIFT ({derived})")),
                    Err(e) => line.push_str(&format!("ENGINE-FAIL ({e})")),
                }
            }
        }
        println!("{line}");
    }
    println!("closed: {closed}/13");
}
