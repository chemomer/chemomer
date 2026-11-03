use chemomer::{Engine, Options};
fn main() {
    census();
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
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
    let c = &fixtures["identity/breadth3-1053"];
    let engine = Engine::with_options(
        c["molfile"].as_str().unwrap(),
        Options::default().with_do_not_add_h(),
    )
    .unwrap();
    let s = engine.structure();
    println!("atoms 17-23 region (working):");
    for (f, t, o) in &s.bonds {
        if [16, 17, 18, 19, 20, 21, 22].contains(f) || [16, 17, 18, 19, 20, 21, 22].contains(t) {
            println!("  bond {}-{} o{}", f + 1, t + 1, o);
        }
    }
    for i in 16..23 {
        let a = &s.atoms[i];
        println!("  atom{} {} q{} h{}", i + 1, a.element, a.charge, a.h_count);
    }
}
fn census() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let cases: Vec<serde_json::Value> = serde_json::from_str(
        &std::fs::read_to_string(format!("{dir}/corpus/fixtures/identity-breadth3.json")).unwrap(),
    )
    .unwrap();
    for c in cases {
        if c["id"] == "identity/breadth3-1053" {
            let engine = chemomer::Engine::with_options(
                c["molfile"].as_str().unwrap(),
                chemomer::Options::default().with_do_not_add_h(),
            )
            .unwrap();
            let s = engine.structure();
            let _counts = chemomer::layers::standardization::adjusted_hydrogens(s);
            // bond valences per atom over the input orders
            let mut bv = vec![0i32; s.atoms.len()];
            for &(f, t, o) in &s.bonds {
                bv[f] += (o as i32).min(3);
                bv[t] += (o as i32).min(3);
            }
            for (i, b) in bv.iter().enumerate() {
                let a = &s.atoms[i];
                if a.element == "C" && a.charge == 0 && *b == 2 {
                    println!("C-II census candidate: atom{} (bonds sum {b})", i + 1);
                }
            }
            println!("total atoms {}", s.atoms.len());
        }
    }
}
