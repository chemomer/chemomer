// Dump atoms/bonds around given indices for the two over-marked fixtures.
use chemomer::{Engine, Options};

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let ids = [
        ("identity/breadth2-1325", vec![13usize, 14]),
        ("identity/breadth2-1960", vec![13, 14, 17, 18]),
    ];
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
    for (id, atoms) in ids {
        let c = &fixtures[id];
        let molfile = c["molfile"].as_str().unwrap();
        let engine = Engine::with_options(molfile, Options::default().with_do_not_add_h()).unwrap();
        let s = engine.structure();
        println!("== {id}");
        for (i, a) in s.atoms.iter().enumerate() {
            if atoms.contains(&(i)) || atoms.contains(&(i - 1)) {
                println!(
                    "  atom{} {} charge={} h={}",
                    i + 1,
                    a.element,
                    a.charge,
                    a.h_count
                );
            }
        }
        for (f, t, o) in &s.bonds {
            if atoms.contains(f) || atoms.contains(t) {
                println!("  bond {}-{} order={}", f + 1, t + 1, o);
            }
        }
    }
}
