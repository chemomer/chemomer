// Dump atoms/bonds around given indices for given fixtures.
use chemomer::{Engine, Options};

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    // fixture id -> atom indices of the bond we mark
    let ids = [
        ("identity/breadth-994", [12usize, 13]),
        ("identity/breadth2-525", [8, 9]),
        ("identity/breadth2-718", [16, 17]),
        ("identity/breadth2-846", [16, 17]),
        ("identity/breadth2-1122", [12, 13]),
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
            if atoms.contains(&i) {
                let neighbors: Vec<String> = s
                    .bonds
                    .iter()
                    .filter_map(|&(f, t, o)| {
                        if f == i {
                            Some(format!("{}(o{o})", t + 1))
                        } else if t == i {
                            Some(format!("{}(o{o})", f + 1))
                        } else {
                            None
                        }
                    })
                    .collect();
                println!(
                    "  atom{} {} charge={} bonds {:?}",
                    i + 1,
                    a.element,
                    a.charge,
                    neighbors
                );
            }
        }
    }
}
