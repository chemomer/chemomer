use chemomer::{Engine, Options};
fn main() {
    dump_aromatic();
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
    let c = &fixtures["identity/breadth-994"];
    let engine = Engine::with_options(
        c["molfile"].as_str().unwrap(),
        Options::default().with_do_not_add_h(),
    )
    .unwrap();
    let s = engine.structure();
    let counts = chemomer::layers::standardization::adjusted_hydrogens(s);
    let _ranks = chemomer::layers::hydrogen::canonical_ranks_for(s);
    let (adj, inv) = chemomer::canonical::graph_of(s);
    let groups = chemomer::canonical::main_canon_groups(s, &counts);
    let (r, sym) = chemomer::canonical::canonical_ranks_and_symmetry(&adj, &inv, &counts, &groups);
    for a in 0..s.atoms.len() {
        let neighbors: Vec<String> = s
            .bonds
            .iter()
            .filter_map(|&(f, t, o)| {
                if f == a {
                    Some(format!("{}(o{o})", t + 1))
                } else if t == a {
                    Some(format!("{}(o{o})", f + 1))
                } else {
                    None
                }
            })
            .collect();
        println!(
            "atom{} rank{} sym{} bonds {:?}",
            a + 1,
            r[a],
            sym[a],
            neighbors
        );
    }
}
// appended: dump aromatic bonds
fn dump_aromatic() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let cases: Vec<serde_json::Value> = serde_json::from_str(
        &std::fs::read_to_string(format!("{dir}/corpus/fixtures/identity-breadth.json")).unwrap(),
    )
    .unwrap();
    for c in cases {
        if c["id"] == "identity/breadth-994" {
            let engine = chemomer::Engine::with_options(
                c["molfile"].as_str().unwrap(),
                chemomer::Options::default().with_do_not_add_h(),
            )
            .unwrap();
            let s = engine.structure();
            let aromatic = chemomer::canonical::aromatic::aromatic_bonds(s);
            println!(
                "aromatic bonds: {:?}",
                aromatic.iter().map(|b| b + 1).collect::<Vec<_>>()
            );
            for (i, &(f, t, o)) in s.bonds.iter().enumerate() {
                println!("bond{} {}-{} o{}", i + 1, f + 1, t + 1, o);
            }
        }
    }
}
