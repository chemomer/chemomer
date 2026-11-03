// Which bonds our double_bond_centers marks on given fixtures.
use chemomer::{Engine, Options};

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let ids: Vec<String> = std::env::args().skip(1).collect();
    let ids = if ids.is_empty() {
        vec![
            "identity/breadth2-1325".to_string(),
            "identity/breadth2-1960".to_string(),
            "identity/breadth3-550".to_string(),
            "identity/breadth-1417".to_string(),
        ]
    } else {
        ids
    };
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
    for id in &ids {
        let c = &fixtures[id];
        let molfile = c["molfile"].as_str().unwrap();
        let engine = Engine::with_options(molfile, Options::default().with_do_not_add_h()).unwrap();
        let s = engine.structure();
        let counts = chemomer::layers::standardization::adjusted_hydrogens(s);
        let ranks = chemomer::layers::hydrogen::canonical_ranks_for(s);
        let (adjacency, invariants) = chemomer::canonical::graph_of(s);
        let groups = chemomer::canonical::main_canon_groups(s, &counts);
        let symmetry = chemomer::canonical::canonical_ranks_and_symmetry(
            &adjacency,
            &invariants,
            &counts,
            &groups,
        )
        .1;
        let centers =
            chemomer::canonical::stereo::center_inputs(s, &counts, &ranks, &symmetry, false, true);
        let marked: Vec<String> = centers
            .iter()
            .filter(|(_, l, _)| l.len() == 1)
            .map(|(a, l, b)| format!("{}-{}(bit{b})", a + 1, l[0] + 1))
            .collect();
        println!("{id}: {} bond marks: {}", marked.len(), marked.join(" "));
    }
}
