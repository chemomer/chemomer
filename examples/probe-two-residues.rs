//! Dump our engine output vs the oracle recording for the two dated
//! residues (identity/breadth3-487, identity/breadth-410) across the
//! DoNotAddH mode.

use std::fs;
use std::path::PathBuf;

use chemomer::Engine;

fn dump(dir: &PathBuf, fixture: &str, id: &str) {
    let path = dir.join("corpus/fixtures").join(fixture);
    let cases: Vec<serde_json::Value> =
        serde_json::from_str(&fs::read_to_string(&path).expect("read")).expect("json");
    for case in cases {
        if case["id"].as_str() != Some(id) {
            continue;
        }
        let molfile = case["molfile"].as_str().unwrap();
        let engine =
            Engine::with_options(molfile, chemomer::Options::default().with_do_not_add_h())
                .unwrap();
        match engine.inchi_string() {
            Ok(s) => println!("{id} identifier:\n  ours   {s}"),
            Err(e) => println!("{id} identifier errored: {e}"),
        }
        match engine.aux_string() {
            Ok(s) => println!("  ours-a {s}"),
            Err(e) => println!("  ours-a errored: {e}"),
        }
        if let Some(rec) = case["donotaddh"].as_object() {
            println!("  oracle  {}", rec["inchi"].as_str().unwrap_or("(none)"));
            println!("  oracle-a {}", rec["aux"].as_str().unwrap_or("(none)"));
        }
    }
}

fn main() {
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    dump(&dir, "identity-breadth3.json", "identity/breadth3-487");
    dump(&dir, "identity-breadth.json", "identity/breadth-410");
}
