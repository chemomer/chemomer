//! Trace the TAUT walk for one fixture: the individualization sequence,
//! the terminal order, and the final winner.

use std::fs;
use std::path::PathBuf;

use chemomer::{Engine, Options};

fn main() {
    let id = std::env::args().nth(1).expect("fixture id");
    let dir = std::env::var("INCHI_TESTS")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("../inchi-tests"));
    for file in [
        "identity.json",
        "identity-breadth.json",
        "identity-breadth2.json",
        "identity-breadth3.json",
    ] {
        let path = dir.join("corpus/fixtures").join(file);
        let Ok(text) = fs::read_to_string(&path) else {
            continue;
        };
        let Ok(cases) = serde_json::from_str::<Vec<serde_json::Value>>(&text) else {
            continue;
        };
        for case in cases {
            if case["id"].as_str() != Some(id.as_str()) {
                continue;
            }
            let molfile = case["molfile"].as_str().unwrap();
            let engine =
                Engine::with_options(molfile, Options::default().with_do_not_add_h()).unwrap();
            let aux = engine.aux_string().unwrap();
            eprintln!("AUX {aux}");
            return;
        }
    }
    panic!("fixture {id} not found");
}
