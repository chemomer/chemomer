//! Debug: engine identifier per probe molfile, for oracle diffing.
use std::fs;

use chemomer::Engine;

fn main() {
    for path in std::env::args().skip(1) {
        let input = fs::read_to_string(&path).expect("read");
        let name = std::path::Path::new(&path)
            .file_stem()
            .unwrap()
            .to_string_lossy();
        match Engine::from_molfile(&input).and_then(|e| e.inchi_string()) {
            Ok(inchi) => println!("{name}\n  {inchi}"),
            Err(e) => println!("{name}\n  ERROR: {}", e.message()),
        }
    }
}
