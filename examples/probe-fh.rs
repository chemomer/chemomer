//! Probe: the fixed-H identifier for a molfile path.
use std::fs;

use chemomer::{Engine, Options};

fn main() {
    let input = fs::read_to_string(std::env::args().nth(1).expect("molfile")).expect("read");
    let engine = Engine::with_options(&input, Options::default().with_fixed_h()).expect("engine");
    println!("{}", engine.inchi_string().expect("identifier"));
}
