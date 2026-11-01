//! Probe: run the engine on a molfile path given as arg; prints the
//! assembled identifier and the aux line.
use std::fs;

use chemomer::Engine;

fn main() {
    let path = std::env::args().nth(1).expect("molfile path");
    let input = fs::read_to_string(path).expect("read");
    let engine = Engine::from_molfile(&input).expect("parses");
    println!("{}", engine.inchi_string().expect("identifier"));
    println!("{}", engine.aux_string().expect("aux"));
}
