//! Probe: the working structure's group census for breadth-305 under dna.
use std::fs;

fn main() {
    let text = fs::read_to_string("/tmp/b305.mol").unwrap();
    let options = chemomer::Options::default().with_do_not_add_h();
    let engine = chemomer::Engine::with_options(&text, options).unwrap();
    println!("identifier: {}", engine.inchi_string().unwrap());
    let aux = engine.aux_string().unwrap();
    println!("aux: {aux}");
}
