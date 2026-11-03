//! Probe: run the engine on a molfile with polymer mode flags.
//! Usage: probe-poly <molfile> [modes...] where modes are the gate's
//! mode words: polymers, polymers105, npzz, fs-none, fs-stars, foldcru.
use std::fs;

use chemomer::options::{FrameShift, PolymersMode};
use chemomer::{Engine, Options};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let input = fs::read_to_string(&args[1]).expect("read");
    let mut options = Options::default();
    for word in &args[2..] {
        match word.as_str() {
            "polymers" => options.polymers = PolymersMode::Modern,
            "polymers105" => options.polymers = PolymersMode::Legacy,
            "npzz" => options.npzz = true,
            "fs-none" => options.frame_shift = FrameShift::None,
            "fs-stars" => options.frame_shift = FrameShift::Opened,
            "foldcru" => options.fold_cru = true,
            _ => {}
        }
    }
    match Engine::with_options(&input, options) {
        Ok(engine) => {
            println!("{}", engine.inchi_string().expect("identifier"));
            println!("{}", engine.aux_string().expect("aux"));
        }
        Err(error) => println!("ERROR: {}", error.message()),
    }
}
