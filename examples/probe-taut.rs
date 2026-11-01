//! Probe: run the engine on a molfile under the experimental
//! tautomerism flags — comma-separated second arg (KET, 15T,
//! PT_06_00, ...); prints the assembled identifier.
use std::fs;

use chemomer::{Engine, Options};

fn main() {
    let path = std::env::args().nth(1).expect("molfile path");
    let flags = std::env::args().nth(2).unwrap_or_default();
    let input = fs::read_to_string(path).expect("read");
    let mut options = Options::default();
    for flag in flags.split(',').filter(|f| !f.is_empty()) {
        options = match flag {
            "KET" => options.with_ket(),
            "15T" => options.with_15t(),
            "PT_22_00" => options.with_pt_22_00(),
            "PT_16_00" => options.with_pt_16_00(),
            "PT_06_00" => options.with_pt_06_00(),
            "PT_39_00" => options.with_pt_39_00(),
            "PT_13_00" => options.with_pt_13_00(),
            "PT_18_00" => options.with_pt_18_00(),
            other => panic!("unknown flag {other}"),
        };
    }
    let engine = Engine::with_options(&input, options).expect("parses");
    println!("{}", engine.inchi_string().expect("identifier"));
}
