// Probe: the 1B × option interplay matrix — one mode word and a
// molfile path on the command line; prints the identifier, the key
// line and the aux line the engine assembles under that option set.
// Development-time tooling for the oneb-interplay unit.
//
//   cargo run --example probe-oneb -- large+fixedh /tmp/oneb-probe/large.mol

use chemomer::options::{FrameShift, PolymersMode};
use chemomer::Options;

fn mode_options(mode: &str) -> Options {
    let mut options = Options::default();
    for part in mode.split('+') {
        match part {
            "polymers" => options.polymers = PolymersMode::Modern,
            "polymers105" => options.polymers = PolymersMode::Legacy,
            "npzz" => options.npzz = true,
            "fs-none" => options.frame_shift = FrameShift::None,
            "fs-stars" => options.frame_shift = FrameShift::Opened,
            "foldcru" => options.fold_cru = true,
            "snon" => options = options.with_no_stereo(),
            "fixedh" => options = options.with_fixed_h(),
            "donotaddh" => options = options.with_do_not_add_h(),
            "recmet" => options = options.with_rec_met(),
            "srel" => options = options.with_relative_stereo(),
            "srac" => options = options.with_racemic_stereo(),
            "sucf" => options = options.with_chiral_flag_stereo(),
            "suu" => options = options.with_include_unknown_stereo(),
            "sluud" => options = options.with_distinct_unknown_undefined(),
            "ket" => options = options.with_ket(),
            "15t" => options = options.with_15t(),
            "pt22" => options = options.with_pt_22_00(),
            "pt16" => options = options.with_pt_16_00(),
            "pt06" => options = options.with_pt_06_00(),
            "pt39" => options = options.with_pt_39_00(),
            "pt13" => options = options.with_pt_13_00(),
            "pt18" => options = options.with_pt_18_00(),
            "satzz" => options = options.with_stereo_at_zz(),
            _ => {}
        }
    }
    options
}

fn main() {
    let mode = std::env::args().nth(1).expect("mode word");
    let path = std::env::args().nth(2).expect("molfile path");
    let molfile = std::fs::read_to_string(&path).expect("molfile");
    let options = mode_options(&mode);
    match chemomer::Engine::with_options(&molfile, options) {
        Ok(engine) => {
            println!("{}", engine.inchi_string().expect("identifier"));
            match engine.aux_string() {
                Ok(aux) => println!("{aux}"),
                Err(error) => println!("aux ERR: {error}"),
            }
            println!("{}", engine.inchikey_output().expect("key output"));
        }
        Err(error) => println!("error: {}", error.message()),
    }
}
