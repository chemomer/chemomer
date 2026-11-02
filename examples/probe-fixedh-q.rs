// Probe: the acetate anion ALONE — standard vs FixedH.
use chemomer::{Engine, Options};

fn show(label: &str, molfile: &str) {
    let std = Engine::from_molfile(molfile)
        .and_then(|e| e.inchi_string())
        .map_err(|e| e.message().to_string());
    let fh = Engine::with_options(molfile, Options::default().with_fixed_h())
        .and_then(|e| e.inchi_string())
        .map_err(|e| e.message().to_string());
    let fh_aux = Engine::with_options(molfile, Options::default().with_fixed_h())
        .and_then(|e| e.aux_string())
        .map_err(|e| e.message().to_string());
    println!("{label}:\n  std: {std:?}\n  fh:  {fh:?}\n  fha: {fh_aux:?}");
}

fn main() {
    let acetate = "acetate\n  chemomer\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   12.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   62.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1\n  2  3  2\n  2  4  1\nM  CHG  1   4  -1\nM  END\n";
    show("acetate anion alone", acetate);
    let acid = "acetic acid\n  chemomer\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   12.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   62.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1\n  2  3  2\n  2  4  1\nM  END\n";
    show("acetic acid alone", acid);
}
