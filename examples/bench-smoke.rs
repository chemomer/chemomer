//! Throwaway smoke timer (delete after use): per-iteration Instant
//! timing; reports min / p50 / p90 per leg — the min is the
//! interference-robust floor on this loaded box.
use std::time::Instant;

const WATER: &str = "water\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
const ALANINE: &str = "l-alanine\n  t\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n   12.0000   22.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0  0\n   72.0000   22.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   67.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  132.0000   17.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  192.0000   62.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  0  0  0  0  0  0  0  0\n  2  3  1  0  0  0  0  0  0  0  0\n  4  5  2  0  0  0  0  0  0  0  0\n  2  4  1  0  0  0  0  0  0  0  0\n  4  6  1  0  0  0  0  0  0  0  0\nM  END\n";

fn time(name: &str, molfile: &str, n: usize) {
    for _ in 0..50 {
        let e = chemomer::Engine::from_molfile(molfile).unwrap();
        let _ = e.inchi_string().unwrap();
    }
    let mut samples: Vec<u32> = Vec::with_capacity(n);
    for _ in 0..n {
        let t = Instant::now();
        let e = chemomer::Engine::from_molfile(molfile).unwrap();
        let _ = e.inchi_string().unwrap();
        samples.push((t.elapsed().as_secs_f64() * 1e9) as u32);
    }
    samples.sort_unstable();
    let pick = |q: f64| samples[((samples.len() as f64 - 1.0) * q) as usize];
    println!(
        "{name}: min {} ns, p50 {} ns, p90 {} ns",
        pick(0.0),
        pick(0.5),
        pick(0.9)
    );
}

fn main() {
    time("water   inchi_string", WATER, 400);
    time("alanine inchi_string", ALANINE, 200);
}
