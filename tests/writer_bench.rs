// Throwaway A/B: wall-clock write over the identity corpus.

#[test]
fn writer_wall_clock_bench() {
    let corpus = std::env::var("INCHI_TESTS").expect("INCHI_TESTS");
    let dir = std::path::Path::new(&corpus).join("corpus/fixtures");
    let mut smiles: Vec<String> = Vec::new();
    for name in [
        "identity.json",
        "identity-breadth.json",
        "identity-metal.json",
        "identity-multicomponent.json",
    ] {
        let Ok(text) = std::fs::read_to_string(dir.join(name)) else {
            continue;
        };
        let cases: Vec<serde_json::Value> = serde_json::from_str(&text).expect("json");
        for case in cases.iter() {
            if let Some(s) = case.get("smiles").and_then(|v| v.as_str()) {
                smiles.push(s.to_string());
            }
        }
    }
    let structures: Vec<_> = smiles
        .iter()
        .filter_map(|s| chemomer::smiles::read(s).ok())
        .collect();
    let n = structures.len();
    let start = std::time::Instant::now();
    let mut bytes = 0usize;
    for structure in &structures {
        bytes += chemomer::smiles::write(structure).map(|s| s.len()).unwrap_or(0);
    }
    let elapsed = start.elapsed();
    println!(
        "WRITER BENCH: {n} rows, {elapsed:?} total, {:.2} ms/row, {bytes} bytes",
        elapsed.as_secs_f64() * 1000.0 / n as f64
    );
}
