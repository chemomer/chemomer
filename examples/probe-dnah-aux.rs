// Probe: the DoNotAddH aux tail — the corpus records donotaddh.aux for
// fixtures carrying the field; no gate compares it yet. Check the
// engine's aux under -DoNotAddH against every recording.
use std::fs;

fn suite_dir() -> std::path::PathBuf {
    std::env::var("INCHI_TESTS")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| std::path::PathBuf::from("../inchi-tests"))
}

fn main() {
    let dir = suite_dir();
    let fixtures = dir.join("corpus/fixtures");
    let mut files: Vec<std::path::PathBuf> = fs::read_dir(&fixtures)
        .expect("corpus fixtures dir")
        .filter_map(|e| {
            let p = e.expect("dir entry").path();
            let n = p.file_name().expect("name").to_string_lossy().to_string();
            (n.starts_with("identity") && n.ends_with(".json")).then_some(p)
        })
        .collect();
    files.sort();
    let mut with_aux = 0;
    let mut pass = 0;
    let mut fail = 0;
    for file in &files {
        let cases: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(file).expect("read")).expect("valid");
        for case in &cases {
            let Some(recorded) = case["donotaddh"]
                .as_object()
                .and_then(|s| s.get("aux"))
                .and_then(|v| v.as_str())
            else {
                continue;
            };
            with_aux += 1;
            let id = case["id"].as_str().expect("id");
            let molfile = case["molfile"].as_str().expect("molfile");
            match chemomer::Engine::with_options(
                molfile,
                chemomer::Options::default().with_do_not_add_h(),
            )
            .and_then(|e| e.aux_string())
            {
                Ok(a) if a == recorded => pass += 1,
                Ok(a) => {
                    fail += 1;
                    println!("{id}:\n  got : {a}\n  want: {recorded}");
                }
                Err(e) => {
                    fail += 1;
                    println!("{id}: ERR {e}");
                }
            }
        }
    }
    println!("\ndonotaddh aux: with recording {with_aux}, pass {pass}, fail {fail}");
}
