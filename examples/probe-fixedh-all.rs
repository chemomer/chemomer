// Probe: every multi/metal fixture's FixedH identifier + aux against
// the recorded oracle, ledgered or not — the decode loop's diff feed.
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
    let ledger: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(dir.join("corpus/aux-triage.json")).expect("aux ledger"),
    )
    .expect("valid ledger");
    let mut id_pass = 0;
    let mut aux_pass = 0;
    let mut id_fail = 0;
    let mut aux_fail = 0;
    for file in files {
        let cases: Vec<serde_json::Value> =
            serde_json::from_str(&fs::read_to_string(&file).expect("read")).expect("valid");
        for case in &cases {
            let Some(fh) = case["fixedh"].as_object() else {
                continue;
            };
            let id = case["id"].as_str().expect("id");
            let molfile = case["molfile"].as_str().expect("molfile");
            let ledgered = ledger.get(id).is_some();
            let engine = match chemomer::Engine::with_options(
                molfile,
                chemomer::Options::default().with_fixed_h(),
            ) {
                Ok(e) => e,
                Err(error) => {
                    println!("{id} [ledgered={ledgered}]: ENGINE ERR {error}");
                    id_fail += 1;
                    continue;
                }
            };
            let assembled = engine.inchi_string();
            let aux = engine.aux_string();
            let rec_inchi = fh.get("inchi").and_then(|v| v.as_str());
            let rec_aux = fh.get("aux").and_then(|v| v.as_str());
            let mut ok_id = true;
            let mut ok_aux = true;
            match (assembled, rec_inchi) {
                (Ok(a), Some(r)) if a == r => id_pass += 1,
                (Ok(a), Some(r)) => {
                    ok_id = false;
                    id_fail += 1;
                    println!("{id} [ledgered={ledgered}]:\n  got : {a}\n  want: {r}");
                }
                (Ok(_), None) => println!("{id} [ledgered={ledgered}]: no recorded inchi"),
                (Err(e), Some(r)) => {
                    ok_id = false;
                    id_fail += 1;
                    println!("{id} [ledgered={ledgered}]: ERR {e} (want {r})");
                }
                (Err(e), None) => {
                    ok_id = false;
                    id_fail += 1;
                    println!("{id} [ledgered={ledgered}]: ERR {e} (no recording)")
                }
            }
            match (aux, rec_aux) {
                (Ok(a), Some(r)) if a == r => aux_pass += 1,
                (Ok(a), Some(r)) => {
                    ok_aux = false;
                    aux_fail += 1;
                    println!("{id} [ledgered={ledgered}]:\n  gotA: {a}\n  wantA: {r}");
                }
                (Ok(_), None) => {}
                (Err(e), Some(r)) => {
                    ok_aux = false;
                    aux_fail += 1;
                    println!("{id} [ledgered={ledgered}]: AUX ERR {e} (want {r})");
                }
                (Err(e), None) => {
                    ok_aux = false;
                    aux_fail += 1;
                    println!("{id} [ledgered={ledgered}]: AUX ERR {e} (no recording)")
                }
            }
            let _ = (ok_id, ok_aux);
        }
    }
    println!("\nfixedh identifier: pass {id_pass}, fail {id_fail}");
    println!("fixedh aux: pass {aux_pass}, fail {aux_fail}");
}
