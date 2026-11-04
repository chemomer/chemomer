// Point checks for the 13 ledgered stereo-walk/TAUT survivors: the
// engine's DoNotAddH identifier+aux against the oracle recordings.
use chemomer::{Engine, Options};

fn main() {
    let dir =
        std::path::Path::new("/Users/mulgogi/src/asciichem/wt/corpus-bns-aux/corpus/fixtures");
    let ids: Vec<&str> = vec![
        "identity/breadth-1231",
        "identity/breadth-27",
        "identity/breadth2-1082",
        "identity/breadth2-1536",
        "identity/breadth2-1549",
        "identity/breadth3-1851",
        "identity/breadth3-333",
        "identity/breadth3-336",
        "identity/breadth3-487",
        "identity/breadth3-501",
        "identity/breadth3-961",
        "identity/breadth-410",
        "identity/breadth2-1937",
    ];
    for file in [
        "identity.json",
        "identity-breadth.json",
        "identity-breadth2.json",
        "identity-breadth3.json",
    ] {
        let path = dir.join(file);
        let Ok(text) = std::fs::read_to_string(&path) else {
            continue;
        };
        let Ok(cases) = serde_json::from_str::<Vec<serde_json::Value>>(&text) else {
            continue;
        };
        for case in &cases {
            let id = case["id"].as_str().unwrap_or("");
            if !ids.contains(&id) {
                continue;
            }
            let molfile = case["molfile"].as_str().unwrap_or("");
            let recorded_aux = case["donotaddh"]["aux"].as_str().unwrap_or("");
            let recorded_inchi = case["donotaddh"]["inchi"].as_str().unwrap_or("");
            let engine = Engine::with_options(molfile, Options::default().with_do_not_add_h());
            match engine {
                Err(e) => println!("{id}: ENGINE ERR {e}"),
                Ok(engine) => {
                    let aux = engine.aux_string().unwrap_or_else(|e| format!("ERR {e}"));
                    let inchi = engine.inchi_string().unwrap_or_else(|e| format!("ERR {e}"));
                    let aux_ok = aux == recorded_aux;
                    let inchi_ok = inchi == recorded_inchi;
                    println!(
                        "{id}: aux {} inchi {}",
                        if aux_ok { "GREEN" } else { "RED" },
                        if inchi_ok { "GREEN" } else { "RED" }
                    );
                    if !aux_ok {
                        println!("  ours:   {aux}");
                        println!("  oracle: {recorded_aux}");
                    }
                    if !inchi_ok {
                        println!("  ours-i:   {inchi}");
                        println!("  oracle-i: {recorded_inchi}");
                    }
                }
            }
        }
    }
}
