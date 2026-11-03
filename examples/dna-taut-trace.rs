fn main() {
    let id = std::env::args().nth(1).expect("id");
    for file in std::fs::read_dir("/Users/mulgogi/src/asciichem/wt/inchi-tests/corpus/fixtures")
        .expect("dir")
    {
        let p = file.expect("e").path();
        let n = p.file_name().unwrap().to_string_lossy().to_string();
        if !n.starts_with("identity") || !n.ends_with(".json") {
            continue;
        }
        let cases: Vec<serde_json::Value> =
            serde_json::from_str(&std::fs::read_to_string(&p).expect("read")).expect("json");
        for case in cases {
            if case["id"].as_str().unwrap() == format!("identity/{id}") {
                let engine = chemomer::Engine::with_options(
                    case["molfile"].as_str().unwrap(),
                    chemomer::Options::default().with_do_not_add_h(),
                )
                .unwrap();
                println!("inchi: {}", engine.inchi_string().unwrap());
                println!("aux: {}", engine.aux_string().unwrap());
                break;
            }
        }
    }
}
