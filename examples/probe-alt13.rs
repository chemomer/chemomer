use chemomer::{Engine, Options};
fn main() {
    let cases = [
        ("propyne", "C3", "InChI=1S/C3/c1-3-2"),
        ("butadyne", "C4", "InChI=1S/C4/c1-3-4-2"),
        ("propynylam", "C3N", "InChI=1S/C3N/c1-2-3-4"),
        ("allene", "C3", "InChI=1S/C3/c1-3-2"),
        ("butenyne", "C4", "InChI=1S/C4/c1-3-4-2"),
    ];
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let _ = dir;
    for (name, _formula, oracle) in cases {
        let path = format!("/tmp/sbprobes/{name}.mol");
        let molfile = std::fs::read_to_string(&path).unwrap();
        for (label, options) in [
            ("dna", Options::default().with_do_not_add_h()),
            ("std", Options::default()),
        ] {
            match Engine::with_options(&molfile, options).and_then(|e| e.inchi_string()) {
                Ok(inchi) => {
                    let mark = if inchi == oracle { "MATCH" } else { "DIFF" };
                    println!("{name} [{label}]: {mark} {inchi}");
                }
                Err(e) => println!("{name} [{label}]: ERR {e}"),
            }
        }
    }
}
