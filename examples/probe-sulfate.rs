// sulfate-2 across the two modes: working structures and aux ranks.
use chemomer::{Engine, Options};

fn describe(tag: &str, engine: &Engine) {
    let s = engine.structure();
    println!("--- {tag}: {} atoms", s.atoms.len());
    for (i, a) in s.atoms.iter().enumerate() {
        println!(
            "  atom{} {} charge={} h={} pinned={}",
            i + 1,
            a.element,
            a.charge,
            a.h_count,
            a.h_pinned
        );
    }
    for (f, t, o) in &s.bonds {
        println!("  bond {}-{} order={}", f + 1, t + 1, o);
    }
    let counts = chemomer::layers::standardization::adjusted_hydrogens(s);
    let ranks = chemomer::layers::hydrogen::canonical_ranks_for(s);
    println!("  h counts: {:?}", counts);
    println!("  identifier ranks: {:?}", ranks);
}

fn main() {
    let dir = std::env::var("INCHI_TESTS").unwrap_or_else(|_| "../inchi-tests".into());
    let mut fixture = None;
    for name in ["identity.json"] {
        let cases: Vec<serde_json::Value> = serde_json::from_str(
            &std::fs::read_to_string(format!("{dir}/corpus/fixtures/{name}")).unwrap(),
        )
        .unwrap();
        for c in cases {
            if c["id"] == "identity/sulfate-2" {
                fixture = Some(c);
            }
        }
    }
    let c = fixture.unwrap();
    let molfile = c["molfile"].as_str().unwrap();
    let std_engine = Engine::from_molfile(molfile).unwrap();
    describe("standard", &std_engine);
    let dna_engine = Engine::with_options(molfile, Options::default().with_do_not_add_h()).unwrap();
    describe("donotaddh", &dna_engine);
    println!("std aux:  {}", std_engine.aux_string().unwrap());
    println!("dna aux:  {}", dna_engine.aux_string().unwrap());
    println!("oracle:   {}", c["aux"].as_str().unwrap());
    println!(
        "oracle dna aux: {}",
        c["donotaddh"]["aux"].as_str().unwrap()
    );
}
