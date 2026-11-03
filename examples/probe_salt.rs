// Salt-bit scan: per fixture, our mobile groups + expected aux digit.
fn main() {
    let dir = std::env::args().nth(1).expect("corpus fixtures dir");
    for entry in std::fs::read_dir(&dir).expect("dir") {
        let path = entry.expect("e").path();
        let name = path.file_name().expect("n").to_string_lossy().to_string();
        if !name.ends_with(".json") { continue; }
        let text = std::fs::read_to_string(&path).expect("read");
        let cases: Vec<serde_json::Value> = serde_json::from_str(&text).expect("json");
        for case in cases {
            let id = case["id"].as_str().expect("id");
            let Some(aux) = case["aux"].as_str() else { continue };
            let Some(rest) = aux.strip_prefix("AuxInfo=1/") else { continue };
            let Some(slash) = rest.find('/') else { continue };
            let expected: i32 = rest[..slash].parse().expect("digit");
            let molfile = case["molfile"].as_str().expect("molfile");
            let mut working = match chemomer::parse_molfile(molfile) { Ok(w) => w, Err(_) => continue };
            chemomer::canonical::metals::materialize_hydrogens(&mut working);
            let counts = chemomer::layers::standardization::adjusted_hydrogens(&working);
            let (adj, inv) = chemomer::canonical::graph_of(&working);
            let ranks = chemomer::canonical::canonical_ranks(&adj, &inv, &counts, &[]);
            let groups = chemomer::canonical::tautomer::mobile_groups(&working, &counts, &ranks);
            let real_h = |g: &chemomer::canonical::tautomer::MobileGroup| g.count > g.charge;
            let taut = groups.iter().any(|g| real_h(g));
            let salt_ours = groups.iter().any(|g| real_h(g) && g.charge > 0);
            let hard = working.protons_hard_moved;
            eprintln!("SALTSCAN\t{}\t{}\t{}\t{}\t{}\t{}", id, expected, taut as i32, salt_ours as i32, hard,
                groups.iter().map(|g| format!("(n{}c{}m{})", g.count, g.charge, g.ranks.len())).collect::<Vec<_>>().join(","));
        }
    }
}
