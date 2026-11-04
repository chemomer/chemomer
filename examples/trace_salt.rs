// Trace the salt-phase probes the aux path sees (INCHI_TRACE_BNS on).
fn main() {
    let ids: Vec<String> = std::env::args().skip(1).collect();
    for file in std::fs::read_dir("/Users/mulgogi/src/asciichem/wt/corpus-bns-aux/corpus/fixtures")
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
            let id = case["id"].as_str().unwrap().to_string();
            if !ids.iter().any(|i| id.ends_with(i)) {
                continue;
            }
            eprintln!("##### {id}");
            let engine = chemomer::Engine::from_molfile(case["molfile"].as_str().unwrap()).unwrap();
            let s = engine.structure();
            let counts = chemomer::layers::standardization::adjusted_hydrogens(s);
            let ranks = chemomer::layers::hydrogen::canonical_ranks_for(s);
            let groups = chemomer::canonical::tautomer::mobile_groups(s, &counts, &ranks);
            let real_h = |g: &chemomer::canonical::tautomer::MobileGroup| g.count > g.charge;
            eprintln!(
                "SUMMARY {id} realgroups={} groups={:?}",
                groups.iter().any(&real_h),
                groups
                    .iter()
                    .map(|g| (g.count, g.charge, g.ranks.clone()))
                    .collect::<Vec<_>>()
            );
        }
    }
}
