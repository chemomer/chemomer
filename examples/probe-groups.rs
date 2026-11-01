//! Probe: print the mobile groups (atoms, elements, counts) for a
//! molfile under tautomer flags (comma-separated second arg).
use std::fs;

use chemomer::{Engine, Options};

fn main() {
    let path = std::env::args().nth(1).expect("molfile path");
    let flags = std::env::args().nth(2).unwrap_or_default();
    let input = fs::read_to_string(path).expect("read");
    let mut options = Options::default();
    for flag in flags.split(',').filter(|f| !f.is_empty()) {
        options = match flag {
            "KET" => options.with_ket(),
            "15T" => options.with_15t(),
            "PT_22_00" => options.with_pt_22_00(),
            "PT_16_00" => options.with_pt_16_00(),
            "PT_06_00" => options.with_pt_06_00(),
            "PT_39_00" => options.with_pt_39_00(),
            "PT_13_00" => options.with_pt_13_00(),
            "PT_18_00" => options.with_pt_18_00(),
            other => panic!("unknown flag {other}"),
        };
    }
    let engine = Engine::with_options(&input, options).expect("parses");
    let s = engine.structure();
    let counts = chemomer::layers::standardization::adjusted_hydrogens(s);
    let (adjacency, invariants) = chemomer::canonical::graph_of(s);
    let mask = chemomer::canonical::tautomer::endpoint_mask(s, &counts);
    let memberships = chemomer::canonical::group_membership(s, &counts);
    let ranks = chemomer::canonical::canonical_ranks(
        &adjacency,
        &invariants,
        &chemomer::canonical::tautomer::without_mobile_h(&counts, &mask),
        &memberships,
    );
    let groups = chemomer::canonical::tautomer::mobile_groups(s, &counts, &ranks);
    for (atom, record) in s.atoms.iter().enumerate() {
        println!(
            "atom {atom} {} H={} charge={} rank={}",
            record.element, counts[atom], record.charge, ranks[atom]
        );
    }
    for group in &groups {
        let members: Vec<String> = group
            .ranks
            .iter()
            .map(|rank| {
                let atom = ranks.iter().position(|&r| r == *rank).unwrap_or(usize::MAX);
                format!("{}(rank {})", s.atoms[atom].element, rank)
            })
            .collect();
        println!(
            "group count={} charge={} {:?}",
            group.count, group.charge, members
        );
    }
    println!("inchi: {}", engine.inchi_string().unwrap());
}
