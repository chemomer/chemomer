// sulfate-2 graph through canonical_ranks_and_symmetry with tracing.
use chemomer::canonical;

fn main() {
    // atoms: 0=O1 1=S2 2=O3 3=O4 4=O5; bonds 1-2 s, 2-3 d, 2-4 s, 2-5 s
    let _adjacency = [vec![1], vec![0, 2, 3, 4], vec![1], vec![1], vec![1]];
    let elements: Vec<String> = ["O", "S", "O", "O", "O"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let invariants = chemomer::canonical::graph_of(&chemomer::Structure::default()).1; // dummy
    let _ = invariants;
    let _ = elements;
    // build invariants via the public helper on a real structure
    let structure = chemomer::Structure {
        atoms: vec![
            chemomer::graph::AtomRecord {
                element: "O".into(),
                h_count: 0,
                charge: 0,
                valence_code: 0,
                isotope: 0,
                h_pinned: false,
            },
            chemomer::graph::AtomRecord {
                element: "S".into(),
                h_count: 0,
                charge: 0,
                valence_code: 0,
                isotope: 0,
                h_pinned: false,
            },
            chemomer::graph::AtomRecord {
                element: "O".into(),
                h_count: 0,
                charge: 0,
                valence_code: 0,
                isotope: 0,
                h_pinned: false,
            },
            chemomer::graph::AtomRecord {
                element: "O".into(),
                h_count: 0,
                charge: 0,
                valence_code: 0,
                isotope: 0,
                h_pinned: false,
            },
            chemomer::graph::AtomRecord {
                element: "O".into(),
                h_count: 0,
                charge: 0,
                valence_code: 0,
                isotope: 0,
                h_pinned: false,
            },
        ],
        bonds: vec![(0, 1, 1), (1, 2, 2), (1, 3, 1), (1, 4, 1)],
        ..Default::default()
    };
    // The real DoNotAddH working structure: pinned zero H, charges.
    let mut structure = structure;
    for (i, atom) in structure.atoms.iter_mut().enumerate() {
        atom.h_pinned = true;
        if i == 0 || i == 4 {
            atom.charge = -1;
        }
    }
    let (adj, inv) = canonical::graph_of(&structure);
    println!("invariants: {:?}", inv);
    let partition = canonical::equitable_partition(&structure);
    println!("equitable ranks: {:?}", partition.ranks());
    let counts = chemomer::layers::standardization::adjusted_hydrogens(&structure);
    let mask = canonical::tautomer::endpoint_mask(&structure, &counts);
    let notaut = canonical::tautomer::without_mobile_h(&counts, &mask);
    let groups = canonical::group_membership(&structure, &counts);
    println!(
        "counts {:?} mask {:?} notaut {:?} groups {:?}",
        counts, mask, notaut, groups
    );
    let (ranks, symmetry) = canonical::canonical_ranks_and_symmetry(&adj, &inv, &notaut, &groups);
    println!("ranks: {:?}", ranks);
    println!("symmetry: {:?}", symmetry);
}
