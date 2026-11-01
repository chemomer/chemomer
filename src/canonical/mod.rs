//! Canonical ranking — a native port of the InChI reference's base
//! canonicalization path (MIT; read at development time, never
//! compiled, vendored, or called). The corpus arbitrates: reproduced
//! rank orders must explain every recorded `/c`.
//!
//! Pipeline: invariants → initial partition → equitable refinement →
//! individualization search over connection tables → `/c` emission
//! from the winning table.
//!
//! The automorphism-orbit pruning of the reference (its `CanonGraph`
//! L10-L13) is a speed optimization only: the winning terminal — the
//! first search-order terminal whose table no later terminal exceeds
//! — is the same without it, because pruned subtrees only ever hold
//! tables equal to an already-found terminal.

pub mod aromatic;
pub mod blossom;
pub(crate) mod bns;
pub mod bond_stereo;
pub mod emit;
pub mod fixedh;
mod invariants;
pub mod metals;
mod partition;
pub mod protons;
mod search;
pub mod stereo;
mod table;
pub mod tautomer;

pub use emit::connection_segment;
pub use partition::Partition;
pub use search::{
    base_symmetry, canonical_ranks, canonical_ranks_and_symmetry,
    canonical_ranks_and_symmetry_stereo, canonical_table,
};
pub use tautomer::group_membership;

use crate::graph::Structure;

/// Atom adjacency and reference invariants for canonicalization.
/// Explicit hydrogen atoms never reach here — they are folded onto
/// their carriers at parse time, so the numbering is heavy-atom-only.
pub fn graph_of(structure: &Structure) -> (Vec<Vec<usize>>, Vec<(u16, u16)>) {
    let n = structure.atoms.len();
    let mut adjacency = vec![Vec::new(); n];
    for &(from, to, _) in &structure.bonds {
        adjacency[from].push(to);
        adjacency[to].push(from);
    }
    let elements: Vec<String> = structure.atoms.iter().map(|a| a.element.clone()).collect();
    let invariants = invariants::atom_invariants(&elements, &adjacency);
    (adjacency, invariants)
}

/// The equitable refinement of the invariant partition — the M2
/// census view over the same machinery the search starts from.
pub fn equitable_partition(structure: &Structure) -> Partition {
    let (adjacency, invariants) = graph_of(structure);
    let mut partition = Partition::initial(&invariants);
    partition.refine(&adjacency);
    partition
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{AtomRecord, Structure};

    fn atom(element: &str) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count: 0,
            h_pinned: false,
            isotope: 0,
            charge: 0,
        }
    }

    fn structure(atoms: Vec<AtomRecord>, bonds: Vec<(usize, usize, u32)>) -> Structure {
        Structure {
            atoms,
            bonds,
            ..Default::default()
        }
    }

    #[test]
    fn benzene_keeps_one_symmetric_class() {
        let atoms = (0..6).map(|_| atom("C")).collect();
        let bonds = (0..6).map(|i| (i, (i + 1) % 6, 4)).collect();
        let partition = equitable_partition(&structure(atoms, bonds));
        assert_eq!(partition.cell_count(), 1);
        assert!(!partition.is_discrete());
    }

    #[test]
    fn ethanol_refines_to_a_total_order() {
        let partition = equitable_partition(&structure(
            vec![atom("C"), atom("C"), atom("O")],
            vec![(0, 1, 1), (1, 2, 1)],
        ));
        assert_eq!(partition.cell_count(), 3);
        assert!(partition.is_discrete());
    }

    #[test]
    fn refinement_is_independent_of_atom_order() {
        // Isopropanol with the atom order rotated: old atom i becomes
        // atom (i + 1) % 4.
        let classes = |partition: &Partition| partition.cell_of_each_atom();
        let rotated = |index: usize| (index + 1) % 4;
        let a = equitable_partition(&structure(
            vec![atom("C"), atom("C"), atom("C"), atom("O")],
            vec![(0, 1, 1), (0, 2, 1), (0, 3, 1)],
        ));
        let b = equitable_partition(&structure(
            vec![atom("O"), atom("C"), atom("C"), atom("C")],
            vec![
                (rotated(0), rotated(1), 1),
                (rotated(0), rotated(2), 1),
                (rotated(0), rotated(3), 1),
            ],
        ));
        let (ca, cb) = (classes(&a), classes(&b));
        for left in 0..4 {
            for right in 0..4 {
                assert_eq!(
                    ca[left] == ca[right],
                    cb[rotated(left)] == cb[rotated(right)],
                    "partition changed for atoms {left}/{right}"
                );
            }
        }
    }
}
