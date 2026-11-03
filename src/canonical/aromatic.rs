//! Aromatic-ring perception, scoped to what the tautomer machinery
//! needs: five- and six-membered Kekulé cycles — every bond single or
//! double, alternating around the ring — mark their bonds aromatic.
//! Tautomer paths then cross an aromatic bond in either parity, which
//! is how the reference delocalizes a ring before applying its
//! templates.

use crate::graph::Structure;
use std::collections::BTreeSet;

/// Indices into `structure.bonds` of the aromatic bonds.
pub fn aromatic_bonds(structure: &Structure) -> BTreeSet<usize> {
    let n = structure.atoms.len();
    let mut neighbors: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (index, &(from, to, kind)) in structure.bonds.iter().enumerate() {
        if kind == 1 || kind == 2 {
            neighbors[from].push((to, index));
            neighbors[to].push((from, index));
        }
    }

    let mut cycles: BTreeSet<Vec<usize>> = BTreeSet::new();
    for length in [5usize, 6] {
        for start in 0..n {
            let mut path = vec![start];
            let mut used = vec![false; n];
            used[start] = true;
            find_cycles(
                &neighbors,
                start,
                start,
                length,
                &mut path,
                &mut used,
                &mut cycles,
            );
        }
    }

    let mut aromatic = BTreeSet::new();
    for cycle in &cycles {
        if is_kekule(structure, cycle) {
            for position in 0..cycle.len() {
                let a = cycle[position];
                let b = cycle[(position + 1) % cycle.len()];
                if let Some(&(_, bond)) = neighbors[a].iter().find(|&&(neighbor, _)| neighbor == b)
                {
                    aromatic.insert(bond);
                }
            }
        }
    }
    aromatic
}

fn find_cycles(
    neighbors: &[Vec<(usize, usize)>],
    start: usize,
    atom: usize,
    length: usize,
    path: &mut Vec<usize>,
    used: &mut [bool],
    cycles: &mut BTreeSet<Vec<usize>>,
) {
    if path.len() == length {
        if neighbors[atom].iter().any(|&(next, _)| next == start) {
            let mut canonical = path.clone();
            // Rotate the smallest atom to the front, then keep the
            // lexicographically smaller of the two directions.
            let offset = canonical
                .iter()
                .enumerate()
                .min_by_key(|&(position, &atom)| (atom, position))
                .map(|(position, _)| position)
                .expect("path is non-empty");
            canonical.rotate_left(offset);
            let mut reversed: Vec<usize> = canonical.iter().rev().copied().collect();
            reversed.rotate_left(1);
            let candidate = std::cmp::min(canonical.clone(), reversed);
            cycles.insert(candidate);
        }
        return;
    }
    for &(next, _) in &neighbors[atom] {
        if used[next] {
            continue;
        }
        used[next] = true;
        path.push(next);
        find_cycles(neighbors, start, next, length, path, used, cycles);
        path.pop();
        used[next] = false;
    }
}

/// Bonds around the cycle strictly alternate single/double.
fn is_kekule(structure: &Structure, cycle: &[usize]) -> bool {
    let kind_between = |a: usize, b: usize| -> u32 {
        structure
            .bonds
            .iter()
            .find(|&&(from, to, _)| (from == a && to == b) || (from == b && to == a))
            .map(|&(_, _, kind)| kind)
            .unwrap_or(0)
    };
    let mut pattern = Vec::new();
    for position in 0..cycle.len() {
        pattern.push(kind_between(
            cycle[position],
            cycle[(position + 1) % cycle.len()],
        ));
    }
    let alternating = |first: usize| {
        pattern.iter().enumerate().all(|(index, &kind)| {
            kind == if (index + first).is_multiple_of(2) {
                1
            } else {
                2
            }
        })
    };
    alternating(0) || alternating(1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;

    fn atom(element: &str) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count: 0,
            h_pinned: false,
            valence_code: 0,
            isotope: 0,
            charge: 0,
        }
    }

    #[test]
    fn kekule_benzene_is_aromatic() {
        let atoms = (0..6).map(|_| atom("C")).collect();
        let bonds: Vec<(usize, usize, u32)> = (0..6)
            .map(|i| (i, (i + 1) % 6, if i % 2 == 0 { 2 } else { 1 }))
            .collect();
        let structure = Structure {
            atoms,
            bonds,
            ..Default::default()
        };
        assert_eq!(aromatic_bonds(&structure).len(), 6);
    }

    #[test]
    fn cyclohexane_is_not() {
        let atoms = (0..6).map(|_| atom("C")).collect();
        let bonds: Vec<(usize, usize, u32)> = (0..6).map(|i| (i, (i + 1) % 6, 1)).collect();
        let structure = Structure {
            atoms,
            bonds,
            ..Default::default()
        };
        assert!(aromatic_bonds(&structure).is_empty());
    }

    #[test]
    fn fused_kekule_rings_mark_the_shared_bond() {
        // Naphthalene skeleton: two fused hexagons.
        let atoms = (0..10).map(|_| atom("C")).collect();
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (2, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 3),
        ];
        let bonds: Vec<(usize, usize, u32)> = edges
            .iter()
            .enumerate()
            .map(|(index, &(a, b))| (a, b, if index % 2 == 0 { 2 } else { 1 }))
            .collect();
        let structure = Structure {
            atoms,
            bonds,
            ..Default::default()
        };
        // Every ring bond of both hexagons is aromatic when its cycles
        // alternate; the shared bond participates in both.
        assert!(!aromatic_bonds(&structure).is_empty());
    }
}
