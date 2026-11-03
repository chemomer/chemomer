//! The linear connection table: for each canonical rank in order, the
//! rank followed by its earlier-ranked neighbors in ascending order.
//! Ordering follows the reference's `CtFullCompare` base layer,
//! including its negative length rule.

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConnectionTable {
    entries: Vec<u32>,
    atoms: usize,
    /// The length of the atom-row prefix: the first `atom_prefix` nodes
    /// are real atoms, any further nodes are the appended t-group
    /// rows — the reference splits its compare there
    /// (`bSplitTautCompare`: layer 0, the atom rows, decides before
    /// the t-group rows are ever consulted).
    atom_prefix_len: usize,
}

impl ConnectionTable {
    /// Builds the table from final canonical ranks (a permutation of
    /// 1..n).
    pub fn build(ranks: &[u32], adjacency: &[Vec<usize>]) -> ConnectionTable {
        Self::build_with_atom_prefix(ranks, adjacency, adjacency.len())
    }

    /// [`ConnectionTable::build`] over a graph whose tail beyond
    /// `atom_prefix` nodes carries the t-group rows.
    pub fn build_with_atom_prefix(
        ranks: &[u32],
        adjacency: &[Vec<usize>],
        atom_prefix: usize,
    ) -> ConnectionTable {
        let atoms = adjacency.len();
        let mut entries = Vec::with_capacity(atoms * 3);
        let mut atom_prefix_len = 0usize;
        for rank in 1..=atoms as u32 {
            let node = ranks
                .iter()
                .position(|&r| r == rank)
                .expect("ranks are a permutation of 1..n");
            entries.push(rank);
            let mut earlier: Vec<u32> = adjacency[node]
                .iter()
                .map(|&b| ranks[b])
                .filter(|&r| r < rank)
                .collect();
            earlier.sort_unstable();
            entries.extend(earlier);
            if rank <= atom_prefix as u32 {
                atom_prefix_len = entries.len();
            }
        }
        ConnectionTable {
            entries,
            atoms,
            atom_prefix_len,
        }
    }

    pub fn entries(&self) -> &[u32] {
        &self.entries
    }

    /// The number of atoms the table numbers.
    pub fn atoms(&self) -> usize {
        self.atoms
    }

    /// The atom-row prefix: the compare's layer 0.
    pub fn atom_entries(&self) -> &[u32] {
        &self.entries[..self.atom_prefix_len]
    }

    /// `CtFullCompare` under the reference's `INCHI_CANON_MIN`
    /// semantics: the search keeps the GREATER result, and greater
    /// means the lexicographically SMALLER connection table. The
    /// reference compares the raw entries, negates the difference
    /// (`INCHI_CANON_MIN` in ichican2.c) and lets the caller replace
    /// its best terminal on a positive result — so a smaller table
    /// wins, and a longer table sorts above a shorter one (the length
    /// difference is negated on entry, then the whole difference is
    /// negated again).
    pub fn compare(&self, other: &ConnectionTable) -> std::cmp::Ordering {
        if self.entries.len() != other.entries.len() {
            return self.entries.len().cmp(&other.entries.len());
        }
        other.entries.cmp(&self.entries)
    }

    /// The layered `bSplitTautCompare` verdict: the atom rows (layer
    /// 0) decide first; the t-group rows compare only when the atoms
    /// tie. Same greater-means-smaller sign rules as [`Self::compare`].
    pub fn compare_layered(&self, other: &ConnectionTable) -> std::cmp::Ordering {
        if self.atom_entries() != other.atom_entries() {
            if self.atom_entries().len() != other.atom_entries().len() {
                return self.atom_entries().len().cmp(&other.atom_entries().len());
            }
            return other.atom_entries().cmp(self.atom_entries());
        }
        self.compare(other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(n: usize, bonds: &[(usize, usize)]) -> Vec<Vec<usize>> {
        let mut adjacency = vec![Vec::new(); n];
        for &(a, b) in bonds {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        adjacency
    }

    #[test]
    fn the_smaller_table_compares_as_greater() {
        // The isobutane skeleton (atom 1 is the hub) under two
        // numberings. Starting at a leaf gives the smaller table.
        let adjacency = graph(4, &[(0, 1), (1, 2), (1, 3)]);
        let leaf_first = ConnectionTable::build(&[2, 1, 3, 4], &adjacency);
        let hub_first = ConnectionTable::build(&[1, 2, 3, 4], &adjacency);
        assert_eq!(leaf_first.entries(), &[1, 2, 1, 3, 1, 4, 1]);
        assert_eq!(hub_first.entries(), &[1, 2, 1, 3, 2, 4, 2]);
        assert_eq!(leaf_first.compare(&hub_first), std::cmp::Ordering::Greater);
        assert_eq!(hub_first.compare(&leaf_first), std::cmp::Ordering::Less);
    }

    #[test]
    fn a_longer_table_compares_as_greater() {
        // The reference negates the length difference on entry and the
        // whole difference again at `done`, so the longer table carries
        // the positive (winning) sign.
        let short = ConnectionTable::build(&[1, 2], &graph(2, &[(0, 1)]));
        let long = ConnectionTable::build(&[1, 2, 3], &graph(3, &[(0, 1), (1, 2)]));
        assert_eq!(long.compare(&short), std::cmp::Ordering::Greater);
        assert_eq!(short.compare(&long), std::cmp::Ordering::Less);
    }
}
