//! The linear connection table: for each canonical rank in order, the
//! rank followed by its earlier-ranked neighbors in ascending order.
//! Ordering follows the reference's `CtFullCompare` base layer,
//! including its negative length rule.

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConnectionTable {
    entries: Vec<u32>,
    /// The hydrogen rows (the reference's `NumH`): empty for the plain
    /// tables, otherwise the non-tautomeric hydrogen count per atom in
    /// rank order followed by each group's mobile hydrogen — the
    /// compare's layer 1 and its layer 2 tail.
    numh: Vec<u32>,
    atoms: usize,
    /// The length of the atom-row prefix: the first `atom_prefix` nodes
    /// are real atoms, any further nodes are the appended t-group
    /// rows — the reference splits its compare there
    /// (`bSplitTautCompare`: layer 0, the atom rows, decides before
    /// the t-group rows are ever consulted).
    atom_prefix_len: usize,
    /// The node count the atom rows number (the atom-only adjacency's
    /// length); the hydrogen prefix slices at it.
    atom_rows: usize,
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
            numh: Vec::new(),
            atoms,
            atom_prefix_len,
            atom_rows: atoms,
        }
    }

    /// The TAUT walk's terminal table, the reference's layered CT
    /// (`CtFullCompare` under `bSplitTautCompare`): layer 0 the atom
    /// rows WITHOUT the t-group edges (the members live in layer 2's
    /// group rows), layer 1 the non-tautomeric hydrogen count per atom
    /// in rank order, layer 2 the group rows (rank + earlier member
    /// ranks) followed by each group's mobile hydrogen. `ranks` covers
    /// every node of `augmented`; `atom_adjacency` numbers the atoms
    /// alone.
    pub fn build_taut(
        ranks: &[u32],
        atom_adjacency: &[Vec<usize>],
        augmented: &[Vec<usize>],
        notaut: &[u32],
        group_h: &[u32],
    ) -> ConnectionTable {
        let atom_count = atom_adjacency.len();
        let total = augmented.len();
        let node_at = |rank: usize| -> usize {
            ranks
                .iter()
                .position(|&r| r == rank as u32)
                .expect("ranks are a permutation of 1..n")
        };
        let mut entries = Vec::new();
        let mut numh = Vec::new();
        // Layer 0: the atom rows over the atom-only adjacency.
        for rank in 1..=atom_count {
            let node = node_at(rank);
            entries.push(rank as u32);
            let mut earlier: Vec<u32> = atom_adjacency[node]
                .iter()
                .map(|&b| ranks[b])
                .filter(|&r| r < rank as u32)
                .collect();
            earlier.sort_unstable();
            entries.extend(earlier);
        }
        let atom_prefix_len = entries.len();
        // Layer 1: the non-tautomeric hydrogen, rank order.
        for rank in 1..=atom_count {
            numh.push(notaut[node_at(rank)]);
        }
        // Layer 2: the group rows over the augmented graph, then each
        // group's mobile hydrogen.
        for rank in atom_count + 1..=total {
            let node = node_at(rank);
            entries.push(rank as u32);
            let mut earlier: Vec<u32> = augmented[node]
                .iter()
                .map(|&b| ranks[b])
                .filter(|&r| r < rank as u32)
                .collect();
            earlier.sort_unstable();
            entries.extend(earlier);
        }
        numh.extend_from_slice(group_h);
        ConnectionTable {
            entries,
            numh,
            atoms: total,
            atom_prefix_len,
            atom_rows: atom_count,
        }
    }

    pub fn entries(&self) -> &[u32] {
        &self.entries
    }

    /// The hydrogen-row length (0 for the plain tables).
    pub fn numh_len(&self) -> usize {
        self.numh.len()
    }

    /// The full-layer identity: every compare layer level. Table-equal
    /// AND hydrogen-equal terminals are automorphic images (the
    /// reference's `CtFullCompare` reading zero — the isomorphism
    /// witness whose transposition joins theta); terminals equal on
    /// the connection rows alone need not be.
    pub fn identical(&self, other: &ConnectionTable) -> bool {
        self.entries == other.entries && self.numh == other.numh
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

    /// The layered `bSplitTautCompare` verdict (`CtFullCompare`): the
    /// whole-CT length first, then the atom rows (layer 0), the
    /// non-tautomeric hydrogen rows (layer 1), the group rows and the
    /// groups' hydrogen (layer 2). Plain tables reduce to the atom
    /// rows alone. Same greater-means-smaller sign rules as
    /// [`Self::compare`].
    pub fn compare_layered(&self, other: &ConnectionTable) -> std::cmp::Ordering {
        // The whole-CT length (endCt covers the atom rows and the group
        // rows; the hydrogen layers ride apart).
        if self.entries.len() != other.entries.len() {
            return self.entries.len().cmp(&other.entries.len());
        }
        // Layer 0: the atom rows.
        let (mine, theirs) = (self.atom_entries(), other.atom_entries());
        if mine != theirs {
            if mine.len() != theirs.len() {
                return mine.len().cmp(&theirs.len());
            }
            return theirs.cmp(mine);
        }
        if self.numh.is_empty() {
            // Plain tables: the atom rows are the whole table.
            return std::cmp::Ordering::Equal;
        }
        // Layer 1: the non-tautomeric hydrogen rows (the equal row
        // count is a graph invariant; the negated-length rule keeps the
        // longer row above).
        let (mine, theirs) = (&self.numh[..self.atom_rows], &other.numh[..other.atom_rows]);
        if mine != theirs {
            if mine.len() != theirs.len() {
                return mine.len().cmp(&theirs.len());
            }
            return theirs.cmp(mine);
        }
        // Layer 2: the group rows, then the groups' mobile hydrogen.
        let (mine, theirs) = (
            &self.entries[self.atom_prefix_len..],
            &other.entries[other.atom_prefix_len..],
        );
        if mine != theirs {
            if mine.len() != theirs.len() {
                return mine.len().cmp(&theirs.len());
            }
            return theirs.cmp(mine);
        }
        let (mine, theirs) = (&self.numh[self.atom_rows..], &other.numh[other.atom_rows..]);
        if mine != theirs {
            if mine.len() != theirs.len() {
                return mine.len().cmp(&theirs.len());
            }
            return theirs.cmp(mine);
        }
        std::cmp::Ordering::Equal
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
