//! Atom partitions: rank-per-atom plus the atoms in rank order.
//! Ported from the reference's `Partition` handling — `SetInitialRanks2`,
//! `DifferentiateRanks3` (equitable refinement), and `PartitionColorVertex`
//! (individualization). Ranks are position-based: a cell's rank is the
//! position just past its end.

#[derive(Clone)]
pub struct Partition {
    rank: Vec<u32>,
    at_number: Vec<usize>,
}

impl Partition {
    /// The invariant-ordered partition. Cells keep their atoms in
    /// atom order; a cell's rank is its end position (1-based).
    pub fn initial(invariants: &[(u16, u16)]) -> Partition {
        let n = invariants.len();
        let mut at_number: Vec<usize> = (0..n).collect();
        at_number.sort_by_key(|&i| (invariants[i], i));
        let mut rank = vec![0u32; n];
        let mut current = n as u32;
        if n > 0 {
            rank[at_number[n - 1]] = current;
        }
        for position in (1..n).rev() {
            if invariants[at_number[position - 1]] != invariants[at_number[position]] {
                current = position as u32;
            }
            rank[at_number[position - 1]] = current;
        }
        Partition { rank, at_number }
    }

    /// The invariant-ordered partition over a pair-valued invariant
    /// (the reference's auxiliary ranking input).
    pub fn initial_pair(invariants: &[(u32, u32)]) -> Partition {
        let n = invariants.len();
        let mut at_number: Vec<usize> = (0..n).collect();
        at_number.sort_by_key(|&i| (invariants[i], i));
        let mut rank = vec![0u32; n];
        let mut current = n as u32;
        if n > 0 {
            rank[at_number[n - 1]] = current;
        }
        for position in (1..n).rev() {
            if invariants[at_number[position - 1]] != invariants[at_number[position]] {
                current = position as u32;
            }
            rank[at_number[position - 1]] = current;
        }
        Partition { rank, at_number }
    }

    /// Equitable refinement: tied atoms are ordered by their sorted
    /// neighbor-rank lists, and cells split at each difference.
    /// Repeats until the partition stops splitting.
    pub fn refine(&mut self, adjacency: &[Vec<usize>]) {
        let n = self.at_number.len();
        loop {
            let neighbors: Vec<Vec<u32>> = adjacency
                .iter()
                .map(|ns| {
                    let mut ranks: Vec<u32> = ns.iter().map(|&b| self.rank[b]).collect();
                    ranks.sort_unstable();
                    ranks
                })
                .collect();
            let mut new_rank = vec![0u32; n];
            let mut position = 0usize;
            let mut expected = 1u32;
            let mut splits = 0usize;
            while position < n {
                let atom = self.at_number[position];
                let mut cell_end = self.rank[atom];
                if expected == cell_end {
                    new_rank[atom] = cell_end;
                    position += 1;
                    expected += 1;
                    continue;
                }
                let start = position;
                let old_end = cell_end;
                self.at_number[start..old_end as usize]
                    .sort_by(|a, b| neighbors[*a].cmp(&neighbors[*b]));
                let mut scan = old_end as usize - 1;
                new_rank[self.at_number[scan]] = old_end;
                while scan > start {
                    if neighbors[self.at_number[scan - 1]] != neighbors[self.at_number[scan]] {
                        cell_end = scan as u32;
                        splits += 1;
                    }
                    scan -= 1;
                    new_rank[self.at_number[scan]] = cell_end;
                }
                position = old_end as usize;
                expected = old_end + 1;
            }
            self.rank = new_rank;
            if splits == 0 {
                return;
            }
        }
    }

    /// Splits `v` off the front of its cell and re-refines
    /// (`PartitionColorVertex`).
    pub fn individualize(&mut self, adjacency: &[Vec<usize>], v: usize) {
        let slot = self
            .at_number
            .iter()
            .position(|&a| a == v)
            .expect("v is in the partition");
        let cell_rank = self.rank[v];
        let mut i = slot;
        let mut j = slot as isize - 1;
        while j >= 0 && self.rank[self.at_number[j as usize]] == cell_rank {
            self.at_number[i] = self.at_number[j as usize];
            i = j as usize;
            j -= 1;
        }
        let new_rank = if i > 0 {
            self.rank[self.at_number[i - 1]] + 1
        } else {
            1
        };
        self.at_number[i] = v;
        self.rank[v] = new_rank;
        self.refine(adjacency);
    }

    /// The raw tie-break pin of the reference's mapping machinery
    /// (`map_an_atom2` / `BreakAllTies`): assign `atom` the given rank
    /// value, re-sort the atom order over the (possibly stale) values,
    /// and run the full `DifferentiateRanks2` — atoms re-ranked from
    /// their (rank, sorted-neighbor-ranks) signature groups, group
    /// ranks the group's top position + 1, repeated until stable.
    /// Unlike [`Partition::refine`] this reassigns EVERY rank from the
    /// sorted order, so a pinned value only steers the sort.
    pub(crate) fn pin_and_refine(&mut self, atom: usize, rank: u32, adjacency: &[Vec<usize>]) {
        let n = self.rank.len();
        self.rank[atom] = rank;
        self.at_number.sort_by_key(|&a| (self.rank[a], a));
        loop {
            // The neighbor signatures: each atom's neighbor ranks,
            // ascending (SortNeighLists2's rank-ordered neighbor
            // lists).
            let neighbors: Vec<Vec<u32>> = adjacency
                .iter()
                .map(|ns| {
                    let mut ranks: Vec<u32> = ns.iter().map(|&b| self.rank[b]).collect();
                    ranks.sort_unstable();
                    ranks
                })
                .collect();
            // SetNewRanksFromNeighLists: atoms sorted by (rank,
            // signature, atom) — a total order — and each signature
            // group assigned its top position + 1.
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by(|&a, &b| {
                (self.rank[a], &neighbors[a])
                    .cmp(&(self.rank[b], &neighbors[b]))
                    .then(a.cmp(&b))
            });
            let mut new_rank = vec![0u32; n];
            let mut current = n as u32;
            if n > 0 {
                new_rank[order[n - 1]] = current;
            }
            for position in (1..n).rev() {
                if (
                    self.rank[order[position - 1]],
                    &neighbors[order[position - 1]],
                ) != (self.rank[order[position]], &neighbors[order[position]])
                {
                    current = position as u32;
                }
                new_rank[order[position - 1]] = current;
            }
            let stable = new_rank == self.rank;
            self.rank = new_rank;
            self.at_number = order;
            if stable {
                return;
            }
        }
    }

    pub fn is_discrete(&self) -> bool {
        self.first_nonsingleton_cell().is_none()
    }

    /// Position span of the first non-singleton cell, if any.
    pub fn first_nonsingleton_cell(&self) -> Option<(usize, usize)> {
        let n = self.at_number.len();
        self.first_nonsingleton_cell_within(n)
    }

    /// Position span of the first non-singleton cell among the first
    /// `bound` positions — the reference's `PartitionGetFirstCell`
    /// over `n` (atoms) while t-group nodes ride the tail: the search
    /// individualizes atoms only, and tied t-group cells never branch.
    pub fn first_nonsingleton_cell_within(&self, bound: usize) -> Option<(usize, usize)> {
        let n = self.at_number.len().min(bound);
        let mut position = 0usize;
        while position < n {
            let r = self.rank[self.at_number[position]];
            let mut end = position;
            while end < n && self.rank[self.at_number[end]] == r {
                end += 1;
            }
            if end - position > 1 {
                return Some((position, end));
            }
            position = end;
        }
        None
    }

    pub fn ranks(&self) -> &[u32] {
        &self.rank
    }

    /// The rank of one atom.
    pub fn rank_of(&self, atom: usize) -> u32 {
        self.rank[atom]
    }

    pub fn at_number(&self) -> &[usize] {
        &self.at_number
    }

    /// The canonical number of each atom, as given by the cell of the
    /// partition the atom sits in (class identity, not final ranks).
    pub fn cell_of_each_atom(&self) -> Vec<usize> {
        self.cell_count();
        let mut cells = vec![0usize; self.rank.len()];
        let mut position = 0usize;
        let n = self.at_number.len();
        let mut class = 0usize;
        while position < n {
            let r = self.rank[self.at_number[position]];
            let mut end = position;
            while end < n && self.rank[self.at_number[end]] == r {
                end += 1;
            }
            for atom in &self.at_number[position..end] {
                cells[*atom] = class;
            }
            class += 1;
            position = end;
        }
        cells
    }

    pub fn cell_count(&self) -> usize {
        let n = self.at_number.len();
        let mut cells = 0usize;
        let mut position = 0usize;
        while position < n {
            let r = self.rank[self.at_number[position]];
            while position < n && self.rank[self.at_number[position]] == r {
                position += 1;
            }
            cells += 1;
        }
        cells
    }
}
