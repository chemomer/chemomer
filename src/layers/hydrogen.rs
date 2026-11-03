//! The `/h` layer: fixed attached hydrogens per canonical rank
//! (counts ascending, ranks ascending within a count, consecutive
//! runs as ranges — the oracle's layout, e.g. `h2H,1H3`), followed by
//! the mobile-H groups (`(H,3,4)`). Counts are the ±1-standardized
//! ones; an anionic O's H rides its group.

use crate::error::Result;
use crate::graph::Structure;
use crate::layers::standardization;
use crate::layers::{order, LayerEmitter};

pub struct HydrogenLayer;

/// The negative charge riding /h groups that hold real mobile
/// hydrogen, summed — the /q layer subtracts it from the fixed
/// accounting. Charge-only groups (no mobile H) stay fixed: the
/// reference never registered them as t-groups.
pub fn mobile_charges(structure: &Structure) -> i32 {
    let counts = standardization::adjusted_hydrogens(structure);
    let ranks = canonical_ranks_for(structure);
    crate::canonical::tautomer::mobile_groups(structure, &counts, &ranks)
        .iter()
        .filter(|group| group.count > group.charge)
        .map(|group| -(group.charge as i32))
        .sum()
}

/// The canonical ranks the /h and /i layers speak: the heavy-atom
/// numbering under the standardized hydrogen counts (explicit H
/// atoms are invisible vertices — see `graph_of`).
pub fn canonical_ranks_for(structure: &Structure) -> Vec<u32> {
    let (adjacency, invariants) = crate::canonical::graph_of(structure);
    let counts = standardization::adjusted_hydrogens(structure);
    let mask = crate::canonical::tautomer::endpoint_mask(structure, &counts);
    let groups = crate::canonical::main_canon_groups(structure, &counts);
    crate::canonical::canonical_ranks(
        &adjacency,
        &invariants,
        &crate::canonical::tautomer::without_mobile_h(&counts, &mask),
        &groups,
    )
}

impl LayerEmitter for HydrogenLayer {
    fn name(&self) -> &'static str {
        "hydrogen"
    }

    fn order(&self) -> u8 {
        order::H
    }

    fn segment_kind(&self) -> Option<char> {
        Some('h')
    }

    fn emit(&self, structure: &Structure) -> Result<Option<String>> {
        standardization::ensure_typable_charges(structure)?;
        let counts = standardization::adjusted_hydrogens(structure);
        let ranks = canonical_ranks_for(structure);
        let groups = crate::canonical::tautomer::mobile_groups(structure, &counts, &ranks);
        let endpoints = crate::canonical::tautomer::endpoint_atoms(&ranks, &groups);
        let mut by_rank = vec![0u32; counts.len()];
        for (atom, count) in counts.iter().enumerate() {
            if !endpoints.contains(&atom) {
                by_rank[ranks[atom] as usize - 1] = *count;
            }
        }
        let mut segment = hydrogen_segment(&by_rank).unwrap_or_default();
        let groups_text = groups
            .iter()
            // A charge-only group (no mobile H) records nothing in
            // /h — the oracle leaves such carboxylates unlisted, the
            // charges folding into the fixed accounting.
            .filter(|group| group.count > group.charge)
            .map(|group| {
                let members: Vec<String> = group.ranks.iter().map(|r| r.to_string()).collect();
                // A charged group writes its mobile hydrogen (the
                // mobile units minus the charges) then the negative
                // charge: (H4-,...), (H2-2,...) — the reference's
                // mobile-group notation.
                let hydrogen = if group.charge > 0 {
                    group.count - group.charge
                } else {
                    group.count
                };
                let hydrogen = if hydrogen > 1 {
                    hydrogen.to_string()
                } else {
                    String::new()
                };
                let charge = match group.charge {
                    0 => String::new(),
                    1 => "-".to_string(),
                    charges => format!("-{charges}"),
                };
                format!("(H{}{},{})", hydrogen, charge, members.join(","))
            })
            .collect::<Vec<_>>()
            .join("");
        if !groups_text.is_empty() && !segment.is_empty() {
            segment.push(',');
        }
        segment.push_str(&groups_text);
        if segment.is_empty() {
            return Ok(None);
        }
        // The segment is the full InChI text: the letter rides along,
        // exactly like `p+1` and `c1-2...`.
        Ok(Some(format!("h{segment}")))
    }
}

fn hydrogen_segment(counts: &[u32]) -> Option<String> {
    let largest = counts.iter().copied().max()?;
    if largest == 0 {
        return None;
    }
    let mut segment = String::new();
    for count in 1..=largest {
        // by_rank is indexed by rank - 1; the segment speaks ranks.
        let ranks: Vec<u32> = (0..counts.len())
            .filter(|&atom| counts[atom] == count)
            .map(|atom| atom as u32 + 1)
            .collect();
        if ranks.is_empty() {
            continue;
        }
        if !segment.is_empty() {
            segment.push(',');
        }
        push_runs(&ranks, &mut segment);
        segment.push('H');
        if count > 1 {
            segment.push_str(&count.to_string());
        }
    }
    Some(segment)
}

fn push_runs(ranks: &[u32], out: &mut String) {
    let mut index = 0;
    while index < ranks.len() {
        let start = ranks[index];
        let mut end = start;
        while index + 1 < ranks.len() && ranks[index + 1] == end + 1 {
            index += 1;
            end = ranks[index];
        }
        if start == end {
            out.push_str(&start.to_string());
        } else {
            out.push_str(&format!("{start}-{end}"));
        }
        index += 1;
        if index < ranks.len() {
            out.push(',');
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;
    use crate::layers::LayerEmitter;

    fn atom(element: &str, h_count: u32, charge: i32) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count,
            h_pinned: false,
            valence_code: 0,
            isotope: 0,
            charge,
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
    fn runs_and_counts_format_like_the_oracle() {
        // by_rank inputs (rank r at index r-1): benzene is h1-6H;
        // acetaldehyde is h2H,1H3 (counts ascending).
        assert_eq!(hydrogen_segment(&[1, 1, 1, 1, 1, 1]).unwrap(), "1-6H");
        assert_eq!(hydrogen_segment(&[3, 1]).unwrap(), "2H,1H3");
        // Multiple runs in one count group (by_rank input).
        assert_eq!(
            hydrogen_segment(&[0, 1, 1, 0, 0, 0, 1, 1, 0, 1]).unwrap(),
            "2-3,7-8,10H"
        );
        assert!(hydrogen_segment(&[0, 0]).is_none());
    }

    #[test]
    fn acetate_carries_its_mobile_hydrogen_in_a_group() {
        // Oracle: InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)/p-1 — the
        // standardized H rides the (H,3,4) group, not a fixed entry.
        let o_anion = atom("O", 0, -1);
        let s = structure(
            vec![atom("C", 0, 0), atom("C", 0, 0), atom("O", 0, 0), o_anion],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        assert_eq!(
            HydrogenLayer.emit(&s).unwrap().as_deref(),
            Some("h1H3,(H,3,4)")
        );
    }

    #[test]
    fn hydrogen_counts_break_connectivity_ties_like_the_oracle() {
        // The 2-methylbutadiene skeleton: the methyl and the two CH2s
        // are connectivity-ambiguous; the terminal compare's hydrogen
        // layer settles them as the oracle records (h4H,1-2H2,3H3 —
        // the CH2s at ranks 1-2, the methyl at rank 3, the hub's CH
        // at rank 4).
        let s = structure(
            vec![
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
            ],
            vec![(0, 1, 1), (1, 2, 2), (3, 4, 2), (1, 3, 1)],
        );
        assert_eq!(
            HydrogenLayer.emit(&s).unwrap().as_deref(),
            Some("h4H,1-2H2,3H3")
        );
    }

    #[test]
    fn ammonium_loses_one_hydrogen() {
        let s = structure(vec![atom("N", 4, 1)], vec![]);
        assert_eq!(HydrogenLayer.emit(&s).unwrap().as_deref(), Some("h1H3"));
    }

    #[test]
    fn quaternary_nitrogen_keeps_zero() {
        let s = structure(
            vec![
                atom("N", 0, 1),
                atom("C", 0, 0),
                atom("C", 0, 0),
                atom("C", 0, 0),
            ],
            vec![(0, 1, 1), (0, 2, 1), (0, 3, 1)],
        );
        let counts = standardization::adjusted_hydrogens(&s);
        assert_eq!(counts[0], 0);
        assert_eq!(counts[1..], vec![3, 3, 3]);
    }
}
