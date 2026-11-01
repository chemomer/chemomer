//! The AuxInfo emitter: the auxiliary information the reference
//! writes beside every identifier — the canonical numbering of the
//! input atoms (`/N:`), the constitutional equivalence classes
//! (`/E:`), the normalized charges (`/CRV:`), and the reversibility
//! echo of the ORIGINAL submitted structure (`/rA`, `/rB`, `/rC`).
//! The canonical fields read the normalized structure the layers
//! speak; the echo fields are a pure function of the parse, explicit
//! hydrogens and all (the reference's `WriteOrigAtoms` /
//! `WriteOrigBonds` / `WriteOrigCoord`, ichiprt1.c).

use crate::error::Result;
use crate::graph::Structure;
use crate::layers::hydrogen::canonical_ranks_for;
use crate::layers::standardization;
use crate::options::Options;

/// The AuxInfo line for a submitted structure (`original`, exactly as
/// parsed) and its normalized form (`working`, the one the layers
/// read). `proton_shift` is the net protons the (DE)PROTONATION
/// standardization moved — the `/p` layer's pinned-path value. The
/// options ride along: `SNon` clears the inverted-stereo block (the
/// oracle's SNon recordings drop `/it:` while keeping the
/// stereo-informed numbering).
pub fn aux_info(
    original: &Structure,
    working: &Structure,
    proton_shift: i32,
    options: &Options,
) -> Result<String> {
    standardization::ensure_typable_charges(working)?;
    let counts = standardization::adjusted_hydrogens(working);
    let ranks = canonical_ranks_for(working);
    let groups = crate::canonical::tautomer::mobile_groups(working, &counts, &ranks);

    // The normalization-type digit (the reference's io.bTautomeric):
    // +1 a tautomeric representation exists — a registered mobile
    // group holding real hydrogen, or the normalization moved protons
    // (the `/p` layer — the reference's FLAG_NORM_CONSIDER_TAUT
    // marker), or an exchangeable isotopic H sits on a heteroatom
    // (the reference's isotopic-endpoint branch). +1 a merged salt
    // group (a real-H group that also carries charge — the (Hn±)
    // super-group). +4 a hard (DE)PROTONATION moved a proton. The
    // salt and hard bits count only when the +1 did (the reference
    // adds them under `if (io.bTautomeric)`); charge-only c-point
    // registrations are not t-groups and count nothing.
    let real_h = |group: &crate::canonical::tautomer::MobileGroup| group.count > group.charge;
    let tautomeric = groups.iter().any(&real_h);
    let salt_merged = groups.iter().any(|group| real_h(group) && group.charge > 0);
    let has_p = proton_shift != 0 || standardization::proton_layer(working).is_some();
    let iso_exchangeable = exchangeable_isotopic_h(working, &counts);
    let mut norm = i32::from(tautomeric || has_p || iso_exchangeable);
    if norm > 0 {
        norm += i32::from(salt_merged);
        norm += 4 * i32::from(working.protons_hard_moved);
    }

    let (adjacency, invariants) = crate::canonical::graph_of(working);
    // The aux speaks the numbering and symmetry of the representation
    // the reference fills AuxInfo from: for a tautomeric structure the
    // NON-TAUT one (mobile H removed from the endpoints, no t-group
    // rows), otherwise the structure as the layers read it. Chiral
    // structures number by the STEREO-informed canonicalization (the
    // reference's nCanonOrdStereo): the terminal compare extends with
    // the stereo rows, and when the inverted-parity pass wins its
    // numbering is adopted (nCompInv2Abs = -1) — the verdict also
    // decides the `/it:` inverted-stereo block.
    let ranks_aux;
    let symmetry;
    let mut inverted_block = String::new();
    // Chiral structures number by the STEREO-informed canonicalization
    // (the reference's nCanonOrdStereo): the terminal compare extends
    // with the stereo rows, and when the inverted-parity pass wins its
    // numbering is adopted (nCompInv2Abs = -1) — the verdict also
    // decides the `/it:` inverted-stereo block. Tautomeric structures
    // rank the NON-TAUT representation (mobile H stripped, no t-group
    // rows); the helper mirrors that input assembly.
    if let Some(verdict) =
        crate::canonical::stereo::informed_verdict(working, &counts, &ranks, tautomeric)
    {
        ranks_aux = verdict.ranks;
        symmetry = verdict.symmetry;
        // SNon withholds the stereo layers and the inverted-stereo
        // fields alike (ichiparm.c:248–253): the block never prints
        // under that mode.
        if verdict.sign != 0 && !options.no_stereo {
            // The inverted-stereo block: the marker `im` ("main =
            // Inv(main)") when the inverted layer is the pairwise
            // parity-flip of the main over the same centers; the
            // single-component corpus never prints the list form.
            let flip_pairs = verdict.main_pairs.len() == verdict.inv_pairs.len()
                && verdict
                    .main_pairs
                    .iter()
                    .zip(&verdict.inv_pairs)
                    .all(|(main, inv)| {
                        main.0 == inv.0
                            && (main.1 + inv.1 == 3
                                || (!matches!(main.1, 1 | 2) && main.1 == inv.1))
                    });
            if flip_pairs {
                inverted_block = "/it:im".to_string();
            }
        }
    } else if tautomeric {
        let mask = crate::canonical::tautomer::endpoint_mask(working, &counts);
        let notaut = crate::canonical::tautomer::without_mobile_h(&counts, &mask);
        let (ranks, classes) =
            crate::canonical::canonical_ranks_and_symmetry(&adjacency, &invariants, &notaut, &[]);
        ranks_aux = ranks;
        symmetry = classes;
    } else {
        let groups_rows = crate::canonical::group_membership(working, &counts);
        let (ranks, classes) = crate::canonical::canonical_ranks_and_symmetry(
            &adjacency,
            &invariants,
            &counts,
            &groups_rows,
        );
        ranks_aux = ranks;
        symmetry = classes;
    }

    let mut out = format!("AuxInfo=1/{norm}/N:");
    // `/N:` — the ORIGINAL atom numbers in canonical order (the
    // reference's nOrigAtNosInCanonOrd reads orig_at_number; the
    // normalized structure's indices shift when explicit hydrogen
    // atoms folded away, so the echo structure supplies the map).
    let original_index: Vec<usize> = original
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.element != "H")
        .map(|(index, _)| index)
        .collect();
    let mut order: Vec<usize> = (0..working.atoms.len()).collect();
    order.sort_by_key(|&atom| ranks_aux[atom]);
    let numbering: Vec<String> = order
        .iter()
        .map(|&atom| (original_index[atom] + 1).to_string())
        .collect();
    out.push_str(&numbering.join(","));

    // `/E:` — the constitutional equivalence classes: members by
    // canonical number, non-singleton classes only, classes ordered
    // by their minimum.
    let mut classes: std::collections::BTreeMap<u32, Vec<u32>> = Default::default();
    for (atom, class) in symmetry.iter().enumerate() {
        classes.entry(*class).or_default().push(ranks_aux[atom]);
    }
    let mut equivalents: Vec<Vec<u32>> = classes
        .into_values()
        .filter(|members| members.len() > 1)
        .map(|mut members| {
            members.sort_unstable();
            members
        })
        .collect();
    equivalents.sort();
    if !equivalents.is_empty() {
        out.push_str("/E:");
        let text: Vec<String> = equivalents
            .iter()
            .map(|members| {
                format!(
                    "({})",
                    members
                        .iter()
                        .map(|m| m.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect();
        out.push_str(&text.join(""));
    }

    // `/gE:` — equivalence of the registered mobile groups: groups
    // whose members fill the same symmetry classes are one class,
    // labelled by the groups' 1-based positions in the identifier.
    // Only t-group-born registrations (real mobile hydrogen) speak —
    // charge-only c-points register no identifier rows.
    let taut_groups: Vec<usize> = groups
        .iter()
        .enumerate()
        .filter(|(_, group)| group.count > group.charge)
        .map(|(index, _)| index)
        .collect();
    if taut_groups.len() > 1 {
        let atom_of_rank = |rank: u32| ranks_aux.iter().position(|&r| r == rank);
        let mut signatures: std::collections::BTreeMap<Vec<u32>, Vec<usize>> = Default::default();
        for index in &taut_groups {
            let group = &groups[*index];
            let mut signature: Vec<u32> = group
                .ranks
                .iter()
                .filter_map(|&rank| atom_of_rank(rank))
                .map(|atom| symmetry[atom])
                .collect();
            signature.sort_unstable();
            signatures.entry(signature).or_default().push(*index);
        }
        let mut equivalents: Vec<Vec<u32>> = signatures
            .into_values()
            .filter(|members| members.len() > 1)
            .map(|mut members| {
                members.sort_unstable();
                members.into_iter().map(|i| i as u32 + 1).collect()
            })
            .collect();
        equivalents.sort();
        if !equivalents.is_empty() {
            out.push_str("/gE:");
            let text: Vec<String> = equivalents
                .iter()
                .map(|members| {
                    format!(
                        "({})",
                        members
                            .iter()
                            .map(|m| m.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    )
                })
                .collect();
            out.push_str(&text.join(""));
        }
    }

    // `/it:` — the inverted-stereo block from the informed search's
    // verdict (empty when the absolute pass won or the rows tied).
    out.push_str(&inverted_block);

    // `/I:/E:m` — the isotopic block: fixed isotopic atoms present and
    // constitutional equivalence spoken — the isotopic numbering adds
    // nothing (the empty /I:) and the isotopic equivalence equals the
    // plain classes (the `m` marker, "same as the main layer").
    let fixed_isotopes = working.atoms.iter().any(|atom| atom.isotope != 0)
        || working
            .isotopic_h
            .iter()
            .any(|(carrier, _)| working.atoms[*carrier].element == "C");
    if fixed_isotopes && out.contains("/E:") {
        out.push_str("/I:/E:m");
    }

    // `/CRV:` — normalized charges and unusual valences, ascending
    // canonical number: `<canon>[±charge][.valence]`. The charges are
    // the reference's norm_at charges: our salt normalization zeroes a
    // bonded fixed-cation/anion pair the reference's norm_at KEEPS
    // (its BNS rides the valence internally), so those two atoms read
    // their submitted charges. The unusual valence reads the
    // BNS-resolved integer bond valence, not the 1.5 aromatic sums.
    let bond_sums = working.bond_order_sums();
    let degrees = working.degrees();
    let kept = kept_ion_pair_charges(original);
    let original_sums = original.bond_order_sums();
    let mut flagged: Vec<(u32, i32, i32)> = working
        .atoms
        .iter()
        .enumerate()
        .map(|(atom, a)| {
            let el = crate::canonical::protons::el_number(&a.element).unwrap_or(0);
            // A kept ion-pair atom reads its SUBMITTED charge and the
            // SUBMITTED bond valence — the reference's norm_at keeps
            // both (its BNS rides the valence without rewriting).
            let kept = kept
                .as_ref()
                .and_then(|map| map.get(&original_index[atom]).copied())
                .unwrap_or(0);
            let (charge, bond_valence) = if a.charge == 0 && kept != 0 {
                (kept, original_sums[original_index[atom]].round() as i32)
            } else {
                (a.charge, bond_sums[atom].round() as i32)
            };
            let unusual = crate::canonical::protons::unusual_el_valence(
                el,
                charge,
                bond_valence,
                counts[atom],
                degrees[atom],
            );
            (ranks_aux[atom], charge, unusual)
        })
        .filter(|&(_, charge, unusual)| charge != 0 || unusual != 0)
        .collect();
    flagged.sort_unstable();
    if !flagged.is_empty() {
        out.push_str("/CRV:");
        let text: Vec<String> = flagged
            .iter()
            .map(|(rank, charge, unusual)| {
                let mut entry = rank.to_string();
                if *charge < 0 {
                    entry.push('-');
                } else if *charge > 0 {
                    entry.push('+');
                }
                if *charge != 0 {
                    entry.push_str(&charge.abs().to_string());
                }
                if *unusual != 0 {
                    entry.push('.');
                    entry.push_str(&unusual.to_string());
                }
                entry
            })
            .collect();
        out.push_str(&text.join(","));
    }

    out.push_str(&reversibility_block(original));
    Ok(out)
}

/// Exchangeable isotopic hydrogen on a heteroatom: the structure
/// carries isotopes AND a non-carbon atom holds hydrogen — the
/// reference's isotopic-endpoint branch (its salt-candidate scan
/// counts the heteroatom and sets `nNumIsotopicEndpoints` above one;
/// isotopes riding carbon — toluene-d3's methyl — never qualify).
fn exchangeable_isotopic_h(structure: &Structure, counts: &[u32]) -> bool {
    let isotopic =
        !structure.isotopic_h.is_empty() || structure.atoms.iter().any(|atom| atom.isotope != 0);
    isotopic
        && structure
            .atoms
            .iter()
            .enumerate()
            .any(|(index, atom)| atom.element != "C" && atom.element != "H" && counts[index] > 0)
}

/// The charges the salt normalization zeroes but the reference's
/// norm_at KEEPS: a directly-bonded fixed-cation/anion pair (the
/// terminal ion pair). CRV reads these submitted charges for the two
/// consumed atoms — and nothing else: everywhere else the normalized
/// structure's charges stand (the hard (DE)PROTONATION sites among
/// them, whose CRV values are the MOVED ones).
fn kept_ion_pair_charges(original: &Structure) -> Option<std::collections::HashMap<usize, i32>> {
    let sums = original.bond_order_sums();
    let degrees = original.degrees();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); original.atoms.len()];
    for &(from, to, _) in &original.bonds {
        neighbors[from].push(to);
        neighbors[to].push(from);
    }
    let mut consumed = vec![false; original.atoms.len()];
    let mut map = std::collections::HashMap::new();
    let anions: Vec<usize> = original
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| {
            atom.charge < 0 && matches!(atom.element.as_str(), "O" | "S" | "Se" | "Te" | "N")
        })
        .map(|(index, _)| index)
        .collect();
    for anion in anions {
        if consumed[anion] {
            continue;
        }
        let found = neighbors[anion].iter().copied().find(|&site| {
            original.atoms[site].charge > 0
                && !consumed[site]
                && degrees[site] == 4
                && !standardization::cation_is_proton_related(
                    &original.atoms[site],
                    degrees[site],
                    sums[site],
                )
        });
        if let Some(site) = found {
            map.insert(anion, original.atoms[anion].charge);
            map.insert(site, original.atoms[site].charge);
            consumed[anion] = true;
            consumed[site] = true;
        }
    }
    (!map.is_empty()).then_some(map)
}

/// The `/rA` `/rB` `/rC` echo of the submitted structure, in the
/// input's own atom order, explicit hydrogens included.
fn reversibility_block(original: &Structure) -> String {
    let n = original.atoms.len();
    // The atom-block prefix letter records the chiral flag; the reader
    // does not carry it yet, so the structure reads non-chiral.
    let mut out = format!("/rA:{n}n");
    for atom in &original.atoms {
        out.push_str(&atom.element);
        match atom.charge {
            0 => {}
            negative if negative < 0 => {
                out.push('-');
                if negative < -1 {
                    out.push_str(&negative.abs().to_string());
                }
            }
            positive => {
                out.push('+');
                if positive > 1 {
                    out.push_str(&positive.to_string());
                }
            }
        }
        let dotted = atom.element.chars().count() == 1;
        if atom.isotope != 0 {
            if dotted {
                out.push('.');
            }
            out.push('i');
            out.push_str(&atom.isotope.to_string());
        }
    }

    // `/rB` — one `;`-terminated entry per atom 2..n: the bonds to
    // lower-numbered neighbors, ascending, each a type char and the
    // neighbor's 1-based number.
    out.push_str("/rB:");
    for j in 1..n {
        let mut entries: Vec<(usize, char)> = Vec::new();
        for (index, &(from, to, kind)) in original.bonds.iter().enumerate() {
            let (_end, other) = if from == j {
                (j, to)
            } else if to == j {
                (j, from)
            } else {
                continue;
            };
            if other >= j {
                continue;
            }
            let stereo = original.stereo_bonds.get(&index).copied().unwrap_or(0);
            // The wedge perspective reads the bond line's FIRST atom
            // (the reference stores +stereo there, -stereo at the
            // other end): `p`/`n` when this end is first, `P`/`N`
            // when the wedge points from the other end.
            let this_end_first = from == j;
            entries.push((other, bond_char(kind, stereo, this_end_first)));
        }
        entries.sort_unstable();
        for (neighbor, character) in entries {
            out.push(character);
            out.push_str(&(neighbor + 1).to_string());
        }
        out.push(';');
    }

    // `/rC` — one `;`-terminated coordinate per atom, trailing zeros
    // trimmed, an all-zero coordinate an empty entry.
    out.push_str("/rC:");
    for (x, y) in &original.positions {
        if *x == 0.0 && *y == 0.0 {
            out.push(';');
        } else {
            out.push_str(&format!("{},{},0;", clean_coord(*x), clean_coord(*y)));
        }
    }
    out
}

/// The bond-type character (the reference's WriteOrigBonds switch,
/// `at[j].bond_stereo` read at the current endpoint with its sign):
/// an up wedge is `p` at the bond-line-first atom, `P` at the other
/// (the wedge points from the neighbor); a dash is `n`/`N` the same
/// way; either is `v`; plain bonds `s`/`d`/`t`/`a`.
fn bond_char(kind: u32, stereo: u8, this_end_first: bool) -> char {
    match kind {
        2 => 'd',
        3 => 't',
        4 => 'a',
        _ => match stereo {
            1 => {
                if this_end_first {
                    'p'
                } else {
                    'P'
                }
            }
            6 => {
                if this_end_first {
                    'n'
                } else {
                    'N'
                }
            }
            4 => 'v',
            _ => 's',
        },
    }
}

/// The reference's `CleanOrigCoord`: fixed four-decimal coordinates
/// with trailing zeros (and the trailing point) trimmed; a leading
/// zero after the sign goes too — the oracle writes `-.7`.
fn clean_coord(value: f64) -> String {
    let mut text = format!("{value:.4}");
    while text.ends_with('0') {
        text.pop();
    }
    if text.ends_with('.') {
        text.pop();
    }
    if let Some(rest) = text.strip_prefix("-0.") {
        return format!("-.{rest}");
    }
    if let Some(rest) = text.strip_prefix("0.") {
        return rest.to_string();
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::AtomRecord;

    fn atom(element: &str) -> AtomRecord {
        AtomRecord {
            element: element.to_string(),
            h_count: 0,
            charge: 0,
            h_pinned: false,
            isotope: 0,
        }
    }

    const WATER_MOL: &str = "water\n  AsciiChem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   12.0000   12.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    #[test]
    fn water_reads_aux_info_with_empty_equivalence() {
        let original = crate::molfile::parse(WATER_MOL).expect("parses");
        let mut working = original.clone();
        working.fold_explicit_hydrogens();
        let aux = aux_info(&original, &working, 0, &Options::default()).expect("aux");
        assert_eq!(aux, "AuxInfo=1/0/N:1/rA:1nO/rB:/rC:12,12,0;");
    }

    #[test]
    fn coord_trimming_matches_the_reference() {
        assert_eq!(clean_coord(12.0), "12");
        assert_eq!(clean_coord(1.4), "1.4");
        assert_eq!(clean_coord(-0.7), "-.7");
        assert_eq!(clean_coord(0.0), "0");
    }

    #[test]
    fn bond_chars_render_the_wedge_perspective() {
        assert_eq!(bond_char(1, 0, true), 's');
        assert_eq!(bond_char(2, 0, true), 'd');
        assert_eq!(bond_char(4, 0, true), 'a');
        // Up: p at the bond-line-first atom, P from the other end.
        assert_eq!(bond_char(1, 1, true), 'p');
        assert_eq!(bond_char(1, 1, false), 'P');
        // Dash: n at the first atom, N from the other end.
        assert_eq!(bond_char(1, 6, true), 'n');
        assert_eq!(bond_char(1, 6, false), 'N');
    }

    #[test]
    fn acetate_reads_its_mobile_group_and_charge() {
        // The recorded oracle line for the acetate battery fixture:
        // the mobile group turns the normalization digit to 1, the
        // anionic oxygen shows in /CRV: by its canonical number.
        let o_anion = AtomRecord {
            element: "O".to_string(),
            h_count: 0,
            charge: -1,
            h_pinned: false,
            isotope: 0,
        };
        let mut working = Structure {
            atoms: vec![atom("C"), atom("C"), atom("O"), o_anion],
            bonds: vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
            ..Default::default()
        };
        working.fold_explicit_hydrogens();
        let original = working.clone();
        let aux = aux_info(&original, &working, 0, &Options::default()).expect("aux");
        assert!(aux.starts_with("AuxInfo=1/1/N:"), "{aux}");
        assert!(aux.contains("/CRV:"), "{aux}");
        assert!(aux.contains("(3,4)"), "{aux}");
    }
}
