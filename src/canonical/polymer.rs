//! The polymer treatment: the reference's validation, CRU caps/ends
//! discovery, frame-shiftability classification, backbone collection
//! and the two graph edits the option modes apply — the v1.05
//! cyclization (`Polymers105`) and the v1.06 frame shift (`Polymers`,
//! the senior-backbone-bond rotation). Everything here speaks the
//! SUBMITTED structure's numbering (the units' atom indices are
//! post-parse indices); the /z layer translates them through the
//! final canonical numbering at emit time.

use crate::error::{err, Result};
use crate::graph::Structure;

/// The SRU closure kind (the reference's CLOSING_SRU_* constants,
/// inpdef.h:54-59).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Closing {
    /// `CLOSING_SRU_NOT_APPLICABLE` — the unit never cyclizes.
    #[default]
    NotApplicable,
    /// `CLOSING_SRU_RING` — ends joined by a new single bond.
    Ring,
    /// `CLOSING_SRU_HIGHER_ORDER_BOND` — ends already adjacent; the
    /// bond order rises.
    HigherOrder,
    /// `CLOSING_SRU_DIRADICAL` — a single atom carries both caps.
    Diradical,
}

/// One treated unit: the reference's per-unit derived state (the
/// `OAD_PolymerUnit` fields the reader leaves unset).
#[derive(Debug, Clone, Default)]
pub struct TreatedUnit {
    /// In-unit ends of the two crossing bonds; `usize::MAX` when the
    /// unit has no crossing bonds.
    pub end1: usize,
    pub end2: usize,
    /// Out-of-unit caps of the two crossing bonds.
    pub cap1: usize,
    pub cap2: usize,
    /// Whether each cap is a star (the Zz pseudo-element).
    pub cap1_is_star: bool,
    pub cap2_is_star: bool,
    /// The closure kind the validation classified.
    pub cyclizable: Closing,
    /// The frame-shiftable backbone bonds (min-ordered 0-based index
    /// pairs), most-senior-first once [`Treated::sort_senior`] ran.
    pub bkbonds: Vec<(usize, usize)>,
    /// Whether the legacy cyclization edit ran on this unit.
    pub cyclized: bool,
    /// The aux `/Z:` echo's crossing-bond endpoint pairs (1-based, in
    /// the drawn order the reader resolved; a frame shift rewrites the
    /// shifted pair to `(end, cap)`).
    pub echo_pairs: Vec<(usize, usize)>,
}

impl TreatedUnit {
    /// Whether both caps are stars — the frame-shift eligibility's
    /// base condition (the reference's cap1_is_undef && cap2_is_undef).
    pub fn bistar(&self) -> bool {
        self.cap1_is_star && self.cap2_is_star
    }
}

/// The treated polymer of a whole structure.
#[derive(Debug, Clone, Default)]
pub struct Treated {
    pub units: Vec<TreatedUnit>,
    /// The polymer-related stars (0-based), in atom order — the
    /// reference's `pd->pzz`.
    pub pzz: Vec<usize>,
    /// The polymer-UNrelated pseudo atoms (0-based): non-empty
    /// rejects the structure unless `NPZz` stands.
    pub zy: Vec<usize>,
}

impl Treated {
    /// Whether any unit data exists.
    pub fn has_units(&self) -> bool {
        !self.units.is_empty()
    }
}

const CONN_HT: u8 = 1;

fn is_zz(element: &str) -> bool {
    element == "Zz"
}

/// The star atoms of the structure (0-based, in atom order).
pub fn star_atoms(structure: &Structure) -> Vec<usize> {
    structure
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| is_zz(&atom.element))
        .map(|(index, _)| index)
        .collect()
}

/// Neighbors of one atom: (neighbor, bond order) pairs.
fn neighbors(structure: &Structure, atom: usize) -> Vec<(usize, u32)> {
    structure
        .bonds
        .iter()
        .filter_map(|&(from, to, order)| {
            if from == atom {
                Some((to, order))
            } else if to == atom {
                Some((from, order))
            } else {
                None
            }
        })
        .collect()
}

/// The bond index between two atoms, if any.
fn bond_index(structure: &Structure, a: usize, b: usize) -> Option<usize> {
    structure
        .bonds
        .iter()
        .position(|&(from, to, _)| (from == a && to == b) || (from == b && to == a))
}

/// The polymer validation and derived state (the reference's
/// `OAD_ValidatePolymerAndPseudoElementData` +
/// `OAD_ValidateAndSortOutPseudoElementAtoms`). `polymers` selects
/// the treatment level; `npzz` admits non-polymeric stars. Without
/// either flag a star rejects the structure and unit data is ignored
/// (the observed default-options behavior: unit data vanishes, the
/// identifier is plain).
pub fn validate(structure: &Structure, polymers: bool, npzz: bool) -> Result<Treated> {
    let stars = star_atoms(structure);
    let Some(polymer) = &structure.polymer else {
        if !stars.is_empty() && !npzz {
            return err("Invalid element(s): Zz");
        }
        return Ok(Treated {
            units: Vec::new(),
            pzz: Vec::new(),
            zy: stars,
        });
    };
    if !polymers {
        // POLYMERS_NO: the unit data is ignored ("Ignore polymer
        // data"); the stars still gate.
        if !stars.is_empty() && !npzz {
            return err("Invalid element(s): Zz");
        }
        return Ok(Treated {
            units: Vec::new(),
            pzz: Vec::new(),
            zy: stars,
        });
    }

    let mut units = Vec::with_capacity(polymer.units.len());
    for unit in &polymer.units {
        let mut treated = TreatedUnit {
            end1: usize::MAX,
            end2: usize::MAX,
            cap1: usize::MAX,
            cap2: usize::MAX,
            cap1_is_star: false,
            cap2_is_star: false,
            cyclizable: Closing::NotApplicable,
            bkbonds: Vec::new(),
            cyclized: false,
            echo_pairs: Vec::new(),
        };
        if unit.alist.is_empty() {
            return err("Empty polymer unit");
        }
        if unit.alist.iter().any(|&a| a >= structure.atoms.len()) {
            return err("Invalid atom number in polymer unit");
        }
        // Crossing bonds: 0 or 2. The reader's blist entries are
        // 0-based BOND indices; the endpoints resolve here.
        let nb = unit.blist.len();
        if nb != 0 && nb != 2 {
            return err("Number of crossing bonds in polymer unit is not 0 or 2");
        }
        if nb == 2 {
            // Ends and caps: the in-unit / out-of-unit split of each
            // crossing bond (FindEndsAndCaps, runichi3.c 2166-2260).
            let mut ends_caps: Vec<Option<(usize, usize)>> = Vec::new();
            for &bond in &unit.blist {
                let Some(&(a, b, _)) = structure.bonds.get(bond) else {
                    return err("Invalid polymer CRU crossing bond");
                };
                treated.echo_pairs.push((a + 1, b + 1));
                let a_in = unit.alist.contains(&a);
                let b_in = unit.alist.contains(&b);
                if a_in && b_in {
                    return err("Caps of polymer unit lie inside it");
                }
                ends_caps.push(if a_in { Some((a, b)) } else { Some((b, a)) });
            }
            if let (Some((end1, cap1)), Some((end2, cap2))) = (ends_caps[0], ends_caps[1]) {
                treated.end1 = end1;
                treated.cap1 = cap1;
                treated.cap1_is_star = is_zz(&structure.atoms[cap1].element);
                treated.end2 = end2;
                treated.cap2 = cap2;
                treated.cap2_is_star = is_zz(&structure.atoms[cap2].element);
                if treated.cap1 == treated.cap2 {
                    return err("Invalid polymer CRU surrounding");
                }
                for cap in [treated.cap1, treated.cap2] {
                    if structure.atoms[cap].element == "H" {
                        return err("H as polymer end group is not supported");
                    }
                }
                // Closure classification (runichi3.c 1882-1920):
                // bistar units only, head-to-tail connections only.
                if treated.bistar() {
                    treated.cyclizable = if unit.alist.len() == 1 {
                        Closing::Diradical
                    } else if unit.alist.len() == 2
                        || neighbors(structure, treated.end1)
                            .iter()
                            .any(|&(nbr, _)| nbr == treated.end2)
                    {
                        Closing::HigherOrder
                    } else {
                        Closing::Ring
                    };
                    if unit.conn != CONN_HT {
                        treated.cyclizable = Closing::NotApplicable;
                    }
                }
            }
        }
        units.push(treated);
    }

    // Pseudo-atom accounting (runichi3.c 4396-4507): strictly
    // univalent stars; caps of bistar units are polymer-related, the
    // rest polymer-unrelated.
    for &star in &stars {
        if neighbors(structure, star).len() > 1 {
            return err("Invalid pseudo element(s) bonding");
        }
    }
    let mut pzz = Vec::new();
    let mut zy = Vec::new();
    for &star in &stars {
        let related = units.iter().any(|unit| {
            unit.bistar()
                && ((unit.cap1_is_star && unit.cap1 == star)
                    || (unit.cap2_is_star && unit.cap2 == star))
        });
        if related {
            pzz.push(star);
        } else {
            zy.push(star);
        }
    }
    if !zy.is_empty() && !npzz {
        return err("Polymer-unrelated pseudoatoms are not allowed");
    }
    if !pzz.is_empty() && !polymers {
        return err("Pseudoelement endgroups are not allowed");
    }

    Ok(Treated { units, pzz, zy })
}

/// The unit header digits the emitters speak after the
/// representation rewrites (runichi3.c 1733-1946): a unit is
/// structure-based when it carries crossing bonds; an all-source
/// polymer reads its SRU as MON, suppresses conn, defaults COP
/// subtype to RAN(2); structure-based SRU/MOD/CRO/MER default a
/// missing conn to EU(3).
pub fn header_digits(
    unit: &crate::graph::PolymerUnit,
    crossing: bool,
    any_structure: bool,
) -> (u8, u8, u8) {
    let (mut ty, mut subtype, mut conn) = (unit.unit_type, unit.subtype, unit.conn);
    if !crossing && !any_structure {
        if ty == 1 {
            ty = 2; // SRU -> MON
        }
        if ty == 3 && subtype == 0 {
            subtype = 2; // COP -> RAN
        }
        conn = 0;
    } else if crossing {
        if matches!(ty, 1 | 4 | 5 | 6) && conn == 0 {
            conn = 3; // EU
        }
        if ty == 3 && subtype == 0 {
            subtype = 2;
        }
    }
    (ty, subtype, conn)
}

/// The frame-shiftable backbone bonds of one unit: the union of bonds
/// over every simple path between the unit's ends within the unit's
/// subgraph (the reference's subgraf pathfinder, strutil.c 5270).
/// Min-ordered pairs in first-discovery order.
pub fn collect_backbone_bonds(
    structure: &Structure,
    alist: &[usize],
    end1: usize,
    end2: usize,
) -> Vec<(usize, usize)> {
    let mut bonds: Vec<(usize, usize)> = Vec::new();
    if end1 == usize::MAX || end2 == usize::MAX {
        return bonds;
    }
    fn walk(
        structure: &Structure,
        alist: &[usize],
        end: usize,
        seen: &mut Vec<usize>,
        bonds: &mut Vec<(usize, usize)>,
    ) {
        let current = *seen.last().expect("non-empty");
        for &(nbr, _) in &neighbors(structure, current) {
            if !alist.contains(&nbr) || seen.contains(&nbr) {
                continue;
            }
            if nbr == end {
                for window in seen.windows(2) {
                    let bond = (window[0].min(window[1]), window[0].max(window[1]));
                    if !bonds.contains(&bond) {
                        bonds.push(bond);
                    }
                }
                let bond = (current.min(end), current.max(end));
                if !bonds.contains(&bond) {
                    bonds.push(bond);
                }
                continue;
            }
            seen.push(nbr);
            walk(structure, alist, end, seen, bonds);
            seen.pop();
        }
    }
    walk(structure, alist, end2, &mut vec![end1], &mut bonds);
    bonds
}

/// Ring-system membership over the whole structure: atom index ->
/// ring-system id for atoms inside a cycle of length > 2 (the
/// reference's `OAD_Polymer_FindRingSystems`, with cyclized closure
/// bonds temporarily broken so artificial rings never count).
pub fn ring_systems(
    structure: &Structure,
    treated: &Treated,
) -> std::collections::HashMap<usize, usize> {
    let mut work = structure.clone();
    for unit in &treated.units {
        if unit.cyclized && unit.end1 != usize::MAX && unit.end2 != usize::MAX {
            if let Some(index) = bond_index(&work, unit.end1, unit.end2) {
                work.bonds.remove(index);
            }
        }
    }
    let n = work.atoms.len();
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(from, to, _) in &work.bonds {
        adjacency[from].push(to);
        adjacency[to].push(from);
    }
    let mut member = std::collections::HashMap::new();
    let mut visited = vec![false; n];
    let mut system = 0;
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut queue = vec![start];
        visited[start] = true;
        let mut nodes = Vec::new();
        let mut edges: Vec<(usize, usize)> = Vec::new();
        while let Some(atom) = queue.pop() {
            nodes.push(atom);
            for &nbr in &adjacency[atom] {
                let edge = (atom.min(nbr), atom.max(nbr));
                if !edges.contains(&edge) {
                    edges.push(edge);
                }
                if !visited[nbr] {
                    visited[nbr] = true;
                    queue.push(nbr);
                }
            }
        }
        // Cycle test: union-find over the component's edges.
        let mut parent: std::collections::HashMap<usize, usize> =
            nodes.iter().map(|&a| (a, a)).collect();
        fn find(parent: &mut std::collections::HashMap<usize, usize>, a: usize) -> usize {
            let p = parent[&a];
            if p == a {
                a
            } else {
                let root = find(parent, p);
                parent.insert(a, root);
                root
            }
        }
        let mut cyclic = false;
        for (a, b) in &edges {
            let (ra, rb) = (find(&mut parent, *a), find(&mut parent, *b));
            if ra == rb {
                cyclic = true;
            } else {
                parent.insert(ra, rb);
            }
        }
        if cyclic && nodes.len() > 2 {
            for node in nodes {
                member.insert(node, system);
            }
            system += 1;
        }
    }
    member
}

/// The legacy cyclization edit (`OAD_Polymer_CyclizeCloseableUnits`):
/// for every cyclizable unit, disconnect the caps and join the ends —
/// a new single bond (Ring), a raised order (HigherOrder) or a
/// triplet-radical reading (Diradical, carried as a pinned zero-H
/// atom — the chemomer has no radical layer). Marks units cyclized.
pub fn cyclize_units(structure: &mut Structure, treated: &mut Treated) {
    for index in 0..treated.units.len() {
        let (end1, end2, cap1, cap2, kind) = {
            let unit = &treated.units[index];
            if unit.cyclizable == Closing::NotApplicable {
                continue;
            }
            (unit.end1, unit.end2, unit.cap1, unit.cap2, unit.cyclizable)
        };
        match kind {
            Closing::Ring => {
                remove_bond(structure, cap1, end1);
                remove_bond(structure, cap2, end2);
                if bond_index(structure, end1, end2).is_none() {
                    structure.bonds.push((end1, end2, 1));
                }
            }
            Closing::HigherOrder => {
                raise_order(structure, end1, end2);
                remove_bond(structure, cap1, end1);
                remove_bond(structure, cap2, end2);
            }
            Closing::Diradical => {
                structure.atoms[end1].h_count = 0;
                structure.atoms[end1].h_pinned = true;
                remove_bond(structure, cap1, end1);
                remove_bond(structure, cap2, end1);
            }
            Closing::NotApplicable => {}
        }
        treated.units[index].cyclized = true;
    }
}

fn remove_bond(structure: &mut Structure, a: usize, b: usize) {
    if let Some(index) = bond_index(structure, a, b) {
        structure.bonds.remove(index);
    }
}

fn raise_order(structure: &mut Structure, a: usize, b: usize) {
    if let Some(index) = bond_index(structure, a, b) {
        structure.bonds[index].2 = (structure.bonds[index].2 + 1).min(3);
    }
}

/// The seniority property tables (the reference's `erank_rule2` /
/// `erank_rule4`, runichi3.c 3312-3348), indexed by REAL periodic
/// number; 0 for anything beyond the table (the pseudo elements,
/// periodic 119+).
const ERANK_RULE2: [u16; 120] = [
    0, 1, 198, 197, 196, 202, 2, 216, 215, 191, 190, 189, 188, 187, 206, 210, 214, 183, 182, 181,
    180, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169, 205, 209, 213, 165, 164, 163, 162,
    161, 160, 159, 158, 157, 156, 155, 154, 153, 152, 151, 204, 208, 212, 147, 146, 145, 144, 143,
    142, 141, 140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 130, 129, 128, 127, 126, 125, 124,
    123, 122, 121, 201, 119, 203, 207, 116, 115, 114, 113, 112, 111, 110, 109, 108, 107, 106, 105,
    104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83,
    82, 81,
];
const ERANK_RULE4: [u16; 120] = [
    0, 1, 198, 197, 196, 202, 2, 211, 215, 191, 190, 189, 188, 187, 206, 210, 214, 183, 182, 181,
    180, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169, 205, 209, 213, 165, 164, 163, 162,
    161, 160, 159, 158, 157, 156, 155, 154, 153, 152, 151, 204, 208, 212, 147, 146, 145, 144, 143,
    142, 141, 140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 130, 129, 128, 127, 126, 125, 124,
    123, 122, 121, 201, 119, 203, 207, 116, 115, 114, 113, 112, 111, 110, 109, 108, 107, 106, 105,
    104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83,
    82, 81,
];

fn erank(table: &[u16; 120], element: &str) -> u16 {
    match crate::canonical::protons::el_number(element) {
        Some(n) if (n as usize) < table.len() => table[n as usize],
        _ => 0,
    }
}

/// The seniority props of one atom (the reference's `OAD_AtProps`).
#[derive(Debug, Clone, Copy, Default)]
pub struct AtomProps {
    pub erank: u16,
    pub ring_erank: u16,
    pub ring_size: usize,
}

/// The seniority props over the whole structure (the reference's
/// `OAD_Polymer_SetAtProps`): rule-4 element ranks, ring systems
/// scored by their max rule-2 erank. Keyed by 0-based atom index.
pub fn atom_props(structure: &Structure, treated: &Treated) -> Vec<AtomProps> {
    let ring = ring_systems(structure, treated);
    let mut props = vec![AtomProps::default(); structure.atoms.len()];
    for (index, atom) in structure.atoms.iter().enumerate() {
        props[index].erank = erank(&ERANK_RULE4, &atom.element);
    }
    let mut max_by_ring: std::collections::HashMap<usize, u16> = Default::default();
    for (&atom, &id) in &ring {
        let rule2 = erank(&ERANK_RULE2, &structure.atoms[atom].element);
        max_by_ring
            .entry(id)
            .and_modify(|m| *m = (*m).max(rule2))
            .or_insert(rule2);
        props[atom].ring_size = 3; // membership (the size > 2 gate)
    }
    for (&atom, &id) in &ring {
        props[atom].ring_erank = max_by_ring[&id];
    }
    props
}

/// The atom-seniority compare (the reference's
/// `OAD_Polymer_CompareRanksOfTwoAtoms`, IUPAC rule 1): -1 when the
/// first atom is senior, +1 when the second is, 0 unresolved.
pub fn compare_ranks(a: usize, b: usize, props: &[AtomProps]) -> i32 {
    const HETEROCYC: i32 = 3;
    const HETEROAT: i32 = 2;
    const CARBOCYC: i32 = 1;
    const CARBOAT: i32 = 0;
    let typ = |atom: usize| -> i32 {
        let p = &props[atom];
        if p.ring_size > 2 {
            if p.ring_erank <= 2 {
                CARBOCYC
            } else {
                HETEROCYC
            }
        } else if p.erank == 2 {
            CARBOAT
        } else {
            HETEROAT
        }
    };
    let (ta, tb) = (typ(a), typ(b));
    if ta == HETEROCYC && tb == HETEROCYC {
        if props[a].ring_erank < props[b].ring_erank {
            return 1;
        }
        if props[a].ring_erank > props[b].ring_erank {
            return -1;
        }
        if props[a].ring_size < props[b].ring_size {
            return 1;
        }
        if props[a].ring_size > props[b].ring_size {
            return -1;
        }
        return 0;
    }
    if ta == HETEROCYC {
        return -1;
    }
    if tb == HETEROCYC {
        return 1;
    }
    if ta == HETEROAT && tb == HETEROAT {
        if props[a].erank < props[b].erank {
            return 1;
        }
        if props[a].erank > props[b].erank {
            return -1;
        }
        return 0;
    }
    if ta == HETEROAT {
        return -1;
    }
    if tb == HETEROAT {
        return 1;
    }
    if ta == CARBOCYC && tb == CARBOCYC {
        if props[a].ring_size < props[b].ring_size {
            return 1;
        }
        if props[a].ring_size > props[b].ring_size {
            return -1;
        }
        return 0;
    }
    if ta == CARBOCYC {
        return -1;
    }
    if tb == CARBOCYC {
        return 1;
    }
    0
}

/// The always-resolving atom compare (the reference's
/// `OAD_Polymer_IsFirstAtomRankLower`): +1 first junior, -1 first
/// senior, by rank then by index.
pub fn is_first_rank_lower(a: usize, b: usize, props: &[AtomProps]) -> i32 {
    let result = compare_ranks(a, b, props);
    if result != 0 {
        return result;
    }
    if a < b {
        1
    } else if a > b {
        -1
    } else {
        0
    }
}

/// The backbone-bond seniority compare (the reference's
/// `OAD_Polymer_CompareBackboneBondsSeniority`): favor the greater
/// max-rank end, then the lesser min-rank end, then indices. -1 = the
/// first bond is senior.
pub fn compare_backbone_bonds(b1: (usize, usize), b2: (usize, usize), props: &[AtomProps]) -> i32 {
    let ext = |bond: (usize, usize)| -> (usize, usize) {
        let (first, second) = bond;
        if is_first_rank_lower(second, first, props) == -1 {
            (second, first)
        } else {
            (first, second)
        }
    };
    let (b1max, b1min) = ext(b1);
    let (b2max, b2min) = ext(b2);
    let cmp = compare_ranks(b1max, b2max, props);
    if cmp == 1 {
        return 1;
    }
    if cmp == -1 {
        return -1;
    }
    let cmp = compare_ranks(b1min, b2min, props);
    if cmp == 1 {
        return -1;
    }
    if cmp == -1 {
        return 1;
    }
    if b1max < b2max {
        return 1;
    }
    if b1max > b2max {
        return -1;
    }
    if b1min < b2min {
        return -1;
    }
    if b1min > b2min {
        return 1;
    }
    0
}

/// Sorts the backbone bonds most-senior-first (the reference's
/// `OAD_PolymerUnit_SortBackboneBondsAndSetSeniors`): the senior bond
/// moves to the front with its senior atom first.
pub fn sort_backbone_bonds(bonds: &[(usize, usize)], props: &[AtomProps]) -> Vec<(usize, usize)> {
    let mut order: Vec<usize> = (0..bonds.len()).collect();
    order.sort_by(|&a, &b| compare_backbone_bonds(bonds[a], bonds[b], props).cmp(&0));
    let mut sorted: Vec<(usize, usize)> = order.iter().map(|&i| bonds[i]).collect();
    if let Some(pair) = sorted.first_mut() {
        if is_first_rank_lower(pair.0, pair.1, props) == 1 {
            *pair = (pair.1, pair.0);
        }
    }
    sorted
}

/// The modern frame-shift edit derivation for one unit
/// (`OAD_Polymer_PrepareFrameShiftEdits`): given the senior backbone
/// bond's ends (in the submitted numbering, senior atom first) and
/// the unit's caps, produce the three bond rewrites — each cap to its
/// senior end (skipped when unchanged), then the senior bond cut and
/// the original ends joined. The vec carries (remove_a, remove_b,
/// add_a, add_b) tuples.
pub fn frame_shift_edits(
    unit: &TreatedUnit,
    end1: usize,
    end2: usize,
) -> Vec<(usize, usize, usize, usize)> {
    let mut edits = Vec::new();
    let old_end1 = unit.end1;
    let old_end2 = unit.end2;
    let same = |a1: usize, a2: usize, b1: usize, b2: usize| {
        (a1 == b1 && a2 == b2) || (a1 == b2 && a2 == b1)
    };
    if !same(old_end1, unit.cap1, end1, unit.cap1) {
        edits.push((old_end1, unit.cap1, end1, unit.cap1));
    }
    if !same(old_end2, unit.cap2, end2, unit.cap2) {
        edits.push((old_end2, unit.cap2, end2, unit.cap2));
    }
    if !same(end1, end2, old_end1, old_end2) {
        edits.push((end1, end2, old_end1, old_end2));
    }
    edits
}

/// Applies one bond rewrite to the structure.
pub fn apply_edit(structure: &mut Structure, edit: &(usize, usize, usize, usize)) {
    let (rm_a, rm_b, add_a, add_b) = *edit;
    remove_bond(structure, rm_a, rm_b);
    structure.bonds.push((add_a, add_b, 1));
}

/// The backbone bonds for every unit (`OAD_Polymer_FindBackbones`):
/// all simple end-to-end paths within each unit, then the delists —
/// intra-ring bonds (the ring perception breaks cyclized closure
/// bonds), higher-order bonds and the bonds the tautomer registration
/// marked (`BOND_TAUTOM` — the cyclization bond itself marks when its
/// closure creates a tautomeric system, the N-oxide corner) — and the
/// ends-bond held when everything delisted.
pub fn fill_backbones(structure: &Structure, poly: &crate::graph::Polymer, treated: &mut Treated) {
    let ring = ring_systems(structure, treated);
    for index in 0..treated.units.len() {
        let state = &mut treated.units[index];
        if state.cyclizable == Closing::NotApplicable {
            continue;
        }
        if state.end1 == usize::MAX || state.end2 == usize::MAX {
            continue;
        }
        let Some(unit) = poly.units.get(index) else {
            continue;
        };
        let mut bonds = collect_backbone_bonds(structure, &unit.alist, state.end1, state.end2);
        bonds.retain(|&(a, b)| match (ring.get(&a), ring.get(&b)) {
            (Some(&x), Some(&y)) => x != y,
            _ => true,
        });
        // Higher-order bonds never carry the frame shift (the
        // reference's DelistHighOrderBackboneBonds namesake filter).
        bonds.retain(|&(a, b)| {
            bond_index(structure, a, b)
                .map(|index| structure.bonds[index].2 <= 1)
                .unwrap_or(true)
        });
        // Delist the tautomeric backbone bonds (the reference's
        // OAD_PolymerUnit_DelistHighOrderBackboneBonds over the
        // normalized composite: a bond the t-group registration made
        // BOND_TAUTOM never carries the frame shift — the cyclization
        // bond itself marks when its closure creates a tautomeric
        // system, the N-oxide corner).
        let counts = crate::layers::standardization::adjusted_hydrogens(structure);
        let marks = crate::canonical::tautomer::tautomeric_bond_marks(structure, &counts);
        bonds.retain(|&(a, b)| match bond_index(structure, a, b) {
            Some(index) => !marks[index],
            None => true,
        });
        if bonds.is_empty() {
            // Hold the ends bond when every candidate delisted (the
            // reference's "single frame-shift bond" fallback).
            bonds.push((state.end1.min(state.end2), state.end1.max(state.end2)));
        }
        state.bkbonds = bonds;
    }
}

/// Whether the FoldCRU attempt fails the way the oracle's does. The
/// reference folds a CRU whose backbone fragments between cuts repeat
/// (runichi2.c `analyze_CRU_folding`): fragments split at every
/// single, non-cyclic backbone bond, signed by the interim's blind
/// equivalence classes, and a class sequence with a repeating period
/// held ≥ 2 whole times prepares fold edits. Every attempt that
/// prepares edits then fails downstream (the edited unit's crossing
/// bonds no longer read 0 or 2), which is the identifier-level error
/// the option produces; CRUs without a repeating pattern fall through
/// and the identifier speaks unedited. Multi-component structures
/// skip silently (the reference's aux-class extraction reads one
/// component — the two-cru recording stays plain).
pub fn cru_folding_fails(
    structure: &Structure,
    poly: &crate::graph::Polymer,
    treated: &Treated,
) -> bool {
    if !single_component(structure) {
        return false;
    }
    // The interim cyclized copy supplies the equivalence classes (the
    // reference reads them from the 105+ aux).
    let mut interim = structure.clone();
    let mut interim_treated = treated.clone();
    cyclize_units(&mut interim, &mut interim_treated);
    fill_backbones(&interim, poly, &mut interim_treated);
    let counts = crate::layers::standardization::adjusted_hydrogens(&interim);
    let (adjacency, invariants) = crate::canonical::graph_of(&interim);
    let mask = crate::canonical::tautomer::endpoint_mask(&interim, &counts);
    let groups = crate::canonical::group_membership(&interim, &counts);
    let (_, symmetry) = crate::canonical::canonical_ranks_and_symmetry(
        &adjacency,
        &invariants,
        &crate::canonical::tautomer::without_mobile_h(&counts, &mask),
        &groups,
    );

    for (index, state) in interim_treated.units.iter().enumerate() {
        let Some(unit) = poly.units.get(index) else {
            continue;
        };
        // Bi-star CRUs with at least two atoms and two crossing bonds.
        if !state.bistar() || unit.alist.len() < 2 || unit.blist.len() < 2 {
            continue;
        }
        // Cuts: the delisted backbone bonds that survive as single
        // bonds of the submitted structure.
        let cuts: Vec<(usize, usize)> = state
            .bkbonds
            .iter()
            .copied()
            .filter(|&(a, b)| {
                bond_index(structure, a, b)
                    .map(|bond| structure.bonds[bond].2 <= 1)
                    .unwrap_or(false)
            })
            .collect();
        if cuts.is_empty() {
            continue;
        }
        let Some(fragments) = backbone_fragments(structure, unit, state, &cuts) else {
            continue;
        };
        // Class-count signatures, classes numbered by first sighting.
        let mut signatures: Vec<Vec<usize>> = Vec::new();
        let mut sequence: Vec<usize> = Vec::with_capacity(fragments.len());
        for fragment in &fragments {
            let mut signature: Vec<usize> = fragment
                .iter()
                .map(|&atom| symmetry[atom] as usize)
                .collect();
            signature.sort_unstable();
            let class = match signatures.iter().position(|known| *known == signature) {
                Some(class) => class,
                None => {
                    signatures.push(signature);
                    signatures.len() - 1
                }
            };
            sequence.push(class);
        }
        // A repeating pattern held at least twice.
        let total = sequence.len();
        for period in 1..total {
            if !total.is_multiple_of(period) || total / period < 2 {
                continue;
            }
            if (0..total).all(|step| sequence[step] == sequence[step % period]) {
                return true;
            }
        }
    }
    false
}

/// The backbone fragments: the connected pieces of the unit's
/// subgraph after the cut bonds come out, ordered along the backbone
/// from end1 (a BFS over the piece adjacency the cuts join).
fn backbone_fragments(
    structure: &Structure,
    unit: &crate::graph::PolymerUnit,
    state: &TreatedUnit,
    cuts: &[(usize, usize)],
) -> Option<Vec<Vec<usize>>> {
    let in_unit = |atom: usize| unit.alist.contains(&atom);
    // Union-find over the unit's atoms; cut bonds never merge.
    let mut parent: std::collections::HashMap<usize, usize> = Default::default();
    let find = |parent: &std::collections::HashMap<usize, usize>, mut atom: usize| -> usize {
        while let Some(&up) = parent.get(&atom) {
            if up == atom {
                break;
            }
            atom = up;
        }
        atom
    };
    for &atom in &unit.alist {
        parent.insert(atom, atom);
    }
    for &(from, to, _) in &structure.bonds {
        if !in_unit(from) || !in_unit(to) || cuts.contains(&(from.min(to), from.max(to))) {
            continue;
        }
        let (ra, rb) = (find(&parent, from), find(&parent, to));
        if ra != rb {
            parent.insert(ra, rb);
        }
    }
    let mut members: std::collections::HashMap<usize, Vec<usize>> = Default::default();
    for &atom in &unit.alist {
        members.entry(find(&parent, atom)).or_default().push(atom);
    }
    let pieces: Vec<Vec<usize>> = members.into_values().collect();
    if pieces.len() != cuts.len() + 1 {
        return None;
    }
    // Order the pieces from the end1 side: walk the cut graph.
    let root_of = |atom: usize| find(&parent, atom);
    let mut order: Vec<Vec<usize>> = Vec::with_capacity(pieces.len());
    let mut seen: Vec<usize> = Vec::new();
    let mut queue = vec![root_of(state.end1)];
    while let Some(root) = queue.pop() {
        if seen.contains(&root) {
            continue;
        }
        seen.push(root);
        let piece = pieces.iter().find(|piece| root_of(piece[0]) == root)?;
        order.push(piece.clone());
        for &(a, b) in cuts {
            if root_of(a) == root {
                queue.push(root_of(b));
            } else if root_of(b) == root {
                queue.push(root_of(a));
            }
        }
    }
    Some(order)
}

/// Whether the structure's atoms form one connected component.
fn single_component(structure: &Structure) -> bool {
    let n = structure.atoms.len();
    if n <= 1 {
        return true;
    }
    let mut seen = vec![false; n];
    let mut queue = vec![0];
    seen[0] = true;
    let mut count = 1;
    while let Some(atom) = queue.pop() {
        for &(from, to, _) in &structure.bonds {
            let other = if from == atom {
                to
            } else if to == atom {
                from
            } else {
                continue;
            };
            if !seen[other] {
                seen[other] = true;
                count += 1;
                queue.push(other);
            }
        }
    }
    count == n
}
