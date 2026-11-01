//! Mobile-H (tautomer) detection — M3, ported from the reference's
//! ichitaut semantics. An endpoint is O/S/Se/Te (endpoint valence 2)
//! or N (3) whose attached hydrogens plus bond-order sum exactly fill
//! the endpoint valence; the fill decides the role — no multiple bond
//! = donor (the H's owner), exactly one double-bond worth = acceptor.
//! Overflowing the endpoint valence disqualifies an atom outright, so
//! a nitro nitrogen (bond sum 4 over N's 3) keeps its aci hydrogens
//! fixed and a sulfone/sulfate sulfur (bond sum ≥ 3 over S's 2)
//! never participates — only thiocarbonyl =S and S-H do.
//! Endpoints link through a shared anchor that holds both a donor and
//! an acceptor (the 1,3 template, the donor's bond to the anchor a
//! real single bond carrying its H), and along alternating
//! single/aromatic/double shift walks whose far terminal is the
//! COMPLEMENTARY role — a phenolic O can never reach another phenolic
//! O across an all-carbon aromatic arc because both terminals are
//! donors. An aromatic run inside a walk is entered from the donor's
//! own bond or a double bond and left into a double bond or the donor
//! terminal: shift paths never jump π systems through an interior
//! chain single. The canonicalization ranks atoms by hydrogens with mobile
//! H removed (the reference's numHNoTautH — see
//! [`without_mobile_h`]), which is what settles terminal-O ordering;
//! a group's hydrogens are the ±1-standardized counts of its members.

use crate::graph::Structure;
use std::collections::{BTreeSet, HashMap};

/// One mobile-H group: the hydrogen count and the canonical ranks of
/// its endpoints (ascending).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MobileGroup {
    /// Mobile units: hydrogens plus negative charges (the reference's
    /// t-group num[0]).
    pub count: u32,
    /// Negative charges among the members (num[1]); the writer emits
    /// them as the group's minus suffix.
    pub charge: u32,
    pub ranks: Vec<u32>,
}

/// Endpoint valence per the reference's `get_endpoint_valence`.
fn endpoint_valence(element: &str) -> Option<u32> {
    match element {
        "O" | "S" | "Se" | "Te" => Some(2),
        "N" => Some(3),
        _ => None,
    }
}

/// Endpoint valence per the reference's `get_endpoint_valence_KET`
/// (util.c:1530): only C and O carry keto-enol endpoints.
fn endpoint_valence_ket(element: &str) -> Option<u32> {
    match element {
        "C" => Some(4),
        "O" => Some(2),
        _ => None,
    }
}

/// Which endpoint predicate a round consults — the reference's
/// `nGetEndpointInfo` variant per experimental flag (ichitaut.c 359
/// base, 452 PT_22_00, 524 PT_16_00, 600 PT_06_00, 677 PT_39_00,
/// 756 PT_13_00, 832 PT_18_00, 916 KET).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Predicate {
    /// The standard heteroatom predicate (O/S/Se/Te 2, N 3).
    Base,
    /// Keto-enol: C 4, O 2, with the connectivity bands.
    Ket,
    /// PT_22_00: C only.
    Pt22,
    /// PT_16_00: the KET table with the same bands.
    Pt16,
    /// PT_06_00: C 4 plus the base table.
    Pt06,
    /// PT_39_00: C 4, N 3.
    Pt39,
    /// PT_13_00: C 4, O/S/Se/Te 2 — the role reads the mobile units,
    /// not the bond arithmetic.
    Pt13,
    /// PT_18_00: O 2, N 3 — no connectivity checks at all.
    Pt18,
}

/// One predicate verdict's full numbers: the role plus the mobile
/// units and the neutral-bonds valence the shift checks read (the
/// reference's `eif.cMobile` / `eif.cNeutralBondsValence`).
#[derive(Clone, Copy, Debug)]
struct EndpointInfo {
    role: Role,
    /// The endpoint's mobile units (H plus an anionic charge).
    mobile: i32,
    /// `endpoint_valence - mobile`: a shift may land here when this
    /// exceeds the heavy-atom valence (the acceptor-side test).
    neutral_bonds_valence: i32,
    /// The endpoint valence itself (the reference's return value);
    /// the KET/PT gates count endpoints by it.
    valence: i32,
}

/// How an endpoint's bond-order sum relates to its heavy-atom valence:
/// all single bonds = [`Role::Donor`], one double-bond worth =
/// [`Role::Acceptor`].
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Role {
    Donor,
    Acceptor,
}

/// The ring perception of the reference's MarkRingSystemsInp: a DFS
/// segmentation at `dfs == low` gives every atom a ring-system id and
/// size — a chain atom is a size-1 system of its own, rings joined by
/// a bridge are separate systems, spiro fusions are one — plus the
/// articulation points the 1,4 template requires.
struct RingSystems {
    system: Vec<usize>,
    size: Vec<usize>,
    cut_vertex: Vec<bool>,
}

fn ring_systems(structure: &Structure) -> RingSystems {
    let n = structure.atoms.len();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(from, to, _) in &structure.bonds {
        neighbors[from].push(to);
        neighbors[to].push(from);
    }
    let mut system = vec![0usize; n];
    let mut size = vec![0usize; n];
    let mut cut_vertex = vec![false; n];

    // Articulation points: a non-root atom whose subtree cannot reach
    // above it; the root only when it has two independent children.
    {
        let mut disc = vec![0u32; n];
        #[allow(clippy::too_many_arguments)]
        fn visit(
            atom: usize,
            parent: Option<usize>,
            parent_edge_used: &mut bool,
            neighbors: &[Vec<usize>],
            disc: &mut Vec<u32>,
            low: &mut Vec<u32>,
            cut_vertex: &mut Vec<bool>,
            timer: &mut u32,
            root: usize,
            root_children: &mut u32,
        ) -> u32 {
            *timer += 1;
            disc[atom] = *timer;
            let mut low_here = *timer;
            for &next in &neighbors[atom] {
                if Some(next) == parent && !*parent_edge_used {
                    *parent_edge_used = true;
                    continue;
                }
                if disc[next] == 0 {
                    *root_children += (atom == root) as u32;
                    let child_low = visit(
                        next,
                        Some(atom),
                        &mut false,
                        neighbors,
                        disc,
                        low,
                        cut_vertex,
                        timer,
                        root,
                        root_children,
                    );
                    low_here = low_here.min(child_low);
                    if child_low >= disc[atom] && (atom != root || *root_children > 1) {
                        cut_vertex[atom] = true;
                    }
                } else {
                    low_here = low_here.min(disc[next]);
                }
            }
            low[atom] = low_here;
            low_here
        }
        for root in 0..n {
            if disc[root] == 0 {
                let mut root_children = 0u32;
                visit(
                    root,
                    None,
                    &mut false,
                    &neighbors,
                    &mut disc,
                    &mut vec![0; n],
                    &mut cut_vertex,
                    &mut 0,
                    root,
                    &mut root_children,
                );
            }
        }
    }

    // Ring-system segmentation: on retracting from an atom whose low
    // equals its own dfs, the DFS stack suffix from it unwinds into a
    // numbered system of that size.
    {
        let mut disc = vec![0u32; n];
        let mut timer = 0u32;
        let mut stack: Vec<usize> = Vec::new();
        let mut next_system = 0usize;
        #[allow(clippy::too_many_arguments)]
        fn visit(
            atom: usize,
            parent: Option<usize>,
            parent_edge_used: &mut bool,
            neighbors: &[Vec<usize>],
            disc: &mut Vec<u32>,
            timer: &mut u32,
            stack: &mut Vec<usize>,
            system: &mut Vec<usize>,
            size: &mut Vec<usize>,
            next_system: &mut usize,
        ) -> u32 {
            *timer += 1;
            disc[atom] = *timer;
            let mut low_here = *timer;
            stack.push(atom);
            for &next in &neighbors[atom] {
                if Some(next) == parent && !*parent_edge_used {
                    *parent_edge_used = true;
                    continue;
                }
                let child_low = if disc[next] == 0 {
                    visit(
                        next,
                        Some(atom),
                        &mut false,
                        neighbors,
                        disc,
                        timer,
                        stack,
                        system,
                        size,
                        next_system,
                    )
                } else {
                    disc[next]
                };
                low_here = low_here.min(child_low);
            }
            if disc[atom] == low_here {
                *next_system += 1;
                let mut count = 0usize;
                while let Some(top) = stack.pop() {
                    system[top] = *next_system;
                    count += 1;
                    if top == atom {
                        break;
                    }
                }
                size[atom] = count;
                let top = system.len();
                for index in 0..top {
                    if system[index] == *next_system {
                        size[index] = count;
                    }
                }
            }
            low_here
        }
        for root in 0..n {
            if disc[root] == 0 {
                visit(
                    root,
                    None,
                    &mut false,
                    &neighbors,
                    &mut disc,
                    &mut timer,
                    &mut stack,
                    &mut system,
                    &mut size,
                    &mut next_system,
                );
            }
        }
    }

    RingSystems {
        system,
        size,
        cut_vertex,
    }
}

/// The reference's SetForbiddenEdges under BNS_PROTECT_FROM_TAUT: a
/// carbonyl carbon carrying two valence-1 oxygens whose bond orders
/// sum to 3 (the acetyl pattern) or a nitro-style nitrogen carrying
/// two valence-1 oxygens summing to 3 or 4 forbids its single-bonded
/// non-oxygen bond — the bond the templates may never cross
/// (ALLOWED_EDGE) and the network may never route.
fn forbidden_edges(
    structure: &Structure,
    neighbors: &[Vec<(usize, u32)>],
    bond_index: &HashMap<(usize, usize), usize>,
) -> BTreeSet<usize> {
    // Bond values as AllocateAndInitBnStruct leaves them: aromatic
    // input reads single.
    let value = |kind: u32| match kind {
        2 => 2,
        3 => 3,
        _ => 1,
    };
    let mut forbidden = BTreeSet::new();
    for center in 0..structure.atoms.len() {
        let element = structure.atoms[center].element.as_str();
        let degree = neighbors[center].len();
        let cbv = neighbors[center]
            .iter()
            .map(|&(_, kind)| value(kind))
            .sum::<i32>();
        let wanted_o_sum: Option<i32> = match (element, degree, cbv) {
            ("C", 3, 4) => Some(3),
            ("N", 3, 4) => Some(3),
            ("N", 3, 5) => Some(4),
            _ => None,
        };
        let Some(wanted_o_sum) = wanted_o_sum else {
            continue;
        };
        let (mut oxy_sum, mut oxy_count) = (0, 0);
        let mut other: Option<(usize, u32)> = None;
        for &(neigh, kind) in &neighbors[center] {
            if structure.atoms[neigh].element == "O" && neighbors[neigh].len() == 1 {
                oxy_sum += value(kind);
                oxy_count += 1;
            } else if other.is_none() {
                other = Some((neigh, kind));
            } else {
                other = None;
                break;
            }
        }
        let Some((neigh, kind)) = other else {
            continue;
        };
        if oxy_count == 2 && oxy_sum == wanted_o_sum && value(kind) == 1 {
            if let Some(&index) = bond_index.get(&(center.min(neigh), center.max(neigh))) {
                forbidden.insert(index);
            }
        }
    }
    forbidden
}

/// Cached per-atom structure facts the endpoint arithmetic needs.
struct Context<'a> {
    structure: &'a Structure,
    hydrogens: &'a [u32],
    neighbors: Vec<Vec<(usize, u32)>>,
    bond_index: HashMap<(usize, usize), usize>,
    /// Perceived-aromatic bonds: the reference's normalization marks
    /// them BOND_ALTERN while the valence arithmetic keeps the input
    /// kekule's integer orders.
    aromatic: BTreeSet<usize>,
    /// Bond-order sums of the INPUT kekule orders, doubled.
    bond_sum2: Vec<i32>,
    /// The MarkRingSystemsInp perception the templates scan.
    rings: RingSystems,
    /// The reference's SetForbiddenEdges (BNS_PROTECT_FROM_TAUT):
    /// acetyl-carbon and nitro-nitrogen bonds the tautomer machinery
    /// may never cross (ALLOWED_EDGE) nor the network route.
    forbidden: BTreeSet<usize>,
    /// Deprotonated endpoints' mobile units (the reference's
    /// `nMobile += (charge == -1)`): on a fully pinned structure —
    /// the (DE)PROTONATION-standardized state — a charge −1 counts
    /// as a mobile unit exactly where no riding hydrogen stands.
    charge_units: Vec<u32>,
    /// Per-bond mobility marks the hard proton ride left (the
    /// reference's ALT12NS / ALTERN|TRIPLE bond types); zero where
    /// the bond carries none.
    ride_marks: Vec<u8>,
}

impl<'a> Context<'a> {
    fn new(structure: &'a Structure, hydrogens: &'a [u32]) -> Self {
        let n = structure.atoms.len();
        let aromatic = crate::canonical::aromatic::aromatic_bonds(structure);
        let mut neighbors = vec![Vec::new(); n];
        let mut bond_index = HashMap::new();
        let mut bond_sum2 = vec![0i32; n];
        let pinned = structure.atoms.iter().all(|atom| atom.h_pinned);
        let charge_units: Vec<u32> = structure
            .atoms
            .iter()
            .map(|atom| (pinned && atom.charge == -1) as u32)
            .collect();
        for (index, &(from, to, kind)) in structure.bonds.iter().enumerate() {
            // The INPUT kekule orders carry the valence arithmetic.
            let order2 = match kind {
                2 => 4,
                3 => 6,
                4 => 3,
                _ => 2,
            };
            neighbors[from].push((to, kind));
            neighbors[to].push((from, kind));
            bond_index.insert((from.min(to), from.max(to)), index);
            bond_sum2[from] += order2;
            bond_sum2[to] += order2;
        }
        let rings = ring_systems(structure);
        let forbidden = forbidden_edges(structure, &neighbors, &bond_index);
        let mut ride_marks = vec![0u8; structure.bonds.len()];
        for (index, &mark) in structure.ride_marks.iter().enumerate() {
            if index < ride_marks.len() {
                ride_marks[index] = mark;
            }
        }
        Context {
            structure,
            hydrogens,
            neighbors,
            bond_index,
            aromatic,
            bond_sum2,
            rings,
            forbidden,
            charge_units,
            ride_marks,
        }
    }

    /// The reference's `ALLOWED_EDGE`: false when the bond carries a
    /// SetForbiddenEdges mark.
    fn allowed_edge(&self, a: usize, b: usize) -> bool {
        self.bond_index
            .get(&(a.min(b), a.max(b)))
            .is_none_or(|index| !self.forbidden.contains(index))
    }

    /// The reference's `nGetEndpointInfo`.
    fn endpoint(&self, state: &Registration, atom: usize) -> Option<Role> {
        let record = &self.structure.atoms[atom];
        let valence = endpoint_valence(&record.element)? as i32;
        let heavy = self.neighbors[atom].len() as i32;
        if valence <= heavy {
            return None;
        }
        if record.charge > 0 {
            // The charge-point branch: a +1 atom inside a c-group
            // whose subtype can hold or donate H is an endpoint.
            if state.c_point_of[atom].is_some() {
                if let Some((_, subtype)) = charge_type(self, state, atom) {
                    if subtype & 4 != 0 {
                        return Some(Role::Acceptor);
                    }
                    if subtype & 8 != 0 {
                        return Some(Role::Donor);
                    }
                }
            }
            return None;
        }
        let mobile = self.mobile(atom) as i32;
        if 2 * mobile + self.bond_sum2[atom] != 2 * valence {
            return None;
        }
        match self.bond_sum2[atom] - 2 * heavy {
            0 => Some(Role::Donor),
            2 => Some(Role::Acceptor),
            _ => None,
        }
    }

    /// The reference's `nGetEndpointInfo` under a flag predicate: the
    /// valence table, the connectivity bands, and the role arithmetic
    /// each flag's function carries (ichitaut.c 452–1080). The base
    /// predicate routes through [`Context::endpoint`] unchanged.
    fn endpoint_by(&self, state: &Registration, atom: usize, pred: Predicate) -> Option<Role> {
        match pred {
            Predicate::Base => self.endpoint(state, atom),
            _ => self
                .endpoint_info_ext(state, atom, pred)
                .map(|info| info.role),
        }
    }

    /// [`Context::endpoint`] plus the shift-check numbers: the mobile
    /// units and the neutral-bonds valence of the base predicate's
    /// verdict.
    fn endpoint_info(&self, state: &Registration, atom: usize) -> Option<EndpointInfo> {
        let role = self.endpoint(state, atom)?;
        let record = &self.structure.atoms[atom];
        let valence = endpoint_valence(&record.element)? as i32;
        let mobile = if record.charge > 0 {
            self.hydrogens[atom] as i32
        } else {
            self.mobile(atom) as i32
        };
        Some(EndpointInfo {
            role,
            mobile,
            neutral_bonds_valence: valence - mobile,
            valence,
        })
    }

    /// The flag predicates' full bodies (the reference's per-flag
    /// `nGetEndpointInfo_*` functions, ichitaut.c 452–1080).
    fn endpoint_info_ext(
        &self,
        state: &Registration,
        atom: usize,
        pred: Predicate,
    ) -> Option<EndpointInfo> {
        let record = &self.structure.atoms[atom];
        let element = record.element.as_str();
        let valence = match pred {
            Predicate::Ket | Predicate::Pt16 => endpoint_valence_ket(element),
            Predicate::Pt22 => (element == "C").then_some(4),
            Predicate::Pt06 => match element {
                "C" => Some(4),
                _ => endpoint_valence(element),
            },
            Predicate::Pt39 => match element {
                "C" => Some(4),
                "N" => Some(3),
                _ => None,
            },
            Predicate::Pt13 => match element {
                "C" => Some(4),
                "O" | "S" | "Se" | "Te" => Some(2),
                _ => None,
            },
            Predicate::Pt18 => match element {
                "O" => Some(2),
                "N" => Some(3),
                _ => None,
            },
            Predicate::Base => endpoint_valence(element),
        }? as i32;
        let heavy = self.neighbors[atom].len() as i32;
        // PT_18_00 carries NO connectivity checks at all — any O/N
        // with charge 0/-1 qualifies (ichitaut.c 832–905).
        if pred == Predicate::Pt18 {
            let raw = self.hydrogens[atom] as i32;
            let role = if raw > 0 { Role::Donor } else { Role::Acceptor };
            return Some(EndpointInfo {
                role,
                mobile: raw,
                neutral_bonds_valence: valence - raw,
                valence,
            });
        }
        // The connectivity bands: below the endpoint valence, and the
        // KET/PT_16_00 extra bands on the C/O valences (the "exclude"
        // checks, ichitaut.c 941–951).
        if valence <= heavy {
            return None;
        }
        if matches!(pred, Predicate::Ket | Predicate::Pt16) {
            if valence == 4 && heavy < 2 {
                return None;
            }
            if valence == 2 && heavy > 1 {
                return None;
            }
        }
        if record.charge > 0 {
            // The charge-point branch: same c-point arithmetic as the
            // base predicate.
            if state.c_point_of[atom].is_some() {
                if let Some((_, subtype)) = charge_type(self, state, atom) {
                    let raw = self.hydrogens[atom] as i32;
                    if subtype & 4 != 0 {
                        return Some(EndpointInfo {
                            role: Role::Acceptor,
                            mobile: raw,
                            neutral_bonds_valence: valence - raw,
                            valence,
                        });
                    }
                    if subtype & 8 != 0 {
                        return Some(EndpointInfo {
                            role: Role::Donor,
                            mobile: raw,
                            neutral_bonds_valence: valence - raw,
                            valence,
                        });
                    }
                }
            }
            return None;
        }
        let mobile = self.mobile(atom) as i32;
        if 2 * mobile + self.bond_sum2[atom] != 2 * valence {
            return None;
        }
        let role = match pred {
            // PT_13_00 reads the mobile units for the role, not the
            // bond arithmetic (ichitaut.c 815–823).
            Predicate::Pt13 => {
                if mobile > 0 {
                    Role::Donor
                } else {
                    Role::Acceptor
                }
            }
            _ => match self.bond_sum2[atom] - 2 * heavy {
                0 => Role::Donor,
                2 => Role::Acceptor,
                _ => return None,
            },
        };
        Some(EndpointInfo {
            role,
            mobile,
            neutral_bonds_valence: valence - mobile,
            valence,
        })
    }

    /// The BNS build's donor wiring (the reference's AddTGroups2BnStruct,
    /// ichi_bns.c 9296–9327): the endpoint info for a registered atom
    /// comes from the FIRST predicate in the flag chain that
    /// classifies it — the base one, then the enabled flags in the
    /// C's order (22, 16, 06, 39, 13, 18, KET). With every flag off
    /// this is the base verdict alone.
    fn bns_role(&self, state: &Registration, atom: usize) -> Option<Role> {
        let flags = self.structure.taut_ext;
        let classify = |pred| self.endpoint_by(state, atom, pred);
        if let Some(role) = classify(Predicate::Base) {
            return Some(role);
        }
        for (enabled, pred) in [
            (flags.pt_22_00, Predicate::Pt22),
            (flags.pt_16_00, Predicate::Pt16),
            (flags.pt_06_00, Predicate::Pt06),
            (flags.pt_39_00, Predicate::Pt39),
            (flags.pt_13_00, Predicate::Pt13),
            (flags.pt_18_00, Predicate::Pt18),
            (flags.keto_enol, Predicate::Ket),
        ] {
            if enabled {
                if let Some(role) = classify(pred) {
                    return Some(role);
                }
            }
        }
        None
    }

    /// Mobile units the reference's `nMobile = num_H + (charge == -1)`.
    /// Raw structures read the legacy riding-hydrogen model (the
    /// stored count already carries the anion's proton); a fully
    /// pinned structure — the (DE)PROTONATION-standardized state —
    /// counts the deprotonated endpoint's charge instead.
    fn mobile(&self, atom: usize) -> u32 {
        self.hydrogens[atom] + self.charge_units[atom]
    }

    /// The reference reads `bond_type & ~BOND_MARK_ALL`, and
    /// SetTautomericBonds REPLACES the marked bond's type with
    /// BOND_TAUTOM — a registered bond reads TAUTOM, its original
    /// order gone. A perceived-aromatic bond reads BOND_ALTERN. A
    /// ride-moved bond reads the mobility mark the ride persisted
    /// (ALT12NS "1 or 2" / ALTERN|TRIPLE "2 or 3") until a
    /// registration replaces it.
    fn bond_type(&self, state: &Registration, a: usize, b: usize) -> BondType {
        let &index = self
            .bond_index
            .get(&(a.min(b), a.max(b)))
            .expect("template walk on a bonded pair");
        if state.marks.contains(&index) {
            return BondType::Tautom;
        }
        if state.altered.contains(&index) || self.aromatic.contains(&index) {
            return BondType::Altern;
        }
        match self.ride_marks[index] {
            9 => return BondType::Alt12ns,
            7 => return BondType::Alt23,
            _ => {}
        }
        match self.structure.bonds[index].2 {
            2 => BondType::Double,
            3 => BondType::Triple,
            4 => BondType::Altern,
            _ => BondType::Single,
        }
    }

    /// The reference's `bIsCenterPointStrict`: all-single with an
    /// endpoint valence above it and a mobile H (or an anionic site),
    /// or exactly one double-bond worth on a strict center element.
    /// The reference's `c_point` disjunct is a charged-BNS context
    /// flag our scope never sets.
    fn is_centerpoint_strict(&self, atom: usize) -> bool {
        let record = &self.structure.atoms[atom];
        let heavy = self.neighbors[atom].len() as i32;
        if self.bond_sum2[atom] == 2 * heavy {
            let ev = endpoint_valence(&record.element).unwrap_or(0) as i32;
            return ev != 0 && ev > heavy && (self.mobile(atom) > 0 || record.charge == -1);
        }
        if self.bond_sum2[atom] == 2 * (heavy + 1) {
            return matches!(record.element.as_str(), "C" | "N" | "P" | "As" | "Sb");
        }
        false
    }

    /// The reference's `is_centerpoint_elem`: elements allowed to
    /// carry the 1,3 template's π bond.
    fn is_centerpoint_elem(&self, atom: usize) -> bool {
        matches!(
            self.structure.atoms[atom].element.as_str(),
            "C" | "N" | "P" | "S" | "I" | "As" | "Sb" | "Se" | "Te" | "Cl" | "Br"
        )
    }

    /// chem_bonds_valence over the input Kekulé orders (doubled sums
    /// halve exactly for single/double/triple input).
    fn max_flow(&self, atom: usize) -> i32 {
        (self.bond_sum2[atom] / 2 - self.neighbors[atom].len() as i32).max(0)
    }

    /// The charge-mode alt path (ALT_PATH_MODE_CHARGE): the moved
    /// atom-atom bonds come back as bond indices for the altered
    /// marks.
    fn charge_alt_path(
        &self,
        state: &Registration,
        charged: usize,
        neutral: usize,
    ) -> (bool, Vec<usize>) {
        let n = self.structure.atoms.len();
        let bond_adjacency: Vec<Vec<(usize, i32, bool)>> = (0..n)
            .map(|atom| {
                self.neighbors[atom]
                    .iter()
                    .map(|&(neigh, kind)| {
                        (
                            neigh,
                            match kind {
                                2 => 1,
                                3 => 2,
                                _ => 0,
                            },
                            !self.allowed_edge(atom, neigh),
                        )
                    })
                    .collect()
            })
            .collect();
        let single_like: Vec<Vec<bool>> = (0..n)
            .map(|atom| {
                self.neighbors[atom]
                    .iter()
                    .map(|&(neigh, _)| {
                        matches!(
                            self.bond_type(state, atom, neigh),
                            BondType::Single
                                | BondType::Altern
                                | BondType::Tautom
                                | BondType::Alt12ns
                        )
                    })
                    .collect()
            })
            .collect();
        let max_flow: Vec<i32> = (0..n).map(|atom| self.max_flow(atom)).collect();
        let donor: Vec<bool> = (0..n)
            .map(|atom| self.bns_role(state, atom) == Some(Role::Donor))
            .collect();
        let net = super::bns::Network::build(&super::bns::Build {
            bond_adjacency: &bond_adjacency,
            max_flow: &max_flow,
            group_of: &state.group_of,
            donor: &donor,
            single_like: &single_like,
            c_point_of: &state.c_point_of,
            c_charged: &state.c_charged,
        });
        let (found, moved) = super::bns::charge_alt_path(&net, charged, neutral);
        let bonds = moved
            .iter()
            .filter_map(|&(a, b)| self.bond_index.get(&(a, b)).copied())
            .collect();
        (found, bonds)
    }

    /// The reference's RegisterCPoints: pair two ungrouped points or
    /// merge/extend their c-groups.
    fn register_c_points(&self, state: &mut Registration, point1: usize, point2: usize) -> bool {
        match (state.c_point_of[point1], state.c_point_of[point2]) {
            (Some(a), Some(b)) if a == b => false,
            (Some(a), Some(b)) => {
                let (keep, gone) = if a < b { (a, b) } else { (b, a) };
                for slot in state.c_point_of.iter_mut() {
                    if *slot == Some(gone) {
                        *slot = Some(keep);
                    }
                }
                true
            }
            (Some(group), None) => {
                state.c_point_of[point2] = Some(group);
                true
            }
            (None, Some(group)) => {
                state.c_point_of[point1] = Some(group);
                true
            }
            (None, None) => {
                let id = state.c_point_of.iter().flatten().max().map_or(0, |g| g + 1);
                state.c_point_of[point1] = Some(id);
                state.c_point_of[point2] = Some(id);
                state.c_charged[point1] = true;
                true
            }
        }
    }

    /// The 2-H salt pair test over a fresh network: donors push,
    /// acceptors pull, grouped atoms act through their t-group vertex.
    fn salt2_alt_path(
        &self,
        state: &Registration,
        first: usize,
        second: usize,
        donors: bool,
    ) -> bool {
        let mut net = self.build_network(state);
        super::bns::exists_2h(&mut net, first, second, donors)
    }

    fn build_network(&self, state: &Registration) -> super::bns::Network {
        let n = self.structure.atoms.len();
        let bond_adjacency: Vec<Vec<(usize, i32, bool)>> = (0..n)
            .map(|atom| {
                self.neighbors[atom]
                    .iter()
                    .map(|&(neigh, kind)| {
                        (
                            neigh,
                            match kind {
                                2 => 1,
                                3 => 2,
                                _ => 0,
                            },
                            !self.allowed_edge(atom, neigh),
                        )
                    })
                    .collect()
            })
            .collect();
        let single_like: Vec<Vec<bool>> = (0..n)
            .map(|atom| {
                self.neighbors[atom]
                    .iter()
                    .map(|&(neigh, _)| {
                        matches!(
                            self.bond_type(state, atom, neigh),
                            BondType::Single
                                | BondType::Altern
                                | BondType::Tautom
                                | BondType::Alt12ns
                        )
                    })
                    .collect()
            })
            .collect();
        let max_flow: Vec<i32> = (0..n).map(|atom| self.max_flow(atom)).collect();
        let donor: Vec<bool> = (0..n)
            .map(|atom| self.bns_role(state, atom) == Some(Role::Donor))
            .collect();
        super::bns::Network::build(&super::bns::Build {
            bond_adjacency: &bond_adjacency,
            max_flow: &max_flow,
            group_of: &state.group_of,
            donor: &donor,
            single_like: &single_like,
            c_point_of: &state.c_point_of,
            c_charged: &state.c_charged,
        })
    }

    /// bExistsAnyAltPath over a fresh network for the current
    /// registration state: every template pairing is confirmed by the
    /// cap-aware balanced search before it registers. `pred` is the
    /// confirming round's mode (the reference's ALT_PATH_MODE_*): the
    /// ungrouped terminals' donor/acceptor precheck consults that
    /// round's predicate (ichi_bns.c 8403–8533).
    fn alt_path_exists(
        &self,
        state: &Registration,
        first: usize,
        second: usize,
        pred: Predicate,
    ) -> bool {
        let result = self.alt_path_exists_traced(state, first, second, pred);
        if std::env::var_os("INCHI_TRACE_BNS").is_some() {
            eprintln!("BNS {first}-{second} -> {result}");
        }
        result
    }

    fn alt_path_exists_traced(
        &self,
        state: &Registration,
        first: usize,
        second: usize,
        pred: Predicate,
    ) -> bool {
        let n = self.structure.atoms.len();
        let bond_adjacency: Vec<Vec<(usize, i32, bool)>> = (0..n)
            .map(|atom| {
                self.neighbors[atom]
                    .iter()
                    .map(|&(neigh, kind)| {
                        (
                            neigh,
                            match kind {
                                2 => 1,
                                3 => 2,
                                _ => 0,
                            },
                            !self.allowed_edge(atom, neigh),
                        )
                    })
                    .collect()
            })
            .collect();
        let single_like: Vec<Vec<bool>> = (0..n)
            .map(|atom| {
                self.neighbors[atom]
                    .iter()
                    .map(|&(neigh, _)| {
                        matches!(
                            self.bond_type(state, atom, neigh),
                            BondType::Single
                                | BondType::Altern
                                | BondType::Tautom
                                | BondType::Alt12ns
                        )
                    })
                    .collect()
            })
            .collect();
        let max_flow: Vec<i32> = (0..n).map(|atom| self.max_flow(atom)).collect();
        let donor: Vec<bool> = (0..n)
            .map(|atom| self.bns_role(state, atom) == Some(Role::Donor))
            .collect();
        let net = super::bns::Network::build(&super::bns::Build {
            bond_adjacency: &bond_adjacency,
            max_flow: &max_flow,
            group_of: &state.group_of,
            donor: &donor,
            single_like: &single_like,
            c_point_of: &state.c_point_of,
            c_charged: &state.c_charged,
        });
        let roles = super::bns::Roles {
            donor_ok: (0..n)
                .map(|atom| {
                    state.registered(atom)
                        || self.endpoint_by(state, atom, pred) == Some(Role::Donor)
                })
                .collect(),
            acceptor_ok: (0..n)
                .map(|atom| {
                    state.registered(atom)
                        || self.endpoint_by(state, atom, pred) == Some(Role::Acceptor)
                })
                .collect(),
        };
        super::bns::any_alt_path(&net, &roles, first, second)
    }
}

/// The atoms whose hydrogens are mobile: members of a mobile group.
/// Endpoint arithmetic alone is not enough — an alcohol oxygen is a
/// donor candidate that never joins a group, and its hydrogen must
/// keep counting in the canonicalization.
/// The mobile charge units per atom: a deprotonated endpoint counts
/// as a mobile unit only on a fully pinned structure — the
/// (DE)PROTONATION-standardized state. Raw structures read the legacy
/// riding-hydrogen model.
fn charge_units_of(structure: &Structure) -> impl Fn(usize) -> u32 + '_ {
    let pinned = structure.atoms.iter().all(|atom| atom.h_pinned);
    move |atom: usize| (pinned && structure.atoms[atom].charge == -1) as u32
}

pub fn endpoint_mask(structure: &Structure, hydrogens: &[u32]) -> Vec<bool> {
    let context = Context::new(structure, hydrogens);
    let state = fresh_state(structure.atoms.len());
    let charge_units = charge_units_of(structure);
    let mut mask = vec![false; structure.atoms.len()];
    for members in components(&context, &state, hydrogens, &charge_units) {
        for atom in members {
            mask[atom] = true;
        }
    }
    mask
}

/// The mobile groups rank-free, as (hydrogen count, member atoms) —
/// the terminal-compare rows inside the canonicalization search (the
/// reference's t-group CT section).
pub fn group_membership(structure: &Structure, hydrogens: &[u32]) -> Vec<(u32, Vec<usize>)> {
    let context = Context::new(structure, hydrogens);
    let state = fresh_state(structure.atoms.len());
    let charge_units = charge_units_of(structure);
    components(&context, &state, hydrogens, &charge_units)
        .into_iter()
        .map(|members| {
            let count = members.iter().map(|&m| hydrogens[m]).sum();
            (count, members)
        })
        .collect()
}

/// Stage-2 hydrogen counts: the given counts with mobile H removed —
/// the reference's `numHNoTautH`, which both the `/c` and `/h`
/// canonicalizations rank by.
pub fn without_mobile_h(hydrogens: &[u32], mask: &[bool]) -> Vec<u32> {
    hydrogens
        .iter()
        .zip(mask)
        .map(|(&count, &endpoint)| if endpoint { 0 } else { count })
        .collect()
}

/// The registration state the templates mutate: group membership per
/// atom (the reference's `at[].endpoint`) and the bonds marked
/// tautomeric (`BOND_TAUTOM`). The marks are the transitive device —
/// registering a group marks all its anchors' bonds, and a
/// tautomeric bond lets its endpoint qualify as either role in
/// later scans.
#[derive(Default)]
struct Registration {
    group_of: Vec<Option<usize>>,
    marks: BTreeSet<usize>,
    /// The charge-group membership (at[].c_point) and which member
    /// carries the moveable charge.
    c_point_of: Vec<Option<usize>>,
    c_charged: Vec<bool>,
    /// Bonds a successful charge move alternated (BOND_ALTERED): the
    /// templates read them as alternating.
    altered: BTreeSet<usize>,
    /// A salt registration succeeded (TG_FLAG_FOUND_SALT_CHARGES_DONE):
    /// later salt gates relax to donor-only.
    salt_charges_found: bool,
    /// A HARD (DE)PROTONATION phase moved a proton (the structure's
    /// FLAG_FORCE_SALT_TAUT): the salt gates relax to donor-only.
    force_salt_taut: bool,
}

impl Registration {
    fn registered(&self, atom: usize) -> bool {
        self.group_of[atom].is_some()
    }

    /// Merges the members into one group; `true` when the state
    /// changed (a member was unregistered or two groups fused).
    fn register(&mut self, members: &[usize]) -> bool {
        let existing: Vec<usize> = members.iter().filter_map(|&m| self.group_of[m]).collect();
        match existing.len() {
            0 => {
                let id = self.group_of.iter().flatten().max().map_or(0, |m| m + 1);
                for &m in members {
                    self.group_of[m] = Some(id);
                }
                true
            }
            1 => {
                let id = existing[0];
                let mut changed = false;
                for &m in members {
                    if self.group_of[m] != Some(id) {
                        self.group_of[m] = Some(id);
                        changed = true;
                    }
                }
                changed
            }
            _ => {
                let keep = *existing.iter().min().unwrap();
                let absorbed: BTreeSet<usize> = existing.into_iter().collect();
                let mut changed = false;
                // Absorb every member of the fused groups, then add
                // the new members themselves.
                for group in self.group_of.iter_mut() {
                    if let Some(id) = *group {
                        if absorbed.contains(&id) && id != keep {
                            *group = Some(keep);
                            changed = true;
                        }
                    }
                }
                for &m in members {
                    if self.group_of[m] != Some(keep) {
                        self.group_of[m] = Some(keep);
                        changed = true;
                    }
                }
                changed
            }
        }
    }
}

/// A blank registration state for the stateless public entry points.
fn fresh_state(n: usize) -> Registration {
    Registration {
        group_of: vec![None; n],
        marks: BTreeSet::new(),
        c_point_of: vec![None; n],
        c_charged: vec![false; n],
        altered: BTreeSet::new(),
        salt_charges_found: false,
        force_salt_taut: false,
    }
}

/// Builds the endpoint linkage: the templates run as registration
/// rounds — anchors, then the shift walks — each registering groups
/// and re-marking bonds, to a fixpoint. The returned adjacency joins
/// every pair of co-members.
fn build_links(context: &Context) -> Vec<BTreeSet<usize>> {
    let n = context.structure.atoms.len();
    let mut state = Registration {
        group_of: vec![None; n],
        marks: BTreeSet::new(),
        c_point_of: vec![None; n],
        c_charged: vec![false; n],
        altered: BTreeSet::new(),
        salt_charges_found: false,
        force_salt_taut: context.structure.protons_hard_moved,
    };
    let mut alternated = alternation_round(context, &mut state);
    for pass in 0..8 {
        let mut changed = registration_round(context, &mut state);
        let (r, mut f, mut w, mut c, mut s0, mut s1, mut a) =
            (changed, false, false, false, false, false, false);
        // The flag passes sit between the base anchors and the walks
        // (the C's block order in MarkTautomerGroups); the walks stay
        // behind the whole call's no-change gate (`!tot_changes`).
        if !changed {
            changed = flag_round(context, &mut state);
            f = changed;
        }
        if !changed {
            changed = walk_round(context, &mut state);
            w = changed;
        }
        if !changed {
            changed = charge_round(context, &mut state);
            c = changed;
        }
        if !changed {
            changed = salt_round(context, &mut state);
            s0 = changed;
        }
        if !changed {
            changed = salt2_round(context, &mut state);
            s1 = changed;
        }
        if !changed {
            changed = alternation_round(context, &mut state);
            a = changed;
        }
        if std::env::var_os("INCHI_TRACE_BNS").is_some() {
            let groups: Vec<usize> = state.group_of.iter().flatten().copied().collect();
            eprintln!(
                "FIXPOINT pass={pass} reg={r} flags={f} walk={w} charge={c} salt0={s0} salt2={s1} merge=false alt={a} groups={} marks={} altered={}",
                groups.len(),
                state.marks.len(),
                state.altered.len()
            );
        }
        if !changed && !alternated {
            break;
        }
        alternated = false;
    }
    // The salt merge is terminal (the reference's MergeSaltTautGroups
    // runs after the whole master loop): its co-membership never
    // feeds back into the charge rounds — a c-point the charge round
    // registers against a merged group is a c-point the reference
    // never registers.
    merge_salt_round(context, &mut state);
    // Adjacency over co-membership (components() then yields the
    // groups; empty ones drop).
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (atom, group) in state.group_of.iter().enumerate() {
        if let Some(group) = group {
            members[*group].push(atom);
        }
    }
    let mut linked: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    for group in members {
        for &x in &group {
            for &y in &group {
                if x != y {
                    linked[x].insert(y);
                    linked[y].insert(x);
                }
            }
        }
    }
    linked
}

/// The reference's GetChargeType + bCanBeACPoint: (type, subtype)
/// for the charge-point candidates — N (0) and P (1) always, plus
/// the ring-bound valence-2 O/S/Se/Te; the subtype bits say which
/// side of a charge move the atom can play.
fn charge_type(context: &Context, state: &Registration, atom: usize) -> Option<(i32, u8)> {
    const CHARGED: u8 = 0;
    const P_DONOR: u8 = 1;
    const P_ACCEPT: u8 = 2;
    const H_ACCEPT: u8 = 4;
    const H_DONOR: u8 = 8;
    const NEUTRAL: u8 = 16;

    let record = &context.structure.atoms[atom];
    // Ion pairs are ignored; a neighbor of the opposite charge that
    // is not a tautomeric endpoint disqualifies.
    if record.charge.abs() == 1 {
        for &(neigh, _) in &context.neighbors[atom] {
            let other = context.structure.atoms[neigh].charge;
            if (other + record.charge).abs() < (other - record.charge).abs()
                && state.group_of[neigh].is_none()
            {
                return None;
            }
        }
    } else if record.charge != 0 {
        return None;
    }

    let element = record.element.as_str();
    let ring_bound = context.rings.size[atom] >= 5 && context.neighbors[atom].len() == 2;
    let ctype: Option<(i32, i32, i32, i32)> = match element {
        "N" => Some((0, 3, 3, 1)),
        "P" => Some((1, 3, 3, 1)),
        "O" if ring_bound => Some((2, 2, 2, 1)),
        "S" if ring_bound => Some((3, 2, 2, 1)),
        "Se" if ring_bound => Some((4, 2, 2, 1)),
        "Te" if ring_bound => Some((5, 2, 2, 1)),
        _ => None,
    };
    let (ctype_id, neutral_valence, neutral_cbv, change_valence) = ctype?;
    let valence = context.neighbors[atom].len() as i32;
    let cbv = context.bond_sum2[atom] / 2;
    let num_h = context.mobile(atom) as i32;
    let endpoint_valence = endpoint_valence(element).unwrap_or(0) as i32;
    let neg = (record.charge == -1) as i32;

    let subtype: u8;
    if record.charge == 1 && valence == cbv && num_h > 0 {
        // Proton-donor candidates (>NH(+)-, >OH(+)): never a point.
        return None;
    } else if record.charge == 1 && valence < cbv {
        let bonds_valence = cbv + num_h;
        let num_bonds = valence + num_h;
        if bonds_valence != neutral_cbv + change_valence {
            return None;
        }
        if num_bonds == neutral_valence {
            subtype = if valence == neutral_valence || endpoint_valence == 0 {
                CHARGED
            } else {
                CHARGED | H_DONOR
            };
        } else if num_bonds == neutral_valence - 1 {
            subtype = if endpoint_valence > 0 {
                if num_h > 0 {
                    CHARGED | H_ACCEPT | P_DONOR
                } else {
                    CHARGED | H_ACCEPT
                }
            } else {
                CHARGED
            };
        } else {
            return None;
        }
    } else if record.charge == 0 || record.charge == -1 {
        let bonds_valence = cbv + num_h + neg;
        let num_bonds = valence + num_h + neg;
        if bonds_valence != neutral_cbv {
            return None;
        }
        if num_bonds == neutral_valence {
            subtype = if valence == neutral_valence || endpoint_valence == 0 {
                NEUTRAL
            } else if valence < neutral_valence {
                if num_h > 0 {
                    NEUTRAL | H_DONOR
                } else {
                    NEUTRAL | H_ACCEPT
                }
            } else {
                return None;
            };
        } else if num_bonds == neutral_valence - 1 {
            subtype = if endpoint_valence > 0 {
                NEUTRAL | H_ACCEPT | P_ACCEPT
            } else {
                NEUTRAL
            };
        } else {
            return None;
        }
    } else {
        return None;
    }
    Some((ctype_id, subtype))
}

/// One charge-group round, the reference's MarkChargeGroups: pair
/// each charged candidate against the neutral candidates of the same
/// type, confirm by the charge-mode balanced network, register the
/// c-points, and record the bonds the move alternated.
fn charge_round(context: &Context, state: &mut Registration) -> bool {
    let n = context.structure.atoms.len();
    let mut candidates: Vec<(usize, i32, u8)> = (0..n)
        .filter_map(|atom| {
            charge_type(context, state, atom).map(|(ctype, subtype)| (atom, ctype, subtype))
        })
        .collect();
    if candidates.len() <= 1 {
        return false;
    }
    candidates.sort_by_key(|&(atom, ctype, subtype)| (ctype, subtype, atom));
    let neutral = |subtype: u8| subtype & 16 != 0;
    let mut changed = false;
    let mut i1 = 0;
    while i1 < candidates.len() {
        while i1 < candidates.len() && neutral(candidates[i1].2) {
            i1 += 1;
        }
        if i1 == candidates.len() {
            break;
        }
        let mut i2 = i1 + 1;
        while i2 < candidates.len()
            && candidates[i2].1 == candidates[i1].1
            && !neutral(candidates[i2].2)
        {
            i2 += 1;
        }
        if i2 == candidates.len() {
            break;
        }
        let mut i3 = i2;
        while i3 < candidates.len() && candidates[i3].1 == candidates[i1].1 {
            i3 += 1;
        }
        if i3 == i2 {
            if i2 < candidates.len() {
                i1 = i3;
                continue;
            }
            break;
        }
        for charged in &candidates[i1..i2] {
            for neutral_candidate in &candidates[i2..i3] {
                let (iat1, _, _) = *charged;
                let (iat2, _, _) = *neutral_candidate;
                if state.c_point_of[iat1].is_some()
                    && state.c_point_of[iat1] == state.c_point_of[iat2]
                {
                    continue;
                }
                let (found, moved) = context.charge_alt_path(state, iat1, iat2);
                if !found {
                    continue;
                }
                for bond in moved {
                    if state.altered.insert(bond) {
                        changed = true;
                    }
                }
                if context.register_c_points(state, iat1, iat2) {
                    changed = true;
                }
                return changed;
            }
        }
        i1 = i3;
    }
    changed
}

/// One BnsTestAndMarkAltBonds pass, the reference's initial bonds
/// normalization: probe every allowed bond at its other achievable
/// flows; when the balanced network can rebalance around the pinned
/// alternative (2*sum == dots, or the trivial dots==0 pin), the bond
/// and everything the probe moved read BOND_ALTERN from then on —
/// wildcard in the pattern arithmetic, either-role in the anchor
/// scan — while the valence arithmetic keeps the input orders.
fn alternation_round(context: &Context, state: &mut Registration) -> bool {
    let n = context.structure.atoms.len();
    let single_like: Vec<Vec<bool>> = (0..n)
        .map(|atom| {
            context.neighbors[atom]
                .iter()
                .map(|&(neigh, _)| {
                    matches!(
                        context.bond_type(state, atom, neigh),
                        BondType::Single | BondType::Altern | BondType::Tautom | BondType::Alt12ns
                    )
                })
                .collect()
        })
        .collect();
    let bond_adjacency: Vec<Vec<(usize, i32, bool)>> = (0..n)
        .map(|atom| {
            context.neighbors[atom]
                .iter()
                .map(|&(neigh, kind)| {
                    (
                        neigh,
                        match kind {
                            2 => 1,
                            3 => 2,
                            _ => 0,
                        },
                        !context.allowed_edge(atom, neigh),
                    )
                })
                .collect()
        })
        .collect();
    let max_flow: Vec<i32> = (0..n).map(|atom| context.max_flow(atom)).collect();
    let donor: Vec<bool> = (0..n)
        .map(|atom| context.bns_role(state, atom) == Some(Role::Donor))
        .collect();
    let net = super::bns::Network::build(&super::bns::Build {
        bond_adjacency: &bond_adjacency,
        max_flow: &max_flow,
        group_of: &state.group_of,
        donor: &donor,
        single_like: &single_like,
        c_point_of: &state.c_point_of,
        c_charged: &state.c_charged,
    });
    let mut changed = false;
    for atom in 0..n {
        for &(neigh, _) in &context.neighbors[atom] {
            if neigh <= atom || !context.allowed_edge(atom, neigh) {
                continue;
            }
            let &index = context.bond_index.get(&(atom, neigh)).expect("bonded pair");
            if state.altered.contains(&index) {
                continue;
            }
            let Some(iedge) = net.edge_between(atom, neigh) else {
                continue;
            };
            let cur = net.cur_flow(iedge);
            let cap = net.edge_cap(iedge);
            let min = net.min_flow_2_check(iedge);
            for test in min..=cap {
                if test == cur {
                    continue;
                }
                let probe = super::bns::probe_altern_flow(&net, iedge, test);
                if std::env::var_os("INCHI_TRACE_BNS").is_some() {
                    eprintln!("ALT {atom}-{neigh} cur={cur} test={test} min={min} cap={cap} -> feasible={} moved={:?}", probe.feasible, probe.moved);
                }
                if probe.feasible {
                    state.altered.insert(index);
                    for &(a, b) in &probe.moved {
                        if let Some(&moved) = context.bond_index.get(&(a, b)) {
                            state.altered.insert(moved);
                        }
                    }
                    changed = true;
                    break;
                }
            }
        }
    }
    changed
}

/// The salt-subtype bits (incomdef.h:97): H-donor, (-)-donor,
/// acceptor.
const SALT_DONOR_H: u8 = 1;
const SALT_DONOR_NEG: u8 = 2;
const SALT_ACCEPTOR: u8 = 4;

/// The reference's GetSaltChargeType: a terminal O/S/Se/Te of standard
/// valence on an unsaturated, uncharged carbon (SALT_WITH_PROTONS=1
/// allows =C(H)-OH / -C(H)=O). Subtypes: H-donor, (-)-donor, the C=O
/// acceptor; t-group members inherit group-level subtypes.
fn salt_charge_type(context: &Context, state: &Registration, atom: usize) -> Option<u8> {
    let record = &context.structure.atoms[atom];
    if context.neighbors[atom].len() != 1 || record.charge < -1 || record.charge > 0 {
        return None;
    }
    if !matches!(record.element.as_str(), "O" | "S" | "Se" | "Te") {
        return None;
    }
    let valence_needed: i32 = match (record.element.as_str(), record.charge) {
        ("O", 0) | ("S", 0) | ("Se", 0) | ("Te", 0) => 2,
        (_, -1) => 1,
        _ => return None,
    };
    // The reference's chem_bonds_valence + num_H — the raw hydrogen
    // count, NOT the mobile units: on a pinned (charged) structure a
    // deprotonated site's charge is not a hydrogen and must not
    // inflate the valence fill.
    if context.bond_sum2[atom] / 2 + context.hydrogens[atom] as i32 != valence_needed {
        return None;
    }
    let (carbon, _) = context.neighbors[atom][0];
    let carbon_record = &context.structure.atoms[carbon];
    if carbon_record.element != "C"
        || carbon_record.charge != 0
        || context.bond_sum2[carbon] / 2 + context.hydrogens[carbon] as i32 != 4
        || context.bond_sum2[carbon] / 2 == context.neighbors[carbon].len() as i32
    {
        return None;
    }
    let mut subtype = 0;
    if let Some(group) = state.group_of[atom] {
        let (mobile, negative) = group_h_and_neg(context, state, group);
        if mobile > negative {
            subtype |= SALT_DONOR_H;
        }
        if negative > 0 {
            subtype |= SALT_DONOR_NEG;
        }
        return Some(subtype | SALT_ACCEPTOR);
    }
    if record.charge == -1 {
        subtype |= SALT_DONOR_NEG;
    }
    if record.charge <= 0 && context.hydrogens[atom] > 0 {
        subtype |= SALT_DONOR_H;
    }
    if record.charge == 0 && context.bond_sum2[atom] / 2 == 2 {
        subtype |= SALT_ACCEPTOR;
    }
    Some(subtype)
}

/// The reference's GetOtherSaltChargeType (bAccept_O=1): a tautomeric
/// endpoint adjacent to a centerpoint-capable neighbor through a bond
/// matching the atom's role.
fn other_salt_charge_type(context: &Context, state: &Registration, atom: usize) -> Option<u8> {
    let role = context.endpoint(state, atom)?;
    let mut subtype = 0;
    let mut centerpoints = 0;
    for &(centerpoint, _kind) in &context.neighbors[atom] {
        let bond = context.bond_type(state, atom, centerpoint);
        let bond_allows = match role {
            Role::Acceptor => {
                matches!(
                    bond,
                    BondType::Double | BondType::Altern | BondType::Tautom | BondType::Alt12ns
                )
            }
            Role::Donor => matches!(
                bond,
                BondType::Single | BondType::Altern | BondType::Tautom | BondType::Alt12ns
            ),
        };
        let center_capable = context.bond_sum2[centerpoint] / 2
            > context.neighbors[centerpoint].len() as i32
            || (context.bond_sum2[centerpoint] / 2 == context.neighbors[centerpoint].len() as i32
                && (state.group_of[centerpoint].is_some()
                    || state.c_point_of[centerpoint].is_some()));
        if bond_allows && center_capable && context.is_centerpoint_elem(centerpoint) {
            centerpoints += 1;
            break;
        }
    }
    if centerpoints == 0 {
        return None;
    }
    if let Some(group) = state.group_of[atom] {
        let (mobile, negative) = group_h_and_neg(context, state, group);
        if mobile > negative {
            subtype |= SALT_DONOR_H;
        }
        if negative > 0 {
            subtype |= SALT_DONOR_NEG;
        }
        return Some(subtype | SALT_ACCEPTOR);
    }
    if role == Role::Acceptor {
        subtype |= SALT_ACCEPTOR;
    }
    if role == Role::Donor {
        if context.structure.atoms[atom].charge == -1 {
            subtype |= SALT_DONOR_NEG;
        }
        if context.hydrogens[atom] > 0 {
            subtype |= SALT_DONOR_H;
        }
    }
    Some(subtype)
}

/// A t-group's (num[0], num[1]): the mobile H count and the negative
/// charges among its members.
fn group_h_and_neg(context: &Context, state: &Registration, group: usize) -> (i32, i32) {
    let (mut mobile, mut negative) = (0, 0);
    for (atom, membership) in state.group_of.iter().enumerate() {
        if *membership == Some(group) {
            mobile += context.mobile(atom) as i32;
            if context.structure.atoms[atom].charge == -1 {
                negative += 1;
            }
        }
    }
    (mobile, negative)
}

/// One salt round, the reference's MarkSaltChargeGroups (salt step 0):
/// pair each H/(-)-donor candidate against each acceptor candidate in
/// both orders; the balanced network confirms (ALT_PATH_MODE_4_SALT —
/// the same setup and role prechecks as the tautomer check), and the
/// pair registers as a t-group.
fn salt_round(context: &Context, state: &mut Registration) -> bool {
    let n = context.structure.atoms.len();
    let candidates: Vec<(usize, u8)> = (0..n)
        .filter_map(|atom| {
            salt_charge_type(context, state, atom)
                .or_else(|| other_salt_charge_type(context, state, atom))
                .map(|subtype| (atom, subtype))
        })
        .collect();
    if std::env::var_os("INCHI_TRACE_BNS").is_some() {
        eprintln!("SALT cands: {:?}", candidates);
    }
    if candidates.len() <= 1 {
        return false;
    }
    let subtype_all = candidates.iter().fold(0u8, |acc, &(_, s)| acc | s);
    if subtype_all & SALT_ACCEPTOR == 0 {
        return false;
    }
    // CHARGED_SALTS_ONLY: without a negative-charge donor (and unless
    // every candidate is a tautomeric endpoint) the round declines —
    // unless FLAG_FORCE_SALT_TAUT (a hard proton move happened), when
    // any donor class suffices.
    let neg_present = subtype_all & SALT_DONOR_NEG != 0;
    let all_other = candidates
        .iter()
        .all(|&(atom, _)| salt_charge_type(context, state, atom).is_none());
    if std::env::var_os("INCHI_TRACE_BNS").is_some() {
        eprintln!("SALT gate: subtype_all={subtype_all} neg={neg_present} all_other={all_other}");
    }
    if state.force_salt_taut {
        if subtype_all & (SALT_DONOR_NEG | SALT_DONOR_H) == 0 {
            return false;
        }
    } else if !neg_present || all_other {
        return false;
    }
    let mut changed = false;
    for i1 in 0..candidates.len() {
        for i2 in (i1 + 1)..candidates.len() {
            let (j1, j2) = (candidates[i1].0, candidates[i2].0);
            if state.group_of[j1].is_some() && state.group_of[j1] == state.group_of[j2] {
                continue;
            }
            for (donor, acceptor) in [(j1, j2), (j2, j1)] {
                let subtype_of = |atom: usize| {
                    candidates
                        .iter()
                        .find(|&&(a, _)| a == atom)
                        .map(|&(_, s)| s)
                        .unwrap_or(0)
                };
                if subtype_of(donor) & (SALT_DONOR_H | SALT_DONOR_NEG) == 0
                    || subtype_of(acceptor) & SALT_ACCEPTOR == 0
                {
                    continue;
                }
                if context.alt_path_exists(state, acceptor, donor, Predicate::Base) {
                    changed |= state.register(&[donor, acceptor]);
                    break;
                }
            }
        }
    }
    changed
}

/// One salt round, the reference's MarkSaltChargeGroups (salt step 0):
/// phase. Each endpoint reaches a centerpoint through a π-class bond;
/// the centerpoint's other endpoint neighbors are classified by bond
/// (single = donor side with mobile H, double = acceptor, tautomeric
/// or alternating = either role), and the pairing is confirmed by the
/// balanced network before any group registers —
/// FindAccessibleEndPoints partitions the admitted endpoints into
/// equivalence classes, and only classes that bridge or extend groups
/// register (the `!nGroupNumber || bDiffGroups > 0` gate). The
/// surviving centerpoint bonds are marked tautomeric regardless.
/// MergeSaltTautGroups: after the salt steps, the acidic-site
/// candidates' groups merge BLIND — one representative per group, all
/// registered into a single t-group — provided the candidate set
/// holds both a negative-charge donor class and an H donor class (or
/// any groups at all once a salt registration has succeeded). The
/// non-taut acceptor-only sites and non-taut nitrogens never join.
fn merge_salt_round(context: &Context, state: &mut Registration) -> bool {
    let n = context.structure.atoms.len();
    let proton_context = crate::canonical::protons::ProtonContext::new(context.structure);
    const SALT_P_DONOR: u8 = 8;
    const SALT_P_ACCEPTOR: u8 = 16;
    const DONOR_ALL: u8 = SALT_DONOR_NEG | SALT_DONOR_H | SALT_P_ACCEPTOR | SALT_P_DONOR;
    const DONOR_NEG2: u8 = SALT_DONOR_NEG | SALT_P_ACCEPTOR;
    const DONOR_H2: u8 = SALT_DONOR_H | SALT_P_DONOR;

    let mut subtype_all = 0u8;
    let mut candidates: Vec<usize> = Vec::new();
    for atom in 0..n {
        let subtype = if let Some(sub) = salt_charge_type(context, state, atom) {
            sub
        } else if let Some(sub) = other_salt_charge_type(context, state, atom) {
            sub
        } else if crate::canonical::protons::has_acidic_hydrogen(&proton_context, atom) {
            SALT_P_DONOR
        } else if crate::canonical::protons::has_acidic_minus(&proton_context, atom) {
            SALT_P_ACCEPTOR
        } else {
            continue;
        };
        if state.registered(atom) {
            // tautomeric atoms always join
        } else if context.structure.atoms[atom].element == "N" {
            continue; // bDoNotMergeNonTautAtom
        }
        if subtype & DONOR_ALL == 0 {
            continue;
        }
        if subtype & SALT_ACCEPTOR != 0 && !state.registered(atom) {
            continue;
        }
        subtype_all |= subtype;
        candidates.push(atom);
    }

    if std::env::var_os("INCHI_TRACE_BNS").is_some() {
        eprintln!(
            "MERGE cands={candidates:?} suball={subtype_all} found={}",
            state.salt_charges_found
        );
    }

    if candidates.len() <= 1 {
        return false;
    }
    // The forced mode (ALLOW_NO_NEGTV_O, FOUND_SALT_CHARGES_DONE, or
    // FLAG_FORCE_SALT_TAUT — a hard proton move happened) declines
    // only when (a NEG2 and an H2 donor both miss) AND no t-group
    // exists; the normal mode demands both donor classes.
    if state.salt_charges_found || state.force_salt_taut {
        if (subtype_all & DONOR_NEG2 == 0 || subtype_all & DONOR_H2 == 0)
            && state.group_of.iter().flatten().next().is_none()
        {
            return false;
        }
    } else if subtype_all & DONOR_NEG2 == 0 || subtype_all & DONOR_H2 == 0 {
        return false;
    }

    // One representative per group (the C's endpoint-change scan);
    // the acceptor-only non-taut sites already filtered above.
    let mut representatives: Vec<usize> = Vec::new();
    let mut seen: Vec<usize> = Vec::new();
    let mut distinct = 0usize;
    let mut previous: Option<usize> = None;
    for &atom in &candidates {
        let group = state.group_of[atom];
        if group != previous || group.is_none() {
            distinct += 1;
        }
        previous = group;
        match group {
            Some(group) => {
                if !seen.contains(&group) {
                    seen.push(group);
                    representatives.push(atom);
                }
            }
            None => representatives.push(atom),
        }
    }
    if distinct <= 1 {
        return false;
    }
    state.register(&representatives)
}

/// MarkSaltChargeGroups2 (salt step 1): the multi-attachment merge.
/// Candidates are acidic sites — genuine salt sites (type 0) and
/// other-type exchangeable sites (type 1) — the first per group kept,
/// tautomeric candidates ranked first. The TST pass confirms both a
/// donor pair and an acceptor pair exist; the CHG pass re-runs the
/// pair tests and marks every endpoint SALT_SELECTED; all selected
/// atoms' groups merge into one t-group.
fn salt2_round(context: &Context, state: &mut Registration) -> bool {
    let n = context.structure.atoms.len();
    let mut candidates: Vec<(usize, u8, u8)> = Vec::new();
    let mut subtype_all = 0u8;
    let mut other_count = 0usize;
    for atom in 0..n {
        let (kind, subtype) = if let Some(sub) = salt_charge_type(context, state, atom) {
            (0u8, sub)
        } else if let Some(sub) = other_salt_charge_type(context, state, atom) {
            other_count += 1;
            (1u8, sub)
        } else {
            continue;
        };
        subtype_all |= subtype;
        candidates.push((atom, kind, subtype));
    }
    // The gate: acceptor present, and a donor unless everything is
    // other-type; FOUND_SALT_CHARGES_DONE or FLAG_FORCE_SALT_TAUT
    // relaxes to donor-only.
    let donor_all = SALT_DONOR_NEG | SALT_DONOR_H;
    let forced = state.salt_charges_found || state.force_salt_taut;
    if candidates.len() <= 1
        || subtype_all & SALT_ACCEPTOR == 0
        || (forced && subtype_all & donor_all == 0)
        || (!forced && (subtype_all & SALT_DONOR_NEG == 0 || other_count == candidates.len()))
    {
        return false;
    }

    // Keep the first candidate per group; tautomeric candidates first.
    let mut seen: Vec<usize> = Vec::new();
    candidates.retain(|&(atom, _, _)| match state.group_of[atom] {
        Some(group) => {
            if seen.contains(&group) {
                false
            } else {
                seen.push(group);
                true
            }
        }
        None => true,
    });
    candidates.sort_by_key(|&(atom, _, _)| state.group_of[atom].is_none());

    let is_donor = |subtype: u8| subtype & donor_all != 0;
    let is_acceptor = |subtype: u8| subtype & SALT_ACCEPTOR != 0;

    // The TST pass: both pair kinds must exist somewhere.
    let mut donor_pair = false;
    let mut acceptor_pair = false;
    'tst: for i in 0..candidates.len() {
        for j in 0..=i {
            if i == j && state.group_of[candidates[i].0].is_none() {
                continue;
            }
            let (_, _, si) = candidates[i];
            let (_, _, sj) = candidates[j];
            if !donor_pair
                && is_donor(si)
                && is_donor(sj)
                && context.salt2_alt_path(state, candidates[i].0, candidates[j].0, true)
            {
                donor_pair = true;
            }
            if !acceptor_pair
                && is_acceptor(si)
                && is_acceptor(sj)
                && context.salt2_alt_path(state, candidates[i].0, candidates[j].0, false)
            {
                acceptor_pair = true;
            }
            if donor_pair && acceptor_pair {
                break 'tst;
            }
        }
    }
    if !(donor_pair && acceptor_pair) {
        return false;
    }

    // The CHG pass: every confirmed pair marks both endpoints
    // SALT_SELECTED.
    let mut selected: Vec<usize> = Vec::new();
    for i in 0..candidates.len() {
        for j in 0..=i {
            if i == j && state.group_of[candidates[i].0].is_none() {
                continue;
            }
            let (_, _, si) = candidates[i];
            let (_, _, sj) = candidates[j];
            if is_donor(si)
                && is_donor(sj)
                && context.salt2_alt_path(state, candidates[i].0, candidates[j].0, true)
                || is_acceptor(si)
                    && is_acceptor(sj)
                    && context.salt2_alt_path(state, candidates[i].0, candidates[j].0, false)
            {
                for &(atom, _, _) in [candidates[i], candidates[j]].iter() {
                    if !selected.contains(&atom) {
                        selected.push(atom);
                    }
                }
            }
        }
    }
    if selected.is_empty() {
        return false;
    }
    let changed = state.register(&selected);
    if changed {
        state.salt_charges_found = true;
    }
    changed
}

fn registration_round(context: &Context, state: &mut Registration) -> bool {
    let n = context.structure.atoms.len();
    let mut changed = false;
    for start in 0..n {
        if context.endpoint(state, start).is_none() {
            continue;
        }
        for &(centerpoint, _) in &context.neighbors[start] {
            if !context.allowed_edge(start, centerpoint) {
                continue;
            }
            let anchor_bond = context.bond_type(state, start, centerpoint);
            if !matches!(
                anchor_bond,
                BondType::Double | BondType::Altern | BondType::Tautom | BondType::Alt12ns
            ) {
                continue;
            }
            if !context.is_centerpoint_elem(centerpoint) {
                continue;
            }
            let mut ends: Vec<usize> = Vec::new();
            let mut groups: Vec<usize> = Vec::new();
            let mut donors = 0u32;
            let mut acceptors = 0u32;
            let mut possible_mobile = 0u32;
            for &(endpoint, _) in &context.neighbors[centerpoint] {
                if !context.allowed_edge(centerpoint, endpoint) {
                    continue;
                }
                let bond = context.bond_type(state, centerpoint, endpoint);
                let non_taut = match bond {
                    BondType::Single | BondType::Double => true,
                    // "2 or 3" — the ride's mark on a demoted triple —
                    // matches no salt role and no template arm.
                    BondType::Triple | BondType::Alt23 => continue,
                    _ => false,
                };
                let Some(arithmetic) = context.endpoint(state, endpoint) else {
                    continue;
                };
                let registered = state.registered(endpoint);
                let mobile = context.mobile(endpoint) > 0 || registered;
                let (donor_side, acceptor_side) = if non_taut {
                    (bond == BondType::Single && mobile, bond == BondType::Double)
                } else {
                    (
                        registered || arithmetic == Role::Donor,
                        registered || arithmetic == Role::Acceptor,
                    )
                };
                if !donor_side && !acceptor_side {
                    continue;
                }
                donors += donor_side as u32;
                acceptors += acceptor_side as u32;
                possible_mobile += mobile as u32;
                ends.push(endpoint);
                groups.push(state.group_of[endpoint].unwrap_or(0));
            }
            if ends.len() <= 1 || possible_mobile == 0 || donors == 0 || acceptors == 0 {
                continue;
            }
            let partitioned = find_accessible_end_points(context, state, &ends, Predicate::Base);
            changed |= register_partition(context, state, centerpoint, &partitioned);
        }
    }
    changed
}

/// RegisterEndPoints + SetTautomericBonds over the FindAccessibleEndPoints
/// survivors: register when the class list introduces a new endpoint
/// or bridges two groups, then mark the surviving centerpoint bonds
/// tautomeric. Shared by the base anchor scan and the flag passes.
fn register_partition(
    context: &Context,
    state: &mut Registration,
    centerpoint: usize,
    partitioned: &Partitioned,
) -> bool {
    let mut changed = false;
    if partitioned.ends.is_empty() {
        return changed;
    }
    // RegisterEndPoints' gate over the surviving endpoints:
    // register when the class list introduces a new endpoint
    // or bridges two groups.
    let mut sentinel = usize::MAX;
    let mut diff: i32 = -1;
    for &atom in &partitioned.ends {
        let group = state.group_of[atom].unwrap_or(0);
        if sentinel != group {
            diff += 1;
            sentinel = group;
        }
    }
    if sentinel == 0 || diff > 0 {
        // Merge each equivalence class into one group,
        // absorbing every current member of the class's old
        // groups (the renumbering merge).
        let mut classes: Vec<(usize, Vec<usize>)> = Vec::new();
        for (position, &atom) in partitioned.ends.iter().enumerate() {
            let equ = partitioned.equ[position];
            let entry = classes.iter_mut().find(|(label, _)| *label == equ);
            let entry = match entry {
                Some(entry) => entry,
                None => {
                    classes.push((equ, Vec::new()));
                    classes.last_mut().unwrap()
                }
            };
            entry.1.push(atom);
            if let Some(group) = state.group_of[atom] {
                for (other, membership) in state.group_of.iter().enumerate() {
                    if membership == &Some(group) && !entry.1.contains(&other) {
                        entry.1.push(other);
                    }
                }
            }
        }
        for (_, members) in &classes {
            changed |= state.register(members);
        }
    }
    // SetTautomericBonds on the surviving centerpoint bonds.
    for &atom in &partitioned.ends {
        let bond = context.bond_index[&(centerpoint.min(atom), centerpoint.max(atom))];
        changed |= state.marks.insert(bond);
    }
    changed
}

/// The endpoint-scan role arithmetic of one flag family — the C's
/// donor/acceptor counting after each candidate endpoint passes its
/// predicate (ichitaut.c 4640–4660 PT_22_00 and friends, 5308–5331
/// PT_13_00, 5487–5508 PT_18_00, 5688–5706 KET).
#[derive(Clone, Copy)]
enum ScanStyle {
    /// The base arithmetic: a real single bond with mobile units (or
    /// registration) donates, a real double bond accepts; on
    /// alternating/tautomeric bonds a registered atom plays either
    /// role and an unregistered one reads the predicate verdict.
    /// `seed_verdict_alt` is the KET quirk: its alt-class arithmetic
    /// reads the SEED's verdict (`eif1`), not the current endpoint's
    /// (`eif2`) — the reference's variable split, ichitaut.c 5667
    /// vs 5693–5705.
    Standard { seed_verdict_alt: bool },
    /// The PT_13_00/PT_18_00 arithmetic: a real bond's role follows
    /// the mobile units alone (donor when any, acceptor when none);
    /// on alt-class bonds the predicate verdict decides, with
    /// `alt_acceptor_by_neutral` picking PT_13_00's neutral-bonds
    /// test over PT_18_00's plain acceptor check.
    Mobile { alt_acceptor_by_neutral: bool },
}

/// The tail gate's counters: every field `None`/false when the flag
/// block counts nothing (PT_22_00/16/06), otherwise the exact
/// populations the C demands (`num_O == 1` everywhere plus the
/// per-flag extras, ichitaut.c 4693, 4861, 5025, 5194, 5362, 5541,
/// 5736).
#[derive(Clone, Copy, Default)]
struct Gate {
    /// Admitted endpoints of valence 2 == 1 (KET/PT_13_00/PT_18_00).
    valence2_equals: Option<u32>,
    /// Admitted endpoints of valence 4: KET demands at least one,
    /// PT_13_00 exactly one.
    valence4: Option<Valence4Gate>,
    /// Admitted endpoints of valence 3 == 1 (PT_18_00).
    valence3_equals: Option<u32>,
    /// PT_39_00's element counts over EVERY centerpoint neighbor at
    /// scan entry: exactly one O and fewer than two N.
    pt39_element_gate: bool,
}

#[derive(Clone, Copy)]
enum Valence4Gate {
    AtLeastOne,
    ExactlyOne,
}

impl Gate {
    fn allows(&self, admitted: &Counts, all_neighbors: &Counts) -> bool {
        if let Some(want) = self.valence2_equals {
            if admitted.valence2 != want {
                return false;
            }
        }
        match self.valence4 {
            Some(Valence4Gate::AtLeastOne) if admitted.valence4 == 0 => return false,
            Some(Valence4Gate::ExactlyOne) if admitted.valence4 != 1 => return false,
            _ => {}
        }
        if let Some(want) = self.valence3_equals {
            if admitted.valence3 != want {
                return false;
            }
        }
        if self.pt39_element_gate && (all_neighbors.oxygen != 1 || all_neighbors.nitrogen >= 2) {
            return false;
        }
        true
    }
}

/// The scan's running populations.
#[derive(Clone, Copy, Default)]
struct Counts {
    oxygen: u32,
    nitrogen: u32,
    valence2: u32,
    valence3: u32,
    valence4: u32,
}

/// One flag pass's anchor-scan configuration — the per-flag block
/// parameters of MarkTautomerGroups (ichitaut.c 4569–5790).
struct AnchorSpec {
    predicate: Predicate,
    /// The seed→centerpoint bond classes (the C's first bond-type
    /// test). The PT_13_00/PT_18_00 blocks also admit real singles;
    /// the BOND_ALT_13 class (single/triple) never occurs in our
    /// structures — no ride mark produces it.
    seed_bonds: &'static [BondType],
    centerpoint_elements: &'static [&'static str],
    /// The centerpoint's heavy-atom valence when the flag pins it
    /// (PT_16_00: 2, PT_39_00: 3, PT_13_00/PT_18_00: 2).
    centerpoint_valence: Option<usize>,
    /// PT_16_00 demands a neutral centerpoint.
    centerpoint_neutral: bool,
    /// KET's centerpoint: a normal tetravalent carbon — no charge,
    /// chem_bonds_valence + num_H == 4 (ichitaut.c 5615–5619).
    ket_normal_carbon: bool,
    /// PT_16_00's filter: another scan endpoint sharing the seed's
    /// element is skipped (ichitaut.c 4805–4806).
    exclude_seed_element: bool,
    /// PT_18_00's skip: a real double bond onto a valence-3 endpoint
    /// with no mobile units never counts (ichitaut.c 5493–5494).
    pt18_double_skip: bool,
    scan: ScanStyle,
    gate: Gate,
}

/// One flag pass: the reference's per-flag anchor blocks. Seeds off
/// the pass's predicate, walks to centerpoints through the pass's
/// bond classes, rescans the centerpoint's neighbors under the pass's
/// arithmetic, and registers through the shared FindAccessibleEndPoints
/// tail with the pass's BNS mode.
fn anchor_round(context: &Context, state: &mut Registration, spec: &AnchorSpec) -> bool {
    let n = context.structure.atoms.len();
    let mut changed = false;
    // The PT_13_00/PT_18_00 blocks file real triples under the
    // non-taut bonds (where they contribute neither role).
    let wide_bonds = matches!(spec.scan, ScanStyle::Mobile { .. });
    for seed in 0..n {
        let Some(seed_info) = (match spec.predicate {
            Predicate::Base => context.endpoint_info(state, seed),
            _ => context.endpoint_info_ext(state, seed, spec.predicate),
        }) else {
            continue;
        };
        for &(centerpoint, _) in &context.neighbors[seed] {
            if !context.allowed_edge(seed, centerpoint) {
                continue;
            }
            let seed_bond = context.bond_type(state, seed, centerpoint);
            if !spec.seed_bonds.contains(&seed_bond) {
                continue;
            }
            let cp_record = &context.structure.atoms[centerpoint];
            if !spec
                .centerpoint_elements
                .contains(&cp_record.element.as_str())
            {
                continue;
            }
            if let Some(want) = spec.centerpoint_valence {
                if context.neighbors[centerpoint].len() != want {
                    continue;
                }
            }
            if spec.centerpoint_neutral && cp_record.charge != 0 {
                continue;
            }
            if spec.ket_normal_carbon
                && (cp_record.charge != 0
                    || context.bond_sum2[centerpoint] / 2 + context.hydrogens[centerpoint] as i32
                        != 4)
            {
                continue;
            }
            let mut ends: Vec<usize> = Vec::new();
            let mut donors = 0u32;
            let mut acceptors = 0u32;
            let mut possible_mobile = 0u32;
            let mut admitted = Counts::default();
            let mut all_neighbors = Counts::default();
            if spec.gate.pt39_element_gate {
                // PT_39_00 counts its O/N population over every
                // neighbor at scan entry, admitted or not
                // (ichitaut.c 5107–5108).
                for &(neighbor, _) in &context.neighbors[centerpoint] {
                    match context.structure.atoms[neighbor].element.as_str() {
                        "O" => all_neighbors.oxygen += 1,
                        "N" => all_neighbors.nitrogen += 1,
                        _ => {}
                    }
                }
            }
            for &(endpoint, _) in &context.neighbors[centerpoint] {
                if !context.allowed_edge(centerpoint, endpoint) {
                    continue;
                }
                let bond = context.bond_type(state, centerpoint, endpoint);
                let non_taut = match bond {
                    BondType::Single | BondType::Double => true,
                    BondType::Triple if wide_bonds => true,
                    BondType::Triple | BondType::Alt23 => continue,
                    _ => false,
                };
                let Some(info) = (match spec.predicate {
                    Predicate::Base => context.endpoint_info(state, endpoint),
                    _ => context.endpoint_info_ext(state, endpoint, spec.predicate),
                }) else {
                    continue;
                };
                if spec.exclude_seed_element
                    && endpoint != seed
                    && context.structure.atoms[endpoint].element
                        == context.structure.atoms[seed].element
                {
                    continue;
                }
                if spec.pt18_double_skip
                    && bond == BondType::Double
                    && info.valence == 3
                    && info.mobile == 0
                {
                    continue;
                }
                let registered = state.registered(endpoint);
                let (donor_side, acceptor_side) = match spec.scan {
                    ScanStyle::Standard { seed_verdict_alt } => {
                        if non_taut {
                            (
                                bond == BondType::Single && (info.mobile > 0 || registered),
                                bond == BondType::Double,
                            )
                        } else {
                            let verdict = if seed_verdict_alt { seed_info } else { info };
                            (
                                registered || verdict.role == Role::Donor,
                                registered
                                    || verdict.neutral_bonds_valence
                                        > context.neighbors[endpoint].len() as i32,
                            )
                        }
                    }
                    ScanStyle::Mobile {
                        alt_acceptor_by_neutral,
                    } => {
                        if non_taut {
                            (info.mobile > 0 || registered, info.mobile == 0)
                        } else {
                            (
                                info.role == Role::Donor,
                                alt_acceptor_by_neutral
                                    && info.neutral_bonds_valence
                                        > context.neighbors[endpoint].len() as i32,
                            )
                        }
                    }
                };
                if !donor_side && !acceptor_side {
                    continue;
                }
                donors += donor_side as u32;
                acceptors += acceptor_side as u32;
                possible_mobile += (info.mobile > 0 || registered) as u32;
                admitted.valence2 += (info.valence == 2) as u32;
                admitted.valence3 += (info.valence == 3) as u32;
                admitted.valence4 += (info.valence == 4) as u32;
                ends.push(endpoint);
            }
            if ends.len() <= 1 || possible_mobile == 0 || donors == 0 || acceptors == 0 {
                continue;
            }
            if !spec.gate.allows(&admitted, &all_neighbors) {
                continue;
            }
            let partitioned = find_accessible_end_points(context, state, &ends, spec.predicate);
            changed |= register_partition(context, state, centerpoint, &partitioned);
        }
    }
    changed
}

/// The flag passes in the reference's MarkTautomerGroups order:
/// PT_22_00, PT_16_00, PT_06_00, PT_39_00, PT_13_00, PT_18_00, KET.
/// Every enabled flag's block runs (the C does not short-circuit);
/// the fixpoint loop re-runs the group until nothing changes.
fn flag_round(context: &Context, state: &mut Registration) -> bool {
    let flags = context.structure.taut_ext;
    static DOUBLE_LIKE: [BondType; 4] = [
        BondType::Double,
        BondType::Altern,
        BondType::Alt12ns,
        BondType::Tautom,
    ];
    static WIDE: [BondType; 5] = [
        BondType::Double,
        BondType::Single,
        BondType::Altern,
        BondType::Alt12ns,
        BondType::Tautom,
    ];
    let mut changed = false;
    if flags.pt_22_00 {
        changed |= anchor_round(
            context,
            state,
            &AnchorSpec {
                predicate: Predicate::Pt22,
                seed_bonds: &DOUBLE_LIKE,
                centerpoint_elements: &["N"],
                centerpoint_valence: None,
                centerpoint_neutral: false,
                ket_normal_carbon: false,
                exclude_seed_element: false,
                pt18_double_skip: false,
                scan: ScanStyle::Standard {
                    seed_verdict_alt: false,
                },
                gate: Gate::default(),
            },
        );
    }
    if flags.pt_16_00 {
        changed |= anchor_round(
            context,
            state,
            &AnchorSpec {
                predicate: Predicate::Pt16,
                seed_bonds: &DOUBLE_LIKE,
                centerpoint_elements: &["N"],
                centerpoint_valence: Some(2),
                centerpoint_neutral: true,
                ket_normal_carbon: false,
                exclude_seed_element: true,
                pt18_double_skip: false,
                scan: ScanStyle::Standard {
                    seed_verdict_alt: false,
                },
                gate: Gate::default(),
            },
        );
    }
    if flags.pt_06_00 {
        changed |= anchor_round(
            context,
            state,
            &AnchorSpec {
                predicate: Predicate::Pt06,
                seed_bonds: &DOUBLE_LIKE,
                centerpoint_elements: &["N", "C", "P"],
                centerpoint_valence: None,
                centerpoint_neutral: false,
                ket_normal_carbon: false,
                exclude_seed_element: false,
                pt18_double_skip: false,
                scan: ScanStyle::Standard {
                    seed_verdict_alt: false,
                },
                gate: Gate::default(),
            },
        );
    }
    if flags.pt_39_00 {
        changed |= anchor_round(
            context,
            state,
            &AnchorSpec {
                predicate: Predicate::Pt39,
                seed_bonds: &DOUBLE_LIKE,
                centerpoint_elements: &["N"],
                centerpoint_valence: Some(3),
                centerpoint_neutral: false,
                ket_normal_carbon: false,
                exclude_seed_element: false,
                pt18_double_skip: false,
                scan: ScanStyle::Standard {
                    seed_verdict_alt: false,
                },
                gate: Gate {
                    pt39_element_gate: true,
                    ..Gate::default()
                },
            },
        );
    }
    if flags.pt_13_00 {
        changed |= anchor_round(
            context,
            state,
            &AnchorSpec {
                predicate: Predicate::Pt13,
                seed_bonds: &WIDE,
                centerpoint_elements: &["C"],
                centerpoint_valence: Some(2),
                centerpoint_neutral: false,
                ket_normal_carbon: false,
                exclude_seed_element: false,
                pt18_double_skip: false,
                scan: ScanStyle::Mobile {
                    alt_acceptor_by_neutral: true,
                },
                gate: Gate {
                    valence2_equals: Some(1),
                    valence4: Some(Valence4Gate::ExactlyOne),
                    ..Gate::default()
                },
            },
        );
    }
    if flags.pt_18_00 {
        changed |= anchor_round(
            context,
            state,
            &AnchorSpec {
                predicate: Predicate::Pt18,
                seed_bonds: &WIDE,
                centerpoint_elements: &["C"],
                centerpoint_valence: Some(2),
                centerpoint_neutral: false,
                ket_normal_carbon: false,
                exclude_seed_element: false,
                pt18_double_skip: true,
                scan: ScanStyle::Mobile {
                    alt_acceptor_by_neutral: false,
                },
                gate: Gate {
                    valence2_equals: Some(1),
                    valence3_equals: Some(1),
                    ..Gate::default()
                },
            },
        );
    }
    if flags.keto_enol {
        changed |= anchor_round(
            context,
            state,
            &AnchorSpec {
                predicate: Predicate::Ket,
                seed_bonds: &DOUBLE_LIKE,
                centerpoint_elements: &["C"],
                centerpoint_valence: None,
                centerpoint_neutral: false,
                ket_normal_carbon: true,
                exclude_seed_element: false,
                pt18_double_skip: false,
                scan: ScanStyle::Standard {
                    seed_verdict_alt: true,
                },
                gate: Gate {
                    valence2_equals: Some(1),
                    valence4: Some(Valence4Gate::AtLeastOne),
                    ..Gate::default()
                },
            },
        );
    }
    changed
}

/// FindAccessibleEndPoints: partition the admitted endpoints by the
/// balanced network's cap-aware connectivity. Endpoints already in
/// groups merge when their representatives confirm an alt path; new
/// endpoints pair among themselves and with group members; classes
/// of one endpoint are dropped (IGNORE_SINGLE_ENDPOINTS). `pred` is
/// the calling round's mode — the confirmations' terminal prechecks
/// consult that round's predicate.
fn find_accessible_end_points(
    context: &Context,
    state: &Registration,
    ends: &[usize],
    pred: Predicate,
) -> Partitioned {
    let mut equ: Vec<usize> = ends
        .iter()
        .map(|&e| state.group_of[e].unwrap_or(0))
        .collect();
    // Representatives of the known groups, first encounter first,
    // with the slot label each carries.
    let mut reps: Vec<(usize, usize)> = Vec::new();
    let mut slots: Vec<usize> = Vec::new();
    for (index, &atom) in ends.iter().enumerate() {
        if equ[index] != 0 && slots.iter().take(reps.len()).all(|&s| s != equ[index]) {
            reps.push((atom, index));
            slots.push(equ[index]);
        }
    }
    let known = reps.len();
    let mut merged = false;
    for i in 0..known {
        for j in (i + 1)..known {
            if context.alt_path_exists(state, reps[i].0, reps[j].0, pred) {
                merged = true;
                let low = slots[i].min(slots[j]);
                let high = slots[i].max(slots[j]);
                for slot in slots.iter_mut() {
                    if *slot == high {
                        *slot = low;
                    }
                }
                for label in equ.iter_mut() {
                    if *label == high {
                        *label = low;
                    }
                }
            }
        }
    }
    // After a merge, only one representative per surviving group
    // label stays.
    let mut members: Vec<(usize, usize)> = Vec::new();
    if merged {
        let mut seen: Vec<usize> = Vec::new();
        for index in 0..reps.len() {
            if seen.contains(&slots[index]) {
                continue;
            }
            seen.push(slots[index]);
            members.push(reps[index]);
        }
        slots = seen;
    } else {
        members = reps.clone();
    }
    for (index, &atom) in ends.iter().enumerate() {
        if equ[index] == 0 {
            members.push((atom, index));
            slots.push(0);
        }
    }
    let mut next_fictive = context.structure.atoms.len() + 1;
    for i in 0..members.len() {
        for j in (i + 1)..members.len() {
            let tested = (slots[i] != slots[j] && (i >= known || j >= known))
                || (slots[i] == 0 && slots[j] == 0);
            if !tested || !context.alt_path_exists(state, members[i].0, members[j].0, pred) {
                continue;
            }
            let (a, b) = (members[i].1, members[j].1);
            if slots[i] != 0 && slots[j] != 0 {
                let low = slots[i].min(slots[j]);
                let high = slots[i].max(slots[j]);
                for slot in slots.iter_mut() {
                    if *slot == high {
                        *slot = low;
                    }
                }
                for label in equ.iter_mut() {
                    if *label == high {
                        *label = low;
                    }
                }
            } else if slots[i] != 0 {
                slots[j] = slots[i];
                equ[b] = slots[i];
            } else if slots[j] != 0 {
                slots[i] = slots[j];
                equ[a] = slots[j];
            } else {
                slots[i] = next_fictive;
                slots[j] = next_fictive;
                equ[a] = next_fictive;
                equ[b] = next_fictive;
                next_fictive += 1;
            }
        }
    }
    let mut partitioned = Partitioned::default();
    for (index, &atom) in ends.iter().enumerate() {
        if equ[index] == 0
            || ends
                .iter()
                .zip(&equ)
                .filter(|(_, label)| *label == &equ[index])
                .count()
                <= 1
        {
            continue;
        }
        partitioned.ends.push(atom);
        partitioned.equ.push(equ[index]);
    }
    partitioned
}

/// The FindAccessibleEndPoints survivors: the endpoints and the
/// equivalence label each carries into RegisterEndPoints (the parallel
/// bond list is derivable — every survivor hangs off the anchor).
#[derive(Default)]
struct Partitioned {
    ends: Vec<usize>,
    equ: Vec<usize>,
}

/// A template check: receives the live registration state, the path
/// atoms, and the path bonds (bond i joins atom i to atom i + 1; the
/// last entry of a ring path is the closure bond).
type DfsCheck<'ck> = &'ck mut dyn FnMut(&mut Registration, &[usize], &[BondType]) -> bool;

/// The shift-walk templates, run only on a pass whose 1,3 phase
/// registered nothing (the reference's `!tot_changes` gate: one
/// MarkTautomerGroups call does the anchors, and only a no-change
/// call falls through to the TAUT_OTHER block). Four template
/// families, in the reference's order: the 1,5 alt path
/// (TAUT_15_NON_RING), 6-member alt rings (TAUT_4PYRIDINOL_RINGS),
/// pyrazole 5-rings (TAUT_PYRAZOLE_RINGS), and the 1,4 tropolone
/// pattern (TAUT_TROPOLONE_5/7).
fn walk_round(context: &Context, state: &mut Registration) -> bool {
    let n = context.structure.atoms.len();
    let mut changed = false;

    // TAUT_15_NON_RING is runtime-gated in the reference
    // (TG_FLAG_1_5_TAUT, off for the standard InChI flow), so the
    // non-ring 1,5 walk never runs here.

    // TAUT_4PYRIDINOL_RINGS: a ring N (endpoint valence 3, two ring
    // bonds) walks 6-member cycles to the para ring atom bearing an
    // exocyclic =O / -OH / =NH / -NH2; BOTH arcs to the N alternate.
    for start in 0..n {
        if context.rings.size[start] < 6
            || context.neighbors[start].len() != 2
            || !context.is_ring_nitrogen(state, start)
        {
            continue;
        }
        let mut path = vec![start];
        let mut bonds = Vec::new();
        let mut on_path = vec![false; n];
        on_path[start] = true;
        changed |= walk_ring_path(
            context,
            state,
            &mut path,
            &mut bonds,
            &mut on_path,
            &[],
            5,
            &mut |state, path, bonds| check_6memb_ring(context, state, path, bonds),
        );
    }

    // TAUT_PYRAZOLE_RINGS: two adjacent ring N's (each two ring
    // bonds, one ring system); a fresh pair needs exactly one mobile
    // H between them; the ring walk closes 5-cycles through both.
    for i1 in 0..n {
        if context.rings.size[i1] < 5
            || context.neighbors[i1].len() != 2
            || !context.is_ring_nitrogen(state, i1)
        {
            continue;
        }
        for &(i2, _) in &context.neighbors[i1] {
            if i2 >= i1 {
                continue;
            }
            if context.rings.system[i2] != context.rings.system[i1] {
                continue;
            }
            if !matches!(
                context.bond_type(state, i1, i2),
                BondType::Single | BondType::Tautom | BondType::Altern | BondType::Alt12ns
            ) {
                continue;
            }
            if context.neighbors[i2].len() != 2 || !context.is_ring_nitrogen(state, i2) {
                continue;
            }
            if !state.registered(i1)
                && !state.registered(i2)
                && context.mobile(i1) + context.mobile(i2) != 1
            {
                continue;
            }
            let mut path = vec![i1, i2];
            let mut bonds = vec![context.bond_type(state, i1, i2)];
            let mut on_path = vec![false; n];
            on_path[i1] = true;
            on_path[i2] = true;
            changed |= walk_ring_path(
                context,
                state,
                &mut path,
                &mut bonds,
                &mut on_path,
                &[],
                4,
                &mut |state, path, bonds| check_pyrazole(context, state, path, bonds),
            );
        }
    }

    // TAUT_TROPOLONE_5/7: two ring cut-vertices Q-R (each three
    // neighbors, strict centerpoints, unregistered) joined by a
    // single-like bond, each carrying an EXOCYCLIC endpoint; the H
    // shifts between them across the ring path (the Q-R bond
    // excluded from the alternation).
    for i1 in 0..n {
        if context.rings.size[i1] < 5
            || !context.is_centerpoint_strict(i1)
            || !context.rings.cut_vertex[i1]
            || context.neighbors[i1].len() != 3
            || state.registered(i1)
        {
            continue;
        }
        for &(i2, _) in &context.neighbors[i1] {
            if context.rings.system[i2] != context.rings.system[i1]
                || !context.is_centerpoint_strict(i2)
                || !context.rings.cut_vertex[i2]
                || context.neighbors[i2].len() != 3
                || state.registered(i2)
            {
                continue;
            }
            if !matches!(
                context.bond_type(state, i1, i2),
                BondType::Single | BondType::Tautom | BondType::Altern | BondType::Alt12ns
            ) {
                continue;
            }
            for &(endpoint1, _) in &context.neighbors[i1] {
                if endpoint1 == i2 {
                    continue;
                }
                let Some(eif1) = context.endpoint(state, endpoint1) else {
                    continue;
                };
                if context.rings.system[endpoint1] == context.rings.system[i1] {
                    continue;
                }
                let b1 = context.bond_type(state, i1, endpoint1);
                if !matches!(
                    b1,
                    BondType::Single
                        | BondType::Double
                        | BondType::Tautom
                        | BondType::Altern
                        | BondType::Alt12ns
                ) {
                    continue;
                }
                for &(endpoint2, _) in &context.neighbors[i2] {
                    if endpoint2 == i1 {
                        continue;
                    }
                    let Some(eif2) = context.endpoint(state, endpoint2) else {
                        continue;
                    };
                    if context.rings.system[endpoint2] == context.rings.system[i2] {
                        continue;
                    }
                    let b2 = context.bond_type(state, i2, endpoint2);
                    if !matches!(
                        b2,
                        BondType::Single
                            | BondType::Double
                            | BondType::Tautom
                            | BondType::Altern
                            | BondType::Alt12ns
                    ) {
                        continue;
                    }
                    let mut shiftable = false;
                    shiftable |= (state.registered(endpoint1) || context.mobile(endpoint1) > 0)
                        && b1 != BondType::Double
                        && (state.registered(endpoint2) || eif2 == Role::Acceptor)
                        && b2 != BondType::Single;
                    shiftable |= (state.registered(endpoint2) || context.mobile(endpoint2) > 0)
                        && b2 != BondType::Double
                        && (state.registered(endpoint1) || eif1 == Role::Acceptor)
                        && b1 != BondType::Single;
                    if !shiftable {
                        continue;
                    }
                    if b1 == b2 && matches!(b1, BondType::Single | BondType::Double) {
                        continue;
                    }
                    for cycle_bonds in [4, 6] {
                        if cycle_bonds == 6 && context.rings.size[i1] < 7 {
                            continue;
                        }
                        let mut path = vec![i1, i2];
                        let mut bonds = vec![context.bond_type(state, i1, i2)];
                        let mut on_path = vec![false; n];
                        on_path[i1] = true;
                        on_path[i2] = true;
                        let mut find = TropoloneFind::default();
                        walk_ring_path(
                            context,
                            state,
                            &mut path,
                            &mut bonds,
                            &mut on_path,
                            &[endpoint1, endpoint2],
                            cycle_bonds,
                            &mut |state, path, bonds| {
                                check_tropolone(
                                    context, state, path, bonds, endpoint1, endpoint2, &mut find,
                                )
                            },
                        );
                        if find.found > 0 {
                            changed |= state.register(&find.endpoints);
                            for &bond in &find.bonds {
                                changed |= state.marks.insert(bond);
                            }
                        }
                    }
                }
            }
        }
    }

    changed
}

/// The reference bond classes the templates distinguish. TAUTOM is a
/// bond a previous registration replaced; ALTERN is an aromatic bond.
/// ALT12NS ("1 or 2", the ride's mark 9) and ALT23 ("2 or 3", mark
/// 7) are the mobility marks a hard proton ride leaves: the first
/// matches either salt role, the second neither.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum BondType {
    Single,
    Double,
    Triple,
    Altern,
    Tautom,
    Alt12ns,
    Alt23,
}

impl<'a> Context<'a> {
    /// The reference's `3 == nGetEndpointInfo(...)`: a classified
    /// nitrogen endpoint.
    fn is_ring_nitrogen(&self, state: &Registration, atom: usize) -> bool {
        self.structure.atoms[atom].element == "N" && self.endpoint(state, atom).is_some()
    }
}

/// `are_alt_bonds`: the path bonds must alternate single/double —
/// ALTERN and TAUTOM match either — and the return is the bond type
/// the final position holds, which the callers read as the pattern
/// bond at the far endpoint. `None` rejects. (The C returns
/// non-single/non-double in the wild-card case, which the callers
/// treat as "no pattern constraint".)
fn are_alt_bonds(bonds: &[BondType]) -> Option<BondType> {
    if bonds.len() < 2 || matches!(bonds[0], BondType::Triple | BondType::Alt23) {
        return None;
    }
    let mut expected = match bonds[0] {
        BondType::Single => Some(BondType::Double),
        BondType::Double => Some(BondType::Single),
        _ => None,
    };
    for &bond in &bonds[1..] {
        let wildcard = matches!(
            bond,
            BondType::Altern | BondType::Tautom | BondType::Alt12ns
        );
        if let Some(want) = expected {
            if bond == want || wildcard {
                expected = Some(if want == BondType::Single {
                    BondType::Double
                } else {
                    BondType::Single
                });
            } else {
                return None;
            }
        } else {
            match bond {
                BondType::Single => expected = Some(BondType::Double),
                BondType::Double => expected = Some(BondType::Single),
                _ if wildcard => {}
                _ => return None,
            }
        }
    }
    Some(match expected {
        Some(BondType::Single) => BondType::Double,
        Some(_) => BondType::Single,
        None => BondType::Altern,
    })
}

/// Mark every marked-able bond of a found path: REPLACE_THE_BOND
/// replaces single, double and alternating bonds with TAUTOM (a
/// triple never; an already-tautom bond stays).
fn mark_path_bonds(
    context: &Context,
    state: &mut Registration,
    path: &[usize],
    bonds: &[BondType],
) -> bool {
    let ring = bonds.len() == path.len();
    let mut marked = false;
    for (step, bond) in path.windows(2).zip(bonds) {
        if matches!(
            *bond,
            BondType::Single | BondType::Double | BondType::Altern | BondType::Alt12ns
        ) {
            if let Some(&index) = context
                .bond_index
                .get(&(step[0].min(step[1]), step[0].max(step[1])))
            {
                marked |= state.marks.insert(index);
            }
        }
    }
    if ring {
        let last = path[path.len() - 1];
        let first = path[0];
        if matches!(
            bonds[bonds.len() - 1],
            BondType::Single | BondType::Double | BondType::Altern | BondType::Alt12ns
        ) {
            if let Some(&index) = context.bond_index.get(&(first.min(last), first.max(last))) {
                marked |= state.marks.insert(index);
            }
        }
    }
    marked
}

/// The reference's DFS_FindTautInARing: paths around a cycle of
/// `cycle_bonds` bonds closing back into the seed's first atom;
/// interior atoms must be strict centerpoints; forbidden atoms (the
/// 1,4 template's exocyclic endpoints) are never touched.
#[allow(clippy::too_many_arguments)]
fn walk_ring_path(
    context: &Context,
    state: &mut Registration,
    path: &mut Vec<usize>,
    bonds: &mut Vec<BondType>,
    on_path: &mut [bool],
    forbidden: &[usize],
    cycle_bonds: usize,
    check: DfsCheck<'_>,
) -> bool {
    let level = path.len() - 1;
    let atom = path[level];
    let mut changed = false;
    for &(next, _) in &context.neighbors[atom] {
        if forbidden.contains(&next) {
            continue;
        }
        let bond = context.bond_type(state, atom, next);
        if on_path[next] {
            if next == path[0] && level == cycle_bonds {
                bonds.push(bond);
                changed |= check(state, path, bonds);
                bonds.pop();
            }
            continue;
        }
        if !context.is_centerpoint_strict(next) {
            continue;
        }
        if level < cycle_bonds {
            on_path[next] = true;
            path.push(next);
            bonds.push(bond);
            changed |= walk_ring_path(
                context,
                state,
                path,
                bonds,
                on_path,
                forbidden,
                cycle_bonds,
                check,
            );
            bonds.pop();
            path.pop();
            on_path[next] = false;
        }
    }
    changed
}

/// Check5MembTautRing (pyrazole): the fresh-pair one-H gate, an
/// alternating 5-cycle, and the closure bond's pattern deciding
/// which N holds the hydrogen.
fn check_pyrazole(
    context: &Context,
    state: &mut Registration,
    path: &[usize],
    bonds: &[BondType],
) -> bool {
    let n1 = path[0];
    let n2 = path[1];
    if context.endpoint(state, n1).is_none() || context.endpoint(state, n2).is_none() {
        return false;
    }
    if !state.registered(n1)
        && !state.registered(n2)
        && context.mobile(n1) + context.mobile(n2) != 1
    {
        return false;
    }
    // The balanced network confirms the fresh pair (the reference
    // skips it only when both N's already share one group).
    if !(state.group_of[n1].is_some() && state.group_of[n1] == state.group_of[n2])
        && !context.alt_path_exists(state, n1, n2, Predicate::Base)
    {
        return false;
    }
    // The reference's path_bonds run from n2 around to the closure
    // INTO n1 (DfsPath[1..]) — the seeded n1-n2 bond is NOT part of
    // the alternation; its pattern type follows from the arc (the two
    // ring arcs alternate). Feeding the seed bond in would demand the
    // arc continue the seed's phase and reject fused cases (the
    // indazole single-single junction) the reference accepts.
    let Some(last) = are_alt_bonds(&bonds[1..]) else {
        return false;
    };
    match last {
        BondType::Single => {
            if (!state.registered(n1) && context.endpoint(state, n1) != Some(Role::Donor))
                || (!state.registered(n2) && context.endpoint(state, n2) != Some(Role::Acceptor))
            {
                return false;
            }
        }
        BondType::Double
            if ((!state.registered(n1) && context.endpoint(state, n1) != Some(Role::Acceptor))
                || (!state.registered(n2) && context.endpoint(state, n2) != Some(Role::Donor))) =>
        {
            return false;
        }
        _ => {}
    }
    let registered = state.register(&[n1, n2]);
    let marked = mark_path_bonds(context, state, path, bonds);
    registered || marked
}

/// Check6MembTautRing (4-pyridinol): the para ring atom (three
/// neighbors) carries an exocyclic endpoint; both arcs to the N
/// alternate; a mobile group must actually shift along one.
fn check_6memb_ring(
    context: &Context,
    state: &mut Registration,
    path: &[usize],
    bonds: &[BondType],
) -> bool {
    let n_at = path[0];
    let mid = path[3];
    if context.neighbors[mid].len() != 3 || !context.rings.cut_vertex[mid] {
        return false;
    }
    let mut o_at = None;
    for &(neighbor, _) in &context.neighbors[mid] {
        if neighbor != path[2] && neighbor != path[4] {
            o_at = Some(neighbor);
            break;
        }
    }
    let Some(o_at) = o_at else {
        return false;
    };
    let exo = context.bond_type(state, mid, o_at);
    if !matches!(
        exo,
        BondType::Single
            | BondType::Double
            | BondType::Tautom
            | BondType::Altern
            | BondType::Alt12ns
    ) {
        return false;
    }
    let Some(eif_o) = context.endpoint(state, o_at) else {
        return false;
    };
    if exo == BondType::Single && !state.registered(o_at) && eif_o != Role::Donor {
        return false;
    }
    let Some(eif_n) = context.endpoint(state, n_at) else {
        return false;
    };
    let mut shiftable = false;
    shiftable |= (state.registered(o_at) || eif_o == Role::Donor)
        && exo != BondType::Double
        && (state.registered(n_at) || eif_n == Role::Acceptor);
    shiftable |= (state.registered(n_at) || eif_n == Role::Donor)
        && (state.registered(o_at) || eif_o == Role::Acceptor)
        && exo != BondType::Single;
    if !shiftable {
        return false;
    }
    let closure = bonds[bonds.len() - 1];
    if !are_alt_bonds(&[exo, bonds[2], bonds[1], bonds[0]]).is_some()
        || !are_alt_bonds(&[exo, bonds[3], bonds[4], closure]).is_some()
    {
        return false;
    }
    // The balanced network confirms the pair (the reference skips it
    // only when both already share one group).
    if !(state.group_of[o_at].is_some() && state.group_of[o_at] == state.group_of[n_at])
        && !context.alt_path_exists(state, n_at, o_at, Predicate::Base)
    {
        return false;
    }
    let registered = state.register(&[n_at, o_at]);
    let marked = mark_path_bonds(context, state, path, bonds);
    registered || marked
}

/// Check7MembTautRing (the 1,4 tropolone pattern, its 5-member
/// form): the exocyclic endpoints hang off the seeded ring path; the
/// alternation runs from one exocyclic bond across the Q-R bond and
/// the ring (the closure bond excluded) to the other. Endpoints and
/// bonds accumulate on the shared [`TropoloneFind`] of one DFS pass —
/// a pair the balanced network rejects still leaves its collection in
/// place, and a later accepted pair registers the whole accumulation
/// (the reference's wrapper arrays).
#[allow(clippy::too_many_arguments)]
fn check_tropolone(
    context: &Context,
    state: &mut Registration,
    path: &[usize],
    bonds: &[BondType],
    endpoint1: usize,
    endpoint2: usize,
    find: &mut TropoloneFind,
) -> bool {
    let i1 = path[0];
    let i2 = path[1];
    let Some(eif1) = context.endpoint(state, endpoint1) else {
        return false;
    };
    let Some(eif2) = context.endpoint(state, endpoint2) else {
        return false;
    };
    // Chain: endpoint2's bond, the ring arc (the walk bonds from i2
    // around to the closure INTO i1 — the seed Q-R bond EXCLUDED, the
    // reference's DfsPath[1..=nLenDfsPath]), endpoint1's bond. The
    // caller already vetted the Q-R bond as single/tautom/altern.
    let mut chain: Vec<BondType> = Vec::with_capacity(bonds.len() + 1);
    chain.push(context.bond_type(state, i2, endpoint2));
    chain.extend_from_slice(&bonds[1..]);
    chain.push(context.bond_type(state, i1, endpoint1));
    let Some(last) = are_alt_bonds(&chain) else {
        return false;
    };
    match last {
        BondType::Single => {
            if (!state.registered(endpoint1) && eif1 != Role::Donor)
                || (!state.registered(endpoint2) && eif2 != Role::Acceptor)
            {
                return false;
            }
        }
        BondType::Double
            if (!state.registered(endpoint1) && eif1 != Role::Acceptor)
                || (!state.registered(endpoint2) && eif2 != Role::Donor) =>
        {
            return false;
        }
        _ => {}
    }
    // AddEndPoints / AddBondsPos: the shared collection dedups;
    // `added` decides whether this pass counts as found.
    let mut added = false;
    for &endpoint in &[endpoint2, endpoint1] {
        if !find.endpoints.contains(&endpoint) {
            find.endpoints.push(endpoint);
            added = true;
        }
    }
    let mut chain_bonds: Vec<(usize, usize)> = vec![(i2, endpoint2)];
    for index in 1..bonds.len() - 1 {
        chain_bonds.push((path[index], path[index + 1]));
    }
    chain_bonds.push((path[bonds.len() - 1], path[0]));
    chain_bonds.push((i1, endpoint1));
    for ((center, other), kind) in chain_bonds.into_iter().zip(&chain) {
        if !matches!(kind, BondType::Single | BondType::Double | BondType::Altern) {
            continue;
        }
        let index = context.bond_index[&(center.min(other), center.max(other))];
        if !find.bonds.contains(&index) {
            find.bonds.push(index);
            added = true;
        }
    }
    if !added {
        return false;
    }
    // A pair already in one group needs no confirmation.
    if state.group_of[endpoint2].is_some() && state.group_of[endpoint2] == state.group_of[endpoint1]
    {
        find.found += 1;
        return true;
    }
    if context.alt_path_exists(state, endpoint2, endpoint1, Predicate::Base) {
        find.found += 1;
        true
    } else {
        // The network said no: this pass contributes nothing, but the
        // collected entries stay for a later accepted pair.
        false
    }
}

/// One DFS pass's shared collection for the tropolone template.
#[derive(Default)]
struct TropoloneFind {
    endpoints: Vec<usize>,
    bonds: Vec<usize>,
    found: i32,
}

/// Groups of endpoints connected by 1,3 anchors or shift walks,
/// filtered to those with at least two members and at least one
/// mobile H. Rank-free — membership depends only on the structure.
fn components(
    context: &Context,
    state: &Registration,
    hydrogens: &[u32],
    charge_units: &impl Fn(usize) -> u32,
) -> Vec<Vec<usize>> {
    let n = context.structure.atoms.len();
    let linked = build_links(context);
    let mut component = vec![usize::MAX; n];
    let mut groups = Vec::new();
    for start in 0..n {
        if context.endpoint(state, start).is_none() || component[start] != usize::MAX {
            continue;
        }
        let mut members = Vec::new();
        let mut queue = vec![start];
        component[start] = start;
        while let Some(atom) = queue.pop() {
            members.push(atom);
            for &next in &linked[atom] {
                if component[next] == usize::MAX {
                    component[next] = start;
                    queue.push(next);
                }
            }
        }
        let count: u32 = members
            .iter()
            .map(|&atom| hydrogens[atom] + charge_units(atom))
            .sum();
        if members.len() >= 2 && count > 0 {
            groups.push(members);
        }
    }
    groups
}

pub fn mobile_groups(structure: &Structure, hydrogens: &[u32], ranks: &[u32]) -> Vec<MobileGroup> {
    let context = Context::new(structure, hydrogens);
    let state = fresh_state(structure.atoms.len());
    let charge_units = charge_units_of(structure);
    let mut groups = Vec::new();
    for members in components(&context, &state, hydrogens, &charge_units) {
        let count: u32 = members
            .iter()
            .map(|&atom| hydrogens[atom] + charge_units(atom))
            .sum();
        let charge: u32 = members.iter().map(|&atom| charge_units(atom)).sum();
        let mut group_ranks: Vec<u32> = members.iter().map(|&atom| ranks[atom]).collect();
        group_ranks.sort_unstable();
        groups.push(MobileGroup {
            count,
            charge,
            ranks: group_ranks,
        });
    }
    groups.sort_by_key(|group| (group.ranks.len(), group.ranks[0]));
    groups
}

/// The atoms whose hydrogens live in mobile groups (excluded from the
/// fixed-H entries).
pub fn endpoint_atoms(ranks: &[u32], groups: &[MobileGroup]) -> BTreeSet<usize> {
    let mut by_rank = std::collections::BTreeMap::new();
    for (atom, &rank) in ranks.iter().enumerate() {
        by_rank.insert(rank, atom);
    }
    let mut endpoints = BTreeSet::new();
    for group in groups {
        for &rank in &group.ranks {
            if let Some(&atom) = by_rank.get(&rank) {
                endpoints.insert(atom);
            }
        }
    }
    endpoints
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
    fn carboxyl_pairs_the_two_oxygens() {
        let s = structure(
            vec![atom("C"), atom("C"), atom("O"), atom("O")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 1], &[1, 2, 3, 4]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].count, 1);
        assert_eq!(groups[0].ranks, vec![3, 4]);
    }

    #[test]
    fn alcohol_oxygen_is_not_an_endpoint() {
        let s = structure(
            vec![atom("C"), atom("C"), atom("O")],
            vec![(0, 1, 1), (1, 2, 1)],
        );
        assert!(endpoint_mask(&s, &[3, 2, 1]).iter().all(|&e| !e));
        assert!(mobile_groups(&s, &[3, 2, 1], &[1, 2, 3]).is_empty());
    }

    #[test]
    fn amide_nitrogen_joins_the_carbonyl_oxygen() {
        let s = structure(
            vec![atom("C"), atom("C"), atom("O"), atom("N")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 2], &[1, 2, 3, 4]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].count, 2);
        assert_eq!(groups[0].ranks, vec![3, 4]);
        let endpoints = endpoint_atoms(&[1, 2, 3, 4], &groups);
        assert_eq!(endpoints, BTreeSet::from([2, 3]));
    }

    #[test]
    fn the_nitro_anion_groups_its_oxygens_but_not_its_nitrogen() {
        // CH3-N(=O)-O(-): the oracle records (H,3,4) over the two
        // oxygens; the nitrogen's bond sum (4) overflows the endpoint
        // valence of N (3), so its hydrogens stay fixed.
        let s = structure(
            vec![atom("C"), atom("N"), atom("O"), atom("O")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 1)],
        );
        let groups = mobile_groups(&s, &[3, 0, 0, 1], &[1, 2, 3, 4]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].ranks, vec![3, 4]);
        assert!(!endpoint_mask(&s, &[3, 0, 0, 1])[1]);
    }

    #[test]
    fn sulfone_sulfur_never_participates_but_thiocarbonyl_does() {
        // R-S(=O)2-R': bond sum 6 overflows S's endpoint valence 2.
        let s = structure(
            vec![atom("C"), atom("S"), atom("O"), atom("O"), atom("C")],
            vec![(0, 1, 1), (1, 2, 2), (1, 3, 2), (1, 4, 1)],
        );
        assert!(!endpoint_mask(&s, &[3, 0, 0, 0, 3])[1]);
        // R-C(=S)-NH2: the thiocarbonyl S is an acceptor endpoint.
        let s = structure(
            vec![atom("C"), atom("S"), atom("N")],
            vec![(0, 1, 2), (0, 2, 1)],
        );
        assert!(endpoint_mask(&s, &[0, 0, 2])[1]);
    }

    #[test]
    fn kekule_double_bond_anchors_the_pyridinol_pair() {
        // 2-hydroxypyridine skeleton in the kekule form the reference
        // normalizes: the OH oxygen (donor) sits on the anchor carbon
        // whose =N- bond (acceptor) is a real input double, so the
        // 1,3 template pairs them.
        let mut atoms: Vec<AtomRecord> = (0..6).map(|_| atom("C")).collect();
        atoms[0] = atom("N");
        atoms.push(atom("O"));
        let s = structure(
            atoms,
            vec![
                (0, 1, 2), // ring N = C1 (the anchor's pi bond)
                (1, 2, 1),
                (2, 3, 2),
                (3, 4, 1),
                (4, 5, 2),
                (5, 0, 1),
                (1, 6, 1), // exocyclic OH on C1
            ],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 0, 0, 0, 1], &[1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].count, 1);
        assert_eq!(groups[0].ranks, vec![1, 7]);
    }

    #[test]
    fn two_phenolic_oxygens_never_merge_across_an_all_carbon_ring() {
        // Catechol skeleton: both O are donors, and a shift walk
        // needs complementary terminals — no group.
        let mut atoms: Vec<AtomRecord> = (0..6).map(|_| atom("C")).collect();
        atoms.push(atom("O"));
        atoms.push(atom("O"));
        let s = structure(
            atoms,
            vec![
                (0, 1, 1),
                (1, 2, 2),
                (2, 3, 1),
                (3, 4, 2),
                (4, 5, 1),
                (5, 0, 2),
                (1, 6, 1),
                (2, 7, 1),
            ],
        );
        let groups = mobile_groups(&s, &[0, 0, 0, 0, 0, 0, 1, 1], &[1, 2, 3, 4, 5, 6, 7, 8]);
        assert!(groups.is_empty());
    }
}
