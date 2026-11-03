//! The reference's balanced network (ichi_bns.c): the flow network
//! the tautomer templates consult through `bExistsAnyAltPath` before
//! registering a pair — atoms carry their valence slack as s-capacity,
//! bonds carry their input order as flow, registered endpoints hang off
//! a t-group vertex carrying the group's mobile hydrogens, and an
//! alternate path exists when the cap-aware mirror search (Kocay &
//! Stone balanced flow, with blossoms) can move one dot of H from the
//! single-bond side to the double-bond side.

use crate::graph::Structure;

pub(crate) const NO_VERTEX: i32 = -2;
const BLOSSOM_BASE: i32 = -1;
const TREE_IN_2: u8 = 1;
const TREE_IN_2BLOSS: u8 = 2;
const TREE_IN_1: u8 = 3;
const MAX_PASSES: usize = 16;
const MAX_BOND_EDGE_CAP: i32 = 2;

const fn prim(v: i32) -> i32 {
    v ^ 1
}

#[derive(Clone)]
struct Edge {
    a: u32,
    b: u32,
    cap: i32,
    flow: i32,
    path: bool,
    forbidden: bool,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Kind {
    Atom,
    CGroup,
    TGroup,
    Temp,
}

#[derive(Clone)]
struct Vertex {
    st_cap: i32,
    st_flow: i32,
    st_path: bool,
    /// Bond edges first (input order), then the c-group edge (if
    /// any), then the t-group edge; a temp vertex holds the single
    /// edge to its host.
    adj: Vec<u32>,
    kind: Kind,
}

/// The inputs one [`Network`] build needs: everything is expressed in
/// structure terms so this module stays free of the tautomer context.
pub(crate) struct Build<'a> {
    pub bond_adjacency: &'a [Vec<(usize, i32, bool)>],
    pub max_flow: &'a [i32],
    pub group_of: &'a [Option<usize>],
    pub donor: &'a [bool],
    /// Per atom, parallel to `bond_adjacency`: the bond reads
    /// single-like to the t-group capacity raise.
    pub single_like: &'a [Vec<bool>],
    /// The charge-group membership (at[].c_point) and which member
    /// carries the moveable charge (its c-group edge flows 0, the
    /// neutral members' flow 1).
    pub c_point_of: &'a [Option<usize>],
    pub c_charged: &'a [bool],
}

/// Which endpoints may play the donor (single-bond / H-owner) or
/// acceptor (double-bond) side of the TAUTOM precheck: registered
/// membership admits either role.
pub(crate) struct Roles {
    pub donor_ok: Vec<bool>,
    pub acceptor_ok: Vec<bool>,
}

#[derive(Clone)]
pub(crate) struct Network {
    verts: Vec<Vertex>,
    edges: Vec<Edge>,
}

impl Network {
    /// AllocateAndInitBnStruct + AddTGroups2BnStruct: atom vertices
    /// with input-order bond edges, then t-group vertices (ascending
    /// endpoint, the group edge appended last, donor edges flowing).
    pub(crate) fn build(spec: &Build<'_>) -> Network {
        let n = spec.max_flow.len();
        let mut verts: Vec<Vertex> = (0..n)
            .map(|i| Vertex {
                st_cap: spec.max_flow[i],
                st_flow: 0,
                st_path: false,
                adj: Vec::new(),
                kind: Kind::Atom,
            })
            .collect();
        let mut edges: Vec<Edge> = Vec::new();
        let mut made: std::collections::HashMap<(usize, usize), u32> =
            std::collections::HashMap::new();
        for i in 0..n {
            for &(neigh, flow, forbidden) in &spec.bond_adjacency[i] {
                if made.contains_key(&(i.min(neigh), i.max(neigh))) {
                    continue;
                }
                let cap = spec.max_flow[i]
                    .min(spec.max_flow[neigh])
                    .min(MAX_BOND_EDGE_CAP);
                let index = edges.len() as u32;
                edges.push(Edge {
                    a: i as u32,
                    b: neigh as u32,
                    cap,
                    flow,
                    path: false,
                    forbidden,
                });
                made.insert((i.min(neigh), i.max(neigh)), index);
                verts[i].st_flow += flow;
                verts[neigh].st_flow += flow;
                verts[i].adj.push(index);
                verts[neigh].adj.push(index);
            }
        }

        // AddCGroups2BnStruct runs before AddTGroups2BnStruct: the
        // c-group vertices sit between the atoms and the t-groups.
        let mut c_groups: Vec<usize> = spec.c_point_of.iter().flatten().copied().collect();
        c_groups.sort_unstable();
        c_groups.dedup();
        for _ in 0..c_groups.len() {
            verts.push(Vertex {
                st_cap: 0,
                st_flow: 0,
                st_path: false,
                adj: Vec::new(),
                kind: Kind::CGroup,
            });
        }
        let c_group_vertex = |gid: usize| n + c_groups.iter().position(|&g| g == gid).unwrap();
        for point in 0..n {
            let Some(gid) = spec.c_point_of.get(point).and_then(|g| *g) else {
                continue;
            };
            let host = c_group_vertex(gid);
            let charged = spec.c_charged[point];
            let index = edges.len() as u32;
            edges.push(Edge {
                a: point as u32,
                b: host as u32,
                cap: 1,
                flow: !charged as i32,
                path: false,
                forbidden: false,
            });
            verts[point].adj.push(index);
            verts[host].adj.push(index);
            if !charged {
                verts[host].st_cap += 1;
                verts[host].st_flow += 1;
                verts[point].st_cap += 1;
                verts[point].st_flow += 1;
            }
        }

        let mut groups: Vec<usize> = spec.group_of.iter().flatten().copied().collect();
        groups.sort_unstable();
        groups.dedup();
        for _ in 0..groups.len() {
            verts.push(Vertex {
                st_cap: 0,
                st_flow: 0,
                st_path: false,
                adj: Vec::new(),
                kind: Kind::TGroup,
            });
        }
        let group_vertex =
            |gid: usize| n + c_groups.len() + groups.iter().position(|&g| g == gid).unwrap();
        for endpoint in 0..n {
            let Some(gid) = spec.group_of.get(endpoint).and_then(|g| *g) else {
                continue;
            };
            let host = group_vertex(gid);
            // Capacity raise: a zero-capacity bond of the endpoint to a
            // neighbor with valence slack becomes movable when it reads
            // single-like.
            for (position, &(neigh, _, _)) in spec.bond_adjacency[endpoint].iter().enumerate() {
                let index = made[&(endpoint.min(neigh), endpoint.max(neigh))] as usize;
                if edges[index].cap == 0
                    && verts[neigh].st_cap >= 1
                    && spec.single_like[endpoint][position]
                {
                    edges[index].cap = 1;
                }
            }
            let donor = spec.donor[endpoint];
            let index = edges.len() as u32;
            edges.push(Edge {
                a: endpoint as u32,
                b: host as u32,
                cap: 1,
                flow: donor as i32,
                path: false,
                forbidden: false,
            });
            verts[endpoint].adj.push(index);
            verts[host].adj.push(index);
            if donor {
                verts[host].st_cap += 1;
                verts[host].st_flow += 1;
                verts[endpoint].st_cap += 1;
                verts[endpoint].st_flow += 1;
            }
        }
        Network { verts, edges }
    }

    /// The edge between two adjacent atoms, from the atom's adjacency.
    pub(crate) fn edge_between(&self, a: usize, b: usize) -> Option<usize> {
        self.verts[a]
            .adj
            .iter()
            .copied()
            .map(|edge| edge as usize)
            .find(|&edge| self.other(edge, a) == b)
    }

    pub(crate) fn cur_flow(&self, edge: usize) -> i32 {
        self.edges[edge].flow
    }

    pub(crate) fn edge_cap(&self, edge: usize) -> i32 {
        self.edges[edge].cap.min(MAX_BOND_EDGE_CAP)
    }

    /// nMinFlow2Check: how far the edge's flow can fall — the flow
    /// minus the smaller of the endpoints' residual capacity from
    /// their OTHER edges, clamped at 0.
    pub(crate) fn min_flow_2_check(&self, edge: usize) -> i32 {
        let edge_rec = &self.edges[edge];
        let f12 = edge_rec.flow;
        if f12 <= 0 {
            return 0;
        }
        let v1 = edge_rec.a as usize;
        let v2 = edge_rec.b as usize;
        let residual_around = |v: usize| -> i32 {
            self.verts[v]
                .adj
                .iter()
                .copied()
                .map(|e| e as usize)
                .filter(|&e| e != edge)
                .map(|e| self.edges[e].cap - self.edges[e].flow)
                .sum()
        };
        f12 - residual_around(v1).min(residual_around(v2)).min(f12)
    }

    fn other(&self, edge: usize, of: usize) -> usize {
        let edge = &self.edges[edge];
        edge.a as usize ^ edge.b as usize ^ of
    }

    /// GetVertexDegree: mirror vertices answer with their structure
    /// vertex degree plus the s/t edge when it has capacity.
    fn degree(&self, u: i32) -> i32 {
        let w = u / 2 - 1;
        if w >= 0 {
            let vert = &self.verts[w as usize];
            if vert.st_cap > 0 {
                vert.adj.len() as i32 + 1
            } else {
                0
            }
        } else {
            self.verts.len() as i32
        }
    }

    /// GetVertexNeighbor: `Some((neighbor, iedge))` with the pseudo
    /// s/t edge encoded as −(neighbor+1), `None` where the reference
    /// reports NO_VERTEX.
    fn neighbor(&self, u: i32, index: i32) -> Option<(i32, i32)> {
        let i = u - 2;
        if i >= 0 {
            if index > 0 {
                let edge = self.verts[(i / 2) as usize].adj[(index - 1) as usize] as usize;
                if self.edges[edge].cap == 0 || self.edges[edge].forbidden {
                    return None;
                }
                let n12 = (self.edges[edge].a ^ self.edges[edge].b) as i32;
                Some(((i ^ (2 * n12 + 1)) + 2, edge as i32))
            } else {
                let st = u & 1;
                Some((st, -(st + 1)))
            }
        } else {
            let w = index as usize;
            if self.verts.get(w)?.st_cap == 0 {
                return None;
            }
            let nb = 2 * w as i32 + 2 + (u & 1);
            Some((nb, -(nb + 1)))
        }
    }

    /// GetEdgePointer: which record backs the mirror edge and whether
    /// the traversal runs against its direction.
    fn edge_pointer(&self, u: i32, v: i32, iuv: i32) -> Option<(bool, usize, bool)> {
        // (st_edge?, record, backward)
        let i = u / 2 - 1;
        let j = v / 2 - 1;
        if i >= 0 {
            if j >= 0 {
                if (u + v) % 2 == 1 {
                    Some((false, iuv as usize, u & 1 == 1))
                } else {
                    None
                }
            } else if v >= 0 && (u + v) % 2 == 0 {
                Some((true, i as usize, v & 1 != 1))
            } else {
                None
            }
        } else if j >= 0 && u >= 0 && (u + v) % 2 == 0 {
            Some((true, j as usize, u & 1 == 1))
        } else {
            None
        }
    }

    fn rescap(&self, u: i32, v: i32, iuv: i32) -> i32 {
        let Some((st, record, backward)) = self.edge_pointer(u, v, iuv) else {
            return -1;
        };
        let (cap, flow) = if st {
            (self.verts[record].st_cap, self.verts[record].st_flow)
        } else {
            (self.edges[record].cap, self.edges[record].flow)
        };
        if backward {
            flow
        } else {
            cap - flow
        }
    }

    /// rescap_mark: the residual, halved on the second traversal of
    /// the same record (an edge shared by the path and its complement).
    fn rescap_mark(&mut self, u: i32, v: i32, iuv: i32) -> i32 {
        let Some((st, record, backward)) = self.edge_pointer(u, v, iuv) else {
            return -1;
        };
        if st {
            let vert = &mut self.verts[record];
            let mut f = if backward {
                vert.st_flow
            } else {
                vert.st_cap - vert.st_flow
            };
            if vert.st_path {
                f /= 2;
            } else {
                vert.st_path = true;
            }
            f
        } else {
            let edge = &mut self.edges[record];
            let mut f = if backward {
                edge.flow
            } else {
                edge.cap - edge.flow
            };
            if edge.path {
                f /= 2;
            } else {
                edge.path = true;
            }
            f
        }
    }

    fn augment_edge(&mut self, u: i32, v: i32, iuv: i32, mut delta: i32) -> i32 {
        let Some((st, record, backward)) = self.edge_pointer(u, v, iuv) else {
            return -1;
        };
        if backward {
            delta = -delta;
        }
        if st {
            let vert = &mut self.verts[record];
            let f = vert.st_flow + delta;
            if delta == 0 {
                vert.st_path = false;
                return f;
            }
            if f < 0 || f > vert.st_cap {
                return -1;
            }
            vert.st_flow = f;
            vert.st_path = false;
            f
        } else {
            let edge = &mut self.edges[record];
            let f = edge.flow + delta;
            if delta == 0 {
                edge.path = false;
                return f;
            }
            if f < 0 || f > edge.cap {
                return -1;
            }
            edge.flow = f;
            edge.path = false;
            f
        }
    }
}

fn tree_mark(tree: &mut [u8], v: i32, mark: u8) {
    let slot = &mut tree[v as usize];
    if *slot < mark {
        *slot = mark;
    }
}

fn s_reachable(tree: &[u8], v: i32) -> bool {
    tree[v as usize] >= TREE_IN_2BLOSS
}

fn find_base(base: &mut [i32], u: i32) -> i32 {
    match base[u as usize] {
        NO_VERTEX => NO_VERTEX,
        BLOSSOM_BASE => u,
        parent => {
            let root = find_base(base, parent);
            base[u as usize] = root;
            root
        }
    }
}

/// Get2ndEdgeVertex over the stored switch edge of `u`.
fn second_vertex(net: &Network, sw_v1: &[i32], sw_ie: &[i32], u: i32) -> i32 {
    let v1 = sw_v1[u as usize];
    let ie = sw_ie[u as usize];
    if ie >= 0 {
        let n12 = (net.edges[ie as usize].a ^ net.edges[ie as usize].b) as i32;
        ((v1 - 2) ^ (2 * n12 + 1)) + 2
    } else if v1 <= 1 {
        -(1 + ie)
    } else {
        v1 % 2
    }
}

fn find_path_to_vertex_s(
    start: i32,
    sw_v1: &[i32],
    base: &mut [i32],
    path: &mut [i32],
) -> Option<usize> {
    let mut i = 0usize;
    path[0] = start;
    let mut cur = start;
    while cur != 0 {
        cur = find_base(base, sw_v1[cur as usize]);
        i += 1;
        if i >= path.len() {
            return None;
        }
        path[i] = cur;
    }
    Some(i)
}

#[allow(clippy::too_many_arguments)]
fn make_blossom(
    net: &mut Network,
    scanq: &mut [i32],
    qsize: &mut i32,
    pu: &mut [i32],
    pv: &mut [i32],
    sw_v1: &mut [i32],
    sw_ie: &mut [i32],
    base: &mut [i32],
    tree: &mut [u8],
    u: i32,
    v: i32,
    iuv: i32,
    b_u: i32,
    b_v: i32,
) -> Option<i32> {
    let len_pu = find_path_to_vertex_s(b_u, sw_v1, base, pu)?;
    let len_pv = find_path_to_vertex_s(b_v, sw_v1, base, pv)?;
    let mut i = len_pu as i32;
    let mut j = len_pv as i32;
    while i >= 0 && j >= 0 && pu[i as usize] == pv[j as usize] {
        i -= 1;
        j -= 1;
    }
    i += 1;
    let mut w = pu[i as usize];
    let (mut z, mut izw) = (sw_v1[w as usize], sw_ie[w as usize]);
    while w != 0 && net.rescap(z, w, izw) >= 2 {
        i += 1;
        w = pu[i as usize];
        z = sw_v1[w as usize];
        izw = sw_ie[w as usize];
    }
    let mut k = i - 1;
    while k >= 0 {
        let vertex = pu[k as usize];
        base[vertex as usize] = w;
        base[prim(vertex) as usize] = w;
        let p = prim(vertex);
        if !s_reachable(tree, p) {
            sw_v1[p as usize] = prim(v);
            sw_ie[p as usize] = iuv;
            *qsize += 1;
            scanq[*qsize as usize] = p;
            tree_mark(tree, p, TREE_IN_2BLOSS);
        }
        k -= 1;
    }
    k = j;
    while k >= 0 {
        let vertex = pv[k as usize];
        base[vertex as usize] = w;
        base[prim(vertex) as usize] = w;
        let p = prim(vertex);
        if !s_reachable(tree, p) {
            sw_v1[p as usize] = u;
            sw_ie[p as usize] = iuv;
            *qsize += 1;
            scanq[*qsize as usize] = p;
            tree_mark(tree, p, TREE_IN_2BLOSS);
        }
        k -= 1;
    }
    let p = prim(w);
    if !s_reachable(tree, p) {
        sw_v1[p as usize] = u;
        sw_ie[p as usize] = iuv;
        *qsize += 1;
        scanq[*qsize as usize] = p;
        tree_mark(tree, p, TREE_IN_2BLOSS);
    }
    Some(w)
}

fn find_path_cap(
    net: &mut Network,
    sw_v1: &[i32],
    sw_ie: &[i32],
    x: i32,
    y: i32,
    delta: i32,
) -> i32 {
    let w = sw_v1[y as usize];
    let z = second_vertex(net, sw_v1, sw_ie, y);
    let iwz = sw_ie[y as usize];
    let cap = net.rescap_mark(w, z, iwz);
    if cap < 0 {
        return cap;
    }
    let mut delta = delta.min(cap);
    if w != x {
        delta = delta.min(find_path_cap(net, sw_v1, sw_ie, x, w, delta));
    }
    if z != y {
        delta = delta.min(find_path_cap(net, sw_v1, sw_ie, prim(y), prim(z), delta));
    }
    delta
}

fn pull_flow(
    net: &mut Network,
    sw_v1: &[i32],
    sw_ie: &[i32],
    x: i32,
    y: i32,
    delta: i32,
    b_reverse: bool,
) -> i32 {
    let w = sw_v1[y as usize];
    let z = second_vertex(net, sw_v1, sw_ie, y);
    let iwz = sw_ie[y as usize];
    let mut ret = 0;
    if b_reverse {
        if z != y {
            ret = pull_flow(net, sw_v1, sw_ie, prim(y), prim(z), delta, !b_reverse);
        }
        if ret >= 0 {
            ret = net.augment_edge(w, z, iwz, delta);
        }
        if w != x && ret >= 0 {
            ret = pull_flow(net, sw_v1, sw_ie, x, w, delta, b_reverse);
        }
    } else {
        if w != x && ret >= 0 {
            ret = pull_flow(net, sw_v1, sw_ie, x, w, delta, b_reverse);
        }
        if ret >= 0 {
            ret = net.augment_edge(w, z, iwz, delta);
        }
        if z != y && ret >= 0 {
            ret = pull_flow(net, sw_v1, sw_ie, prim(y), prim(z), delta, !b_reverse);
        }
    }
    ret
}

/// BalancedNetworkSearch: one augmenting pass over the mirror network;
/// the return is the path capacity, 0 when exhausted, negative on the
/// reference's parameter errors (treated as absence by the caller).
fn balanced_network_search(net: &mut Network) -> i32 {
    let size = 2 * net.verts.len() + 2;
    let mut base = vec![NO_VERTEX; size];
    let mut sw_v1 = vec![NO_VERTEX; size];
    let mut sw_ie = vec![NO_VERTEX; size];
    let mut tree = vec![0u8; size];
    let mut scanq = vec![0i32; size];
    let mut pu = vec![0i32; size];
    let mut pv = vec![0i32; size];

    let mut qsize: i32 = 0;
    scanq[0] = 0;
    base[1] = 0;
    base[0] = BLOSSOM_BASE;
    tree[0] = TREE_IN_1;
    let mut k: i32 = 0;
    loop {
        let u = scanq[k as usize];
        let mut b_u = find_base(&mut base, u);
        let degree = net.degree(u);
        for i in 0..degree {
            let Some((v, iuv)) = net.neighbor(u, i) else {
                continue;
            };
            let v1u = sw_v1[u as usize];
            if !(v1u != v || second_vertex(net, &sw_v1, &sw_ie, u) != u) {
                continue;
            }
            let r = net.rescap(u, v, iuv);
            if r <= 0 {
                if r < 0 {
                    return 0;
                }
                continue;
            }
            let b_v = find_base(&mut base, v);
            if b_v == NO_VERTEX {
                qsize += 1;
                scanq[qsize as usize] = v;
                tree_mark(&mut tree, v, TREE_IN_1);
                tree_mark(&mut tree, prim(v), TREE_IN_2);
                sw_v1[v as usize] = u;
                sw_ie[v as usize] = iuv;
                base[prim(v) as usize] = v;
                base[v as usize] = BLOSSOM_BASE;
            } else {
                let pu_edge = sw_v1[prim(u) as usize];
                let t_reachable = s_reachable(&tree, prim(v))
                    && !(pu_edge == prim(v)
                        && second_vertex(net, &sw_v1, &sw_ie, prim(u)) == prim(u))
                    && b_u != b_v;
                if t_reachable {
                    let Some(w) = make_blossom(
                        net, &mut scanq, &mut qsize, &mut pu, &mut pv, &mut sw_v1, &mut sw_ie,
                        &mut base, &mut tree, u, v, iuv, b_u, b_v,
                    ) else {
                        return 0;
                    };
                    b_u = w;
                    if prim(w) == 1 {
                        let delta = find_path_cap(net, &sw_v1, &sw_ie, 0, 1, 10000);
                        if delta < 0 {
                            return 0;
                        }
                        let ret = pull_flow(net, &sw_v1, &sw_ie, 0, 1, delta, false);
                        return if ret < 0 { 0 } else { delta };
                    }
                }
            }
        }
        k += 1;
        if k > qsize {
            break;
        }
    }
    0
}

/// bSetFlowToCheckOneBond: pin the edge at `flow` by absorbing the
/// difference from the endpoints' residuals and neighbors; the return
/// is the dot balance the operation creates (negative = failure).
fn set_flow_to_check_one_bond(net: &mut Network, iedge: usize, flow: i32) -> Option<i32> {
    let f12 = net.edges[iedge].flow;
    let v1 = net.edges[iedge].a as usize;
    let v2 = net.edges[iedge].b as usize;
    let mut ndots = 0;
    if f12 < flow {
        if net.verts[v1].st_cap < flow || net.verts[v2].st_cap < flow {
            return None;
        }
        if net.verts[v1].st_flow < f12 || net.verts[v2].st_flow < f12 {
            return None;
        }
        let (mut d1, mut d2) = (flow - f12, flow - f12);
        if f12 > 0 {
            for v in [v1, v2] {
                net.verts[v].st_cap -= f12;
                net.verts[v].st_flow -= f12;
            }
        }
        net.edges[iedge].flow = 0;
        net.edges[iedge].cap = 0;

        let mut rescap = net.verts[v1].st_cap - net.verts[v1].st_flow;
        while rescap > 0 && d1 > 0 {
            rescap -= 1;
            d1 -= 1;
            net.verts[v1].st_cap -= 1;
            ndots -= 1;
        }
        let mut rescap = net.verts[v2].st_cap - net.verts[v2].st_flow;
        while rescap > 0 && d2 > 0 {
            rescap -= 1;
            d2 -= 1;
            net.verts[v2].st_cap -= 1;
            ndots -= 1;
        }

        for &(host, delta) in &[(v1, d1), (v2, d2)] {
            let mut delta = delta;
            let adjacency: Vec<usize> = net.verts[host]
                .adj
                .iter()
                .copied()
                .map(|e| e as usize)
                .collect();
            for &edge in &adjacency {
                if delta == 0 {
                    break;
                }
                if edge == iedge || net.edges[edge].forbidden {
                    continue;
                }
                let mut f = net.edges[edge].flow;
                let neighbor = net.other(edge, host);
                while f > 0 && delta > 0 {
                    f -= 1;
                    delta -= 1;
                    net.edges[edge].flow -= 1;
                    net.verts[neighbor].st_flow -= 1;
                    net.verts[host].st_cap -= 1;
                    net.verts[host].st_flow -= 1;
                    ndots += 1;
                }
            }
            if delta > 0 {
                return None;
            }
        }
    } else {
        if net.verts[v1].st_flow < f12
            || net.verts[v2].st_flow < f12
            || net.verts[v1].st_cap < flow
            || net.verts[v2].st_cap < flow
        {
            return None;
        }
        let delta = f12 - flow;
        for v in [v1, v2] {
            net.verts[v].st_flow -= f12;
            net.verts[v].st_cap -= flow;
        }
        net.edges[iedge].flow = 0;
        net.edges[iedge].cap = 0;
        ndots = 2 * delta;
    }
    Some(ndots)
}

/// bAddNewVertex: attach a temp vertex carrying one dot to `host`.
fn add_new_vertex(net: &mut Network, host: usize, ndots: &mut i32) {
    add_new_vertex_cap(net, host, 1, ndots);
}

/// bAddStCapToAVertex: one more dot on `v1` and the raised bond
/// capacities toward every non-group neighbor.
fn add_st_cap(net: &mut Network, v1: usize, skip: usize, ndots: &mut i32) {
    net.verts[v1].st_cap += 1;
    *ndots += 1;
    // The reference raises no edge caps from a group vertex — the
    // whole increment block sits under `if (!(pVert1->type & ANY_GROUP))`.
    if net.verts[v1].kind != Kind::Atom {
        return;
    }
    let cap = net.verts[v1].st_cap;
    let adj = net.verts[v1].adj.clone();
    for edge in adj {
        let edge = edge as usize;
        let other = net.other(edge, v1);
        if other == skip || net.verts[other].kind != Kind::Atom {
            continue;
        }
        net.edges[edge].cap = net.verts[other].st_cap.min(cap).min(MAX_BOND_EDGE_CAP);
    }
}

/// One BnsTestAndMarkAltBonds probe: pin the edge at `test_flow`,
/// run the augmentation, and report feasibility (2*sum == dots, or
/// the trivially-pinned dots==0) with the bonds whose flow moved.
pub(crate) struct FlowProbe {
    pub feasible: bool,
    pub moved: Vec<(usize, usize)>,
}

pub(crate) fn probe_altern_flow(net: &Network, iedge: usize, test_flow: i32) -> FlowProbe {
    let before: Vec<i32> = net.edges.iter().map(|edge| edge.flow).collect();
    let mut work = net.clone();
    let dots = match set_flow_to_check_one_bond(&mut work, iedge, test_flow) {
        Some(dots) => dots,
        None => {
            return FlowProbe {
                feasible: false,
                moved: Vec::new(),
            }
        }
    };
    if dots == 0 {
        return FlowProbe {
            feasible: true,
            moved: Vec::new(),
        };
    }
    if dots < 0 {
        return FlowProbe {
            feasible: false,
            moved: Vec::new(),
        };
    }
    let mut sum = 0;
    let mut feasible = false;
    for _ in 0..MAX_PASSES {
        let delta = balanced_network_search(&mut work);
        if delta <= 0 {
            feasible = delta == 0 && 2 * sum == dots;
            break;
        }
        sum += delta;
    }
    let moved = if feasible {
        (0..net.edges.len())
            .filter(|&index| work.edges[index].flow != before[index])
            .map(|index| {
                let edge = &net.edges[index];
                (edge.a.min(edge.b) as usize, edge.a.max(edge.b) as usize)
            })
            .collect()
    } else {
        Vec::new()
    };
    FlowProbe { feasible, moved }
}

/// The t-group or c-group edge of an atom: GetEdgeToGroupVertex
/// scans the adjacency from the end for a neighbor of the exact type.
fn group_edge_of(net: &Network, atom: usize, kind: Kind) -> Option<usize> {
    let vert = &net.verts[atom];
    if vert.kind != Kind::Atom {
        return None;
    }
    vert.adj
        .iter()
        .rev()
        .copied()
        .map(|edge| edge as usize)
        .find(|&edge| net.verts[net.other(edge, atom)].kind == kind)
}

/// The t-group vertex acting for `atom`, when the atom is an
/// endpoint (GetGroupVertex).
fn group_vertex_of(net: &Network, atom: usize, kind: Kind) -> Option<usize> {
    group_edge_of(net, atom, kind).map(|edge| net.other(edge, atom))
}

/// bAddNewVertex with an explicit dot capacity: a fresh `cap`-sized
/// sink edge on `host` and the host's st-cap rebalance.
fn add_new_vertex_cap(net: &mut Network, host: usize, cap: i32, ndots: &mut i32) {
    let new = net.verts.len();
    net.verts.push(Vertex {
        st_cap: cap,
        st_flow: 0,
        st_path: false,
        adj: Vec::new(),
        kind: Kind::Temp,
    });
    let index = net.edges.len() as u32;
    net.edges.push(Edge {
        a: host as u32,
        b: new as u32,
        cap,
        flow: 0,
        path: false,
        forbidden: false,
    });
    net.verts[host].adj.push(index);
    net.verts[new].adj.push(index);
    *ndots += 1;
    *ndots -= net.verts[host].st_cap - net.verts[host].st_flow;
    if net.verts[host].st_cap < net.verts[host].st_flow {
        net.verts[host].st_cap = net.verts[host].st_flow;
    }
    *ndots += net.verts[host].st_cap - net.verts[host].st_flow;
}

/// bExistsAltPath's 2-H setup (ALT_PATH_MODE_{ADD2H,REM2H}_{TST,CHG}):
/// REM2H endpoints push (st-cap increments, or a dot vertex when they
/// act through a t-group), ADD2H endpoints pull (dot vertices, or
/// st-cap when grouped); same-group endpoints share one 2-unit dot.
/// An alternate path exists when the augmentation covers the dots.
pub(crate) fn exists_2h(net: &mut Network, first: usize, second: usize, donors: bool) -> bool {
    // Adjacency precheck per endpoint: donors need a neighbor with st
    // capacity (or the partner over a non-maxed bond); acceptors need
    // a neighbor carrying flow.
    let capable = |net: &Network, v: usize, partner: usize| -> bool {
        let mut n = 0;
        let mut adjacent = 0;
        for &edge in &net.verts[v].adj {
            let edge = edge as usize;
            let other = net.other(edge, v);
            if donors {
                n += (net.verts[other].st_cap > 0) as i32;
                adjacent += (other == partner && net.edges[edge].cap < MAX_BOND_EDGE_CAP) as i32;
            } else {
                n += (net.edges[edge].flow > 0) as i32;
            }
        }
        n > 0 || adjacent > 0
    };
    if !capable(net, first, second) || !capable(net, second, first) {
        return false;
    }

    let mut ndots = 0;
    if std::env::var_os("INCHI_TRACE_NET").is_some() {
        for (atom, vert) in net.verts.iter().enumerate() {
            if vert.kind != Kind::Atom {
                continue;
            }
            let edges: Vec<String> = vert
                .adj
                .iter()
                .map(|&e| {
                    format!(
                        "{}:{}/{}",
                        net.other(e as usize, atom),
                        net.edges[e as usize].cap,
                        net.edges[e as usize].flow
                    )
                })
                .collect();
            eprintln!(
                "NET atom={atom} st={}/{} adj={}",
                vert.st_cap,
                vert.st_flow,
                edges.join(",")
            );
        }
    }
    let v1t = group_vertex_of(net, first, Kind::TGroup);
    let v2t = group_vertex_of(net, second, Kind::TGroup);
    let v1act = v1t.unwrap_or(first);
    let v2act = v2t.unwrap_or(second);
    let cap = if v1t.is_some() && v1t == v2t { 2 } else { 1 };

    // (bDonors != 0) == (vt != NO_VERTEX) picks the dot side.
    let side = |grouped: bool| donors == grouped;
    let mut set_v1 = false;
    let mut set_v2 = false;
    if side(v1t.is_some()) {
        add_new_vertex_cap(net, v1act, cap, &mut ndots);
        set_v1 = true;
    }
    if side(v2t.is_some()) {
        if cap == 1 {
            add_new_vertex_cap(net, v2act, 1, &mut ndots);
        }
        set_v2 = true;
    }
    if !set_v1 {
        add_st_cap(net, v1act, v2act, &mut ndots);
    }
    if !set_v2 {
        add_st_cap(net, v2act, v1act, &mut ndots);
    }
    if ndots < 0 || ndots % 2 != 0 {
        return false;
    }

    let mut sum = 0;
    for _ in 0..MAX_PASSES {
        let delta = balanced_network_search(net);
        if delta <= 0 {
            let ok = delta == 0 && 2 * sum >= ndots;
            if std::env::var_os("INCHI_TRACE_BNS").is_some() {
                eprintln!("SALT2H {first}-{second} donors={donors} dots={ndots} sum={sum} delta={delta} -> {ok}");
            }
            return ok;
        }
        sum += delta;
    }
    let ok = 2 * sum >= ndots;
    if std::env::var_os("INCHI_TRACE_BNS").is_some() {
        eprintln!(
            "SALT2H {first}-{second} donors={donors} dots={ndots} sum={sum} passes-out -> {ok}"
        );
    }
    ok
}

/// bExistsAltPath under ALT_PATH_MODE_TAUTOM: the role precheck, the
/// setup that pins the two endpoints' group flows and hangs the H dot,
/// then up to MAX_PASSES augmenting searches — an alternate path
/// exists when the eliminations cover the planted dots.
fn exists(net: &mut Network, roles: &Roles, double_side: usize, single_side: usize) -> bool {
    if !roles.donor_ok[single_side] || !roles.acceptor_ok[double_side] {
        return false;
    }
    if net.verts[double_side].st_flow == 0 {
        return false;
    }
    let mut ndots = 0;
    // Prohibit charge movement: pin each side's c-group edge at flow 1.
    if let Some(edge) = group_edge_of(net, single_side, Kind::CGroup) {
        let Some(ret) = set_flow_to_check_one_bond(net, edge, 1) else {
            return false;
        };
        ndots += ret;
    }
    if let Some(edge) = group_edge_of(net, double_side, Kind::CGroup) {
        let Some(ret) = set_flow_to_check_one_bond(net, edge, 1) else {
            return false;
        };
        ndots += ret;
    }
    if let Some(edge) = group_edge_of(net, single_side, Kind::TGroup) {
        let Some(ret) = set_flow_to_check_one_bond(net, edge, 1) else {
            return false;
        };
        ndots += ret;
    }
    if let Some(edge) = group_edge_of(net, double_side, Kind::TGroup) {
        let Some(ret) = set_flow_to_check_one_bond(net, edge, 0) else {
            return false;
        };
        ndots += ret;
    }
    add_new_vertex(net, double_side, &mut ndots);
    add_st_cap(net, single_side, double_side, &mut ndots);
    if ndots < 0 || ndots % 2 != 0 {
        return false;
    }

    let mut sum = 0;
    for _ in 0..MAX_PASSES {
        let delta = balanced_network_search(net);
        if delta <= 0 {
            return delta == 0 && 2 * sum >= ndots;
        }
        sum += delta;
    }
    2 * sum >= ndots
}

/// bExistsAltPath under ALT_PATH_MODE_REM_PROTON as the instrumented
/// oracle runs it on the acid-removal side (HardRemoveAcidicProtons):
/// BOTH endpoints are the virtual t-group stars, not atoms — the
/// H-donor star and the acidic =O/N-O star built over the whole
/// structure. One dot vertex lands on the donor-side star
/// (bAddNewVertex, cap 1) and one st-cap raise lands on the
/// acceptor-side star (bAddStCapToAVertex); the removal is viable
/// when the balanced augmentation covers the planted dots — the
/// donor proton can ride the skeleton onto the acidic oxygen, which
/// is where the reference then strips it as an acidic proton.
/// The REM_PROTON search with the alternation it moved: the atom-atom
/// bonds whose flow changed, signed by the π delta the reference's
/// BNS_EF_SAVE_ALL|BNS_EF_ALTR_NS ride keeps (positive = the bond
/// gained a π unit, negative = lost one). Star-member edges are
/// dropped — only real bonds shift orders.
pub(crate) fn rem_proton_route_moved(
    net: &Network,
    donor: usize,
    acceptor: usize,
) -> (bool, Vec<(usize, usize, i32)>) {
    let donor_star = group_vertex_of(net, donor, Kind::TGroup);
    let acid_star = group_vertex_of(net, acceptor, Kind::TGroup);
    let (Some(donor_star), Some(acid_star)) = (donor_star, acid_star) else {
        return (false, Vec::new());
    };
    if donor_star == acid_star {
        return (false, Vec::new());
    }
    // The reference's prechecks: the dot side must carry flow, and a
    // group vertex whose edges all read doubled flow is refused
    // (so far a group edge moves at most one unit per member).
    if net.verts[donor_star].st_flow == 0 {
        return (false, Vec::new());
    }
    if net.verts[acid_star].st_flow == 2 * net.verts[acid_star].adj.len() as i32 {
        return (false, Vec::new());
    }
    let before: Vec<i32> = net.edges.iter().map(|edge| edge.flow).collect();
    let mut work = net.clone();
    let mut ndots = 0;
    add_new_vertex_cap(&mut work, donor_star, 1, &mut ndots);
    add_st_cap(&mut work, acid_star, donor_star, &mut ndots);
    if ndots < 0 || ndots % 2 != 0 {
        return (false, Vec::new());
    }
    let mut sum = 0;
    let ok = loop {
        let delta = balanced_network_search(&mut work);
        if delta <= 0 {
            break delta == 0 && 2 * sum >= ndots;
        }
        sum += delta;
    };
    let moved: Vec<(usize, usize, i32)> = (0..net.edges.len())
        .filter(|&index| work.edges[index].flow != before[index])
        .filter_map(|index| {
            let edge = &net.edges[index];
            let (a, b) = (edge.a as usize, edge.b as usize);
            // Atom-atom bonds only; the group/star edges are the
            // virtual plumbing.
            (net.verts[a].kind == Kind::Atom && net.verts[b].kind == Kind::Atom)
                .then(|| (a.min(b), a.max(b), work.edges[index].flow - before[index]))
        })
        .collect();
    if std::env::var_os("INCHI_TRACE_PROTON").is_some() {
        eprintln!(
            "REM_PROTON donor={donor} acceptor={acceptor} dots={ndots} sum={sum} -> {ok} moved={moved:?}"
        );
    }
    (ok, moved)
}

/// The ride result over the two virtual stars: the H deltas on the
/// t-group star edges, the charge deltas on the c-group star edges,
/// and the atom-atom bond flow deltas (the pi alternation the ride
/// keeps — BNS_EF_SAVE_ALL | BNS_EF_ALTR_NS).
#[derive(Debug, Default)]
pub(crate) struct StarRide {
    pub route: bool,
    pub tg: Vec<(usize, i32)>,
    pub cg: Vec<(usize, i32)>,
    pub bonds: Vec<(usize, usize, i32)>,
}

/// The HardRemoveHplusNP network inputs: the skeleton at the standard
/// flow convention plus the three stars in the reference's creation
/// order — the (+) charge group, the (-) charge group, the PR_HARD
/// H-group (vertex ids follow that order).
pub(crate) struct HnpSpec {
    /// Atom-atom bonds in input order: (a, b, flow) with the flow the
    /// bond order reads (single 0, double 1, triple 2).
    pub bonds: Vec<(usize, usize, i32)>,
    /// Per atom: MAX_AT_FLOW — the st capacity.
    pub max_flow: Vec<i32>,
    /// The (+) group members: (atom, carries-the-plus).
    pub plus: Vec<(usize, bool)>,
    /// The (-) group members (charge 0 or -1 under the NEG classes).
    pub minus: Vec<(usize, bool)>,
    /// The H-group members: (atom, edge cap, edge flow).
    pub h_star: Vec<(usize, i32, i32)>,
}

/// AllocateAndInitBnStruct + CreateCGroupInBnStruct(+) +
/// CreateCGroupInBnStruct(-) + CreateTGroupInBnStruct, in that order,
/// each with its member-edge bookkeeping and the zero-capacity bond
/// raises the reference runs at creation time. Returns the network
/// and the (t-group star, plus star) vertex pair.
pub(crate) fn build_hnp(spec: &HnpSpec) -> (Network, usize, usize) {
    let n = spec.max_flow.len();
    let mut verts: Vec<Vertex> = (0..n)
        .map(|i| Vertex {
            st_cap: spec.max_flow[i],
            st_flow: 0,
            st_path: false,
            adj: Vec::new(),
            kind: Kind::Atom,
        })
        .collect();
    let mut edges: Vec<Edge> = Vec::new();
    let mut made: std::collections::HashMap<(usize, usize), u32> = std::collections::HashMap::new();
    for i in 0..n {
        for &(a, b, flow) in &spec.bonds {
            if a != i && b != i {
                continue;
            }
            let neigh = if a == i { b } else { a };
            if made.contains_key(&(i.min(neigh), i.max(neigh))) {
                continue;
            }
            let cap = spec.max_flow[i]
                .min(spec.max_flow[neigh])
                .min(MAX_BOND_EDGE_CAP);
            let index = edges.len() as u32;
            edges.push(Edge {
                a: i as u32,
                b: neigh as u32,
                cap,
                flow,
                path: false,
                forbidden: false,
            });
            made.insert((i.min(neigh), i.max(neigh)), index);
            verts[i].st_flow += flow;
            verts[neigh].st_flow += flow;
            verts[i].adj.push(index);
            verts[neigh].adj.push(index);
        }
    }
    // The two charge-group stars, then the H-group star (the
    // reference's creation order fixes the vertex ids).
    for kind in [Kind::CGroup, Kind::CGroup, Kind::TGroup] {
        verts.push(Vertex {
            st_cap: 0,
            st_flow: 0,
            st_path: false,
            adj: Vec::new(),
            kind,
        });
    }
    let plus_star = n;
    let minus_star = n + 1;
    let tg_star = n + 2;
    // CreateCGroupInBnStruct: member edges cap 1; the edge flows 1
    // unless the member carries the star's own charge — a (+) star's
    // charged member flows 0, a (-) star's charged member flows 1 —
    // and the st books count the flow on both ends, then the member's
    // zero-capacity bonds raise.
    let connect_c_group =
        |verts: &mut Vec<Vertex>, edges: &mut Vec<Edge>, star: usize, members: &[(usize, bool)]| {
            for &(atom, charged) in members {
                let index = edges.len() as u32;
                let flow = i32::from(!charged);
                edges.push(Edge {
                    a: atom as u32,
                    b: star as u32,
                    cap: 1,
                    flow,
                    path: false,
                    forbidden: false,
                });
                verts[atom].adj.push(index);
                verts[star].adj.push(index);
                if flow > 0 {
                    verts[star].st_cap += 1;
                    verts[star].st_flow += 1;
                    verts[atom].st_cap += 1;
                    verts[atom].st_flow += 1;
                }
                raise_zero_cap_bonds(verts, edges, atom);
            }
        };
    connect_c_group(&mut verts, &mut edges, plus_star, &spec.plus);
    connect_c_group(&mut verts, &mut edges, minus_star, &spec.minus);
    // CreateTGroupInBnStruct: member edges carry the per-member
    // capacity and the capped hydrogen flow.
    for &(atom, cap, flow) in &spec.h_star {
        let index = edges.len() as u32;
        edges.push(Edge {
            a: atom as u32,
            b: tg_star as u32,
            cap,
            flow,
            path: false,
            forbidden: false,
        });
        verts[atom].adj.push(index);
        verts[tg_star].adj.push(index);
        verts[tg_star].st_cap += flow;
        verts[tg_star].st_flow += flow;
        verts[atom].st_cap += flow;
        verts[atom].st_flow += flow;
    }
    for &(atom, _, _) in &spec.h_star {
        raise_zero_cap_bonds(&mut verts, &mut edges, atom);
    }
    (Network { verts, edges }, tg_star, plus_star)
}

/// The FIX_CPOINT_BOND_CAP / CreateTGroupInBnStruct bond raise: a
/// member's zero-capacity bond opens to min(neighbor st-cap, member
/// st-cap, cap bound) over the member's current st-cap.
fn raise_zero_cap_bonds(verts: &mut [Vertex], edges: &mut [Edge], atom: usize) {
    let member_cap = verts[atom].st_cap;
    let bond_edges: Vec<usize> = verts[atom]
        .adj
        .iter()
        .rev()
        .copied()
        .map(|edge| edge as usize)
        .filter(|&edge| {
            let edge = &edges[edge];
            let other = if edge.a as usize == atom {
                edge.b as usize
            } else {
                edge.a as usize
            };
            verts[other].kind == Kind::Atom
        })
        .collect();
    for edge in bond_edges {
        if edges[edge].cap == 0 {
            let edge_ref = &edges[edge];
            let neigh = if edge_ref.a as usize == atom {
                edge_ref.b as usize
            } else {
                edge_ref.a as usize
            };
            let raised = verts[neigh].st_cap.min(member_cap).min(MAX_BOND_EDGE_CAP);
            edges[edge].cap = raised;
        }
    }
}

/// bExistsAltPath under ALT_PATH_MODE_REM_PROTON between the freshly
/// created HardRemoveHplusNP stars (the PR_HARD t-group as the donor
/// side, the (+) charge group as the acceptor side): one dot vertex
/// on the t-group star, one st-cap raise on the charge-group star,
/// and every flow change the augmentation made reports per site.
pub(crate) fn star_ride(net: &Network, donor_star: usize, acceptor_star: usize) -> StarRide {
    let mut ride = StarRide::default();
    if net.verts[donor_star].st_flow == 0 {
        return ride;
    }
    if net.verts[acceptor_star].st_flow == 2 * net.verts[acceptor_star].adj.len() as i32 {
        return ride;
    }
    let before: Vec<i32> = net.edges.iter().map(|edge| edge.flow).collect();
    let mut work = net.clone();
    let mut ndots = 0;
    add_new_vertex_cap(&mut work, donor_star, 1, &mut ndots);
    add_st_cap(&mut work, acceptor_star, donor_star, &mut ndots);
    if ndots < 0 || ndots % 2 != 0 {
        return ride;
    }
    let mut sum = 0;
    let ok = loop {
        let delta = balanced_network_search(&mut work);
        if delta <= 0 {
            break delta == 0 && 2 * sum >= ndots;
        }
        sum += delta;
    };
    ride.route = ok;
    if !ok {
        return ride;
    }
    for (index, edge) in net.edges.iter().enumerate() {
        let delta = work.edges[index].flow - before[index];
        if delta == 0 {
            continue;
        }
        let (a, b) = (edge.a as usize, edge.b as usize);
        match (net.verts[a].kind, net.verts[b].kind) {
            (Kind::Atom, Kind::Atom) => ride.bonds.push((a.min(b), a.max(b), delta)),
            (Kind::Atom, star) => {
                if star == Kind::TGroup {
                    ride.tg.push((a, delta));
                } else if star == Kind::CGroup {
                    ride.cg.push((a, delta));
                }
            }
            (star, Kind::Atom) => {
                if star == Kind::TGroup {
                    ride.tg.push((b, delta));
                } else if star == Kind::CGroup {
                    ride.cg.push((b, delta));
                }
            }
            _ => {}
        }
    }
    ride
}

/// bExistsAltPath under ALT_PATH_MODE_CHARGE: no role precheck, the
/// generic setup (the H dot on the charged side, an extra st-cap on
/// the neutral side), and the changed atom-atom bonds reported for
/// the BOND_ALTERED marks. Success mirrors the reference's return:
/// bit 0 = found, and `3` (bits 0|1) when bonds moved.
pub(crate) fn charge_alt_path(
    net: &Network,
    charged: usize,
    neutral: usize,
) -> (bool, Vec<(usize, usize)>) {
    let before: Vec<i32> = net.edges.iter().map(|edge| edge.flow).collect();
    let mut work = net.clone();
    let found = exists_generic(&mut work, charged, neutral);
    let moved: Vec<(usize, usize)> = (0..net.edges.len())
        .filter(|&index| work.edges[index].flow != before[index])
        .map(|index| {
            let edge = &net.edges[index];
            (edge.a.min(edge.b) as usize, edge.a.max(edge.b) as usize)
        })
        .collect();
    (found, moved)
}

/// The generic-mode existence check: no role precheck, the temp dot
/// on `double_side`, the extra st-cap on `single_side`.
fn exists_generic(net: &mut Network, double_side: usize, single_side: usize) -> bool {
    if net.verts[double_side].st_flow == 0 {
        return false;
    }
    let mut ndots = 0;
    add_new_vertex(net, double_side, &mut ndots);
    add_st_cap(net, single_side, double_side, &mut ndots);
    if ndots < 0 || ndots % 2 != 0 {
        return false;
    }
    let mut sum = 0;
    for _ in 0..MAX_PASSES {
        let delta = balanced_network_search(net);
        if delta <= 0 {
            return delta == 0 && 2 * sum >= ndots;
        }
        sum += delta;
    }
    2 * sum >= ndots
}

/// bExistsAnyAltPath: either direction of the (double, single) pair.
pub(crate) fn any_alt_path(net: &Network, roles: &Roles, first: usize, second: usize) -> bool {
    let mut forward = net.clone();
    if exists(&mut forward, roles, first, second) {
        return true;
    }
    let mut backward = net.clone();
    exists(&mut backward, roles, second, first)
}

// ─── The resolved integer bond valence ──────────────────────────────

/// The reference's read-back of the BNS flow solution
/// (SetAtomRadAndChemValFromVertexCapFlow): after normalization every
/// atom's `chem_bonds_valence` settles at `valence + st_edge.flow` —
/// an integer. Real single/double/triple orders stand unchanged; the
/// aromatic bonds the input left alternating resolve one π per matched
/// bond, where the resolution is a maximum matching restricted to the
/// atoms whose preferred registered valence still lacks one unit after
/// their σ bonds and real multiples (witnessed by the instrumented
/// oracle: furan's O resolves π-free at cbv 2 while its four carbons
/// take one π each, naphthalene's bridgeheads match once for cbv 4
/// against the outer carbons' cbv 3, and a terminal imine N keeps its
/// submitted double — cbv 2 — because neither endpoint sits on an
/// aromatic bond). `counts` are the structure's per-atom H counts
/// (zero under DoNotAddH, the read state otherwise).
pub(crate) fn resolved_chem_valence(structure: &Structure, counts: &[u32]) -> Vec<i32> {
    let n = structure.atoms.len();
    let mut degree = vec![0i32; n];
    let mut real_pi = vec![0i32; n];
    let mut aromatic: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(from, to, order) in &structure.bonds {
        degree[from] += 1;
        degree[to] += 1;
        if order > 3 {
            aromatic[from].push(to);
            aromatic[to].push(from);
        } else {
            let unit = order as i32 - 1;
            real_pi[from] += unit;
            real_pi[to] += unit;
        }
    }
    let eligible: Vec<bool> = (0..n)
        .map(|atom| {
            let record = &structure.atoms[atom];
            let Some(el) = crate::canonical::protons::el_number(&record.element) else {
                return false;
            };
            let preferred = i32::from(crate::canonical::protons::el_valences(el, record.charge)[0]);
            if preferred == 0 {
                return false;
            }
            preferred - counts[atom] as i32 - degree[atom] - real_pi[atom] >= 1
        })
        .collect();
    let mut partner: Vec<Option<usize>> = vec![None; n];
    for atom in 0..n {
        if eligible[atom] {
            let mut seen = vec![false; n];
            kuhn_augment(atom, &aromatic, &eligible, &mut partner, &mut seen);
        }
    }
    (0..n)
        .map(|atom| degree[atom] + real_pi[atom] + i32::from(partner[atom].is_some()))
        .collect()
}

/// One Kuhn augmentation round from `atom` over the aromatic subgraph.
fn kuhn_augment(
    atom: usize,
    aromatic: &[Vec<usize>],
    eligible: &[bool],
    partner: &mut [Option<usize>],
    seen: &mut [bool],
) -> bool {
    for &neigh in &aromatic[atom] {
        if !eligible[neigh] || seen[neigh] {
            continue;
        }
        seen[neigh] = true;
        let free = match partner[neigh] {
            None => true,
            Some(prev) => kuhn_augment(prev, aromatic, eligible, partner, seen),
        };
        if free {
            partner[neigh] = Some(atom);
            partner[atom] = Some(neigh);
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    fn net(
        atoms: usize,
        bonds: &[(usize, usize, u32)],
        max_flow: &[i32],
        group_of: &[Option<usize>],
        donor: &[bool],
    ) -> Network {
        let mut bond_adjacency: Vec<Vec<(usize, i32, bool)>> = vec![Vec::new(); atoms];
        let mut single_like: Vec<Vec<bool>> = vec![Vec::new(); atoms];
        let flow = |kind: u32| match kind {
            2 => 1,
            3 => 2,
            _ => 0,
        };
        for &(a, b, kind) in bonds {
            bond_adjacency[a].push((b, flow(kind), false));
            bond_adjacency[b].push((a, flow(kind), false));
            single_like[a].push(true);
            single_like[b].push(true);
        }
        Network::build(&Build {
            bond_adjacency: &bond_adjacency,
            max_flow,
            group_of,
            donor,
            single_like: &single_like,
            c_point_of: &[None; 0],
            c_charged: &[false; 0],
        })
    }

    #[test]
    fn acetic_acid_oxygens_have_an_alternate_path() {
        // O=C-OH: the carbonyl O (acceptor, flow from the double
        // bond) and the OH oxygen (donor) can trade the H along the
        // carbon.
        let net = net(
            3,
            &[(0, 1, 2), (0, 2, 1)],
            &[1, 1, 1],
            &[None, None, None],
            &[false, false, false],
        );
        let roles = Roles {
            donor_ok: vec![false, false, true],
            acceptor_ok: vec![false, true, false],
        };
        assert!(any_alt_path(&net, &roles, 1, 2));
    }

    #[test]
    fn two_registered_groups_admit_a_two_unit_donor_path() {
        // Two carboxyl groups joined through a carbon chain, both
        // registered: the 2-H donor test between their hydroxyl
        // oxygens must find the paired augmentation (the salt-step-1
        // shape on b53: every grouped pair routes through the group
        // vertices).
        let net = net(
            7,
            &[
                (0, 1, 2),
                (0, 2, 1),
                (2, 3, 1),
                (3, 4, 2),
                (3, 5, 1),
                (5, 6, 1),
            ],
            &[1, 1, 1, 1, 1, 1, 1],
            &[Some(0), None, Some(0), None, Some(1), Some(1), None],
            &[false, false, true, false, false, true, false],
        );
        assert!(exists_2h(&mut { net }, 2, 5, true));
    }

    #[test]
    fn a_pair_without_a_donor_side_is_rejected_before_the_search() {
        // Two unregistered donors (a catechol's two OH oxygens): the
        // acceptor side of the precheck kills both directions.
        let net = net(
            3,
            &[(0, 1, 1), (0, 2, 1)],
            &[1, 1, 1],
            &[None, None, None],
            &[false, false, false],
        );
        let roles = Roles {
            donor_ok: vec![false, true, true],
            acceptor_ok: vec![false, false, false],
        };
        assert!(!any_alt_path(&net, &roles, 1, 2));
    }

    #[test]
    fn an_empty_network_never_finds_a_path() {
        // No bonds, no groups: both directions come back empty.
        let net = net(2, &[], &[0, 0], &[None, None], &[false, false]);
        let roles = Roles {
            donor_ok: vec![false, false],
            acceptor_ok: vec![false, false],
        };
        assert!(!any_alt_path(&net, &roles, 0, 1));
    }

    #[test]
    fn the_rem_proton_route_rides_the_donor_proton_to_the_acid_star() {
        // N-H-C(=O)-C-N+: the H-donor star (the amide N) and the
        // acidic =O star hang off the skeleton; the planted dots
        // travel s -> dot -> donor star -> N - C(=O) - acid star -> t,
        // exactly the fixture identity/breadth-580 shape.
        let net = net(
            5,
            &[(0, 1, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1)],
            &[0, 1, 1, 0, 0],
            &[Some(0), None, Some(1), None, None],
            &[true, false, false, false, false],
        );
        assert!(rem_proton_route_moved(&net, 0, 2).0);
    }

    #[test]
    fn a_rem_proton_route_needs_an_acidic_star_member() {
        // Same skeleton with the carbonyl oxygen unregistered: no
        // acceptor-side star, no route — the N-H stays.
        let net = net(
            5,
            &[(0, 1, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1)],
            &[0, 1, 1, 0, 0],
            &[Some(0), None, None, None, None],
            &[true, false, false, false, false],
        );
        assert!(!rem_proton_route_moved(&net, 0, 2).0);
    }

    #[test]
    fn a_flowless_donor_star_refuses_the_route_before_the_search() {
        // The reference's precheck: the dot side must carry flow — a
        // registered donor without its hydrogen plants no dots.
        let net = net(
            5,
            &[(0, 1, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1)],
            &[0, 1, 1, 0, 0],
            &[Some(0), None, Some(1), None, None],
            &[false, false, false, false, false],
        );
        assert!(!rem_proton_route_moved(&net, 0, 2).0);
    }
}
