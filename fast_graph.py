"""
fast_graph.py

Array-based (CSR) graph routines that accelerate the lazy-constraint
separation callback for the diameter-bounded compact districting MIP.

Node labels must be integers 0..n-1 (true for the JSON instances). These
routines replace the per-callback networkx calls (to_directed, diameter,
single_source_dijkstra with the 'sep-weight' trick, node_boundary,
strongly_connected_components) that dominate the original callback's runtime.

CSRGraph builds the neighbor structure ONCE (cache it on the Gurobi model and
reuse across every callback) and keeps three independent scratch areas so that
two BFS distance maps and a separation BFS can be live simultaneously without
reallocating or materializing per-node lists.

Because the inputs are PLANAR (county/tract adjacency graphs), m = O(n), so each
(capped) BFS is O(n) and every separation is near-linear.

All separation routines are validated against the original networkx logic in
test_fast_graph_all.py: every cut produced is a valid + minimal length-s
a,b-separator, and the SAME set of violating districts is flagged, so the MIP
accepts exactly the same integer solutions (identical optimum) as the reference
callbacks (compactness_callback_ref / continuous_ONLY_callback_ref).
"""


class CSRGraph:
    __slots__ = ("n", "start", "idx",
                 "_s1", "_d1", "_q1", "_g1", "_t1",
                 "_s2", "_d2", "_q2", "_g2", "_t2",
                 "_s3", "_d3", "_q3", "_g3",
                 "_member", "_blocked")

    def __init__(self, n, edges):
        deg = [0] * n
        elist = []
        for u, v in edges:
            if u == v:
                continue
            deg[u] += 1
            deg[v] += 1
            elist.append((u, v))
        start = [0] * (n + 1)
        for i in range(n):
            start[i + 1] = start[i] + deg[i]
        idx = [0] * start[n]
        pos = list(start[:n])
        for u, v in elist:
            idx[pos[u]] = v
            pos[u] += 1
            idx[pos[v]] = u
            pos[v] += 1
        self.n = n
        self.start = start
        self.idx = idx
        # three scratch areas (stamp / dist / queue / generation / tail)
        self._s1 = [0] * n; self._d1 = [0] * n; self._q1 = [0] * n; self._g1 = 0; self._t1 = 0
        self._s2 = [0] * n; self._d2 = [0] * n; self._q2 = [0] * n; self._g2 = 0; self._t2 = 0
        self._s3 = [0] * n; self._d3 = [0] * n; self._q3 = [0] * n; self._g3 = 0
        self._member = [False] * n
        self._blocked = [False] * n

    # ---- BFS into scratch area 1 (optionally capped at depth `cap`) ------- #
    def bfs1(self, src, allowed=None, cap=None):
        self._g1 += 1
        g = self._g1
        stamp = self._s1; dist = self._d1; q = self._q1
        start = self.start; idx = self.idx
        if cap is None:
            cap = 1 << 30
        stamp[src] = g; dist[src] = 0
        head = 0; q[0] = src; tail = 1
        while head < tail:
            u = q[head]; head += 1
            du = dist[u]
            if du >= cap:
                continue
            du1 = du + 1
            base = start[u]; end = start[u + 1]
            if allowed is None:
                for p in range(base, end):
                    w = idx[p]
                    if stamp[w] != g:
                        stamp[w] = g; dist[w] = du1; q[tail] = w; tail += 1
            else:
                for p in range(base, end):
                    w = idx[p]
                    if stamp[w] != g and allowed[w]:
                        stamp[w] = g; dist[w] = du1; q[tail] = w; tail += 1
        self._t1 = tail
        return g

    # ---- BFS into scratch area 2 (capped) -------------------------------- #
    def bfs2(self, src, cap=None):
        self._g2 += 1
        g = self._g2
        stamp = self._s2; dist = self._d2; q = self._q2
        start = self.start; idx = self.idx
        if cap is None:
            cap = 1 << 30
        stamp[src] = g; dist[src] = 0
        head = 0; q[0] = src; tail = 1
        while head < tail:
            u = q[head]; head += 1
            du = dist[u]
            if du >= cap:
                continue
            du1 = du + 1
            for p in range(start[u], start[u + 1]):
                w = idx[p]
                if stamp[w] != g:
                    stamp[w] = g; dist[w] = du1; q[tail] = w; tail += 1
        self._t2 = tail
        return g


def graph_diameter(csr):
    """Exact diameter of a CONNECTED graph: the largest eccentricity over all
    sources (deepest BFS level from any node). Pure array BFS -- far faster than
    networkx's nx.diameter on the large tract instances. Returns None if the graph
    is disconnected."""
    n = csr.n
    best = 0
    for src in range(n):
        csr.bfs1(src)
        if csr._t1 != n:
            return None  # disconnected: no finite diameter
        d = csr._d1; q = csr._q1
        ecc = 0
        for i in range(csr._t1):
            dw = d[q[i]]
            if dw > ecc:
                ecc = dw
        if ecc > best:
            best = ecc
    return best


def components_within(csr, nodes, member):
    """Connected components of the subgraph induced by `member` (boolean list).
    `nodes` lists the member node ids. Uses scratch area 1."""
    start = csr.start; idx = csr.idx
    csr._g1 += 1; g = csr._g1
    seen = csr._s1; q = csr._q1
    comps = []
    for s in nodes:
        if seen[s] == g:
            continue
        seen[s] = g
        comp = [s]
        head = 0; q[0] = s; tail = 1
        while head < tail:
            u = q[head]; head += 1
            for p in range(start[u], start[u + 1]):
                w = idx[p]
                if member[w] and seen[w] != g:
                    seen[w] = g; q[tail] = w; tail += 1; comp.append(w)
        comps.append(comp)
    return comps


def fischetti_separator(csr, component, member_comp, b):
    """Length-1 (contiguity) separator: nodes adjacent to `component` and
    reachable from b without entering the component's boundary. Replicates
    aux_functions.find_fischetti_separator. Uses scratch areas 1 (boundary)
    and 2 (visited)."""
    n = csr.n; start = csr.start; idx = csr.idx
    csr._g1 += 1; gb = csr._g1; boundary = csr._s1
    for u in component:
        for p in range(start[u], start[u + 1]):
            w = idx[p]
            if not member_comp[w]:
                boundary[w] = gb
    csr._g2 += 1; gv = csr._g2; visited = csr._s2; q = csr._q2
    visited[b] = gv; q[0] = b; head = 0; tail = 1
    while head < tail:
        u = q[head]; head += 1
        if boundary[u] == gb:
            continue
        for p in range(start[u], start[u + 1]):
            w = idx[p]
            if visited[w] != gv:
                visited[w] = gv; q[tail] = w; tail += 1
    return [i for i in range(n) if boundary[i] == gb and visited[i] == gv]


def separates(csr, a, b, s):
    """True iff dist_{G - blocked}(a,b) > s, where blocked = csr._blocked.
    a,b assumed not blocked. Uses scratch area 3 (so callers may keep BFS
    distance maps live in areas 1 and 2)."""
    start = csr.start; idx = csr.idx; blocked = csr._blocked
    csr._g3 += 1; g = csr._g3
    stamp = csr._s3; dist = csr._d3; q = csr._q3
    stamp[a] = g; dist[a] = 0; head = 0; q[0] = a; tail = 1
    while head < tail:
        u = q[head]; head += 1
        du = dist[u]
        if du >= s:
            continue
        du1 = du + 1
        for p in range(start[u], start[u + 1]):
            w = idx[p]
            if stamp[w] != g and not blocked[w]:
                if w == b:
                    return False
                stamp[w] = g; dist[w] = du1; q[tail] = w; tail += 1
    return True


def minimalize(csr, a, b, s, cand):
    """Greedily reduce candidate set `cand` to a MINIMAL length-s a,b-separator.
    Assumes `cand` is already a valid separator. Uses csr._blocked (left
    all-False on return)."""
    blocked = csr._blocked
    for c in cand:
        blocked[c] = True
    minC = []
    for c in cand:
        blocked[c] = False
        if not separates(csr, a, b, s):
            blocked[c] = True
            minC.append(c)
    for c in cand:
        blocked[c] = False
    return minC


def _length_s_candidates(csr, a, b, s, member_vj):
    """Candidate separator nodes for the diameter branch: non-V_j nodes v with
    dist_G(a,v)+dist_G(v,b) <= s. Capped BFS from a and b keeps this O(reach)."""
    ga = csr.bfs1(a, cap=s)            # area 1: distances from a (<= s)
    tail_a = csr._t1
    qa = csr._q1
    reached_a = qa[:tail_a]            # snapshot before area-2 bfs (separate q)
    gb = csr.bfs2(b, cap=s)            # area 2: distances from b (<= s)
    s1 = csr._s1; d1 = csr._d1; s2 = csr._s2; d2 = csr._d2
    cand = []
    for v in reached_a:
        if member_vj[v]:
            continue
        if s2[v] == gb and d1[v] + d2[v] <= s:
            cand.append(v)
    return cand


def _fischetti_for_component(csr, comp, b):
    """Fischetti separator for a single component (marks comp membership in a
    scratch boolean, then calls fischetti_separator)."""
    mc = csr._blocked  # borrow blocked as comp-membership marker (all-False in/out)
    for v in comp:
        mc[v] = True
    C = fischetti_separator(csr, comp, mc, b)
    for v in comp:
        mc[v] = False
    return C


def separate_compactness(csr, s, selector, xval, yval, minority_range,
                         majority_range, k):
    """Return the list of (a, b, minC) separators to add as lazy cuts for the
    given integer solution. Mirrors the original callback's branching:
      * disconnected district  -> Fischetti separator, minimalized to length s
      * connected, diameter > s -> length-s separator per violating pair
    Every minC is a minimal length-s a,b-separator inside V \\ V_j."""
    n = csr.n
    member = csr._member
    minority_set = set(minority_range)
    scan = list(minority_range) if selector == 'minority' else list(range(k))

    cuts = []
    for j in scan:
        if selector == 'minority':
            Vj = [v for v in range(n) if xval[v, j] > 0.5]
        elif selector == 'majority':
            Vj = [v for v in range(n) if yval[v, j] > 0.5]
        else:
            # 'both': a node in district j is assigned via x[v,j] (if j is a
            # majority-minority district) or y[v,j] (if non-MM). Read the UNION
            # so the district is identified correctly no matter which index ends
            # up being MM -- required when the z-ordering symmetry-breaking is
            # absent (e.g., bare runs). Guarded for reduced models where x/y may
            # be defined over fewer indices.
            Vj = [v for v in range(n)
                  if ((v, j) in xval and xval[v, j] > 0.5)
                  or ((v, j) in yval and yval[v, j] > 0.5)]

        if len(Vj) <= 1:
            continue

        for v in Vj:
            member[v] = True

        comps = components_within(csr, Vj, member)

        if len(comps) > 1:
            # disconnected -> contiguity cut (Fischetti, then length-s minimal)
            smallest = min(comps, key=len)
            b = smallest[0]
            for comp in comps:
                if comp is smallest:
                    continue
                a = comp[0]
                # Fischetti separator relative to THIS component:
                cand = _fischetti_for_component(csr, comp, b)
                minC = minimalize(csr, a, b, s, cand)
                cuts.append((a, b, minC))
        else:
            # connected -> diameter check via within-district BFS
            for a in Vj:
                gsub = csr.bfs1(a, allowed=member)
                s1 = csr._s1; d1 = csr._d1
                viol = [bb for bb in Vj if a < bb and s1[bb] == gsub and d1[bb] > s]
                if not viol:
                    continue
                for b in viol:
                    cand = _length_s_candidates(csr, a, b, s, member)
                    minC = minimalize(csr, a, b, s, cand)
                    cuts.append((a, b, minC))

        for v in Vj:
            member[v] = False

    return cuts


def separate_contiguity(csr, selector, xval, yval, minority_range,
                        majority_range, k):
    """Contiguity-only separation (no diameter bound): for each disconnected
    scanned district, return (a, b, C) where C is a Fischetti a,b-separator.
    Mirrors aux_functions.continuous_ONLY_callback label selection."""
    n = csr.n
    member = csr._member
    minority_set = set(minority_range)
    labels = []
    if selector != 'majority':
        labels += list(minority_range)
    if selector != 'minority':
        labels += list(majority_range)
    if selector == 'both':
        labels = list(majority_range)

    cuts = []
    for j in labels:
        if selector == 'minority':
            Vj = [v for v in range(n) if xval[v, j] > 0.5]
        elif selector == 'majority':
            Vj = [v for v in range(n) if yval[v, j] > 0.5]
        else:
            if j in minority_set:
                Vj = [v for v in range(n) if xval[v, j] > 0.5]
            else:
                Vj = [v for v in range(n) if yval[v, j] > 0.5]
        if len(Vj) <= 1:
            continue
        for v in Vj:
            member[v] = True
        comps = components_within(csr, Vj, member)
        if len(comps) > 1:
            smallest = min(comps, key=len)
            b = smallest[0]
            for comp in comps:
                if comp is smallest:
                    continue
                a = comp[0]
                C = _fischetti_for_component(csr, comp, b)
                cuts.append((a, b, C))
        for v in Vj:
            member[v] = False
    return cuts


# ---------------------------------------------------------------------- #
# Power-graph helpers (avoid the expensive nx.power / nx.complement).     #
# For PLANAR inputs m = O(n), so each capped BFS is O(n); building G^s    #
# edges this way is much faster than materializing nx.power(G, s).        #
# Validated against networkx in test_fast_graph_all.py.                   #
# ---------------------------------------------------------------------- #
def power_edges(csr, s):
    """Edges of G^s: all pairs {u, v}, u < v, with 1 <= dist_G(u,v) <= s."""
    n = csr.n
    out = []
    for u in range(n):
        csr.bfs1(u, cap=s)
        tail = csr._t1; q = csr._q1; d = csr._d1
        for i in range(tail):
            w = q[i]
            if w > u and d[w] >= 1:
                out.append((u, w))
    return out


def greedy_power_independent_set(csr, s, k):
    """A greedy maximal independent set in the power graph G^s.

    Returns (certified, vertices). If `certified` is True, `vertices` is an
    independent set in G^s with MORE than k members, which proves alpha(G^s) > k
    WITHOUT building the (possibly huge) max-independent-set MIP -- this happens
    for every small-s step of the lower-bound binary search, where alpha is large.
    Certifying ">k" costs at most k+1 capped BFS calls.

    If `certified` is False, the greedy maximal set has size <= k; this is a lower
    bound on alpha and is INCONCLUSIVE for the ">k" test, so the caller must solve
    the exact MIP. (greedy is sound, not optimal: a small greedy set does not imply
    alpha <= k.)
    """
    n = csr.n
    blocked = bytearray(n)
    chosen = []
    for u in range(n):
        if blocked[u]:
            continue
        chosen.append(u)
        # block u and all vertices within distance s of u (its closed G^s-ball)
        csr.bfs1(u, cap=s)
        q = csr._q1
        for i in range(csr._t1):
            blocked[q[i]] = 1
        if len(chosen) > k:
            return True, chosen
    return False, chosen


def far_pairs(csr, s, nodes):
    """Pairs {u, v} (u < v) among `nodes` with dist_G(u,v) > s or disconnected;
    i.e., the conflict pairs of complement(G^s) restricted to `nodes`."""
    snodes = sorted(nodes)
    out = []
    for u in snodes:
        g = csr.bfs1(u)
        s1 = csr._s1; d1 = csr._d1
        for v in snodes:
            if v <= u:
                continue
            if not (s1[v] == g and d1[v] <= s):
                out.append((u, v))
    return out


def check_feasible_plan(csr, label, k, s, pop=None, L=None, U=None,
                        minority_pop=None, vap=None, f=None, min_minority=None):
    """Efficiently verify that `label` (label[v] in 0..k-1) is a feasible
    districting plan, using only array BFS (no networkx):
      * every node assigned to a valid, non-empty district (a partition);
      * every district is CONNECTED;
      * every district has DIAMETER <= s;
      * if pop/L/U given: every district population in [L, U];
      * if minority_pop/vap/f/min_minority given: at least `min_minority`
        districts are minority-majority (minority_pop >= f * vap).
    Returns (ok: bool, reason: str). Diameter is checked exactly via a capped
    BFS from each node (capped at depth s, so each BFS is O(district size)).
    Validated against networkx in test_feasible.py.
    """
    n = csr.n
    districts = [[] for _ in range(k)]
    for v in range(n):
        j = label[v]
        if j is None or j < 0 or j >= k:
            return False, "node %d has invalid label %r" % (v, j)
        districts[j].append(v)

    member = csr._member
    minority_count = 0
    for j in range(k):
        nodes = districts[j]
        if not nodes:
            return False, "district %d is empty" % j

        if pop is not None:
            p = 0
            for v in nodes:
                p += pop[v]
            if p < L or p > U:
                return False, "district %d population %d outside [%d, %d]" % (j, p, L, U)

        for v in nodes:
            member[v] = True
        nd = len(nodes)

        # connectivity: one BFS from nodes[0] must reach all district nodes
        csr.bfs1(nodes[0], allowed=member)
        if csr._t1 != nd:
            for v in nodes:
                member[v] = False
            return False, "district %d is disconnected" % j

        # diameter <= s: a capped BFS from each node must reach all nd nodes
        bad = None
        for a in nodes:
            csr.bfs1(a, allowed=member, cap=s)
            if csr._t1 != nd:
                bad = a
                break
        if bad is not None:
            for v in nodes:
                member[v] = False
            return False, "district %d has diameter > %d (from node %d)" % (j, s, bad)

        for v in nodes:
            member[v] = False

        if minority_pop is not None and vap is not None and f is not None:
            mv = 0
            vv = 0
            for v in nodes:
                mv += minority_pop[v]
                vv += vap[v]
            if mv >= f * vv:
                minority_count += 1

    if min_minority is not None and minority_count < min_minority:
        return False, "only %d minority-majority districts (< %d)" % (minority_count, min_minority)

    return True, "feasible"


# ---------------------------------------------------------------------- #
# Combinatorial fixing pre-screens (NO MIP).                              #
# A district D with diam(D) <= s containing v satisfies D subseteq        #
# B_s(v) (radius-s ball of v in the induced graph). Hence two SOUND       #
# necessary conditions for v to lie in a feasible minority-majority       #
# district (pop in [L,U], minority_pop >= f*VAP):                         #
#   (1) pop(B_s(v)) >= L                                                   #
#   (2) max_{S in B_s(v), v in S, pop(S)>=L} sum(minority_i - f*VAP_i) >=0 #
# (2) is evaluated by its LP relaxation (drops connectivity, the diameter, #
# and the upper bound U), so it is a sound over-estimate: value < 0 =>     #
# fix v. Validated against brute force in test_prefilter.py (0 false      #
# fixings). Use this to skip the per-vertex feasibility MIP.              #
# ---------------------------------------------------------------------- #
def _ball_nodes(csr, v, s, allowed):
    csr.bfs1(v, allowed=allowed, cap=s)
    return list(csr._q1[:csr._t1])


def _minority_surplus_ub(nodes, v, pop, minority, vap, f, L):
    surplus = {u: minority[u] - f * vap[u] for u in nodes}
    S_surp = surplus[v]
    S_pop = pop[v]
    for u in nodes:
        if u != v and surplus[u] > 0:
            S_surp += surplus[u]
            S_pop += pop[u]
    if S_pop >= L:
        return S_surp, True
    deficit = L - S_pop
    rest = [u for u in nodes if u != v and surplus[u] <= 0 and pop[u] > 0]
    rest.sort(key=lambda u: surplus[u] / pop[u], reverse=True)
    for u in rest:
        if deficit <= 0:
            break
        take = min(1.0, deficit / pop[u])
        S_surp += take * surplus[u]
        deficit -= take * pop[u]
    if deficit > 1e-9:
        return S_surp, False
    return S_surp, True


def _minority_surplus_ub_under_U(nodes, v, pop, minority, vap, f, U):
    """Upper bound on the minority surplus sum(minority - f*vap) achievable by a
    set T with v in T and pop(T) <= U. Negative-surplus nodes never help, so the
    bound is surplus[v] plus a fractional (LP) knapsack of the positive-surplus
    nodes under the residual population budget U - pop[v]. Returns None if v alone
    already exceeds U (then no district of population <= U can contain v)."""
    budget = U - pop[v]
    if budget < -1e-9:
        return None
    total = minority[v] - f * vap[v]
    pos = []
    for u in nodes:
        if u == v:
            continue
        su = minority[u] - f * vap[u]
        if su > 0:
            if pop[u] <= 0:
                total += su          # zero-population positive node: free surplus
            else:
                pos.append((su, pop[u]))
    # fill the budget with the highest surplus-per-population nodes first
    pos.sort(key=lambda sp: sp[0] / sp[1], reverse=True)
    for su, pp in pos:
        if budget <= 1e-9:
            break
        take = min(1.0, budget / pp)
        total += take * su
        budget -= take * pp
    return total


def can_fix_minority(csr, v, s, allowed, pop, L, U, minority, vap, f):
    """Sound: returns True only if v provably cannot belong to any feasible
    minority-majority district with diameter <= s (so X[v, *] may be fixed to 0
    without solving the per-vertex MIP). Three cheap, sound certificates, each a
    relaxation of the per-vertex feasibility MIP:
      (1) the radius-s ball cannot even supply the minimum population L;
      (2) max surplus subject to pop >= L is negative (drops the U cap);
      (3) max surplus subject to pop <= U is negative, or v alone exceeds U
          (drops the L floor).
    If any relaxation is infeasible, the integer MIP is too."""
    nodes = _ball_nodes(csr, v, s, allowed)
    if sum(pop[u] for u in nodes) < L:
        return True
    ub, pop_ok = _minority_surplus_ub(nodes, v, pop, minority, vap, f, L)
    if not pop_ok:
        return True
    if ub < -1e-9:
        return True
    ub_U = _minority_surplus_ub_under_U(nodes, v, pop, minority, vap, f, U)
    if ub_U is None or ub_U < -1e-9:
        return True
    return False
