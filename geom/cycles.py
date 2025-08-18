# geom/cycles.py
from collections import deque

def connected_components(G):
    seen = set()
    comps = []
    for start in range(len(G.nodes)):
        if start in seen:
            continue
        if start not in G.adj:
            seen.add(start)
            comps.append([start])
            continue
        q = deque([start])
        seen.add(start)
        comp = [start]
        while q:
            u = q.popleft()
            for eid in G.adj.get(u, []):
                a, b = G.edge_nodes(eid)
                v = b if a == u else a
                if v not in seen:
                    seen.add(v)
                    q.append(v)
                    comp.append(v)
        comps.append(comp)
    return comps

def is_degree2_component(G, comp_nodes):
    for nid in comp_nodes:
        if G.degree(nid) != 2:
            return False
    return True

def _component_edges(G, comp_nodes):
    comp = set(comp_nodes)
    eids = []
    for nid in comp_nodes:
        for eid in G.adj.get(nid, []):
            a, b = G.edge_nodes(eid)
            if a in comp and b in comp:
                eids.append(eid)
    return eids

def _next_neighbor(G, prev_node, curr_node):
    """В degree=2 у узла ровно два соседа; берём тот, который != prev."""
    out = []
    for eid in G.adj[curr_node]:
        a, b = G.edge_nodes(eid)
        v = b if a == curr_node else a
        out.append((eid, v))
    # out = [(eid1, v1), (eid2, v2)]
    if not out:
        return None
    if len(out) == 1:
        return out[0]  # деградация, но обычно не бывает в degree=2
    (e1, v1), (e2, v2) = out
    if v1 != prev_node:
        return (e1, v1)
    else:
        return (e2, v2)

# geom/cycles.py
from collections import defaultdict, deque

def find_simple_cycles(G):
    """
    Возвращает список простых циклов как списки node_id в порядке обхода.
    Работает на общем неориентированном графе (узлы могут иметь любую степень).
    Алгоритм: строим остов (DFS/BFS) по компонентам; для каждого нетростового ребра (u,v)
    восстанавливаем путь u→v по остову и получаем фундаментальный цикл.
    Дедупликуем циклы канонизацией (минимальная ротация + выбор направления).
    """
    n = len(G.nodes)
    adj = defaultdict(list)
    for eid, (u, v) in enumerate(G.edges):
        adj[u].append(v)
        adj[v].append(u)

    seen_cycles = set()
    out_cycles = []

    seen_node = set()
    for root in range(n):
        if root in seen_node or root not in G.adj:
            continue

        # BFS-остов для компоненты
        parent = {root: None}
        depth  = {root: 0}
        order  = []
        q = deque([root])
        seen_node.add(root)

        while q:
            u = q.popleft()
            order.append(u)
            for v in adj[u]:
                if v not in parent:
                    parent[v] = u
                    depth[v]  = depth[u] + 1
                    seen_node.add(v)
                    q.append(v)

        # множество остовных рёбер (для быстрого теста)
        tree_edges = set()
        for v in parent:
            u = parent[v]
            if u is not None:
                a, b = (u, v) if u < v else (v, u)
                tree_edges.add((a, b))

        # хелпер: путь между двумя вершинами по остову
        def path_between(a, b):
            pa, pb = [a], [b]
            ua, ub = a, b
            # выровнять глубины
            while depth[ua] > depth[ub]:
                ua = parent[ua]; pa.append(ua)
            while depth[ub] > depth[ua]:
                ub = parent[ub]; pb.append(ub)
            # подниматься до LCA
            while ua != ub:
                ua = parent[ua]; pa.append(ua)
                ub = parent[ub]; pb.append(ub)
            lca = ua
            # цикл: a..LCA + reverse(b..LCA без повторения LCA)
            return pa + list(reversed(pb[:-1]))

        # обойти нетростовые рёбра и собрать фундаментальные циклы
        visited_back = set()
        for u in parent.keys():
            for v in adj[u]:
                a, b = (u, v) if u < v else (v, u)
                if (a, b) in tree_edges:
                    continue  # остовное ребро
                if (a, b) in visited_back:
                    continue
                visited_back.add((a, b))

                # соберём цикл
                cyc = path_between(u, v)
                if len(cyc) < 3:
                    continue

                # канонизация для дедупликации
                # — сдвиг до минимального id
                k = len(cyc)
                min_pos = min(range(k), key=lambda i: cyc[i])
                rot1 = cyc[min_pos:] + cyc[:min_pos]
                # — и обратное направление
                rcyc = list(reversed(cyc))
                min_pos_r = min(range(k), key=lambda i: rcyc[i])
                rot2 = rcyc[min_pos_r:] + rcyc[:min_pos_r]
                canon = tuple(rot2) if tuple(rot2) < tuple(rot1) else tuple(rot1)

                if canon not in seen_cycles:
                    seen_cycles.add(canon)
                    out_cycles.append(list(canon))

    return out_cycles

# --- ВСТАВЬ НИЖЕ В КОНЕЦ geom/cycles.py ---
import math
from collections import defaultdict

def _angle(p, q):
    return math.atan2(q[1]-p[1], q[0]-p[0])

def find_planar_faces(G, include_outer=False, right_hand=True):
    """
    Возвращает список границ граней (каждая — список node_id по кругу).
    Алгоритм: half-edge обход. У каждого узла соседи сортируются по углу.
    include_outer=False — внешнюю грань удаляем (по максимальной |площади|).
    right_hand=True — правило правой руки; False — левой.
    """
    # 1) угловой порядок соседей
    nbrs = {}
    for u in range(len(G.nodes)):
        if u not in G.adj:
            continue
        P = G.nodes[u]
        seen = set()
        neigh = []
        for eid in G.adj[u]:
            a, b = G.edge_nodes(eid)
            v = b if a == u else a
            if v in seen:
                continue
            seen.add(v)
            ang = _angle(P, G.nodes[v])
            neigh.append((ang, v))
        neigh.sort()  # CCW
        nbrs[u] = [v for _, v in neigh]

    # 2) полурёбра
    visited = set()  # направленные ребра (u,v)

    def next_right(u, v):
        """Из u->v поворачиваем 'вправо' вокруг v (берём соседа перед u в CCW-списке v)."""
        L = nbrs.get(v)
        if not L: return None
        i = L.index(u)
        return L[i-1] if right_hand else L[(i+1) % len(L)]

    def face_walk(u0, v0):
        u, v = u0, v0
        cyc = [u]
        while True:
            visited.add((u, v))
            cyc.append(v)
            w = next_right(u, v)
            if w is None:
                return None
            u, v = v, w
            if (u, v) == (u0, v0):
                break
            if (u, v) in visited:
                break
        return cyc[:-1] if cyc and cyc[0] == cyc[-1] else cyc

    faces = []
    for u in nbrs:
        for v in nbrs[u]:
            if (u, v) in visited:
                continue
            cyc = face_walk(u, v)
            if not cyc or len(cyc) < 3:
                continue
            faces.append(cyc)

    # 3) площадь и фильтр внешней
    def area(cyc):
        pts = [G.nodes[i] for i in cyc]
        s = 0.0
        for i in range(len(pts)):
            x1, y1 = pts[i]
            x2, y2 = pts[(i+1) % len(pts)]
            s += x1*y2 - x2*y1
        return 0.5*s

    if not faces:
        return []

    # внешняя грань — с максимальной |площади|
    if not include_outer:
        outer_idx = max(range(len(faces)), key=lambda i: abs(area(faces[i])))
        faces.pop(outer_idx)

    # нормализуем направление (CCW)
    out = []
    for cyc in faces:
        if area(cyc) < 0:
            cyc = list(reversed(cyc))
        out.append(cyc)
    return out
