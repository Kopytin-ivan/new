# geom/outer.py
from typing import List, Tuple, Dict, Iterable
import math

from geom.graph import Graph
from geom.cycles import find_planar_faces
from geom.io import save_segments

# Если у тебя уже есть эта функция из прошлого шага — оставь её и импортируй здесь:
from geom.export_dxf import save_dxf_polyline  # пишет ОДНУ закрытую LWPOLYLINE в DXF (Meters)
from geom.planarize import planarize_graph

from collections import defaultdict


Point   = Tuple[float, float]
Segment = Tuple[Point, Point]

# ---------- вспомогалки ----------

def _poly_area(node_ids: List[int], G: Graph) -> float:
    """Подписанная площадь по формуле Шоена."""
    s = 0.0
    k = len(node_ids)
    for i in range(k):
        x1, y1 = G.nodes[node_ids[i]]
        x2, y2 = G.nodes[node_ids[(i + 1) % k]]
        s += x1 * y2 - x2 * y1
    return 0.5 * s

def _cycle_to_segments(node_ids: List[int], G: Graph) -> List[Segment]:
    out: List[Segment] = []
    k = len(node_ids)
    for i in range(k):
        a = node_ids[i]
        b = node_ids[(i + 1) % k]
        out.append((G.nodes[a], G.nodes[b]))
    return out

def _cycle_to_points(node_ids: List[int], G: Graph) -> List[Point]:
    return [G.nodes[nid] for nid in node_ids]

def prune_leaves(G: Graph, min_len: float = 0.0) -> Graph:
    """
    Итеративно срезает листья (degree=1) и, опц., короткие рёбра < min_len.
    Работает на копии.
    """
    H = G.clone()
    while True:
        # степень для каждой вершины
        deg = {nid: H.degree(nid) for nid in range(len(H.nodes))}
        to_remove = set()
        for ei, (a, b) in enumerate(H.edges):
            if a == -1 or b == -1:
                continue
            if min_len > 0.0:
                xa, ya = H.nodes[a]; xb, yb = H.nodes[b]
                if math.hypot(xa - xb, ya - yb) < min_len:
                    to_remove.add(ei); continue
            if deg.get(a, 0) == 1 or deg.get(b, 0) == 1:
                to_remove.add(ei)
        if not to_remove:
            break
        for ei in to_remove:
            H.remove_edge(ei)
    return H

# ---------- основной пайплайн «внешки» ----------

def extract_outer_cycle_nodes(G: Graph,
                              drop_leaves: bool = True,
                              leaf_len_mm: float = 0.0,
                              eps_snap_m: float = 0.002) -> Tuple[List[int], Graph]:
    """
    Планаризует G -> H, (опц.) срезает листья в H и возвращает:
      (список id узлов цикла внешней грани в H, сам граф H)
    """
    # 1) планаризация (режем пересечения и T-стыки)
    H = planarize_graph(G, eps=eps_snap_m)

    # 2) (опц.) срез листьев/шипов уже в H
    if drop_leaves:
        H = prune_leaves(H, min_len=(leaf_len_mm/1000.0) if leaf_len_mm > 0 else 0.0)

    # 3) поиск граней на H
    faces = find_planar_faces(H, include_outer=True, right_hand=False)
    if not faces:
        return [], H

    # 4) внешняя — с максимальной |площадью| на H
    outer_nodes = max(faces, key=lambda cyc: abs(_poly_area(cyc, H)))
    return outer_nodes, H

def save_outer_from_graph(G: Graph,
                          out_json_path: str,
                          out_dxf_path: str,
                          drop_leaves: bool = True,
                          leaf_len_mm: float = 0.0,
                          eps_snap_m: float = 0.002) -> Dict:
    """
    Сохраняет внешний контур:
    - JSON: список сегментов
    - DXF: одна закрытая LWPOLYLINE (Meters)
    Возвращает метаданные.
    """
    cyc, H = extract_outer_cycle_nodes(G, drop_leaves=drop_leaves,
                                       leaf_len_mm=leaf_len_mm,
                                       eps_snap_m=eps_snap_m)
    if not cyc:
        save_segments(out_json_path, [], params={"meta": {"faces": 0}})
        return {"faces": 0}

    # Сегменты и точки берём из H (не из G!)
    outline = _cycle_to_segments(cyc, H)
    pts     = _cycle_to_points(cyc, H)

    area = _poly_area(cyc, H)
    meta = {
        "nodes": len(H.nodes),
        "edges": sum(1 for u, v in H.edges if u != -1 and v != -1),
        "faces": 1,
        "signed_area": area,
        "orientation": "CCW" if area > 0 else "CW"
    }

    save_segments(out_json_path, outline, params={"meta": meta, "source": "graph_planarized"})
    save_dxf_polyline(pts, out_dxf_path, layer="OUTER", color=7, lineweight=25,
                      insunits="Meters", closed=True)
    return meta


# ---------- утилиты площади/габарита ----------
def _bbox(nodes):
    xs = [x for x, _ in nodes]
    ys = [y for _, y in nodes]
    return (min(xs), min(ys), max(xs), max(ys))

def _area_abs(nodes, cyc):
    s = 0.0
    for i in range(len(cyc)):
        x1, y1 = nodes[cyc[i]]
        x2, y2 = nodes[cyc[(i+1) % len(cyc)]]
        s += x1*y2 - x2*y1
    return abs(0.5*s)

# ---------- fallback: прямой «правый» обход внешки ----------
def _angle(p, q):
    return math.atan2(q[1]-p[1], q[0]-p[0])

def _build_nbrs_with_eids(G: Graph):
    """Возвращает: nbrs[u] = [(v, eid)] в CCW-порядке; pos[u][(v,eid)] -> idx."""
    nbrs, pos = {}, {}
    for u in range(len(G.nodes)):
        eids = G.adj.get(u, [])
        if not eids:
            continue
        Pu = G.nodes[u]
        L = []
        for eid in eids:
            a, b = G.edge_nodes(eid)
            if a == -1 or b == -1:
                continue
            v = b if a == u else a
            L.append((_angle(Pu, G.nodes[v]), v, eid))
        if not L:
            continue
        L.sort()  # CCW
        nbrs[u] = [(v, eid) for _, v, eid in L]
        pos[u]  = {(v, eid): i for i, (v, eid) in enumerate(nbrs[u])}
    return nbrs, pos

def _next_right_halfedge(nbrs, pos, u, v, eid, right_hand=True):
    L = nbrs.get(v)
    if not L: 
        return None
    i = pos[v].get((u, eid))
    if i is None: 
        return None
    j = (i - 1) % len(L) if right_hand else (i + 1) % len(L)
    v2, eid2 = L[j]
    return (v, v2, eid2)

def trace_outer_by_giftwrap(G: Graph, right_hand=False, max_steps_factor=20):
    """
    Возвращает один цикл внешней границы правым обходом.
    Работает даже на разреженных графах, где face-обход даёт «обрубки».
    """
    nbrs, pos = _build_nbrs_with_eids(G)
    if not nbrs:
        return []

    # старт — вершина с минимальным y (при равенстве — с минимальным x)
    start = min(range(len(G.nodes)), key=lambda i: (G.nodes[i][1], G.nodes[i][0]))
    if start not in nbrs:
        # возьмём ближайшую с рёбрами
        candidates = [u for u in nbrs.keys()]
        start = min(candidates, key=lambda i: (G.nodes[i][1], G.nodes[i][0]))

    # из стартовой берём полуребро с минимальным углом к +X (самое "вправо")
    Pu = G.nodes[start]
    best = None
    best_ang = +1e9
    for v, eid in nbrs[start]:
        ang = _angle(Pu, G.nodes[v])  # [-pi, pi]
        if ang < best_ang:
            best_ang = ang
            best = (start, v, eid)
    if best is None:
        return []

    u0, v0, e0 = best
    u, v, eid = u0, v0, e0
    cyc = [u0]
    visited_steps = 0
    max_steps = max_steps_factor * max(1, sum(1 for (a,b) in G.edges if a!=-1 and b!=-1))

    while visited_steps < max_steps:
        visited_steps += 1
        cyc.append(v)
        nxt = _next_right_halfedge(nbrs, pos, u, v, eid, right_hand=right_hand)
        if nxt is None:
            break
        u, v, eid = nxt
        if (u, v, eid) == (u0, v0, e0):
            break

    # замкнутся должны в ту же полурёберу
    if len(cyc) >= 3 and cyc[0] == cyc[-1]:
        cyc = cyc[:-1]
    return cyc if len(cyc) >= 3 else []

# ---------- «лучший из двух»: faces → иначе giftwrap ----------
def extract_outer_cycle_best(G: Graph, eps_snap_m: float = 0.002):
    """
    1) Планаризация → faces (включая outer) → берём цикл с макс |площади|.
    2) Если цикл аномально мал по отношению к bbox сцены — fallback: giftwrap.
    Возвращает (cyc_nodes, Graph_used_for_coords)
    """
    H = planarize_graph(G, eps=eps_snap_m)

    faces = find_planar_faces(H, include_outer=True, right_hand=False)
    if faces:
        outer = max(faces, key=lambda cyc: _area_abs(H.nodes, cyc))
        # sanity-check: цикл должен покрывать значимую долю bbox
        bx0, by0, bx1, by1 = _bbox(H.nodes)
        bbox_area = max((bx1 - bx0) * (by1 - by0), 1e-12)
        outer_area = _area_abs(H.nodes, outer)
        if outer_area >= 0.2 * bbox_area:
            return outer, H

    # fallback: прямой правый обход на исходном графе (без планаризации)
    cyc2 = trace_outer_by_giftwrap(G, right_hand=False)
    if cyc2:
        return cyc2, G

    # не получилось — вернём пусто (пусть вызывающий обработает)
    return [], H

def save_outer_best(G: Graph, out_json_path: str, out_dxf_path: str,
                    eps_snap_m: float = 0.002) -> Dict:
    cyc, X = extract_outer_cycle_best(G, eps_snap_m=eps_snap_m)
    if not cyc:
        save_segments(out_json_path, [], params={"meta": {"faces": 0}})
        return {"faces": 0}

    outline = _cycle_to_segments(cyc, X)
    pts     = _cycle_to_points(cyc, X)

    area = _poly_area(cyc, X)
    meta = {
        "nodes": len(X.nodes),
        "edges": sum(1 for u, v in X.edges if u != -1 and v != -1),
        "faces": 1,
        "signed_area": area,
        "orientation": "CCW" if area > 0 else "CW",
    }

    save_segments(out_json_path, outline, params={"meta": meta, "source": "outer_best"})
    save_dxf_polyline(pts, out_dxf_path, layer="OUTER", color=7, lineweight=25,
                      insunits="Meters", closed=True)
    return meta


def _build_ccw_nbrs_from_edges(nodes, undirected_edges):
    """
    undirected_edges: iterable[(u,v)] без дублей (u!=v), неориентированные.
    Возвращает:
      nbrs[u] = [v1, v2, ...]  (CCW по углу вокруг u)
      pos[(u,v)] = index в nbrs[u]
    """
    # собираем неориентированный набор
    st = set()
    for (u, v) in undirected_edges:
        if u == v:
            continue
        a, b = (u, v) if u < v else (v, u)
        st.add((a, b))

    # делаем ориентированный список (оба направления)
    dir_adj = defaultdict(list)
    import math as _m
    for (a, b) in st:
        Pa = nodes[a]; Pb = nodes[b]
        ang_ab = _m.atan2(Pb[1]-Pa[1], Pb[0]-Pa[0])
        ang_ba = _m.atan2(Pa[1]-Pb[1], Pa[0]-Pb[0])
        dir_adj[a].append((ang_ab, b))
        dir_adj[b].append((ang_ba, a))

    # сортировка CCW и позиции
    nbrs = {}
    pos  = {}
    for u, L in dir_adj.items():
        L.sort()                     # по углу CCW
        nbrs[u] = [v for _, v in L]
        pos.update({(u, v): i for i, v in enumerate(nbrs[u])})
    return nbrs, pos

def _walk_cycles_on_edges(nodes, boundary_edges):
    """
    boundary_edges — неориентированные ребра, образующие (возможно несколько) циклов.
    Обходит все циклы по правилу "самый правый после твина".
    Возвращает список циклов (списки node_id), ориентированные CCW.
    """
    nbrs, pos = _build_ccw_nbrs_from_edges(nodes, boundary_edges)
    visited = set()  # ориентированные полурёбра (u,v)

    def next_right(u, v):
        if v not in nbrs:
            return None
        i = pos.get((v, u))
        if i is None:
            return None
        j = (i + 1) % len(nbrs[v])
        return (v, nbrs[v][j])

    cycles = []
    for u in list(nbrs.keys()):
        for v in nbrs[u]:
            if (u, v) in visited:
                continue
            # старт обхода
            u0, v0 = u, v
            cyc = [u0]
            while True:
                visited.add((u, v))
                cyc.append(v)
                nxt = next_right(u, v)
                if nxt is None:
                    cyc = []  # несвязность — цикл не замкнулся
                    break
                u, v = nxt
                if (u, v) == (u0, v0):
                    break
                if (u, v) in visited:
                    cyc = []  # зашли в уже посещённое полуребро — не цикл
                    break
            if len(cyc) >= 3:
                if cyc[0] == cyc[-1]:
                    cyc = cyc[:-1]
                # нормализуем CCW
                if _area_abs(nodes, cyc) < 0:
                    cyc = list(reversed(cyc))
                cycles.append(cyc)
    return cycles

def extract_outer_via_faces_union(G: Graph, eps_snap_m: float = 0.002):
    """
    1) Планаризуем граф → H
    2) Находим ТОЛЬКО ограниченные лица (include_outer=False)
    3) Считаем для каждой неориентированной кромки (u,v), сколько раз она встречается в границах лиц.
       Кромки с count==1 образуют границу объединения внутренних лиц.
    4) Собираем из них циклы и берём максимальный по |площади| — внешний контур.
    Возвращает (cyc_nodes, Graph H, все_циклы).
    """
    H = planarize_graph(G, eps=eps_snap_m)

    faces = find_planar_faces(H, include_outer=False, right_hand=False)  # только внутренние
    if not faces:
        return [], H, []

    # считаем использование кромок (u,v) независимо от ориентации
    edge_use = defaultdict(int)
    for cyc in faces:
        for i in range(len(cyc)):
            u = cyc[i]
            v = cyc[(i+1) % len(cyc)]
            a, b = (u, v) if u < v else (v, u)
            edge_use[(a, b)] += 1

    # кромки на границе объединения — те, где счётчик == 1
    boundary_edges = [(u, v) for (u, v), c in edge_use.items() if c == 1]
    if not boundary_edges:
        return [], H, []

    # из этих ребер строим циклы
    cycles = _walk_cycles_on_edges(H.nodes, boundary_edges)
    if not cycles:
        return [], H, []

    # выбираем цикл с максимальной площадью — внешний
    best = max(cycles, key=lambda cyc: _area_abs(H.nodes, cyc))
    return best, H, cycles

def save_outer_via_faces_union(G: Graph, out_json_path: str, out_dxf_path: str,
                               eps_snap_m: float = 0.002) -> Dict:
    cyc, H, _ = extract_outer_via_faces_union(G, eps_snap_m=eps_snap_m)
    if not cyc:
        save_segments(out_json_path, [], params={"meta": {"faces": 0}})
        return {"faces": 0}

    outline = _cycle_to_segments(cyc, H)
    pts     = _cycle_to_points(cyc, H)

    area = _poly_area(cyc, H)
    meta = {
        "nodes": len(H.nodes),
        "edges": sum(1 for u, v in H.edges if u != -1 and v != -1),
        "faces": 1,
        "signed_area": area,
        "orientation": "CCW" if area > 0 else "CW",
        "method": "faces_union"
    }

    save_segments(out_json_path, outline, params={"meta": meta, "source": "faces_union"})
    save_dxf_polyline(pts, out_dxf_path, layer="OUTER", color=7, lineweight=25,
                      insunits="Meters", closed=True)
    return meta