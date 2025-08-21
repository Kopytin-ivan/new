# geom/outer.py
from typing import List, Tuple, Dict, Iterable
import math

from geom.graph import Graph
from geom.cycles import find_planar_faces
from geom.io import save_segments

# Если у тебя уже есть эта функция из прошлого шага — оставь её и импортируй здесь:
from geom.export_dxf import save_dxf_polyline  # пишет ОДНУ закрытую LWPOLYLINE в DXF (Meters)
from geom.planarize import planarize_graph



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
    faces = find_planar_faces(H, include_outer=True, right_hand=True)
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
