from typing import List, Tuple, Dict, Set
import math
from geom.graph import Graph
from geom.grid import UniformGrid

Point = Tuple[float, float]

def _dot(ax, ay, bx, by): return ax*bx + ay*by
def _cross(ax, ay, bx, by): return ax*by - ay*bx

def _orient(a: Point, b: Point, c: Point) -> float:
    return _cross(b[0]-a[0], b[1]-a[1], c[0]-a[0], c[1]-a[1])

def _between(a: Point, b: Point, c: Point, eps: float) -> bool:
    # c лежит на отрезке ab (с допуском)
    if abs(_orient(a,b,c)) > eps: 
        return False
    minx, maxx = (a[0], b[0]) if a[0] <= b[0] else (b[0], a[0])
    miny, maxy = (a[1], b[1]) if a[1] <= b[1] else (b[1], a[1])
    return (minx - eps <= c[0] <= maxx + eps) and (miny - eps <= c[1] <= maxy + eps)

def _seg_intersection(a: Point, b: Point, c: Point, d: Point, eps: float):
    """
    ('proper', P)  — пересечение в интерьерах обоих
    ('t_on_cd', P) — a или b лежит на cd
    ('u_on_ab', P) — c или d лежит на ab
    ('none', None) — иначе
    Коллинеарные перекрытия здесь не обрабатываем (их лучше схлопнуть на подготовке).
    """
    # bbox-тест
    if (max(a[0], b[0]) + eps < min(c[0], d[0]) or
        max(c[0], d[0]) + eps < min(a[0], b[0]) or
        max(a[1], b[1]) + eps < min(c[1], d[1]) or
        max(c[1], d[1]) + eps < min(a[1], b[1])):
        return ('none', None)

    o1 = _orient(a, b, c)
    o2 = _orient(a, b, d)
    o3 = _orient(c, d, a)
    o4 = _orient(c, d, b)

    # proper cross — строгие разнонаправленные знаки
    if (o1 * o2 < 0) and (o3 * o4 < 0):
        ax, ay = a; bx, by = b; cx, cy = c; dx, dy = d
        r_x, r_y = bx-ax, by-ay
        s_x, s_y = dx-cx, dy-cy
        denom = _cross(r_x, r_y, s_x, s_y)
        if abs(denom) < 1e-18:
            return ('none', None)
        t = _cross(cx-ax, cy-ay, s_x, s_y) / denom
        P = (ax + t*r_x, ay + t*r_y)
        return ('proper', P)

    # T-стыки (конец одного на интерьере другого)
    for p in (a, b):
        if _between(c, d, p, eps):
            return ('t_on_cd', p)
    for p in (c, d):
        if _between(a, b, p, eps):
            return ('u_on_ab', p)

    return ('none', None)

def planarize_graph(G: Graph, eps: float = 0.002, grid_cell: float = None, max_passes: int = 4) -> Graph:
    """
    Итеративная планаризация: режем пересечения и T-стыки, повторяем до стационарности.
    eps — в метрах (0.002 = 2 мм). Возвращает новый граф.
    """
    H = G.clone()

    for _pass in range(max_passes):
        # 1) грид
        if grid_cell is None:
            grid_cell_use = max(eps * 50.0, 0.02)  # ≥ 2 см
        else:
            grid_cell_use = grid_cell
        grid = UniformGrid(grid_cell_use)

        alive = []
        for ei, (u, v) in enumerate(H.edges):
            if u == -1 or v == -1:
                continue
            a, b = H.nodes[u], H.nodes[v]
            grid.insert_segment(ei, a, b, pad=0.0)
            alive.append(ei)
        alive_set: Set[int] = set(alive)

        # 2) поиск кандидатов и накопление точек разреза
        cut_points: Dict[int, List[Point]] = {}
        seen_pairs: Set[Tuple[int, int]] = set()

        def acc(ei: int, P: Point):
            cut_points.setdefault(ei, []).append(P)

        changed = False
        for _, eids in grid.cells.items():
            m = len(eids)
            if m < 2:
                continue
            for i in range(m):
                ei = eids[i]
                if ei not in alive_set: 
                    continue
                u1, v1 = H.edges[ei]
                if u1 == -1 or v1 == -1:
                    continue
                a = H.nodes[u1]; b = H.nodes[v1]

                for j in range(i+1, m):
                    ej = eids[j]
                    if ej not in alive_set or ei == ej:
                        continue
                    pair = (min(ei, ej), max(ei, ej))
                    if pair in seen_pairs:
                        continue
                    seen_pairs.add(pair)

                    u2, v2 = H.edges[ej]
                    if u2 == -1 or v2 == -1:
                        continue
                    c = H.nodes[u2]; d = H.nodes[v2]

                    # общая вершина — уже планарно
                    if len({u1, v1, u2, v2}) < 4:
                        continue

                    kind, P = _seg_intersection(a, b, c, d, eps)
                    if kind == 'proper':
                        acc(ei, P); acc(ej, P); changed = True
                    elif kind == 't_on_cd':
                        acc(ej, P); changed = True
                    elif kind == 'u_on_ab':
                        acc(ei, P); changed = True

        # 3) если нечего резать — выходим
        if not changed:
            break

        # 4) режем по отсортированным точкам вдоль каждого ребра
        for ei, pts in list(cut_points.items()):
            if ei >= len(H.edges):
                continue
            u, v = H.edges[ei]
            if u == -1 or v == -1:
                continue

            ax, ay = H.nodes[u]; bx, by = H.nodes[v]
            vx, vy = (bx - ax, by - ay)
            L2 = vx*vx + vy*vy
            if L2 <= 1e-18:
                continue

            def _t(p: Point) -> float:
                return _dot(p[0]-ax, p[1]-ay, vx, vy) / L2

            pts_on = [p for p in pts if _between((ax, ay), (bx, by), p, eps)]
            if not pts_on:
                continue
            pts_sorted = sorted(pts_on, key=_t)

            current_eid = ei
            for P in pts_sorted:
                ucur, vcur = H.edges[current_eid]
                if ucur == -1 or vcur == -1:
                    break
                # если близко к концам — пропускаем
                if (math.hypot(H.nodes[ucur][0] - P[0], H.nodes[ucur][1] - P[1]) <= eps or
                    math.hypot(H.nodes[vcur][0] - P[0], H.nodes[vcur][1] - P[1]) <= eps):
                    continue
                _, e1, e2 = H.split_edge(current_eid, P)
                current_eid = e2  # продолжаем резку правой части

    return H
