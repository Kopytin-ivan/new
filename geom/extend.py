# geom/extend.py
import math
from geom.templates import make_template
from geom.cycles import find_planar_faces

def _len2(ax, ay): return ax*ax + ay*ay
def _len(ax, ay): return math.hypot(ax, ay)
def _dot(ax, ay, bx, by): return ax*bx + ay*by
def _cross(ax, ay, bx, by): return ax*by - ay*bx

def _unit(dx, dy):
    L = _len(dx, dy)
    return (0.0, 0.0) if L <= 1e-18 else (dx/L, dy/L)

def _ray_segment_intersection(P, D, A, B, ang_tol_deg=2.0, col_perp_tol=10.0, eps=1e-12):
    """
    Пересечение луча P + t*D (t>=0) с отрезком [A,B].
    Возвращает (t, u, X) или None.
    Коллинеарный случай: если [A,B] почти на оси луча (поперечное отклонение <= col_perp_tol),
    берём ближайшую вперёд точку (A/B).
    """
    rx, ry = D
    sx, sy = (B[0]-A[0], B[1]-A[1])
    denom = _cross(rx, ry, sx, sy)

    if abs(denom) > eps:
        APx, APy = (A[0]-P[0], A[1]-P[1])
        t = _cross(APx, APy, sx, sy) / denom
        u = _cross(APx, APy, rx, ry) / denom
        if t < 0.0 or u < 0.0 or u > 1.0:
            return None
        X = (P[0] + t*rx, P[1] + t*ry)
        return t, u, X

    # почти параллельно → смотрим поперечное отклонение концов сегмента от оси луча
    def _perp(Z):
        vx, vy = (Z[0]-P[0], Z[1]-P[1])
        s = _dot(vx, vy, rx, ry)
        per2 = max(0.0, vx*vx + vy*vy - s*s)
        return s, math.sqrt(per2)

    sA, perA = _perp(A)
    sB, perB = _perp(B)
    if perA <= col_perp_tol and perB <= col_perp_tol:
        cand = [s for s in (sA, sB) if s >= 0.0]
        if not cand:
            return None
        t = min(cand)
        X = (P[0] + t*rx, P[1] + t*ry)
        vv = sx*sx + sy*sy
        u = 0.0 if vv <= eps else max(0.0, min(1.0, _dot(X[0]-A[0], X[1]-A[1], sx, sy) / vv))
        return t, u, X

    return None

def _is_isolated_tail(G, tail) -> bool:
    """Истинно, если рассматриваемый хвост принадлежит отрезку, у которого оба конца degree=1."""
    eid = tail["eid"]
    u, v = G.edges[eid]
    if u == -1 or v == -1:
        return False
    return (G.degree(u) == 1) and (G.degree(v) == 1)



def _nearest_point_on_segment(P, A, B):
    vx, vy = B[0]-A[0], B[1]-A[1]
    wx, wy = P[0]-A[0], P[1]-A[1]
    vv = vx*vx + vy*vy
    if vv <= 1e-18: return _len(P[0]-A[0], P[1]-A[1]), A, 0.0
    u = max(0.0, min(1.0, (wx*vx + wy*vy)/vv))
    X = (A[0] + u*vx, A[1] + u*vy)
    return _len(P[0]-X[0], P[1]-X[1]), X, u

def _face_contains_directed_edge(face_nodes, a, b):
    k = len(face_nodes)
    for i in range(k):
        if face_nodes[i] == a and face_nodes[(i+1)%k] == b:
            return True
    return False

def _decompose_along(P, D, X):
    """PX = s*D + n. Возвращает (s — вдоль, per — поперёк оси луча)."""
    vx, vy = X[0]-P[0], X[1]-P[1]
    s = _dot(vx, vy, D[0], D[1])
    per2 = max(0.0, vx*vx + vy*vy - s*s)
    return s, math.sqrt(per2)

import os, json  # если хочешь лог

def _bbox(P, Q):
    x1, y1 = P; x2, y2 = Q
    return (min(x1,x2), min(y1,y2), max(x1,x2), max(y1,y2))

def _seg_intersection_proper(P, Q, A, B, eps=1e-9):
    def orient(a,b,c): return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    o1 = orient(P,Q,A); o2 = orient(P,Q,B)
    o3 = orient(A,B,P); o4 = orient(A,B,Q)
    minx1, miny1, maxx1, maxy1 = _bbox(P,Q)
    minx2, miny2, maxx2, maxy2 = _bbox(A,B)
    if maxx1+eps < minx2 or maxx2+eps < minx1 or maxy1+eps < miny2 or maxy2+eps < miny1:
        return False
    return (o1*o2 < -eps) and (o3*o4 < -eps)


def _has_blocking_hit(G, grid, A, B, ignore_eids, eps):
    # проверяем пересечения на всём отрезке A–B
    def orient(p,q,r): return (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0])
    def bbox(p,q): return (min(p[0],q[0]), min(p[1],q[1]), max(p[0],q[0]), max(p[1],q[1]))
    minx1,miny1,maxx1,maxy1 = bbox(A,B)
    cx, cy = (A[0]+B[0])*0.5, (A[1]+B[1])*0.5
    r = 0.5*math.hypot(B[0]-A[0], B[1]-A[1]) + 5*eps
    for eid in grid.nearby_segments_by_point(cx, cy, r):
        if eid in ignore_eids: 
            continue
        u,v = G.edges[eid]
        if u == -1 or v == -1: 
            continue
        C, D = G.nodes[u], G.nodes[v]
        minx2,miny2,maxx2,maxy2 = bbox(C,D)
        if maxx1+eps < minx2 or maxx2+eps < minx1 or maxy1+eps < miny2 or maxy2+eps < miny1:
            continue
        o1 = orient(A,B,C); o2 = orient(A,B,D)
        o3 = orient(C,D,A); o4 = orient(C,D,B)
        if o1*o2 < -eps and o3*o4 < -eps:
            return True
    return False



def _gather_candidates(G, grid, tail, params):
    P = tail["P"]; D = tail["D"]
    prev = tail["prev"]; end = tail["end"]

    max_extend = params["MAX_EXTEND"]
    iso_mult   = params.get("ISOLATED_EXT_MULT", 2.0)
    # если хвост на изолированном сегменте — разрешим до 2×MAX_EXTEND
    ext_limit  = max_extend * (iso_mult if _is_isolated_tail(G, tail) else 1.0)

    ang_tol    = params["ANGLE_TOL"]
    near_perp_max    = params.get("NEAR_PERP_MAX", 10.0)
    near_forward_min = params.get("NEAR_FORWARD_MIN", 0.0)
    eps = params.get("EPS_SNAP", 2.0)

    # радиус поиска берём по расширенному лимиту
    search_r = params.get("SEARCH_RADIUS") or (ext_limit * 1.25)
    eids = [e for e in grid.nearby_segments_by_point(P[0], P[1], search_r)
            if e != tail["eid"] and G.edges[e][0] != -1]

    cand = []

    # --- RAY (по направлению хвоста) ---
    best_ray = None
    for eid in eids:
        A = G.nodes[G.edges[eid][0]]; B = G.nodes[G.edges[eid][1]]
        hit = _ray_segment_intersection(P, D, A, B, ang_tol_deg=ang_tol, col_perp_tol=near_perp_max)
        if not hit:
            continue
        t, u, X = hit
        if t < 1e-9 or t > ext_limit:     # ← ИСПОЛЬЗУЕМ ext_limit
            continue
        # видимость только для новой добавки P→T
        target = G.edges[eid][1] if (u >= 0.5 or abs(u-1.0) < 1e-9) else G.edges[eid][0]
        T = X if (1e-9 < u < 1-1e-9) else G.nodes[target]
        if _has_blocking_hit(G, grid, P, T, ignore_eids={tail["eid"], eid}, eps=eps):
            continue
        if (best_ray is None) or (t < best_ray["dist"]):
            best_ray = {"mode":"ray","target_eid":eid,"u":u,"X":X,"dist":t}
    if best_ray:
        cand.append(best_ray)

    # --- NEAREST (узкий коридор вдоль D) ---
    best_near = None
    for eid in eids:
        A = G.nodes[G.edges[eid][0]]; B = G.nodes[G.edges[eid][1]]
        _, X, u = _nearest_point_on_segment(P, A, B)
        s, per = _decompose_along(P, D, X)
        if s < near_forward_min or s > ext_limit or per > near_perp_max:   # ← ext_limit
            continue

        # оценим излом для выбора (как раньше)
        vx, vy = X[0]-P[0], X[1]-P[1]
        Lv = math.hypot(vx, vy) or 1.0
        ux, uy = vx/Lv, vy/Lv
        dot = max(-1.0, min(1.0, D[0]*ux + D[1]*uy))
        kink_deg = math.degrees(math.acos(dot))

        T = X
        if _has_blocking_hit(G, grid, P, T, ignore_eids={tail["eid"], eid}, eps=eps):
            continue
        if (best_near is None) or (s < best_near["dist"]):
            best_near = {"mode":"nearest","target_eid":eid,"u":u,"X":X,"dist":s, "kink_deg":kink_deg}
    if best_near:
        cand.append(best_near)

    return cand



def _make_tail(P, prev):
    """Направление 'наружу' из конца prev->P."""
    return _unit(P[0]-prev[0], P[1]-prev[1])

def find_tails(G):
    tails = []
    for eid, (u, v) in enumerate(G.edges):
        if u == -1: continue
        # конец в u
        if G.degree(u) == 1:
            P = G.nodes[u]; prev = v
            tails.append({"eid":eid, "end":u, "prev":prev, "P":P, "D":_make_tail(P, G.nodes[prev])})
        if G.degree(v) == 1:
            P = G.nodes[v]; prev = u
            tails.append({"eid":eid, "end":v, "prev":prev, "P":P, "D":_make_tail(P, G.nodes[prev])})
    return tails

def _try_template_for_choice(G, end_nid, choice, tdb, params):
    """Проверяем: если добавить это ребро, получается ли грань, похожая на шаблон из базы?"""
    # 1) клон графа и временное добавление
    G2 = G.clone()
    eid = choice["target_eid"]; u = choice["u"]; X = choice["X"]
    a, b = G2.edges[eid]
    if u > 1e-9 and u < 1-1e-9:
        mid, e1, e2 = G2.split_edge(eid, X)
        target = mid
    else:
        target = b if u >= 0.5 else a
    # подключаем конец
    G2.add_edge(end_nid, target)

    # 2) находим грани и выбираем те, где есть directed edge end_nid->target
    faces = find_planar_faces(G2, include_outer=False, right_hand=True)
    faces_hit = [cyc for cyc in faces if _face_contains_directed_edge(cyc, end_nid, target)]
    if not faces_hit:
        return False

    # 3) строим шаблоны и спрашиваем память
    ang_tol = params.get("TEMPLATE_ANG_TOL", 2)
    len_tol = params.get("TEMPLATE_LEN_TOL", 5)
    for cyc in faces_hit:
        tpl = make_template(G2, cyc)
        if not tpl: 
            continue
        if tdb.has_similar(tpl, ang_tol=ang_tol, len_tol=len_tol):
            return True
    return False

def choose_by_templates(G, grid, tail, tdb, params):
    """Сначала пытаемся найти вариант доводки, который даёт 'похожую' грань из памяти."""
    cand = _gather_candidates(G, grid, tail, params)
    for choice in cand:
        if _try_template_for_choice(G, tail["end"], choice, tdb, params):
            return choice  # этот вариант порождает грань, похожую на известную
    return None

def choose_rule_based(G, grid, tail, params):
    """Fallback: правило 'ray vs nearest' c коэффициентом и расширенным лимитом для изолированных хвостов."""
    cand = _gather_candidates(G, grid, tail, params)
    if not cand:
        return None

    max_extend = params["MAX_EXTEND"]
    iso_mult   = params.get("ISOLATED_EXT_MULT", 2.0)
    ext_limit  = max_extend * (iso_mult if _is_isolated_tail(G, tail) else 1.0)

    ratio        = params.get("DIR_TO_NEAR_RATIO", 2.0)
    near_kinkmax = params.get("NEAR_MAX_KINK_DEG", 12)

    ray  = next((c for c in cand if c["mode"] == "ray"), None)
    near = next((c for c in cand if c["mode"] == "nearest"), None)

    ray_ok  = bool(ray  and ray["dist"]  <= ext_limit)
    near_ok = bool(near and near["dist"] <= ext_limit)

    if ray_ok and not near_ok:
        return ray
    if near_ok and not ray_ok:
        return near  # луча нет/слишком далеко — берём ближайший

    if not ray_ok and not near_ok:
        return None

    # Оба валидны: если nearest делает заметный излом — выбираем луч
    if near.get("kink_deg", 0.0) > near_kinkmax:
        return ray

    # Иначе стандартное сравнение по ratio
    if ray["dist"] <= ratio * near["dist"]:
        return ray
    return near


def apply_choice(G, grid, params, tail, choice):
    eid_target = choice["target_eid"]
    u          = choice["u"]
    X          = choice["X"]
    mode       = choice["mode"]
    dist       = choice["dist"]

    a, b = G.edges[eid_target]
    if 1e-9 < u < 1.0 - 1e-9:
        mid, e1, e2 = G.split_edge(eid_target, X)
        target = mid
    else:
        target = b if u >= 0.5 else a

    prev = tail["prev"]; end = tail["end"]
    P = G.nodes[end]
    T = G.nodes[target]
    eps =  params.get("EPS_SNAP", 2.0)

    if mode == "ray":
        # безопасно «дотянуть по прямой»: новый отрезок prev→target коллинеарен старому
        # (дополнительно убеждаемся, что по пути нет блокеров)
        if _has_blocking_hit(G, grid, P, T, ignore_eids={tail["eid"], eid_target}, eps=eps):
            return None
        G.remove_edge(tail["eid"])
        eid_new = G.add_edge(prev, target)
        return {"mode": mode, "tail_eid_old": tail["eid"], "tail_eid_new": eid_new,
                "target_nid": target, "dist": dist}

    else:  # mode == "nearest"
        # НЕ поворачиваем всю длинную грань! Добавляем КОРОТКИЙ коннектор end→target.
        if _has_blocking_hit(G, grid, P, T, ignore_eids={tail["eid"], eid_target}, eps=eps):
            return None  # даже короткий отрезок чем-то перекрыт
        # оставляем prev→end как есть и добавляем end→target
        eid_conn = G.add_edge(end, target)
        return {"mode": mode, "connector_eid": eid_conn, "end": end,
                "target_nid": target, "dist": dist}


def close_tails_smart(G, grid, params, templates_db, iter_max=5):
    """
    Итеративно доводит хвосты.
    1) Сначала пробует варианты, которые дают грань «похожую» на шаблон из базы.
    2) Если не нашлось — правило «луч vs ближайший».
    Возвращает СПИСОК операций (каждая — dict из apply_choice).
    """
    ops = []
    for _ in range(iter_max):
        tails = find_tails(G)
        did = 0
        for t in tails:
            choice = choose_by_templates(G, grid, t, templates_db, params)
            if not choice:
                choice = choose_rule_based(G, grid, t, params)
            if not choice:
                continue
            op = apply_choice(G, grid, params, t, choice)
            if not op:
                continue
            ops.append(op)
            did += 1

        if did == 0:
            break
        # После серии изменений обновим сетку, чтобы следующие итерации видели новые рёбра
        try:
            from geom.grid import build_grid_from_graph
            grid = build_grid_from_graph(G, grid.cell, pad=params.get("EPS_SNAP", 0.0))
        except Exception:
            # Если по какой-то причине обновить сетку не удалось, продолжаем со старой
            pass
    return ops

