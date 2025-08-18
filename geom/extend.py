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

def _ray_segment_intersection(P, D, A, B, ang_tol_deg=2.0, eps=1e-12):
    """
    Пересечение луча P + t*D (t>=0) с отрезком [A,B].
    Возвращает (t, u, X) где t>=0 — вдоль луча (при условии, что D единичный),
    0<=u<=1 — параметр на отрезке, X — точка. Обрабатывает и коллинеарный случай.
    """
    rx, ry = D
    sx, sy = (B[0]-A[0], B[1]-A[1])
    denom = _cross(rx, ry, sx, sy)

    # --- НЕ совпадающие направления: обычное пересечение прямых ---
    if abs(denom) > eps:
        APx, APy = (A[0]-P[0], A[1]-P[1])
        t = _cross(APx, APy, sx, sy) / denom
        u = _cross(APx, APy, rx, ry) / denom
        if t < 0.0 or u < 0.0 or u > 1.0:
            return None
        X = (P[0] + t*rx, P[1] + t*ry)
        return t, u, X

    # --- Почти параллельны: проверяем, не КОЛЛИНЕАРНЫ ли ---
    # Вектор от P до A должен быть тоже параллелен D (крест. произведение ~ 0)
    if abs(_cross(A[0]-P[0], A[1]-P[1], rx, ry)) > 1e-9:
        # параллельны, но не на одной прямой — пересечения нет
        return None

    # Коллинеарный случай: проектируем концы сегмента на ось луча и берём ближайший ВПЕРЁД
    # D у нас единичный, поэтому параметр t = проекция на D.
    sA = _dot(A[0]-P[0], A[1]-P[1], rx, ry)
    sB = _dot(B[0]-P[0], B[1]-P[1], rx, ry)
    candidates = [s for s in (sA, sB) if s >= 0.0]       # только вперёд по лучу
    if not candidates:
        return None
    t = min(candidates)                                  # ближайшая вперёд точка
    X = (P[0] + t*rx, P[1] + t*ry)
    # параметр u на [A,B]
    vv = sx*sx + sy*sy
    if vv <= eps:
        u = 0.0
    else:
        u = max(0.0, min(1.0, _dot(X[0]-A[0], X[1]-A[1], sx, sy) / vv))
    return t, u, X


def _decompose_along(P, D, X):
    vx, vy = X[0]-P[0], X[1]-P[1]
    s = _dot(vx, vy, D[0], D[1])              # вдоль
    perp2 = max(0.0, vx*vx + vy*vy - s*s)
    per = math.sqrt(perp2)
    return s, per


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
    """Разложение PX: s — вперёд по D (может быть <0), per — поперёк оси луча."""
    vx, vy = X[0]-P[0], X[1]-P[1]
    s = _dot(vx, vy, D[0], D[1])
    perp2 = max(0.0, vx*vx + vy*vy - s*s)
    return s, math.sqrt(perp2)

def _gather_candidates(G, grid, tail, params):
    """
    Возвращает кандидатов-целей с ЖЁСТКИМ отсевом по MAX_EXTEND.
    Для 'nearest' учитываем только точки почти по направлению хвоста (узкий коридор).
    """
    P = tail["P"]; D = tail["D"]
    max_extend = params["MAX_EXTEND"]
    ang_tol    = params["ANGLE_TOL"]
    near_perp_max    = params.get("NEAR_PERP_MAX",  max(0.01, params["EPS_SNAP"]*5))
    near_forward_min = params.get("NEAR_FORWARD_MIN", 0.0)

    # Радиус поиска — не шире лимита
    search_r = params.get("SEARCH_RADIUS")
    if not search_r:
        search_r = max_extend * 1.25  # чуть «с запасом», но не как R_QUERY=60

    eids = [e for e in grid.nearby_segments_by_point(P[0], P[1], search_r)
            if e != tail["eid"] and G.edges[e][0] != -1]

    cand = []

    # 1) Лучом (по направлению). t — это расстояние вдоль луча (D — единичный)
    best_ray = None
    for eid in eids:
        A = G.nodes[G.edges[eid][0]]; B = G.nodes[G.edges[eid][1]]
        hit = _ray_segment_intersection(P, D, A, B, ang_tol_deg=ang_tol)
        if not hit:
            continue
        t, u, X = hit
        if t < 1e-9 or t > max_extend:
            continue
        if (best_ray is None) or (t < best_ray["dist"]):
            best_ray = {"mode":"ray","target_eid":eid,"u":u,"X":X,"dist":t}
    if best_ray:
        cand.append(best_ray)

    # 2) К ближайшему — но только в узком коридоре вдоль D и в пределах max_extend
    best_near = None
    for eid in eids:
        A = G.nodes[G.edges[eid][0]]; B = G.nodes[G.edges[eid][1]]
        _, X, u = _nearest_point_on_segment(P, A, B)
        s, per = _decompose_along(P, D, X)
        if s < near_forward_min or s > max_extend:
            continue          # позади / слишком далеко вперёд
        if per > near_perp_max:
            continue          # слишком сбоку
        if (best_near is None) or (s < best_near["dist"]):
            best_near = {"mode":"nearest","target_eid":eid,"u":u,"X":X,"dist":s}
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
    """Fallback: правило 'ray vs nearest' c коэффициентом и MAX_EXTEND."""
    cand = _gather_candidates(G, grid, tail, params)
    if not cand: return None
    max_extend = params["MAX_EXTEND"]
    ratio = params.get("DIR_TO_NEAR_RATIO", 2.0)
    ray = next((c for c in cand if c["mode"]=="ray"), None)
    near = next((c for c in cand if c["mode"]=="nearest"), None)

    if ray and not near:
        return ray if ray["dist"] <= max_extend else None
    if near and not ray:
        return near if near["dist"] <= max_extend else None
    if not ray and not near:
        return None

    if (ray["dist"] <= ratio * near["dist"]) and (ray["dist"] <= max_extend):
        return ray
    return near if near["dist"] <= max_extend else None

def apply_choice(G, tail, choice):
    """
    Вносит изменения в G, ПЕРЕПОДКЛЮЧАЯ исходный хвостовой отрезок на цель.
    Никаких «добавочных» ребер: меняем сам отрезок tail['eid'].
    Возвращает dict-лог операции (пригодится для отладки).
    """
    eid_target = choice["target_eid"]
    u          = choice["u"]
    X          = choice["X"]
    mode       = choice["mode"]
    dist       = choice["dist"]

    # 1) подготовить целевую вершину (если попали в середину сегмента — разрезать)
    a, b = G.edges[eid_target]
    if u > 1e-9 and u < 1.0 - 1e-9:
        mid, e1, e2 = G.split_edge(eid_target, X)
        target_nid = mid
        target_xy  = G.nodes[mid]
        did_split  = True
        created_id = mid
        created_xy = target_xy
    else:
        target_nid = b if u >= 0.5 else a
        target_xy  = G.nodes[target_nid]
        did_split  = False
        created_id = None
        created_xy = None

    # 2) переподключить исходный «хвостовой» сегмент
    #    tail['eid'] соединяет (prev, end). Меняем его на (prev, target_nid).
    eid_tail = tail["eid"]
    prev     = tail["prev"]
    end      = tail["end"]

    # удалить старое ребро и добавить новое (переподключение)
    G.remove_edge(eid_tail)
    eid_new = G.add_edge(prev, target_nid)

    # (опц.) старый узел end станет изолированным — можно оставить, не мешает

    return {
        "tail_eid_old": eid_tail,
        "tail_prev": prev,
        "tail_end": end,
        "tail_end_xy": G.nodes[end],
        "mode": mode,
        "dist": dist,
        "target_eid": eid_target,
        "u_on_target": u,
        "split": did_split,
        "created_node_id": created_id,
        "created_xy": created_xy,
        "target_nid": target_nid,
        "target_xy": target_xy,
        "tail_eid_new": eid_new,
    }


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
            op = apply_choice(G, t, choice)
            ops.append(op)
            did += 1
        if did == 0:
            break
    return ops

