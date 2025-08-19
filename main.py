# main.py
import json
import os

from geom.io import load_segments
from geom.units import ensure_meters, quantize_mm_round   # или quantize_mm_trunc
from geom.snap import snap_points
from geom.graph import Graph
from geom.grid import pick_grid_params, build_grid_from_graph
from geom.config import PARAMS
from geom.prof import Prof

# faces/шаблоны/доводка
from geom.cycles import find_planar_faces
from geom.templates import TemplatesDB, make_template
from geom.extend import close_tails_smart

from geom.export_dxf import save_dxf_lines


# создаём профайлер
prof = Prof(enabled=True)

def stage0_1_init(path_in, unit_scale=0.001):
    # 0) загрузка и перевод единиц
    with prof.section("stage0: load+scale"):
        segs_in = load_segments(path_in)
        segs = ensure_meters(segs_in, scale=unit_scale)

    # 0.1) квантуем до 1 мм (ускоряет снап/грид)
    with prof.section("stage0: quantize_mm"):
        segs = quantize_mm_round(segs)
        # если хочешь “обрезать” дробную часть вместо округления — используй quantize_mm_trunc

    # 1) снап вершин
    with prof.section("stage1: snap_points"):
        nodes, edges = snap_points(segs, PARAMS["EPS_SNAP"])
        G = Graph(nodes, edges)

    # 2) автоподбор параметров сетки
    with prof.section("stage1: pick_grid_params"):
        grid_cell, r_query = pick_grid_params(G.nodes, PARAMS["EPS_SNAP"])
        if PARAMS["GRID_CELL"] is None:
            PARAMS["GRID_CELL"] = grid_cell
        if PARAMS["R_QUERY"] is None:
            PARAMS["R_QUERY"] = r_query
        # если хочешь фиксированные значения — раскомментируй:
        # PARAMS["GRID_CELL"] = 20.0
        # PARAMS["R_QUERY"]   = 60.0

    # 3) строим равномерную сетку
    with prof.section("stage1: build_grid"):
        grid = build_grid_from_graph(G, PARAMS["GRID_CELL"], pad=PARAMS["EPS_SNAP"])

    params = {
        "EPS_SNAP":   PARAMS["EPS_SNAP"],
        "GRID_CELL":  PARAMS["GRID_CELL"],
        "R_QUERY":    PARAMS["R_QUERY"],
        "MAX_EXTEND": PARAMS["MAX_EXTEND"],
        "ANGLE_TOL":  PARAMS["ANGLE_TOL"],
        "LEN_TOL":    PARAMS["LEN_TOL"],
    }
    return G, grid, params


if __name__ == "__main__":
    os.makedirs("output", exist_ok=True)

    # ---- Этап 0–1: загрузка → снап → сетка 
    with prof.section("STAGE0_1 total"):
        G, grid, params = stage0_1_init(r"data\\linesln.json", unit_scale=0.001)  

    # --- проверка масштаба
    xs = [x for x, y in G.nodes]
    ys = [y for x, y in G.nodes]
    bb = (min(xs), min(ys), max(xs), max(ys))
    dx = bb[2] - bb[0]
    dy = bb[3] - bb[1]
    print(f"BBox: {bb}")
    print(f"Size: {dx:.3f} x {dy:.3f} (в текущих внутренних единицах)")

    print("nodes:", len(G.nodes), "edges:", len(G.edges))
    print("params:", params)

    # sanity-check по сетке
    nid = 0
    x, y = G.nodes[nid]
    cand = grid.nearby_segments_by_point(x, y, params["R_QUERY"])
    print("candidates near node0:", len(cand))

    # --- Память шаблонов: подгружаем прошлые + засеваем текущими замкнутыми (если есть)
    tmem = TemplatesDB()
    tmem.load(os.path.join("output", "closed_templates.json"))  # ок, если ещё нет файла

    faces0 = find_planar_faces(G, include_outer=False, right_hand=True)
    for cyc in faces0:
        tpl = make_template(G, cyc)
        if tpl:
            tmem.add(tpl)

    # --- «Умное доведение» хвостов: сначала по памяти, потом по правилу
    with prof.section("stage1.9: smart extend tails"):
        extended = close_tails_smart(G, grid, PARAMS, templates_db=tmem, iter_max=5)
        print("tails extended:", extended)

    # --- После доводок пересобираем сетку (геометрия изменилась)
    with prof.section("stage1.95: rebuild grid after extend"):
        grid = build_grid_from_graph(G, PARAMS["GRID_CELL"], pad=PARAMS["EPS_SNAP"])

    # === Этап 2: поиск замкнутых граней + запись (и все, и уникальные)
    with prof.section("STAGE2 total"):

        # 2.1 Найти грани (faces)
        with prof.section("stage2: find_faces"):
            cycles = find_planar_faces(G, include_outer=False, right_hand=True)

        # 2.2 Генерация шаблонов, потоковая запись «всех», параллельно считаем «уникальные»
        with prof.section("stage2: make_templates"):
            os.makedirs("output", exist_ok=True)

            write_all = True  # ← если не нужно писать «все», поставь False
            all_path = os.path.join("output", "closed_templates_all.ndjson")
            fout = open(all_path, "w", encoding="utf-8") if write_all else None

            total_all = 0
            try:
                for cyc in cycles:
                    tpl = make_template(G, cyc)
                    if not tpl:
                        continue
                    # пополняем память уникальных
                    tmem.add(tpl)

                    # потоковая запись «всех» упрощённых записей
                    if fout:
                        rec = {
                            "n": tpl["n"],
                            "angles": tpl["angles"],
                            "percents": tpl["percents"],
                            "K": tpl["K"],
                        }
                        fout.write(json.dumps(rec, ensure_ascii=False) + "\n")
                        total_all += 1
                        if (total_all % 1000) == 0:
                            fout.flush()
            finally:
                if fout:
                    fout.close()

        # 2.3 Сохранить уникальные (память)
        with prof.section("stage2: save_templates"):
            tmem.save(os.path.join("output", "closed_templates.json"))

    # Итоги
    print(f"closed shapes found: {len(cycles)} | unique templates: {tmem.size()} | all saved: {total_all}")
    for i, rec in enumerate(tmem.to_list()[:5], 1):
        print(f"{i}. {rec['text']}  (K={rec['K']}, count={rec['count']})")

    # ---- отчёт по времени
    prof.report_console()
    prof.dump_csv("output/timings.csv")
    print("timings saved to output/timings.csv")


        # --- После доводок пересобираем сетку (геометрия изменилась)
    with prof.section("stage1.95: rebuild grid after extend"):
        grid = build_grid_from_graph(G, PARAMS["GRID_CELL"], pad=PARAMS["EPS_SNAP"])

    # --- ВЫГРУЗИТЬ ВСЕ СЕГМЕНТЫ ПОСЛЕ ДОВОДКИ (округление до 5 знаков)
    out_segs = []
    for eid, (u, v) in enumerate(G.edges):
        if u == -1 or v == -1:
            continue
        x1, y1 = G.nodes[u]
        x2, y2 = G.nodes[v]
        out_segs.append([[round(x1, 5), round(y1, 5)], [round(x2, 5), round(y2, 5)]])

    # имя файла — как просишь: «polylineOut (округление до 5 знаков).json»
    out_path = os.path.join("output", "polylineOut (округление до 5 знаков).json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump({"segments": out_segs}, f, ensure_ascii=False, indent=2)
    print("saved:", out_path)
    dxf_segs = [((float(a[0]), float(a[1])), (float(b[0]), float(b[1]))) for a, b in out_segs]
    out_dxf = os.path.join("output", "polylineOut.dxf")
    save_dxf_lines(dxf_segs, out_dxf, layer="OUTLINE", color=7, lineweight=25, insunits="Meters")


    print("saved:", out_dxf)