# main.py
import json
import os

from geom.io import load_segments
from geom.units import ensure_meters, quantize_mm_round   # или quantize_mm_trunc
from geom.snap import snap_points
from geom.graph import Graph
from geom.grid import pick_grid_params, build_grid_from_graph
from geom.config import PARAMS
from geom.cycles import find_simple_cycles
from geom.templates import make_template, TemplatesDB
from geom.prof import Prof

# создаём профайлер
prof = Prof(enabled=True)

def stage0_1_init(path_in, unit_scale=1.0):
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
        PARAMS["GRID_CELL"] = 20.0
        PARAMS["R_QUERY"]   = 60.0

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

    # ---- Этап 0–1
    with prof.section("STAGE0_1 total"):
        G, grid, params = stage0_1_init(r"data\\linesln.json")  # ← поменяй имя, если у тебя другой файл

    print("nodes:", len(G.nodes), "edges:", len(G.edges))
    print("params:", params)

    # sanity-check по сетке
    nid = 0
    x, y = G.nodes[nid]
    cand = grid.nearby_segments_by_point(x, y, params["R_QUERY"])
    print("candidates near node0:", len(cand))

# === Этап 2: поиск замкнутых фигур + запись (и все, и уникальные) ===
with prof.section("STAGE2 total"):

    # 2.1 Найти циклы
    with prof.section("stage2: find_cycles"):
        from geom.cycles import find_simple_cycles as find_cycles  # можно оставить find_simple_cycles
        cycles = find_cycles(G)  # список: [[n0,n1,...], ...]

    # 2.2 Генерация шаблонов, потоковая запись «всех», параллельно считаем «уникальные»
    with prof.section("stage2: make_templates"):
        from geom.templates import make_template, TemplatesDB
        os.makedirs("output", exist_ok=True)
        tmem = TemplatesDB()                                 # для уникальных (как раньше)

        write_all = True                                     # ← если не нужно писать «все», поставь False
        all_path = os.path.join("output", "closed_templates_all.ndjson")
        fout = open(all_path, "w", encoding="utf-8") if write_all else None

        total_all = 0
        try:
            for cyc in cycles:
                tpl = make_template(G, cyc)
                if not tpl:
                    continue
                tmem.add(tpl)                                # учёт уникальных
                if fout:                                     # потоковая запись «всех»
                    # можно писать облегчённую запись:
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

    # 2.3 Сохранить уникальные
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
