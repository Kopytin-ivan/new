# main.py
import json
import os

from geom.io import load_segments
from geom.units import ensure_meters, quantize_mm_round
from geom.snap import snap_points
from geom.graph import Graph
from geom.grid import pick_grid_params, build_grid_from_graph
from geom.config import PARAMS
from geom.prof import Prof

from geom.templates import TemplatesDB, make_template     # для close_tails_smart (можно пустую память)
from geom.cycles import find_planar_faces  # только чтобы засеять память из уже замкнутых
from geom.extend import close_tails_smart, connect_closed_islands_to_host
from geom.export_dxf import save_dxf_lines


from pathlib import Path
from geom.outer import save_outer_from_graph
prof = Prof(enabled=True)

def stage0_1_init(path_in, unit_scale=0.001):
    # 0) загрузка и масштаб (0.001 = мм → м)
    with prof.section("stage0: load+scale"):
        segs_in = load_segments(path_in)
        segs = ensure_meters(segs_in, scale=unit_scale)

    # 0.1) округление до мм (ускоряет снап/грид)
    with prof.section("stage0: quantize_mm"):
        segs = quantize_mm_round(segs)

    # 1) снап вершин
    with prof.section("stage1: snap_points"):
        nodes, edges = snap_points(segs, PARAMS["EPS_SNAP"])
        G = Graph(nodes, edges)

    # 2) параметры сетки
    with prof.section("stage1: pick_grid_params"):
        grid_cell, r_query = pick_grid_params(G.nodes, PARAMS["EPS_SNAP"])
        if PARAMS["GRID_CELL"] is None:
            PARAMS["GRID_CELL"] = grid_cell
        if PARAMS["R_QUERY"] is None:
            PARAMS["R_QUERY"] = r_query
        # при желании зафиксируй вручную:
        # PARAMS["GRID_CELL"] = 20.0
        # PARAMS["R_QUERY"]   = 60.0

    # 3) построение сетки
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

    # === Этап 0–1: загрузка → снап → сетка
    with prof.section("STAGE0_1 total"):
        G, grid, params = stage0_1_init(r"data\\linesln.json", unit_scale=0.001)

    print("nodes:", len(G.nodes), "edges:", len(G.edges))
    print("params:", params)
    if G.nodes:
        x0, y0 = G.nodes[0]
        print("candidates near node0:", len(grid.nearby_segments_by_point(x0, y0, params["R_QUERY"])))

    # === Память шаблонов (минимально: подгрузим, затем засеем текущими замкнутыми)
    tmem = TemplatesDB()
    tmem.load(os.path.join("output", "closed_templates.json"))  # ок, если файла нет
    faces0 = find_planar_faces(G, include_outer=False, right_hand=True)
    for cyc in faces0:
        tpl = make_template(G, cyc)
        if tpl:
            tmem.add(tpl)

    # === 1) Умное доведение хвостов (по памяти → правило)
    with prof.section("stage1.9: smart extend tails"):
        ops1 = close_tails_smart(G, grid, PARAMS, templates_db=tmem, iter_max=5)
        print("tails extended:", len(ops1))

    # перестроим сетку после правок
    with prof.section("stage1.95: rebuild grid after tails"):
        grid = build_grid_from_graph(G, PARAMS["GRID_CELL"], pad=PARAMS["EPS_SNAP"])

    # === 2) Замкнутые острова: два параллельных мостика к host с минимальной суммой доводок
    with prof.section("stage1.96: bridge closed islands → host"):
        ops2 = connect_closed_islands_to_host(G, grid, PARAMS, rebuild_grid_each=True)
        print("closed-island bridges:", len(ops2))

    # финальная перестройка сетки
    with prof.section("stage1.97: rebuild grid after bridges"):
        grid = build_grid_from_graph(G, PARAMS["GRID_CELL"], pad=PARAMS["EPS_SNAP"])

    # === Экспорт ВСЕХ текущих сегментов (уже с доводкой и мостиками)
    out_segs = []
    for eid, (u, v) in enumerate(G.edges):
        if u == -1 or v == -1:
            continue
        x1, y1 = G.nodes[u]
        x2, y2 = G.nodes[v]
        out_segs.append([[round(x1, 5), round(y1, 5)], [round(x2, 5), round(y2, 5)]])

    out_json = os.path.join("output", "polylineOut (округление до 5 знаков).json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump({"segments": out_segs}, f, ensure_ascii=False, indent=2)
    print("saved:", out_json)

    dxf_segs = [((float(a[0]), float(a[1])), (float(b[0]), float(b[1]))) for a, b in out_segs]
    out_dxf = os.path.join("output", "polylineOut.dxf")
    save_dxf_lines(dxf_segs, out_dxf, layer="OUTLINE", color=7, lineweight=25, insunits="Meters")
    print("saved:", out_dxf)

    # === ВНЕШНИЙ КОНТУР ИЗ ФИНАЛЬНОГО ГРАФА (без переснэпа/перечтения)
    out_stem = Path(out_json).stem  # 'polylineOut (округление до 5 знаков)'
    out_outer_json = Path("output") / f"{out_stem}_outer.json"
    out_outer_dxf  = Path("output") / f"{out_stem}_outer.dxf"

    meta_outer = save_outer_from_graph(
        G,
        out_json_path=str(out_outer_json),
        out_dxf_path=str(out_outer_dxf),
        drop_leaves=True,       # срезать degree=1
        leaf_len_mm=0.0         # можно 2–5 мм, если мешают короткие «шипы»
    )
    print("outer meta:", meta_outer)
    print("saved:", out_outer_json)
    print("saved:", out_outer_dxf)


    # профилирование
    prof.report_console()
    prof.dump_csv("output/timings.csv")
    print("timings saved to output/timings.csv")
