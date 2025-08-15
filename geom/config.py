# config.py
EPS_SNAP   = 2   # м, 1 мм
GRID_CELL  = None    # подберём автоматически (см. ниже)
R_QUERY    = None    # подберём как 3 * GRID_CELL
MAX_EXTEND = 10    # м, 10 см
ANGLE_TOL  = 2       # градусы
LEN_TOL    = 5       # проценты
PARAMS = {
    "EPS_SNAP": 0.002,      # 1 мм (если ваши единицы — метры)
    "MAX_EXTEND": 0.05,     # 5 см максимум продления (под ваш кейс)
    "ANGLE_TOL": 2,         # ±2°
    "LEN_TOL": 5,           # ±5%

    # Эти два ниже пересчитаем автоматически после загрузки
    "GRID_CELL": None,
    "R_QUERY":   None,
}
def tune_spatial_grid(nodes_xy, params):
    # Оценим характерный шаг узлов: медиана расстояния до ближайшего соседа
    import numpy as np
    if len(nodes_xy) < 2:
        cell = 1.0; r = 3.0
    else:
        # Быстро и без тяжёлого KD-дерева: возьмём подвыборку
        sample = np.array(nodes_xy[::max(1, len(nodes_xy)//5000)])
        # «Грубый» шаг как медиана по осям (дешёво, но стабильно)
        dx = np.median(np.abs(np.diff(np.sort(sample[:,0])))) or 0.05
        dy = np.median(np.abs(np.diff(np.sort(sample[:,1])))) or 0.05
        dnn = max(min(dx, dy), params["EPS_SNAP"]*5)  # не меньше 5×EPS
        cell = max(0.02, min(5.0, 3.0*dnn))          # 2 см…5 м безопасные рамки
        r    = 3.0*cell
    params["GRID_CELL"] = cell
    params["R_QUERY"]   = r
