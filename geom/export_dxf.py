# geom/export_dxf.py
# Экспорт отрезков в DXF. Если ezdxf установлен — используем его;
# иначе пишем минимальный ASCII-DXF вручную.
from typing import List, Tuple

Point = Tuple[float, float]
Segment = Tuple[Point, Point]

_INSUNITS = {
    "Unitless": 0, "Inches": 1, "Feet": 2, "Miles": 3, "Millimeters": 4,
    "Centimeters": 5, "Meters": 6, "Kilometers": 7
}

def save_dxf_lines(
    segments: List[Segment],
    path: str,
    layer: str = "OUTLINE",
    color: int = 7,            # 7=белый/чёрный
    lineweight: int = 25,      # 0.25 мм → код 25 (DXF group 370)
    insunits: str = "Meters"   # твои координаты — в метрах
) -> None:
    try:
        import ezdxf  # type: ignore
        _save_with_ezdxf(segments, path, layer, color, lineweight, insunits)
    except Exception:
        _save_plain_ascii_dxf(segments, path, layer, color, lineweight, insunits)


def _save_with_ezdxf(
    segments: List[Segment], path: str, layer: str, color: int, lw: int, insunits: str
) -> None:
    import ezdxf  # type: ignore
    doc = ezdxf.new(setup=True)  # AC1027 by default
    # Единицы файла — метры
    doc.header["$INSUNITS"] = _INSUNITS.get(insunits, 6)
    # Настроим слой
    if layer not in doc.layers:
        doc.layers.add(name=layer, color=color, lineweight=lw)
    msp = doc.modelspace()
    for (p1, p2) in segments:
        msp.add_line((p1[0], p1[1], 0.0), (p2[0], p2[1], 0.0), dxfattribs={"layer": layer, "color": color, "lineweight": lw})
    doc.saveas(path)


def _save_plain_ascii_dxf(
    segments: List[Segment], path: str, layer: str, color: int, lw: int, insunits: str
) -> None:
    """Минимальный ASCII-DXF без внешних зависимостей (LINE в секции ENTITIES)."""
    iu = _INSUNITS.get(insunits, 6)
    lines = []
    push = lines.append
    # HEADER
    push("0"); push("SECTION")
    push("2"); push("HEADER")
    push("9"); push("$INSUNITS")
    push("70"); push(str(iu))
    push("0"); push("ENDSEC")

    # TABLES (минимум, без определения слоёв — слои создадутся на лету)
    push("0"); push("SECTION")
    push("2"); push("TABLES")
    push("0"); push("ENDSEC")

    # ENTITIES
    push("0"); push("SECTION")
    push("2"); push("ENTITIES")
    for (p1, p2) in segments:
        push("0"); push("LINE")
        push("8"); push(layer)                 # layer
        push("62"); push(str(color))           # color (ACI)
        push("370"); push(str(lw))             # lineweight (0.01 mm units)
        push("10"); push(f"{p1[0]:.12f}")
        push("20"); push(f"{p1[1]:.12f}")
        push("30"); push("0.0")
        push("11"); push(f"{p2[0]:.12f}")
        push("21"); push(f"{p2[1]:.12f}")
        push("31"); push("0.0")
    push("0"); push("ENDSEC")
    push("0"); push("EOF")

    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write("\n".join(lines))
