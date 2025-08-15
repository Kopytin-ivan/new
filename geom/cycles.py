# geom/cycles.py
from collections import deque

def connected_components(G):
    seen = set()
    comps = []
    for start in range(len(G.nodes)):
        if start in seen:
            continue
        if start not in G.adj:
            seen.add(start)
            comps.append([start])
            continue
        q = deque([start])
        seen.add(start)
        comp = [start]
        while q:
            u = q.popleft()
            for eid in G.adj.get(u, []):
                a, b = G.edge_nodes(eid)
                v = b if a == u else a
                if v not in seen:
                    seen.add(v)
                    q.append(v)
                    comp.append(v)
        comps.append(comp)
    return comps

def is_degree2_component(G, comp_nodes):
    for nid in comp_nodes:
        if G.degree(nid) != 2:
            return False
    return True

def _component_edges(G, comp_nodes):
    comp = set(comp_nodes)
    eids = []
    for nid in comp_nodes:
        for eid in G.adj.get(nid, []):
            a, b = G.edge_nodes(eid)
            if a in comp and b in comp:
                eids.append(eid)
    return eids

def _next_neighbor(G, prev_node, curr_node):
    """В degree=2 у узла ровно два соседа; берём тот, который != prev."""
    out = []
    for eid in G.adj[curr_node]:
        a, b = G.edge_nodes(eid)
        v = b if a == curr_node else a
        out.append((eid, v))
    # out = [(eid1, v1), (eid2, v2)]
    if not out:
        return None
    if len(out) == 1:
        return out[0]  # деградация, но обычно не бывает в degree=2
    (e1, v1), (e2, v2) = out
    if v1 != prev_node:
        return (e1, v1)
    else:
        return (e2, v2)

def find_simple_cycles(G):
    cycles = []
    for comp in connected_components(G):
        if not is_degree2_component(G, comp):
            continue
        # boolean-маркировка рёбер компоненты
        comp_edges = _component_edges(G, comp)
        if not comp_edges:
            continue
        visited = {eid: False for eid in comp_edges}

        # стартуем из узлов компоненты, вытягиваем орбиты
        for start in comp:
            # возьмём первое неиспользованное ребро у старта
            e0 = None
            for eid in G.adj[start]:
                if eid in visited and not visited[eid]:
                    e0 = eid
                    break
            if e0 is None:
                continue

            # первый шаг: старт → сосед
            a, b = G.edge_nodes(e0)
            curr = start
            nxt = b if a == curr else a
            visited[e0] = True
            cycle_nodes = [curr, nxt]
            prev = curr
            curr = nxt

            # идём, пока не вернёмся в start
            while True:
                step = _next_neighbor(G, prev, curr)
                if step is None:
                    cycle_nodes = None
                    break
                eid_next, v_next = step
                if eid_next in visited:
                    if visited[eid_next]:
                        # если ребро уже использовали — значит вернулись на пройденный путь:
                        # либо цикл замкнулся, либо что-то пошло не так.
                        if v_next == cycle_nodes[0]:
                            # замыкание
                            break
                        else:
                            cycle_nodes = None
                            break
                    visited[eid_next] = True
                else:
                    # ребро вне компоненты (не должно быть)
                    cycle_nodes = None
                    break
                cycle_nodes.append(v_next)
                prev, curr = curr, v_next
                if curr == cycle_nodes[0]:
                    break

            if cycle_nodes and len(cycle_nodes) >= 3:
                # убрать последний дубликат, если есть
                if cycle_nodes[0] == cycle_nodes[-1]:
                    cycle_nodes = cycle_nodes[:-1]
                # быстрый фильтр самоповторов
                if len(set(cycle_nodes)) == len(cycle_nodes):
                    cycles.append(cycle_nodes)

    return cycles
