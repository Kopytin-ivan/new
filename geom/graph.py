# graph.py
from collections import defaultdict

class Graph:
    def __init__(self, nodes, edges):
        self.nodes = nodes              # id -> (x,y)
        self.edges = list(edges)        # list[(u,v)]
        self.adj   = defaultdict(list)  # id -> [edge_id]
        for eid,(u,v) in enumerate(self.edges):
            self.adj[u].append(eid)
            self.adj[v].append(eid)

    def degree(self, nid):
        return len(self.adj[nid])

    def edge_nodes(self, eid):
        return self.edges[eid]
