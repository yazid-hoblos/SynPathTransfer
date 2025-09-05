import requests
from collections import deque, defaultdict

KEGG_BASE = "http://rest.kegg.jp"

def get_compound_id(name):
    """Return KEGG compound ID for a given name."""
    url = f"{KEGG_BASE}/find/compound/{name}"
    r = requests.get(url)
    if not r.ok or not r.text.strip():
        return None
    return r.text.split("\t")[0].replace("cpd:", "").strip()  # e.g., C01013

def get_reactions_for_compound(cid):
    """Return all KEGG reactions involving the compound."""
    url = f"{KEGG_BASE}/link/reaction/cpd:{cid}"
    r = requests.get(url)
    if not r.ok or not r.text.strip():
        return []
    return [line.split("\t")[1].replace("rn:", "").strip()
            for line in r.text.strip().split("\n")]

def get_reaction_equation(rid):
    """Parse substrates and products from KEGG reaction EQUATION line."""
    url = f"{KEGG_BASE}/get/rn:{rid}"
    r = requests.get(url)
    if not r.ok:
        return [], [], False
    for line in r.text.split("\n"):
        if line.startswith("EQUATION"):
            eq = line.replace("EQUATION", "").strip()
            reversible = "<=>" in eq
            if reversible:
                left, right = eq.split("<=>")
            elif "=>" in eq:
                left, right = eq.split("=>")
            else:
                return [], [], False
            subs = [x.strip() for x in left.split("+")]
            prods = [x.strip() for x in right.split("+")]
            return subs, prods, reversible
    return [], [], False

def build_graph(target_cid):
    """Build a directed graph with substrates -> products edges.
       Reversible reactions add edges in both directions."""
    graph = defaultdict(set)
    visited_reactions = set()
    queue = deque([target_cid])

    while queue:
        cid = queue.popleft()
        reactions = get_reactions_for_compound(cid)
        for rid in reactions:
            if rid in visited_reactions:
                continue
            visited_reactions.add(rid)
            subs, prods, reversible = get_reaction_equation(rid)
            for s in subs:
                for p in prods:
                    graph[s].add(p)
                    if reversible:
                        graph[p].add(s)
            # Enqueue new compounds for further graph expansion
            for c in subs + prods:
                if c not in graph:
                    queue.append(c)
    # draw graph for debugging
    import matplotlib.pyplot as plt
    import networkx as nx
    G = nx.DiGraph()
    for src, dsts in graph.items():
        for dst in dsts:
            G.add_edge(src, dst)
    plt.figure(figsize=(12, 12))
    nx.draw(G, with_labels=True, node_size=500, font_size=8)
    plt.show()
    return graph

def reverse_bfs(graph, target_cid, steps=3):
    """Traverse reverse edges to find precursors within N steps."""
    # Build reverse graph
    reverse_graph = defaultdict(set)
    for src, dsts in graph.items():
        for dst in dsts:
            reverse_graph[dst].add(src)

    precursors = set()
    queue = deque([(target_cid, 0)])
    visited = {target_cid}

    while queue:
        cid, depth = queue.popleft()
        if depth >= steps:
            continue
        for precursor in reverse_graph.get(cid, []):
            if precursor not in visited:
                precursors.add(precursor)
                visited.add(precursor)
                queue.append((precursor, depth + 1))
    return precursors

def get_compound_name(cid):
    """Return the primary human-readable name of a KEGG compound."""
    url = f"{KEGG_BASE}/get/cpd:{cid}"
    r = requests.get(url)
    if not r.ok:
        return cid
    for line in r.text.split("\n"):
        if line.startswith("NAME"):
            return line.replace("NAME", "").strip().split(";")[0]
    return cid

# --- Main function ---
def find_precursors(target_name, steps=3):
    target_cid = get_compound_id(target_name)
    if not target_cid:
        return []

    graph = build_graph(target_cid)
    precursors = reverse_bfs(graph, target_cid, steps)
    # Convert to human-readable names
    return [(cid, get_compound_name(cid)) for cid in precursors]

# --- Example usage ---
if __name__ == "__main__":
    target = "3-hydroxypropionate"
    precursors = find_precursors(target, steps=2)
    print(f"Compounds that can form {target} in ≤ 3 steps:")
    for cid, name in precursors:
        print(f"{cid} → {name}")
