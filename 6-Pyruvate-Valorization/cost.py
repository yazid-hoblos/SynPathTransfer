
import re
import time
import requests
from collections import defaultdict

BASE = "https://rest.kegg.jp"

# --- KEGG compound IDs we need ---
CID = {
    "ATP":"C00002","ADP":"C00008","AMP":"C00020",
    "GTP":"C00044","GDP":"C00035","UTP":"C00075","CTP":"C00063",
    "Pi":"C00009","PPi":"C00013",
    "NAD+":"C00003","NADH":"C00004","NADP+":"C00006","NADPH":"C00005",
    "FAD":"C00016","FADH2":"C01352",
    "O2":"C00007","CO2":"C00011","H2O":"C00001","H+":"C00080",
}

TRIVIAL = {CID["H2O"], CID["H+"], CID["Pi"]}

# --- Weights for RELATIVE score (can be tweaked) ---
REL_W = dict(
    alpha=1.0,   # ATP equivalents
    beta=1.0,    # redox (in ATP-eq via P/O)
    gamma=0.3,   # O2 consumption
    delta=0.25,  # CO2 released
    eps=0.20,    # complexity
    zeta=1.0     # precedent term
)

# --- Helpers to parse KEGG entries ---
CID_TOKEN   = re.compile(r'(C\d{5})')
COEFF_TOKEN = re.compile(r'^\s*(\d+(?:\.\d+)?)\s+(\S.*)$')

# Functions to access the reactions using their RID
def kegg_get(rid):
    r = requests.get(f"{BASE}/get/{rid}", timeout=60)
    r.raise_for_status()
    return r.text

def parse_entry(raw):
    """Extract ENTRY, EQUATION, PATHWAY lines from KEGG reaction entry."""
    data = {"ENTRY": None, "EQUATION": None, "PATHWAY": []}
    key = None
    for line in raw.splitlines():
        if not line.strip(): 
            continue
        if re.match(r'^[A-Z]{2,}\s{2,}', line):
            key = line[:12].strip()
            val = line[12:].rstrip()
        else:
            val = line.rstrip()
        if key == "ENTRY":
            data["ENTRY"] = val.split()[0]  # RXXXXX
        elif key == "EQUATION":
            data["EQUATION"] = (data["EQUATION"] + " " + val) if data["EQUATION"] else val
        elif key == "PATHWAY" and val:
            data["PATHWAY"].append(val.split()[0])  # mapXXXXX
    return data

# Functions to manipulate the reactions data
def parse_side(side_text):
    """Turn '2 C00002 + C00003' into {C00002:2, C00003:1}."""
    out = defaultdict(float)
    for term in side_text.split('+'):
        term = term.strip()
        if not term: 
            continue
        m = COEFF_TOKEN.match(term)
        if m:
            coeff = float(m.group(1)); rest = m.group(2)
        else:
            coeff = 1.0; rest = term
        m2 = CID_TOKEN.search(rest)
        if m2:
            out[m2.group(1)] += coeff
    return dict(out)

def parse_equation(eq):
    if '<=>' in eq:
        L, R = eq.split('<=>', 1)
    elif '=>' in eq:
        L, R = eq.split('=>', 1)
    else:
        return {}, {}
    return parse_side(L), parse_side(R)

def reverse_equation(L, R):
    return R, L

# --- Feature extraction for the cost formula ---
def net(dleft, dright, cids):
    """moles on left minus right for the set of cids."""
    return sum(dleft.get(c,0.0) for c in cids) - sum(dright.get(c,0.0) for c in cids)

def complexity(L, R):
    """unique non-trivial species count across both sides."""
    cleanL = {c for c in L if c not in TRIVIAL}
    cleanR = {c for c in R if c not in TRIVIAL}
    return len(cleanL | cleanR)

def features_from_entry(entry_text, direction=+1):
    ent = parse_entry(entry_text)
    L, R = parse_equation(ent["EQUATION"] or "")
    if direction == -1:
        L, R = reverse_equation(L, R)
    pcount = len(set(ent["PATHWAY"]))

    # 1) ΔATPeq
    triphos = {CID["ATP"], CID["GTP"], CID["UTP"], CID["CTP"]}
    recover = {CID["ADP"], CID["AMP"], CID["GDP"]}
    dATPeq = net(L, R, triphos) - 0.9*net(L, R, recover)

    # 2) ΔREDOX_ATP (ATP-equivalents via P/O: 2.5, 1.5)
    red_cons = 2.5*net(L, R, {CID["NADH"]}) + 2.5*net(L, R, {CID["NADPH"]}) + 1.5*net(L, R, {CID["FADH2"]})
    red_prod = 2.5*net(R, L, {CID["NAD+"], CID["NADP+"]}) + 1.5*net(R, L, {CID["FAD"]})
    dREDOX_ATP = red_cons - red_prod

    # 3) O2 consumption (≥0)
    O2cons = max(0.0, net(L, R, {CID["O2"]}))

    # 4) CO2 released (≥0)
    CO2rel = max(0.0, net(R, L, {CID["CO2"]}))  # right-left

    # 5) complexity
    cx = complexity(L, R)

    # 6) precedent
    precedent = 1.0/(1.0 + pcount)

    return dict(
        dATPeq=dATPeq,
        dREDOX_ATP=dREDOX_ATP,
        O2cons=O2cons,
        CO2rel=CO2rel,
        complexity=cx,
        precedent=precedent
    )

def cost_relative_from_features(F, W=REL_W):
    """Exact formula from your image (relative score)."""
    return (W["alpha"]*F["dATPeq"]
            + W["beta"] *F["dREDOX_ATP"]
            + W["gamma"]*F["O2cons"]
            + W["delta"]*F["CO2rel"]
            + W["eps"]  *F["complexity"]
            + W["zeta"] *F["precedent"])

# --- Public convenience functions ---

def reaction_cost_relative(rid, direction=+1, weights=REL_W):
    """Relative cost for a single reaction ID (KEGG)."""
    txt = kegg_get(rid)
    F = features_from_entry(txt, direction=direction)
    return cost_relative_from_features(F, weights), F

def subpathway_cost_relative(steps, weights=REL_W):
    """
    steps = [(reaction_id, direction), ...] with direction in {+1, -1}
    Returns (total_cost, per_step_details)
    """
    total = 0.0
    details = []
    # simple caching so we don't refetch the same R-id twice
    cache = {}
    for rid, direction in steps:
        if rid not in cache:
            cache[rid] = kegg_get(rid)
            time.sleep(0.2)  # be polite to KEGG
        F = features_from_entry(cache[rid], direction=direction)
        c = cost_relative_from_features(F, weights)
        total += c
        details.append({"reaction_id": rid, "direction": direction, "cost": c, **F})
    return total, details

