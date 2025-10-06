
# Pyruvate Valorization Module — SynPathTransfer

> **What this file is:** a single, self-contained GitHub README that precisely documents *what we built* for the Pyruvate Valorization Module: a KEGG-driven sub-pathway enumerator plus a relative energetic cost model to rank candidate routes starting from pyruvate (or any seed compound).

---

## What We Did
- Implemented sub-pathway enumeration over KEGG pathways by traversing reaction adjacencies derived from parsed reaction equations.
- Implemented a relative cost model that scores each reaction (and sums along a path) using ATP equivalents, redox burden (via P/O), O₂ usage, CO₂ release, a “complexity” proxy, and a small “precedent” prior.
- Selected the lowest-cost (most favorable) sub-pathway and visualized it on KEGG’s map UI with highlighted reactions and the seed compound.
- Wrote two Python modules:
  - `search.py` — graph exploration and visualization
  - `cost.py` — feature extraction from KEGG entries and linear cost scoring

## Why This Module Exists
The broader SynPathTransfer pipeline explores transplanting the 3-hydroxypropionate (3HP) cycle into cyanobacteria as a strategy to bypass RuBisCO inefficiencies. A common by-product is excess pyruvate (KEGG: `C00022`). Our module automatically discovers and ranks downstream routes that valorize pyruvate toward useful products (e.g., lactate, ethanol, butanol, isoprene, PHB), favoring routes that are energetically/redox-wise more attractive and operationally simpler.

---

## Repository Layout
```
synpathtransfer/
└── pyruvate_valorization/
    ├── cost.py        # relative cost model + KEGG feature extraction
    ├── search.py      # KEGG graph traversal, path enumeration, visualization
    ├── examples/
    │   └── demo.py    # minimal runnable example
    └── README.md      # this file
```
> Only `requests` is required at runtime; no heavy scientific stack is necessary.

---

## How It Works (Design Overview)

### Reaction-Centric Graph in a KEGG Map
- **Nodes:** KEGG reactions (RIDs like `R00200`).
- **Edges (Rᵢ → Rⱼ):** If **products(Rᵢ)** intersect **substrates(Rⱼ)** after removing common cofactors to avoid inflated connectivity.
- **Directionality:** From the equation arrow in KEGG entries:
  - `A => B` = irreversible (we record sign `-1` for direction penalties)
  - `A <=> B` = reversible (sign `+1`)

**Cofactors filtered by default:** `ATP, ADP, AMP, Pi, PPi, NAD+, NADH, NADP+, NADPH, H2O, H+, CO2, CoA, NHx`.

### Sub-Pathway Enumeration
1. Pull all reactions belonging to the selected KEGG map.
2. Preload their equations (lower API latency).
3. For each candidate start reaction (consuming the seed compound), depth-first traverse via reaction–reaction adjacencies.
4. Record each reachable reaction as `(RID, sign)`; gather all subpaths with length > 1.

### Relative Cost Model (Per-Reaction, Summed Over Path)
We extract features from `get/reaction:{RID}` and apply a linear score:

- **ΔATPeq** — triphosphate usage minus partial credit for ADP/AMP “recovery”
- **ΔREDOX_ATP** — NADH/NADPH/FADH₂ converted to ATP-equivalents using **P/O ≈ 2.5/2.5/1.5**
- **O₂ consumption** (≥ 0)
- **CO₂ release** (≥ 0)
- **Complexity** — count of unique non-trivial metabolites across both sides
- **Precedent** — `1/(1 + p)` where *p* = number of KEGG pathways listing this reaction

**Default weights:**
```
alpha=1.0   # ATP equivalents
beta=1.0    # redox equivalents (ATP-eq)
gamma=0.3   # O2 consumption penalty
delta=0.25  # CO2 release penalty
eps=0.20    # complexity proxy
zeta=1.0    # precedent prior
```

**Total path score** = sum of per-reaction costs along the subpath.

---

## Installation

```bash
# Create a fresh environment (optional but recommended)
python -m venv .venv && source .venv/bin/activate  # on Windows: .venv\Scripts\activate

# Install dependency
pip install requests
```

---

## Quickstart

```python
# examples/demo.py
from search import visualize_best_subpathway

# KEGG map for pyruvate-related metabolism (example)
map_id = "map00720"
seed_cpd = "C00022"  # pyruvate

visualize_best_subpathway(map_id, seed_cpd)
```
Run:
```bash
python examples/demo.py
```
This opens KEGG with the selected reaction chain highlighted and the seed compound marked.

---

## API (Minimal Surface)

### `cost.py`
- `reaction_cost_relative(rid: str, direction: int = +1, weights: dict = REL_W) -> (float, dict)`

  Compute cost and return `(score, feature_dict)` for a single KEGG reaction.
- `subpathway_cost_relative(steps: list[tuple[str,int]], weights: dict = REL_W) -> (float, list[dict])`

  Sum costs along a sequence `[(RID, dir), ...]` and return `(total, per_step_details)`.

### `search.py`
- `find_all_subpathways(pathway: str, compound: str) -> list[list[tuple[str,int]]]`

  Enumerate subpaths within a KEGG map starting from reactions that consume the seed compound.
- `best_pw(subpaths: list[list[tuple[str,int]]]) -> str`

  Return KEGG URL segment for the minimal-cost subpath (internal helper).
- `visualize_best_subpathway(pathway: str, compound: str) -> None`

  Open the best subpath overlay on KEGG’s interactive pathway map.

---

## Example: Programmatic Scoring

```python
from cost import subpathway_cost_relative

steps = [("R00200", +1), ("R00703", -1), ("R01015", +1)]  # example RIDs
total, details = subpathway_cost_relative(steps)
print("Total score:", total)
for d in details:
    print(d["reaction_id"], d["direction"], d["cost"], d["dATPeq"], d["dREDOX_ATP"])
```

---

## Notes, Assumptions & Limitations
- **Direction vs. physiology:** Arrow direction in KEGG equations may not match in vivo flux; our model uses the text arrows unless overridden.
- **Cofactor list:** Fixed list; extend if your map is cofactor-heavy.
- **Equation parsing edge cases:** A few entries can be missing/irregular; we default to empty sets and skip.
- **Graph simplification:** Reaction-centric adjacency can over-connect via hub metabolites even after cofactor filtering; tune filters as needed.
- **Scoring realism:** P/O conversions are approximate. Toxicity and explicit path-length terms are future work.
- **Rate limiting:** We cache reactions and throttle (sleep ~0.2 s) to be polite to KEGG; large maps will still take time.
- **Implementation footnote:** Initialize your best-path search with `min_cost = float('inf')` so positive-cost paths are considered.

---

## Extending the Module
- Add **toxicity** (lookup tables) and **explicit length penalty**.
- Integrate **ΔG′m** estimates for thermodynamic feasibility.
- Support **multi-seed** and **multi-objective** optimization.
- Move to **KGML** topology when pathway layout/compartment info is needed.
- Persist a **local cache** of KEGG entries for reproducibility.

---

## License & Attribution
- This module queries **KEGG REST** endpoints. Please respect KEGG’s terms of use and rate limits.
- © 2025 SynPathTransfer contributors. See project LICENSE for details.
