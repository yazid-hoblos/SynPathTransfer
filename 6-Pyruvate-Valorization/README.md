# Pyruvate Valorization Module — SynPathTransfer

> **What this file is:** a single, self-contained GitHub README that precisely documents what we built for the **Pyruvate Valorization Module**: a KEGG-driven sub-pathway enumerator plus a **relative energetic cost model** to rank candidate routes starting from pyruvate (or any seed compound).

---

## 1) What We Did (general outlook on the module)

We implemented a small, composable toolkit that:

- **Enumerates reaction sub-pathways** inside a KEGG pathway, starting from a given compound (e.g., pyruvate, `C00022`), using reaction equations from KEGG REST.
- **Scores each sub-pathway** with a **relative energetic cost** that approximates ATP usage, redox balance (via P/O ratios), O₂ consumption, CO₂ release, reaction “complexity,” and literature precedent.
- **Selects and visualizes** the best-scoring route on KEGG’s pathway map.

This is packaged as two Python modules—`search.py` (graph traversal over reactions) and `cost.py` (feature extraction and scoring)—plus a tiny demo snippet you can adapt.

---

## 2) Why This Module Exists

In SynPathTransfer, the 3-HP cycle can generate **excess pyruvate**. To turn that surplus into a useful product (lactate, ethanol, butanol, isoprene, PHB, etc.), we can explore **downstream routes** and rank them by plausible **bioenergetic/physiological cost**. Our module automates this exploration from KEGG, surfaces candidate sub-pathways, and provides a **relative** score to prioritize engineering options.

---

## 3) Repository Layout

```
synpathtransfer/
├─ pyruvate_valorization/
│  ├─ cost.py     # KEGG parsing + feature extraction + relative cost model
│  ├─ search.py   # Reaction-centric traversal and sub-pathway enumeration
│  └─ demo.py     # (optional) Example usage wiring the two together
└─ README.md      # this file
```

### a) `cost`
- Fetches KEGG reaction entries and parses **ENTRY**, **EQUATION**, and **PATHWAY** fields.
- Computes per-reaction features (ATP-eq, redox ATP-eq via P/O, O₂ use, CO₂ release, complexity, precedent).
- Aggregates per-step features into a **relative** path score.

### b) `demo`
- Shows typical calls to compute costs for a reaction and for a sub-pathway; and how to visualize the best route.

### c) `search`
- Builds a reaction-adjacency from KEGG equations (with simple cofactor filtering).
- Enumerates sub-pathways with a DFS-style traversal.
- Hands candidate paths to `cost` for ranking; opens the best in a browser on KEGG.

---

## 3) How Does It Work?

> Below is a function-by-function tour aligned with the code included in this module.

### a) `search.py` — sub-pathway enumeration

- **`find_reactions_from_compound(map: str, compound: str) -> dict`**  
  Uses KEGG REST to intersect:
  - all reactions in a **pathway** (`/link/reaction/{map}`), with
  - all reactions **involving** the **compound** (`/link/reaction/{compound}`),  
  then filters out cases where the compound only appears on the **product** side (we only keep reactions where it’s a substrate/start).

- **`clean_equation(eq: str) -> str`**  
  Removes common **cofactors** (ATP/ADP/AMP, Pi/PPi, NAD(P)H, H₂O, H⁺, CO₂, CoA, NH, etc.) to reduce false connectivity via hub cofactors.

- **`parse_equation(equation: str) -> (set, set, int)`**  
  Splits an equation into **substrate set**, **product set**, and **sign**: `-1` for irreversible (`=>`) and `+1` for reversible (`<=>`).

- **`load_equations(all_reactions: list) -> dict`**  
  Bulk-loads equations for all reactions via `/list/reaction:{RID}` and returns a mapping:
  ```python
  equations[RID] = (substrates_set, products_set, sign)
  ```
  Missing/irregular entries degrade gracefully to empty sets.

- **`find_subpathways_from_reac(start_reac, all_reactions, equations) -> list`**  
  Reaction-centric DFS using a **stack**. Two reactions are considered adjacent if **products** of one overlap **substrates** of the other (and vice-versa), after cofactor filtering. Returns a “pathway” list:
  ```python
  [(reaction_id, sign), ...]
  ```

- **`find_all_subpathways(pathway: str, compound: str) -> list[list[tuple]]`**  
  For every reaction in the pathway graph, launches the local traversal to capture subpaths. Keeps any with length > 1.

- **`best_pw(subpaths) -> str`**  
  Calls into `cost.subpathway_cost_relative(...)` for each candidate, picks the **minimum** total cost, and returns a **slash-joined** list of reaction IDs for KEGG visualization.  
  > **Important:** initialize `min_cost = float('inf')` so positive-cost paths are considered (see Notes §4).

- **`visualize_best_subpathway(pathway: str, compound: str)`**  
  Computes the best path and opens:
  ```
  https://www.kegg.jp/kegg-bin/show_pathway?{pathway}/{R1}/{R2}/.../{compound}%20%23ff0000
  ```
  highlighting the selected reactions and the seed compound (red).

### b) `cost.py` — KEGG parsing & relative cost model

- **`kegg_get(rid) -> str`**  
  Downloads a KEGG **reaction** entry (`/get/{rid}`).

- **`parse_entry(raw) -> dict`**  
  Extracts `ENTRY` (RID), `EQUATION` (raw text), and `PATHWAY` (list of maps). We use `PATHWAY` count to model a crude **literature precedent** term.

- **`parse_side(side_text) -> dict[Cxxxx: stoich]`** & **`parse_equation(eq)`**  
  Tokenize LHS/RHS, handle optional coefficients, return two dicts `{Cxxxx: stoich}`. Supports both `<=>` and `=>`.

- **`features_from_entry(entry_text, direction=+1) -> dict`**  
  Computes per-reaction features used by the score:
  - **ΔATPeq** — Net ATP-equivalent triphosphate hydrolysis (ATP/GTP/UTP/CTP minus partial recovery of ADP/AMP/GDP); we include it to penalize reactions that directly drain cellular energy.
  - **ΔREDOX_ATP** — Net change in NAD(P)H/FADH₂ converted to ATP-equivalents via P/O ratios; we use it to capture the opportunity cost or benefit of consuming or generating reducing power.
  - **O₂ consumption** — Moles of molecular oxygen required by the reaction; we track it because O₂ dependence can be limiting in some hosts and adds physiological burden.
  - **CO₂ release** — Moles of carbon released as CO₂; we penalize it because decarboxylation lowers carbon yield from pyruvate.
  - **Complexity** — Count of distinct non-trivial metabolites on either side of the equation; we use it as a proxy for implementation burden (more parts, transport, regulation).
  - **Precedent** — `1/(1 + pcount)` where `pcount` is the number of KEGG pathway maps containing the reaction; this favors well-charted, commonly used steps and slightly penalizes rare ones.


- **`cost_relative_from_features(F, W)`**  
  Linear combination:
  ```
  cost = α·ΔATPeq + β·ΔREDOX_ATP + γ·O2 + δ·CO2 + ε·complexity + ζ·precedent
  ```
  with default weights that can be tweaked:
  ```
  α=1.0, β=1.0, γ=0.3, δ=0.25, ε=0.20, ζ=1.0
  ```

- **`reaction_cost_relative(rid, direction=+1, weights=REL_W)`**  
  Returns `(score, features)` for a single reaction.

- **`subpathway_cost_relative(steps, weights=REL_W)`**  
  For `steps = [(RID, direction), ...]`, fetches each reaction once (polite **0.2 s sleep**) and sums per-step costs. Returns:
  ```python
  total_cost, [
    {"reaction_id": RID, "direction": ±1, "cost": c, **features}, ...
  ]
  ```

### c) `demo.py` — minimal usage

```python
# demo.py
from search import visualize_best_subpathway, find_all_subpathways, best_pw
from cost import subpathway_cost_relative, reaction_cost_relative

# Example: start from pyruvate (C00022) inside a KEGG pathway map, e.g. map00720
PATHWAY = "map00720"
PYRUVATE = "C00022"

# 1) Visualize best route directly on KEGG
visualize_best_subpathway(PATHWAY, PYRUVATE)

# 2) Or enumerate & score yourself
subs = find_all_subpathways(PATHWAY, PYRUVATE)
best_reac_string = best_pw(subs)  # "R00001/R01234/..."
print("Best subpathway (KEGG string):", best_reac_string)

# 3) Inspect detailed costs for one candidate subpath
if subs:
    total, details = subpathway_cost_relative(subs[0])
    print("Total relative cost:", total)
    for step in details:
        print(step["reaction_id"], step["cost"], step["dATPeq"], step["dREDOX_ATP"])
```

---

## 4) Notes, Assumptions & Limitations

- **Direction vs. physiology:** Arrow direction in KEGG equations may not match in vivo flux; our model uses the text arrows unless overridden.  
- **Cofactor list:** Fixed list; extend if your map is cofactor-heavy.  
- **Equation parsing edge cases:** A few entries can be missing/irregular; we default to empty sets and skip.  
- **Graph simplification:** Reaction-centric adjacency can over-connect via hub metabolites even after cofactor filtering; tune filters as needed.  
- **Scoring realism:** P/O conversions are approximate. Toxicity and explicit path-length terms are future work.  
- **Rate limiting:** We cache reactions and throttle (sleep ~0.2 s) to be polite to KEGG; large maps will still take time.  
- **Implementation footnote:** Initialize your best-path search with `min_cost = float('inf')` so positive-cost paths are considered.

---

## 5) Extending the Module

- **Add toxicity** (lookup tables) and **explicit length penalty**.  
- Integrate **ΔG′m** estimates for thermodynamic feasibility.  
- Support **multi-seed** and **multi-objective** optimization.  
- Move to **KGML topology** when pathway layout/compartment info is needed.  
- Persist a **local cache** of KEGG entries for reproducibility.

---

## 6) License & Attribution

This module queries **KEGG REST** endpoints. Please respect KEGG’s terms of use and rate limits.  
© 2025 SynPathTransfer contributors. See project LICENSE for details.
