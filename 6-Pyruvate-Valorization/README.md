
# Pyruvate Valorization Module — SynPathTransfer

> **What this file is:** a single, self-contained GitHub README that precisely documents *what we built* for the Pyruvate Valorization Module: a KEGG-driven sub-pathway enumerator plus a relative energetic cost model to rank candidate routes starting from pyruvate (or any seed compound).

---

## TL;DR (What We Did)

- Implemented **sub-pathway enumeration** over KEGG pathways by traversing reaction adjacencies derived from parsed reaction equations.
- Implemented a **relative cost model** that scores each reaction (and sums along a path) using ATP equivalents, redox burden (via P/O), O₂ usage, CO₂ release, a “complexity” proxy, and a small “precedent” prior.
- **Selected** the **lowest-cost** (most favorable) sub-pathway and **visualized** it on KEGG’s map UI with highlighted reactions and the seed compound.
- Wrote two Python modules:
  - `search.py` — graph exploration and visualization
  - `cost.py` — feature extraction from KEGG entries and linear cost scoring

---

## Why This Module Exists

The broader **SynPathTransfer** pipeline explores transplanting the **3-hydroxypropionate (3HP) cycle** into cyanobacteria as a strategy to bypass RuBisCO inefficiencies. A common by-product is excess **pyruvate (KEGG: C00022)**. Our module **automatically discovers and ranks downstream routes** that valorize pyruvate toward useful products (e.g., lactate, ethanol, butanol, isoprene, PHB), favoring routes that are energetically/redox-wise more attractive and operationally simpler.

---

## Repository Layout

