# hormone-sim

**Population pharmacokinetic / pharmacodynamic (PK–PD) simulation of gender‑affirming
hormone therapy (GAHT).**

A reproducible, literature‑parameterized modeling toolkit for the *endocrinology of
cross‑sex hormone therapy* — how serum **estradiol (E2)** and **testosterone (T)**
concentrations evolve over time under different administration routes, doses, and
anti‑androgen regimens. The intent is research, education, and study *design* — the kind
of pre‑clinical modeling that precedes and informs properly supervised clinical work.

> ⚠️ **This is a population model, not a patient.** It produces *population‑typical*
> concentration curves calibrated to published steady‑state values and guideline target
> ranges. It does **not** produce individualized dosing. Real GAHT is titrated by a
> qualified clinician against *measured* serum levels, renal/hepatic function, and
> cardiovascular risk. Nothing here should be used to self‑medicate or to dose a real
> person. See [`docs/governance.md`](docs/governance.md).

---

## What it models

| Component | Content |
|---|---|
| **PK — estradiol** | 1‑compartment model with first‑order (or zero‑order, for patches) absorption and first‑order elimination. Route differences live entirely in the *input function* and bioavailability; clearance and volume are shared molecule properties. |
| **Routes** | Oral (micronized), transdermal (patch/gel, zero‑order), intramuscular & subcutaneous esters (depot, slow absorption). |
| **PD — testosterone** | Suppression of the hypothalamic‑pituitary‑gonadal (HPG) axis by estradiol negative feedback, plus additive anti‑androgen effects (spironolactone, cyproterone acetate, GnRH agonists, bicalutamide). |
| **Targets** | Endocrine Society / WPATH ranges: E2 ≈ 100–200 pg/mL; T < 50 ng/dL (female range). |
| **Population variability** | Monte‑Carlo sampling of clearance, volume, and sensitivity to produce confidence bands, not single lines. |

The model is deliberately **simplified and transparent**: every parameter is a named,
citable constant (see [`docs/references.md`](docs/references.md)), and the few that are
*calibrated* rather than directly measured are flagged in [`src/hormone_sim/params.py`](src/hormone_sim/params.py).

---

## Install

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

Requires Python ≥ 3.10, `numpy`, `scipy`, `matplotlib` (and `pytest` to run tests).

## Quick start

```bash
# Simulate a transdermal 50 µg/day patch for 12 weeks, plot E2 + T
hormone-sim run --regimen configs/tf_transdermal.yaml --weeks 12 --out outputs/

# Dose–response: what steady-state E2 does each transdermal dose give?
hormone-sim scan --route transdermal --doses 25 50 75 100 --weeks 16 --out outputs/

# Population bands: 200 subjects, transdermal 50 µg/day
hormone-sim pop --regimen configs/tf_transdermal.yaml --n 200 --weeks 12 --out outputs/
```

Programmatic use:

```python
from hormone_sim import Regimen, simulate
import pandas as pd

r = Regimen.from_yaml("configs/tf_transdermal.yaml")
sim = simulate(r, weeks=12)          # returns a SimResult
sim.plot()                            # E2 + T vs time, with target bands
print(sim.steady_state_e2_pg_ml())
```

## Project layout

```
hormone-sim/
├── src/hormone_sim/     # the package (pk, pd, params, regimens, population, sim, cli)
├── configs/             # YAML regimen definitions (editable, reproducible inputs)
├── scripts/             # dose-scan / study-design drivers
├── tests/               # unit + integration tests (assert physiological sanity)
├── docs/                # references, limitations, governance/ethics
└── outputs/             # generated plots + CSVs (gitignored)
```

## Reproducibility

- Regimens are **YAML files** — the input is the artifact; the simulation is deterministic
  given a seed (`--seed`).
- All parameters carry a citation and a `calibrated: true/false` flag.
- `scripts/run_scans.py` reproduces the study figures from the README.
- Run `pytest -q` before claiming results; the numerical invariants are asserted in
  [`tests/`](tests/).

## Citation & licensing

Modeling framework: [MIT License](LICENSE). Clinical target ranges and parameter sources:
[docs/references.md](docs/references.md) — cite the primary guidelines when you publish
work that builds on this.
