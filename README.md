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
| **Estrogenic activity** | A registry of real phytoestrogens (genistein, coumestrol, …) plus a constructor for *hypothetical* molecules, each characterized by relative estrogenic activity (REA) vs estradiol. Folds any compound's plasma curve into **estradiol‑equivalents (EE)** so breast‑growth creams, soy supplements, and estradiol sit on one scale. |
| **Inverse design** | Given a target EE (e.g. the 100–200 pg/mL window), what dose does a candidate/hypothetical molecule need? "Care yet to be discovered" explored as a property space. |
| **Breast development** | Cumulative‑EE‑driven saturating‑growth model mapped to Tanner stages (population‑typical; genetics dominate the ceiling). |
| **Galactagogues & lactation** | D2‑antagonist registry (metoclopramide, domperidone, sulpiride, amisulpride) with PK + prolactin Emax/ED50; Newman–Goldfarb protocol with **two‑compartment priming** (E2‑driven ductal scaffold + P‑driven lobuloalveolar units), **folded‑in progesterone PK** (micronized P's ~3× peak/mean swing; LA integrates exposure), **P‑withdrawal lactogenesis gate**, **E2 taper vs hard drop**, and the **estrogen brake** on maintained HRT. End‑to‑end bands (drug *and* gland uncertainty); equieffective `compare`; low‑dose `sweep`; `risk` exposure framing; **`personalize`** — importance‑sampling inference from a participant's own prolactin/supply measurements, with a posterior‑predictive titration table. **No controlled trials in trans women** — calibrated from postpartum cis literature, large uncertainty flagged. |

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
hormone-sim pop --route transdermal --dose 50 --n 200 --weeks 12 --out outputs/

# Estrogenic activity: how estrogenic is genistein? EE at 500 pg/mL plasma?
hormone-sim compound
hormone-sim compound --name genistein --plasma-pgml 500

# Inverse design: what dose does a hypothetical molecule (REA=0.5) need?
hormone-sim hypothetical --name future-drug --rea 0.5 --target-ee 150

# Breast development: Tanner stage over 2 years on transdermal 50 µg/day
hormone-sim breast --regimen configs/tf_transdermal.yaml --weeks 104

# Galactagogues: registry + prolactin dose-response
hormone-sim galactagogues
hormone-sim galactagogues --name domperidone

# Lactation induction: Newman-Goldfarb protocol, 7 mo prep, domperidone 40 mg/day
hormone-sim lactation --galactagogue domperidone --dose 40 \
    --prep-weeks 30 --sessions 8 --out outputs/

# Same, with Monte-Carlo sensitivity bands on the milk parameters
hormone-sim lactation --galactagogue domperidone --dose 40 --bands --n 300

# Interrogate progesterone timing: P last 20 wk vs P-free protocol
hormone-sim lactation --prep-weeks 30 --p-start-weeks 10
hormone-sim lactation --prep-weeks 30 --p-start-weeks 30   # P-free

# Head-to-head: agents dose-matched to the SAME prolactin target
hormone-sim compare --agents domperidone metoclopramide --target-prolactin 130

# Risk-exposure framing for a planned course (duration/age/cumulative tiers)
hormone-sim risk --agent metoclopramide --dose 15 --weeks 12 --age 30

# Personal inference from a participant's measurements (study layer):
# week 2 post-drop on 10 mg/d, PRL 72, pumping yield 95 mL/day, etc.
hormone-sim personalize --agent metoclopramide --dose 10 \
    --obs '2:prl=72,milk=95' --obs '6:prl=78,milk=140' \
    --candidate-doses 10 15 20

# Study tools (for clinician-supervised titration studies):
hormone-sim ingest --csv examples/participant_log.csv --agent metoclopramide \
    --milk-unit oz --out outputs/personal_posterior.png   # log -> posterior + figure
hormone-sim aims --baseline '0:all=0' --visit '8:lips_perioral=2'  # AIMS TD monitoring
hormone-sim protocol --agent metoclopramide --out docs/protocol-draft.md  # IRB skeleton
```

The example participant log lives in [`examples/`](examples/).

## Trans-masculine / masculinization domain (v0.8)

The mirror domain is now built: testosterone PK (`androgens.py` — injectable
esters with the peak–trough story quantified, gel for stable levels; oral
undecanoate deliberately unregistered), reversed-direction HPG suppression
(ovarian output falls as T rises), masculinization outcomes
(`masculinization.py` — cycle-cessation anchor ~5 months; voice/hair/body/bottom
saturating clocks with genetics-dominant ceilings), and erythrocytosis risk
framing (`erythrocytosis.py` — the case where the risk IS lab-monitorable,
dose/route-dependent, and reversible: the mirror image of the TD asymmetry).

```bash
hormone-sim androgens
hormone-sim masc --agent testosterone_cypionate --dose 100 --interval-days 7 \
    --weeks 104 --n 150 --out outputs/   # outcomes + bands + figure
hormone-sim hct --agent testosterone_cypionate --dose 100 --interval-days 7
```

Study-layer generalization (ingest/personalize on T logs) is the next open
item — see [`docs/roadmap-trans-masc.md`](docs/roadmap-trans-masc.md).

## Intersex condition packs (v0.9)

Parameter overlays that adapt the models to documented intersex endocrinology:
CAIS, PAIS (framework-only — individual receptor sensitivity is not
predictable from population data and the pack says so), CAH (21-hydroxylase),
Klinefelter, and Turner. Each pack changes the endogenous baselines, carries
citations and an anchor-quality rating, and states its limits in plain
language. The module docstring names why this exists: intersex people —
including intersex trans people — have historically been denied both bodily
autonomy (non-consensual normalization procedures) and information about
their own endocrinology. Legibility is the point.

```bash
hormone-sim intersex --name cah_21ohd           # a pack, with its limits
hormone-sim masc --pack cah_21ohd --dose 100 --interval-days 7
hormone-sim ingest-t --csv examples/t_log.csv --pack cah_21ohd
```

The T-log study layer now mirrors the lactation one: `ingest-t` reads a
trans-masculine log (serum T + CBCs, mixed date/weeks axes, fail-closed
validation) and runs importance-sampling personalization with a posterior
predictive over candidate (dose, interval) regimens — trough/mean/peak plus
Hct plateau, because for T the schedule changes the *swing*, not just the
level. See [`examples/t_log.csv`](examples/t_log.csv).

Programmatic use:

```python
from hormone_sim import Regimen, simulate
from hormone_sim.plot import plot_single

r = Regimen.from_yaml("configs/tf_transdermal.yaml")
sim = simulate(r, weeks=12)          # returns a SimResult
plot_single(sim, out_path="outputs/demo.png")   # E2 + T vs time, with target bands
print(sim.steady_state_e2_pgml())
```

Estrogenic-activity example:

```python
from hormone_sim import get_compound, ee_of
import numpy as np

gen = get_compound("genistein")
# EE at a plasma concentration of 500 pg/mL:
ee = ee_of(gen, np.array([500.0])).ee_curve_pgml[0]
print(f"genistein at 500 pg/mL -> {ee:.3f} pg/mL estradiol-equivalent")
```
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
