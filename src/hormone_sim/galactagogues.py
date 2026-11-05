"""Galactagogue pharmacology: D2-antagonist lactation-induction agents.

Prolactin secretion is under tonic inhibition by tuberoinfundibular dopamine
(TIDA) neurons acting on D2 receptors on lactotrophs. D2 antagonists remove
that brake, raising prolactin. This module models the **pharmacology of the
agents** — their PK (one-compartment, repeated oral dosing) and their
**prolactin dose-response** (Emax/ED50 model on daily dose).

Registry covers the agents actually used for lactation induction
(pharmacologic galactagogues). Herbal galactagogues (fenugreek, blessed thistle,
shatavari) are NOT here — they have no reliable D2-antagonist mechanism and no
quantifiable potency; see ``docs/references.md`` for why they are excluded
rather than hand-waved in.

Evidence caveat (important): the PK and prolactin dose-response values below
come from **postpartum cis women and medication-induced hyperprolactinemia**
literature. There are no controlled trials in trans women. The model
extrapolates to a trans-feminine physiology (on estradiol, no pregnancy
priming) and is explicit about that.

Safety flags carried in the registry:
  * metoclopramide — boxed warning: tardive dyskinesia; also EPS, depression.
  * domperidone    — not FDA-approved (QT prolongation / arrhythmia risk);
                     approved in EU/Canada/UK for other indications.
  * sulpiride/amisulpride — antipsychotics; off-label galactagogues in some
                     countries; significant prolactin elevation.

This is a pharmacology model, not medical advice. Real lactation induction is
run with a clinician monitoring prolactin, cardiac status (domperidone), and
neuropsychiatric side effects (metoclopramide).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np


# --- Prolactin baseline (population-typical) ---------------------------------
# Adult AMAB baseline ~2-18 ng/mL (mean ~8). Estradiol is a lactotroph
# secretagogue: on trans-feminine HRT (E2 ~100-200 pg/mL) baseline runs
# modestly higher, roughly 1.2-1.5x. Use ~12 ng/mL as the on-HRT baseline.
PROLACTIN_BASELINE_NG_ML = 12.0
# Pregnancy (term) reaches ~200-400 ng/mL; full lactation ~100-300 ng/mL.
PROLACTIN_LACTATION_REFERENCE_NG_ML = 150.0


@dataclass(frozen=True)
class Galactagogue:
    """A D2-antagonist galactagogue's PK and prolactin dose-response.

    PK (one-compartment, oral, repeated dosing):
      * ``half_life_h``      — apparent elimination half-life (h)
      * ``bioavailability``  — fraction of oral dose reaching systemic circulation
      * ``vd_l_kg``          — volume of distribution (L/kg)

    Prolactin effect (Emax/ED50 on **daily dose**, for transparency):
      * ``prolactin_emax_ng_ml`` — maximum prolactin *rise* above baseline at
        saturating dose (ng/mL).
      * ``prolactin_ed50_mg_day`` — daily oral dose achieving 50% of Emax.

    Dose-response is given in dose (mg/day) rather than concentration because
    that is what a protocol specifies and what is measurable. The PK layer
    converts dose -> plasma level -> time-varying effect.
    """

    name: str
    half_life_h: float
    bioavailability: float
    vd_l_kg: float
    prolactin_emax_ng_ml: float
    prolactin_ed50_mg_day: float
    citation: str
    safety_note: str = ""
    approved_us: bool = True   # domperidone is not FDA-approved

    # --- PK helpers ---------------------------------------------------------
    def clearance_l_per_h(self, body_weight_kg: float = 70.0) -> float:
        """Implied clearance (L/h) = ln2 * Vd / t1/2."""
        return (0.6931471805599453 * self.vd_l_kg * body_weight_kg) / self.half_life_h

    def peak_conc_ng_ml(self, dose_mg: float, body_weight_kg: float = 70.0) -> float:
        """Single-dose peak plasma concentration (ng/mL) — bolus approximation."""
        vd_ml = self.vd_l_kg * body_weight_kg * 1000.0
        return (dose_mg * 1e6 * self.bioavailability) / vd_ml


# --- Registry -----------------------------------------------------------------
# PK values are population-typical literature estimates (see references.md);
# prolactin Emax/ED50 are calibrated to observed responses at standard doses.
GALACTAGOGUES: Dict[str, Galactagogue] = {
    "metoclopramide": Galactagogue(
        name="metoclopramide",
        half_life_h=5.5,
        bioavailability=0.80,
        vd_l_kg=2.5,
        prolactin_emax_ng_ml=150.0,
        prolactin_ed50_mg_day=25.0,
        citation="Metoclopramide PK: t1/2 ~5-6 h, F ~80%, Vd ~2.5 L/kg "
        "(StatPearls; Sankar et al. 2017 review of galactagogues). "
        "Prolactin response calibrated to 10 mg TID/QID -> ~30-80 ng/mL rise.",
        safety_note="Boxed warning: tardive dyskinesia (often irreversible). "
        "Also EPS, depression. Shortest safe course.",
        approved_us=True,
    ),
    "domperidone": Galactagogue(
        name="domperidone",
        half_life_h=15.0,
        bioavailability=0.15,
        vd_l_kg=3.0,
        prolactin_emax_ng_ml=200.0,
        prolactin_ed50_mg_day=30.0,
        citation="Domperidone PK: t1/2 ~14-16 h, F ~13-17% (hepatic first-pass + "
        "CYP3A4 degradation), Vd ~3 L/kg. Prolactin response calibrated to "
        "10-20 mg TID/QID -> ~50-150 ng/mL rise (Newman-Goldfarb protocol "
        "literature; OSHERC 2018 review).",
        safety_note="NOT FDA-approved (QT prolongation / ventricular "
        "arrhythmia risk; CYP3A4 interactions). Approved in EU/Canada/UK for "
        "other indications. Baseline ECG and potassium are recommended.",
        approved_us=False,
    ),
    "sulpiride": Galactagogue(
        name="sulpiride",
        half_life_h=7.5,
        bioavailability=0.35,
        vd_l_kg=2.7,
        prolactin_emax_ng_ml=250.0,
        prolactin_ed50_mg_day=75.0,
        citation="Sulpiride PK: t1/2 ~7-8 h, F ~35%, Vd ~2.7 L/kg. Antipsychotic "
        "with potent prolactin elevation; used as a galactagogue in some "
        "countries (off-label).",
        safety_note="Antipsychotic: sedation, EPS, weight gain. Not available "
        "in the US. Used off-label elsewhere.",
        approved_us=False,
    ),
    "amisulpride": Galactagogue(
        name="amisulpride",
        half_life_h=12.0,
        bioavailability=0.48,
        vd_l_kg=5.3,
        prolactin_emax_ng_ml=280.0,
        prolactin_ed50_mg_day=100.0,
        citation="Amisulpride PK: t1/2 ~12 h, F ~48%, Vd ~5.3 L/kg. Potent "
        "hyperprolactinemic agent; off-label galactagogue use in some countries.",
        safety_note="Antipsychotic: EPS, sedation, QT prolongation at higher "
        "doses. Not approved for lactation.",
        approved_us=False,
    ),
}


def get_galactagogue(name: str) -> Galactagogue:
    """Look up a galactagogue by name (raises KeyError if unknown)."""
    if name not in GALACTAGOGUES:
        raise KeyError(f"Unknown galactagogue {name!r}. Known: {sorted(GALACTAGOGUES)}")
    return GALACTAGOGUES[name]


# --- PK: repeated oral dosing, one-compartment --------------------------------
def simulate_dosing_pk(
    g: Galactagogue,
    dose_mg: float,
    dosing_interval_h: float = 8.0,
    duration_h: float = 24.0 * 30.0,
    dt_h: float = 0.1,
    body_weight_kg: float = 70.0,
) -> np.ndarray:
    """Simulate plasma concentration (ng/mL) under repeated oral dosing.

    One-compartment, first-order elimination, bolus doses every
    ``dosing_interval_h`` (oral F applied at each dose). Returns the
    concentration curve over ``duration_h``.
    """
    k_el = np.log(2) / g.half_life_h
    vd_ml = g.vd_l_kg * body_weight_kg * 1000.0

    n = int(round(duration_h / dt_h)) + 1
    A = np.zeros(n)   # amount in body, ng

    # Dose events every dosing_interval_h (first dose at t=0).
    dose_interval_steps = max(1, int(round(dosing_interval_h / dt_h)))
    dose_ug = dose_mg * 1e6 * g.bioavailability  # ng

    for i in range(1, n):
        A[i] = A[i - 1] * np.exp(-k_el * dt_h)
        if (i - 1) % dose_interval_steps == 0:
            A[i] += dose_ug

    return A / vd_ml


# --- Prolactin pharmacodynamics ------------------------------------------------
def prolactin_response(
    g: Galactagogue,
    dose_mg_day: float,
) -> float:
    """Steady-state prolactin (ng/mL) at a given daily dose (Emax/ED50 model).

    P_ss = baseline + Emax * dose / (ED50 + dose)

    This is the *dose-based* steady-state estimate (the PK layer gives the
    time-varying version). Baseline is the trans-feminine on-HRT value.
    """
    d = max(0.0, float(dose_mg_day))
    rise = g.prolactin_emax_ng_ml * (d / (g.prolactin_ed50_mg_day + d))
    return PROLACTIN_BASELINE_NG_ML + rise


def prolactin_curve(
    g: Galactagogue,
    dose_mg_day: float,
    dosing_interval_h: float = 8.0,
    duration_h: float = 24.0 * 30.0,
    dt_h: float = 0.1,
    body_weight_kg: float = 70.0,
) -> np.ndarray:
    """Time-varying prolactin (ng/mL) from the drug's PK curve.

    Uses a rapid-equilibrium approximation: prolactin (t1/2 ~20-30 min) tracks
    the drug concentration quasi-instantaneously, so the prolactin curve is the
    dose-response evaluated at each drug concentration. The daily dose is
    distributed across the dosing interval to get the concentration.
    """
    per_dose_mg = dose_mg_day * dosing_interval_h / 24.0
    conc = simulate_dosing_pk(
        g, per_dose_mg, dosing_interval_h=dosing_interval_h,
        duration_h=duration_h, dt_h=dt_h, body_weight_kg=body_weight_kg,
    )
    # Map concentration -> prolactin rise via a Hill curve anchored so that the
    # time-averaged concentration reproduces the dose-based steady state.
    avg_conc = float(conc[conc.size // 2:].mean())  # steady-state average
    if avg_conc <= 0:
        return np.full_like(conc, PROLACTIN_BASELINE_NG_ML)
    # EC50 in concentration units: the average concentration that gives 50% Emax
    # given the dose->concentration relationship.
    ss_prolactin = prolactin_response(g, dose_mg_day)
    target_rise = ss_prolactin - PROLACTIN_BASELINE_NG_ML
    if target_rise >= g.prolactin_emax_ng_ml:
        # Saturating; cap at Emax.
        return np.full_like(conc, PROLACTIN_BASELINE_NG_ML + g.prolactin_emax_ng_ml)
    # Rise/concentration relation: rise = Emax * C/(EC50 + C) -> EC50 from anchor.
    # At C = avg_conc, rise = target_rise.
    # target_rise/Emax = avg_conc/(EC50 + avg_conc) -> solve for EC50.
    frac = target_rise / g.prolactin_emax_ng_ml
    if frac <= 0 or frac >= 1:
        return np.full_like(conc, PROLACTIN_BASELINE_NG_ML + g.prolactin_emax_ng_ml)
    ec50_conc = avg_conc * (1.0 - frac) / frac
    if ec50_conc <= 0:
        return np.full_like(conc, PROLACTIN_BASELINE_NG_ML)
    rise = g.prolactin_emax_ng_ml * (conc / (ec50_conc + conc))
    return PROLACTIN_BASELINE_NG_ML + rise


def list_galactagogues() -> str:
    """Human-readable registry summary."""
    lines = [
        "Pharmacologic galactagogues (D2 antagonists):",
        f"{'name':>15} | {'t1/2 h':>6} | {'F':>5} | {'Vd L/kg':>7} | "
        f"{'Emax ng/mL':>10} | {'ED50 mg/d':>9} | {'US-approved':>11}",
        "-" * 84,
    ]
    for name in sorted(GALACTAGOGUES):
        g = GALACTAGOGUES[name]
        lines.append(
            f"{name:>15} | {g.half_life_h:>6.1f} | {g.bioavailability:>5.2f} | "
            f"{g.vd_l_kg:>7.1f} | {g.prolactin_emax_ng_ml:>10.0f} | "
            f"{g.prolactin_ed50_mg_day:>9.0f} | {str(g.approved_us):>11}"
        )
    lines.append("")
    lines.append("Safety: " + "; ".join(
        f"{g.name}: {g.safety_note}" for g in GALACTAGOGUES.values()
    ))
    return "\n".join(lines)
