"""Breast development model — a cumulative-exposure PD outcome.

Breast growth on feminizing hormone therapy is driven primarily by **cumulative
estrogenic exposure** (estradiol-equivalents integrated over time), modulated by
genetics, age-at-start, and body fat. This module models glandular tissue volume
as a saturating function of cumulative EE exposure, mapped to a Tanner-stage-like
scale.

Model (phenomenological, population-typical):
    dV/dt = k_growth * EE(t) * (1 - V/V_max)
    V(0) = V_baseline
where V is a normalized glandular-volume index (0–1), EE is the
estradiol-equivalent exposure, ``k_growth`` sets the rate, and ``V_max`` the
genetic ceiling. The output is mapped to approximate Tanner stages.

This is deliberately a *population-typical* estimate. Real breast development
varies enormously between individuals (genetics dominate the ceiling), depends on
age and adiposity, and is not predictable from hormone levels alone. It is a
modeling/education tool, not a prediction for a specific person.

Epistemic status: breast-tissue pharmacodynamics is far less characterized than
plasma PK. Treat k_growth/V_max as calibrated-with-large-uncertainty.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .params import TARGET_E2_PGML


@dataclass(frozen=True)
class BreastModel:
    """Breast development parameters (calibrated, population-typical).

    ``k_growth``      — growth rate constant per pg/mL-EE per week (saturating).
                        Calibrated so ~2 yr at EE=150 pg/mL reaches ~Tanner 3
                        (V~0.45), with slower development at lower exposure.
    ``v_max``         — genetic ceiling of normalized glandular volume (0–1).
                        Population-typical ~0.75; reflects that most trans women
                        do not reach full Tanner 5.
    ``v_baseline``    — starting normalized volume (pre-HRT, pre-puberty index ~0).
    ``weeks_to_tanner`` — mapping from V to approximate Tanner stage.
    """

    k_growth_per_ee_week: float = 7e-5
    v_max: float = 0.75
    v_baseline: float = 0.0

    k_growth_source = "Calibrated (large uncertainty): population-typical breast " \
        "growth on HRT reaches Tanner 2–4 over 1–3 yr; see WPATH SOC-8 (2022), " \
        "Grimley et al. (2017). NOT individualized."

    def tanner_stage(self, v: float) -> int:
        """Map normalized volume index -> approximate Tanner stage (1–5)."""
        # Map V (0..v_max) to Tanner. V_max=0.75 -> Tanner ~4.
        frac = max(0.0, min(1.0, v / self.v_max))
        if frac < 0.05:
            return 1
        if frac < 0.25:
            return 2
        if frac < 0.50:
            return 3
        if frac < 0.80:
            return 4
        return 5


@dataclass
class BreastResult:
    """Breast development outcome over time."""

    weeks: np.ndarray
    volume: np.ndarray        # normalized glandular volume index (0..v_max)
    tanner: np.ndarray        # approximate Tanner stage per time point
    ee_curve: np.ndarray      # estradiol-equivalent exposure (pg/mL)

    def final_tanner(self) -> int:
        return int(self.tanner[-1])

    def final_volume_frac(self) -> float:
        """Fraction of genetic ceiling reached (0–1)."""
        return float(self.volume[-1] / BreastModel().v_max)

    def report(self) -> str:
        m = BreastModel()
        return (
            f"Breast development (population-typical, EE-driven):\n"
            f"  Final normalized volume: {self.volume[-1]:.2f} "
            f"({self.final_volume_frac()*100:.0f}% of genetic ceiling)\n"
            f"  Approximate final Tanner stage: {self.final_tanner()}\n"
            f"  NOTE: genetics dominate the ceiling; this is a population estimate, "
            f"not an individual prediction."
        )


def simulate_breast_development(
    ee_curve_pgml: np.ndarray,
    weeks: float,
    model: Optional[BreastModel] = None,
    dt_weeks: float = 0.25,
) -> BreastResult:
    """Integrate breast tissue growth driven by an EE exposure curve.

    ``ee_curve_pgml`` is the estradiol-equivalent curve (pg/mL) over the same
    time axis as ``weeks``. Uses explicit Euler in weeks.
    """
    model = model or BreastModel()
    ee = np.asarray(ee_curve_pgml, dtype=float)

    # Build a week-indexed EE series (interpolate onto weekly grid).
    n_weeks = int(round(weeks / dt_weeks)) + 1
    # Map ee (hourly) onto weekly grid: assume ee indexed by weeks already if len matches.
    if ee.size == n_weeks:
        ee_weekly = ee
    else:
        # ee is on an arbitrary grid; sample weekly by index proportional to time.
        ee_weekly = np.array([
            ee[min(ee.size - 1, int(round(w / dt_weeks * (ee.size - 1) / (weeks / dt_weeks))))]
            for w in range(n_weeks)
        ])

    volume = np.empty(n_weeks)
    tanner = np.empty(n_weeks, dtype=int)
    volume[0] = model.v_baseline
    tanner[0] = model.tanner_stage(volume[0])

    for i in range(1, n_weeks):
        # growth proportional to EE and to remaining capacity.
        dV = model.k_growth_per_ee_week * ee_weekly[i] * (model.v_max - volume[i - 1]) * dt_weeks
        volume[i] = volume[i - 1] + dV
        if volume[i] > model.v_max:
            volume[i] = model.v_max
        tanner[i] = model.tanner_stage(volume[i])

    weeks_axis = np.arange(n_weeks) * dt_weeks
    return BreastResult(weeks=weeks_axis, volume=volume, tanner=tanner, ee_curve=ee)
