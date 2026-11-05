"""Erythrocytosis risk framing for testosterone therapy (trans-masculine).

The deliberate mirror image of the metoclopramide/TD framing in
:mod:`hormone_sim.risk_exposure` — and the case where the asymmetry argument
lands the **other** way:

* Testosterone raises hematocrit (via hepcidin suppression → increased iron
  availability → erythropoiesis). This is the *dominant* lab-tracked risk of
  trans-masculine GAHT.
* Unlike TD: it is **acute-to-subacute, fully lab-monitorable (a CBC), and
  reversible** on dose reduction or pause. There IS a biomarker. Risk is
  dose- and route-dependent (injectable esters > gel, consistent with the
  peak-exposure story). This is the "certain bodies, certain conditions"
  heterogeneity case where monitoring genuinely works — the same structure
  as the user's estradiol/VTE experience, on the masculinization side.
* Common clinical practice (flagged as practice, not law): check CBC at
  baseline, ~3 and 6 months, then annually; evaluate a Hct climbing above
  ~52%; consider holding/reducing T above ~54% (also thrombosis-avoidance
  measures). Exact thresholds are clinician judgment; this module outputs
  *projected trajectory + tiering*, not a decision.

Model: saturating first-order rise toward a T-dependent target.

    Hct(t) = baseline + rise_max * drive(T_mean) * (1 - exp(-t / tau))

with drive saturating in mean T, rise_max ~13 points (female-typical ~40 →
male-typical ceiling ~48-53), tau ~10 weeks (RBC-mass response scale).
Calibrated with moderate uncertainty; the *shape* (rise over ~3-6 months,
plateau, reversible) is well described in the literature.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from .androgens import Androgen, simulate_androgen


@dataclass(frozen=True)
class HctModel:
    """Hematocrit response parameters (moderate uncertainty)."""

    baseline_hct: float = 40.0        # female-typical
    max_rise_points: float = 13.0     # saturating ceiling of T-driven rise
    km_ngdl: float = 400.0            # mean-T level for half-maximal drive
    tau_wk: float = 10.0              # RBC-mass response time constant

    def drive(self, t_mean_ngdl: float) -> float:
        t = max(0.0, float(t_mean_ngdl))
        return t / (self.km_ngdl + t)

    def target_hct(self, t_mean_ngdl: float) -> float:
        return self.baseline_hct + self.max_rise_points * self.drive(t_mean_ngdl)


@dataclass
class HctResult:
    weeks: np.ndarray
    hct: np.ndarray
    t_mean_ngdl: np.ndarray
    model: HctModel

    def steady_state_hct(self, tail_weeks: float = 8.0) -> float:
        n_tail = max(1, int(round(tail_weeks / (self.weeks[1] - self.weeks[0]))))
        return float(self.hct[-n_tail:].mean())

    def weeks_to_hct(self, level: float) -> Optional[float]:
        hits = np.nonzero(self.hct >= level)[0]
        return float(self.weeks[hits[0]]) if hits.size else None


def simulate_hct(
    androgen: Androgen,
    dose_mg: float,
    interval_days: float,
    weeks: float = 52.0,
    dt_weeks: float = 0.25,
    model: Optional[HctModel] = None,
    body_weight_kg: float = 70.0,
) -> HctResult:
    """Project the hematocrit trajectory under a T regimen."""
    model = model or HctModel()
    pk = simulate_androgen(androgen, dose_mg, interval_days, weeks,
                           body_weight_kg=body_weight_kg)
    step_h = dt_weeks * 7 * 24
    n = int(round(weeks * 7 * 24 / step_h)) + 1
    idx = np.clip((np.arange(n) * step_h / pk["hours"][-1]
                   * (pk["hours"].size - 1)).astype(int), 0, pk["hours"].size - 1)
    t_curve = pk["t_ngdl"][idx]
    win = max(1, int(round(2.0 / dt_weeks)))
    t_mean = np.array([t_curve[max(0, i - win + 1):i + 1].mean() for i in range(n)])

    hct = np.empty(n)
    hct[0] = model.baseline_hct
    k = 1.0 / model.tau_wk
    for i in range(1, n):
        target = model.target_hct(t_mean[i])
        hct[i] = hct[i - 1] + k * (target - hct[i - 1]) * dt_weeks

    weeks_axis = np.arange(n) * dt_weeks
    return HctResult(weeks=weeks_axis, hct=hct, t_mean_ngdl=t_mean, model=model)


def hct_risk_report(
    androgen: Androgen,
    dose_mg: float,
    interval_days: float,
    body_weight_kg: float = 70.0,
) -> str:
    """Structured risk framing: projected trajectory + tiers + lab schedule.

    The mirror-image note to the TD tool: this risk IS biomarker-monitored,
    dose/route-dependent, and reversible — the user's 'certain bodies'
    monitoring argument applies with full force here.
    """
    res = simulate_hct(androgen, dose_mg, interval_days,
                       body_weight_kg=body_weight_kg)
    ss = res.steady_state_hct()
    route_note = ("injectable ester: larger peak exposure drives a higher "
                  "plateau than gel at the same weekly dose"
                  if androgen.input_kind == "depot"
                  else "gel: near-flat exposure, lowest Hct plateau of the routes")

    if ss < 48:
        tier = "low (projected plateau within common reference range)"
    elif ss < 52:
        tier = "elevated (plateau above female-typical; monitor per schedule)"
    else:
        tier = "high (projected plateau near/above the ~52% evaluation zone)"

    lines = [
        f"Erythrocytosis risk framing — {androgen.name} {dose_mg:g} mg every "
        f"{interval_days:g} d",
        f"  Projected Hct plateau: ~{ss:.0f}% "
        f"(baseline {res.model.baseline_hct:.0f}%)",
        f"  Time to plateau: ~{res.model.tau_wk * 3:.0f} weeks (3x time constant)",
        f"  Tier: {tier}",
        f"  Route note: {route_note}",
        "",
        "  Why this is the 'good' risk case (vs the TD asymmetry):",
        "    - biomarker exists: a CBC is the monitor; no exam required",
        "    - dose- and route-dependent: schedule changes move it",
        "    - reversible: dose reduction / pause brings Hct back down",
        "",
        "  Commonly used monitoring cadence (practice, not a prescription):",
        "    - baseline CBC",
        "    - ~3 and ~6 months (the rise window)",
        "    - 12 months, then annually",
        "    - evaluate a Hct climbing above ~52%; holding/reducing T above "
        "~54% is common practice",
        "",
        "  This is a trajectory-projection tool, not a dosing decision; "
        "the CBC is the decision input.",
    ]
    return "\n".join(lines)
