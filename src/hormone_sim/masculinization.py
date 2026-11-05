"""Masculinization PD outcomes: cumulative-T-exposure models, trans-masculine.

The mirror of :mod:`hormone_sim.breast_development`, with the same honesty
rules. Outcomes are **population-typical, cumulative-exposure-driven,
saturating** — and genetics dominate every ceiling. This model is about
comfort and access (a body that stops feeling like a stranger), NOT about
"passing" as judged by others; the ceilings vary enormously and that is the
truth, not a limitation to hide.

**Calibration anchor:** menstrual cessation is the one masculinization
outcome with reasonably consistent reported timelines (most trans-masculine
people on stable T cease cycling within ~1-6 months; median ~5 months is
used here, flagged as calibrated with wide uncertainty). It anchors the
T-exposure → hazard mapping; the other outcomes are scaled from the same
drive function and carry correspondingly larger uncertainty.

Model shape (all saturating first-order in an *effective exposure* clock):

    drive(t)      = T_mean(t) / (T_mean(t) + Km)      (saturating in mean T)
    outcome(t)    = ceiling * (1 - exp(-eff_time(t) / tau))
    eff_time accrues at rate drive(t) — so higher troughs move every clock
    faster, which is the mechanistic version of why split-dose schedules and
    gels "feel" steadier than big-bolus biweekly plans.

Cessation is a hazard, not a clock: dP/dt = λ(t) · (1 − P), with
λ(t) = λ_max · drive(t) — cumulative probability of having stopped cycling.

Population mode samples genetics (ceilings, CV 0.15 — the dominant
uncertainty), drive sensitivity, and outcome time constants.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from .androgens import Androgen, simulate_androgen


@dataclass(frozen=True)
class MasculinizationModel:
    """PD outcome parameters (calibrated, LARGE uncertainty except where noted)."""

    # Cycle cessation hazard scale (per month at full drive). Anchored so
    # typical T (mean ~800-1000 ng/dL) gives median cessation ~5 months.
    cycle_lambda_max_per_mo: float = 0.19
    cycle_km_ngdl: float = 300.0

    # Saturating clocks (weeks at FULL drive). Genetics move the ceilings,
    # not the clocks.
    voice_tau_wk: float = 40.0
    hair_tau_wk: float = 150.0
    body_comp_tau_wk: float = 80.0
    bottom_tau_wk: float = 26.0

    # Population ceilings (fraction of full typical effect; genetics).
    voice_ceiling: float = 0.85      # rarely reaches cis-male typical fully
    hair_ceiling: float = 0.85
    body_comp_ceiling: float = 0.95
    bottom_ceiling: float = 1.00

    def drive(self, t_mean_ngdl: float) -> float:
        t = max(0.0, float(t_mean_ngdl))
        return t / (self.cycle_km_ngdl + t)


@dataclass
class MascResult:
    """Masculinization outcomes over time for one simulated body."""

    weeks: np.ndarray
    t_mean_ngdl: np.ndarray          # trailing mean T (the drive signal)
    t_ngdl: np.ndarray               # actual (peak-trough) T curve
    e2_pgml: np.ndarray
    cycle_ceased_prob: np.ndarray    # cumulative probability cycles have stopped
    voice: np.ndarray                # 0-1 index
    hair: np.ndarray
    body_comp: np.ndarray
    bottom: np.ndarray
    model: MasculinizationModel

    def report(self, weeks_out: Tuple[float, ...] = (13.0, 26.0, 52.0, 104.0)) -> str:
        m = self.model
        dt = self.weeks[1] - self.weeks[0] if self.weeks.size > 1 else 1.0
        def at(arr: np.ndarray, w: float) -> float:
            return float(arr[min(int(round(w / dt)), arr.size - 1)])
        lines = [
            f"Masculinization outcomes (population-typical; genetics dominate):",
            f"  T exposure: mean {self.t_mean_ngdl[-1]:.0f} ng/dL "
            f"(curve range {self.t_ngdl.min():.0f}-{self.t_ngdl.max():.0f})",
            f"  Cycle cessation: {self.cycle_ceased_prob[-1]*100:.0f}% by week "
            f"{self.weeks[-1]:.0f} "
            f"({self.weeks_to_frac(0.5) or 'n/r'} wk to 50% likelihood)",
            f"  {'outcome':<12} {'wk13':>6} {'wk26':>6} {'wk52':>6} {'wk104':>6}   "
            f"(index 0-1 of typical full effect)",
        ]
        for name, arr, ceil in (("voice", self.voice, m.voice_ceiling),
                                ("hair", self.hair, m.hair_ceiling),
                                ("body comp", self.body_comp, m.body_comp_ceiling),
                                ("bottom", self.bottom, m.bottom_ceiling)):
            cells = " ".join(f"{at(arr, w):>6.2f}" for w in weeks_out)
            lines.append(f"  {name:<12} {cells}   ceiling {ceil:.2f}")
        lines.append("")
        lines.append(
            "  Ceilings are genetics-dominant and vary enormously; this is a "
            "population\n  picture for comfort/expectation-setting, not an "
            "individual prediction.")
        return "\n".join(lines)

    def weeks_to_frac(self, frac: float) -> Optional[float]:
        """First week the cycle-cessation probability reaches ``frac``."""
        hits = np.nonzero(self.cycle_ceased_prob >= frac)[0]
        return float(self.weeks[hits[0]]) if hits.size else None


def simulate_masculinization(
    androgen: Androgen,
    dose_mg: float,
    interval_days: float,
    weeks: float = 104.0,
    dt_weeks: float = 0.25,
    model: Optional[MasculinizationModel] = None,
    body_weight_kg: float = 70.0,
    e2_ec50_ngdl: float = 170.0,
    baseline_t_ngdl: Optional[float] = None,
    baseline_e2_pgml: Optional[float] = None,
) -> MascResult:
    """Run the T PK and integrate the masculinization outcomes on it.

    Baseline overrides support the intersex packs (e.g. CAH's adrenal-driven
    endogenous T, Turner's near-zero gonadal output)."""
    model = model or MasculinizationModel()
    pk = simulate_androgen(androgen, dose_mg, interval_days, weeks,
                           body_weight_kg=body_weight_kg, e2_ec50_ngdl=e2_ec50_ngdl,
                           baseline_t_ngdl=baseline_t_ngdl,
                           baseline_e2_pgml=baseline_e2_pgml)

    hours = pk["hours"]
    # Weekly-grid resample of the PK curves (dt_weeks).
    step_h = dt_weeks * 7 * 24
    n = int(round(weeks * 7 * 24 / step_h)) + 1
    wk_idx = np.clip((np.arange(n) * step_h / hours[-1] * (hours.size - 1)).astype(int),
                     0, hours.size - 1)
    t_curve = pk["t_ngdl"][wk_idx]
    e2 = pk["e2_pgml"][wk_idx]

    # Trailing mean T over ~2 weeks as the drive signal.
    win = max(1, int(round(2.0 / dt_weeks)))
    t_mean = np.array([t_curve[max(0, i - win + 1):i + 1].mean() for i in range(n)])

    months_per_week = 12.0 / 52.0
    cycle = np.empty(n)
    voice = np.empty(n)
    hair = np.empty(n)
    body = np.empty(n)
    bottom = np.empty(n)
    surv = 1.0
    voice_e = hair_e = body_e = bottom_e = 0.0
    d_wk = dt_weeks

    for i in range(n):
        f = model.drive(t_mean[i])
        if i > 0:
            d_months = d_wk * months_per_week
            surv *= np.exp(-model.cycle_lambda_max_per_mo * f * d_months)
            voice_e += f * d_wk
            hair_e += f * d_wk
            body_e += f * d_wk
            bottom_e += f * d_wk
        cycle[i] = 1.0 - surv
        voice[i] = model.voice_ceiling * (1.0 - np.exp(-voice_e / model.voice_tau_wk))
        hair[i] = model.hair_ceiling * (1.0 - np.exp(-hair_e / model.hair_tau_wk))
        body[i] = model.body_comp_ceiling * (1.0 - np.exp(-body_e / model.body_comp_tau_wk))
        bottom[i] = model.bottom_ceiling * (1.0 - np.exp(-bottom_e / model.bottom_tau_wk))

    weeks_axis = np.arange(n) * dt_weeks
    return MascResult(
        weeks=weeks_axis, t_mean_ngdl=t_mean, t_ngdl=t_curve, e2_pgml=e2,
        cycle_ceased_prob=cycle, voice=voice, hair=hair, body_comp=body,
        bottom=bottom, model=model,
    )


# --- Population bands --------------------------------------------------------------
@dataclass
class MascPopulationResult:
    protocol_note: str
    weeks: np.ndarray
    t_ngdl: np.ndarray        # (n, t) actual curves (peak-trough visible)
    cycle_ceased_prob: np.ndarray
    voice: np.ndarray
    hair: np.ndarray
    body_comp: np.ndarray
    bottom: np.ndarray
    n: int

    def cycle_ceased_band(self, q: float = 0.8) -> Tuple[float, float, float]:
        """(p10, median, p90) of cessation probability at the horizon."""
        tail = self.cycle_ceased_prob[:, -1]
        return (float(np.percentile(tail, (1 - q) / 2 * 100)),
                float(np.median(tail)),
                float(np.percentile(tail, (1 + q) / 2 * 100)))

    def outcome_band(self, name: str, q: float = 0.8) -> Tuple[float, float, float]:
        arr = getattr(self, name)
        tail = arr[:, -1]
        return (float(np.percentile(tail, (1 - q) / 2 * 100)),
                float(np.median(tail)),
                float(np.percentile(tail, (1 + q) / 2 * 100)))

    def report(self) -> str:
        cp = self.cycle_ceased_band()
        lines = [
            f"Population bands (n={self.n}) — {self.protocol_note}",
            f"  Cycle cessation by horizon: p10 {cp[0]*100:.0f}% / median "
            f"{cp[1]*100:.0f}% / p90 {cp[2]*100:.0f}%",
        ]
        for name in ("voice", "hair", "body_comp", "bottom"):
            lo, med, hi = self.outcome_band(name)
            lines.append(f"  {name:<10} p10 {lo:.2f} / median {med:.2f} / "
                         f"p90 {hi:.2f}")
        lines.append(
            "  Genetics dominate the ceilings: the p10-p90 spread IS the "
            "point,\n  not noise to reduce.")
        return "\n".join(lines)


def simulate_masc_population(
    androgen: Androgen,
    dose_mg: float,
    interval_days: float,
    weeks: float = 104.0,
    n: int = 200,
    seed: int = 42,
    dt_weeks: float = 0.5,
) -> MascPopulationResult:
    """Monte-Carlo: PK + drive sensitivity + genetics (ceilings/taus)."""
    rng = np.random.default_rng(seed)

    t_stack, cyc_stack = [], []
    v_stack, h_stack, b_stack, bo_stack = [], [], [], []
    weeks_axis = None
    for _ in range(n):
        m = MasculinizationModel(
            cycle_lambda_max_per_mo=max(1e-3, 0.19 * rng.lognormal(0.0, 0.30)),
            cycle_km_ngdl=max(50.0, 300.0 * rng.lognormal(0.0, 0.30)),
            voice_tau_wk=max(5.0, 40.0 * rng.lognormal(0.0, 0.40)),
            hair_tau_wk=max(20.0, 150.0 * rng.lognormal(0.0, 0.40)),
            body_comp_tau_wk=max(10.0, 80.0 * rng.lognormal(0.0, 0.40)),
            bottom_tau_wk=max(5.0, 26.0 * rng.lognormal(0.0, 0.40)),
            voice_ceiling=float(np.clip(rng.normal(0.85, 0.15), 0.3, 1.0)),
            hair_ceiling=float(np.clip(rng.normal(0.85, 0.15), 0.3, 1.0)),
            body_comp_ceiling=float(np.clip(rng.normal(0.95, 0.10), 0.4, 1.0)),
            bottom_ceiling=1.0,
        )
        res = simulate_masculinization(
            androgen, dose_mg, interval_days, weeks=weeks, dt_weeks=dt_weeks,
            model=m,
        )
        if weeks_axis is None:
            weeks_axis = res.weeks
        t_stack.append(res.t_ngdl)
        cyc_stack.append(res.cycle_ceased_prob)
        v_stack.append(res.voice)
        h_stack.append(res.hair)
        b_stack.append(res.body_comp)
        bo_stack.append(res.bottom)

    return MascPopulationResult(
        protocol_note=f"{androgen.name} {dose_mg:g} mg every {interval_days:g} d",
        weeks=weeks_axis,
        t_ngdl=np.vstack(t_stack),
        cycle_ceased_prob=np.vstack(cyc_stack),
        voice=np.vstack(v_stack),
        hair=np.vstack(h_stack),
        body_comp=np.vstack(b_stack),
        bottom=np.vstack(bo_stack),
        n=n,
    )
