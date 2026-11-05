"""Sensitivity bands on the lactation model's milk parameters.

Every number in the lactation model is calibrated with large uncertainty and
extrapolated from postpartum cis literature — a single deterministic line
therefore overstates what anyone actually knows. This module replaces the one
line with a **plausible envelope**: Monte-Carlo sampling of the uncertain
parameters, reporting percentile bands of milk supply over time.

Sampled parameters (lognormal around population defaults, CV in parentheses):
    * ``milk_vmax_ml_day``    (0.25) — ceiling for a fully primed gland
    * ``milk_km_prolactin_ng_ml`` (0.30) — prolactin sensitivity of milk output
    * ``maturation_tau_weeks`` (0.60) — lactogenesis ramp; LEAST constrained
      parameter in the model, deliberately given the widest band
    * ``ductal_tau_weeks``    (0.35) — E2-driven ductal growth
    * ``la_tau_weeks``        (0.35) — P-driven lobuloalveolar growth
    * ``involution_rate_per_week`` (0.50) — post-drop decay when unstimulated

NOT sampled: prolactin dose-response (separate uncertainty, see
``galactagogues``), protocol structure. The bands answer "given the protocol,
what range of supply is plausible?" — they do not include the additional
uncertainty in the drug's prolactin effect itself.

Interpretation guidance: the bands are **only as good as the assumed CVs**,
which are themselves judgment calls informed by the spread in the literature.
The honest reading is the *shape* (steep ramp, plateau at a fraction of Vmax)
and the *order of magnitude*, not the precise percentiles.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .lactation import (
    LactationModel,
    LactationProtocol,
    LactationResult,
    simulate_lactation,
)


@dataclass(frozen=True)
class SensitivitySpec:
    """Which parameters vary, and by how much (lognormal CV).

    ``maturation_prolactin_threshold_ng_ml`` is included because it is the
    parameter that governs whether a LOW galactagogue dose works at all: the
    gland starts maturing only above a prolactin threshold, and where that
    threshold sits between individuals is unknown. Without sampling it, the
    low-dose regime shows a fake sharp cliff (everyone matures above exactly
    40 ng/mL, nobody below); with it, low doses show what is actually claimed:
    a heterogeneous zone where some modeled bodies respond and others don't.

    ``galactagogue_cvs`` sample the DRUG side of the chain — each modeled
    body gets its own prolactin Emax/ED50 and baseline — so the bands are
    END-TO-END: drug response *and* gland response, not gland-only.
    """

    cvs: Dict[str, float] = field(default_factory=lambda: {
        "milk_vmax_ml_day": 0.25,
        "milk_km_prolactin_ng_ml": 0.30,
        "maturation_tau_weeks": 0.60,
        "maturation_prolactin_threshold_ng_ml": 0.40,
        "ductal_tau_weeks": 0.35,
        "la_tau_weeks": 0.35,
        "involution_rate_per_week": 0.50,
    })

    galactagogue_cvs: Dict[str, float] = field(default_factory=lambda: {
        "prolactin_emax_ng_ml": 0.30,
        "prolactin_ed50_mg_day": 0.30,
        "baseline": 0.20,
    })

    def sample(self, rng: np.random.Generator) -> LactationModel:
        """Sample one individual's model parameters from the population."""
        base = LactationModel()
        kwargs = {}
        for name, cv in self.cvs.items():
            val = float(getattr(base, name)) * float(
                rng.lognormal(mean=0.0, sigma=cv)
            )
            kwargs[name] = max(val, 1e-6)
        return LactationModel(**kwargs)


@dataclass
class LactationPopulationResult:
    """Aggregate of many individual lactation simulations."""

    protocol: LactationProtocol
    weeks: np.ndarray
    milk: np.ndarray           # shape (n, t) mL/day

    def band(self, q: float = 0.8, hormone: str = "milk") -> Tuple:
        """(lower, median, upper) milk curves for the central ``q`` mass."""
        lo_pct = (1.0 - q) / 2.0 * 100.0
        hi_pct = (1.0 + q) / 2.0 * 100.0
        lo = np.percentile(self.milk, lo_pct, axis=0)
        hi = np.percentile(self.milk, hi_pct, axis=0)
        med = np.median(self.milk, axis=0)
        return lo, med, hi

    def steady_state(self, tail_weeks: float = 8.0) -> Tuple[float, float, float]:
        """(p10, median, p90) of milk over the last ``tail_weeks``."""
        dt = self.weeks[1] - self.weeks[0] if self.weeks.size > 1 else 1.0
        n_tail = max(1, int(round(tail_weeks / dt)))
        tail = self.milk[:, -n_tail:].mean(axis=1)
        return (float(np.percentile(tail, 10)),
                float(np.median(tail)),
                float(np.percentile(tail, 90)))

    def peak(self) -> Tuple[float, float, float]:
        peaks = self.milk.max(axis=1)
        return (float(np.percentile(peaks, 10)),
                float(np.median(peaks)),
                float(np.percentile(peaks, 90)))

    def weeks_to_supply(self, threshold_ml_day: float = 50.0,
                        lactation_start: Optional[float] = None) -> Tuple[Optional[float], Optional[float], Optional[float]]:
        """(p10, median, p90) weeks after the drop to reach ``threshold``."""
        start = lactation_start if lactation_start is not None else self.protocol.prep_weeks
        dt = self.weeks[1] - self.weeks[0]
        out = []
        for i in range(self.milk.shape[0]):
            reached = None
            for j, w in enumerate(self.weeks):
                if w >= start and self.milk[i, j] >= threshold_ml_day:
                    reached = float(w - start)
                    break
            out.append(reached)
        vals = [v for v in out if v is not None]
        frac = len(vals) / len(out)
        if not vals:
            return (None, None, None)
        p10, med, p90 = (float(np.percentile(vals, 10)),
                         float(np.median(vals)),
                         float(np.percentile(vals, 90)))
        # If fewer than half reached it, the median is 'not reached'.
        if frac < 0.5:
            med = None
        if frac < 0.1:
            p10 = None
        if frac < 0.9:
            p90 = None
        return (p10, med, p90)

    def report(self) -> str:
        ss10, ss50, ss90 = self.steady_state()
        pk10, pk50, pk90 = self.peak()
        w10, w50, w90 = self.weeks_to_supply()
        def fmt(v, unit=""):
            return f"{v:.0f}{unit}" if v is not None else "n/r"
        lines = [
            f"Lactation sensitivity bands (n={self.milk.shape[0]} samples):",
            f"  Steady-state milk (last 8 wk): "
            f"p10 {fmt(ss10)} / median {fmt(ss50)} / p90 {fmt(ss90)} mL/day",
            f"  Peak milk: "
            f"p10 {fmt(pk10)} / median {fmt(pk50)} / p90 {fmt(pk90)} mL/day",
            f"  Weeks to >50 mL/day after the drop: "
            f"p10 {w10 if w10 is None else f'{w10:.1f}'} / "
            f"median {w50 if w50 is None else f'{w50:.1f}'} / "
            f"p90 {w90 if w90 is None else f'{w90:.1f}'}  (n/r = not reached)",
            "",
            "  Caveats: bands reflect only the sampled parameter CVs (judgment",
            "  calls, not measured population spread); the drug's prolactin",
            "  dose-response uncertainty is NOT included; no controlled trials",
            "  in trans women. Read the shape and order of magnitude.",
        ]
        return "\n".join(lines)


def simulate_lactation_population(
    protocol: LactationProtocol,
    n: int = 200,
    seed: int = 42,
    spec: Optional[SensitivitySpec] = None,
    total_weeks: Optional[float] = None,
    dt_weeks: float = 0.25,
    include_drug_uncertainty: bool = True,
) -> LactationPopulationResult:
    """Monte-Carlo the lactation model's uncertain parameters.

    By default the sampled parameters include the DRUG side (each modeled
    body's prolactin Emax/ED50 and baseline, per ``spec.galactagogue_cvs``),
    making the envelope end-to-end: drug response *and* gland response.
    Pass ``include_drug_uncertainty=False`` to hold the agent's dose-response
    at its population value (gland-only bands, the pre-0.5 behavior).
    """
    from .galactagogues import get_galactagogue, PROLACTIN_BASELINE_NG_ML

    spec = spec or SensitivitySpec()
    rng = np.random.default_rng(seed)
    g = get_galactagogue(protocol.galactagogue)

    stack: List[np.ndarray] = []
    weeks = None
    for _ in range(n):
        m = spec.sample(rng)

        p_drug = None
        if include_drug_uncertainty:
            emax = g.prolactin_emax_ng_ml * float(
                rng.lognormal(0.0, spec.galactagogue_cvs["prolactin_emax_ng_ml"]))
            ed50 = max(1e-3, g.prolactin_ed50_mg_day * float(
                rng.lognormal(0.0, spec.galactagogue_cvs["prolactin_ed50_mg_day"])))
            base = PROLACTIN_BASELINE_NG_ML * float(
                rng.lognormal(0.0, spec.galactagogue_cvs["baseline"]))
            d = max(0.0, float(protocol.galactagogue_dose_mg_day))
            p_drug = base + emax * d / (ed50 + d)

        res = simulate_lactation(protocol, total_weeks=total_weeks,
                                 dt_weeks=dt_weeks, model=m,
                                 individual_prolactin_ng_ml=p_drug)
        if weeks is None:
            weeks = res.weeks
        stack.append(res.milk_ml_day)

    return LactationPopulationResult(
        protocol=protocol,
        weeks=weeks,
        milk=np.vstack(stack),
    )


# --- Low-dose sweep: "useful for which bodies?" ------------------------------
@dataclass
class DoseSweepRow:
    """One dose's outcome, with between-body heterogeneity made explicit."""

    agent: str
    dose_mg_day: float
    prolactin_ng_ml: float
    milk_p10: float
    milk_p50: float
    milk_p90: float
    frac_reaching_threshold: float   # fraction of sampled bodies >= threshold

    def report(self) -> str:
        return (
            f"{self.dose_mg_day:>6.1f} mg/d | PRL {self.prolactin_ng_ml:>5.0f} | "
            f"milk p10 {self.milk_p10:>5.0f} / p50 {self.milk_p50:>5.0f} / "
            f"p90 {self.milk_p90:>5.0f} mL/d | "
            f"{self.frac_reaching_threshold*100:>5.1f}% of modeled bodies "
            f"reach the threshold"
        )


def dose_sweep_population(
    agent: str,
    doses: List[float],
    prep_weeks: float = 30.0,
    progesterone_start_weeks: Optional[float] = None,
    pumping_sessions_per_day: float = 8.0,
    n: int = 150,
    seed: int = 42,
    threshold_ml_day: float = 100.0,
    spec: Optional[SensitivitySpec] = None,
) -> List[DoseSweepRow]:
    """Sweep galactagogue doses and quantify between-body heterogeneity.

    For each dose: run the population simulation and report the milk-supply
    envelope plus **what fraction of modeled bodies reach a meaningful
    supply**. This is the quantified version of the heterogeneity question:
    at low doses the answer is not yes/no but *which fraction*, driven by the
    (unknown, sampled) spread in prolactin sensitivity and maturation
    threshold. Uses the SAME seed for every dose so rows are paired — the
    same "bodies" are evaluated at each dose.
    """
    from .galactagogues import get_galactagogue, prolactin_response

    g = get_galactagogue(agent)
    spec = spec or SensitivitySpec()
    rows: List[DoseSweepRow] = []

    for dose in doses:
        protocol = LactationProtocol(
            galactagogue=agent,
            galactagogue_dose_mg_day=float(dose),
            prep_weeks=prep_weeks,
            progesterone_start_weeks=progesterone_start_weeks,
            pumping_sessions_per_day=pumping_sessions_per_day,
        )
        # Same seed per dose -> paired samples across doses.
        pop = simulate_lactation_population(
            protocol, n=n, seed=seed, spec=spec,
        )
        p10, p50, p90 = pop.steady_state()
        tail_dt = pop.weeks[1] - pop.weeks[0]
        n_tail = max(1, int(round(8.0 / tail_dt)))
        tail = pop.milk[:, -n_tail:].mean(axis=1)
        frac = float(np.mean(tail >= threshold_ml_day))
        rows.append(DoseSweepRow(
            agent=g.name,
            dose_mg_day=float(dose),
            prolactin_ng_ml=prolactin_response(g, dose),
            milk_p10=p10, milk_p50=p50, milk_p90=p90,
            frac_reaching_threshold=frac,
        ))
    return rows


def dose_sweep_table(rows: List[DoseSweepRow],
                     threshold_ml_day: float = 100.0) -> str:
    """Human-readable low-dose sweep table."""
    agent = rows[0].agent if rows else "?"
    lines = [
        f"Low-dose sweep — {agent} (paired samples; n per dose = "
        f"same bodies at each dose):",
        f"{'dose':>10} | {'PRL ng/mL':>9} | {'milk p10/p50/p90 mL/d':>24} | "
        f"{'% bodies >=' + f'{threshold_ml_day:.0f}':>16}",
        "-" * 70,
    ]
    for r in rows:
        lines.append(r.report())
    lines.append("")
    lines.append(
        "  Reading this: the p50 column is the 'average body'; the p10-p90"
    )
    lines.append(
        "  spread and the % column are the heterogeneity question. Low doses"
    )
    lines.append(
        "  live in the marginal zone — useful for some modeled physiologies,"
    )
    lines.append(
        "  not others, which is why titration-to-effect (not fixed dosing)"
    )
    lines.append(
        "  is the rational protocol."
    )
    return "\n".join(lines)
