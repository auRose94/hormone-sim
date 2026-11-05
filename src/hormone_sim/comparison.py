"""Head-to-head comparison of galactagogues at equieffective prolactin doses.

The practical question isn't "which drug raises prolactin more" — any D2
antagonist can be dose-adjusted to hit the same prolactin level. The question
is: **at equieffective prolactin, what differs?** The answer separates into:

1. **Pharmacokinetics** — how fast the drug reaches steady state (metoclopramide
   t½ ~5.5 h: 1-2 days; domperidone t½ ~15 h: 3-5 days), and how lumpy the
   prolactin curve is between doses (a function of t½ vs dosing interval).
2. **Pharmacodynamics** — dose required to get there (potency at the receptor
   and in the gut), and the shape of the Emax curve (how much headroom is left
   before the ceiling).
3. **Safety profile** — metoclopramide crosses the blood-brain barrier
   (CNS effects: EPS, depression, boxed-warning tardive dyskinesia);
   domperidone is largely peripherally restricted (fewer CNS effects) but
   carries QT-prolongation risk and is not FDA-approved.

This module finds, for each agent, the daily dose that reaches a target
prolactin level (inverting the Emax model), then simulates the full
Newman-Goldfarb protocol with each and reports supply, onset, and the
PK/safety tradeoffs side by side.

Same caveat as everywhere else: no controlled trials in trans women; the
comparison is mechanistic extrapolation, not clinical guidance.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from .galactagogues import (
    GALACTAGOGUES,
    Galactagogue,
    get_galactagogue,
    prolactin_response,
    PROLACTIN_BASELINE_NG_ML,
)
from .lactation import LactationProtocol, LactationResult, simulate_lactation


def equieffective_dose(g: Galactagogue, target_prolactin_ng_ml: float) -> float:
    """Daily oral dose (mg/day) reaching ``target_prolactin_ng_ml``.

    Inverts the Emax model: rise = P_target - baseline;
    dose = ED50 * rise / (Emax - rise). Raises ValueError if the target
    exceeds the agent's ceiling (baseline + Emax).
    """
    rise = target_prolactin_ng_ml - PROLACTIN_BASELINE_NG_ML
    if rise <= 0:
        return 0.0
    if rise >= g.prolactin_emax_ng_ml:
        raise ValueError(
            f"{g.name}: target {target_prolactin_ng_ml:g} ng/mL exceeds its "
            f"ceiling ({PROLACTIN_BASELINE_NG_ML + g.prolactin_emax_ng_ml:g} ng/mL); "
            "no dose can reach it."
        )
    return g.prolactin_ed50_mg_day * rise / (g.prolactin_emax_ng_ml - rise)


@dataclass
class AgentComparison:
    """One agent's row in a head-to-head comparison."""

    name: str
    dose_mg_day: float
    prolactin_ng_ml: float
    result: LactationResult
    half_life_h: float
    approved_us: bool
    safety_note: str

    # Convenience passthroughs
    @property
    def ss_milk(self) -> float:
        return self.result.steady_state_milk()

    @property
    def peak_milk(self) -> float:
        return self.result.peak_milk()

    @property
    def weeks_to_supply(self) -> Optional[float]:
        return self.result.weeks_to_first_milk()


def compare_galactagogues(
    agents: List[str],
    target_prolactin_ng_ml: float = 130.0,
    prep_weeks: float = 30.0,
    pumping_sessions_per_day: float = 8.0,
    progesterone_start_weeks: Optional[float] = None,
    total_weeks: Optional[float] = None,
    model=None,
) -> List[AgentComparison]:
    """Compare agents at equieffective prolactin over the same protocol.

    Every agent is dose-adjusted to reach ``target_prolactin_ng_ml``, then run
    through an identical Newman-Goldfarb protocol. Returns rows sorted by
    steady-state milk (desc).
    """
    protocol_kwargs = dict(
        prep_weeks=prep_weeks,
        progesterone_start_weeks=progesterone_start_weeks,
        pumping_sessions_per_day=pumping_sessions_per_day,
    )
    rows: List[AgentComparison] = []
    for name in agents:
        g = get_galactagogue(name)
        dose = equieffective_dose(g, target_prolactin_ng_ml)
        protocol = LactationProtocol(
            galactagogue=name,
            galactagogue_dose_mg_day=dose,
            **protocol_kwargs,
        )
        res = simulate_lactation(protocol, total_weeks=total_weeks, model=model)
        rows.append(AgentComparison(
            name=g.name,
            dose_mg_day=dose,
            prolactin_ng_ml=prolactin_response(g, dose),
            result=res,
            half_life_h=g.half_life_h,
            approved_us=g.approved_us,
            safety_note=g.safety_note,
        ))
    rows.sort(key=lambda r: r.ss_milk, reverse=True)
    return rows


def comparison_table(rows: List[AgentComparison]) -> str:
    """Human-readable side-by-side table."""
    lines = [
        "Equieffective comparison (same target prolactin, same protocol):",
        f"{'agent':>15} | {'dose mg/d':>9} | {'PRL ng/mL':>9} | {'SS milk':>8} | "
        f"{'peak':>6} | {'>50 mL/wk':>9} | {'t1/2 h':>6} | {'US':>4}",
        "-" * 84,
    ]
    for r in rows:
        wk = r.weeks_to_supply
        wk_txt = f"{wk:.1f}" if wk is not None else "never"
        lines.append(
            f"{r.name:>15} | {r.dose_mg_day:>9.1f} | {r.prolactin_ng_ml:>9.0f} | "
            f"{r.ss_milk:>8.0f} | {r.peak_milk:>6.0f} | {wk_txt:>9} | "
            f"{r.half_life_h:>6.1f} | {'yes' if r.approved_us else 'no':>4}"
        )
    return "\n".join(lines)


def tradeoff_notes(rows: List[AgentComparison]) -> str:
    """The pharmacology that the table can't show: onset, lumpy-ness, safety."""
    lines = ["PK / safety tradeoffs at equieffective prolactin:"]
    for r in rows:
        g = get_galactagogue(r.name)
        # Time to steady state ~ 4-5 half-lives.
        ss_days = 4.5 * g.half_life_h / 24.0
        lines.append(
            f"  {r.name}: steady state in ~{ss_days:.1f} d; dose {r.dose_mg_day:.0f} mg/day. "
            f"{r.safety_note}"
        )
    lines.append(
        "  NOTE: at equal *average* prolactin, a shorter t½ drug gives a lumpier"
    )
    lines.append(
        "  prolactin curve between doses (more peak-trough swing at the same"
    )
    lines.append(
        "  daily dose split); the lactation-phase model here uses steady-state"
    )
    lines.append(
        "  averages, so within-day fluctuation is a qualitative caveat, not"
    )
    lines.append("  quantified output.")
    return "\n".join(lines)
