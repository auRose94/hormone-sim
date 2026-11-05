"""Population (Monte-Carlo) simulation.

Samples between-subject variability in PK/PD parameters and runs many
individuals to produce confidence bands (e.g. the 10th–90th percentile envelope
of steady-state estradiol), not a single line. This is the part of the model
that answers "what fraction of a population reaches target?" — a study-design
question, not a per-patient prediction.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from .params import (
    EstradiolPK,
    AndrogenPD,
    Variability,
    e2_pgml_to_pmol,
    t_ngdl_to_nmol,
)
from .sim import simulate, SimResult


@dataclass
class PopulationResult:
    """Aggregate of many individual simulations."""

    individuals: List[SimResult]
    weeks: float
    hours: np.ndarray
    e2_pgml: np.ndarray      # shape (n, t)
    t_ngdl: np.ndarray       # shape (n, t)

    def percentile_band(self, q: float = 0.8, hormone: str = "e2") -> tuple:
        """Return (lower, median, upper) curves for the given central mass.

        ``q`` is the central fraction (e.g. 0.8 -> the middle 80%).
        ``hormone`` selects "e2" (estradiol) or "t" (testosterone).
        """
        series = self.e2_pgml if hormone == "e2" else self.t_ngdl
        lo_pct = (1.0 - q) / 2.0 * 100.0
        hi_pct = (1.0 + q) / 2.0 * 100.0
        lo = np.percentile(series, lo_pct, axis=0)
        hi = np.percentile(series, hi_pct, axis=0)
        med = np.median(series, axis=0)
        return lo, med, hi

    def steady_state_e2_mean(self) -> float:
        return float(self.e2_pgml[:, -1].mean())

    def fraction_e2_in_target(self) -> float:
        """Fraction of individuals whose steady-state E2 is in the target band."""
        from .params import TARGET_E2_PGML
        lo, hi = TARGET_E2_PGML
        last = self.e2_pgml[:, -1]
        return float(np.mean((last >= lo) & (last <= hi)))

    def fraction_t_suppressed(self, upper: float = 50.0) -> float:
        last = self.t_ngdl[:, -1]
        return float(np.mean(last <= upper))


def run_population(
    build_regimen,
    weeks: float = 12.0,
    n: int = 200,
    seed: int = 42,
    dt_hours: float = 0.5,
) -> PopulationResult:
    """Run ``n`` individuals with sampled parameters.

    ``build_regimen`` is a zero-arg callable returning a :class:`Regimen` (so
    the caller can vary the regimen per run). Variability is sampled from
    :class:`Variability` (log-normal for positive quantities).
    """
    rng = np.random.default_rng(seed)
    vb = Variability()

    individuals: List[SimResult] = []
    e2_stack: List[np.ndarray] = []
    t_stack: List[np.ndarray] = []

    for _ in range(n):
        # Sample individual PK/PD parameters around the population mean on the
        # model's tunable knobs (half-life, steady-state factor).
        half_life = rng.normal(1.0, vb.half_life_cv) * EstradiolPK.half_life_h_val
        half_life = max(half_life, 1.0)
        ssf = rng.normal(1.0, vb.steady_state_cv) * EstradiolPK.steady_state_factor
        ssf = max(ssf, 0.1)
        ec50 = rng.normal(1.0, vb.e2_ec50_cv) * AndrogenPD.e2_ec50_pgml
        ec50 = max(ec50, 10.0)
        baseline_t = rng.normal(1.0, vb.baseline_t_cv) * 300.0
        baseline_t = max(baseline_t, 50.0)

        pk = EstradiolPK(half_life_h_val=half_life, steady_state_factor=ssf)
        pd = AndrogenPD(e2_ec50_pgml=ec50)

        sim = simulate(build_regimen(), weeks=weeks, pk=pk, pd=pd, dt_hours=dt_hours)
        # Override baseline testosterone for this individual (PD starts from baseline).
        sim.t_ngdl = baseline_t / 300.0 * sim.t_ngdl
        individuals.append(sim)
        e2_stack.append(sim.e2_pgml)
        t_stack.append(sim.t_ngdl)

    hours = individuals[0].hours
    return PopulationResult(
        individuals=individuals,
        weeks=weeks,
        hours=hours,
        e2_pgml=np.vstack(e2_stack),
        t_ngdl=np.vstack(t_stack),
    )
