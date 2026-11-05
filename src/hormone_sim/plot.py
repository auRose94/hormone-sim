"""Plotting for simulation results.

Produces the standard figure: estradiol (top) and testosterone (bottom) vs time
in weeks, with guideline target bands shaded. matplotlib is imported lazily so
the package works headless without it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from .params import TARGET_E2_PGML, TARGET_T_NGDL
from .sim import SimResult
from .population import PopulationResult


def _weeks_axis(hours: np.ndarray) -> np.ndarray:
    return hours / (7.0 * 24.0)


def plot_single(result: SimResult, out_path: Optional[str] = None) -> "object":
    """Plot a single :class:`SimResult`. Returns the matplotlib Figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    weeks = _weeks_axis(result.hours)
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True, constrained_layout=True
    )

    # Estradiol.
    ax1.plot(weeks, result.e2_pgml, color="#c23b6a", lw=1.8, label="Estradiol")
    lo, hi = TARGET_E2_PGML
    ax1.axhspan(lo, hi, color="#c23b6a", alpha=0.12, label=f"E2 target {lo:.0f}-{hi:.0f} pg/mL")
    ax1.axhline(hi, color="#c23b6a", ls=":", lw=0.8)
    ax1.axhline(lo, color="#c23b6a", ls=":", lw=0.8)
    ax1.set_ylabel("Estradiol (pg/mL)")
    ax1.set_title(f"Estradiol & Testosterone over time — {result.regimen.name}")

    # Testosterone.
    ax2.plot(weeks, result.t_ngdl, color="#2a6fb0", lw=1.8, label="Testosterone")
    ax2.axhspan(0, TARGET_T_NGDL[1], color="#2a6fb0", alpha=0.12,
                label=f"T target <{TARGET_T_NGDL[1]:.0f} ng/dL")
    ax2.axhline(TARGET_T_NGDL[1], color="#2a6fb0", ls=":", lw=0.8)
    ax2.set_ylabel("Testosterone (ng/dL)")
    ax2.set_xlabel("Time (weeks)")

    for ax in (ax1, ax2):
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle(result.regimen.summary().splitlines()[0], fontsize=11)
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=130)
    return fig


def plot_population(pop: PopulationResult, out_path: Optional[str] = None) -> "object":
    """Plot population bands (median + 80% interval) for E2 and T."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    weeks = _weeks_axis(pop.hours)
    e2_lo, e2_med, e2_hi = pop.percentile_band(0.8, hormone="e2")
    t_lo, t_med, t_hi = pop.percentile_band(0.8, hormone="t")

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True, constrained_layout=True
    )

    ax1.plot(weeks, e2_med, color="#c23b6a", lw=2, label="median")
    ax1.fill_between(weeks, e2_lo, e2_hi, color="#c23b6a", alpha=0.2,
                     label="80% interval")
    lo, hi = TARGET_E2_PGML
    ax1.axhspan(lo, hi, color="#c23b6a", alpha=0.12)
    ax1.set_ylabel("Estradiol (pg/mL)")
    ax1.set_title(
        f"Population bands (n={len(pop.individuals)}, "
        f"E2 in target {pop.fraction_e2_in_target()*100:.0f}%, "
        f"T suppressed {pop.fraction_t_suppressed()*100:.0f}%) — "
        f"{pop.individuals[0].regimen.name}"
    )

    ax2.plot(weeks, t_med, color="#2a6fb0", lw=2, label="median")
    ax2.fill_between(weeks, t_lo, t_hi, color="#2a6fb0", alpha=0.2, label="80% interval")
    ax2.axhspan(0, TARGET_T_NGDL[1], color="#2a6fb0", alpha=0.12)
    ax2.set_ylabel("Testosterone (ng/dL)")
    ax2.set_xlabel("Time (weeks)")

    for ax in (ax1, ax2):
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.3)

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=130)
    return fig
