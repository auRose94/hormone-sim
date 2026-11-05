"""High-level simulation orchestration.

Ties the PK (estradiol) and PD (testosterone) models together for a full
regimen over a specified horizon, returning a :class:`SimResult` with the time
series and summary statistics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .params import (
    EstradiolPK,
    AndrogenPD,
    TARGET_E2_PGML,
    TARGET_T_NGDL,
    e2_pgml_to_pmol,
    t_ngdl_to_nmol,
)
from .pk import build_schedule, simulate_pk
from .pd import suppress_testosterone
from .regimens import Regimen


@dataclass
class SimResult:
    """A single simulation's output."""

    regimen: Regimen
    weeks: float
    hours: np.ndarray          # time axis (hours)
    e2_pgml: np.ndarray        # estradiol (pg/mL)
    t_ngdl: np.ndarray         # testosterone (ng/dL)
    seed: Optional[int] = None

    @property
    def e2_pmol(self) -> np.ndarray:
        return e2_pgml_to_pmol(self.e2_pgml)

    @property
    def t_nmol(self) -> np.ndarray:
        return t_ngdl_to_nmol(self.t_ngdl)

    # --- summary statistics -------------------------------------------------
    def steady_state_e2_pgml(self, tail_frac: float = 0.25) -> float:
        """Mean estradiol over the last ``tail_frac`` of the horizon."""
        n = self.e2_pgml.size
        cut = int(n * (1.0 - tail_frac))
        return float(self.e2_pgml[cut:].mean())

    def steady_state_t_ngdl(self, tail_frac: float = 0.25) -> float:
        n = self.t_ngdl.size
        cut = int(n * (1.0 - tail_frac))
        return float(self.t_ngdl[cut:].mean())

    def e2_in_target(self) -> bool:
        lo, hi = TARGET_E2_PGML
        return lo <= self.steady_state_e2_pgml() <= hi

    def t_suppressed(self) -> bool:
        return self.steady_state_t_ngdl() <= TARGET_T_NGDL[1]

    def report(self) -> str:
        e2 = self.steady_state_e2_pgml()
        t = self.steady_state_t_ngdl()
        lines = [
            self.regimen.summary(),
            f"  Horizon: {self.weeks:g} weeks ({self.hours[-1]:.0f} h)",
            f"  Steady-state E2: {e2:.1f} pg/mL ({e2_pgml_to_pmol(e2):.0f} pmol/L)"
            f"  [{'in target' if self.e2_in_target() else 'OUT of target'}]",
            f"  Steady-state T:  {t:.1f} ng/dL ({t_ngdl_to_nmol(t):.1f} nmol/L)"
            f"  [{'suppressed' if self.t_suppressed() else 'NOT suppressed'}]",
        ]
        return "\n".join(lines)


def simulate(
    regimen: Regimen,
    weeks: float = 12.0,
    pk: Optional[EstradiolPK] = None,
    pd: Optional[AndrogenPD] = None,
    dt_hours: float = 0.5,
) -> SimResult:
    """Run a single population-typical simulation.

    ``weeks`` is the simulation horizon. Returns a :class:`SimResult`.
    """
    pk = pk or EstradiolPK()
    pd = pd or AndrogenPD()

    hours = weeks * 7.0 * 24.0
    route = regimen.estradiol.route
    spec = regimen.estradiol

    schedule = build_schedule(
        daily_rate_ug_day=spec.daily_rate_ug_day,
        route=spec.route_params(),
        period_hours=hours,
        injection_ug=spec.injection_ug,
        injection_interval_hours=spec.injection_interval_hours,
    )

    e2 = simulate_pk(schedule, pk, hours, dt_hours=dt_hours, route=route)
    t = suppress_testosterone(e2, pd, window_hours=72.0, dt_hours=dt_hours,
                              agents=regimen.anti_androgen)

    time_axis = np.arange(e2.size) * dt_hours
    return SimResult(
        regimen=regimen,
        weeks=weeks,
        hours=time_axis,
        e2_pgml=e2,
        t_ngdl=t,
    )
