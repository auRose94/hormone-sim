"""Estradiol one-compartment pharmacokinetics.

The model is intentionally transparent: a single well whose mass changes by an
explicit **input function** minus **first-order elimination**. All route
differences live in the input function and the depot absorption rate —
clearance and volume are molecule properties shared across routes.

Two-state model:

    A  = estradiol mass in the body (µg)        C = A / Vd  (pg/mL == µg/L)
    D  = mass in a "depot" site (µg): the gut lumen (oral) or the oil/water
         depot (injectable ester). The bolus lands here, then is absorbed
         first-order into the body.

Dynamics (per hour):
    dA/dt = production  +  k_abs * D  -  Cl * C  +  body_input(t)
    dD/dt = bolus(t)    -  k_abs * D

``production`` is a constant endogenous estradiol production rate calibrated so
that, with **zero exogenous input**, the body sits at ``BASELINE_E2_PGML`` (the
untreated AMAB baseline). It equals ``Cl * baseline_C`` so the baseline is a true
steady state: at C = baseline, production exactly balances clearance. This makes
"off therapy" = baseline rather than collapsing to zero.

Units note: 1 pg/mL estradiol = 1 µg/L, so mass(µg) = C(pg/mL) * Vd(L).

This is a *population-typical* model; individual variability is layered on in
:mod:`hormone_sim.population`. Nothing here is patient-specific dosing.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .params import EstradiolPK, RouteParams, BASELINE_E2_PGML


@dataclass
class DoseEvent:
    """A single discrete dosing event: ``amount_ug`` delivered into the depot at
    ``t_hours``."""

    t_hours: float
    amount_ug: float


@dataclass
class DosingSchedule:
    """A regimen's timeline, expanded into a continuous flux and discrete events.

    * ``continuous_rate_ug_h`` — zero-order input (µg/h) while on (transdermal).
    * ``events`` — discrete bolus doses (µg) placed into the depot.
    * ``k_abs_h`` — first-order depot->body absorption rate (1/h).
    """

    continuous_rate_ug_h: float = 0.0
    events: List[DoseEvent] = field(default_factory=list)
    k_abs_h: float = 0.3


def build_schedule(
    daily_rate_ug_day: float,
    route: RouteParams,
    period_hours: float,
    injection_ug: Optional[float] = None,
    injection_interval_hours: Optional[float] = None,
) -> DosingSchedule:
    """Construct a :class:`DosingSchedule` for a regimen over ``period_hours``.

    * transdermal (``input_kind == 'zero_order'``): continuous input at the
      daily rate, spread across the day.
    * injectable / oral (``input_kind == 'first_order'``): discrete depot events
      every ``injection_interval_hours`` (oral defaults to daily; injectables to
      weekly unless overridden).
    """
    if daily_rate_ug_day < 0:
        raise ValueError("daily_rate_ug_day must be >= 0")

    # Calibrated delivery factor: folds bioavailability/first-pass/ester mass
    # into the *effective* delivered dose for this route.
    df = route.delivery_factor

    if route.input_kind == "zero_order":
        rate = daily_rate_ug_day * df / 24.0  # µg/h effective systemic flux
        return DosingSchedule(continuous_rate_ug_h=rate, k_abs_h=0.0)

    # first-order routes: discrete depot events.
    # Injections require an explicit per-event dose (injection_ug). Oral is
    # dosed daily, so a daily_rate is split into one event per day.
    if route.name in ("im_injection", "subq_injection"):
        if injection_ug is None:
            raise ValueError(
                f"Route {route.name!r} needs injection_ug (a per-event dose)."
            )
        interval = injection_interval_hours or 24.0 * 7.0  # default weekly
        per_event = injection_ug * df
    elif daily_rate_ug_day > 0:
        # oral: one discrete dose per day
        interval = 24.0
        per_event = daily_rate_ug_day * df
    else:
        raise ValueError(
            f"Route {route.name!r} needs injection_ug or a positive daily_rate."
        )

    events = [
        DoseEvent(t_hours=float(t), amount_ug=per_event)
        for t in np.arange(0.0, period_hours, interval)
    ]
    return DosingSchedule(
        continuous_rate_ug_h=0.0,
        events=events,
        k_abs_h=route.absorption_rate_h,
    )


def simulate_pk(
    schedule: DosingSchedule,
    pk: EstradiolPK,
    period_hours: float,
    dt_hours: float = 0.1,
    route: str = "transdermal",
    baseline_e2_pgml: Optional[float] = None,
) -> np.ndarray:
    """Integrate the estradiol PK model; return concentration (pg/mL) over time.

    Explicit Euler on a fine grid (``dt_hours``, default 0.1 h). The grid is
    fine enough that Euler error is negligible versus the model's inherent
    simplification. Starts at baseline estradiol.

    ``baseline_e2_pgml`` overrides the population endogenous baseline
    (intersex packs: e.g. Turner or post-gonadectomy CAIS sit near-adrenal
    levels with no gonadal production).
    """
    n = int(round(period_hours / dt_hours)) + 1
    A = np.empty(n)   # body mass (µg)
    D = np.empty(n)   # depot mass (µg)

    # Initial condition: baseline estradiol mass (pg/mL * Vd(L) = µg).
    A[0] = BASELINE_E2_PGML * pk.v_d_l
    D[0] = 0.0

    # Effective clearance (oral gets extra first-pass).
    cl = pk.effective_clearance(route)

    # Endogenous production rate calibrated so zero input -> baseline steady state.
    baseline_C = BASELINE_E2_PGML if baseline_e2_pgml is None else float(baseline_e2_pgml)
    production = cl * baseline_C  # µg/h (balances clearance at baseline C)

    # Dose->level calibration: scale the delivered dose (already delivery_factor-
    # adjusted in build_schedule) so steady-state E2 = baseline +
    # steady_state_factor * delivered_daily_dose.
    scale = pk.input_scale()

    # Precompute bolus events by grid index (already delivery_factor-adjusted;
    # now apply the dose->level calibration scale).
    event_map: Dict[int, float] = {}
    for ev in schedule.events:
        idx = int(round(ev.t_hours / dt_hours))
        if idx < n:
            event_map[idx] = event_map.get(idx, 0.0) + ev.amount_ug * scale

    for i in range(1, n):
        C = A[i - 1] / pk.v_d_l
        elimination = cl * C * dt_hours
        absorption = schedule.k_abs_h * D[i - 1] * dt_hours
        body_input = schedule.continuous_rate_ug_h * scale * dt_hours
        bolus = event_map.get(i, 0.0)

        A[i] = A[i - 1] - elimination + absorption + body_input + production * dt_hours
        if A[i] < 0.0:
            A[i] = 0.0
        D[i] = D[i - 1] + bolus - absorption
        if D[i] < 0.0:
            D[i] = 0.0

    return A / pk.v_d_l
