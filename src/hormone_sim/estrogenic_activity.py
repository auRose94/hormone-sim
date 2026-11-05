"""Estrogenic activity: folding a compound's PK into estradiol-equivalents.

The central quantity is the **estradiol-equivalent (EE)** concentration:

    EE(t) = C_compound(t) [pg/mL] * potency(compound)

where ``potency`` (from :mod:`hormone_sim.compounds`) is the receptor-weighted
relative potency per pg/mL vs estradiol. An EE of 100 pg/mL means "the total
estrogenic receptor activation equals that of 100 pg/mL estradiol," regardless
of whether it came from estradiol itself, genistein, or a hypothetical drug.

This is what lets us compare a breast-growth cream, a soy supplement, and
estradiol on one scale — and what makes **inverse design** possible: given a
target EE (e.g. the 100–200 pg/mL feminizing window), what PK (rate,
potency, delivery) does a candidate molecule need?

Epistemic status: REA/potency values are assay-dependent literature central
estimates with wide uncertainty. Treat EE as an order-of-magnitude, not a
precise prediction. This is a research model, not dosing advice.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from .compounds import Compound, get_compound, MW_ESTRADIOL


@dataclass
class EstrogenicResult:
    """The estrogenic activity of a compound over time.

    ``ee_pgml`` is the estradiol-equivalent concentration curve. ``e2_curve`` is
    the compound's own plasma concentration curve (pg/mL).
    """

    compound_name: str
    potency: float
    hours: np.ndarray
    e2_curve_pgml: np.ndarray   # the compound's own plasma concentration
    ee_curve_pgml: np.ndarray   # estradiol-equivalents

    def steady_state_ee_pgml(self, tail_frac: float = 0.25) -> float:
        n = self.ee_curve_pgml.size
        cut = int(n * (1.0 - tail_frac))
        return float(self.ee_curve_pgml[cut:].mean())

    def e2_curve_pmol(self) -> np.ndarray:
        return self.e2_curve_pgml * 3.671

    def report(self) -> str:
        ee = self.steady_state_ee_pgml()
        return (
            f"Compound: {self.compound_name}\n"
            f"  potency (per pg/mL vs E2): {self.potency:.3g}\n"
            f"  steady-state compound level: {np.mean(self.e2_curve_pgml):.3g} pg/mL\n"
            f"  steady-state EE: {ee:.4g} pg/mL "
            f"({ee * 3.671:.3g} pmol/L estradiol-equivalent)"
        )


def ee_of(compound: Compound, e2_curve_pgml: np.ndarray) -> EstrogenicResult:
    """Compute estradiol-equivalents for a compound's plasma concentration curve."""
    ee = np.asarray(e2_curve_pgml, dtype=float) * compound.potency
    return EstrogenicResult(
        compound_name=compound.name,
        potency=compound.potency,
        hours=None,  # caller attaches time axis
        e2_curve_pgml=np.asarray(e2_curve_pgml, dtype=float),
        ee_curve_pgml=ee,
    )


# --- inverse design: what does a hypothetical molecule need? -----------------
@dataclass
class TargetProfile:
    """A target feminizing estrogenic exposure.

    ``ee_target_pgml`` is the desired steady-state estradiol-equivalent (e.g. 150
    pg/mL for the center of the feminizing window). ``ee_window`` is the
    acceptable band.
    """

    ee_target_pgml: float = 150.0
    ee_window: tuple = (100.0, 200.0)


def required_daily_rate_ug_day(
    compound: Compound,
    target_ee_pgml: float,
    route: str = "oral",
) -> float:
    """Inverse design: daily dose (µg/day) of ``compound`` needed to reach
    ``target_ee_pgml`` estradiol-equivalents, given a route's delivery factor.

    Uses the estradiol model's dose→level mapping as the reference: estradiol
    reaches ~150 pg/mL at 50 µg/day transdermal (delivery_factor 1.0). For a
    compound, the delivered dose must be scaled up by 1/potency to achieve the
    same EE. This is the *property-space* answer — the dose a hypothetical
    molecule with this potency would need.

    Caveat: assumes linear PK and that the compound's plasma curve shape matches
    estradiol's for the route (reasonable first approximation; see module docstring).
    """
    from .params import ROUTES, EstradiolPK

    ref_ee = 150.0          # estradiol reference EE at 50 ug/day transdermal
    ref_dose = 50.0         # µg/day
    df = ROUTES[route].delivery_factor
    # estradiol EE per delivered µg/day (transdermal): 150 pg/mL / 50 µg/day = 3.
    e2_ee_per_ug = ref_ee / ref_dose
    # compound EE per delivered µg/day = e2_ee_per_ug * potency * df/df_ref.
    # df_ref for transdermal = 1.0; for other routes it's the route's factor.
    # We want delivered_ug * e2_ee_per_ug * potency * df = target_ee.
    if compound.potency <= 0:
        raise ValueError("compound potency must be > 0")
    required = target_ee_pgml / (e2_ee_per_ug * compound.potency * df)
    return required


def design_summary(compound: Compound, target: TargetProfile, routes: Optional[List[str]] = None) -> str:
    """Human-readable inverse-design summary across routes."""
    routes = routes or ["transdermal", "oral", "im_injection"]
    lines = [
        f"Inverse design: {compound.name}",
        f"  REA={compound.rea:g}, MW={compound.mw:g}, potency={compound.potency:.3g}",
        f"  Target EE: {target.ee_target_pgml:g} pg/mL "
        f"(window {target.ee_window[0]:.0f}-{target.ee_window[1]:.0f})",
        f"  Required daily dose (µg/day) by route:",
    ]
    for route in routes:
        try:
            d = required_daily_rate_ug_day(compound, target.ee_target_pgml, route)
            lines.append(f"    {route:>15}: {format_dose(d)}")
        except Exception as e:  # noqa: BLE001
            lines.append(f"    {route:>15}: n/a ({e})")
    return "\n".join(lines)


def format_dose(ug_day: float) -> str:
    """Format a daily dose (µg/day) into a human-friendly unit."""
    if ug_day < 1.0:
        return f"{ug_day * 1000:.2f} µg/day"
    if ug_day < 1000.0:
        return f"{ug_day:.1f} µg/day"
    return f"{ug_day / 1000:.2f} mg/day"
