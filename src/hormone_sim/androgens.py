"""Androgen pharmacokinetics: testosterone agents for masculinizing GAHT.

Same architecture as the estradiol/galactagogue layer: registry-style PK,
one-compartment core, route differences living in the input function.

The headline structural fact for injectable esters is the **peak–trough
swing**: cypionate/enanthate depots release over days-to-weeks, so a
200 mg/2-week schedule and a 100 mg/weekly schedule deliver the *same
average* dose but very different curves — the two-week schedule rides peaks
well above the male reference range and spends the back half of the interval
in decline. Patient-reported symptom patterns (irritability, fatigue) are
commonly tied to troughs; the model tracks peak, trough, and mean separately
so that pattern can be reasoned about. The gel, by contrast, is near-flat —
the masculinization analog of "transdermal estradiol for stable levels."

Units: serum testosterone in **ng/dL** (clinical convention). Internally the
integrator works in µg/L (1 ng/dL = 0.01 µg/L... to be exact: 1 ng/dL =
10 pg/mL = 0.01 µg/L × — no: 1 ng/dL = 1 ng / 0.1 L = 10 ng/L = 0.01 µg/L ×
10 → 1 ng/dL = 0.01 µg/L is wrong by 10×; 1 ng/dL = 10 ng/L = 0.01 µg/L.
Check: 500 ng/dL = 5 µg/L. Yes: ng/dL = 100 × µg/L. The converter lives in
:func:`ugl_to_ngdl`.)

Baselines are **trans-masculine pre-T**: ovarian testosterone ~15–70 ng/dL
(mean ~30), which the HPG suppression in :mod:`hormone_sim.masculinization`
then drives down as exogenous T rises — the reversed-direction mirror of the
estradiol model's testicular suppression.

Evidence anchors: Endocrine Society CPG (Grimley et al. 2017) target range
300–1000 ng/dL trough-adjusted; ester half-lives from standard references
(cypionate ~8 d; enanthate ~4.5–10 d; undecanoate IM ~20–40 d, wide).
Oral testosterone undecanoate is deliberately NOT registered — its
food-dependent, highly variable absorption would need data this package
doesn't have; excluding it is the honest move (same rule as the herbal
galactagogues).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np


# --- Unit helpers -------------------------------------------------------------
# 1 ng/dL = 1 ng / 0.1 L = 10 ng/L; 1 µg/L = 1000 ng/L = 100 ng/dL.
# Check: 500 ng/dL = 5 µg/L. ✓
def ugl_to_ngdl(ug_l: float) -> float:
    """µg/L -> ng/dL (1 µg/L = 100 ng/dL)."""
    return ug_l * 100.0


def ngdl_to_ugl(ngdl: float) -> float:
    """ng/dL -> µg/L."""
    return ngdl / 100.0


# --- Pre-T baselines (population-typical, trans-masculine) ---------------------
BASELINE_T_NGDL = 30.0        # ovarian range ~15-70
BASELINE_E2_PGML = 60.0       # cycling mean; phases span ~30-400 (flagged wide)


@dataclass(frozen=True)
class Androgen:
    """A testosterone agent's PK.

    ``input_kind``: ``"depot"`` (IM/SubQ esters — first-order release from
    the depot, like the estradiol injectables) or ``"zero_order"`` (gel —
    near-constant daily flux).

    ``half_life_h`` is the *molecule's* apparent elimination half-life (fast);
    ``depot_half_life_days`` is the ester release half-life (slow, and what
    dominates the observed curve for injectables).
    """

    name: str
    input_kind: str
    half_life_h: float                 # testosterone molecule (apparent)
    vd_l_kg: float
    bioavailability: float             # fraction of nominal dose reaching systemic
    depot_half_life_days: Optional[float]  # esters only
    typical_dose_mg: float             # per administration
    typical_interval_days: float
    citation: str
    note: str = ""

    def clearance_l_per_h(self, body_weight_kg: float = 70.0) -> float:
        return (0.6931471805599453 * self.vd_l_kg * body_weight_kg
                / self.half_life_h)

    def depot_rate_h(self) -> float:
        """First-order depot release constant (1/h) for esters."""
        if self.depot_half_life_days is None:
            return 0.0
        return 0.6931471805599453 / (self.depot_half_life_days * 24.0)


# --- Registry --------------------------------------------------------------------
ANDROGENS: Dict[str, Androgen] = {
    "testosterone_cypionate": Androgen(
        name="testosterone_cypionate",
        input_kind="depot",
        half_life_h=2.0,          # T molecule, apparent; ester release dominates
        vd_l_kg=2.0,
        bioavailability=0.90,
        depot_half_life_days=8.0,
        typical_dose_mg=100.0,
        typical_interval_days=7.0,
        citation="T cypionate: depot release t1/2 ~8 d; commonly 100 mg/wk "
        "SubQ or 200 mg/2wk IM (Endocrine Society CPG, Grimley et al. 2017; "
        "WPATH SOC-8 2022). T molecule apparent t1/2 and Vd calibrated "
        "(fast relative to depot).",
        note="The peak-trough agent: weekly and biweekly schedules deliver "
        "the same average but very different swings. Troughs are where "
        "patient-reported symptom patterns commonly land.",
    ),
    "testosterone_enanthate": Androgen(
        name="testosterone_enanthate",
        input_kind="depot",
        half_life_h=2.0,
        vd_l_kg=2.0,
        bioavailability=0.90,
        depot_half_life_days=7.0,
        typical_dose_mg=100.0,
        typical_interval_days=7.0,
        citation="T enanthate: depot t1/2 ~4.5-10 d (midpoint 7 used, "
        "flagged wide); dosing as for cypionate (Grimley et al. 2017).",
        note="Slightly faster release than cypionate; near-interchangeable "
        "in practice at the same schedule.",
    ),
    "testosterone_undecanoate_im": Androgen(
        name="testosterone_undecanoate_im",
        input_kind="depot",
        half_life_h=2.0,
        vd_l_kg=2.0,
        bioavailability=0.90,
        depot_half_life_days=30.0,   # wide literature spread; flagged
        typical_dose_mg=750.0,
        typical_interval_days=70.0,  # ~10 weeks after loading
        citation="T undecanoate IM: very long depot (loading then ~10-wk "
        "intervals); t1/2 ~20-40 d in castor oil (WIDE — calibrated midpoint "
        "used). Grimley et al. 2017.",
        note="Flat-at-scale: months-long release. Clinic-administered; "
        "large single doses.",
    ),
    "testosterone_gel": Androgen(
        name="testosterone_gel",
        input_kind="zero_order",
        half_life_h=2.0,
        vd_l_kg=2.0,
        bioavailability=0.10,        # ~10% of nominal gel dose systemic
        depot_half_life_days=None,
        typical_dose_mg=50.0,        # nominal mg/day (e.g. 5 g of 1% gel)
        typical_interval_days=1.0,
        citation="T gel: ~10% of nominal dose absorbed systemically, "
        "daily application, near-flat levels (Grimley et al. 2017).",
        note="The 'stable levels' route: minimal swing, lower peak exposure, "
        "transfer precautions (skin-to-skin contact with others).",
    ),
}


def get_androgen(name: str) -> Androgen:
    if name not in ANDROGENS:
        raise KeyError(
            f"Unknown androgen {name!r}. Known: {sorted(ANDROGENS)}. "
            "Oral testosterone undecanoate is deliberately unregistered "
            "(food-dependent variable absorption; no honest anchor)."
        )
    return ANDROGENS[name]


# --- PK: coupled T curve with HPG-suppressed ovarian production ------------------
def simulate_androgen(
    a: Androgen,
    dose_mg: float,
    interval_days: float,
    weeks: float,
    dt_hours: float = 1.0,
    e2_ec50_ngdl: float = 170.0,
    e2_hill_n: float = 2.0,
    body_weight_kg: float = 70.0,
) -> Dict[str, np.ndarray]:
    """Simulate serum testosterone (ng/dL) under repeated dosing.

    The exogenous agent enters via depot release (esters) or continuous flux
    (gel). Endogenous **ovarian** production is suppressed dynamically by
    the testosterone level through HPG negative feedback (Hill function) —
    this is the reversed-direction mirror of the trans-feminine model, and
    it is why T monotherapy stops cycles.

    Returns dict with keys: ``hours``, ``t_ngdl``, ``e2_pgml``,
    ``ovarian_frac`` (0-1 surviving ovarian contribution).
    """
    cl = a.clearance_l_per_h(body_weight_kg)
    vd_ugl = a.vd_l_kg * body_weight_kg            # Vd in L; C(µg/L) = A(µg)/Vd

    # Ovarian production calibrated so zero-dose -> baseline T.
    base_c_ugl = ngdl_to_ugl(BASELINE_T_NGDL)
    p_ova_base = cl * base_c_ugl                    # µg/h

    n = int(round(weeks * 7 * 24 / dt_hours)) + 1
    hours = np.arange(n) * dt_hours
    A = np.zeros(n)
    D = np.zeros(n)

    k_el = 0.6931471805599453 / a.half_life_h
    k_dep = a.depot_rate_h() if a.input_kind == "depot" else 0.0

    # Dose events (depot) or continuous flux (gel). mg -> µg is ×1e3.
    interval_h = interval_days * 24.0
    dose_events: Dict[int, float] = {}
    if a.input_kind == "depot":
        t = 0.0
        while t < weeks * 7 * 24:
            dose_events[int(round(t / dt_hours))] = dose_mg * 1e3 * a.bioavailability  # µg
            t += interval_h
        gel_flux = 0.0
    else:
        gel_flux = dose_mg * 1e3 * a.bioavailability / 24.0  # µg/h

    # Suppression uses a trailing mean of T (days-scale, like the E2 model).
    suppress_window_h = 72.0
    win = max(1, int(round(suppress_window_h / dt_hours)))
    ec50_ugl = ngdl_to_ugl(e2_ec50_ngdl)

    ovarian_frac = np.empty(n)
    for i in range(1, n):
        c_ugl = A[i - 1] / vd_ugl
        # Trailing-mean suppression of ovarian production.
        lo = max(0, i - win)
        c_avg = A[lo:i].mean() / vd_ugl
        denom = ec50_ugl ** e2_hill_n + c_avg ** e2_hill_n
        s = (c_avg ** e2_hill_n / denom) if denom > 0 else 0.0
        p_ova = p_ova_base * (1.0 - s)
        ovarian_frac[i] = 1.0 - s

        elim = cl * c_ugl * dt_hours
        depot_release = k_dep * D[i - 1] * dt_hours
        bolus = dose_events.get(i, 0.0)

        A[i] = A[i - 1] - elim + depot_release + gel_flux * dt_hours + p_ova * dt_hours
        if A[i] < 0.0:
            A[i] = 0.0
        D[i] = D[i - 1] + bolus - depot_release
        if D[i] < 0.0:
            D[i] = 0.0
    ovarian_frac[0] = 1.0

    t_ngdl = ugl_to_ngdl(A / vd_ugl)
    # Estradiol: ovarian E2 falls with the same suppression; aromatization
    # from the (suppressed) total T adds ~0.05 pg/mL per ng/dL.
    e2 = BASELINE_E2_PGML * ovarian_frac + 0.05 * t_ngdl

    return {
        "hours": hours,
        "t_ngdl": t_ngdl,
        "e2_pgml": e2,
        "ovarian_frac": ovarian_frac,
    }
