"""Testosterone suppression pharmacodynamics.

Model the hypothalamic-pituitary-gonadal (HPG) axis suppression driven by
estradiol negative feedback, plus additive pharmacologic anti-androgen effects.

    T(t) = T_baseline * (1 - S_e2(t)) * (1 - S_aa)

* ``S_e2`` — estradiol-driven HPG suppression, a Hill function of the *recent
  average* estradiol level. Using a moving average (not instantaneous) reflects
  that gonadal suppression tracks sustained exposure, not momentary peaks.
* ``S_aa`` — anti-androgen effect, the combined fractional suppression from the
  agents in the regimen (each contributes Emax-weighted, capped so the total
  stays in [0, 1]).

This is a phenomenological, population-typical model. It captures the *direction
and magnitude* of suppression seen clinically; it is not a mechanistic
gonadotropin model and is not patient-specific.
"""

from __future__ import annotations

import numpy as np

from .params import AndrogenPD, BASELINE_T_NGDL


def estradiol_suppression(avg_e2_pgml: np.ndarray, pd: AndrogenPD) -> np.ndarray:
    """Hill-function suppression fraction S_e2 from average estradiol (pg/mL).

    S_e2 = C^n / (EC50^n + C^n),  with Emax = 1 (complete suppression).
    At C=0 -> S=0 (no suppression); at high C -> S->1.
    """
    ec50 = pd.e2_ec50_pgml
    n = pd.hill_n
    c = np.asarray(avg_e2_pgml, dtype=float)
    return (c ** n) / (ec50 ** n + c ** n)


def anti_androgen_suppression(agents: dict) -> float:
    """Combined fractional anti-androgen suppression from a regimen.

    ``agents`` maps agent name -> dose (relative units; see regimens). The
    contribution of each agent is its Emax scaled by a simple dose-response
    (dose / (dose + 1), normalized so dose=1 ~ typical clinical dose). Combined
    multiplicatively on the *un-suppressed* fraction so the total is bounded.
    """
    pd = AndrogenPD()
    emax = pd.anti_androgen_emax
    unsuppressed = 1.0
    for name, dose in (agents or {}).items():
        if name not in emax:
            raise KeyError(f"Unknown anti-androgen: {name!r}")
        dose = float(dose)
        # dose-response: fraction of Emax achieved at this (relative) dose.
        frac = dose / (dose + 1.0)
        unsuppressed *= (1.0 - emax[name] * frac)
    return max(0.0, min(1.0, 1.0 - unsuppressed))


def suppress_testosterone(
    e2_curve_pgml: np.ndarray,
    pd: AndrogenPD,
    window_hours: float = 72.0,
    dt_hours: float = 0.1,
    baseline_t_ngdl: float = BASELINE_T_NGDL,
    agents: dict | None = None,
) -> np.ndarray:
    """Return testosterone (ng/dL) over time given an estradiol curve.

    Uses a trailing moving average of estradiol over ``window_hours`` as the
    driver of HPG suppression, and applies the anti-androgen effect as a
    constant multiplier.
    """
    e2 = np.asarray(e2_curve_pgml, dtype=float)
    n = e2.size
    t = np.empty(n)

    # Precompute trailing moving average of estradiol.
    win = max(1, int(round(window_hours / dt_hours)))
    kernel = np.ones(win) / win
    # 'same' with trailing alignment: shift so average is over [i-win, i].
    conv = np.convolve(e2, kernel, mode="full")
    # conv index for trailing window ending at i is i + win - 1.
    avg_e2 = np.empty(n)
    for i in range(n):
        lo = max(0, i - win + 1)
        avg_e2[i] = e2[lo : i + 1].mean()

    s_e2 = estradiol_suppression(avg_e2, pd)
    s_aa = anti_androgen_suppression(agents)

    # T = baseline * (1 - S_e2) * (1 - S_aa)
    t = baseline_t_ngdl * (1.0 - s_e2) * (1.0 - s_aa)
    return np.clip(t, 0.0, None)
