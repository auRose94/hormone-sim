"""Tests for the trans-masculine personalized figure (AGENTS.md §5a item 3).

Covers the weighted band API on ``AndrogenPersonalizedResult`` (the mirror of
the lactation result's ``band()``) and the posterior-vs-prior figure written
by ``ingest-t --out``.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from hormone_sim import (
    AndrogenObservation,
    personalize_androgen,
    plot_personal_androgen,
)


def _obs(**kw):
    base = dict(weeks_on_t=6.0, dose_mg=100.0, interval_days=7.0,
                testosterone_ngdl=880.0, hct_pct=None)
    base.update(kw)
    return AndrogenObservation(**base)


@pytest.fixture(scope="module")
def result():
    return personalize_androgen(
        "testosterone_cypionate",
        [_obs(), _obs(weeks_on_t=10.0, testosterone_ngdl=None, hct_pct=45.5)],
        n=60, seed=3,
    )


# --- weighted bands -----------------------------------------------------------------

def test_band_ordering_and_shape(result):
    lo, med, hi = result.band(0.8, which="t")
    assert lo.shape == med.shape == hi.shape == result.weeks.shape
    assert np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))
    assert np.all(lo <= med) and np.all(med <= hi)
    hlo, hmed, hhi = result.band(0.8, which="hct")
    assert np.all(hlo <= hmed) and np.all(hmed <= hhi)


def test_band_q_extremes_contain_interior(result):
    lo80, _, hi80 = result.band(0.8, which="t")
    lo98, _, hi98 = result.band(0.98, which="t")
    assert np.all(lo98 <= lo80) and np.all(hi98 >= hi80)


def test_band_rejects_unknown_curve(result):
    with pytest.raises(ValueError):
        result.band(0.8, which="voice")


def test_posterior_median_moves_toward_observation(result):
    """A constraining T observation moves the weighted median toward the
    measured value relative to the prior median. (Weighted bands SHIFT with
    the likelihood; they do not nest inside the prior band — the posterior
    can extend past the prior envelope on the side the data favors.)"""
    week_idx = int(round(6.0 / (result.weeks[1] - result.weeks[0])))
    measured = 880.0
    prior_med = float(np.median(result.t_ngdl[:, week_idx]))
    post_med = result.weighted_percentile(result.t_ngdl[:, week_idx], 50.0)
    assert abs(post_med - measured) < abs(prior_med - measured)


# --- figure --------------------------------------------------------------------------

def test_plot_personal_androgen_writes_png(result, tmp_path):
    out = tmp_path / "t_posterior.png"
    fig = plot_personal_androgen(
        result,
        [_obs(), _obs(weeks_on_t=10.0, testosterone_ngdl=None, hct_pct=45.5)],
        out_path=str(out))
    assert out.exists() and out.stat().st_size > 0
    plt.close(fig)


def test_plot_personal_androgen_no_observations_on_axis(result, tmp_path):
    """Programmatic use with no overlay points must not crash."""
    out = tmp_path / "t_posterior_bare.png"
    fig = plot_personal_androgen(result, [], out_path=str(out))
    assert out.exists()
    plt.close(fig)


def test_plot_personal_androgen_returns_fig_without_path(result):
    fig = plot_personal_androgen(result, [])
    assert fig is not None
    plt.close(fig)
