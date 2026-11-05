#!/usr/bin/env python3
"""Reproduce the study figures/scans from the README.

Runs a dose-response scan for each route and a population band plot for the
first-line transdermal regimen, writing PNG + CSV into ``outputs/``.

This is the "one script reproduces the study" driver. Run it after ``pip install``.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make the package importable when run from the repo root without install.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from hormone_sim import simulate, run_population  # noqa: E402
from hormone_sim.regimens import Regimen, EstradiolSpec  # noqa: E402
from hormone_sim.plot import plot_single, plot_population  # noqa: E402
from hormone_sim.params import e2_pgml_to_pmol, t_ngdl_to_nmol, TARGET_E2_PGML  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(exist_ok=True)

ROUTES = ["transdermal", "oral", "im_injection"]
# transdermal & oral doses are ug/day; injection doses are ug PER INJECTION (weekly).
DOSES = {
    "transdermal": [25, 50, 75, 100],
    "oral": [1000, 2000, 4000, 8000],
    "im_injection": [10000, 20000, 25000],  # 10, 20, 25 mg weekly
}


def scan_route(route: str, doses: list[int], weeks: float = 12.0) -> None:
    print(f"\n=== Dose-response: {route} ===")
    print(f"{'Dose':>10} | {'SS E2 pg/mL':>12} | {'E2 pmol/L':>10} | "
          f"{'E2 target':>9} | {'SS T ng/dL':>11} | {'T supp':>7}")
    for dose in doses:
        if route == "im_injection":
            # dose is per-injection (ug), weekly.
            spec = EstradiolSpec(
                route=route, daily_rate_ug_day=0.0,
                injection_ug=dose, injection_interval_hours=168.0,
            )
        else:
            spec = EstradiolSpec(
                route=route, daily_rate_ug_day=dose,
            )
        r = Regimen(name=f"{route}_{dose}", population="trans_feminine",
                    estradiol=spec, anti_androgen={"spironolactone": 1.0})
        sim = simulate(r, weeks=weeks)
        e2 = sim.steady_state_e2_pgml()
        t = sim.steady_state_t_ngdl()
        in_target = "yes" if TARGET_E2_PGML[0] <= e2 <= TARGET_E2_PGML[1] else "no"
        print(f"{dose:>10g} | {e2:>12.1f} | {e2_pgml_to_pmol(e2):>10.0f} | "
              f"{in_target:>9} | {t:>11.1f} | {'yes' if t <= 50 else 'no':>7}")
        plot_single(sim, out_path=str(OUT / f"{route}_{dose}.png"))


def population_plot() -> None:
    print("\n=== Population bands: transdermal 50 ug/day ===")
    spec = EstradiolSpec(route="transdermal", daily_rate_ug_day=50)
    r = Regimen(name="pop_transdermal_50", population="trans_feminine",
                estradiol=spec, anti_androgen={"spironolactone": 1.0})
    pop = run_population(lambda: r, weeks=12, n=200, seed=42)
    plot_population(pop, out_path=str(OUT / "population_bands.png"))
    print(f"  E2 in target: {pop.fraction_e2_in_target()*100:.1f}%")
    print(f"  T suppressed: {pop.fraction_t_suppressed()*100:.1f}%")


if __name__ == "__main__":
    for route in ROUTES:
        doses = DOSES.get(route, [50])
        scan_route(route, doses)
    population_plot()
    print(f"\nDone. Outputs in {OUT}")
