"""Command-line interface for hormone-sim.

Subcommands:
  run    Simulate one regimen (YAML) and report + plot.
  scan   Sweep a dose range for a route; report steady-state E2 vs dose.
  pop    Monte-Carlo population simulation; report target-achievement rates.

This CLI is a research/education tool. It produces population-typical curves,
not individualized dosing.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from .regimens import Regimen
from .sim import simulate
from .population import run_population
from .plot import plot_single, plot_population
from .params import (
    ROUTES,
    e2_pgml_to_pmol,
    t_ngdl_to_nmol,
    TARGET_E2_PGML,
    TARGET_T_NGDL,
)


def _safe_slug(name: str) -> str:
    """Turn a regimen name into a filesystem-safe slug (no path separators)."""
    return "".join(c if c.isalnum() or c in "_-." else "_" for c in name)


def _build_regimen(
    route: str,
    dose_ug_day: float,
    anti_androgen: dict | None,
    name: str,
    injection_ug: float | None = None,
    injection_interval: float | None = None,
) -> Regimen:
    from .regimens import EstradiolSpec
    return Regimen(
        name=name,
        population="trans_feminine",
        estradiol=EstradiolSpec(
            route=route,
            daily_rate_ug_day=dose_ug_day,
            injection_ug=injection_ug,
            injection_interval_hours=injection_interval,
        ),
        anti_androgen=anti_androgen or {},
    )


def cmd_run(args: argparse.Namespace) -> int:
    regimen = Regimen.from_yaml(args.regimen)
    sim = simulate(regimen, weeks=args.weeks)
    print(sim.report())
    if args.out:
        out = Path(args.out)
        slug = _safe_slug(regimen.name)
        png = out / f"{slug}.png"
        plot_single(sim, out_path=str(png))
        csv = out / f"{slug}.csv"
        import csv as _csv
        with csv.open("w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["hours", "weeks", "e2_pgml", "e2_pmol", "t_ngdl", "t_nmol"])
            for i in range(sim.hours.size):
                w.writerow([
                    round(sim.hours[i], 2),
                    round(sim.hours[i] / 168.0, 3),
                    round(sim.e2_pgml[i], 3),
                    round(e2_pgml_to_pmol(sim.e2_pgml[i]), 1),
                    round(sim.t_ngdl[i], 3),
                    round(t_ngdl_to_nmol(sim.t_ngdl[i]), 3),
                ])
        print(f"\nWrote: {png}")
        print(f"Wrote: {csv}")
    return 0


def cmd_scan(args: argparse.Namespace) -> int:
    print(f"Route: {args.route}\n{'Dose (ug/day)':>14} | "
          f"{'SS E2 pg/mL':>12} | {'E2 pmol/L':>10} | {'E2 target':>9} | "
          f"{'SS T ng/dL':>11} | {'T supp':>7}")
    print("-" * 80)
    rows = []
    for dose in args.doses:
        aa = {"spironolactone": 1.0} if args.with_antiandrogen else {}
        regimen = _build_regimen(args.route, dose, aa, f"{args.route}_{dose}ug")
        sim = simulate(regimen, weeks=args.weeks)
        e2 = sim.steady_state_e2_pgml()
        t = sim.steady_state_t_ngdl()
        in_target = "yes" if TARGET_E2_PGML[0] <= e2 <= TARGET_E2_PGML[1] else "no"
        suppressed = "yes" if t <= TARGET_T_NGDL[1] else "no"
        print(f"{dose:>14g} | {e2:>12.1f} | {e2_pgml_to_pmol(e2):>10.0f} | "
              f"{in_target:>9} | {t:>11.1f} | {suppressed:>7}")
        rows.append({
            "route": args.route, "dose_ug_day": dose,
            "ss_e2_pgml": round(e2, 2), "ss_e2_pmol": round(e2_pgml_to_pmol(e2), 1),
            "ss_t_ngdl": round(t, 2), "e2_in_target": in_target,
            "t_suppressed": suppressed,
        })
    if args.out:
        import csv
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote: {out}")
    return 0


def cmd_pop(args: argparse.Namespace) -> int:
    from .regimens import EstradiolSpec
    regimen = _build_regimen(args.route, args.dose, {"spironolactone": 1.0},
                             f"pop_{args.route}_{args.dose}ug")
    pop = run_population(
        lambda: regimen, weeks=args.weeks, n=args.n, seed=args.seed
    )
    print(f"Population: {args.n} subjects, {args.weeks:g} weeks, "
          f"route={args.route}, dose={args.dose:g} ug/day")
    print(f"  Mean steady-state E2: {pop.steady_state_e2_mean():.1f} pg/mL "
          f"({e2_pgml_to_pmol(pop.steady_state_e2_mean()):.0f} pmol/L)")
    print(f"  Fraction E2 in target ({TARGET_E2_PGML[0]:.0f}-{TARGET_E2_PGML[1]:.0f} "
          f"pg/mL): {pop.fraction_e2_in_target()*100:.1f}%")
    print(f"  Fraction T suppressed (<{TARGET_T_NGDL[1]:.0f} ng/dL): "
          f"{pop.fraction_t_suppressed()*100:.1f}%")
    if args.out:
        out = Path(args.out)
        png = out / f"pop_{args.route}_{args.dose}ug.png"
        plot_population(pop, out_path=str(png))
        print(f"Wrote: {png}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hormone-sim",
        description="Population PK/PD simulation of gender-affirming hormone therapy. "
        "Research/education tool — population curves, not individualized dosing.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    pr = sub.add_parser("run", help="simulate one regimen (YAML) and report/plot")
    pr.add_argument("--regimen", required=True, help="path to regimen YAML")
    pr.add_argument("--weeks", type=float, default=12.0)
    pr.add_argument("--out", help="output directory for PNG + CSV")
    pr.set_defaults(func=cmd_run)

    ps = sub.add_parser("scan", help="sweep a dose range for a route")
    ps.add_argument("--route", required=True, choices=sorted(ROUTES.keys()))
    ps.add_argument("--doses", nargs="+", type=float, required=True)
    ps.add_argument("--weeks", type=float, default=12.0)
    ps.add_argument("--with-antiandrogen", action="store_true",
                    help="add spironolactone 1.0 (relative units)")
    ps.add_argument("--out", help="write scan results CSV")
    ps.set_defaults(func=cmd_scan)

    pp = sub.add_parser("pop", help="Monte-Carlo population simulation")
    pp.add_argument("--route", required=True, choices=sorted(ROUTES.keys()))
    pp.add_argument("--dose", type=float, required=True, help="ug/day")
    pp.add_argument("--weeks", type=float, default=12.0)
    pp.add_argument("--n", type=int, default=200)
    pp.add_argument("--seed", type=int, default=42)
    pp.add_argument("--out", help="output directory for PNG")
    pp.set_defaults(func=cmd_pop)

    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
