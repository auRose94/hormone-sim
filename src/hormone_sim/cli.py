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
from .compounds import (
    PHYTOESTROGENS,
    get_compound,
    hypothetical,
    MW_ESTRADIOL,
)
from .estrogenic_activity import (
    ee_of,
    required_daily_rate_ug_day,
    design_summary,
    TargetProfile,
)
from .breast_development import simulate_breast_development


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


def cmd_compound(args: argparse.Namespace) -> int:
    """List known compounds, or report one compound's potency + EE at a dose."""
    if args.name:
        c = get_compound(args.name)
        print(f"Compound: {c.name}")
        print(f"  REA (vs E2=1):        {c.rea:.3g}")
        print(f"  MW:                   {c.mw:g} g/mol")
        print(f"  free fraction:        {c.free_fraction:g}")
        print(f"  potency (per pg/mL):  {c.potency:.3g}")
        print(f"  calibrated:           {c.calibrated}")
        print(f"  note:                 {c.note}")
        print(f"  source:               {c.citation}")
        # What EE does a given plasma level produce?
        if args.plasma_pgml is not None:
            res = ee_of(c, np.array([args.plasma_pgml]))
            print(f"\n  At {args.plasma_pgml:g} pg/mL plasma -> EE {res.steady_state_ee_pgml():.4g} pg/mL")
        return 0

    print("Known estrogenic compounds (REA = relative estrogenic activity vs estradiol):")
    print(f"{'name':>14} | {'REA':>9} | {'MW':>8} | {'potency/pg':>11}")
    print("-" * 52)
    for name in sorted(PHYTOESTROGENS):
        c = PHYTOESTROGENS[name]
        print(f"{name:>14} | {c.rea:>9.1e} | {c.mw:>8.1f} | {c.potency:>11.2e}")
    print("\nTip: 'hormone-sim compound <name>' for details; use --hypothetical to")
    print("explore a hypothetical molecule. See docs/references.md for provenance.")
    return 0


def cmd_hypothetical(args: argparse.Namespace) -> int:
    """Inverse design for a hypothetical molecule (care yet to be discovered)."""
    c = hypothetical(
        name=args.name,
        rea=args.rea,
        mw=args.mw,
        free_fraction=args.free_fraction,
        note="hypothetical (user-specified target product profile)",
    )
    target = TargetProfile(ee_target_pgml=args.target_ee)
    print(design_summary(c, target, routes=["transdermal", "oral", "im_injection"]))
    print(f"\n  (REA={c.rea:.3g}, MW={c.mw:g}, potency={c.potency:.3g}; "
          "linear-PK property-space estimate, not a discovered drug.)")
    return 0


def cmd_breast(args: argparse.Namespace) -> int:
    """Simulate population-typical breast development from an EE exposure."""
    # Build an EE curve. If a regimen is given, use its E2 as EE; otherwise a
    # flat EE at the target level.
    if args.regimen:
        r = Regimen.from_yaml(args.regimen)
        sim = simulate(r, weeks=args.weeks)
        ee_curve = sim.e2_pgml
        label = r.name
    else:
        # flat EE exposure (steady state) for a target level.
        ee_curve = np.full(int(args.weeks / 0.25) + 1, args.flat_ee)
        label = f"flat EE {args.flat_ee:g} pg/mL"
    res = simulate_breast_development(ee_curve, weeks=args.weeks)
    print(f"Breast development: {label}, {args.weeks:g} weeks")
    print(res.report())
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

    pc = sub.add_parser("compound", help="list compounds or report one's potency/EE")
    pc.add_argument("--name", help="compound name (e.g. genistein); omit to list all")
    pc.add_argument("--plasma-pgml", type=float,
                    help="plasma concentration to compute EE at (for --name)")
    pc.set_defaults(func=cmd_compound)

    ph = sub.add_parser("hypothetical", help="inverse-design a hypothetical molecule")
    ph.add_argument("--name", required=True, help="label for the hypothetical compound")
    ph.add_argument("--rea", type=float, required=True,
                    help="relative estrogenic activity vs estradiol (=1)")
    ph.add_argument("--mw", type=float, default=MW_ESTRADIOL,
                    help="molecular weight (g/mol), default = estradiol")
    ph.add_argument("--free-fraction", type=float, default=0.02,
                    help="plasma free fraction, default = estradiol (0.02)")
    ph.add_argument("--target-ee", type=float, default=150.0,
                    help="target estradiol-equivalent (pg/mL)")
    ph.set_defaults(func=cmd_hypothetical)

    pb = sub.add_parser("breast", help="simulate population-typical breast development")
    pb.add_argument("--regimen", help="regimen YAML; use its E2 as the EE exposure")
    pb.add_argument("--flat-ee", type=float, default=150.0,
                    help="flat EE exposure (pg/mL) if no --regimen given")
    pb.add_argument("--weeks", type=float, default=104.0, help="horizon (default 2 yr)")
    pb.set_defaults(func=cmd_breast)

    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
