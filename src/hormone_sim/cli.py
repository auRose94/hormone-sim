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
from .galactagogues import (
    GALACTAGOGUES,
    get_galactagogue,
    prolactin_response,
    prolactin_curve,
    list_galactagogues,
)
from .lactation import LactationProtocol, simulate_lactation
from .sensitivity import (
    dose_sweep_population,
    dose_sweep_table,
    simulate_lactation_population,
)
from .risk_exposure import exposure_report


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


def cmd_galactagogues(args: argparse.Namespace) -> int:
    """List pharmacologic galactagogues, or detail one + its prolactin response."""
    if args.name:
        g = get_galactagogue(args.name)
        print(f"Galactagogue: {g.name}")
        print(f"  half-life:       {g.half_life_h:g} h")
        print(f"  bioavailability: {g.bioavailability:g}")
        print(f"  Vd:              {g.vd_l_kg:g} L/kg")
        print(f"  prolactin Emax:  +{g.prolactin_emax_ng_ml:g} ng/mL")
        print(f"  prolactin ED50:  {g.prolactin_ed50_mg_day:g} mg/day")
        print(f"  US-approved:     {g.approved_us}")
        print(f"  safety:          {g.safety_note}")
        print(f"  source:          {g.citation}")
        print()
        print(f"  Prolactin (ng/mL) vs daily dose:")
        print(f"  {'10 mg/d':>9} | {'30 mg/d':>9} | {'60 mg/d':>9} | {'90 mg/d':>9}")
        row = [prolactin_response(g, d) for d in (10, 30, 60, 90)]
        print("  " + " | ".join(f"{v:>9.0f}" for v in row))
        print()
        print("  PK note: domperidone accumulates over ~3-5 days (t1/2 ~15 h);")
        print("  metoclopramide reaches steady state in ~1-2 days (t1/2 ~5.5 h).")
        return 0

    print(list_galactagogues())
    return 0


def cmd_lactation(args: argparse.Namespace) -> int:
    """Simulate a Newman-Goldfarb-style lactation induction protocol."""
    protocol = LactationProtocol(
        galactagogue=args.galactagogue,
        galactagogue_dose_mg_day=args.dose,
        prep_weeks=args.prep_weeks,
        progesterone_start_weeks=args.p_start_weeks,
        pumping_sessions_per_day=args.sessions,
    )
    result = simulate_lactation(
        protocol,
        total_weeks=args.prep_weeks + args.lactation_weeks,
        dt_weeks=0.25,
    )
    print(result.report())
    print()
    print("  Post-drop milk ramp (mL/day):")
    drop = protocol.prep_weeks
    for wk in (0, 1, 2, 4, 6, 8, 12, 26, 52):
        i = int(round((drop + wk) / 0.25))
        if i < result.milk_ml_day.size:
            print(f"    +{wk:>2} wk: {result.milk_ml_day[i]:>6.1f}")

    if getattr(args, "bands", False):
        from .sensitivity import simulate_lactation_population
        pop = simulate_lactation_population(
            protocol, n=args.n, seed=args.seed,
            total_weeks=args.prep_weeks + args.lactation_weeks,
        )
        print()
        print(pop.report())
        result_for_plot = pop
    else:
        result_for_plot = None

    if args.out:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from .galactagogues import PROLACTIN_BASELINE_NG_ML

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 9), sharex=True,
                                            constrained_layout=True)
        ax1.plot(result.weeks, result.ductal, color="#8e44ad", lw=1.6,
                 label="ductal (E2-driven)")
        ax1.plot(result.weeks, result.lobuloalveolar, color="#c0392b", lw=1.6,
                 label="lobuloalveolar (P-driven)")
        ax1.plot(result.weeks, result.priming, color="#555", lw=1.2, ls="--",
                 label="capacity = sqrt(ductal·LA)")
        ax1.plot(result.weeks, result.maturation, color="#e67e22", lw=1.6,
                 label="lactogenesis maturation")
        # Mark the progesterone window.
        p_on = protocol.p_start()
        if p_on < protocol.prep_weeks:
            ax1.axvspan(p_on, protocol.prep_weeks, color="#c0392b", alpha=0.08)
        ax1.set_ylabel("index (0-1)")
        ax1.legend(fontsize=7, loc="upper left")
        ax1.grid(alpha=0.3)

        ax2.plot(result.weeks, result.prolactin, color="#c0392b", lw=1.8)
        ax2.axhline(PROLACTIN_BASELINE_NG_ML, color="#999", ls=":", lw=0.8,
                    label="baseline")
        ax2.set_ylabel("Prolactin (ng/mL)")
        ax2.legend(fontsize=8)
        ax2.grid(alpha=0.3)

        if result_for_plot is not None:
            lo, med, hi = result_for_plot.band(0.8)
            ax3.fill_between(result.weeks, lo, hi, color="#2980b9", alpha=0.25,
                             label="80% interval")
            ax3.plot(result.weeks, med, color="#2980b9", lw=1.8, label="median")
            ax3.legend(fontsize=8)
        else:
            ax3.fill_between(result.weeks, 0, result.milk_ml_day,
                             color="#2980b9", alpha=0.35)
            ax3.plot(result.weeks, result.milk_ml_day, color="#2980b9", lw=1.8)
        ax3.set_ylabel("Milk (mL/day)")
        ax3.set_xlabel("Weeks")
        ax3.grid(alpha=0.3)

        drop = protocol.prep_weeks
        for ax in fig.axes:
            ax.axvline(drop, color="#666", ls="--", lw=0.8)

        title = (f"Lactation induction — {protocol.galactagogue} "
                 f"{protocol.galactagogue_dose_mg_day:g} mg/day, "
                 f"{protocol.prep_weeks:g} wk prep, "
                 f"{protocol.pumping_sessions_per_day:g} sessions/day")
        if result_for_plot is not None:
            title += f" [bands n={result_for_plot.milk.shape[0]}]"
        fig.suptitle(title, fontsize=11)
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        png = out / f"lactation_{protocol.galactagogue}_{int(protocol.galactagogue_dose_mg_day)}mg.png"
        fig.savefig(png, dpi=130)
        plt.close(fig)
        print(f"\nWrote: {png}")
    return 0


def cmd_compare(args: argparse.Namespace) -> int:
    """Head-to-head galactagogue comparison at equieffective prolactin."""
    from .comparison import compare_galactagogues, comparison_table, tradeoff_notes
    rows = compare_galactagogues(
        agents=args.agents,
        target_prolactin_ng_ml=args.target_prolactin,
        prep_weeks=args.prep_weeks,
        pumping_sessions_per_day=args.sessions,
        progesterone_start_weeks=args.p_start_weeks,
    )
    print(comparison_table(rows))
    print()
    print(tradeoff_notes(rows))
    if args.out:
        import csv
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["agent", "dose_mg_day", "prolactin_ng_ml",
                        "ss_milk_ml_day", "peak_milk_ml_day",
                        "weeks_to_50ml", "half_life_h", "us_approved"])
            for r in rows:
                w.writerow([r.name, round(r.dose_mg_day, 2),
                            round(r.prolactin_ng_ml, 1),
                            round(r.ss_milk, 1), round(r.peak_milk, 1),
                            (round(r.weeks_to_supply, 2)
                             if r.weeks_to_supply is not None else "never"),
                            r.half_life_h, r.approved_us])
        print(f"\nWrote: {out}")
    return 0


def cmd_sweep(args: argparse.Namespace) -> int:
    """Low-dose sweep with between-body heterogeneity quantified."""
    rows = dose_sweep_population(
        agent=args.agent,
        doses=args.doses,
        prep_weeks=args.prep_weeks,
        progesterone_start_weeks=args.p_start_weeks,
        pumping_sessions_per_day=args.sessions,
        n=args.n,
        seed=args.seed,
        threshold_ml_day=args.threshold,
    )
    print(dose_sweep_table(rows, threshold_ml_day=args.threshold))
    if args.out:
        import csv
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["agent", "dose_mg_day", "prolactin_ng_ml",
                        "milk_p10", "milk_p50", "milk_p90",
                        "frac_reaching_threshold"])
            for r in rows:
                w.writerow([r.agent, r.dose_mg_day,
                            round(r.prolactin_ng_ml, 1),
                            round(r.milk_p10, 1), round(r.milk_p50, 1),
                            round(r.milk_p90, 1),
                            round(r.frac_reaching_threshold, 4)])
        print(f"\nWrote: {out}")
    return 0


def cmd_risk(args: argparse.Namespace) -> int:
    """Structured risk-exposure framing for a planned galactagogue course."""
    prof = exposure_report(
        agent=args.agent,
        dose_mg_day=args.dose,
        weeks=args.weeks,
        age=args.age,
        concurrent_antipsychotics=args.antipsychotics,
    )
    print(prof.report())
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

    pg = sub.add_parser("galactagogues",
                        help="list pharmacologic galactagogues or detail one")
    pg.add_argument("--name", help="galactagogue name; omit to list all")
    pg.set_defaults(func=cmd_galactagogues)

    pl = sub.add_parser("lactation",
                        help="simulate Newman-Goldfarb-style lactation induction")
    pl.add_argument("--galactagogue", default="domperidone",
                    choices=sorted(GALACTAGOGUES.keys()))
    pl.add_argument("--dose", type=float, default=40.0,
                    help="galactagogue daily dose (mg/day)")
    pl.add_argument("--prep-weeks", type=float, default=30.0,
                    help="E2+anti-androgen preparation phase (weeks)")
    pl.add_argument("--p-start-weeks", type=float, default=None,
                    help="week progesterone starts (default: prep minus 12; "
                    "set >= prep-weeks for a P-free protocol)")
    pl.add_argument("--sessions", type=float, default=8.0,
                    help="pumping/suckling sessions per day")
    pl.add_argument("--lactation-weeks", type=float, default=52.0,
                    help="lactation phase to simulate after the drop (weeks)")
    pl.add_argument("--bands", action="store_true",
                    help="Monte-Carlo sensitivity bands on the milk parameters")
    pl.add_argument("--n", type=int, default=200,
                    help="sample count for --bands")
    pl.add_argument("--seed", type=int, default=42, help="seed for --bands")
    pl.add_argument("--out", help="output directory for the figure")
    pl.set_defaults(func=cmd_lactation)

    pcm = sub.add_parser("compare",
                         help="compare galactagogues at equieffective prolactin")
    pcm.add_argument("--agents", nargs="+",
                     default=["domperidone", "metoclopramide"],
                     choices=sorted(GALACTAGOGUES.keys()))
    pcm.add_argument("--target-prolactin", type=float, default=130.0,
                     help="target prolactin (ng/mL) all agents are dose-matched to")
    pcm.add_argument("--prep-weeks", type=float, default=30.0)
    pcm.add_argument("--p-start-weeks", type=float, default=None)
    pcm.add_argument("--sessions", type=float, default=8.0)
    pcm.add_argument("--out", help="write comparison CSV")
    pcm.set_defaults(func=cmd_compare)

    psw = sub.add_parser("sweep",
                         help="low-dose sweep; quantifies between-body heterogeneity")
    psw.add_argument("--agent", default="metoclopramide",
                     choices=sorted(GALACTAGOGUES.keys()))
    psw.add_argument("--doses", nargs="+", type=float,
                     default=[5, 10, 15, 20, 30, 45])
    psw.add_argument("--prep-weeks", type=float, default=30.0)
    psw.add_argument("--p-start-weeks", type=float, default=None)
    psw.add_argument("--sessions", type=float, default=8.0)
    psw.add_argument("--n", type=int, default=150, help="modeled bodies per dose")
    psw.add_argument("--seed", type=int, default=42)
    psw.add_argument("--threshold", type=float, default=100.0,
                     help="'meaningful supply' threshold (mL/day)")
    psw.add_argument("--out", help="write sweep CSV")
    psw.set_defaults(func=cmd_sweep)

    prk = sub.add_parser("risk",
                         help="structured risk-exposure framing for a planned course")
    prk.add_argument("--agent", default="metoclopramide",
                     choices=sorted(GALACTAGOGUES.keys()))
    prk.add_argument("--dose", type=float, default=10.0, help="mg/day")
    prk.add_argument("--weeks", type=float, default=8.0, help="course length")
    prk.add_argument("--age", type=int, default=30,
                     help="age (strongest known TD risk modifier)")
    prk.add_argument("--antipsychotics", action="store_true",
                     help="flag concurrent dopamine-antagonist use")
    prk.set_defaults(func=cmd_risk)

    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
