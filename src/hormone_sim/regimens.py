"""Regimen definitions and YAML loading.

A *regimen* is the editable input to a simulation: route, estradiol dose,
injectable interval, and any anti-androgen agents. Regimens are YAML files so
the input is the reproducible artifact (see ``configs/``).

Example (``configs/tf_transdermal.yaml``)::

    name: transdermal 50 ug/day
    population: trans_feminine
    estradiol:
      route: transdermal
      daily_rate_ug_day: 50
    anti_androgen:
      spironolactone: 1.0      # relative dose units (~100-200 mg/d)
    notes: >
      Typical first-line transdermal regimen; avoids hepatic first-pass.

This is a *population* modeling input, not a prescription.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

import yaml

from .params import ROUTES, RouteParams


@dataclass
class EstradiolSpec:
    route: str
    daily_rate_ug_day: float = 0.0
    injection_ug: Optional[float] = None
    injection_interval_hours: Optional[float] = None

    def route_params(self) -> RouteParams:
        if self.route not in ROUTES:
            raise KeyError(
                f"Unknown route {self.route!r}. Valid: {sorted(ROUTES)}"
            )
        return ROUTES[self.route]


@dataclass
class Regimen:
    name: str
    population: str
    estradiol: EstradiolSpec
    anti_androgen: Dict[str, float] = field(default_factory=dict)
    notes: str = ""

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Regimen":
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(p)
        with p.open() as fh:
            data = yaml.safe_load(fh)
        if not data:
            raise ValueError(f"Empty regimen file: {p}")
        e = data.get("estradiol", {})
        if not e.get("route"):
            raise ValueError("regimen.estradiol.route is required")
        return cls(
            name=data.get("name", p.stem),
            population=data.get("population", "unspecified"),
            estradiol=EstradiolSpec(
                route=e["route"],
                daily_rate_ug_day=float(e.get("daily_rate_ug_day", 0.0)),
                injection_ug=(
                    float(e["injection_ug"]) if e.get("injection_ug") else None
                ),
                injection_interval_hours=(
                    float(e["injection_interval_hours"])
                    if e.get("injection_interval_hours")
                    else None
                ),
            ),
            anti_androgen={
                str(k): float(v)
                for k, v in (data.get("anti_androgen") or {}).items()
            },
            notes=data.get("notes", ""),
        )

    def summary(self) -> str:
        lines = [
            f"Regimen: {self.name}  (population: {self.population})",
            f"  Estradiol: {self.estradiol.daily_rate_ug_day:g} ug/day, "
            f"route={self.estradiol.route}",
        ]
        if self.anti_androgen:
            aa = ", ".join(f"{k}={v:g}" for k, v in self.anti_androgen.items())
            lines.append(f"  Anti-androgen: {aa}")
        if self.notes:
            lines.append(f"  Notes: {self.notes}")
        return "\n".join(lines)
