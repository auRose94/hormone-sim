"""hormone-sim: population PK/PD simulation of gender-affirming hormone therapy.

A research/education toolkit. Produces population-typical concentration curves,
NOT individualized dosing. See docs/governance.md for scope and ethics.
"""

from .sim import simulate, SimResult
from .regimens import Regimen, EstradiolSpec
from .population import run_population, PopulationResult
from .params import (
    EstradiolPK,
    AndrogenPD,
    Variability,
    TARGET_E2_PGML,
    TARGET_T_NGDL,
    e2_pgml_to_pmol,
    t_ngdl_to_nmol,
)

__version__ = "0.1.0"

__all__ = [
    "simulate",
    "SimResult",
    "Regimen",
    "EstradiolSpec",
    "run_population",
    "PopulationResult",
    "EstradiolPK",
    "AndrogenPD",
    "Variability",
    "TARGET_E2_PGML",
    "TARGET_T_NGDL",
    "e2_pgml_to_pmol",
    "t_ngdl_to_nmol",
    "__version__",
]
