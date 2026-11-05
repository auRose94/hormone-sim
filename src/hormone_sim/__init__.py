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
from .compounds import (
    Compound,
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
    EstrogenicResult,
)
from .breast_development import (
    BreastModel,
    BreastResult,
    simulate_breast_development,
)

__version__ = "0.2.0"

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
    # estrogenic activity / compounds
    "Compound",
    "PHYTOESTROGENS",
    "get_compound",
    "hypothetical",
    "MW_ESTRADIOL",
    "ee_of",
    "required_daily_rate_ug_day",
    "design_summary",
    "TargetProfile",
    "EstrogenicResult",
    # breast development
    "BreastModel",
    "BreastResult",
    "simulate_breast_development",
    "__version__",
]
