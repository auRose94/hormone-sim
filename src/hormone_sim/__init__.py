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
from .galactagogues import (
    GALACTAGOGUES,
    Galactagogue,
    get_galactagogue,
    simulate_dosing_pk,
    prolactin_response,
    prolactin_curve,
    list_galactagogues,
    PROLACTIN_BASELINE_NG_ML,
)
from .lactation import (
    LactationModel,
    LactationProtocol,
    LactationResult,
    simulate_lactation,
)
from .comparison import (
    AgentComparison,
    compare_galactagogues,
    equieffective_dose,
    comparison_table,
    tradeoff_notes,
)
from .sensitivity import (
    SensitivitySpec,
    LactationPopulationResult,
    simulate_lactation_population,
)

__version__ = "0.4.0"

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
    # galactagogues / lactation
    "GALACTAGOGUES",
    "Galactagogue",
    "get_galactagogue",
    "simulate_dosing_pk",
    "prolactin_response",
    "prolactin_curve",
    "list_galactagogues",
    "PROLACTIN_BASELINE_NG_ML",
    "LactationModel",
    "LactationProtocol",
    "LactationResult",
    "simulate_lactation",
    # comparison / sensitivity
    "AgentComparison",
    "compare_galactagogues",
    "equieffective_dose",
    "comparison_table",
    "tradeoff_notes",
    "SensitivitySpec",
    "LactationPopulationResult",
    "simulate_lactation_population",
    "__version__",
]
