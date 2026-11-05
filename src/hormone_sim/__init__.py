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
    DoseSweepRow,
    dose_sweep_population,
    dose_sweep_table,
)
from .risk_exposure import exposure_report, ExposureProfile
from .progestogens import (
    PROGESTOGENS,
    Progestogen,
    get_progestogen,
    custom_progestogen,
)
from .personalization import (
    Observation,
    PersonalizedResult,
    personalize,
    posterior_predictive,
)
from .ingest import (
    load_observations,
    personalize_from_csv,
    plot_personal,
    IngestError,
)
from .aims import (
    AIMS_ITEMS,
    AIMSVisit,
    AIMSRecord,
    MONITORING_SCHEDULE,
    parse_visit,
)
from .protocol import ProtocolConfig, generate_protocol
from .androgens import (
    ANDROGENS,
    Androgen,
    get_androgen,
    simulate_androgen,
    ugl_to_ngdl,
    ngdl_to_ugl,
    BASELINE_T_NGDL,
    BASELINE_E2_PGML,
)
from .masculinization import (
    MasculinizationModel,
    MascResult,
    simulate_masculinization,
    MascPopulationResult,
    simulate_masc_population,
)
from .erythrocytosis import (
    HctModel,
    HctResult,
    simulate_hct,
    hct_risk_report,
)
from .intersex import (
    PACKS as INTERSEX_PACKS,
    ConditionPack,
    get_pack,
    pack_summary,
)
from .personalization_androgen import (
    AndrogenObservation,
    AndrogenPersonalizedResult,
    personalize_androgen,
    posterior_predictive_androgen,
    load_androgen_observations,
)

__version__ = "0.9.0"

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
    "DoseSweepRow",
    "dose_sweep_population",
    "dose_sweep_table",
    "exposure_report",
    "ExposureProfile",
    "PROGESTOGENS",
    "Progestogen",
    "get_progestogen",
    "custom_progestogen",
    "Observation",
    "PersonalizedResult",
    "personalize",
    "posterior_predictive",
    "load_observations",
    "personalize_from_csv",
    "plot_personal",
    "IngestError",
    "AIMS_ITEMS",
    "AIMSVisit",
    "AIMSRecord",
    "MONITORING_SCHEDULE",
    "parse_visit",
    "ProtocolConfig",
    "generate_protocol",
    "ANDROGENS",
    "Androgen",
    "get_androgen",
    "simulate_androgen",
    "ugl_to_ngdl",
    "ngdl_to_ugl",
    "BASELINE_T_NGDL",
    "BASELINE_E2_PGML",
    "MasculinizationModel",
    "MascResult",
    "simulate_masculinization",
    "MascPopulationResult",
    "simulate_masc_population",
    "HctModel",
    "HctResult",
    "simulate_hct",
    "hct_risk_report",
    "INTERSEX_PACKS",
    "ConditionPack",
    "get_pack",
    "pack_summary",
    "AndrogenObservation",
    "AndrogenPersonalizedResult",
    "personalize_androgen",
    "posterior_predictive_androgen",
    "load_androgen_observations",
    "__version__",
]
