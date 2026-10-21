"""Project-wide constants: paths, study selection, endpoint panels, matching levels."""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# The raw collections are ~9 GB; point VCG_DATA_DIR somewhere with space if the
# repository itself lives on a synced or small volume.
DATA_DIR = Path(os.environ.get("VCG_DATA_DIR", ROOT / "data"))
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"

RAW_FILES = {
    "clin_chem": "202602_Clinical_Chemistry_IAD.xlsx",
    "hematology": "202602_Hematology_IAD.xlsx",
    "organ_weight": "202602_Organ_Weight_IAD.xlsx",
    "terminal_bw": "202602_Terminal_Bodyweight_IAD.xlsx",
    "inlife_bw": "202507_In_Life_Bodyweight_IAD.csv",
    "histopath": "202602_Histopathology_IAD.csv",
    "catalog": "cebs_study_catalog.csv",
}

# --- study selection --------------------------------------------------------
STRAINS = {"Rat": "F344", "Mouse": "B6C3F1/N"}          # species -> strain
STUDY_TYPES = {"subchronic": "3-month toxicology design", "chronic": "Chronic"}
# the in-life collection uses its own study type labels
STUDY_TYPE_ALIASES = {
    "Chronic toxicology study design": "Chronic",
    "3-month toxicology design": "3-month toxicology design",
}
CONTROL_ROLES = ("Vehicle control", "Naive control", "Untreated control")  # in-life file says "Untreated"
TREATED_ROLE = "Treated"
# animals considered for terminal measurements / incidence denominators
TERMINAL_REMOVALS = ("Scheduled Sacrifice", "Scheduled removal", "Scheduled removal, terminal")
MISSING_VEHICLE = "(none recorded)"   # feed/water studies often record no vehicle; treated as its own category
INTERIM_REMOVAL = "Scheduled removal, interim"

# terminal windows in days on study, per study design
TERMINAL_WINDOW_DAYS = {"subchronic": (80, 105), "chronic": (690, 770)}
# in-life body weight: nominal time points (days) and tolerance
INLIFE_TIMEPOINTS = {"subchronic": {"initial": (0, 0, 3), "week 13": (91, 84, 98)},
                     "chronic": {"initial": (0, 0, 3), "week 52": (364, 350, 380), "week 104": (728, 700, 745)}}

# --- endpoint panels (standard DTT assay names) -----------------------------
ORGAN_WEIGHTS = ["Liver weight", "Kidney weight, right", "Heart weight", "Lung weight",
                 "Thymus weight", "Testis weight, right"]
CLIN_CHEM = ["Alanine aminotransferase level of activity assay",
             "Sorbitol dehydrogenase level of activity assay",
             "Alkaline phosphatase level of activity assay",
             "Urea nitrogen assay", "Bile acid concentration assay",
             "Albumin concentration assay", "Creatinine concentration assay",
             "Protein concentration assay", "Creatine kinase level of activity assay"]
HEMATOLOGY = ["Hemoglobin concentration assay", "Erythrocyte count assay", "Hematocrit assay",
              "Mean cell hemoglobin concentration assay", "Mean cell volume assay",
              "Mean cell hemoglobin assay", "Platelet count assay", "Leukocyte count assay",
              "Lymphocyte count assay", "Neutrophil count assay", "Reticulocyte count assay",
              "Eosinophil count assay", "Monocyte count assay"]
SHORT_NAMES = {
    "Alanine aminotransferase level of activity assay": "ALT",
    "Sorbitol dehydrogenase level of activity assay": "SDH",
    "Alkaline phosphatase level of activity assay": "ALP",
    "Urea nitrogen assay": "Urea nitrogen",
    "Bile acid concentration assay": "Bile acids",
    "Albumin concentration assay": "Albumin",
    "Creatinine concentration assay": "Creatinine",
    "Protein concentration assay": "Total protein",
    "Creatine kinase level of activity assay": "Creatine kinase",
    "Hemoglobin concentration assay": "Hemoglobin",
    "Erythrocyte count assay": "Erythrocytes",
    "Hematocrit assay": "Hematocrit",
    "Mean cell hemoglobin concentration assay": "MCHC",
    "Mean cell volume assay": "MCV",
    "Mean cell hemoglobin assay": "MCH",
    "Platelet count assay": "Platelets",
    "Leukocyte count assay": "Leukocytes",
    "Lymphocyte count assay": "Lymphocytes",
    "Neutrophil count assay": "Neutrophils",
    "Reticulocyte count assay": "Reticulocytes",
    "Eosinophil count assay": "Eosinophils",
    "Monocyte count assay": "Monocytes",
    "Terminal body weight assay": "Terminal body weight",
}
# minimum number of affected animals (all groups pooled) for a lesion to be tested
MIN_AFFECTED_FOR_INCIDENCE = 3

# --- matching ---------------------------------------------------------------
# Nested strictness levels: each adds one criterion to the previous one.  The
# base level is always species + strain + sex + study design.
MATCH_LEVELS: dict[str, tuple[str, ...]] = {
    "L0 design": (),
    "L1 +route": ("route",),
    "L2 +vehicle": ("route", "vehicle"),
    "L3 +laboratory": ("route", "vehicle", "laboratory"),
    "L4 +period (3y)": ("route", "vehicle", "laboratory", "period"),
}
PERIOD_YEARS = 3                      # |start year difference| allowed at the strictest level
DRIFT_LEVEL = "L2 +vehicle"           # level used for the drift (time gap) analysis
DRIFT_BANDS = [(0, 2), (3, 5), (6, 10), (11, 20), (21, 60)]
N_DRAWS = 20                          # virtual control groups drawn per configuration
ALPHA = 0.05
SEED = 20261021

# --- laboratory name normalisation -----------------------------------------
# The catalog and the histopathology table spell the same laboratory several
# ways (legacy upper-case codes versus full names).
LAB_ALIASES = {
    "MASON": "TSI Mason Research", "TSI MASON RESEARCH": "TSI Mason Research",
    "BATTELLE COLUMBUS": "Battelle Columbus", "BATTELLE COLUMBUS LABORATORY": "Battelle Columbus",
    "BATTELLE NORTHWEST": "Battelle Northwest",
    "SOUTHERN": "Southern Research", "SOUTHERN RESEARCH INSTITUTE": "Southern Research",
    "SOUTHERN RESEARCH": "Southern Research",
    "HAZLETON": "Hazleton", "HAZLETON LABORATORIES AMERICA, INC.": "Hazleton",
    "HAZLETON LABS, ROCKVILLE, GOV'T SERVICES, INC": "Hazleton",
    "IITRI": "IIT Research Institute", "I. I. T. RESEARCH INSTITUTE": "IIT Research Institute",
    "LITTON BIONETICS": "Litton Bionetics",
    "NATIONAL CENTER FOR TOXICOLOGICAL RESEARCH": "NCTR",
    "SRI": "SRI International", "SRI INTERNATIONAL": "SRI International",
    "MICROBIOLOGICAL": "Microbiological Associates",
    "MICROBIOLOGICAL ASSOCIATES": "Microbiological Associates",
    "PAPANICOLAOU": "Papanicolaou Cancer Research Institute",
    "PHYSIOLOGICAL": "Physiological Research Laboratories",
    "GULF SOUTH": "Gulf South Research Institute",
    "FCRF": "Frederick Cancer Research Facility",
    "SPRINGBORN": "Springborn Laboratories",
    "BIORELIANCE INC.": "BioReliance",
    "INTEGRATED LABORATORY SYSTEMS": "Integrated Laboratory Systems",
    "LOVELACE INHALATION TOXICOLOGY RESEARCH INSTITUTE": "Lovelace",
    "MIDWEST RESEARCH INSTITUTE (INHALATION)": "Midwest Research Institute",
    "MIDWEST RESEARCH INSTITUTE": "Midwest Research Institute",
    "EG&G MASON RESEARCH INSTITUTE": "TSI Mason Research",
    "NATL INST. FOR ENVIRONMENTAL HEALTH SCIENCES": "NIEHS",
}


def normalise_lab(name: object) -> str | None:
    if name is None or (isinstance(name, float) and name != name):
        return None
    key = str(name).strip().upper()
    if not key:
        return None
    return LAB_ALIASES.get(key, str(name).strip())
