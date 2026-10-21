# Data

No data is committed to this repository. `download_ntp.py` fetches everything
from the U.S. National Toxicology Program (NTP) / Division of Translational
Toxicology (DTT) into `data/raw/` (or `$VCG_DATA_DIR/raw/`):

| File | Content | Size |
|---|---|---|
| `202602_Clinical_Chemistry_IAD.xlsx` | individual-animal clinical chemistry | 74 MB |
| `202602_Hematology_IAD.xlsx` | individual-animal hematology | 150 MB |
| `202602_Organ_Weight_IAD.xlsx` | individual-animal organ weights | 66 MB |
| `202602_Terminal_Bodyweight_IAD.xlsx` | individual-animal terminal body weights | 77 MB |
| `202507_In_Life_Bodyweight_IAD.csv` | individual-animal in-life body weights (all weighings) | 5.8 GB |
| `202602_Histopathology_IAD.csv` | individual-animal histopathology findings | 2.5 GB |
| `DTT_Integrated_Data_Page_Definitions.xlsx` | column definitions for the above | 39 KB |
| `cebs_study_catalog.csv` | one row per CEBS study: laboratory, start year, route, vehicle | 2.5 MB |

The first seven are the versioned "DTT Integrated High Value Datasets"
published through CEBS (https://cebs.niehs.nih.gov/cebs/paper/16015,
DOI 10.22427/NTP-DATA-110-021-001-000-6). Their SHA-256 checksums are pinned in
the script. The study catalog is exported through the same endpoint the
"Download" button of the CEBS Guided Search uses
(https://cebs-ext.niehs.nih.gov/datasets/search/cebs-study-catalog); it is not
versioned, so the export date is recorded in `download_manifest.json`.

`--skip-large` skips the two CSV files. The pipeline then runs without the
in-life body weight and histopathology endpoints and says so in the log.

The NTP Historical Controls database (https://ntp.niehs.nih.gov/data/controls)
was reviewed but is **not** used: it publishes group-level summary tables as PDF,
whereas this project needs individual animals. The CEBS individual-animal
collections above contain the same control animals.

Processed tables are written to `data/processed/` (parquet) and are likewise
not committed.
