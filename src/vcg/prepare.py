"""Turn the raw NTP/DTT collections into four tidy analysis tables.

Outputs (all parquet, in ``data/processed/``):

* ``studies``       one row per CEBS accession: species, strain, design, route,
                    vehicle, laboratory, start year, test article;
* ``measurements``  one row per animal x continuous endpoint (terminal value);
* ``animals``       roster of animals per study/sex/dose with time on study and
                    removal reason (denominators for incidence endpoints);
* ``lesions``       one row per animal x histopathology finding;
* ``covariates``    pre-dosing body weight (``bw0``) per animal, used for matching
                    and as an ANCOVA covariate, never as an endpoint.

Only F344 rats and B6C3F1/N mice from simple dose-response 3-month and
chronic (2-year) studies are kept (see ``config``).
"""
from __future__ import annotations

import logging
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as C

log = logging.getLogger(__name__)

COMMON_COLUMNS = {
    "NIEHS CEBS Accession Number": "study",
    "Unique CEBS Subject Identifier": "animal",
    "Species": "species", "Strain": "strain", "Sex": "sex", "Study Type": "study_type",
    "Dose": "dose", "Treatment Role": "role", "Assay Name": "endpoint",
    "Assay Result": "value", "Assay Unit": "unit",
    "Exposure Duration": "duration", "Exposure Duration Unit": "duration_unit",
    "Removal Reason": "removal", "Route": "route", "Vehicle": "vehicle",
    "Study Start Year": "year", "DTXSID": "dtxsid", "Study Test Article": "test_article",
    "Dose Response Study Flag": "dose_response",
}
INLIFE_COLUMNS = {
    "Study Accession Number": "study", "Unique Subject Identifier": "animal",
    "Species": "species", "Strain": "strain", "Sex": "sex", "Study Type": "study_type",
    "Dose": "dose", "Treatment Role": "role", "Assay Name": "endpoint",
    "Standard Value": "value", "Standard Unit": "unit",
    "Time in Study": "duration", "Time in Study Unit": "duration_unit",
    "Removal Reason": "removal", "Study Start Year": "year", "DTXSID": "dtxsid",
    "Study Test Article": "test_article", "Is Study Dose Response": "dose_response",
}
HISTO_COLUMNS = {
    "NIEHS CEBS Accession Number": "study", "Unique CEBS Subject Identifier": "animal",
    "Species": "species", "Strain": "strain", "Sex": "sex", "Study Type": "study_type",
    "Dose": "dose", "Treatment Role": "role", "Tissue": "tissue", "Morphology": "morphology",
    "Neoplastic Flag": "neoplastic", "Malignant Flag": "malignant", "Severity": "severity",
    "Exposure duration": "duration", "Exposure duration unit": "duration_unit",
    "Removal Reason": "removal", "Route": "route", "Vehicle": "vehicle",
    "Study Start Year": "year", "DTXSID": "dtxsid", "Study Test Article": "test_article",
    "Dose Response Study Flag": "dose_response", "Facility": "laboratory",
}
CATALOG_COLUMNS = {
    "NIEHS CEBS Accession Number": "study", "Study Type": "study_type",
    "Study Start Year": "year", "Study Test Article": "test_article", "DTXSID": "dtxsid",
    "Species": "species", "Strain": "strain", "Route": "route", "Vehicle": "vehicle",
    "Laboratory": "laboratory", "Publication": "publication",
}


# ----------------------------------------------------------------------------
# raw readers
# ----------------------------------------------------------------------------

def read_xlsx(path: Path) -> pd.DataFrame:
    """Read a single-sheet DTT spreadsheet quickly with openpyxl in read-only mode."""
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=True)
    ws = wb.worksheets[0]
    rows = ws.iter_rows(values_only=True)
    header = [str(h) for h in next(rows)]
    df = pd.DataFrame(list(rows), columns=header)
    wb.close()
    return df


def _cached_xlsx(name: str, raw_dir: Path, cache_dir: Path) -> pd.DataFrame:
    cache = cache_dir / f"{name}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    t0 = time.time()
    df = read_xlsx(raw_dir / C.RAW_FILES[name])
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype("string")
    cache_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache, index=False)
    log.info("read %s: %s rows in %.0f s", name, len(df), time.time() - t0)
    return df


def _duration_days(duration: pd.Series, unit: pd.Series) -> pd.Series:
    d = pd.to_numeric(duration, errors="coerce")
    u = unit.astype("string").str.strip().str.lower()
    days = pd.Series(np.nan, index=d.index, dtype=float)
    days[u == "day"] = d[u == "day"]
    days[u == "week"] = d[u == "week"] * 7
    days[u == "month"] = d[u == "month"] * 30.4
    return days


def _strain_matches(df: pd.DataFrame) -> np.ndarray:
    """True where the row's strain is the one selected for its species (NA-safe)."""
    expected = df["species"].astype("string").map(C.STRAINS)
    return (expected == df["strain"].astype("string")).fillna(False).to_numpy(dtype=bool)


def _design(study_type: pd.Series) -> pd.Series:
    st = study_type.astype("string").replace(C.STUDY_TYPE_ALIASES)
    inv = {v: k for k, v in C.STUDY_TYPES.items()}
    return st.map(inv)


def _select(df: pd.DataFrame) -> pd.DataFrame:
    """Keep the species/strain/design of interest, simple dose-response designs,
    control and treated animals with a numeric dose."""
    df = df.copy()
    df["design"] = _design(df["study_type"])
    keep = df["design"].notna()
    keep &= _strain_matches(df)
    keep &= df["dose_response"].astype("string").eq("T").fillna(False).to_numpy(dtype=bool)
    keep &= df["role"].isin((*C.CONTROL_ROLES, C.TREATED_ROLE))
    df = df.loc[keep].copy()
    # identifiers are integers in the spreadsheets and strings in the CSVs
    df["animal"] = df["animal"].astype("string")
    df["study"] = df["study"].astype("string")
    df["dose"] = pd.to_numeric(df["dose"], errors="coerce")
    df = df[df["dose"].notna()]
    df["is_control"] = df["role"].isin(C.CONTROL_ROLES) & (df["dose"] == 0)
    df = df[df["is_control"] | (df["role"] == C.TREATED_ROLE)]
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    return df


def _single_control_group(df: pd.DataFrame) -> pd.DataFrame:
    """If a study/sex has both a vehicle and an untreated control group, keep the vehicle one."""
    ctrl = df[df["is_control"]]
    roles = ctrl.groupby(["study", "sex"])["role"].agg(lambda s: set(s))
    both = roles[roles.apply(lambda s: len(s) > 1)].index
    if len(both) == 0:
        return df
    drop = df.set_index(["study", "sex"]).index.isin(both) & (df["role"] == "Naive control").to_numpy()
    return df.loc[~drop]


# ----------------------------------------------------------------------------
# continuous endpoints
# ----------------------------------------------------------------------------

def _terminal_measurements(df: pd.DataFrame, domain: str, endpoints: list[str] | None) -> pd.DataFrame:
    df = _select(df.rename(columns=COMMON_COLUMNS))
    if endpoints is not None:
        df = df[df["endpoint"].isin(endpoints)]
    df["days"] = _duration_days(df["duration"], df["duration_unit"])
    lo_hi = df["design"].map(C.TERMINAL_WINDOW_DAYS)
    df = df[df["days"].between(lo_hi.str[0], lo_hi.str[1])]
    if domain == "terminal_bw":
        df = df[df["removal"].isin(C.TERMINAL_REMOVALS)]
    else:
        df = df[~df["removal"].eq(C.INTERIM_REMOVAL)]
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df[df["value"].notna()]
    # one unit per endpoint: keep the modal one
    modal = df.groupby("endpoint")["unit"].agg(lambda s: s.mode().iat[0] if len(s.mode()) else pd.NA)
    df = df[df["unit"].eq(df["endpoint"].map(modal))]
    # one value per animal x endpoint: the latest measurement in the window
    df = df.sort_values("days").groupby(["study", "sex", "endpoint", "animal"], as_index=False).last()
    df["domain"] = domain
    return df


def _read_inlife(path: Path, cache_dir: Path) -> pd.DataFrame:
    """Chunked pass over the multi-GB in-life CSV keeping body mass rows of the
    species/strains of interest; cached because it takes minutes."""
    cache = cache_dir / "inlife_bw.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    parts = []
    t0 = time.time()
    for i, chunk in enumerate(pd.read_csv(path, usecols=list(INLIFE_COLUMNS), chunksize=1_000_000, dtype=str,
                                          encoding="utf-8", encoding_errors="replace",
                                          on_bad_lines="skip")):
        chunk = chunk.rename(columns=INLIFE_COLUMNS)
        chunk = chunk[_strain_matches(chunk) & chunk["endpoint"].eq("Body mass assay").to_numpy()]
        if len(chunk):
            parts.append(chunk)
        if i % 5 == 0:
            log.info("  in-life chunk %d (%.0f s), kept %d rows so far", i, time.time() - t0,
                     sum(len(p) for p in parts))
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=list(INLIFE_COLUMNS.values()))
    for col in df.columns:
        df[col] = df[col].astype("string")
    cache_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache, index=False)
    return df


def _inlife_measurements(raw: pd.DataFrame, studies: pd.DataFrame) -> pd.DataFrame:
    """Body mass at nominal time points plus body weight gain."""
    df = raw[raw["study"].isin(set(studies.index))]
    if df.empty:
        return pd.DataFrame()
    df = _select(df)
    df["days"] = _duration_days(df["duration"], df["duration_unit"])
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df[df["value"].notna() & df["days"].notna() & df["unit"].str.lower().eq("g")]
    out = []
    for design, points in C.INLIFE_TIMEPOINTS.items():
        d = df[df["design"] == design]
        for label, (nominal, lo, hi) in points.items():
            sel = d[d["days"].between(lo, hi)].copy()
            sel["dist"] = (sel["days"] - nominal).abs()
            sel = sel.sort_values("dist").groupby(["study", "sex", "animal"], as_index=False).first()
            sel["endpoint"] = f"Body weight, {label}"
            out.append(sel)
    df = pd.concat(out, ignore_index=True)
    # body weight gain from the initial weighing to the last nominal point of the design
    gains = []
    for design, points in C.INLIFE_TIMEPOINTS.items():
        last = list(points)[-1]
        d = df[df["design"] == design]
        init = d[d["endpoint"] == "Body weight, initial"].set_index(["study", "sex", "animal"])["value"]
        fin = d[d["endpoint"] == f"Body weight, {last}"].set_index(["study", "sex", "animal"])
        common = fin.index.intersection(init.index)
        g = fin.loc[common].copy()
        g["value"] = g["value"] - init.loc[common].to_numpy()
        g["endpoint"] = f"Body weight gain, initial to {last}"
        gains.append(g.reset_index())
    df = pd.concat([df, *gains], ignore_index=True)
    df["domain"] = "inlife_bw"
    return df


# ----------------------------------------------------------------------------
# study metadata
# ----------------------------------------------------------------------------

def _study_table(catalog: pd.DataFrame, frames: list[pd.DataFrame]) -> pd.DataFrame:
    """One row per accession.  Catalog values first, table values as fall-back."""
    cat = catalog.rename(columns=CATALOG_COLUMNS)[list(CATALOG_COLUMNS.values())].copy()
    cat["year"] = pd.to_numeric(cat["year"], errors="coerce")
    cat = cat.drop_duplicates("study").set_index("study")
    cols = ["species", "strain", "study_type", "route", "vehicle", "year", "dtxsid", "test_article"]
    tab = pd.concat([f[["study", *[c for c in cols if c in f.columns],
                        *(["laboratory"] if "laboratory" in f.columns else [])]] for f in frames],
                    ignore_index=True)
    tab = tab.groupby("study").agg(lambda s: s.dropna().mode().iat[0] if s.dropna().size else pd.NA)
    studies = tab.reindex(tab.index.union(cat.index))
    for col in [*cols, "laboratory", "publication"]:
        if col in cat.columns:
            studies[col] = cat[col].reindex(studies.index).fillna(studies.get(col, pd.Series(dtype=object)))
    studies["laboratory"] = studies["laboratory"].map(C.normalise_lab)
    studies["vehicle"] = studies["vehicle"].fillna(C.MISSING_VEHICLE)
    studies["design"] = _design(studies["study_type"].astype("string").replace(C.STUDY_TYPE_ALIASES))
    studies["year"] = pd.to_numeric(studies["year"], errors="coerce")
    studies = studies[studies["design"].notna().to_numpy() & _strain_matches(studies)]
    studies.index.name = "study"
    return studies


# ----------------------------------------------------------------------------
# incidence: roster and lesions
# ----------------------------------------------------------------------------

def _read_histopath(path: Path, cache_dir: Path) -> pd.DataFrame:
    cache = cache_dir / "histopath.parquet"
    if cache.exists():
        return _select(pd.read_parquet(cache))
    parts = []
    t0 = time.time()
    for i, chunk in enumerate(pd.read_csv(path, usecols=list(HISTO_COLUMNS), chunksize=1_000_000,
                                          dtype=str, encoding="utf-8", encoding_errors="replace",
                                          on_bad_lines="skip")):
        chunk = chunk.rename(columns=HISTO_COLUMNS)
        chunk = chunk[_strain_matches(chunk) & chunk["study_type"].isin(C.STUDY_TYPES.values()).to_numpy()]
        parts.append(chunk)
        log.info("  histopathology chunk %d (%.0f s)", i, time.time() - t0)
    df = pd.concat(parts, ignore_index=True)
    for col in df.columns:
        df[col] = df[col].astype("string")
    cache_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache, index=False)
    return _select(df)


def _roster(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Animals per study/sex/dose with time on study (days) and removal reason."""
    cols = ["study", "sex", "dose", "role", "is_control", "animal", "removal", "days", "design"]
    pieces = []
    for prio, f in enumerate(frames):
        p = f[cols].copy()
        p["prio"] = prio
        pieces.append(p)
    r = pd.concat(pieces, ignore_index=True)
    r = r.sort_values(["prio", "days"], ascending=[True, False])
    r = r.groupby(["study", "sex", "animal"], as_index=False).first()
    r = r[~r["removal"].eq(C.INTERIM_REMOVAL)]
    r["terminal"] = r["removal"].isin(C.TERMINAL_REMOVALS)
    # terminal sacrifice time per study/sex: the typical time of scheduled sacrifice
    tmax = (r[r["terminal"]].groupby(["study", "sex"])["days"].median()
            .rename("t_max").reset_index())
    r = r.merge(tmax, on=["study", "sex"], how="left")
    default = r["design"].map({"subchronic": 91.0, "chronic": 728.0})
    r["t_max"] = r["t_max"].fillna(default)
    # animals without a recorded time: terminal animals get the sacrifice time,
    # others (rare) the same, which gives them full poly-3 weight
    r["days"] = r["days"].fillna(r["t_max"])
    return r.drop(columns="prio")


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------

def _load_one(args: tuple[str, Path, Path]) -> tuple[str, pd.DataFrame]:
    name, raw_dir, cache_dir = args
    return name, _cached_xlsx(name, raw_dir, cache_dir)


def build(raw_dir: Path = C.RAW_DIR, out_dir: Path = C.PROCESSED_DIR, workers: int = 4) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / "cache"
    names = ["clin_chem", "hematology", "organ_weight", "terminal_bw"]
    missing = [n for n in names if not (raw_dir / C.RAW_FILES[n]).exists() and not (cache_dir / f"{n}.parquet").exists()]
    if missing:
        raise FileNotFoundError(f"raw files missing: {missing}; run data/download_ntp.py first")

    with ProcessPoolExecutor(max_workers=min(workers, len(names))) as pool:
        raw = dict(pool.map(_load_one, [(n, raw_dir, cache_dir) for n in names]))

    panels = {"clin_chem": C.CLIN_CHEM, "hematology": C.HEMATOLOGY,
              "organ_weight": C.ORGAN_WEIGHTS, "terminal_bw": None}
    frames = {n: _terminal_measurements(raw[n], n, panels[n]) for n in names}
    for n, f in frames.items():
        log.info("%s: %d terminal values, %d studies", n, len(f), f["study"].nunique())

    catalog = pd.read_csv(raw_dir / C.RAW_FILES["catalog"], encoding="utf-8-sig")

    histo_path = raw_dir / C.RAW_FILES["histopath"]
    histo = None
    if histo_path.exists():
        histo = _read_histopath(histo_path, cache_dir)
        histo["days"] = _duration_days(histo["duration"], histo["duration_unit"])
        log.info("histopathology: %d findings, %d studies", len(histo), histo["study"].nunique())
    else:
        log.warning("histopathology CSV not found; incidence endpoints will be skipped")

    studies = _study_table(catalog, [*frames.values(), *([histo] if histo is not None else [])])

    inlife_path = raw_dir / C.RAW_FILES["inlife_bw"]
    if inlife_path.exists():
        frames["inlife_bw"] = _inlife_measurements(_read_inlife(inlife_path, cache_dir), studies)
        log.info("in-life body weight: %d values, %d studies", len(frames["inlife_bw"]),
                 frames["inlife_bw"]["study"].nunique())
    else:
        log.warning("in-life body weight CSV not found; body weight trajectory endpoints skipped")

    keep = ["study", "sex", "dose", "role", "is_control", "animal", "domain", "endpoint", "value", "unit", "days"]
    meas = pd.concat([f[keep] for f in frames.values() if len(f)], ignore_index=True)
    meas = _single_control_group(meas)
    meas = meas[meas["study"].isin(studies.index)]
    meas["endpoint"] = meas["endpoint"].replace(C.SHORT_NAMES)
    # The pre-dosing body weight is a matching covariate, not an endpoint: it is
    # measured before any treatment, so a "treatment effect" on it can only be chance
    # or a mismatch between the control groups.
    is_bw0 = meas["endpoint"].eq(C.INITIAL_BW_ENDPOINT)
    covariates = (meas.loc[is_bw0, ["study", "sex", "animal", "value"]]
                  .rename(columns={"value": "bw0"})
                  .drop_duplicates(["study", "sex", "animal"]))
    meas = meas.loc[~is_bw0]

    roster_frames = [f for f in ([histo] if histo is not None else []) + [frames["terminal_bw"]]
                     if f is not None and len(f)]
    for f in roster_frames:
        f["design"] = f["design"] if "design" in f else _design(f["study_type"])
    roster = _roster(roster_frames)
    roster = _single_control_group(roster)
    roster = roster[roster["study"].isin(studies.index)]

    if histo is not None:
        lesions = histo[["study", "sex", "animal", "tissue", "morphology", "neoplastic", "malignant", "severity"]]
        lesions = lesions.drop_duplicates(["study", "sex", "animal", "tissue", "morphology"])
        lesions = lesions.merge(roster[["study", "sex", "animal"]], on=["study", "sex", "animal"])
    else:
        lesions = pd.DataFrame(columns=["study", "sex", "animal", "tissue", "morphology", "neoplastic", "malignant", "severity"])

    studies.reset_index().to_parquet(out_dir / "studies.parquet", index=False)
    meas.to_parquet(out_dir / "measurements.parquet", index=False)
    roster.to_parquet(out_dir / "animals.parquet", index=False)
    lesions.to_parquet(out_dir / "lesions.parquet", index=False)
    covariates.to_parquet(out_dir / "covariates.parquet", index=False)
    log.info("written: %d studies, %d measurements, %d animals, %d lesion records, %d pre-dosing weights",
             len(studies), len(meas), len(roster), len(lesions), len(covariates))


TABLES = ("studies", "measurements", "animals", "lesions", "covariates")


def load_processed(dir_: Path = C.PROCESSED_DIR) -> dict[str, pd.DataFrame]:
    """Load the processed tables; ``covariates`` holds the pre-dosing body weight (bw0)."""
    return {name: pd.read_parquet(dir_ / f"{name}.parquet") for name in TABLES}
