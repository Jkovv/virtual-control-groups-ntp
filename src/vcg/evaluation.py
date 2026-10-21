"""Compare treated-versus-control conclusions under concurrent and virtual controls.

For every study/sex ("unit"), endpoint and dose group we record whether the
comparison with the control is significant (alpha = 0.05) when the control is

* the concurrent control group (the reference conclusion),
* a virtual control group drawn from matched historical controls, at each
  matching level and for each of ``N_DRAWS`` random draws,
* a partial VCG (half concurrent, half virtual),
* a VCG restricted to historical studies in a given time-gap band (drift).

The raw per-draw records are aggregated by ``summarise`` into agreement, false
positive and false negative rates.  "False positive" means the VCG comparison is
significant while the concurrent one is not; "false negative" the reverse.  Both
are relative to the concurrent control, never to the truth.
"""
from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import config as C
from . import matching as M
from . import stats as S

log = logging.getLogger(__name__)

RECORD_COLUMNS = ["study", "sex", "species", "design", "domain", "endpoint", "dose_rank", "dose",
                  "mode", "level", "band", "draw", "n_control", "pool_animals", "pool_studies",
                  "p_conc", "p_vcg", "sig_conc", "sig_vcg", "dir_conc", "dir_vcg", "p_ctrl_vs_vcg"]


@dataclass
class Unit:
    """Everything the worker needs for one study/sex: its own data and the pools."""
    target: M.Target
    continuous: pd.DataFrame        # measurements of this unit (all endpoints)
    roster: pd.DataFrame            # animals of this unit with time on study
    lesions: pd.DataFrame           # findings of this unit
    hist_controls: pd.DataFrame     # candidate historical control measurements (all endpoints)
    hist_roster: pd.DataFrame       # candidate historical control animals
    hist_lesions: pd.DataFrame      # findings of those animals
    n_draws: int = C.N_DRAWS
    seed: int = C.SEED


# ----------------------------------------------------------------------------
# continuous endpoints
# ----------------------------------------------------------------------------

def _dunnett_flags(groups: list[np.ndarray], control: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    res = S.dunnett(groups, control)
    return res.pvalue, res.pvalue < C.ALPHA, np.sign(res.statistic)


def _welch_p(a: np.ndarray, b: np.ndarray) -> float:
    from scipy import stats
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return float(stats.ttest_ind(a, b, equal_var=False).pvalue)


def evaluate_continuous(unit: Unit) -> list[list]:
    rng = np.random.default_rng([unit.seed, M.stable_key(unit.target.study, unit.target.sex)])
    out: list[list] = []
    tgt = unit.target
    if unit.hist_controls.empty:
        return out
    level_masks = {lvl: M.pool_mask(unit.hist_controls, tgt, crit) for lvl, crit in C.MATCH_LEVELS.items()}
    drift_base = level_masks[C.DRIFT_LEVEL]
    band_masks = {f"{lo}-{hi}": drift_base & M.year_gap_mask(unit.hist_controls, tgt, (lo, hi))
                  for lo, hi in C.DRIFT_BANDS}
    configs = [("vcg", lvl, "", m) for lvl, m in level_masks.items()]
    configs += [("partial", lvl, "", m) for lvl, m in level_masks.items()]
    configs += [("vcg", C.DRIFT_LEVEL, band, m) for band, m in band_masks.items()]
    hist_endpoint = unit.hist_controls["endpoint"].to_numpy()
    hist_values = unit.hist_controls["value"].to_numpy(dtype=float)
    hist_study = unit.hist_controls["study"].to_numpy()

    for (domain, endpoint), d in unit.continuous.groupby(["domain", "endpoint"], sort=False):
        ctrl = d.loc[d["is_control"], "value"].to_numpy(dtype=float)
        treated = d[~d["is_control"]]
        doses = [x for x in sorted(treated["dose"].unique()) if (treated["dose"] == x).sum() >= 2]
        groups = [treated.loc[treated["dose"] == x, "value"].to_numpy(dtype=float) for x in doses]
        n0 = len(ctrl)
        if n0 < 3 or not groups:
            continue
        ref = S.dunnett(groups, ctrl)
        sig_conc, dir_conc = ref.pvalue < C.ALPHA, np.sign(ref.statistic)
        base = [tgt.study, tgt.sex, tgt.species, tgt.design, domain, endpoint]
        in_endpoint = hist_endpoint == endpoint
        R = unit.n_draws

        for mode, level, band, mask in configs:
            sel = mask & in_endpoint
            pool = hist_values[sel]
            n_virtual = n0 if mode == "vcg" else n0 - n0 // 2
            # require two groups' worth of animals for a full VCG so that draws differ
            if len(pool) < (2 * n0 if mode == "vcg" else n_virtual):
                continue
            pool_s = len(np.unique(hist_study[sel]))
            draws = np.stack([rng.choice(pool, size=n_virtual, replace=False) for _ in range(R)])
            if mode == "partial":
                kept = np.stack([rng.choice(ctrl, size=n0 - n_virtual, replace=False) for _ in range(R)])
                draws = np.concatenate([kept, draws], axis=1)
            res = S.dunnett_many(groups, draws)
            p_cv = S.welch_p_many(ctrl, draws) if mode == "vcg" else np.full(R, np.nan)
            sig, direction = res.pvalue < C.ALPHA, np.sign(res.statistic)
            for k in range(R):
                for r, x in enumerate(doses):
                    out.append(base + [r + 1, x, mode, level, band, k, n0, int(len(pool)), int(pool_s),
                                       ref.pvalue[r], res.pvalue[k, r], bool(sig_conc[r]), bool(sig[k, r]),
                                       int(dir_conc[r]), int(direction[k, r]), float(p_cv[k])])
    return out


# ----------------------------------------------------------------------------
# incidence endpoints
# ----------------------------------------------------------------------------

def _incidence_matrix(roster: pd.DataFrame, lesions: pd.DataFrame, findings: list[tuple[str, str]]) -> np.ndarray:
    """animals x findings boolean matrix aligned with ``roster`` rows."""
    key = roster["animal"].to_numpy()
    pos = {a: i for i, a in enumerate(key)}
    mat = np.zeros((len(roster), len(findings)), dtype=bool)
    col = {f: j for j, f in enumerate(findings)}
    sub = lesions[lesions["animal"].isin(pos)]
    for animal, tissue, morph in sub[["animal", "tissue", "morphology"]].itertuples(index=False):
        j = col.get((tissue, morph))
        if j is not None:
            mat[pos[animal], j] = True
    return mat


def evaluate_incidence(unit: Unit) -> list[list]:
    rng = np.random.default_rng([unit.seed, M.stable_key(unit.target.study, unit.target.sex), 7])
    out: list[list] = []
    tgt = unit.target
    roster = unit.roster.reset_index(drop=True)
    if roster.empty or unit.lesions.empty:
        return out
    neoplastic = tgt.design == "chronic"
    les = unit.lesions[unit.lesions["neoplastic"].eq("Neoplastic") == neoplastic]
    counts = les.groupby(["tissue", "morphology"])["animal"].nunique()
    findings = [f for f, n in counts.items() if n >= C.MIN_AFFECTED_FOR_INCIDENCE]
    if not findings:
        return out
    domain = "neoplasm" if neoplastic else "nonneoplastic lesion"

    mat = _incidence_matrix(roster, les, findings)
    ctrl_rows = roster["is_control"].to_numpy()
    n0 = int(ctrl_rows.sum())
    if n0 < 3:
        return out
    doses = sorted(roster.loc[~ctrl_rows, "dose"].unique())
    dose_rows = [(roster["dose"] == x).to_numpy() & ~ctrl_rows for x in doses]
    dose_rows = [(x, m) for x, m in zip(doses, dose_rows) if m.sum() >= 2]
    if not dose_rows:
        return out
    t_days = roster["days"].to_numpy(dtype=float)
    t_max = roster["t_max"].to_numpy(dtype=float)

    hist = unit.hist_roster.reset_index(drop=True)
    if hist.empty:
        return out
    hmat = _incidence_matrix(hist, unit.hist_lesions, findings)
    h_days = hist["days"].to_numpy(dtype=float)
    h_tmax = hist["t_max"].to_numpy(dtype=float)
    level_masks = {lvl: M.pool_mask(hist, tgt, crit) for lvl, crit in C.MATCH_LEVELS.items()}
    drift_base = level_masks[C.DRIFT_LEVEL]
    band_masks = {f"{lo}-{hi}": drift_base & M.year_gap_mask(hist, tgt, (lo, hi)) for lo, hi in C.DRIFT_BANDS}

    def test(aff_t, days_t, tmax_t, aff_c, days_c, tmax_c):
        if neoplastic:
            w_t = S.poly_k_weights(days_t, aff_t, tmax_t)
            w_c = S.poly_k_weights(days_c, aff_c, tmax_c)
            z, p = S.poly3_test(aff_t, w_t, aff_c, w_c)
            return p, np.sign(z) if np.isfinite(z) else 0
        p = S.fisher_exact_2x2(aff_t.sum(), aff_t.size, aff_c.sum(), aff_c.size)[0]
        return p, np.sign(aff_t.mean() - aff_c.mean())

    for j, (tissue, morph) in enumerate(findings):
        endpoint = f"{tissue}: {morph}"
        aff_all = mat[:, j]
        ref = [test(aff_all[m], t_days[m], t_max[m], aff_all[ctrl_rows], t_days[ctrl_rows], t_max[ctrl_rows])
               for _, m in dose_rows]
        base = [tgt.study, tgt.sex, tgt.species, tgt.design, domain, endpoint]
        ctrl_idx = np.flatnonzero(ctrl_rows)

        def configs():
            for lvl, m in level_masks.items():
                yield "vcg", lvl, "", m
                yield "partial", lvl, "", m
            for band, m in band_masks.items():
                yield "vcg", C.DRIFT_LEVEL, band, m

        for mode, level, band, mask in configs():
            pool_idx = np.flatnonzero(mask)
            n_virtual = n0 if mode == "vcg" else n0 - n0 // 2
            if len(pool_idx) < max(n_virtual, 2 * n0 if mode == "vcg" else n_virtual):
                continue
            pool_s = int(hist.loc[pool_idx, "study"].nunique())
            for k in range(unit.n_draws):
                pick = rng.choice(pool_idx, size=n_virtual, replace=False)
                aff_c = hmat[pick, j]
                days_c, tmax_c = h_days[pick], h_tmax[pick]
                if mode == "partial":
                    keep = rng.choice(ctrl_idx, size=n0 - n_virtual, replace=False)
                    aff_c = np.concatenate([aff_all[keep], aff_c])
                    days_c = np.concatenate([t_days[keep], days_c])
                    tmax_c = np.concatenate([t_max[keep], tmax_c])
                for r, (x, m) in enumerate(dose_rows):
                    p, direction = test(aff_all[m], t_days[m], t_max[m], aff_c, days_c, tmax_c)
                    p_ref, dir_ref = ref[r]
                    out.append(base + [r + 1, x, mode, level, band, k, n0, int(len(pool_idx)), pool_s,
                                       p_ref, p, bool(p_ref < C.ALPHA), bool(p < C.ALPHA),
                                       int(dir_ref), int(direction), np.nan])
    return out


# ----------------------------------------------------------------------------
# orchestration
# ----------------------------------------------------------------------------

def _controls_with_metadata(frame: pd.DataFrame, studies: pd.DataFrame) -> pd.DataFrame:
    cols = ["species", "strain", "design", "route", "vehicle", "laboratory", "year", "dtxsid"]
    meta = studies.set_index("study")[cols]
    base = frame.loc[frame["is_control"], [c for c in frame.columns if c not in cols]]
    ctrl = base.merge(meta, left_on="study", right_index=True, how="inner")
    return ctrl.reset_index(drop=True)


def _by_key(frame: pd.DataFrame, keys: list[str]) -> dict:
    return {k: g for k, g in frame.groupby(keys, sort=False)} if len(frame) else {}


def build_units(data: dict[str, pd.DataFrame], n_draws: int = C.N_DRAWS, seed: int = C.SEED,
                limit: int | None = None) -> list[Unit]:
    studies, meas, roster, lesions = (data[k] for k in ("studies", "measurements", "animals", "lesions"))
    hist_controls = _controls_with_metadata(meas, studies)
    hist_roster = _controls_with_metadata(roster, studies)
    hist_lesions = lesions[lesions["animal"].isin(set(hist_roster["animal"]))]

    meas_by = _by_key(meas, ["study", "sex"])
    roster_by = _by_key(roster, ["study", "sex"])
    lesions_by = _by_key(lesions, ["study", "sex"])
    # candidate pools share species/sex/design with the target; levels only narrow further
    hc_by = _by_key(hist_controls, ["species", "sex", "design"])
    hr_by = _by_key(hist_roster, ["species", "sex", "design"])
    hl_by = {k: hist_lesions[hist_lesions["animal"].isin(set(g["animal"]))] for k, g in hr_by.items()}

    meta = studies.set_index("study")
    pairs = pd.concat([meas[["study", "sex"]], roster[["study", "sex"]]]).drop_duplicates()
    pairs = pairs[pairs["study"].isin(meta.index)].sort_values(["study", "sex"])
    if limit:
        pairs = pairs.head(limit)
    empty = pd.DataFrame()
    units: list[Unit] = []
    for study, sex in pairs.itertuples(index=False):
        row = meta.loc[study]
        tgt = M.target_from_row(pd.Series({"study": study, **row.to_dict()}), sex)
        key = (tgt.species, sex, tgt.design)
        hc = hc_by.get(key, empty)
        hr = hr_by.get(key, empty)
        if len(hc):
            hc = hc[hc["study"] != study]
        if len(hr):
            hr = hr[hr["study"] != study]
        l = lesions_by.get((study, sex), empty)
        hl = hl_by.get(key, empty) if len(l) else empty
        units.append(Unit(tgt, meas_by.get((study, sex), empty), roster_by.get((study, sex), empty), l,
                          hc, hr, hl, n_draws, seed))
    return units


_CATEGORICAL = ["study", "sex", "species", "design", "domain", "endpoint", "mode", "level", "band"]
_SMALL_INT = ["dose_rank", "draw", "n_control", "pool_studies", "dir_conc", "dir_vcg"]
_FLOAT32 = ["dose", "p_conc", "p_vcg", "p_ctrl_vs_vcg"]


def _compact(rows: list[list]) -> pd.DataFrame:
    """Turn a unit's records into a memory-friendly frame (20M+ rows fit in a few GB)."""
    df = pd.DataFrame(rows, columns=RECORD_COLUMNS)
    for c in _CATEGORICAL:
        df[c] = df[c].astype("category")
    for c in _SMALL_INT:
        df[c] = df[c].astype("int16")
    df["pool_animals"] = df["pool_animals"].astype("int32")
    for c in _FLOAT32:
        df[c] = df[c].astype("float32")
    return df


def _run_unit(unit: Unit) -> pd.DataFrame:
    rows = evaluate_continuous(unit) if len(unit.continuous) else []
    rows += evaluate_incidence(unit) if len(unit.roster) and len(unit.lesions) else []
    return _compact(rows)


def run(data: dict[str, pd.DataFrame], workers: int = 1, n_draws: int = C.N_DRAWS,
        seed: int = C.SEED, limit: int | None = None) -> pd.DataFrame:
    units = build_units(data, n_draws=n_draws, seed=seed, limit=limit)
    log.info("evaluating %d study/sex units with %d workers", len(units), workers)
    parts: list[pd.DataFrame] = []
    n = 0
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for i, part in enumerate(pool.map(_run_unit, units, chunksize=4)):
                parts.append(part)
                n += len(part)
                if i % 50 == 0:
                    log.info("  %d / %d units done, %d records", i + 1, len(units), n)
    else:
        for i, u in enumerate(units):
            part = _run_unit(u)
            parts.append(part)
            n += len(part)
            if i % 50 == 0:
                log.info("  %d / %d units done, %d records", i + 1, len(units), n)
    out = pd.concat(parts, ignore_index=True)
    for c in _CATEGORICAL:
        if out[c].dtype != "category":          # concat widens mixed categories to object
            out[c] = out[c].astype("category")
    return out


# ----------------------------------------------------------------------------
# summaries
# ----------------------------------------------------------------------------

def endpoint_level(records: pd.DataFrame) -> pd.DataFrame:
    """Collapse dose groups: an endpoint is 'significant' if any dose group is."""
    keys = ["study", "sex", "species", "design", "domain", "endpoint", "mode", "level", "band", "draw"]
    agg = records.groupby(keys, sort=False, observed=True).agg(
        sig_conc=("sig_conc", "any"), sig_vcg=("sig_vcg", "any"),
        n_control=("n_control", "first"), pool_animals=("pool_animals", "first"),
        pool_studies=("pool_studies", "first"), p_ctrl_vs_vcg=("p_ctrl_vs_vcg", "first")).reset_index()
    return agg


def rates(df: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """Agreement, false positive and false negative rates of VCG conclusions
    relative to the concurrent control, averaged over draws within endpoint first."""
    d = df.copy()
    d["agree"] = (d["sig_conc"] == d["sig_vcg"]).astype(float)
    d["fp"] = (~d["sig_conc"] & d["sig_vcg"]).astype(float)
    d["fn"] = (d["sig_conc"] & ~d["sig_vcg"]).astype(float)
    unit_keys = ["study", "sex", "endpoint", "mode", "level", "band"]
    per_endpoint = d.groupby([*unit_keys, *[b for b in by if b not in unit_keys]], sort=False, observed=True).agg(
        agree=("agree", "mean"), fp=("fp", "mean"), fn=("fn", "mean"),
        sig_conc=("sig_conc", "mean"), sig_vcg=("sig_vcg", "mean")).reset_index()
    g = per_endpoint.groupby(by, sort=False, observed=True)
    out = g.agg(n_endpoints=("agree", "size"), agreement=("agree", "mean"),
                false_positive=("fp", "mean"), false_negative=("fn", "mean"),
                conc_sig_rate=("sig_conc", "mean"), vcg_sig_rate=("sig_vcg", "mean")).reset_index()
    # Chance-corrected agreement (Cohen's kappa) and the agreement a constant
    # "never significant" verdict would achieve; both put the raw agreement in context
    pe = out["conc_sig_rate"] * out["vcg_sig_rate"] + (1 - out["conc_sig_rate"]) * (1 - out["vcg_sig_rate"])
    out["kappa"] = (out["agreement"] - pe) / (1 - pe).where((1 - pe) > 1e-9, np.nan)
    out["trivial_agreement"] = np.maximum(out["conc_sig_rate"], 1 - out["conc_sig_rate"])
    # conditional rates: FP among endpoints not significant with concurrent control, FN among significant
    cond = per_endpoint.assign(nonsig=lambda x: 1 - x["sig_conc"]).groupby(by, sort=False, observed=True).apply(
        lambda x: pd.Series({
            "fp_given_conc_nonsig": x["fp"].sum() / max(x["nonsig"].sum(), 1e-9),
            "fn_given_conc_sig": x["fn"].sum() / max(x["sig_conc"].sum(), 1e-9)}),
        include_groups=False).reset_index()
    out = out.merge(cond, on=by)
    out["sensitivity"] = 1 - out["fn_given_conc_sig"]      # concurrent-significant findings also found with the VCG
    out["specificity"] = 1 - out["fp_given_conc_nonsig"]   # concurrent-nonsignificant endpoints also quiet with the VCG
    return out


def summarise(records: pd.DataFrame) -> dict[str, pd.DataFrame]:
    ep = endpoint_level(records)
    main = ep[(ep["band"] == "")]
    tables = {
        "agreement_by_level": rates(main, ["mode", "level"]),
        "agreement_by_domain_level": rates(main, ["mode", "level", "domain"]),
        "agreement_by_species_design": rates(main, ["mode", "level", "species", "design"]),
        "agreement_by_endpoint": rates(main[(main["mode"] == "vcg")], ["level", "domain", "endpoint"]),
        "agreement_by_drift": rates(ep[ep["band"] != ""], ["mode", "level", "band"]),
        "agreement_by_drift_domain": rates(ep[ep["band"] != ""], ["mode", "level", "band", "domain"]),
        "dose_group_level": rates(records[records["band"] == ""], ["mode", "level"]),
    }
    # coverage: how many study/sex/endpoint units could be served at each level
    universe = main[["study", "sex", "endpoint", "domain"]].drop_duplicates().shape[0]
    cov = (main[main["mode"] == "vcg"].groupby("level", observed=True)[["study", "sex", "endpoint"]]
           .apply(lambda x: x.drop_duplicates().shape[0]).rename("n_units").reset_index())
    cov["coverage"] = cov["n_units"] / universe
    tables["coverage_by_level"] = cov
    # control comparability: how often the concurrent controls differ from the VCG itself
    comp = main[(main["mode"] == "vcg") & main["p_ctrl_vs_vcg"].notna()].copy()
    comp["differs"] = comp["p_ctrl_vs_vcg"] < C.ALPHA
    tables["control_vs_vcg"] = comp.groupby(["level", "domain"], observed=True)["differs"].mean().rename(
        "share_concurrent_differs_from_vcg").reset_index()
    return tables
