"""Pre-specified comparison of ways to build and analyse a virtual control group.

See ``docs/analysis_plan.md``.  All methods share the matching level
``config.METHOD_LEVEL``; they differ in which historical animals are used and
how the comparison is analysed:

=================  ===========================================================
``baseline``       current method: individual animals sampled uniformly
``bw_matched``     individual animals sampled with pre-dosing body weight weights
``bw_ancova``      ``bw_matched`` animals, ANCOVA on pre-dosing body weight
``whole_group``    intact historical control groups
``hierarchical``   random-effects MAP prior (full) / robust MAP (partial)
``baseline_trend`` Cochran-Armitage trend test with the ``baseline`` VCG
``tarone``         Tarone's trend test with a beta-binomial historical prior
=================  ===========================================================

Each worker returns one row per study x sex x endpoint x method x mode, with the
draws already averaged.  ``compare`` turns those rows into metric tables and
paired, study-cluster bootstrap comparisons with ``baseline``.
"""
from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from . import config as C
from . import evaluation as E
from . import matching as M
from . import stats as S

log = logging.getLogger(__name__)

ENDPOINT_COLUMNS = ["set", "study", "sex", "species", "design", "domain", "endpoint", "method", "mode",
                    "n_draws", "n_control", "pool_animals", "pool_groups",
                    "sig_conc", "sig_vcg", "agree", "fp", "fn", "ctrl_differs"]
METHOD_ORDER = ["baseline", "bw_matched", "bw_ancova", "whole_group", "hierarchical", "baseline_trend", "tarone"]
COMPARATOR = {"bw_matched": "baseline", "bw_ancova": "baseline", "whole_group": "baseline",
              "hierarchical": "baseline", "tarone": "baseline_trend"}
# Post-hoc, not in the analysis plan: is a partial-mode gain due to the borrowed
# information, or only to dropping the virtual animals?
POSTHOC_PAIRS = [("half_concurrent", "baseline"), ("hierarchical", "half_concurrent"),
                 ("half_concurrent_trend", "baseline_trend"), ("tarone", "half_concurrent_trend")]
SETS = ("development", "holdout")


def set_of(year: float) -> str | None:
    """Development or held-out set by the target's start year; None if unknown."""
    if year is None or not np.isfinite(year):
        return None
    return "holdout" if year >= C.CUTOFF_YEAR else "development"


def _row(base: list, method: str, mode: str, sig_ref: bool, sig_draws: np.ndarray,
         p_check: np.ndarray | None, n_control: int, pool_animals: int, pool_groups: int) -> list:
    sig_draws = np.asarray(sig_draws, dtype=bool)
    if p_check is not None and np.isfinite(p_check).any():
        differs = float(np.nanmean(np.asarray(p_check, dtype=float) < C.ALPHA))
    else:
        differs = np.nan
    return base + [method, mode, int(sig_draws.size), int(n_control), int(pool_animals), int(pool_groups),
                   float(sig_ref), float(sig_draws.mean()), float((sig_draws == sig_ref).mean()),
                   float((sig_draws & (not sig_ref)).mean()), float((~sig_draws & sig_ref).mean()), differs]


def _rng(unit: E.Unit, *tags: object) -> np.random.Generator:
    return np.random.default_rng([unit.seed, M.stable_key(unit.target.study, unit.target.sex, *tags)])


# ----------------------------------------------------------------------------
# continuous endpoints
# ----------------------------------------------------------------------------

def _two_sided_normal_p(z: float) -> float:
    from scipy import stats
    return float(2 * stats.norm.sf(abs(z)))


def evaluate_continuous(unit: E.Unit, set_label: str) -> list[list]:
    tgt = unit.target
    rows: list[list] = []
    hc = unit.hist_controls
    if hc.empty or unit.continuous.empty:
        return rows
    R = unit.n_draws
    for (domain, endpoint), d in unit.continuous.groupby(["domain", "endpoint"], sort=False):
        is_ctrl = d["is_control"].to_numpy(dtype=bool)
        ctrl = d.loc[is_ctrl, "value"].to_numpy(dtype=float)
        ctrl_bw = d.loc[is_ctrl, "bw0"].to_numpy(dtype=float)
        treated = d.loc[~is_ctrl]
        doses = [x for x in sorted(treated["dose"].unique()) if (treated["dose"] == x).sum() >= 2]
        groups = [treated.loc[treated["dose"] == x, "value"].to_numpy(dtype=float) for x in doses]
        groups_bw = [treated.loc[treated["dose"] == x, "bw0"].to_numpy(dtype=float) for x in doses]
        n0 = ctrl.size
        if n0 < 3 or not groups:
            continue
        h = hc[hc["endpoint"] == endpoint]
        if h.empty:
            continue
        hv = h["value"].to_numpy(dtype=float)
        hb = h["bw0"].to_numpy(dtype=float)
        hg = h["study"].to_numpy()
        n_groups = len(np.unique(hg))
        base = [set_label, tgt.study, tgt.sex, tgt.species, tgt.design, domain, endpoint]
        rng = _rng(unit, endpoint)
        sig_ref = bool((S.dunnett(groups, ctrl).pvalue < C.ALPHA).any())
        n_keep = n0 // 2                 # partial: concurrent animals kept
        n_virt = n0 - n_keep             # partial: virtual animals added

        def dunnett_rows(method: str, controls: list[np.ndarray] | np.ndarray, mode: str, pool_n: int,
                         pool_g: int, check: bool) -> None:
            if isinstance(controls, np.ndarray):          # equal sizes: one vectorised call
                res = S.dunnett_many(groups, controls)
                sig = (res.pvalue < C.ALPHA).any(axis=1)
                p_cv = S.welch_p_many(ctrl, controls) if check else None
            else:
                sig = np.array([(S.dunnett(groups, c).pvalue < C.ALPHA).any() for c in controls])
                p_cv = np.array([S.welch_p_many(ctrl, c[None, :])[0] for c in controls]) if check else None
            rows.append(_row(base, method, mode, sig_ref, sig, p_cv, n0, pool_n, pool_g))

        # --- baseline: uniform individual sampling --------------------------------
        if hv.size >= 2 * n0:
            draws = np.stack([hv[rng.choice(hv.size, n0, replace=False)] for _ in range(R)])
            dunnett_rows("baseline", draws, "vcg", hv.size, n_groups, check=True)
        if hv.size >= n_virt:
            draws = np.stack([np.r_[rng.choice(ctrl, n_keep, replace=False),
                                    hv[rng.choice(hv.size, n_virt, replace=False)]] for _ in range(R)])
            dunnett_rows("baseline", draws, "partial", hv.size, n_groups, check=False)

        # --- body-weight matched sampling ----------------------------------------
        w = M.bw_weights(hb, d["bw0"].to_numpy(dtype=float))
        n_elig = int((w > 0).sum()) if w is not None else 0
        g_elig = len(np.unique(hg[w > 0])) if n_elig else 0
        if w is not None and n_elig >= 2 * n0:
            draws = np.stack([hv[M.draw_weighted(w, n0, rng)] for _ in range(R)])
            dunnett_rows("bw_matched", draws, "vcg", n_elig, g_elig, check=True)
        if w is not None and n_elig >= n_virt:
            draws = np.stack([np.r_[rng.choice(ctrl, n_keep, replace=False),
                                    hv[M.draw_weighted(w, n_virt, rng)]] for _ in range(R)])
            dunnett_rows("bw_matched", draws, "partial", n_elig, g_elig, check=False)

        # --- ANCOVA on pre-dosing body weight (weight endpoints only) --------------
        if domain in C.ANCOVA_DOMAINS and w is not None:
            rows.extend(_ancova_rows(base, rng, R, ctrl, ctrl_bw, groups, groups_bw, hv, hb, w, hg))

        # --- whole historical control groups -----------------------------------
        if n_groups >= C.MIN_HIST_GROUPS and hv.size >= 2 * n0:
            draws = [hv[M.draw_whole_groups(hg, n0, rng)] for _ in range(R)]
            dunnett_rows("whole_group", draws, "vcg", hv.size, n_groups, check=True)
        if n_groups >= C.MIN_HIST_GROUPS and hv.size >= n_virt:
            draws = [np.r_[rng.choice(ctrl, n_keep, replace=False), hv[M.draw_whole_groups(hg, n_virt, rng)]]
                     for _ in range(R)]
            dunnett_rows("whole_group", draws, "partial", hv.size, n_groups, check=False)

        # --- hierarchical borrowing --------------------------------------------
        hist_groups = [hv[hg == g] for g in np.unique(hg)]
        hist_groups = [g for g in hist_groups if g.size >= 2]
        if len(hist_groups) >= C.MIN_HIST_GROUPS:
            re, s2_hist = S.historical_random_effects(hist_groups)
            res = S.map_predictive_test(groups, re, C.HIER_MAX_DF)
            if np.isfinite(res.pvalue).all():
                z = (ctrl.mean() - re.mu) / np.sqrt(ctrl.var(ddof=1) / n0 + re.predictive_var)
                rows.append(_row(base, "hierarchical", "vcg", sig_ref, np.array([(res.pvalue < C.ALPHA).any()]),
                                 np.array([_two_sided_normal_p(z)]), n0, hv.size, len(hist_groups)))
            sig = []
            for _ in range(R):
                kept = rng.choice(ctrl, n_keep, replace=False)
                r = S.robust_map_test(groups, kept, re, s2_hist, C.ROBUST_MAP_WEIGHT, C.HIER_MAX_DF)
                sig.append((r.pvalue < C.ALPHA).any())
            rows.append(_row(base, "hierarchical", "partial", sig_ref, np.array(sig), None, n0, hv.size,
                             len(hist_groups)))
    return rows


def _ancova_rows(base: list, rng: np.random.Generator, R: int, ctrl: np.ndarray, ctrl_bw: np.ndarray,
                 groups: list[np.ndarray], groups_bw: list[np.ndarray], hv: np.ndarray, hb: np.ndarray,
                 w: np.ndarray, hg: np.ndarray) -> list[list]:
    """ANCOVA variant: animals without a pre-dosing weight are left out of both
    the reference analysis and the VCG analysis."""
    ok = np.isfinite(ctrl_bw)
    c, cb = ctrl[ok], ctrl_bw[ok]
    gs = [(g[np.isfinite(b)], b[np.isfinite(b)]) for g, b in zip(groups, groups_bw)]
    gs = [(g, b) for g, b in gs if g.size >= 2]
    n0 = c.size
    if n0 < 3 or not gs:
        return []
    ys, xs = [g for g, _ in gs], [b for _, b in gs]
    sig_ref = bool((S.ancova_dunnett(ys, xs, c, cb).pvalue < C.ALPHA).any())
    eligible = w > 0
    n_elig = int(eligible.sum())
    g_elig = len(np.unique(hg[eligible])) if n_elig else 0
    out = []
    if n_elig >= 2 * n0:
        sig, p_cv = [], []
        for _ in range(R):
            idx = M.draw_weighted(w, n0, rng)
            sig.append((S.ancova_dunnett(ys, xs, hv[idx], hb[idx]).pvalue < C.ALPHA).any())
            p_cv.append(S.ancova_dunnett([hv[idx]], [hb[idx]], c, cb).pvalue[0])
        out.append(_row(base, "bw_ancova", "vcg", sig_ref, np.array(sig), np.array(p_cv), n0, n_elig, g_elig))
    n_keep = n0 // 2
    n_virt = n0 - n_keep
    if n_elig >= n_virt:
        sig = []
        for _ in range(R):
            keep = rng.choice(n0, n_keep, replace=False)
            idx = M.draw_weighted(w, n_virt, rng)
            sig.append((S.ancova_dunnett(ys, xs, np.r_[c[keep], hv[idx]], np.r_[cb[keep], hb[idx]]).pvalue
                        < C.ALPHA).any())
        out.append(_row(base, "bw_ancova", "partial", sig_ref, np.array(sig), None, n0, n_elig, g_elig))
    return out


# ----------------------------------------------------------------------------
# incidence endpoints
# ----------------------------------------------------------------------------

def evaluate_incidence(unit: E.Unit, set_label: str) -> list[list]:
    tgt = unit.target
    rows: list[list] = []
    roster = unit.roster.reset_index(drop=True)
    hist = unit.hist_roster.reset_index(drop=True)
    if roster.empty or unit.lesions.empty or hist.empty:
        return rows
    neoplastic = tgt.design == "chronic"
    les = unit.lesions[unit.lesions["neoplastic"].eq("Neoplastic") == neoplastic]
    counts = les.groupby(["tissue", "morphology"])["animal"].nunique()
    findings = [f for f, n in counts.items() if n >= C.MIN_AFFECTED_FOR_INCIDENCE]
    if not findings:
        return rows
    domain = "neoplasm" if neoplastic else "nonneoplastic lesion"
    mat = E._incidence_matrix(roster, les, findings)
    is_ctrl = roster["is_control"].to_numpy(dtype=bool)
    ctrl_idx = np.flatnonzero(is_ctrl)
    n0 = ctrl_idx.size
    if n0 < 3:
        return rows
    doses = sorted(roster.loc[~is_ctrl, "dose"].unique())
    dose_rows = [(roster["dose"] == x).to_numpy() & ~is_ctrl for x in doses]
    dose_rows = [m for m in dose_rows if m.sum() >= 2]
    if not dose_rows:
        return rows
    t_days = roster["days"].to_numpy(dtype=float)
    t_max = roster["t_max"].to_numpy(dtype=float)
    hmat = E._incidence_matrix(hist, unit.hist_lesions, findings)
    h_days = hist["days"].to_numpy(dtype=float)
    h_tmax = hist["t_max"].to_numpy(dtype=float)
    hg = hist["study"].to_numpy()
    g_labels, g_inv = np.unique(hg, return_inverse=True)
    n_groups = g_labels.size
    rng = _rng(unit, "incidence")
    R = unit.n_draws
    n_keep = n0 // 2
    n_virt = n0 - n_keep

    # draw the virtual controls once per method/mode: a VCG is one set of animals
    # used for every finding of the study
    draws: dict[tuple[str, str], list[tuple[np.ndarray, np.ndarray]]] = {}

    def add(method: str, mode: str, pick) -> None:
        out = []
        for _ in range(R):
            keep = rng.choice(ctrl_idx, n_keep, replace=False) if mode == "partial" else ctrl_idx[:0]
            out.append((keep, pick(n0 if mode == "vcg" else n_virt)))
        draws[(method, mode)] = out

    n_hist = hist.shape[0]
    for mode in ("vcg", "partial"):
        need = 2 * n0 if mode == "vcg" else n_virt
        if n_hist >= need:
            add("baseline", mode, lambda n: rng.choice(n_hist, n, replace=False))
    w = M.bw_weights(hist["bw0"].to_numpy(dtype=float), roster["bw0"].to_numpy(dtype=float))
    n_elig = int((w > 0).sum()) if w is not None else 0
    for mode in ("vcg", "partial"):
        need = 2 * n0 if mode == "vcg" else n_virt
        if w is not None and n_elig >= need:
            add("bw_matched", mode, lambda n: M.draw_weighted(w, n, rng))
    for mode in ("vcg", "partial"):
        need = 2 * n0 if mode == "vcg" else n_virt
        if n_groups >= C.MIN_HIST_GROUPS and n_hist >= need:
            add("whole_group", mode, lambda n: M.draw_whole_groups(hg, n, rng))
    pool_info = {"baseline": (n_hist, n_groups), "whole_group": (n_hist, n_groups),
                 "bw_matched": (n_elig, len(np.unique(hg[w > 0])) if n_elig else 0)}

    def pairwise_any(aff_c, days_c, tmax_c, aff_all) -> bool:
        for m in dose_rows:
            if neoplastic:
                z, p = S.poly3_test(aff_all[m], S.poly_k_weights(t_days[m], aff_all[m], t_max[m]),
                                    aff_c, S.poly_k_weights(days_c, aff_c, tmax_c))
            else:
                p = S.fisher_exact_2x2(aff_all[m].sum(), m.sum(), aff_c.sum(), aff_c.size)[0]
            if p < C.ALPHA:
                return True
        return False

    scores = np.arange(len(dose_rows) + 1, dtype=float)
    n_treated = np.array([m.sum() for m in dose_rows], dtype=float)
    for j, (tissue, morph) in enumerate(findings):
        base = [set_label, tgt.study, tgt.sex, tgt.species, tgt.design, domain, f"{tissue}: {morph}"]
        aff = mat[:, j]
        haff = hmat[:, j]
        sig_ref = pairwise_any(aff[ctrl_idx], t_days[ctrl_idx], t_max[ctrl_idx], aff)
        for (method, mode), ds in draws.items():
            sig = []
            for keep, pick in ds:
                aff_c = np.r_[aff[keep], haff[pick]]
                days_c = np.r_[t_days[keep], h_days[pick]]
                tmax_c = np.r_[t_max[keep], h_tmax[pick]]
                sig.append(pairwise_any(aff_c, days_c, tmax_c, aff))
            pn, pg = pool_info[method]
            rows.append(_row(base, method, mode, sig_ref, np.array(sig), None, n0, pn, pg))

        if not neoplastic:
            continue
        # --- trend tests: Cochran-Armitage reference, baseline VCG, Tarone --------
        x_treated = np.array([aff[m].sum() for m in dose_rows], dtype=float)
        _, p_ref = S.cochran_armitage(np.r_[aff[ctrl_idx].sum(), x_treated], np.r_[n0, n_treated], scores)
        sig_trend_ref = bool(p_ref < C.ALPHA)
        for mode in ("vcg", "partial"):
            ds = draws.get(("baseline", mode))
            if ds is None:
                continue
            sig = []
            for keep, pick in ds:
                xc = aff[keep].sum() + haff[pick].sum()
                sig.append(S.cochran_armitage(np.r_[xc, x_treated], np.r_[keep.size + pick.size, n_treated],
                                              scores)[1] < C.ALPHA)
            rows.append(_row(base, "baseline_trend", mode, sig_trend_ref, np.array(sig), None, n0, n_hist, n_groups))
        if n_groups >= C.MIN_HIST_GROUPS:
            x_h = np.bincount(g_inv, weights=haff.astype(float), minlength=n_groups)
            n_h = np.bincount(g_inv, minlength=n_groups).astype(float)
            a, b = S.fit_beta_binomial(x_h, n_h)
            p_full = S.tarone_trend(np.r_[0.0, x_treated], np.r_[0.0, n_treated], scores, a, b)[1]
            rows.append(_row(base, "tarone", "vcg", sig_trend_ref, np.array([p_full < C.ALPHA]), None,
                             n0, n_hist, n_groups))
            sig = []
            for _ in range(R):
                keep = rng.choice(ctrl_idx, n_keep, replace=False)
                p = S.tarone_trend(np.r_[aff[keep].sum(), x_treated], np.r_[n_keep, n_treated], scores, a, b)[1]
                sig.append(p < C.ALPHA)
            rows.append(_row(base, "tarone", "partial", sig_trend_ref, np.array(sig), None, n0, n_hist, n_groups))
    return rows


# ----------------------------------------------------------------------------
# post-hoc control: the kept concurrent half alone, without any borrowing
# ----------------------------------------------------------------------------

def evaluate_half_concurrent(unit: E.Unit, set_label: str) -> list[list]:
    """Post-hoc (not in the analysis plan): analyse the partial design using only
    the kept half of the concurrent controls, without any virtual animals.

    This separates two effects of the partial VCG methods: removing half of the
    concurrent animals, and what replaces them.  Separate random streams keep the
    pre-specified results unchanged.  Methods: ``half_concurrent`` (Dunnett,
    Fisher, poly-3) and ``half_concurrent_trend`` (Cochran-Armitage).
    """
    tgt = unit.target
    rows: list[list] = []
    R = unit.n_draws
    if not unit.continuous.empty and not unit.hist_controls.empty:
        hist_endpoints = set(unit.hist_controls["endpoint"])
        for (domain, endpoint), d in unit.continuous.groupby(["domain", "endpoint"], sort=False):
            if endpoint not in hist_endpoints:
                continue
            is_ctrl = d["is_control"].to_numpy(dtype=bool)
            ctrl = d.loc[is_ctrl, "value"].to_numpy(dtype=float)
            treated = d.loc[~is_ctrl]
            doses = [x for x in sorted(treated["dose"].unique()) if (treated["dose"] == x).sum() >= 2]
            groups = [treated.loc[treated["dose"] == x, "value"].to_numpy(dtype=float) for x in doses]
            n0 = ctrl.size
            if n0 < 3 or not groups or n0 // 2 < 2:
                continue
            rng = _rng(unit, endpoint, "half_concurrent")
            sig_ref = bool((S.dunnett(groups, ctrl).pvalue < C.ALPHA).any())
            kept = np.stack([rng.choice(ctrl, n0 // 2, replace=False) for _ in range(R)])
            sig = (S.dunnett_many(groups, kept).pvalue < C.ALPHA).any(axis=1)
            base = [set_label, tgt.study, tgt.sex, tgt.species, tgt.design, domain, endpoint]
            rows.append(_row(base, "half_concurrent", "partial", sig_ref, sig, None, n0, 0, 0))

    roster = unit.roster.reset_index(drop=True)
    if roster.empty or unit.lesions.empty or unit.hist_roster.empty:
        return rows
    neoplastic = tgt.design == "chronic"
    les = unit.lesions[unit.lesions["neoplastic"].eq("Neoplastic") == neoplastic]
    counts = les.groupby(["tissue", "morphology"])["animal"].nunique()
    findings = [f for f, n in counts.items() if n >= C.MIN_AFFECTED_FOR_INCIDENCE]
    is_ctrl = roster["is_control"].to_numpy(dtype=bool)
    ctrl_idx = np.flatnonzero(is_ctrl)
    n0 = ctrl_idx.size
    doses = sorted(roster.loc[~is_ctrl, "dose"].unique())
    dose_rows = [m for m in ((roster["dose"] == x).to_numpy() & ~is_ctrl for x in doses) if m.sum() >= 2]
    if not findings or n0 < 3 or not dose_rows:
        return rows
    domain = "neoplasm" if neoplastic else "nonneoplastic lesion"
    mat = E._incidence_matrix(roster, les, findings)
    t_days = roster["days"].to_numpy(dtype=float)
    t_max = roster["t_max"].to_numpy(dtype=float)
    rng = _rng(unit, "incidence", "half_concurrent")
    keeps = [rng.choice(ctrl_idx, n0 // 2, replace=False) for _ in range(R)]
    scores = np.arange(len(dose_rows) + 1, dtype=float)
    n_treated = np.array([m.sum() for m in dose_rows], dtype=float)

    def pairwise_any(c_idx, aff) -> bool:
        for m in dose_rows:
            if neoplastic:
                _, p = S.poly3_test(aff[m], S.poly_k_weights(t_days[m], aff[m], t_max[m]),
                                    aff[c_idx], S.poly_k_weights(t_days[c_idx], aff[c_idx], t_max[c_idx]))
            else:
                p = S.fisher_exact_2x2(aff[m].sum(), m.sum(), aff[c_idx].sum(), c_idx.size)[0]
            if p < C.ALPHA:
                return True
        return False

    for j, (tissue, morph) in enumerate(findings):
        aff = mat[:, j]
        base = [set_label, tgt.study, tgt.sex, tgt.species, tgt.design, domain, f"{tissue}: {morph}"]
        sig = np.array([pairwise_any(k, aff) for k in keeps])
        rows.append(_row(base, "half_concurrent", "partial", pairwise_any(ctrl_idx, aff), sig, None, n0, 0, 0))
        if neoplastic:
            x_t = np.array([aff[m].sum() for m in dose_rows], dtype=float)
            ref = S.cochran_armitage(np.r_[aff[ctrl_idx].sum(), x_t], np.r_[n0, n_treated], scores)[1] < C.ALPHA
            sig = np.array([S.cochran_armitage(np.r_[aff[k].sum(), x_t], np.r_[k.size, n_treated], scores)[1]
                            < C.ALPHA for k in keeps])
            rows.append(_row(base, "half_concurrent_trend", "partial", bool(ref), sig, None, n0, 0, 0))
    return rows


# ----------------------------------------------------------------------------
# orchestration
# ----------------------------------------------------------------------------

def _with_bw(frame: pd.DataFrame, covariates: pd.DataFrame) -> pd.DataFrame:
    return frame.merge(covariates[["study", "sex", "animal", "bw0"]], on=["study", "sex", "animal"], how="left")


def build_units(data: dict[str, pd.DataFrame], n_draws: int = C.N_DRAWS, seed: int = C.SEED,
                limit: int | None = None) -> list[tuple[E.Unit, str]]:
    """Units of the current evaluation, with pre-dosing weights attached and the
    historical pools narrowed to ``METHOD_LEVEL`` once, up front."""
    cov = data["covariates"]
    data = {**data, "measurements": _with_bw(data["measurements"], cov), "animals": _with_bw(data["animals"], cov)}
    criteria = C.MATCH_LEVELS[C.METHOD_LEVEL]
    out = []
    for unit in E.build_units(data, n_draws=n_draws, seed=seed, limit=limit):
        label = set_of(unit.target.year)
        if label is None:
            continue
        if len(unit.hist_controls):
            unit.hist_controls = unit.hist_controls.loc[M.pool_mask(unit.hist_controls, unit.target, criteria)]
        if len(unit.hist_roster):
            unit.hist_roster = unit.hist_roster.loc[M.pool_mask(unit.hist_roster, unit.target, criteria)]
            if len(unit.hist_lesions):
                unit.hist_lesions = unit.hist_lesions[unit.hist_lesions["animal"].isin(set(unit.hist_roster["animal"]))]
        out.append((unit, label))
    return out


def _run_unit(item: tuple[E.Unit, str]) -> list[list]:
    unit, label = item
    return evaluate_continuous(unit, label) + evaluate_incidence(unit, label)


def _run_unit_posthoc(item: tuple[E.Unit, str]) -> list[list]:
    unit, label = item
    return evaluate_half_concurrent(unit, label)


def run(data: dict[str, pd.DataFrame], workers: int = 1, n_draws: int = C.N_DRAWS, seed: int = C.SEED,
        limit: int | None = None, posthoc: bool = False) -> pd.DataFrame:
    """Evaluate the pre-specified methods or, with ``posthoc``, only the
    half-concurrent control analysis."""
    worker = _run_unit_posthoc if posthoc else _run_unit
    items = build_units(data, n_draws=n_draws, seed=seed, limit=limit)
    log.info("method comparison: %d study/sex units (%d held-out) with %d workers", len(items),
             sum(lbl == "holdout" for _, lbl in items), workers)
    rows: list[list] = []
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for i, part in enumerate(pool.map(worker, items, chunksize=4)):
                rows.extend(part)
                if i % 100 == 0:
                    log.info("  %d / %d units done, %d endpoint rows", i + 1, len(items), len(rows))
    else:
        for i, item in enumerate(items):
            rows.extend(worker(item))
    df = pd.DataFrame(rows, columns=ENDPOINT_COLUMNS)
    for c in ["set", "study", "sex", "species", "design", "domain", "endpoint", "method", "mode"]:
        df[c] = df[c].astype("category")
    return df


# ----------------------------------------------------------------------------
# between-study variance of control groups
# ----------------------------------------------------------------------------

def implied_false_positive_rate(icc: np.ndarray, n: int = 10, sampling: str = "individual") -> np.ndarray:
    """False positive rate of one two-sided t test of a treated group against a
    virtual control group of the same size ``n``, when there is no treatment
    effect but a between-study share of variance ``icc`` = tau^2 / (tau^2 + sigma^2).

    The treated group carries its own study's offset (variance tau^2).

    * ``individual``: the VCG's ``n`` animals come from different studies, so its
      offset averages out but its animals vary by sigma^2 + tau^2.  The statistic
      has variance (2 sigma^2 + (n + 1) tau^2) / (2 sigma^2 + tau^2).
    * ``whole_group``: the VCG is one other study's control group with its own
      offset.  The statistic has variance 1 + n tau^2 / sigma^2.

    The statistic is taken as sqrt(var) times a t variable with 2n - 2 degrees
    of freedom, so the rate is exactly 5% when icc = 0.
    """
    from scipy import stats

    icc = np.asarray(icc, dtype=float)
    ratio = icc / (1 - icc)                      # tau^2 / sigma^2
    if sampling == "individual":
        var_z = (2 + (n + 1) * ratio) / (2 + ratio)
    elif sampling == "whole_group":
        var_z = 1 + n * ratio
    else:
        raise ValueError(sampling)
    df = 2 * n - 2
    return 2 * stats.t.sf(stats.t.ppf(0.975, df) / np.sqrt(var_z), df)


def between_study_variance(data: dict[str, pd.DataFrame], level: str = C.METHOD_LEVEL) -> pd.DataFrame:
    """Between-study share of control-group variance per endpoint, within strata
    of studies that match each other at ``level``.

    For each stratum (species, sex, design and the level's criteria) and endpoint
    with at least ``MIN_HIST_GROUPS`` control groups, DerSimonian-Laird gives
    tau^2 and the pooled within-group variance gives sigma^2.
    """
    criteria = [c for c in C.MATCH_LEVELS[level] if c != "period"]
    ctrl = E._controls_with_metadata(data["measurements"], data["studies"])
    keys = ["species", "sex", "design", *criteria, "domain", "endpoint"]
    rows = []
    for key, g in ctrl.groupby(keys, observed=True, dropna=True, sort=False):
        groups = [v.to_numpy(dtype=float) for _, v in g.groupby("study", observed=True)["value"]]
        groups = [v for v in groups if v.size >= 2]
        if len(groups) < C.MIN_HIST_GROUPS:
            continue
        re, s2 = S.historical_random_effects(groups)
        rows.append(dict(zip(keys, key), n_groups=len(groups), tau2=re.tau2, sigma2=s2,
                         icc=re.tau2 / (re.tau2 + s2) if re.tau2 + s2 > 0 else np.nan))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# summaries
# ----------------------------------------------------------------------------

def _metrics_from_sums(n, agree, pc, pv, fp, fn):
    """Metrics from per-group sums of per-endpoint (draw-averaged) quantities."""
    n = np.asarray(n, dtype=float)
    po, rc, rv = agree / n, pc / n, pv / n
    pe = rc * rv + (1 - rc) * (1 - rv)
    with np.errstate(divide="ignore", invalid="ignore"):
        kappa = np.where(1 - pe > 1e-12, (po - pe) / (1 - pe), np.nan)
        sens = np.where(pc > 0, 1 - fn / pc, np.nan)
        spec = np.where(n - pc > 0, 1 - fp / (n - pc), np.nan)
    return po, kappa, sens, spec, fp / n, fn / n, np.maximum(rc, 1 - rc)


def metrics(df: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    g = df.groupby(by, observed=True, sort=False)
    s = g.agg(n_units=("agree", "size"), agree=("agree", "sum"), pc=("sig_conc", "sum"), pv=("sig_vcg", "sum"),
              fp=("fp", "sum"), fn=("fn", "sum"), ctrl_differs=("ctrl_differs", "mean")).reset_index()
    po, kappa, sens, spec, fpr, fnr, trivial = _metrics_from_sums(s["n_units"], s["agree"], s["pc"], s["pv"],
                                                                  s["fp"], s["fn"])
    out = s[by + ["n_units"]].copy()
    out["agreement"], out["trivial_agreement"], out["kappa"] = po, trivial, kappa
    out["sensitivity"], out["specificity"] = sens, spec
    out["false_positive"], out["false_negative"] = fpr, fnr
    out["ctrl_differs"] = s["ctrl_differs"]
    return out


def _paired(df: pd.DataFrame, method: str, comparator: str) -> pd.DataFrame:
    keys = ["set", "study", "sex", "domain", "endpoint", "mode"]
    cols = ["agree", "sig_conc", "sig_vcg", "fp", "fn"]
    a = df[df["method"] == method][keys + cols]
    b = df[df["method"] == comparator][keys + cols]
    return a.merge(b, on=keys, suffixes=("_m", "_c"))


def _bootstrap_rows(p: pd.DataFrame, rng: np.random.Generator, resamples: int) -> dict:
    """Kappa difference (method - comparator) with a study-cluster bootstrap CI."""
    sums = p.groupby("study", observed=True).agg(
        n=("agree_m", "size"),
        **{f"{c}_{s}": (f"{c}_{s}", "sum") for s in ("m", "c") for c in ("agree", "sig_conc", "sig_vcg", "fp", "fn")})
    w = rng.multinomial(len(sums), np.full(len(sums), 1 / len(sums)), size=resamples)   # (B, studies)
    tot = {k: w @ sums[k].to_numpy(dtype=float) for k in sums.columns}
    point = {k: float(sums[k].sum()) for k in sums.columns}

    def kappas(t):
        km = _metrics_from_sums(t["n"], t["agree_m"], t["sig_conc_m"], t["sig_vcg_m"], t["fp_m"], t["fn_m"])
        kc = _metrics_from_sums(t["n"], t["agree_c"], t["sig_conc_c"], t["sig_vcg_c"], t["fp_c"], t["fn_c"])
        return km, kc

    (pm, km, sm, spm, fpm, fnm, _), (pc, kc, sc, spc, fpc, fnc, _) = kappas({k: np.array(v) for k, v in point.items()})
    (_, kmb, *_), (_, kcb, *_) = kappas(tot)
    diff = kmb - kcb
    lo, hi = np.nanpercentile(diff, [2.5, 97.5])
    return {"n_paired": int(point["n"]), "n_studies": len(sums),
            "kappa_method": float(km), "kappa_baseline": float(kc), "kappa_diff": float(km - kc),
            "ci_low": float(lo), "ci_high": float(hi),
            "agreement_method": float(pm), "agreement_baseline": float(pc),
            "sensitivity_method": float(sm), "sensitivity_baseline": float(sc),
            "specificity_method": float(spm), "specificity_baseline": float(spc),
            "false_positive_method": float(fpm), "false_positive_baseline": float(fpc)}


def verdict(diff: float, lo: float, hi: float) -> str:
    if diff >= C.HELPED_MIN_KAPPA_GAIN and lo > 0:
        return "helped"
    if diff < 0 and hi < 0:
        return "hurt"
    return "no clear difference"


def compare(df: pd.DataFrame, resamples: int = C.BOOTSTRAP_RESAMPLES, seed: int = C.SEED) -> dict[str, pd.DataFrame]:
    """Metric tables per method and paired comparisons with the comparator."""
    rng = np.random.default_rng(seed)
    tables = {
        "methods_by_domain": metrics(df, ["set", "method", "mode", "domain"]),
        "methods_overall": metrics(df[df["method"] != "baseline_trend"], ["set", "method", "mode"]),
    }
    tables["methods_vs_baseline"] = _comparisons(df, list(COMPARATOR.items()), rng, resamples)
    if df["method"].isin(["half_concurrent", "half_concurrent_trend"]).any():
        posthoc = _comparisons(df[df["mode"] == "partial"], POSTHOC_PAIRS, rng, resamples)
        posthoc.insert(0, "analysis", "post-hoc")
        tables["posthoc_half_concurrent"] = posthoc
    return tables


def _comparisons(df: pd.DataFrame, pairs: list[tuple[str, str]], rng: np.random.Generator,
                 resamples: int) -> pd.DataFrame:
    rows = []
    for method, comparator in pairs:
        p_all = _paired(df, method, comparator)
        for (set_label, mode), p in p_all.groupby(["set", "mode"], observed=True):
            parts = [(dom, q) for dom, q in p.groupby("domain", observed=True)] + [("all", p)]
            for dom, q in parts:
                if len(q) < 20 or q["study"].nunique() < 5:
                    continue
                r = _bootstrap_rows(q, rng, resamples)
                rows.append({"set": set_label, "method": method, "mode": mode, "domain": dom,
                             "comparator": comparator, **r,
                             "verdict": verdict(r["kappa_diff"], r["ci_low"], r["ci_high"])})
    return pd.DataFrame(rows)
