"""Select historical control animals for a target study and draw virtual control groups.

A *pool* is the set of control animals from other studies that satisfy a
matching level (see ``config.MATCH_LEVELS``).  A *virtual control group* (VCG)
is a random sample from the pool of the same size as the concurrent control
group it stands in for.  A *partial* VCG keeps a random half of the concurrent
controls and fills the rest from the pool.
"""
from __future__ import annotations

import zlib
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import config as C

MATCH_COLUMNS = ("species", "strain", "sex", "design", "route", "vehicle", "laboratory", "year", "dtxsid")


@dataclass(frozen=True)
class Target:
    """Attributes of the study/sex whose concurrent controls are replaced."""
    study: str
    sex: str
    species: str
    strain: str
    design: str
    route: object
    vehicle: object
    laboratory: object
    year: float
    dtxsid: object


def stable_key(*parts: object) -> int:
    """Deterministic integer from strings, for seeding.  Python's ``hash`` is
    salted per process, so it cannot be used for reproducible random streams."""
    return zlib.crc32("|".join(map(str, parts)).encode())


def target_from_row(row: pd.Series, sex: str) -> Target:
    return Target(study=row["study"], sex=sex, species=row["species"], strain=row["strain"],
                  design=row["design"], route=row["route"], vehicle=row["vehicle"],
                  laboratory=row["laboratory"], year=row["year"], dtxsid=row["dtxsid"])


def _eq(series: pd.Series, value: object) -> np.ndarray:
    """Equality that treats missing values as never matching."""
    if value is None or (isinstance(value, float) and np.isnan(value)) or value is pd.NA:
        return np.zeros(len(series), dtype=bool)
    return (series == value).to_numpy(dtype=bool, na_value=False)


def pool_mask(controls: pd.DataFrame, target: Target, criteria: tuple[str, ...],
              period_years: int = C.PERIOD_YEARS) -> np.ndarray:
    """Boolean mask over ``controls`` (historical control animals with the
    MATCH_COLUMNS) selecting the animals eligible for the target at a level."""
    m = _eq(controls["species"], target.species) & _eq(controls["strain"], target.strain)
    m &= _eq(controls["sex"], target.sex) & _eq(controls["design"], target.design)
    # never borrow from the same study, nor from any study of the same test article
    m &= ~_eq(controls["study"], target.study)
    m &= ~_eq(controls["dtxsid"], target.dtxsid)
    for crit in criteria:
        if crit == "period":
            gap = (controls["year"] - target.year).abs().to_numpy(dtype=float, na_value=np.nan)
            m &= np.nan_to_num(gap, nan=np.inf) <= period_years
        else:
            m &= _eq(controls[crit], getattr(target, crit))
    return m


def year_gap_mask(controls: pd.DataFrame, target: Target, band: tuple[int, int]) -> np.ndarray:
    gap = (controls["year"] - target.year).abs().to_numpy(dtype=float, na_value=np.nan)
    lo, hi = band
    return (gap >= lo) & (gap <= hi)


def draw_vcg(pool_values: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    """Sample ``n`` historical animals without replacement."""
    return rng.choice(pool_values, size=n, replace=False)


def draw_partial(concurrent: np.ndarray, pool_values: np.ndarray, n_virtual: int,
                 rng: np.random.Generator) -> np.ndarray:
    """Keep a random subset of the concurrent controls and top up from the pool."""
    keep = rng.choice(concurrent, size=len(concurrent) - n_virtual, replace=False)
    return np.concatenate([keep, rng.choice(pool_values, size=n_virtual, replace=False)])


def bw_weights(pool_bw: np.ndarray, target_bw: np.ndarray, window_sd: float = C.BW_WINDOW_SD,
               min_target: int = C.BW_MIN_TARGET_ANIMALS) -> np.ndarray | None:
    """Sampling weights that match the target's pre-dosing body weight distribution.

    Following the weighted-sampling approach of Gurjanov et al. (ALTEX 2024), a
    pool animal's weight is the normal density of its pre-dosing body weight
    under the target study's pre-dosing distribution (mean and SD of all target
    animals, before any treatment).  Animals outside +/- ``window_sd`` target SDs,
    or without a pre-dosing weight, get weight 0.  Returns None when the target
    has fewer than ``min_target`` pre-dosing weights.
    """
    tb = np.asarray(target_bw, dtype=float)
    tb = tb[np.isfinite(tb)]
    if tb.size < min_target:
        return None
    mu, sd = tb.mean(), tb.std(ddof=1)
    if not sd > 0:
        return None
    z = (np.asarray(pool_bw, dtype=float) - mu) / sd
    w = np.exp(-0.5 * z**2)
    return np.where(np.isfinite(z) & (np.abs(z) <= window_sd), w, 0.0)


def draw_weighted(weights: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    """Indices of ``n`` pool animals drawn without replacement with probability
    proportional to ``weights``."""
    return rng.choice(weights.size, size=n, replace=False, p=weights / weights.sum())


def draw_whole_groups(group_ids: np.ndarray, n_needed: int, rng: np.random.Generator) -> np.ndarray:
    """Indices of all animals of randomly ordered whole historical control groups,
    adding groups until at least ``n_needed`` animals are included.

    Keeping groups intact preserves the within-study structure (shared housing,
    diet batch, necropsy day, analyser run) that sampling individual animals
    from many studies averages away.
    """
    labels, inverse, counts = np.unique(group_ids, return_inverse=True, return_counts=True)
    order = rng.permutation(labels.size)
    take = np.searchsorted(np.cumsum(counts[order]), n_needed) + 1
    chosen = np.zeros(labels.size, dtype=bool)
    chosen[order[:take]] = True
    return np.flatnonzero(chosen[inverse])


def describe_pool(controls: pd.DataFrame, mask: np.ndarray) -> dict:
    sub = controls.loc[mask]
    return {"n_animals": int(mask.sum()), "n_studies": int(sub["study"].nunique()),
            "n_labs": int(sub["laboratory"].nunique())}
