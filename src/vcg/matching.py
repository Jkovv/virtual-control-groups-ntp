"""Select historical control animals for a target study and draw virtual control groups.

A *pool* is the set of control animals from other studies that satisfy a
matching level (see ``config.MATCH_LEVELS``).  A *virtual control group* (VCG)
is a random sample from the pool of the same size as the concurrent control
group it stands in for.  A *partial* VCG keeps a random half of the concurrent
controls and fills the rest from the pool.
"""
from __future__ import annotations

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


def describe_pool(controls: pd.DataFrame, mask: np.ndarray) -> dict:
    sub = controls.loc[mask]
    return {"n_animals": int(mask.sum()), "n_studies": int(sub["study"].nunique()),
            "n_labs": int(sub["laboratory"].nunique())}
