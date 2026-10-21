import numpy as np
import pandas as pd

from vcg import config as C
from vcg import matching as M


def _controls():
    rows = []
    for study, dtx, route, veh, lab, year in [
        ("S1", "D1", "Gavage", "Corn oil", "Lab A", 2000),
        ("S2", "D2", "Gavage", "Corn oil", "Lab A", 2001),
        ("S3", "D3", "Gavage", "Water", "Lab A", 2002),
        ("S4", "D4", "Gavage", "Corn oil", "Lab B", 2000),
        ("S5", "D5", "Dosed feed", "Corn oil", "Lab A", 2000),
        ("S6", "D1", "Gavage", "Corn oil", "Lab A", 2000),   # same chemical as S1
        ("S7", "D7", "Gavage", "Corn oil", "Lab A", 2010),   # far away in time
    ]:
        for i in range(10):
            rows.append(dict(study=study, animal=f"{study}-{i}", species="Rat", strain="F344", sex="Male",
                             design="subchronic", route=route, vehicle=veh, laboratory=lab, year=year,
                             dtxsid=dtx, value=float(i)))
    return pd.DataFrame(rows)


TARGET = M.Target(study="S1", sex="Male", species="Rat", strain="F344", design="subchronic",
                  route="Gavage", vehicle="Corn oil", laboratory="Lab A", year=2000, dtxsid="D1")


def test_levels_are_nested_and_exclude_same_study_and_chemical():
    ctrl = _controls()
    sizes = []
    for lvl, crit in C.MATCH_LEVELS.items():
        m = M.pool_mask(ctrl, TARGET, crit)
        studies = set(ctrl.loc[m, "study"])
        assert "S1" not in studies and "S6" not in studies, lvl
        sizes.append(m.sum())
    assert sizes == sorted(sizes, reverse=True)
    assert set(ctrl.loc[M.pool_mask(ctrl, TARGET, ()), "study"]) == {"S2", "S3", "S4", "S5", "S7"}
    assert set(ctrl.loc[M.pool_mask(ctrl, TARGET, ("route",)), "study"]) == {"S2", "S3", "S4", "S7"}
    assert set(ctrl.loc[M.pool_mask(ctrl, TARGET, ("route", "vehicle")), "study"]) == {"S2", "S4", "S7"}
    assert set(ctrl.loc[M.pool_mask(ctrl, TARGET, ("route", "vehicle", "laboratory")), "study"]) == {"S2", "S7"}
    assert set(ctrl.loc[M.pool_mask(ctrl, TARGET, ("route", "vehicle", "laboratory", "period")), "study"]) == {"S2"}


def test_missing_attribute_never_matches():
    ctrl = _controls()
    ctrl.loc[ctrl["study"] == "S2", "laboratory"] = None
    m = M.pool_mask(ctrl, TARGET, ("route", "vehicle", "laboratory"))
    assert set(ctrl.loc[m, "study"]) == {"S7"}


def test_year_gap_band():
    ctrl = _controls()
    m = M.year_gap_mask(ctrl, TARGET, (6, 20))
    assert set(ctrl.loc[m, "study"]) == {"S7"}


def test_draws_have_requested_size_and_partial_keeps_half():
    rng = np.random.default_rng(0)
    pool = np.arange(100, 200, dtype=float)
    concurrent = np.arange(10, dtype=float)
    vcg = M.draw_vcg(pool, 10, rng)
    assert vcg.size == 10 and len(set(vcg)) == 10 and vcg.min() >= 100
    partial = M.draw_partial(concurrent, pool, 5, rng)
    assert partial.size == 10
    assert ((partial < 10).sum(), (partial >= 100).sum()) == (5, 5)
