import numpy as np
import pytest
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


def test_bw_weights_follow_the_target_distribution():
    target = np.array([150.0, 152, 148, 151, 149, 150, 153, 147])
    pool = np.array([150.0, 153, 140, 170, np.nan])
    w = M.bw_weights(pool, target, window_sd=2.0)
    sd = target.std(ddof=1)
    assert w[0] == 1.0                                   # at the target mean
    assert w[1] == pytest.approx(np.exp(-0.5 * ((153 - 150) / sd) ** 2))
    assert w[2] == 0.0 and w[3] == 0.0                   # outside +/- 2 SD
    assert w[4] == 0.0                                   # no pre-dosing weight
    assert M.bw_weights(pool, target[:3]) is None        # too few target weights


def test_weighted_draw_prefers_matching_animals():
    rng = np.random.default_rng(0)
    w = np.r_[np.full(50, 1.0), np.full(50, 0.05)]
    picks = np.concatenate([M.draw_weighted(w, 10, rng) for _ in range(200)])
    assert (picks < 50).mean() > 0.85
    assert len(set(M.draw_weighted(w, 10, rng))) == 10   # without replacement


def test_whole_groups_are_kept_intact():
    rng = np.random.default_rng(1)
    groups = np.repeat(np.array(["A", "B", "C", "D"]), [10, 8, 12, 10])
    for _ in range(50):
        idx = M.draw_whole_groups(groups, 10, rng)
        chosen = set(groups[idx])
        assert all((groups == g).sum() == (groups[idx] == g).sum() for g in chosen)
        assert idx.size >= 10
        # groups stop being added once n_needed is reached: without the last
        # (hence without the largest) group the draw would have been too small
        assert idx.size - max((groups == g).sum() for g in chosen) < 10


def test_stable_key_is_deterministic():
    assert M.stable_key("002-01015-0016-0000-6", "Male") == M.stable_key("002-01015-0016-0000-6", "Male")
    assert M.stable_key("a", "Male") != M.stable_key("a", "Female")
