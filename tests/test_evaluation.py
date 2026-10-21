"""End-to-end checks of the evaluation on small synthetic data sets."""
import numpy as np
import pandas as pd
import pytest

from vcg import evaluation as E
from vcg import matching as M


def _unit(effect: float, pool_shift: float = 0.0, n_hist_studies: int = 6, seed: int = 0) -> E.Unit:
    rng = np.random.default_rng(seed)
    target = M.Target(study="T", sex="Male", species="Rat", strain="F344", design="subchronic",
                      route="Gavage", vehicle="Corn oil", laboratory="Lab A", year=2000, dtxsid="DT")
    rows = []
    for dose in [0, 1, 2, 3]:
        for i in range(10):
            rows.append(dict(study="T", sex="Male", dose=float(dose), role="Vehicle control" if dose == 0 else "Treated",
                             is_control=dose == 0, animal=f"T{dose}{i}", domain="organ_weight",
                             endpoint="Liver weight", value=rng.normal(10 + effect * dose), unit="g", days=91.0))
    cont = pd.DataFrame(rows)
    hist = []
    for s in range(n_hist_studies):
        for i in range(10):
            hist.append(dict(study=f"H{s}", sex="Male", dose=0.0, role="Vehicle control", is_control=True,
                             animal=f"H{s}{i}", domain="organ_weight", endpoint="Liver weight",
                             value=rng.normal(10 + pool_shift), unit="g", days=91.0, species="Rat", strain="F344",
                             design="subchronic", route="Gavage", vehicle="Corn oil", laboratory="Lab A",
                             year=2000 + s % 3, dtxsid=f"DH{s}"))
    hist = pd.DataFrame(hist)
    empty = pd.DataFrame()
    return E.Unit(target, cont, empty, empty, hist, empty, empty, n_draws=5, seed=1)


def test_records_have_expected_shape_and_modes():
    rows = E.evaluate_continuous(_unit(effect=0.0))
    df = pd.DataFrame(rows, columns=E.RECORD_COLUMNS)
    assert set(df["mode"]) == {"vcg", "partial"}
    assert set(df["dose_rank"]) == {1, 2, 3}
    assert df["draw"].max() == 4
    assert (df["n_control"] == 10).all()
    # concurrent p-values are identical across configurations for the same dose group
    assert df.groupby("dose_rank")["p_conc"].nunique().eq(1).all()


def test_strong_effect_is_found_with_both_controls():
    df = pd.DataFrame(E.evaluate_continuous(_unit(effect=3.0)), columns=E.RECORD_COLUMNS)
    assert df["sig_conc"].all()
    assert df["sig_vcg"].mean() > 0.95


def test_shifted_pool_creates_false_positives():
    df = pd.DataFrame(E.evaluate_continuous(_unit(effect=0.0, pool_shift=1.0)), columns=E.RECORD_COLUMNS)
    ep = E.endpoint_level(df)
    full = ep[ep["mode"] == "vcg"]
    assert not full["sig_conc"].any()
    assert full["sig_vcg"].mean() > 0.5           # the mismatched pool makes doses look different
    partial = ep[ep["mode"] == "partial"]
    assert partial["sig_vcg"].mean() < full["sig_vcg"].mean()  # keeping half the concurrent animals dampens it


def test_rates_definitions():
    df = pd.DataFrame(E.evaluate_continuous(_unit(effect=0.0, pool_shift=2.5)), columns=E.RECORD_COLUMNS)
    tables = E.summarise(df)
    t = tables["agreement_by_level"]
    row = t[(t["mode"] == "vcg") & (t["level"] == "L0 design")].iloc[0]
    assert row["agreement"] + row["false_positive"] + row["false_negative"] == pytest.approx(1.0)
    assert row["false_negative"] == 0.0                   # nothing was significant with the concurrent control
    assert row["fp_given_conc_nonsig"] == pytest.approx(row["false_positive"])
    assert "coverage_by_level" in tables and tables["coverage_by_level"]["coverage"].between(0, 1).all()


def test_small_pool_is_skipped_for_full_vcg_only():
    unit = _unit(effect=0.0, n_hist_studies=1)       # 10 animals: < 2 x n_control, but enough for 5 virtual
    df = pd.DataFrame(E.evaluate_continuous(unit), columns=E.RECORD_COLUMNS)
    assert set(df["mode"]) == {"partial"}
    assert E.evaluate_continuous(_unit(effect=0.0, n_hist_studies=0)) == []
