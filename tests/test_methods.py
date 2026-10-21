"""Method comparison on synthetic studies where the mechanism is known."""
import numpy as np
import pandas as pd
import pytest

from vcg import evaluation as E
from vcg import matching as M
from vcg import methods as MM

TARGET = M.Target(study="T", sex="Male", species="Rat", strain="F344", design="subchronic",
                  route="Gavage", vehicle="Corn oil", laboratory="Lab A", year=1995, dtxsid="DT")


def _unit(seed=0, tau=0.0, bw_shift=0.0, slope=0.0, effect=0.0, n_hist=12, n_draws=10):
    """A 3-month study (control + 3 dose groups of 10) and ``n_hist`` historical
    control groups whose means vary between studies with SD ``tau``.  The
    endpoint depends on pre-dosing body weight with ``slope``; the target's
    animals start ``bw_shift`` grams heavier than the historical ones."""
    rng = np.random.default_rng(seed)
    rows = []
    target_offset = rng.normal(0, tau)
    for dose in range(4):
        for i in range(10):
            bw0 = rng.normal(150 + bw_shift, 8)
            rows.append(dict(study="T", sex="Male", dose=float(dose), role="Vehicle control" if dose == 0 else "Treated",
                             is_control=dose == 0, animal=f"T{dose}-{i}", domain="organ_weight",
                             endpoint="Liver weight", value=10 + target_offset + slope * (bw0 - 150) + effect * dose
                             + rng.normal(0, 0.5), unit="g", days=91.0, bw0=bw0))
    hist = []
    for s in range(n_hist):
        offset = rng.normal(0, tau)
        for i in range(10):
            bw0 = rng.normal(150, 8)
            hist.append(dict(study=f"H{s}", sex="Male", dose=0.0, role="Vehicle control", is_control=True,
                             animal=f"H{s}-{i}", domain="organ_weight", endpoint="Liver weight",
                             value=10 + offset + slope * (bw0 - 150) + rng.normal(0, 0.5), unit="g", days=91.0,
                             bw0=bw0, species="Rat", strain="F344", design="subchronic", route="Gavage",
                             vehicle="Corn oil", laboratory="Lab A", year=1990, dtxsid=f"D{s}"))
    empty = pd.DataFrame()
    return E.Unit(TARGET, pd.DataFrame(rows), empty, empty, pd.DataFrame(hist), empty, empty, n_draws=n_draws, seed=1)


def _frame(rows):
    return pd.DataFrame(rows, columns=MM.ENDPOINT_COLUMNS)


def test_all_continuous_methods_and_modes_are_evaluated():
    df = _frame(MM.evaluate_continuous(_unit(), "holdout"))
    got = set(zip(df["method"], df["mode"]))
    for method in ("baseline", "bw_matched", "bw_ancova", "whole_group", "hierarchical"):
        assert (method, "vcg") in got and (method, "partial") in got
    assert df.loc[(df["method"] == "hierarchical") & (df["mode"] == "vcg"), "n_draws"].item() == 1
    assert (df["agree"] + df["fp"] + df["fn"]).round(9).eq(1).all()
    assert df["set"].eq("holdout").all()


def _vcg_fp_rate(method: str, reps: int = 40, **kw) -> float:
    rates = []
    for seed in range(reps):
        df = _frame(MM.evaluate_continuous(_unit(seed=seed, **kw), "holdout"))
        row = df[(df["method"] == method) & (df["mode"] == "vcg")]
        if len(row) and row["sig_conc"].item() == 0:
            rates.append(row["sig_vcg"].item())
    return float(np.mean(rates))


def test_between_study_variance_creates_false_positives_that_the_hierarchical_model_absorbs():
    # no treatment effect, between-study SD twice the within-study SD
    fp_baseline = _vcg_fp_rate("baseline", tau=1.0)
    fp_hier = _vcg_fp_rate("hierarchical", tau=1.0)
    assert fp_baseline > 0.4
    assert fp_hier < fp_baseline / 2


def test_without_between_study_variance_baseline_is_well_calibrated():
    assert _vcg_fp_rate("baseline", tau=0.0) < 0.2


def test_ancova_absorbs_a_pre_dosing_weight_mismatch():
    # target animals start 12 g heavier; the endpoint scales with starting weight
    kw = dict(tau=0.0, bw_shift=12.0, slope=0.08)
    fp_baseline = _vcg_fp_rate("baseline", **kw)
    fp_ancova = _vcg_fp_rate("bw_ancova", **kw)
    assert fp_baseline > 0.5
    assert fp_ancova < 0.25


def test_missing_pre_dosing_weights_switch_off_the_bw_methods():
    unit = _unit()
    unit.continuous["bw0"] = np.nan
    df = _frame(MM.evaluate_continuous(unit, "development"))
    assert not df["method"].isin(["bw_matched", "bw_ancova"]).any()
    assert df["method"].isin(["baseline", "whole_group", "hierarchical"]).any()


def test_too_few_historical_groups_switch_off_group_methods():
    df = _frame(MM.evaluate_continuous(_unit(n_hist=2), "holdout"))
    assert not df["method"].isin(["whole_group", "hierarchical"]).any()


def test_set_assignment_uses_the_cutoff_year():
    assert MM.set_of(1989.0) == "development"
    assert MM.set_of(1990.0) == "holdout"
    assert MM.set_of(float("nan")) is None


def test_verdict_rule():
    assert MM.verdict(0.08, 0.02, 0.14) == "helped"
    assert MM.verdict(0.03, 0.01, 0.05) == "no clear difference"    # gain below 0.05
    assert MM.verdict(0.08, -0.01, 0.17) == "no clear difference"   # interval includes 0
    assert MM.verdict(-0.04, -0.07, -0.01) == "hurt"


def test_compare_pairs_methods_on_the_same_endpoints():
    frames = [_frame(MM.evaluate_continuous(_unit(seed=s, tau=1.0), "holdout")) for s in range(25)]
    df = pd.concat(frames, ignore_index=True)
    df["study"] = np.repeat([f"T{s}" for s in range(25)], [len(f) for f in frames])
    out = MM.compare(df, resamples=200)
    vs = out["methods_vs_baseline"]
    row = vs[(vs["method"] == "hierarchical") & (vs["mode"] == "vcg") & (vs["domain"] == "all")].iloc[0]
    assert row["n_paired"] == 25 and row["n_studies"] == 25
    assert row["ci_low"] <= row["kappa_diff"] <= row["ci_high"]
    assert row["false_positive_method"] < row["false_positive_baseline"]
    m = out["methods_by_domain"]
    assert m["agreement"].between(0, 1).all()
    assert (m["false_positive"] + m["false_negative"] + m["agreement"]).round(9).eq(1).all()


def _chronic_unit(seed=0, hist_rate=0.04, treated_rates=(0.15, 0.30, 0.45)):
    rng = np.random.default_rng(seed)
    target = M.Target(study="T", sex="Male", species="Mouse", strain="B6C3F1/N", design="chronic",
                      route="Gavage", vehicle="Corn oil", laboratory="Lab A", year=1995, dtxsid="DT")
    roster, lesions = [], []
    for dose in range(4):
        rate = hist_rate if dose == 0 else treated_rates[dose - 1]
        for i in range(50):
            a = f"T{dose}-{i}"
            roster.append(dict(study="T", sex="Male", dose=float(dose), is_control=dose == 0, animal=a,
                               days=730.0, t_max=730.0, bw0=np.nan))
            if rng.random() < rate:
                lesions.append(dict(study="T", sex="Male", animal=a, tissue="Liver",
                                    morphology="Hepatocellular adenoma", neoplastic="Neoplastic"))
    hist, hles = [], []
    for s in range(8):
        for i in range(50):
            a = f"H{s}-{i}"
            hist.append(dict(study=f"H{s}", sex="Male", dose=0.0, is_control=True, animal=a, days=730.0,
                             t_max=730.0, bw0=np.nan, laboratory="Lab A", year=1990, dtxsid=f"D{s}"))
            if rng.random() < hist_rate:
                hles.append(dict(study=f"H{s}", sex="Male", animal=a, tissue="Liver",
                                 morphology="Hepatocellular adenoma", neoplastic="Neoplastic"))
    return E.Unit(target, pd.DataFrame(), pd.DataFrame(roster), pd.DataFrame(lesions), pd.DataFrame(),
                  pd.DataFrame(hist), pd.DataFrame(hles), n_draws=10, seed=1)


def test_neoplasm_methods_include_trend_tests():
    df = _frame(MM.evaluate_incidence(_chronic_unit(), "holdout"))
    got = set(zip(df["method"], df["mode"]))
    for method in ("baseline", "whole_group", "baseline_trend", "tarone"):
        assert (method, "vcg") in got and (method, "partial") in got
    assert ("bw_matched", "vcg") not in got                  # no pre-dosing weights
    tarone = df[(df["method"] == "tarone") & (df["mode"] == "vcg")].iloc[0]
    assert tarone["n_draws"] == 1 and tarone["domain"] == "neoplasm"
    # a strong dose-related increase is found by the reference and by Tarone's test
    assert tarone["sig_conc"] == 1 and tarone["sig_vcg"] == 1


@pytest.mark.parametrize("sampling", ["individual", "whole_group"])
def test_implied_false_positive_rate_matches_simulation(sampling):
    from scipy import stats
    assert MM.implied_false_positive_rate(0.0, sampling=sampling) == pytest.approx(0.05, abs=0.002)
    rng = np.random.default_rng(0)
    sims, n, tau, icc = 20000, 10, 0.5, 0.2               # tau^2 = 0.25, sigma^2 = 1
    treated = rng.normal(rng.normal(0, tau, (sims, 1)), 1.0, (sims, n))
    if sampling == "individual":
        control = rng.normal(rng.normal(0, tau, (sims, n)), 1.0)        # one study per animal
    else:
        control = rng.normal(rng.normal(0, tau, (sims, 1)), 1.0, (sims, n))
    p = stats.ttest_ind(treated, control, axis=1).pvalue
    assert (p < 0.05).mean() == pytest.approx(MM.implied_false_positive_rate(icc, n, sampling), abs=0.02)


def test_whole_groups_are_predicted_to_be_worse_than_individual_sampling():
    icc = np.array([0.05, 0.1, 0.3])
    assert (MM.implied_false_positive_rate(icc, sampling="whole_group")
            > MM.implied_false_positive_rate(icc, sampling="individual")).all()


def test_half_concurrent_control_uses_only_concurrent_animals():
    df = _frame(MM.evaluate_half_concurrent(_unit(tau=1.0), "holdout"))
    assert set(df["method"]) == {"half_concurrent"} and set(df["mode"]) == {"partial"}
    assert (df["pool_animals"] == 0).all()
    inc = _frame(MM.evaluate_half_concurrent(_chronic_unit(), "holdout"))
    assert set(inc["method"]) == {"half_concurrent", "half_concurrent_trend"}


def test_split_half_reference_is_disjoint_from_the_tested_half():
    """The two halves must never share an animal: with a control group whose
    values are all distinct, the recorded halves cannot overlap."""
    rng = np.random.default_rng(0)
    unit = _unit(seed=3, n_draws=30)
    ctrl = unit.continuous.loc[unit.continuous["is_control"], "value"]
    unit.continuous.loc[ctrl.index, "value"] = np.arange(len(ctrl)) * 100.0   # unique, far apart
    df = _frame(MM.evaluate_reliability(unit, "holdout"))
    assert "split_half" in set(df["method"])
    assert set(df["mode"]) == {"disjoint"}
    row = df[df["method"] == "split_half"].iloc[0]
    assert row["n_draws"] == 30 and row["pool_animals"] == 0
    assert (row["agree"] + row["fp"] + row["fn"]) == pytest.approx(1.0)


def test_split_half_is_a_ceiling_not_a_perfect_score():
    """Two halves of the same control group disagree sometimes, and they agree
    better than a control borrowed from other studies does."""
    units = [_unit(seed=s, tau=1.0, effect=0.35) for s in range(30)]
    rel = pd.concat([_frame(MM.evaluate_reliability(u, "holdout")) for u in units], ignore_index=True)
    pre = pd.concat([_frame(MM.evaluate_continuous(u, "holdout")) for u in units], ignore_index=True)
    ceiling = rel[rel["method"] == "split_half"]["agree"].mean()
    borrowed = pre[(pre["method"] == "baseline") & (pre["mode"] == "vcg")]["agree"].mean()
    assert 0.5 < ceiling < 1.0
    assert ceiling > borrowed


def test_disjoint_mode_methods_share_the_reference():
    df = _frame(MM.evaluate_reliability(_unit(tau=1.0), "holdout"))
    got = set(df["method"])
    assert {"split_half", "baseline", "hierarchical"} <= got
    # every method in a unit/endpoint is scored against the same reference draws
    assert df.groupby("endpoint")["sig_conc"].nunique().eq(1).all()


def test_reliability_covers_incidence_and_trend_tests():
    df = _frame(MM.evaluate_reliability(_chronic_unit(), "holdout"))
    got = set(df["method"])
    assert {"split_half", "split_half_trend", "baseline", "baseline_trend", "tarone"} <= got
    assert set(df["domain"]) == {"neoplasm"}


def test_reliability_needs_at_least_six_control_animals():
    unit = _unit()
    keep = unit.continuous[~unit.continuous["is_control"]]
    small = unit.continuous[unit.continuous["is_control"]].head(5)
    unit.continuous = pd.concat([small, keep], ignore_index=True)
    assert MM.evaluate_reliability(unit, "holdout") == []


def test_half_size_vcg_matches_the_ceiling_group_size():
    """``vcg_half`` must hold exactly as many animals as the reference half, so
    that comparing it with ``split_half`` is like for like."""
    df = _frame(MM.evaluate_reliability(_unit(tau=1.0), "holdout"))
    assert {"vcg_half", "split_half"} <= set(df["method"])
    inc = _frame(MM.evaluate_reliability(_chronic_unit(), "holdout"))
    assert {"vcg_half", "vcg_half_trend"} <= set(inc["method"])
    # same reference draws for both, so the comparison is paired
    both = df[df["method"].isin(["vcg_half", "split_half"])]
    assert both.groupby("endpoint")["sig_conc"].nunique().eq(1).all()


def test_half_size_vcg_agrees_less_than_a_real_half_when_studies_differ():
    units = [_unit(seed=s, tau=1.0, effect=0.35) for s in range(30)]
    df = pd.concat([_frame(MM.evaluate_reliability(u, "holdout")) for u in units], ignore_index=True)
    ceiling = df[df["method"] == "split_half"]["agree"].mean()
    virtual = df[df["method"] == "vcg_half"]["agree"].mean()
    assert virtual < ceiling
    # with no between-study variance the two are close
    flat = pd.concat([_frame(MM.evaluate_reliability(_unit(seed=s, tau=0.0, effect=0.35), "holdout"))
                      for s in range(30)], ignore_index=True)
    assert flat[flat["method"] == "vcg_half"]["agree"].mean() == pytest.approx(
        flat[flat["method"] == "split_half"]["agree"].mean(), abs=0.08)


def test_half_size_vcg_is_size_matched_to_the_reference_half():
    """``vcg_half`` holds exactly as many animals as the reference half, so that
    comparing it with ``split_half`` is like for like."""
    df = _frame(MM.evaluate_reliability(_unit(tau=1.0), "holdout"))
    assert {"vcg_half", "split_half"} <= set(df["method"])
    inc = _frame(MM.evaluate_reliability(_chronic_unit(), "holdout"))
    assert {"vcg_half", "vcg_half_trend"} <= set(inc["method"])
    both = df[df["method"].isin(["vcg_half", "split_half"])]
    assert both.groupby("endpoint")["sig_conc"].nunique().eq(1).all()   # same reference draws


def test_half_size_vcg_trails_a_real_half_only_when_studies_differ():
    units = [_unit(seed=s, tau=1.0, effect=0.35) for s in range(30)]
    df = pd.concat([_frame(MM.evaluate_reliability(u, "holdout")) for u in units], ignore_index=True)
    assert df[df["method"] == "vcg_half"]["agree"].mean() < df[df["method"] == "split_half"]["agree"].mean()
    flat = pd.concat([_frame(MM.evaluate_reliability(_unit(seed=s, tau=0.0, effect=0.35), "holdout"))
                      for s in range(30)], ignore_index=True)
    assert flat[flat["method"] == "vcg_half"]["agree"].mean() == pytest.approx(
        flat[flat["method"] == "split_half"]["agree"].mean(), abs=0.08)
