import numpy as np
import pytest
from scipy import stats as sps

from vcg import stats as S


@pytest.mark.parametrize("seed", range(6))
def test_dunnett_matches_scipy(seed):
    rng = np.random.default_rng(seed)
    k = rng.integers(1, 5)
    control = rng.normal(size=rng.integers(5, 40))
    groups = [rng.normal(loc=rng.normal(0, 0.7), size=rng.integers(5, 40)) for _ in range(k)]
    mine = S.dunnett(groups, control)
    ref = sps.dunnett(*groups, control=control, random_state=0)
    assert np.allclose(mine.statistic, ref.statistic)
    assert np.allclose(mine.pvalue, ref.pvalue, atol=2e-3)


def test_dunnett_detects_large_shift():
    rng = np.random.default_rng(1)
    control = rng.normal(size=10)
    res = S.dunnett([rng.normal(size=10), rng.normal(loc=3, size=10)], control)
    assert res.pvalue[0] > 0.05 and res.pvalue[1] < 0.001
    assert res.statistic[1] > 0


def test_dunnett_degenerate_groups_give_nan():
    res = S.dunnett([np.array([1.0])], np.array([0.0, 1.0]))
    assert np.isnan(res.pvalue).all()


def test_jonckheere_trend():
    rng = np.random.default_rng(2)
    groups = [rng.normal(loc=i, size=12) for i in range(4)]
    z, p = S.jonckheere_terpstra(groups)
    assert z > 3 and p < 0.01
    _, p_null = S.jonckheere_terpstra([rng.normal(size=12) for _ in range(4)])
    assert p_null > 0.05


def test_fisher_matches_scipy():
    tables = [(3, 10, 0, 10), (0, 10, 0, 10), (10, 10, 2, 10), (5, 50, 5, 50), (7, 12, 1, 9)]
    for x1, n1, x0, n0 in tables:
        mine = S.fisher_exact_2x2(x1, n1, x0, n0)[0]
        ref = sps.fisher_exact([[x1, n1 - x1], [x0, n0 - x0]])[1]
        assert mine == pytest.approx(ref, abs=1e-9)


def test_poly3_weights():
    t = np.array([730.0, 365.0, 365.0])
    affected = np.array([False, True, False])
    w = S.poly_k_weights(t, affected, 730.0)
    assert w[0] == 1.0            # survived to terminal sacrifice
    assert w[1] == 1.0            # tumour-bearing animals always count fully
    assert w[2] == pytest.approx(0.5**3)


def test_poly3_direction_and_null():
    t = np.full(50, 730.0)
    control = np.zeros(50, bool); control[:2] = True
    treated = np.zeros(50, bool); treated[:15] = True
    w = S.poly_k_weights(t, control, 730.0)
    z, p = S.poly3_test(treated, S.poly_k_weights(t, treated, 730.0), control, w)
    assert z > 2 and p < 0.01
    z0, p0 = S.poly3_test(control, w, control, w)
    assert z0 == 0 and p0 == 0.5


def test_dunnett_many_matches_single_calls():
    rng = np.random.default_rng(3)
    groups = [rng.normal(size=10) for _ in range(3)]
    controls = rng.normal(size=(6, 10))
    many = S.dunnett_many(groups, controls)
    assert many.pvalue.shape == (6, 3)
    for i in range(6):
        one = S.dunnett(groups, controls[i])
        assert np.allclose(one.statistic, many.statistic[i])
        assert np.allclose(one.pvalue, many.pvalue[i], atol=1e-9)


def test_interpolated_sf_matches_exact_quadrature():
    lam = np.sqrt(np.array([0.5, 0.4, 0.6]))
    t = np.linspace(0.1, 5.5, 40)
    assert np.allclose(S.dunnett_sf(t, lam, 36), S.dunnett_sf(t, lam, 36, exact=True), atol=2e-4)


def test_welch_many_matches_scipy():
    rng = np.random.default_rng(4)
    a = rng.normal(size=12)
    rows = rng.normal(loc=0.5, size=(5, 9))
    ref = [sps.ttest_ind(a, r, equal_var=False).pvalue for r in rows]
    assert np.allclose(S.welch_p_many(a, rows), ref)
