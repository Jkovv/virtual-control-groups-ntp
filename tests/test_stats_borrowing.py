"""Reference-value tests for the statistics added for the method comparison:
ANCOVA, DerSimonian-Laird, (robust) MAP priors, Cochran-Armitage, Tarone and
the beta-binomial fit.  Each is checked against an independent source: a
published example, statsmodels, scipy, or numerical integration."""
import numpy as np
import pytest
import statsmodels.api as sm
from scipy import integrate
from scipy import stats as sps
from statsmodels.stats.meta_analysis import combine_effects

from vcg import stats as S


# --- Cochran-Armitage and Tarone ---------------------------------------------

AGRESTI_X = np.array([48, 38, 5, 1, 1])                       # malformations
AGRESTI_N = np.array([17114, 14502, 793, 127, 38])            # infants
AGRESTI_SCORES = np.array([0.0, 0.5, 1.5, 4.0, 7.0])          # drinks per day


def test_cochran_armitage_matches_agresti_example():
    # Agresti (2002), Categorical Data Analysis, 2nd ed., section 3.4.4, Table 3.7:
    # linear-by-linear statistic M^2 = 6.57, two-sided P = 0.01.  M^2 = (N - 1)/N * z_CA^2.
    z, p = S.cochran_armitage(AGRESTI_X, AGRESTI_N, AGRESTI_SCORES)
    n_total = AGRESTI_N.sum()
    assert z**2 * (n_total - 1) / n_total == pytest.approx(6.57, abs=0.005)
    assert 2 * p == pytest.approx(0.0104, abs=0.0005)


def test_tarone_without_history_is_cochran_armitage():
    x, n, d = np.array([1, 2, 4, 9]), np.array([50, 50, 50, 50]), np.arange(4.0)
    assert S.tarone_trend(x, n, d, 0.0, 0.0) == S.cochran_armitage(x, n, d)


def test_tarone_equals_cochran_armitage_with_historical_pseudo_group():
    x, n, d = np.array([2, 3, 6, 9]), np.array([50, 50, 50, 50]), np.arange(4.0)
    a, b = 3.0, 97.0
    with_prior = S.tarone_trend(x, n, d, a, b)
    pseudo = S.cochran_armitage(np.r_[a, x], np.r_[a + b, n], np.r_[0.0, d])
    assert with_prior == pytest.approx(pseudo)


def test_tarone_historical_information_sharpens_the_test():
    # a low historical rate makes the same treated incidences more surprising
    x, n, d = np.array([1, 2, 3, 6]), np.array([50, 50, 50, 50]), np.arange(4.0)
    _, p_ca = S.cochran_armitage(x, n, d)
    _, p_low = S.tarone_trend(x, n, d, a=1.0, b=199.0)
    _, p_high = S.tarone_trend(x, n, d, a=20.0, b=180.0)
    assert p_low < p_ca < p_high


def test_tarone_full_replacement_of_control():
    # the control may be absent (n = 0); the prior alone acts as control
    z, p = S.tarone_trend([0, 5, 10], [0, 50, 50], [0, 1, 2], a=2.0, b=98.0)
    assert z > 3 and p < 0.001


# --- beta-binomial ------------------------------------------------------------

def test_betabinom_loglik_matches_scipy():
    x, n, a, b = np.array([0, 3, 7, 1]), np.array([50, 49, 50, 48]), 1.7, 30.0
    ref = sps.betabinom.logpmf(x, n, a, b).sum()
    assert S.betabinom_loglik(a, b, x, n) == pytest.approx(ref)


def test_beta_binomial_fit_is_the_likelihood_maximum():
    rng = np.random.default_rng(0)
    p = rng.beta(2.0, 38.0, size=300)
    n = np.full(300, 50)
    x = rng.binomial(n, p)
    a, b = S.fit_beta_binomial(x, n)
    assert a / (a + b) == pytest.approx(0.05, rel=0.1)
    assert 25 < a + b < 70                       # true a + b = 40
    # no point on a grid around the estimate has a higher likelihood
    best = S.betabinom_loglik(a, b, x, n)
    for fa in (0.9, 1.1):
        for fb in (0.9, 1.1):
            assert S.betabinom_loglik(a * fa, b * fb, x, n) <= best + 1e-6


def test_beta_binomial_without_overdispersion_hits_the_pooling_cap():
    rng = np.random.default_rng(1)
    n = np.full(200, 50)
    x = rng.binomial(n, 0.1)
    a, b = S.fit_beta_binomial(x, n)
    assert a + b > 0.5 * n.sum()                 # near the cap = total historical animals
    assert a / (a + b) == pytest.approx(x.sum() / n.sum(), rel=0.02)


# --- DerSimonian-Laird and MAP -------------------------------------------------

def test_dersimonian_laird_matches_statsmodels():
    y = np.array([0.10, 0.50, 0.30, 0.90, -0.20, 0.40])
    v = np.array([0.04, 0.05, 0.03, 0.06, 0.05, 0.02])
    ref = combine_effects(y, v, method_re="dl")
    re = S.dersimonian_laird(y, v)
    assert re.mu == pytest.approx(ref.mean_effect_re)
    assert re.se_mu == pytest.approx(ref.sd_eff_w_re)
    assert re.tau2 == pytest.approx(ref.tau2)


def test_dersimonian_laird_homogeneous_groups_give_zero_tau():
    re = S.dersimonian_laird(np.array([1.0, 1.0, 1.0]), np.array([0.1, 0.2, 0.3]))
    assert re.tau2 == 0.0 and re.mu == pytest.approx(1.0)


def test_historical_random_effects_recovers_between_study_variance():
    rng = np.random.default_rng(2)
    groups = [rng.normal(loc=rng.normal(10, 1.0), scale=2.0, size=10) for _ in range(400)]
    re, s2 = S.historical_random_effects(groups)
    assert s2 == pytest.approx(4.0, rel=0.05)
    assert re.tau2 == pytest.approx(1.0, rel=0.25)
    assert re.mu == pytest.approx(10.0, abs=0.15)


def test_map_predictive_test_inflates_variance_with_tau():
    rng = np.random.default_rng(3)
    groups = [rng.normal(11.0, 1.0, size=10) for _ in range(3)]
    narrow = S.map_predictive_test(groups, S.RandomEffects(10.0, 0.05, 0.0, 30))
    wide = S.map_predictive_test(groups, S.RandomEffects(10.0, 0.05, 1.0, 30))
    assert (narrow.pvalue < wide.pvalue).all()


def test_shared_control_contrasts_reduce_to_dunnett():
    # with V_c = s^2/n0 and V_i = s^2/n_i the shared-control test is Dunnett's
    rng = np.random.default_rng(4)
    control = rng.normal(size=12)
    groups = [rng.normal(0.8, size=10), rng.normal(0.2, size=11)]
    ref = S.dunnett(groups, control)
    s2, df = S._pooled_within_variance([control, *groups])
    ns = np.array([10.0, 11.0])
    diff = np.array([g.mean() for g in groups]) - control.mean()
    t, p = S.shared_control_contrasts(diff, s2 / ns, s2 / 12, df)
    assert np.allclose(t, ref.statistic)
    assert np.allclose(p, ref.pvalue, atol=2e-4)


# --- robust MAP posterior ---------------------------------------------------------

def _numeric_posterior_moments(prior: S.NormalMixture, ybar: float, var_ybar: float) -> tuple[float, float]:
    def density(theta):
        prior_pdf = sum(w * sps.norm.pdf(theta, m, np.sqrt(v))
                        for w, m, v in zip(prior.weights, prior.means, prior.variances))
        return prior_pdf * sps.norm.pdf(ybar, theta, np.sqrt(var_ybar))
    lo, hi = -60, 60
    z = integrate.quad(density, lo, hi, limit=400)[0]
    m1 = integrate.quad(lambda t: t * density(t), lo, hi, limit=400)[0] / z
    m2 = integrate.quad(lambda t: t * t * density(t), lo, hi, limit=400)[0] / z
    return m1, m2 - m1**2


@pytest.mark.parametrize("ybar", [0.0, 1.5, 6.0])
def test_mixture_posterior_matches_numerical_integration(ybar):
    prior = S.robust_map_prior(S.RandomEffects(0.0, 0.2, 0.3, 20), unit_var=4.0, vague_weight=0.2)
    post = S.normal_mixture_posterior(prior, ybar, 0.5)
    mean, var = _numeric_posterior_moments(prior, ybar, 0.5)
    assert post.mean == pytest.approx(mean, abs=1e-6)
    assert post.var == pytest.approx(var, rel=1e-5)


def test_robust_prior_discounts_conflicting_data_less():
    re = S.RandomEffects(0.0, 0.1, 0.1, 20)
    informative = S.robust_map_prior(re, unit_var=4.0, vague_weight=0.0)
    robust = S.robust_map_prior(re, unit_var=4.0, vague_weight=0.2)
    # concurrent data far from the prior: the robust posterior follows the data
    far_inf = S.normal_mixture_posterior(informative, 4.0, 0.4).mean
    far_rob = S.normal_mixture_posterior(robust, 4.0, 0.4).mean
    assert far_rob > 3.0 > far_inf


# --- ANCOVA --------------------------------------------------------------------------

def test_ancova_t_statistics_match_statsmodels_ols():
    rng = np.random.default_rng(5)
    sizes = [10, 9, 10, 8]
    xs = [rng.normal(150, 10, size=n) for n in sizes]
    ys = [20 + 0.1 * x + shift + rng.normal(0, 1, size=x.size) for x, shift in zip(xs, [0, 0.3, 1.0, 2.0])]
    res = S.ancova_dunnett(ys[1:], xs[1:], ys[0], xs[0])
    group = np.concatenate([np.full(n, i) for i, n in enumerate(sizes)])
    design = np.column_stack([np.ones(group.size)] + [(group == i).astype(float) for i in (1, 2, 3)]
                             + [np.concatenate(xs)])
    ols = sm.OLS(np.concatenate(ys), design).fit()
    assert np.allclose(res.statistic, ols.tvalues[1:4])
    assert res.df == int(ols.df_resid)


def test_ancova_removes_a_covariate_imbalance():
    # treated animals start heavier but there is no treatment effect: Dunnett
    # rejects often, ANCOVA keeps its nominal 5% false-positive rate
    rng = np.random.default_rng(6)
    rej_dunnett = rej_ancova = 0
    reps = 400
    for _ in range(reps):
        x0, x1 = rng.normal(150, 8, 10), rng.normal(160, 8, 10)
        y0, y1 = 0.2 * x0 + rng.normal(0, 1, 10), 0.2 * x1 + rng.normal(0, 1, 10)
        rej_dunnett += S.dunnett([y1], y0).pvalue[0] < 0.05
        rej_ancova += S.ancova_dunnett([y1], [x1], y0, x0).pvalue[0] < 0.05
    assert rej_dunnett / reps > 0.5
    assert rej_ancova / reps == pytest.approx(0.05, abs=0.03)
