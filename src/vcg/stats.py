"""Statistical tests used to compare treated groups with a control group.

Each endpoint type gets the test a toxicologist would reach for:

* continuous endpoints (body weight, organ weights, clinical pathology):
  Dunnett's many-to-one t-test, two-sided, plus the Jonckheere-Terpstra trend
  test that NTP uses to flag dose-related trends;
* lesion incidences in short studies: Fisher's exact test, per dose group;
* neoplasm incidences in chronic studies: the poly-3 test (Bailer & Portier
  1988) with the Bieler & Williams (1993) variance, as used in NTP technical
  reports.

The Dunnett p-value is computed from the exact multivariate-t distribution
with the product correlation structure implied by the group sizes.  Because the
correlations are of the form rho_ij = lambda_i * lambda_j, the k-dimensional
integral collapses to a two-dimensional one, which we evaluate with Gauss
quadrature.  This is a few hundred times faster than ``scipy.stats.dunnett``
and agrees with it to about 1e-3 (see tests), which is what makes resampling
thousands of virtual control groups affordable.  The survival function is also
tabulated per correlation structure, since every virtual control group of a
given endpoint shares the same group sizes.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy import special, stats

# ----------------------------------------------------------------------------
# Dunnett's test
# ----------------------------------------------------------------------------

_GH_NODES, _GH_WEIGHTS = np.polynomial.hermite_e.hermegauss(48)  # probabilists' Hermite
_GH_WEIGHTS = _GH_WEIGHTS / np.sqrt(2 * np.pi)                   # integrate against N(0,1)
# grid on which the adjusted survival function is tabulated per (lambdas, df)
_T_GRID = np.concatenate([np.linspace(0.0, 6.0, 241), np.linspace(6.05, 15.0, 60)])


@lru_cache(maxsize=None)
def _chi_nodes(df: int, n_nodes: int = 40) -> tuple[np.ndarray, np.ndarray]:
    """Quadrature nodes/weights for s = sqrt(chi2_df / df) on its density.

    Generalised Gauss-Laguerre quadrature in x = s^2 * df / 2 integrates
    x^(df/2-1) e^-x exactly for polynomial integrands, which is the chi-square
    kernel; we then map back to s.
    """
    alpha = df / 2.0 - 1.0
    x, w = special.roots_genlaguerre(n_nodes, alpha)
    w = w / special.gamma(df / 2.0)          # normalise the Gamma(df/2, 1) density
    s = np.sqrt(2.0 * x / df)
    return s, w


def _dunnett_sf_quadrature(t: np.ndarray, lam: np.ndarray, df: int) -> np.ndarray:
    denom = np.sqrt(1.0 - lam**2)
    s, ws = _chi_nodes(int(df))
    z, wz = _GH_NODES, _GH_WEIGHTS
    ts = t[:, None, None] * s[None, :, None]            # (t, s, 1)
    prod = np.ones((t.size, s.size, z.size))
    for li, di in zip(lam, denom):
        upper = (ts + li * z[None, None, :]) / di
        lower = (-ts + li * z[None, None, :]) / di
        prod *= special.ndtr(upper) - special.ndtr(lower)
    cdf = np.einsum("tsz,s,z->t", prod, ws, wz)
    return np.clip(1.0 - cdf, 0.0, 1.0)


@lru_cache(maxsize=8192)
def _dunnett_table(lam_key: tuple[float, ...], df: int) -> np.ndarray:
    return _dunnett_sf_quadrature(_T_GRID, np.array(lam_key), df)


def dunnett_sf(t: np.ndarray, lambdas: np.ndarray, df: int, exact: bool = False) -> np.ndarray:
    """P(max_i |T_i| > t) for Dunnett's statistics with correlations lambda_i*lambda_j.

    ``lambdas`` has one entry per treated group, lambda_i = sqrt(n_i/(n_0+n_i)).
    ``t`` may be an array; a two-sided adjusted p-value is returned per entry.

    By default the survival function is tabulated once per (lambdas, df) on a
    fine grid and interpolated linearly (error below 1e-4), which makes the
    thousands of resampled tests per endpoint cheap.  ``exact=True`` evaluates
    the quadrature directly.
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    lam = np.asarray(lambdas, dtype=float)
    if exact:
        return _dunnett_sf_quadrature(t, lam, int(df))
    table = _dunnett_table(tuple(np.round(lam, 6)), int(df))
    return np.interp(t, _T_GRID, table, right=0.0)


@dataclass(frozen=True)
class DunnettResult:
    statistic: np.ndarray   # t statistic per treated group (treated minus control)
    pvalue: np.ndarray      # two-sided, adjusted for the number of comparisons
    df: int


def dunnett(groups: list[np.ndarray], control: np.ndarray) -> DunnettResult:
    """Dunnett's test of each group against ``control`` (two-sided)."""
    control = np.asarray(control, dtype=float)
    groups = [np.asarray(g, dtype=float) for g in groups]
    n0 = control.size
    ns = np.array([g.size for g in groups])
    k = len(groups)
    df = int(n0 + ns.sum() - k - 1)
    if df < 1 or n0 < 2 or (ns < 2).any():
        return DunnettResult(np.full(k, np.nan), np.full(k, np.nan), df)
    ss = ((control - control.mean()) ** 2).sum() + sum(((g - g.mean()) ** 2).sum() for g in groups)
    mse = ss / df
    if mse <= 0:
        return DunnettResult(np.zeros(k), np.ones(k), df)
    means = np.array([g.mean() for g in groups])
    tstat = (means - control.mean()) / np.sqrt(mse * (1.0 / ns + 1.0 / n0))
    lambdas = np.sqrt(ns / (n0 + ns))
    pvals = dunnett_sf(np.abs(tstat), lambdas, df)
    return DunnettResult(tstat, pvals, df)


def dunnett_many(groups: list[np.ndarray], controls: np.ndarray) -> DunnettResult:
    """Dunnett's test of fixed treated ``groups`` against each row of ``controls``.

    ``controls`` is an (R, n0) matrix of R alternative control groups of equal
    size (for example R virtual control groups).  Because every row has the same
    size, all R x k statistics share one correlation structure and one degrees of
    freedom, so the adjusted p-values come from a single quadrature call.
    Returns arrays of shape (R, k).
    """
    controls = np.atleast_2d(np.asarray(controls, dtype=float))
    groups = [np.asarray(g, dtype=float) for g in groups]
    r, n0 = controls.shape
    ns = np.array([g.size for g in groups])
    k = len(groups)
    df = int(n0 + ns.sum() - k - 1)
    if df < 1 or n0 < 2 or (ns < 2).any():
        return DunnettResult(np.full((r, k), np.nan), np.full((r, k), np.nan), df)
    ss_groups = sum(((g - g.mean()) ** 2).sum() for g in groups)
    ss_ctrl = controls.var(axis=1, ddof=0) * n0                      # (R,)
    mse = (ss_groups + ss_ctrl) / df                                 # (R,)
    means = np.array([g.mean() for g in groups])                     # (k,)
    diff = means[None, :] - controls.mean(axis=1)[:, None]           # (R, k)
    se = np.sqrt(np.maximum(mse, 1e-300)[:, None] * (1.0 / ns + 1.0 / n0)[None, :])
    tstat = diff / se
    lambdas = np.sqrt(ns / (n0 + ns))
    pvals = dunnett_sf(np.abs(tstat).ravel(), lambdas, df).reshape(r, k)
    pvals[mse <= 0] = 1.0
    return DunnettResult(tstat, pvals, df)


def welch_p_many(a: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """Two-sided Welch t-test p-values of sample ``a`` against each row of ``rows``."""
    a = np.asarray(a, dtype=float)
    rows = np.atleast_2d(np.asarray(rows, dtype=float))
    na, nb = a.size, rows.shape[1]
    if na < 2 or nb < 2:
        return np.full(rows.shape[0], np.nan)
    va = a.var(ddof=1) / na
    vb = rows.var(axis=1, ddof=1) / nb
    denom = va + vb
    t = (a.mean() - rows.mean(axis=1)) / np.sqrt(np.where(denom > 0, denom, np.nan))
    df = denom**2 / (va**2 / (na - 1) + vb**2 / (nb - 1))
    p = 2 * special.stdtr(df, -np.abs(t))
    return np.where(np.isfinite(p), p, 1.0)


# ----------------------------------------------------------------------------
# Jonckheere-Terpstra trend test (normal approximation, ties handled)
# ----------------------------------------------------------------------------

def jonckheere_terpstra(groups: list[np.ndarray]) -> tuple[float, float]:
    """Jonckheere-Terpstra test for a monotone trend across ordered groups.

    ``groups`` are ordered from control upwards.  Returns (z, two-sided p).
    """
    groups = [np.asarray(g, dtype=float) for g in groups if len(g)]
    if len(groups) < 2:
        return np.nan, np.nan
    j = 0.0
    for a in range(len(groups) - 1):
        for b in range(a + 1, len(groups)):
            diff = groups[b][:, None] - groups[a][None, :]
            j += (diff > 0).sum() + 0.5 * (diff == 0).sum()
    ns = np.array([len(g) for g in groups], dtype=float)
    n = ns.sum()
    allv = np.concatenate(groups)
    _, tie_counts = np.unique(allv, return_counts=True)
    t = tie_counts.astype(float)
    mean = (n**2 - (ns**2).sum()) / 4.0
    var = (n * (n - 1) * (2 * n + 5) - (ns * (ns - 1) * (2 * ns + 5)).sum()
           - (t * (t - 1) * (2 * t + 5)).sum()) / 72.0
    var += ((ns * (ns - 1) * (ns - 2)).sum() * (t * (t - 1) * (t - 2)).sum()) / (36.0 * n * (n - 1) * (n - 2))
    var += ((ns * (ns - 1)).sum() * (t * (t - 1)).sum()) / (8.0 * n * (n - 1))
    if var <= 0:
        return 0.0, 1.0
    z = (j - mean) / np.sqrt(var)
    return float(z), float(2 * stats.norm.sf(abs(z)))


# ----------------------------------------------------------------------------
# Fisher's exact test (two-sided, vectorised over tables)
# ----------------------------------------------------------------------------

def fisher_exact_2x2(x1: np.ndarray, n1: np.ndarray, x0: np.ndarray, n0: np.ndarray) -> np.ndarray:
    """Two-sided Fisher exact p-values for many 2x2 tables at once.

    Table per element: treated ``x1`` of ``n1`` affected versus control ``x0`` of
    ``n0``.  Uses the "sum of probabilities not larger than the observed one"
    definition, matching ``scipy.stats.fisher_exact``.
    """
    x1, n1, x0, n0 = (np.atleast_1d(np.asarray(a, dtype=int)) for a in (x1, n1, x0, n0))
    out = np.ones(x1.shape, dtype=float)
    for i in range(x1.size):
        m, nn, k = n1[i] + n0[i], n1[i], x1[i] + x0[i]
        if m == 0 or k == 0 or k == m:
            continue
        support = np.arange(max(0, k - n0[i]), min(k, nn) + 1)
        pmf = stats.hypergeom.pmf(support, m, nn, k)
        p_obs = stats.hypergeom.pmf(x1[i], m, nn, k)
        out[i] = min(1.0, pmf[pmf <= p_obs * (1 + 1e-7)].sum())
    return out


# ----------------------------------------------------------------------------
# Poly-3 test for neoplasm incidence
# ----------------------------------------------------------------------------

def poly_k_weights(time_on_study: np.ndarray, affected: np.ndarray, t_max: np.ndarray | float,
                   k: float = 3.0) -> np.ndarray:
    """Poly-k weights: 1 for affected animals or animals reaching the terminal
    sacrifice, (t/t_max)^k otherwise."""
    t = np.asarray(time_on_study, dtype=float)
    tmax = np.broadcast_to(np.asarray(t_max, dtype=float), t.shape)
    w = np.clip(t / np.where(tmax > 0, tmax, 1.0), 0.0, 1.0) ** k
    return np.where(np.asarray(affected, dtype=bool) | (t >= tmax), 1.0, w)


def _poly_k_group(affected: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    """Return (adjusted n, incidence rate, Bieler-Williams variance of the rate)."""
    x = np.asarray(affected, dtype=float)
    w = np.asarray(weights, dtype=float)
    n = x.size
    n_adj = w.sum()
    if n < 2 or n_adj <= 0:
        return n_adj, np.nan, np.nan
    p = x.sum() / n_adj
    r = x - p * w                       # ratio-estimator residuals
    var = n * ((r - r.mean()) ** 2).sum() / (n - 1) / n_adj**2
    return n_adj, p, var


def poly3_test(aff_treated: np.ndarray, w_treated: np.ndarray,
               aff_control: np.ndarray, w_control: np.ndarray) -> tuple[float, float]:
    """Poly-3 pairwise comparison of a treated group against control.

    Returns (z, one-sided p for an *increase* in incidence), the form reported in
    NTP technical reports.  A continuity correction is not applied.
    """
    n1, p1, v1 = _poly_k_group(aff_treated, w_treated)
    n0, p0, v0 = _poly_k_group(aff_control, w_control)
    if not np.isfinite(p1) or not np.isfinite(p0):
        return np.nan, np.nan
    if (p1 == 0 and p0 == 0) or (p1 == 1 and p0 == 1):
        return 0.0, 1.0
    var = v1 + v0
    if var <= 0:
        return (np.inf if p1 > p0 else -np.inf if p1 < p0 else 0.0), (0.0 if p1 > p0 else 1.0)
    z = (p1 - p0) / np.sqrt(var)
    return float(z), float(stats.norm.sf(z))
