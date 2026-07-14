"""
Tail approximations: how a permutation test stops needing so many permutations.

Counting exceedances cannot report a p-value below ``1/n_perm``. That floor -- not
the statistic, and not the number of elements -- is what forces a permutation test
to run many thousands of permutations, and it caps FDR too, since FDR is computed
from the uncorrected p-values. Two fits remove it.

**Gamma**, for the FWE-corrected p-values. The null of the *maximum* statistic is
fitted from its first three moments rather than counted, so a corrected p-value is
not floored and converges with far fewer permutations.

**Generalised Pareto**, for the uncorrected p-values. A GPD is fitted to the tail
of each element's own permutation distribution. Its *shape* is pooled across
elements: they all carry the same statistic under the same design, so they differ
in scale, not in shape, and fitting the shape per element would spend a hundred
tail values on a number they all share.

Winkler et al., Faster permutation inference in brain imaging, NeuroImage
141:502-516, 2016.
"""

import numpy as np
from scipy import stats

__all__ = [
    "gamma_pvalue",
    "pareto_pvalue",
    "gpd_fit_pwm",
    "gpd_sf",
]


# ----------------------------------------------------------------------------
# Gamma: the null of the maximum statistic
# ----------------------------------------------------------------------------
def gamma_pvalue(stat, null_max, floor=None):
    """FWE-corrected p-values from a Gamma fit to the maximum distribution.

    Parameters
    ----------
    stat : ndarray
        Observed statistic, one value per element.
    null_max : ndarray, shape (n_perm,)
        The maximum of the statistic over the image, one per permutation.
    floor : float, optional
        Smallest p-value to report. Defaults to ``1/n_perm``, which is what
        counting could have resolved, so the fit is never asked to extrapolate
        further than the permutations can support.

    Returns
    -------
    ndarray
        p-values, same shape as ``stat``.

    Notes
    -----
    A Gamma is fitted by matching the first three moments, which is what makes
    this cheap and stable: with shape ``k``, scale ``theta`` and location ``mu``,

        mean = mu + k*theta,   var = k*theta**2,   skew = 2/sqrt(k)

    so the three moments invert in closed form. A permutation distribution of a
    maximum is right-skewed, which is exactly the regime a Gamma is for; if the
    sample skewness comes out non-positive the fit has nothing to hold on to and
    the empirical count is returned instead.
    """
    stat = np.asarray(stat, dtype=float)
    null_max = np.asarray(null_max, dtype=float).ravel()
    n_perm = null_max.size

    if floor is None:
        floor = 1.0 / n_perm

    mean = null_max.mean()
    var = null_max.var(ddof=1)
    skew = stats.skew(null_max, bias=False)

    if not np.isfinite(skew) or skew <= 0 or var <= 0:
        # nothing for a Gamma to fit; fall back to counting
        return _count_sf(stat, null_max)

    k = 4.0 / skew**2
    theta = np.sqrt(var / k)
    mu = mean - k * theta

    p = stats.gamma.sf(stat, a=k, loc=mu, scale=theta)
    return np.clip(p, floor, 1.0)


def _count_sf(stat, null):
    """P(null >= stat), by counting."""
    null = np.sort(np.asarray(null, dtype=float))
    # number of null values >= each stat
    n_ge = null.size - np.searchsorted(null, stat, side="left")
    return n_ge / null.size


# ----------------------------------------------------------------------------
# Generalised Pareto: the tail of each element's own distribution
# ----------------------------------------------------------------------------
def gpd_fit_pwm(y):
    """Generalised Pareto parameters by probability-weighted moments.

    Parameters
    ----------
    y : ndarray, shape (m, n)
        Exceedances above a threshold, ascending down each column, one column per
        element.

    Returns
    -------
    k, sigma : ndarray, shape (n,)
        Shape and scale, parameterised as in Hosking & Wallis (1987):

            F(y) = 1 - (1 - k*y/sigma)**(1/k)

    Notes
    -----
    Closed form, hence every element is fitted at once rather than one at a time.
    """
    y = np.asarray(y, dtype=float)
    m = y.shape[0]
    j = np.arange(1, m + 1)[:, None]

    a0 = y.mean(axis=0)
    a1 = (y * ((m - j) / (m * (m - 1)))).sum(axis=0)

    d = a0 - 2 * a1
    with np.errstate(divide="ignore", invalid="ignore"):
        k = a0 / d - 2
        sigma = 2 * a0 * a1 / d
    return k, sigma


def gpd_sf(z, k, sigma):
    """Survival function of the Generalised Pareto."""
    z = np.asarray(z, dtype=float)
    k = np.broadcast_to(np.asarray(k, dtype=float), z.shape)
    sigma = np.broadcast_to(np.asarray(sigma, dtype=float), z.shape)

    out = np.zeros(z.shape, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        lin = np.abs(k) < 1e-8              # the k -> 0 limit is exponential
        out[lin] = np.exp(-z[lin] / sigma[lin])

        q = 1 - k * z / sigma
        q = np.where(q < 0, 0.0, q)         # beyond the end point of a k > 0 fit
        nl = ~lin
        out[nl] = q[nl] ** (1.0 / k[nl])

    return out


def _gpd_anderson(y, k, sigma):
    """Anderson-Darling statistic of a fitted GPD against the tail it was fitted to.

    Used only to rank candidate tail sizes against each other, so what matters is
    that it grows when the fit describes the observed tail less well.
    """
    m = y.shape[0]
    j = np.arange(1, m + 1)[:, None]

    F = 1 - gpd_sf(y, k, sigma)
    F = np.clip(F, 1e-12, 1 - 1e-12)

    return -m - ((2 * j - 1) * (np.log(F) + np.log(1 - F[::-1]))).sum(axis=0) / m


def pareto_pvalue(
    stat,
    null,
    *,
    n_exc_min=25,
    n_tail=100,
    n_shape=20000,
    ad_max=5.0,
):
    """Uncorrected p-values, resolved below the ``1/n_perm`` floor of counting.

    Parameters
    ----------
    stat : ndarray, shape (n_elements,)
        Observed statistic. Must already be an *upper* tail: flip the sign of the
        elements whose effect is negative before calling, so that both sides are
        the same question.
    null : ndarray, shape (n_perm, n_elements)
        The permutation distribution of each element, same sign convention. Only
        its upper tail is used, so it is enough to pass the largest ``n_tail``
        values per element.
    n_exc_min : int, default 25
        Elements whose statistic was exceeded at least this often are left to
        plain counting: the count is then already precise enough (25 exceedances
        is a relative standard error of 20%) and there is nothing for a fit to add.
    n_tail : int, default 100
        Largest tail the fit may use.
    n_shape : int, default 20000
        How many elements the shape is pooled over.
    ad_max : float, default 5.0
        A tail that no candidate fit describes this well is counted, not
        extrapolated. Deliberately loose: it is there to catch a fit that has gone
        wrong, not to test a hypothesis about the tail.

    Returns
    -------
    ndarray, shape (n_elements,)
        p-values.
    """
    stat = np.asarray(stat, dtype=float).ravel()
    null = np.asarray(null, dtype=float)
    n_perm, n_elem = null.shape

    if stat.size != n_elem:
        raise ValueError(
            "stat has %d elements but null describes %d" % (stat.size, n_elem)
        )

    cnt = (null >= stat[None, :]).sum(axis=0)
    p = cnt / n_perm

    K = min(n_tail, n_perm - 1)
    if K < 30:
        return p

    # only the elements counting cannot resolve are worth fitting
    need = (cnt < n_exc_min) & np.isfinite(stat)
    if not need.any():
        return p

    # descending tails: the fitted elements, and a sample of ALL of them for the
    # shape. The shape is pooled over the whole image on purpose -- the elements
    # being fitted were selected for carrying a large statistic, and selecting on
    # that would bias the tails it is estimated from.
    G = -np.sort(-null, axis=0)[:K]                   # (K, n_elem), descending
    Gn = G[:, need]
    ob = stat[need]

    step = max(1, n_elem // n_shape)
    Gs = G[:, ::step]

    best_ad = np.full(ob.shape, np.inf)
    best_p = p[need].copy()

    cands = sorted({int(round(K * f)) for f in np.arange(0.9, 0.29, -0.1)})
    cands = [m for m in cands if 20 <= m <= K - 1]

    for m in cands:
        # pooled shape, from the sampled elements
        us = (Gs[m - 1] + Gs[m]) / 2
        ys = Gs[:m][::-1] - us                       # ascending exceedances
        ks, _ = gpd_fit_pwm(ys)
        ks = ks[np.isfinite(ks)]
        if ks.size == 0:
            continue
        k_pool = float(np.median(ks))

        u = (Gn[m - 1] + Gn[m]) / 2
        y = Gn[:m][::-1] - u

        # the pooled shape, and the exponential limit as a safety net: a GPD with
        # k > 0 has a finite upper end point at u + sigma/k, and where the observed
        # statistic lies beyond it the fit returns exactly zero. The observation is
        # itself proof that the end point is wrong, so such a fit is rejected
        # rather than believed. The exponential has infinite support and always
        # answers.
        variants = []
        if np.isfinite(k_pool) and k_pool > -1:
            variants.append(np.full(ob.shape, k_pool))
        variants.append(np.zeros(ob.shape))

        for kk in variants:
            # given the shape, the scale follows from the mean of the exceedances,
            # because the GPD has mean sigma/(1 + k)
            sigma = y.mean(axis=0) * (1 + kk)

            with np.errstate(invalid="ignore", divide="ignore"):
                pp = (m / n_perm) * gpd_sf(ob - u, kk, sigma)
                ad = _gpd_anderson(y, kk, sigma)

            ok = (
                (sigma > 0)
                & (ob >= u)
                & np.isfinite(pp)
                & (pp > 0)
                & (pp <= 1)
                & np.isfinite(ad)
            )
            upd = ok & (ad < best_ad)
            best_ad = np.where(upd, ad, best_ad)
            best_p = np.where(upd, pp, best_p)

    trust = np.isfinite(best_ad) & (best_ad < ad_max)
    out = p.copy()
    out[need] = np.where(trust, best_p, p[need])
    return out
