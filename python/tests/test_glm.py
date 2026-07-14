"""
Does the accelerated GLM return what the unaccelerated one returns?

The whole point of :class:`tfce.glm.PermutedGLM` is that it never forms the
permuted dataset and never forms the residual matrix. That is an algebraic
rearrangement, not an approximation, so it must agree with an honest refit to
floating-point precision. The reference here is the naive thing: permute the data,
call ``lstsq``, compute the residuals, form the t.
"""

import numpy as np
import pytest
from scipy import stats

from tfce.glm import PermutedGLM, partition_design


def naive_t(Y, X, c, perm=None):
    """The obvious implementation. Slow, and deliberately not clever."""
    if perm is not None:
        Y = Y[:, perm]

    beta, *_ = np.linalg.lstsq(X, Y.T, rcond=None)      # (r, n_elem)
    resid = Y.T - X @ beta
    df = X.shape[0] - np.linalg.matrix_rank(X)
    res_ms = (resid**2).sum(axis=0) / df

    pinvX = np.linalg.pinv(X)
    var_c = c @ (pinvX @ pinvX.T) @ c
    return (c @ beta) / np.sqrt(res_ms * var_c)


def naive_f(Y, X, C, perm=None):
    if perm is not None:
        Y = Y[:, perm]

    beta, *_ = np.linalg.lstsq(X, Y.T, rcond=None)
    resid = Y.T - X @ beta
    df = X.shape[0] - np.linalg.matrix_rank(X)
    res_ms = (resid**2).sum(axis=0) / df

    M = C.T @ np.linalg.pinv(X.T @ X) @ C
    cb = C.T @ beta                                     # (q, n_elem)
    ess = np.einsum("ij,ik,kj->j", cb, np.linalg.pinv(M), cb)
    return (ess / np.linalg.matrix_rank(C)) / res_ms


@pytest.fixture
def data():
    rng = np.random.default_rng(7)
    n, n_elem = 40, 300
    g = np.repeat([1.0, 0.0], n // 2)
    age = rng.standard_normal(n)
    tiv = rng.standard_normal(n)
    X = np.column_stack([g, 1 - g, age, tiv])
    Y = rng.standard_normal((n_elem, n)) + 0.3 * g
    return Y, X, rng


def test_t_matches_a_naive_refit(data):
    Y, X, rng = data
    c = np.array([1.0, -1.0, 0.0, 0.0])
    model = PermutedGLM(Y, X, c)

    for _ in range(10):
        perm = rng.permutation(X.shape[0])
        got = model.fit(perm)
        want = naive_t(Y, X, c, perm)
        assert np.allclose(got, want, rtol=1e-9, atol=1e-9)


def test_unpermuted_fit_matches(data):
    Y, X, _ = data
    c = np.array([1.0, -1.0, 0.0, 0.0])
    model = PermutedGLM(Y, X, c)
    assert np.allclose(model.fit(None), naive_t(Y, X, c), rtol=1e-9, atol=1e-9)


def test_f_matches_a_naive_refit():
    rng = np.random.default_rng(8)
    n, n_elem = 45, 200
    grp = np.arange(n) % 3
    X = np.column_stack([(grp == k).astype(float) for k in range(3)])
    X = np.column_stack([X, rng.standard_normal(n)])          # a nuisance column
    Y = rng.standard_normal((n_elem, n))

    C = np.array([[1.0, 0.0], [-1.0, 1.0], [0.0, -1.0], [0.0, 0.0]])
    model = PermutedGLM(Y, X, C)

    for _ in range(6):
        perm = rng.permutation(n)
        assert np.allclose(
            model.fit(perm), naive_f(Y, X, C, perm), rtol=1e-8, atol=1e-8
        )


def test_sign_flipping_matches_a_naive_refit():
    """A one-sample design is permuted by flipping signs, not by shuffling rows."""
    rng = np.random.default_rng(9)
    n, n_elem = 30, 200
    X = np.ones((n, 1))
    Y = rng.standard_normal((n_elem, n)) + 0.2
    c = np.array([1.0])

    model = PermutedGLM(Y, X, c)

    for _ in range(10):
        signs = rng.choice([-1.0, 1.0], size=n)
        got = model.fit_signs(signs)
        want = naive_t(Y * signs[None, :], X, c)
        assert np.allclose(got, want, rtol=1e-9, atol=1e-9)


def test_t_is_calibrated_under_the_null():
    """Sanity: the statistic must actually be a t under the null."""
    rng = np.random.default_rng(21)
    n, n_elem = 30, 5000
    g = np.repeat([1.0, 0.0], n // 2)
    X = np.column_stack([g, 1 - g])
    Y = rng.standard_normal((n_elem, n))

    t = PermutedGLM(Y, X, np.array([1.0, -1.0])).fit(None)
    p = stats.kstest(t, "t", args=(n - 2,)).pvalue
    assert p > 0.01, f"the statistic is not t-distributed under the null (p={p})"


def test_partition_design():
    X = np.zeros((10, 4))
    idx_i, idx_n = partition_design(X, np.array([1.0, -1.0, 0.0, 0.0]))
    assert list(idx_i) == [0, 1]
    assert list(idx_n) == [2, 3]
