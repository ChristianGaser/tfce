# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2020-2026, Christian Gaser. See LICENSE.
"""
The worked examples in docs/permutation.md, actually run.

Documentation that is not executed rots, and a worked example that does not work is
worse than none at all -- it is a trap for the one reader who trusts it. So the two
pipelines that page describes are run here, on data with a planted effect, and are
required to find it while controlling the false-positive rate.

If you change docs/permutation.md, change this too.
"""

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

import tfce
from tfce.glm import PermutedGLM
from tfce.tails import gamma_pvalue, pareto_pvalue

SHAPE = (14, 14, 14)
N_PERM = 200
N_TAIL = 50


def make_data(rng, n_subj, effect=0.0, blob=None):
    Y = np.stack(
        [gaussian_filter(rng.standard_normal(SHAPE), 1.8) for _ in range(n_subj)]
    )
    Y /= Y.std()
    if effect and blob is not None:
        Y[:, blob] += effect
    return Y


def blob_mask():
    m = np.zeros(SHAPE, dtype=bool)
    m[5:9, 5:9, 5:9] = True
    return m


def run_permutation_test(model, permute, rng, n_elem):
    """Exactly the loop docs/permutation.md describes."""
    t0 = model.fit(None).reshape(SHAPE)
    tfce0 = tfce.tfce(t0, E=0.5, H=2.0)

    sgn = np.sign(tfce0.ravel())
    sgn[sgn == 0] = 1
    obs = np.abs(tfce0).ravel()

    null_max = np.empty(N_PERM)
    tail = np.full((N_TAIL, n_elem), -np.inf)
    cnt = np.zeros(n_elem)

    for p in range(N_PERM):
        tp = tfce.tfce(permute(rng).reshape(SHAPE), E=0.5, H=2.0)

        null_max[p] = np.abs(tp).max()

        v = tp.ravel() * sgn

        # exact, over every permutation -- the truncated tail cannot answer this
        cnt += v >= obs

        tail[0] = v
        tail.sort(axis=0)

    p_fwe = gamma_pvalue(obs, null_max).reshape(SHAPE)
    p_unc = pareto_pvalue(obs, tail, cnt, N_PERM).reshape(SHAPE)

    return tfce0, p_fwe, p_unc


class TestOneSample:
    """The sign-flipping example."""

    def _model(self, Y):
        n_subj = Y.shape[0]
        flat = Y.reshape(n_subj, -1).T          # (n_elements, n_subjects)
        X = np.ones((n_subj, 1))
        c = np.array([1.0])
        return PermutedGLM(flat, X, c), flat.shape[0]

    def test_detects_a_real_effect(self):
        rng = np.random.default_rng(0)
        blob = blob_mask()
        Y = make_data(rng, 20, effect=1.0, blob=blob)

        model, n_elem = self._model(Y)
        permute = lambda r: model.fit_signs(r.choice([-1.0, 1.0], size=20))

        _, p_fwe, p_unc = run_permutation_test(model, permute, rng, n_elem)

        assert p_fwe[blob].min() < 0.05, "the planted effect was not detected"
        assert np.all(p_fwe >= 0) and np.all(p_fwe <= 1)
        assert np.all(p_unc > 0), "an uncorrected p-value of zero is unusable"

    def test_controls_false_positives_under_the_null(self):
        rng = np.random.default_rng(1)
        Y = make_data(rng, 20)                  # no effect anywhere

        model, n_elem = self._model(Y)
        permute = lambda r: model.fit_signs(r.choice([-1.0, 1.0], size=20))

        _, p_fwe, p_unc = run_permutation_test(model, permute, rng, n_elem)

        # FWE: nothing anywhere should survive, most of the time
        assert p_fwe.min() > 0.01

        # The uncorrected p-value is one-sided and conditioned on the sign of the
        # observed effect, so under the null it is uniform on (0, 0.5] -- NOT on
        # (0, 1]. It can never approach 1: a large positive statistic cannot also
        # be unusually small. So the upper tail to check is 0.45, and ~10% of the
        # elements should clear it. This is the check docs/permutation.md prints.
        assert p_unc.max() <= 0.6, (
            f"p_unc reaches {p_unc.max():.2f}; a sign-conditioned one-sided "
            "p-value cannot exceed ~0.5, so the construction is wrong"
        )

        frac_high = np.mean(p_unc > 0.45)
        assert 0.05 < frac_high < 0.18, (
            f"{frac_high:.1%} of uncorrected p exceed 0.45, expected ~10% — "
            "the permutation null is the wrong width"
        )


class TestTwoSample:
    """The group-permutation example, with a nuisance covariate."""

    def test_detects_a_real_effect_with_a_nuisance_covariate(self):
        rng = np.random.default_rng(2)
        n_per, blob = 10, blob_mask()
        n_subj = 2 * n_per

        Y = make_data(rng, n_subj)
        g = np.repeat([1.0, 0.0], n_per)
        Y[g == 1] += 1.2 * blob                 # effect in group A only

        flat = Y.reshape(n_subj, -1).T
        age = rng.standard_normal(n_subj)
        X = np.column_stack([g, 1 - g, age])
        c = np.array([1.0, -1.0, 0.0])

        model = PermutedGLM(flat, X, c)
        permute = lambda r: model.fit(r.permutation(n_subj))

        _, p_fwe, p_unc = run_permutation_test(model, permute, rng, flat.shape[0])

        assert p_fwe[blob].min() < 0.05
        assert np.all(p_unc > 0)


def test_pareto_reaches_below_the_counting_floor():
    """The claim the docs make: p below 1/n_perm, which counting cannot reach."""
    rng = np.random.default_rng(3)
    blob = blob_mask()
    Y = make_data(rng, 20, effect=1.5, blob=blob)

    n_subj = 20
    flat = Y.reshape(n_subj, -1).T
    model = PermutedGLM(flat, np.ones((n_subj, 1)), np.array([1.0]))
    permute = lambda r: model.fit_signs(r.choice([-1.0, 1.0], size=n_subj))

    _, _, p_unc = run_permutation_test(model, permute, rng, flat.shape[0])

    floor = 1.0 / N_PERM
    assert p_unc.min() < floor, (
        f"smallest uncorrected p is {p_unc.min():.2e}, which counting could have "
        f"reached on its own (floor {floor:.2e})"
    )
