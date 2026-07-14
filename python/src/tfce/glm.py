# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2020-2026, Christian Gaser. See LICENSE.
"""
The permuted GLM.

The point of this module is one piece of algebra. A permutation test refits the
same model thousands of times, and the obvious way to do it -- build the permuted
dataset, then regress -- costs an ``n_elements x n x n`` multiply per permutation
and is by far the dominant term. It is also unnecessary.

Two facts remove it. For any design, ``R = I - X @ pinv(X)`` is a projector, so the
residual sum of squares is

    ResSS = ||y||**2 - beta.T @ (X.T @ X) @ beta

and never needs the residual matrix. And under Freedman-Lane the nuisance columns
are part of the full design, so the full-model residual former already annihilates
them (``R @ Rz == R``) and ``Rz`` drops out of the residuals altogether:

    ResSS = ||R @ Rz @ P @ y||**2 = ||P @ y||**2 - beta.T @ (X.T @ X) @ beta

with ``beta = pinv(X) @ P @ y``. Since a permutation matrix only moves the rows it
touches, ``||P @ y||**2`` is the same for every permutation and is computed once.
The whole permutation therefore collapses into the small ``rank(X) x n`` matrix
``pinv(X) @ P``, and the permuted dataset is never formed at all.
"""

import numpy as np

__all__ = [
    "PermutedGLM",
    "partition_design",
]


def partition_design(design, contrast):
    """Split a design into the columns a contrast tests and the rest.

    Returns ``(idx_interest, idx_nuisance)``. This is the Guttman partition: the
    columns the contrast loads on are the effect of interest, everything else is
    nuisance.
    """
    contrast = np.atleast_2d(np.asarray(contrast, dtype=float))
    if contrast.shape[0] == 1 and contrast.shape[1] != design.shape[1]:
        contrast = contrast.T

    loaded = np.any(contrast != 0, axis=tuple(range(contrast.ndim - 1)))
    idx_interest = np.flatnonzero(loaded)
    idx_nuisance = np.flatnonzero(~loaded)
    return idx_interest, idx_nuisance


class PermutedGLM:
    """Fit ``Y = X @ B`` and re-fit it under permutations of the data.

    Parameters
    ----------
    Y : ndarray, shape (n_elements, n_samples)
        Data. One row per element, one column per observation -- the transpose of
        the way a design matrix is usually written, because the elements are the
        long axis and this keeps them contiguous.
    design : ndarray, shape (n_samples, n_regressors)
    contrast : ndarray, shape (n_regressors,) or (n_regressors, q)
        A vector gives a t-contrast, a matrix an F-contrast.

    Notes
    -----
    The permutation enters through the data (Freedman-Lane), so the design -- and
    therefore ``pinv(X)``, ``X.T @ X`` and the error degrees of freedom -- is the
    same for every permutation and is factorised once, here.
    """

    def __init__(self, Y, design, contrast):
        self.Y = np.asarray(Y, dtype=np.float64)
        self.X = np.asarray(design, dtype=np.float64)
        self.c = np.asarray(contrast, dtype=np.float64)

        if self.c.ndim == 1:
            self.c = self.c[:, None]
        self.is_t = self.c.shape[1] == 1

        n, r = self.X.shape
        if self.Y.shape[1] != n:
            raise ValueError(
                "Y has %d columns but the design has %d rows"
                % (self.Y.shape[1], n)
            )

        self.pinvX = np.linalg.pinv(self.X)
        self.XtX = self.X.T @ self.X
        self.rank = np.linalg.matrix_rank(self.X)
        self.df = n - self.rank

        if self.df <= 0:
            raise ValueError("the design leaves no residual degrees of freedom")

        # permutation-invariant: a permutation matrix only moves rows around
        self.ssy = np.einsum("ij,ij->i", self.Y, self.Y)

        # for an F-contrast, an orthonormal basis of the contrast space leaves the
        # ESS unchanged and makes the middle matrix invertible, which turns the
        # pinv of the textbook formula into an ordinary solve
        if not self.is_t:
            self.Q = np.linalg.qr(self.c)[0]
            self.eidf = np.linalg.matrix_rank(self.c)

    def fit(self, perm=None):
        """Statistic under one permutation.

        Parameters
        ----------
        perm : ndarray of int, or None
            The permutation, read the way numpy reads a fancy index: observation
            ``j`` of the permuted data is observation ``perm[j]`` of the original,
            i.e. the permuted data is ``Y[:, perm]``. ``None`` is the identity.
            A sign-flip is not a permutation of rows and goes through
            :meth:`fit_signs`.

        Returns
        -------
        ndarray, shape (n_elements,)

        Notes
        -----
        The permuted data is never formed. Writing out

            beta = pinv(X) @ Y[:, perm].T

        and substituting ``k = perm[j]`` moves the permutation off the data and
        onto the columns of ``pinv(X)``, where it becomes the *inverse*
        permutation -- that inverse is the whole trick, and getting it the wrong
        way round is silent, because a wrong permutation is still a permutation
        and the statistic still looks perfectly reasonable.
        """
        if perm is None:
            G = self.pinvX
        else:
            perm = np.asarray(perm)
            inv = np.argsort(perm)                # k -> where it came from
            G = self.pinvX[:, inv]

        # a permutation only reorders the observations, so ||y||**2 is untouched
        beta = self.Y @ G.T                       # (n_elem, r): the only big GEMM
        return self._statistic(beta, self.ssy)

    def fit_signs(self, signs):
        """Statistic under one sign-flip, for a one-sample design."""
        signs = np.asarray(signs, dtype=np.float64)
        G = self.pinvX * signs[None, :]
        beta = self.Y @ G.T
        return self._statistic(beta, self.ssy)     # sign-flips are orthogonal too

    def _statistic(self, beta, ssy):
        ss_fit = np.einsum("ij,jk,ik->i", beta, self.XtX, beta)
        res_ss = np.maximum(ssy - ss_fit, 0.0)     # the subtraction can cancel
        res_ms = res_ss / self.df

        if self.is_t:
            c = self.c[:, 0]
            var_c = c @ (self.pinvX @ self.pinvX.T) @ c
            num = beta @ c
            with np.errstate(divide="ignore", invalid="ignore"):
                return num / np.sqrt(res_ms * var_c + np.finfo(float).eps)

        # F: ESS = (Q'b)' inv(Q' inv(XtX) Q) (Q'b), with Q an orthonormal basis
        M = self.Q.T @ np.linalg.pinv(self.XtX) @ self.Q
        Minv = np.linalg.pinv(M)
        cb = beta @ self.Q                          # (n_elem, q)
        ess = np.einsum("ij,jk,ik->i", cb, Minv, cb)
        with np.errstate(divide="ignore", invalid="ignore"):
            return (ess / self.eidf) / res_ms
