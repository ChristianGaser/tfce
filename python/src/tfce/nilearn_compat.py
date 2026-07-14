# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2020-2026, Christian Gaser. See LICENSE.
"""
Drop-in replacement for nilearn's TFCE.

nilearn computes TFCE by stepping the height over a grid and summing
(``nilearn.mass_univariate._utils.calculate_tfce``, 100 steps by default,
following fslmaths). That costs a step size, a discretisation error that depends
on it, and -- because fslmaths does not multiply the partial scores by ``dh`` --
values whose scale depends on the range of the map they came from.

:func:`calculate_tfce` here has the same signature and returns the same thing,
computed exactly: the max-tree integrates each piece of the extent function in
closed form, so there is no ``dh`` and nothing to tune. Monkey-patching it in is
enough to give ``permuted_ols(..., tfce=True)`` an exact transform::

    from nilearn.mass_univariate import _utils
    import tfce.nilearn_compat as tc
    _utils.calculate_tfce = tc.calculate_tfce

Note the connectivity. nilearn builds its neighbourhood with
``generate_binary_structure(3, 1)``, which is 6-connectivity -- faces only. The
MATLAB TFCE toolbox and fslmaths use 26. ``calculate_tfce`` below reads the
neighbourhood out of the ``bin_struct`` it is handed, so it reproduces whichever
one the caller meant.
"""

import numpy as np

from .core import tfce

__all__ = ["calculate_tfce", "connectivity_from_bin_struct"]


def connectivity_from_bin_struct(bin_struct):
    """Read 6 / 18 / 26 out of a scipy ``generate_binary_structure`` array.

    ``generate_binary_structure(3, 1)`` is 6, ``(3, 2)`` is 18 and ``(3, 3)`` is
    26. Rather than guess, this counts the neighbours the structure actually
    marks, so any equivalent array works.
    """
    bs = np.asarray(bin_struct, dtype=bool)
    if bs.shape != (3, 3, 3):
        raise ValueError("bin_struct must be (3, 3, 3), got %r" % (bs.shape,))

    n = int(bs.sum()) - int(bs[1, 1, 1])   # the centre is not its own neighbour
    if n not in (6, 18, 26):
        raise ValueError(
            "bin_struct marks %d neighbours; expected 6, 18 or 26" % n
        )
    return n


def calculate_tfce(arr4d, bin_struct, E=0.5, H=2, dh="auto", two_sided_test=True):
    """Exact TFCE, with nilearn's ``calculate_tfce`` signature.

    Parameters
    ----------
    arr4d : ndarray, shape (X, Y, Z, R)
        Unthresholded statistic maps, one per regressor.
    bin_struct : ndarray, shape (3, 3, 3)
        Connectivity, as nilearn passes it.
    E, H : float
        Extent and height weights.
    dh : ignored
        Accepted so that this is a drop-in, and ignored because there is no step
        size: the integral is evaluated exactly. Passing one is not an error, it
        simply has nothing to change.
    two_sided_test : bool, default True

    Returns
    -------
    ndarray, shape (X, Y, Z, R)
    """
    arr4d = np.asarray(arr4d)
    if arr4d.ndim != 4:
        raise ValueError("arr4d must be 4D, got %r" % (arr4d.shape,))

    conn = connectivity_from_bin_struct(bin_struct)

    return tfce(
        arr4d,
        connectivity=conn,
        E=E,
        H=H,
        two_sided=bool(two_sided_test),
        n_jobs=-1,
    )
