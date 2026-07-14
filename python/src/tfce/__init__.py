# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2020-2026, Christian Gaser. See LICENSE.
"""Exact TFCE and permutation inference.

The core is deliberately framework-free: arrays in, arrays out, numpy and scipy
and nothing else. The pieces are

* :mod:`tfce.core`    -- the exact TFCE transform (max-tree, no step size)
* :mod:`tfce.tails`   -- Gamma and Generalised Pareto tail approximations, which
                         are what let a permutation test report a p-value below
                         1/n_perm and so run far fewer permutations
* :mod:`tfce.glm`     -- the permuted GLM
* :mod:`tfce.nilearn_compat` -- drop-in for nilearn's ``calculate_tfce``
"""

from .core import TfceError, adjacency_from_faces, tfce
from . import glm, tails

try:  # the version lives in pyproject.toml and nowhere else
    from importlib.metadata import PackageNotFoundError, version

    __version__ = version("tfce")
except (ImportError, PackageNotFoundError):  # running from a source tree
    __version__ = "unknown"

__all__ = [
    "tfce",
    "adjacency_from_faces",
    "TfceError",
    "glm",
    "tails",
    "__version__",
]
