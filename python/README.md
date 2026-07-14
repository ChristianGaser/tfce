# tfce (Python)

Exact threshold-free cluster enhancement, and the permutation machinery around it,
as a plain Python package.

This is **not a reimplementation**. The transform is the same C max-tree the
MATLAB toolbox runs — compiled straight from the repository root into a Cython
extension — so the two give bit-identical answers, and there is only ever one
implementation of TFCE in this repository to get right.

```python
import numpy as np
import tfce

# volume: (nx, ny, nz) or (nx, ny, nz, n_permutations)
t = tfce.tfce(stat_map, E=0.5, H=2.0)

# surface: faces are 1-based, as GIFTI stores them
adj = tfce.adjacency_from_faces(faces, n_vertices)
t = tfce.tfce(surf_map, adjacency=adj, E=1.0, H=2.0)

# a block of permutations, one per thread
t = tfce.tfce(perms, E=0.5, H=2.0, n_jobs=-1)
```

## Why exact matters

The TFCE of an element is an integral,

```
TFCE(v) = ∫ e_v(h)^E · h^H dh
```

and implementations normally approximate it by stepping `h` over a grid and
summing. That costs a step size, a discretisation error that depends on it, and an
accuracy parameter the caller has to guess. This one builds the **max-tree** (the
component tree) with union-find, and because the extent function `e_v(h)` is
piecewise constant it integrates each piece in closed form. There is nothing to
tune.

The difference is not academic. Against nilearn's stepped transform on a 60×72×60
volume:

| n_steps | max error vs exact |
| --- | --- |
| 50 | 3.4% |
| 100 *(nilearn's default)* | **1.7%** |
| 200 | 0.8% |
| 400 | 0.4% |

The error halves every time the steps double — first order, exactly as a step-size
approximation must. The exact transform has no such term, and is also faster:

| | one volume | 16 permutations |
| --- | --- | --- |
| nilearn (100 steps) | 0.38 s | 9.1 s (569 ms/perm) |
| exact max-tree | **0.03 s** (14×) | **0.13 s** (8 ms/perm, **72×**) |

## Using it inside nilearn

nilearn already has `permuted_ols(..., tfce=True)`. Its TFCE is
`nilearn.mass_univariate._utils.calculate_tfce`, and `tfce.nilearn_compat`
provides a drop-in with the same signature:

```python
from nilearn.mass_univariate import _utils
import tfce.nilearn_compat as tc

_utils.calculate_tfce = tc.calculate_tfce   # now exact, and much faster
```

One thing to know: nilearn builds its neighbourhood with
`generate_binary_structure(3, 1)`, which is **6-connectivity**. The MATLAB toolbox
and fslmaths use **26**. The drop-in reads the neighbourhood out of the
`bin_struct` it is handed, so it reproduces whichever one the caller meant, and
`tfce.tfce(..., connectivity=6|18|26)` lets you say so directly.

## What's here

| module | what it is |
| --- | --- |
| `tfce.core` | the exact transform: arrays in, arrays out |
| `tfce.tails` | Gamma and Generalised Pareto tail approximations |
| `tfce.glm` | the permuted GLM, which never forms the permuted data |
| `tfce.nilearn_compat` | drop-in for `calculate_tfce` |

The core is deliberately framework-free — numpy and scipy, nothing else. No image
objects, no file I/O, no design parsing. Those belong in a layer above, so that
the core can be vendored anywhere (nilearn included) without dragging a dependency
tree behind it. `nibabel` is an optional extra, wanted only by I/O.

### Fewer permutations

Counting exceedances cannot report a p-value below `1/n_perm`. That floor — not
the statistic — is what forces a permutation test to run many thousands of
permutations, and it caps FDR too, since FDR is computed from the uncorrected
p-values. `tfce.tails` removes it: a **Gamma** fit to the null of the maximum for
the FWE-corrected p-values, and a **Generalised Pareto** fit to each element's own
tail for the uncorrected ones, with the *shape* pooled across elements (they all
carry the same statistic under the same design, so they differ in scale, not in
shape).

## Install

```bash
pip install -e .          # needs a C compiler; the C core sits at ../
pytest                    # 49 checks
```

## Status

Early. The transform, the tails and the permuted GLM are here and tested. The
design layer — turning a **BIDS Stats Model** into a design matrix, contrasts and,
crucially, an *exchangeability* structure — is not, and is the real work: BIDS-SM
specifies `X`, `Formula`, `Contrasts` and `GroupBy`, but has **no** concept of
exchangeability blocks or variance groups, so what may be permuted with what has
to be defined on top of it rather than read out of it.
