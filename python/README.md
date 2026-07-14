# tfce

**Exact threshold-free cluster enhancement, and the permutation inference around it.**

[![PyPI](https://img.shields.io/pypi/v/tfce.svg)](https://pypi.org/project/tfce/)
[![Python](https://img.shields.io/pypi/pyversions/tfce.svg)](https://pypi.org/project/tfce/)
[![License](https://img.shields.io/pypi/l/tfce.svg)](https://github.com/ChristianGaser/tfce/blob/master/python/LICENSE)

TFCE combines focal effects of large height with broad effects of large extent, and needs **no
cluster-forming threshold** - the arbitrary choice that cluster-based inference forces on you, and
that the result can depend on heavily.

```bash
pip install tfce
```

Wheels for Linux, macOS and Windows. No compiler needed.

---

## Exact, not stepped

The TFCE of an element is an integral:

```text
TFCE(v) = ∫ e_v(h)^E · h^H dh
```

over the extent `e_v(h)` of the cluster containing `v` at height `h`. Implementations normally
approximate it by stepping `h` over a grid and summing. That costs a step size `dh`, a discretisation
error that depends on it, and an accuracy parameter you have to guess.

This one doesn't. It builds the **max-tree** (the component tree) with union-find, and because
`e_v(h)` is piecewise constant, integrates each piece in closed form. The answer is the integral, not
a sample of it. **There is nothing to tune.**

The difference is not academic. Against a stepped implementation on a 60×72×60 volume:

| steps | error vs exact |
| --- | --- |
| 50 | 3.4% |
| 100 | **1.7%** |
| 200 | 0.8% |
| 400 | 0.4% |

The error halves every time the steps double - first order, exactly as a step-size approximation must
behave. The exact transform has no such term.

## Quick start

```python
import numpy as np
import tfce

# a volume: (nx, ny, nz)
t = tfce.tfce(stat_map, E=0.5, H=2.0)

# a surface: faces are 1-based, as GIFTI stores them
adj = tfce.adjacency_from_faces(faces, n_vertices)
t = tfce.tfce(surf_map, adjacency=adj, E=1.0, H=2.0)

# a block of permutations, one per thread
#   (nx, ny, nz, n_perm) -> (nx, ny, nz, n_perm)
t = tfce.tfce(perms, E=0.5, H=2.0, n_jobs=-1)
```

Volumes take `connectivity=6 | 18 | 26` (26 is the default, and what fslmaths and the MATLAB toolbox
use). Surfaces take the mesh, so the neighbourhood is whatever the mesh says it is.

Everything is arrays in, arrays out. No image objects, no file I/O, no design parsing - those belong
in a layer above, so the core can be dropped anywhere.

## Speed

Permutations are independent, so they are TFCE'd a block at a time, one per thread, with the GIL
released. On the same 60×72×60 volume:

| | one map | 16 permutations |
| --- | --- | --- |
| stepped (100 steps) | 0.38 s | 9.1 s - 569 ms/perm |
| **tfce** | **0.03 s** (14×) | **0.13 s** - 8 ms/perm (**72×**) |

## Using it with nilearn

nilearn has `permuted_ols(..., tfce=True)`, and its TFCE is a stepped approximation. `tfce` ships a
drop-in with the same signature:

```python
from nilearn.mass_univariate import _utils
import tfce.nilearn_compat as tc

_utils.calculate_tfce = tc.calculate_tfce   # now exact, and much faster
```

It reads the neighbourhood out of the `bin_struct` it is handed, so it reproduces whichever
connectivity the caller meant.

## Fewer permutations

Counting exceedances cannot report a p-value below `1/n_perm`. That floor - not the statistic - is
what forces a permutation test to run many thousands of permutations, and it caps FDR too, since FDR
is computed from the uncorrected p-values. `tfce.tails` removes it:

- **Gamma** fit to the null of the *maximum*, for FWE-corrected p-values.
- **Generalised Pareto** fit to each element's own tail, for the uncorrected ones - with the *shape*
  pooled across elements, since they all carry the same statistic under the same design and so differ
  in scale, not in shape.

From 1000 permutations this recovers `p ≈ 1e-4` with the median right, where plain counting returns
zero for 91% of elements and tells you nothing at all.

## What's in the box

| | |
| --- | --- |
| `tfce.core` | the exact transform - arrays in, arrays out |
| `tfce.tails` | Gamma and Generalised Pareto tail approximations |
| `tfce.glm` | the permuted GLM, which never forms the permuted data |
| `tfce.nilearn_compat` | drop-in for nilearn's `calculate_tfce` |

Dependencies are deliberately thin: **numpy and scipy**, nothing else. `nibabel` is an optional extra,
wanted only by I/O.

## Trust

This is not a reimplementation. It is the same C core that the
[MATLAB TFCE toolbox](https://github.com/ChristianGaser/tfce) has been running for years, with a
Cython binding instead of a MEX one - so the two give **bit-identical** results, and the test suite
holds them to it. It ships a validation suite of its own: the exactness of the max-tree is established
against an *independent* stepped implementation, which must converge onto it at first order.

## Citing

If you use this, please cite the method:

> Smith SM, Nichols TE (2009). *Threshold-free cluster enhancement: addressing problems of smoothing,
> threshold dependence and localisation in cluster inference.* NeuroImage 44:83–98.
> [doi:10.1016/j.neuroimage.2008.03.061](https://doi.org/10.1016/j.neuroimage.2008.03.061)

The tail approximations are from:

> Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM (2016). *Faster permutation inference in
> brain imaging.* NeuroImage 141:502–516.
> [doi:10.1016/j.neuroimage.2016.05.068](https://doi.org/10.1016/j.neuroimage.2016.05.068)

## Status

Early, but the parts that are here are tested and are the parts that matter. The transform, the tail
approximations and the permuted GLM are done. Not yet done: the **design layer** - turning a
[BIDS Stats Model](https://bids-standard.github.io/stats-models/) into a design matrix, contrasts and,
crucially, an *exchangeability* structure. That last one is the real work, because BIDS-SM specifies
`X`, `Formula`, `Contrasts` and `GroupBy` but has **no** concept of exchangeability blocks, so what may
be permuted with what has to be defined on top of it rather than read out of it.

## Documentation

| | |
| --- | --- |
| [Installation](https://github.com/ChristianGaser/tfce/blob/master/python/docs/installation.md) | Wheels, building from source, extras |
| [Usage](https://github.com/ChristianGaser/tfce/blob/master/python/docs/usage.md) | Volumes, surfaces, batches, connectivity, `E` and `H` |
| [A permutation test, end to end](https://github.com/ChristianGaser/tfce/blob/master/python/docs/permutation.md) | The whole thing, as a worked example you can run |
| [Theory](https://github.com/ChristianGaser/tfce/blob/master/python/docs/theory.md) | Why the integral is exact, and why you need fewer permutations |
| [Using it with nilearn](https://github.com/ChristianGaser/tfce/blob/master/python/docs/nilearn.md) | The drop-in - and a bug in nilearn's stepped TFCE |
| [API reference](https://github.com/ChristianGaser/tfce/blob/master/python/docs/api.md) | Every public function |
| [Benchmarks](https://github.com/ChristianGaser/tfce/blob/master/python/docs/benchmarks.md) | Accuracy and speed, with scripts to reproduce |
| [Validation](https://github.com/ChristianGaser/tfce/blob/master/python/docs/validation.md) | What the test suite establishes, and how |

## Links

- **Source, issues, MATLAB toolbox:** <https://github.com/ChristianGaser/tfce>
- **Licence:** BSD-3-Clause. Permissive on purpose: `nilearn` is BSD-3 and
  `nipreps` is Apache-2.0, and neither can take on a GPL dependency. Nothing here
  is derived from SPM, so nothing here has to be GPL. (The MATLAB/SPM toolbox in
  the same repository *is* GPL, because it genuinely is derived from SPM - see
  [LICENSE.md](https://github.com/ChristianGaser/tfce/blob/master/LICENSE.md).)

Developed by Christian Gaser, Structural Brain Mapping Group, Jena University Hospital.
