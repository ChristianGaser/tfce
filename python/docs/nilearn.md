# Using it with nilearn

nilearn already has TFCE: `permuted_ols(..., tfce=True)`, backed by
`nilearn.mass_univariate._utils.calculate_tfce`. That implementation **steps** the height over a grid.
This package gives you the same function with the same signature, computed exactly.

```python
from nilearn.mass_univariate import _utils
import tfce.nilearn_compat as tc

_utils.calculate_tfce = tc.calculate_tfce
```

That is the whole integration. Everything downstream - `permuted_ols`, the maskers, the outputs -
carries on unchanged, and now gets an exact transform, 14× faster on one map and 72× on a block of
permutations.

`nilearn` is **not** a dependency of this package. `tfce.nilearn_compat` imports it lazily, so the
package installs and works without it.

---

## Mind the connectivity

nilearn builds its neighbourhood with `generate_binary_structure(3, 1)`, which is **6-connectivity** -
faces only. The MATLAB TFCE toolbox and `fslmaths` use **26**.

The drop-in reads the neighbourhood out of the `bin_struct` it is handed, so it reproduces whichever
one the caller meant. But if you are moving between the two ecosystems, know that **their default
answers differ**, and not by a little: a looser neighbourhood merges clusters, so extents grow and
TFCE values grow with them.

```python
tfce.tfce(x, connectivity=6)    # what nilearn does
tfce.tfce(x, connectivity=26)   # what fslmaths and the MATLAB toolbox do
```

---

## A bug worth knowing about

nilearn's stepped implementation has a defect that is easy to miss, and it is the strongest reason to
replace it.

Look at `_return_score_threshs`: with `dh="auto"` it takes **100 steps between 0 and `max(|map|)`** -
so `dh = max(|map|)/100`, recomputed **for every map**. And `calculate_tfce` deliberately does *not*
multiply the partial scores by `dh` (it follows `fslmaths` here, and says so).

Each of those choices is defensible. Together they are not:

```text
TFCE_nilearn(map) ≈ (100 / max|map|) · TFCE_true(map)
```

Omitting `dh` is only a harmless constant **if `dh` is a constant**. It isn't - it depends on the
map's own maximum, and in a permutation test every permutation has a different maximum. So every
permutation's TFCE is rescaled by a different factor, and the permutation null becomes a mixture of
differently-scaled maxima.

The test stays *valid* - the statistic is still a deterministic function of the map, so exchangeability
is untouched - but it is no longer TFCE, and the rescaling is noise.

**It costs power.** One-sample sign-flip test, planted blob, FWE on the max statistic, same data and
same sign-flips for both:

| | false positives | power |
| --- | --- | --- |
| exact | 0.013 | **0.920** |
| nilearn (stepped) | 0.007 | 0.853 |

Both control errors; the exact transform detects the effect about **7 points more often**. (150
simulations, so that gap is roughly 2 standard errors - suggestive rather than settled.)

---

## Accuracy of the stepped transform

Independently of the rescaling, the step size costs accuracy. On a 60×72×60 volume, against the exact
integral:

| steps | error |
| --- | --- |
| 50 | 3.4% |
| **100** (nilearn's default) | **1.7%** |
| 200 | 0.8% |
| 400 | 0.4% |

The error halves every time the steps double - first order, exactly as a Riemann sum must behave. The
exact transform has no such term at all.

---

## Could this go into nilearn itself?

Not as a compiled extension: **nilearn ships no compiled code**, by policy, and a Cython dependency
would be rejected on those grounds before anyone read it.

Two routes exist:

1. **As an optional backend** - nilearn uses `tfce` if it is installed, and falls back to its own
   stepped implementation otherwise. Small, and needs no compiled code in nilearn.
2. **A pure-Python exact max-tree**, which turns out to be possible using only what nilearn already
   depends on: `scipy.sparse.csgraph.minimum_spanning_tree` for the merge order, a single-linkage
   dendrogram, and a vectorised pointer-doubling accumulation. Exact, and roughly at parity with
   nilearn's current speed - but it needs the C to be genuinely fast.

Either way, the licence is not an obstacle: this package is **BSD-3-Clause**, same as nilearn.
