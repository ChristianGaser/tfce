# TFCE for SPM (MATLAB)

The SPM12 toolbox: non-parametric permutation inference with threshold-free cluster enhancement, for
3D volume and surface data.

TFCE combines focal effects of large height with broad effects of large extent, and needs **no
cluster-forming threshold** — the arbitrary choice that cluster-based inference forces on you, and
that the result can depend on heavily. It is also fairly robust to the non-stationarity that is
common in VBM data.

The toolbox runs on **any existing second-level SPM design**. Point it at an `SPM.mat` from a
parametric analysis you have already estimated, and it re-does the inference non-parametrically.

> There is also a **Python package** with the same C core and bit-identical results:
> [`pip install tfce`](../python/README.md). It is not an SPM toolbox — it is a library, and it drops
> into [nilearn](https://nilearn.github.io).

**Licence:** GPL-2.0-or-later (see [../COPYING](../COPYING)), because this toolbox is built on SPM and
carries code from SnPM and PALM. The C core in [../c](../c) and the Python package are BSD-3-Clause —
see [../LICENSE.md](../LICENSE.md).

---

## Contents

- [Download](#download)
- [Quick start](#quick-start)
- [Features](#features)
- [Settings](#settings)
- [Output files](#output-files)
- [Validation](#validation)
- [References](#references)

---

## Download

**[Download the latest release](https://github.com/ChristianGaser/tfce/releases/latest)** - grab the
`tfce_*.zip` asset and unpack it into your SPM `toolbox` directory.

[All releases](https://github.com/ChristianGaser/tfce/releases) ·
[Source](https://github.com/ChristianGaser/tfce)

The compiled mex files ship with the release. To rebuild them, run `compile` from the toolbox
directory.

There is also a **Python package**: `pip install tfce`. It is the same C core with a different
caller, so it gives bit-identical results, and it drops into
[nilearn](https://nilearn.github.io) as an exact replacement for its stepped TFCE. See
[../python/README.md](../python/README.md).

## Quick start

1. Estimate a second-level model in SPM as usual (volume or surface), and define your contrasts in
   the contrast manager.
2. Open the SPM batch editor and choose **TFCE → Estimate TFCE**.
3. Select the `SPM.mat`, pick the contrast(s), and run.

That is the whole workflow. Everything else has a sensible default.

Results are viewed with the toolbox's own results interface, which reads the `TFCE_log_p*` images the
same way SPM reads its own.

---

## Features

### Inference

- **Exact TFCE, with no step size.** The TFCE integral is
  `TFCE(v) = ∫ e(h)^E · h^H dh`. Most implementations approximate it by sampling `h` on a grid and
  summing, which introduces a discretisation error that depends on a step size `dh` you have to
  choose. This toolbox does not. It builds a **max-tree** (component tree) with union-find, and
  because the extent function `e(h)` is piecewise constant, integrates each piece in closed form. The
  result is the exact integral, at `O(N log N)`, with no precision parameter to get wrong.
- **Permutation inference throughout.** No parametric assumptions about the distribution of the
  statistic, and no random field theory.
- **FWE, FDR and uncorrected p-values**, for both the TFCE map and the raw t/F statistic.
- **t- and F-contrasts**, on volumes and on surfaces.

### Designs it understands

- **Any second-level SPM design**: one-sample, two-sample, paired, correlation, one-way and factorial
  ANOVA, repeated-measures ANOVA, ANCOVA, interactions.
- **Three ways of handling nuisance variables**: Draper-Stoneman, Smith (the default when nuisance
  variables are present), and Freedman-Lane. See [Winkler et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.01.060).
- **Exchangeability blocks** for longitudinal and repeated-measures designs, so that data are only
  ever permuted within the blocks they may be permuted within.
- **Sign-flipping** for one-sample designs, permutation for everything else - chosen automatically.
- **Voxel- and vertex-wise covariates**: a covariate given as one image per subject, so that the
  design matrix differs at every element.
- **Surfaces and TBSS**: surface meshes use the mesh adjacency rather than a voxel grid, and TBSS data
  use the 2D-optimised TFCE parameters (`E=1`, `H=2` instead of `E=0.5`, `H=2`).
- **Small volume correction** by supplying an additional mask.

### Fewer permutations

The number of permutations you need is set by the smallest p-value you want to resolve. Counting
exceedances cannot report a p-value below `1/n_perm` - that floor, not the statistic, is what forces a
permutation test to run many thousands of permutations. Three things remove it:

- **Gamma approximation of the maximum distribution** → FWE-corrected p-values. The distribution of
  the maximum statistic is fitted from its first three moments rather than counted, so corrected
  p-values are not floored at `1/n_perm` and converge with far fewer permutations.
- **Generalised Pareto tail approximation of the element-wise distributions** → *uncorrected* p-values,
  and with them FDR, which is computed from them. A GPD is fitted to the tail of each element's
  permutation distribution. Its *shape* is pooled across elements - they all carry the same statistic
  under the same design, so they differ in scale, not in shape, and fitting the shape per element
  spends a hundred tail values on a number they all share. From 1000 permutations this recovers
  p ≈ 1e-4 with the median right, where plain counting returns zero for 91% of elements and tells you
  nothing at all.
- **Sequential stopping** (optional, off by default). Stops as soon as it is certain that *nothing* in
  the image can become significant. The largest value in the image is watched, because no element
  rests on fewer exceedances than it does. A null image reaches that almost at once and stops at the
  floor of 500 permutations; an image with a real effect, or one sitting anywhere near alpha, never
  reaches it and runs the full count. On a null analysis this is a ~6.5x saving with the decision
  unchanged.

### Faster permutations

- **The permuted GLM never forms the permuted data.** Under Freedman-Lane the nuisance columns are
  part of the full design, so the residual-forming matrix already annihilates them and the
  `n_vox × n × n` multiply that builds the permuted dataset disappears entirely. The residual sum of
  squares needs no residual matrix either. The whole permutation collapses into one small projection.
- **Permutations are TFCE'd in parallel**, a block at a time, one permutation per thread. They are
  independent of one another, so this is where the parallelism belongs.
- **The max-tree sorts by radix, not by comparison.** Every value it sorts is a positive float, and
  for positive IEEE-754 floats the bit pattern read as an integer is monotonic - so the sort key needs
  no transformation and four counting passes replace `n log n` comparisons. This is the single
  largest cost of the transform.
- **The voxel-wise covariate path takes no pseudoinverse per voxel.** The design differs from element
  to element in only a few columns, so the normal equations are assembled for every element at once
  from a handful of matrix products and batch-solved. This is ~17x faster, and brings the voxel-wise
  path to roughly the cost of an ordinary design, where it used to be an order of magnitude slower.
- **Half the permutations, for free**, in balanced two-sample designs, where the mirrored permutation
  provably gives the negated statistic. Correctly refused when the design is not balanced.
- **Multiple `SPM.mat` files in parallel** across cores (`nproc`).

### Safety nets

A permutation test is only as good as its exchangeability assumptions, and a wrong design gives a
confidently wrong answer rather than an error. The toolbox therefore checks itself as it goes:

- **The permutations are verified before use**: data may only ever be exchanged inside their own
  exchangeability block. A wrong block structure invalidates the whole test and cannot be diagnosed
  from the results afterwards.
- **The permutation null is checked for width.** Under a valid permutation the uncorrected p-values are
  uniform on (0, 0.5] — one-sided, because they are conditioned on the sign of the observed effect —
  so ~5% of them land in the upper tail. True effects only produce *small* p-values and cannot inflate
  that upper tail, which makes it a clean check. A null that is too narrow makes the test
  anti-conservative, and is reported.
- **The parametric and the non-parametric statistic are compared.** A large discrepancy points to a
  misspecified permutation scheme and is reported at the command line.

---

## Settings

Everything below has a working default; in most analyses only the `SPM.mat` and the contrast need
choosing.

| Setting | Default | What it does |
| --- | --- | --- |
| **SPM.mat** | - | The second-level design to re-analyse. Contrasts must already be defined in the contrast manager. |
| **Additional mask** | the analysis mask | Restrict the analysis to a region of interest (small volume/surface correction). |
| **Results title** | automatic | Heading on the results page. |
| **Contrast index** | `Inf` | Index of the contrast(s), as numbered in the contrast manager. `Inf` selects them interactively. Define the contrasts themselves in the contrast manager, not here. |
| **Number of permutations** | 5000 | 5000–10000 gives reliable estimates. If the design admits fewer than this, all of them are used and the test is exact. Thanks to the tail approximations, the smallest reportable p-value is no longer tied to this number. |
| **Permutation method for nuisance variables** | Smith | Smith when nuisance variables are present, Draper-Stoneman when there are none. Only worth changing if a large discrepancy between the parametric and non-parametric statistic is reported. Freedman-Lane is available but should be used with care. |
| **Stop early if nothing can become significant** | no | Sequential stopping (see above). Off by default: it changes how many permutations an analysis runs, and a fixed count is easier to compare across analyses. Worth switching on when screening many contrasts most of which are expected to be null. |
| **TBSS data** | no | 2D-optimised TFCE parameters for TBSS/DTI data. |
| **Exchangeability blocks** | - | For longitudinal and repeated-measures designs. |
| **Number of processes** | number of cores | Parallelisation across several `SPM.mat` files. |

### TFCE parameters E and H

`E` (extent) and `H` (height) set how focal and broad effects are weighted against one another.
They are **not** user settings. They are fixed by the kind of data:

| Data | E | H |
| --- | --- | --- |
| 3D volume | 0.5 | 2 |
| Surface | 1 | 2 |
| TBSS (2D) | 1 | 2 |

These are the values [Smith and Nichols (2009)](https://doi.org/10.1016/j.neuroimage.2008.03.061)
established empirically as giving good power, and the surface/TBSS values follow the same reasoning
for data of lower effective dimension. Earlier versions of the toolbox exposed `E` as a setting; it no
longer is, because changing it changes what "significant" means in a way that is not comparable
across studies and is easy to misuse.

---

## Output files

One set per contrast `nnnn`. All p-value images are `-log10(p)`, in SPM's usual convention, and carry
the sign of the effect.

| File | Contents |
| --- | --- |
| `T_nnnn` / `F_nnnn` | The unpermuted statistic. |
| `TFCE_nnnn` | The TFCE of that statistic. |
| `TFCE_log_p_nnnn` | Uncorrected p, TFCE. |
| `TFCE_log_pFWE_nnnn` | FWE-corrected p, TFCE. |
| `TFCE_log_pFDR_nnnn` | FDR-corrected p, TFCE. |
| `T_log_p_nnnn`, `T_log_pFWE_nnnn`, `T_log_pFDR_nnnn` | The same three, for the raw statistic rather than its TFCE. |
| `T_nnnn.txt` | The number of permutations actually run. |

For most purposes `TFCE_log_pFWE_nnnn` is the map to look at.

---

## Validation

The toolbox ships a validation suite that establishes, rather than assumes, that the things above are
true. From the `matlab/` directory:

```matlab
addpath /path/to/spm12
spm('defaults','fmri')
compile              % builds the mex files against the shared C core in ../c
validation/run_all
```

It runs 94 checks and prints a `[PASS]`/`[FAIL]` line for each. Among them: that the max-tree really
is the exact TFCE integral (an independent dh-stepping implementation converges onto it at first
order); that the Gamma and Pareto approximations are calibrated against counting, and that the Pareto
fit is never *less* accurate than the counting it replaces; that all three nuisance methods control
the false-positive rate, including when the nuisance variable is correlated with the effect of
interest; that the accelerated GLM returns exactly the statistic the unaccelerated one returns; and
that sequential stopping reaches the same decision as running every permutation.

See [validation/README.md](validation/README.md) for what each check establishes and why it is
written the way it is.

---

## References

**TFCE**
Smith SM, Nichols TE (2009). *Threshold-free cluster enhancement: addressing problems of smoothing,
threshold dependence and localisation in cluster inference.* NeuroImage 44:83–98.
[doi:10.1016/j.neuroimage.2008.03.061](https://doi.org/10.1016/j.neuroimage.2008.03.061)

**Permutation inference and nuisance variables**
Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE (2014). *Permutation inference for the
general linear model.* NeuroImage 92:381–397.
[doi:10.1016/j.neuroimage.2014.01.060](https://doi.org/10.1016/j.neuroimage.2014.01.060)

**Tail approximations and faster permutation inference**
Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM (2016). *Faster permutation inference in brain
imaging.* NeuroImage 141:502–516.
[doi:10.1016/j.neuroimage.2016.05.068](https://doi.org/10.1016/j.neuroimage.2016.05.068)

**Sequential Monte Carlo tests**
Besag J, Clifford P (1991). *Sequential Monte Carlo p-values.* Biometrika 78:301–304.
Gandy A (2009). *Sequential implementation of Monte Carlo tests with uniformly bounded resampling
risk.* JASA 104:1504–1511.
