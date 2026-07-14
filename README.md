# TFCE: Threshold-Free Cluster Enhancement

Non-parametric permutation inference with threshold-free cluster enhancement, for 3D volume and
surface data - as an **SPM toolbox** and as a **Python package**.

TFCE combines focal effects of large height with broad effects of large extent, and needs **no
cluster-forming threshold** - the arbitrary choice that cluster-based inference forces on you, and
that the result can depend on heavily. It is also fairly robust to the non-stationarity that is
common in VBM data.

The TFCE of an element is an integral,

```text
TFCE(v) = ∫ e_v(h)^E · h^H dh
```

over the extent `e_v(h)` of the cluster containing `v` at height `h`. Implementations normally
approximate it by stepping `h` over a grid and summing. **This one does not.** It builds the
**max-tree** (the component tree) with union-find, and because the extent function is piecewise
constant, integrates each piece in closed form. There is no step size `dh`, no discretisation error,
and no accuracy parameter to get wrong.

---

## Two implementations, one core

Both call the **same C**, so they give bit-identical results. There is exactly one implementation of
TFCE in this repository to get right.

| | | |
| --- | --- | --- |
| **[MATLAB / SPM toolbox](matlab/README.md)** | An SPM12 toolbox. Point it at an `SPM.mat` from any second-level design you have already estimated, and it re-does the inference non-parametrically. | [**Download**](https://github.com/ChristianGaser/tfce/releases/latest) · [docs](matlab/README.md) |
| **[Python package](python/README.md)** | A library, not a pipeline. Arrays in, arrays out. Drops into [nilearn](https://nilearn.github.io) as an exact replacement for its stepped TFCE. | `pip install tfce` · [docs](python/docs/) |

**Which one do I want?**

- You already work in **SPM or CAT12**, and you have an `SPM.mat` → the **MATLAB toolbox**.
- You work in **Python**, or want TFCE inside **nilearn / nipreps / fitlins**, or just want the
  transform as a function → the **Python package**.

---

## Repository layout

| | |
| --- | --- |
| [`c/`](c/) | The TFCE core, plain C. Shared by both bindings and owned by neither. No dependency beyond libc and pthreads. |
| [`matlab/`](matlab/) | The SPM toolbox: m-files, mex glue, prebuilt binaries, the [validation suite](matlab/validation/) and the [help pages](matlab/html/). |
| [`python/`](python/) | The pip package: a Cython binding to the same core, plus the permutation machinery in numpy. |

An installed SPM toolbox is a single flat folder, so `make zip` collapses `matlab/` and the headers
from `c/` into one. `compile.m` looks for the core beside itself first and in `../c` second, which is
what lets it build from either layout.

---

## Licensing

The two halves have different licences, and the split follows exactly what each is **derived from**:

| | Licence | Why |
| --- | --- | --- |
| [`c/`](c/), [`python/`](python/) | **BSD-3-Clause** | Original work. No SPM, SnPM or PALM code - nothing obliges it to be GPL, and permissive projects (nilearn is BSD-3, nipreps is Apache-2.0) can depend on it freely. |
| [`matlab/`](matlab/) | **GPL-2.0-or-later** | Genuinely derived from GPL code: `snpm_P_FDR.m` from SnPM, `cat_spm_results_ui.m` from SPM, the `palm_*` subfunctions from PALM, and every `.m` calls SPM. |

BSD-3-Clause is GPLv2-compatible, so the same C core is compiled into the GPL mex-files and the BSD
Python wheel alike. See **[LICENSE.md](LICENSE.md)**.

---

## Validation

Both implementations ship a validation suite, because a permutation test that is subtly wrong gives a
confidently wrong answer rather than an error.

- **MATLAB** - 105 checks: [`matlab/validation/`](matlab/validation/)
- **Python** - 53 checks: [`python/tests/`](python/tests/)

Among them: that the max-tree really is the *exact* TFCE integral (an independent dh-stepping
implementation must converge onto it at first order); that the Gamma and Pareto tail approximations
are calibrated against counting; that all three nuisance methods control the false-positive rate; and
that the Python binding is **bit-identical** to the MATLAB mex.

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

---

Developed by Christian Gaser, Structural Brain Mapping Group, Departments of Psychiatry and
Neurology, Jena University Hospital. <https://neuro-jena.github.io>
