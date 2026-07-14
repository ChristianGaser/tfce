# TFCE — repository guide

Threshold-free cluster enhancement and the permutation inference around it, in two
implementations that share one C core.

The distinguishing claim, and the thing to protect in any change: **the TFCE integral
is evaluated exactly.** Every other implementation (fslmaths, nilearn, PALM) sums the
integrand over a ladder of thresholds and inherits a discretisation error — 1.7% at
nilearn's default. Here a max-tree over a union-find gives the integral in closed form.
There is no `dh`, no step count, nothing to tune.

---

## Layout

| | license | what it is |
| --- | --- | --- |
| `c/` | **BSD-3** | The core. Header-only, plain C, no MATLAB and no Python in it. |
| `matlab/` | **GPL-2.0-or-later** | The SPM toolbox: m-files, mex glue, prebuilt binaries. |
| `python/` | **BSD-3** | The pip package `tfce`. |
| `Makefile` | | Releases the **MATLAB toolbox only**. |
| `.github/workflows/` | | `compile_mex.yml` (MATLAB, manual) and `python_wheels.yml` (PyPI). |

`c/` is owned by neither binding. Both compile the same headers, so a change there
changes both — and the Python tests assert bit-identical agreement with MATLAB
(`python/tests/test_matlab_parity.py`), which is what keeps that honest.

### `c/`
- `tfce_maxtree.h` — the algorithm. Union-find over voxels/vertices sorted descending.
  `tfce_val` is `float`. Sorting is a **radix sort on the raw bit patterns**: every value
  is a positive IEEE-754 float, and those are monotonic read as `uint32`, so no key
  transform is needed. That sort, not single precision, is where the speed came from.
- `tfce_batch.h` — threaded batch driver, one permutation per thread.
- `tfce_threads.h` — pthreads / Win32 shim.
- `tfce_capi.h` / `.c` — flat C API (`tfce_volume`, `tfce_mesh`, `*_batch`,
  `tfce_adjacency_from_faces`, `tfce_free`). This is what Cython calls.

### `matlab/`
- `tfce_estimate_stat.m` (~4100 lines) — the whole permutation engine. Everything
  interesting lives in its subfunctions; see **Config flags** below.
- `CHANGES.txt` — the toolbox changelog, newest first, `Added:`/`Changed:`/`Fixed:`
  style, back to the 2010 SVN import. Add a `TFCE x.y | rNNN | date` block on top for
  each release. Shipped in the toolbox and uploaded beside the zip by `make scp`.
- `tfceMex_maxtree.c`, `tfceMex_maxtree_batch.c`, `tfceMex_resss.c` — mex glue only.
- `tbx_cfg_tfce.m` — the SPM batch GUI. New user-facing options go here **and** in
  `tfce_estimate_stat.m`.
- `compile.m` — finds the core beside itself (installed flat toolbox) *or* at `../c`
  (repo). Both layouts must keep working.
- `validation/` — **105 checks**, `run_all.m`.
- `html/` — generated from `html/tfce.txt` by `make doc`. Never edit `tfce.html`.

### `python/`
- `src/tfce/core.py` — `tfce()`, `adjacency_from_faces()`, `TfceError`.
- `src/tfce/glm.py` — `PermutedGLM`, `partition_design`.
- `src/tfce/tails.py` — `gamma_pvalue`, `pareto_pvalue`, `gpd_fit_pwm`, `gpd_sf`.
- `src/tfce/nilearn_compat.py` — drop-in for `nilearn.mass_univariate._utils`.
- `src/tfce/_maxtree.pyx` — Cython binding. Releases the GIL around every core call.
- `tests/` — **53 tests**.
- `docs/` — 8 pages, the real documentation. `README.md` is the PyPI front page.

---

## Building and testing

```bash
# MATLAB
cd matlab && matlab -batch "compile(1); cd validation; run_all"     # 105 checks

# Python
cd python && pip install -e '.[dev]' && pytest                      # 53 tests

# MATLAB release
make zip           # or: make install / make doc / make update / make scp
```

`setup.py` **vendors `c/` into `src/tfce/_c/` at build time** — setuptools rejects source
paths outside the package root. So `python/` alone is not buildable; the whole repo must
be checked out. That is why `python_wheels.yml` checks out everything, not just `python/`.

Python deps are deliberately thin: `numpy`, `scipy`. `nibabel` is `[io]`; `nilearn` is
`[dev]` only, so the core can be depended on *by* nilearn without a dependency cycle.
`test_nilearn_compat.py` uses `importorskip` and correctly skips in a bare install.

---

## Config flags in `tfce_estimate_stat.m`

Near the top of the file. These are the knobs that change results, not just speed:

| flag | default | |
| --- | --- | --- |
| `use_gamma_tail_approximation` | `true` | Gamma fit to the max distribution → FWE p below `1/n_perm`. |
| `use_pareto_tail_approximation` | `true` | Generalised Pareto fit to the upper tail → uncorrected p below the counting floor. ~5× fewer permutations. |
| `n_tail` | `100` | Tail values kept per element. |
| `n_exc_min` | `25` | Exceeded this often → keep the count, the fit adds nothing. |
| `tail_memory_limit` | `4` GB | Above this the Pareto tail buffer is refused. |
| `tfce_block_size` | `0` | Threads for the batch mex. `0` = auto. |
| `use_sequential_stopping` | **off**, from the GUI | Besag–Clifford + Gandy margin. |
| `n_exceed_stop`, `n_sigma_stop`, `n_perm_min` | `20`, `3`, `500` | |

---

## Things that have bitten us, and must not be undone

These are all fixed. Each one was subtle, each one is easy to reintroduce.

**Surface faces are not doubles.** GIFTI `g.faces` is `int32`. The mex once read it with
`mxGetPr()`, reinterpreting int32 as double → garbage vertex indices → heap corruption on
*every* surface analysis. `faces_to_int()` handles double/single/int16/uint16/int32/uint32,
and `build_adjacency()` validates every index is in `[1, N]`. Keep both.

**Mesh degree is unbounded.** Neighbours were copied into a fixed `nbuf[26]`; a vertex with
more than 26 neighbours overflowed the stack. `neigh_get()` now returns a *pointer into the
CSR row* for meshes. Do not reintroduce a fixed buffer.

**The Pareto finite-endpoint trap.** A PWM fit can land on `k > 0`, giving the distribution
a finite upper endpoint — and if the observation sits beyond it, `p = 0`. This hit 14–20% of
elements. Fits the observation contradicts are rejected and fall back to the exponential
(`k = 0`).

**Sequential stopping needs Gandy's margin.** Plain Besag–Clifford exceedance counting
disagreed with the full run 40% of the time at `p ≈ alpha`. The 3-sigma margin
(`n_sigma_stop`) is not optional.

**The permuted GLM never forms the permuted data.** Under Freedman–Lane `R·Rz = R`, so `Rz`
drops out, and `‖Pset·y‖²` is permutation-invariant. `calc_GLM_fast` exploits this — but
only when `W = I`, so there is a guard that checks `R*Rz ≈ R` and falls back otherwise.
Do not remove the guard.

**Python permutations are inverted relative to the index.** `PermutedGLM.fit(perm)` uses
`np.argsort(perm)` because the permuted data is `Y[:, perm]` and the design must move the
*other* way. Tests catch this; believe them.

**Normal equations square the condition number.** `calc_GLM_voxelwise` accumulates in double
and tests rank with a Schur complement — an LU solve of a singular system returns huge
*finite* values, so `~isfinite` does not catch degeneracy.

---

## Conventions

- **Two licenses, on purpose.** `c/` and `python/` are BSD-3 so permissive projects
  (nilearn is BSD-3, nipreps is Apache-2.0) can take them without friction; BSD-3 is also
  GPLv2-compatible, so the MATLAB toolbox stays GPL. Nothing in `c/` or `python/` is derived
  from SPM/SnPM/PALM. **Keep it that way** — copying SPM code into either would relicense it.
  Every BSD file carries an SPDX header.
- **Never hand-edit generated files**: `matlab/Contents.m` (from `Contents_info.txt`),
  `matlab/INSTALL.txt` (from `INSTALL_info.txt`), `matlab/html/tfce.html` (from `tfce.txt`).
  `make update` regenerates them, so version/revision stamps in them are stale between
  releases — expected, not a bug.
- **Check counts appear in prose** in `README.md`, `python/docs/installation.md` and
  `python/docs/validation.md`. Adding a test means updating them; they have drifted before.
- The root `README.md` is a **hub**. It points at the two implementations and does not
  document either.
- `python/docs/benchmarks.md` numbers are **measured, not estimated**, and the file carries
  the scripts that reproduce them. Do not edit a number without rerunning it.
- Don't expect thread scaling. The max-tree is memory-bound: it plateaus at ~2× from six
  threads and gets *worse* beyond eight.

---

## Loose ends

- `methods_validation_draft.md` (untracked), `SPM-Designs-Permutation.rtf` and `spm/`
  (tracked) are working notes at the repo root, not shipped artefacts.
- Not done: non-macOS wheel builds have never actually run; PyPI trusted publishing needs
  the `pypi` environment configured on GitHub; the nilearn PR; the BIDS Stats Models
  exchangeability layer.
