# Changelog

All notable changes to the `tfce` Python package. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

This is the Python package only. The MATLAB/SPM toolbox in the same repository is versioned
separately; its history is in [`matlab/CHANGES.txt`](../matlab/CHANGES.txt).

## [Unreleased]

## [0.1.0] - 2026-07-14

First release. It is a new package, but not new code: the C core it compiles has been in use in
the MATLAB/SPM toolbox for years, and the test suite asserts that this package reproduces that
toolbox **bit for bit**.

### Added

- `tfce.tfce()` - exact threshold-free cluster enhancement, for volumes (6-, 18- or
  26-connectivity) and for surface meshes. The integral is evaluated in closed form by a max-tree,
  so there is **no step size**: nothing to tune, and no discretisation error. Accepts a stack of
  maps and transforms them across threads, releasing the GIL.
- `tfce.adjacency_from_faces()` - CSR adjacency of a triangulated surface, with every face index
  checked to name a vertex that exists.
- `tfce.glm.PermutedGLM` - a GLM that re-fits under permutations and sign-flips **without ever
  forming the permuted data**, along with `partition_design()` for the Guttman partition into
  effects of interest and nuisance.
- `tfce.tails` - `gamma_pvalue()` for FWE-corrected p-values from the distribution of the maximum,
  and `pareto_pvalue()` for uncorrected p-values, both resolving p below the `1/n_perm` floor of
  counting. Roughly a fifth of the permutations for the same accuracy.
- `tfce.nilearn_compat.calculate_tfce()` - a drop-in for nilearn's internal TFCE with the same
  signature. `dh` is accepted and ignored, because there is no step size.
- Wheels for Linux (x86_64, aarch64), macOS (Intel, Apple silicon) and Windows, on Python
  3.9–3.13. No compiler required.
- Documentation in [`docs/`](docs/): installation, usage, permutation testing, theory, the nilearn
  comparison, the API reference, benchmarks and validation.
- 53 tests, including bit-exact parity against reference output from the MATLAB toolbox.

### Notes

- **License**: BSD-3-Clause. This package and the C core it compiles contain no SPM, SnPM or PALM
  code and do not inherit their GPL. The MATLAB toolbox in the same repository remains
  GPL-2.0-or-later.
- **Accuracy**: against a stepped TFCE, the error halves each time the step count doubles (1.7% at
  nilearn's default of 100 steps, 0.8% at 200) - a discretisation error converging onto the exact
  transform this package computes.
- **Speed**: 14× a stepped implementation on one map, 72× on a block of 16 permutations. The
  transform is memory-bound, not compute-bound, so it plateaus at about 2× from six threads and
  gets *worse* beyond eight. `n_jobs=-1` sits on that plateau.

[Unreleased]: https://github.com/ChristianGaser/tfce/compare/v1.3...HEAD
[0.1.0]: https://github.com/ChristianGaser/tfce/releases/tag/v1.3
