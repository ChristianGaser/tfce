# Installation

```bash
pip install tfce
```

That is all. Wheels are built for Linux (x86-64, aarch64), macOS (Intel and Apple silicon) and
Windows, for CPython 3.9–3.13, so **no compiler is needed**.

## What it pulls in

Only **numpy** and **scipy**. Nothing else.

That is deliberate. The core is meant to be depended on — by nilearn, by a pipeline, by a one-off
script — and a library that drags a tree of dependencies behind it is a library people work around
instead of using.

## Optional extras

```bash
pip install "tfce[io]"     # + nibabel, for reading NIfTI and GIFTI
pip install "tfce[test]"   # + pytest
pip install "tfce[dev]"    # + nibabel, nilearn, cython, pytest
```

`nilearn` is **not** a dependency, even for `tfce.nilearn_compat`. That module imports it lazily, so
the package works without it — and so nilearn could one day depend on `tfce` without a cycle.

## Building from source

The source distribution carries the C core, so it builds anywhere with a C compiler:

```bash
pip install --no-binary tfce tfce
```

From a git checkout:

```bash
git clone https://github.com/ChristianGaser/tfce
cd tfce/python
pip install -e ".[dev]"
pytest                       # 49 checks
```

The C lives in [`c/`](../../c) at the repository root, shared with the MATLAB toolbox, and
`setup.py` vendors it into the package at build time. So an edit to the C is picked up by the next
build with nothing to remember, and the sdist is still self-contained.

## Checking it worked

```python
import numpy as np
import tfce

x = np.zeros((9, 9, 9))
x[4, 4, 4] = 3.0
x[4, 5, 4] = 2.0

print(tfce.__version__)
print(tfce.tfce(x).max())     # ~10.1
```

## Requirements

| | |
| --- | --- |
| Python | ≥ 3.9 |
| numpy | ≥ 1.22 |
| scipy | ≥ 1.8 |
| nibabel *(optional)* | ≥ 4.0 |

## Threads

The batched transform runs one permutation per thread and **releases the GIL** for the whole call, so
`n_jobs` does what it says even inside a thread pool. It uses pthreads (Win32 threads on Windows)
directly — not OpenMP, not joblib — so there is no thread-pool interaction to reason about and no
environment variable to set.
