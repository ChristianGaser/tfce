# Benchmarks

All numbers below are measured, not estimated, on an Apple M2 (8 cores). Reproduce them with the
scripts at the bottom.

---

## Accuracy

Against a stepped implementation (nilearn's), on a 60×72×60 volume:

| steps | max error vs exact |
| --- | --- |
| 50 | 3.4% |
| **100** — nilearn's default | **1.7%** |
| 200 | 0.8% |
| 400 | 0.4% |
| 800 | 0.2% |

The error **halves every time the steps double** — first order, exactly as a Riemann sum must behave.
That is the signature of a discretisation error, and it is what proves the max-tree is not merely
another approximation: the stepped transform converges *onto* it.

The exact transform has no such term. There is no step size to increase.

---

## Speed

Same volume, one map:

| | time |
| --- | --- |
| stepped (100 steps) | 0.38 s |
| **exact max-tree** | **0.03 s** — **14×** |

A block of 16 permutations, which is what a permutation test actually does:

| | total | per permutation |
| --- | --- | --- |
| stepped (100 steps) | 9.1 s | 569 ms |
| **exact, threaded** | **0.13 s** | **8 ms** — **72×** |

The gap widens with the batch because permutations are independent, so they are handed out one per
thread and never lock.

### Do not expect it to scale with core count

The max-tree is bound by **memory**, not by arithmetic. On an 8-core machine:

| threads | s / permutation | speed-up |
| --- | --- | --- |
| 1 | 0.128 | 1.0× |
| 2 | 0.121 | 1.1× |
| 4 | 0.070 | 1.8× |
| **6** | **0.063** | **2.0×** |
| 8 | 0.063 | 2.0× |
| 12 | 0.077 | 1.7× |
| 16 | 0.075 | 1.7× |

It plateaus at about **2×** from six threads onwards, and gets *worse* beyond eight — eight threads
each streaming ~7 MB in and 7 MB out saturate the memory system long before they saturate the cores.
`n_jobs=-1` sits on that plateau, which is why it is the right default and why asking for more is
counter-productive.

---

## Where the speed comes from

Not from where you would guess. Measured, step by step, on a 91×109×91 volume:

| | s / permutation |
| --- | --- |
| original: `double`, `qsort` | 0.128 |
| `single`, `qsort` | 0.116 |
| **`single` + radix sort** | **0.081** |
| batched, 8 threads | **0.051** |

**Single precision on its own bought almost nothing** (0.128 → 0.116). The transform is bound by
*random access* — the sort and the union-find pointer chasing — not by streaming bandwidth, so halving
the size of arrays you are already missing cache on does not help.

What it *enabled* was the sort. Every value the max-tree sorts is a **positive float**, and for
positive IEEE-754 floats the bit pattern read as a `uint32` increases monotonically with the value. So
the sort key needs no transformation at all, and four counting passes replace `n log n` comparisons
through a function pointer. **That** is where the time was.

---

## Reproducing

```python
# accuracy and speed against nilearn
import time
import numpy as np
from scipy.ndimage import gaussian_filter, generate_binary_structure
from nilearn.mass_univariate._utils import calculate_tfce as stepped
from tfce.nilearn_compat import calculate_tfce as exact

rng = np.random.default_rng(0)
arr = gaussian_filter(rng.standard_normal((60, 72, 60)), 2.0)
arr = (arr / arr.std())[..., None]
bs = generate_binary_structure(3, 1)          # 6-connectivity, as nilearn uses

ref = exact(arr, bs, E=0.5, H=2)

for n_steps in (50, 100, 200, 400):
    dh = np.abs(arr).max() / n_steps
    approx = stepped(arr, bs, E=0.5, H=2, dh=dh)
    # nilearn follows fslmaths and does not multiply by dh, so rescale to compare
    scale = np.sum(ref * approx) / np.sum(approx * approx)
    err = np.abs(scale * approx - ref).max() / np.abs(ref).max()
    print(f"{n_steps:4d} steps: {err:.4f}")

t0 = time.perf_counter(); stepped(arr, bs, E=0.5, H=2); t_s = time.perf_counter() - t0
t0 = time.perf_counter(); exact(arr, bs, E=0.5, H=2);   t_e = time.perf_counter() - t0
print(f"stepped {t_s:.3f}s   exact {t_e:.3f}s   ({t_s / t_e:.1f}x)")
```

```python
# thread scaling
import time
import numpy as np
import tfce
from scipy.ndimage import gaussian_filter

rng = np.random.default_rng(0)
B = 16
vols = np.stack(
    [gaussian_filter(rng.standard_normal((91, 109, 91)), 2.0) for _ in range(B)],
    axis=-1,
)

for nt in (1, 2, 4, 6, 8, 12, 16):
    t0 = time.perf_counter()
    tfce.tfce(vols, n_jobs=nt)
    dt = time.perf_counter() - t0
    print(f"{nt:2d} threads: {dt / B * 1000:6.1f} ms/perm")
```

---

## A note on what to time

Timing the transform alone flatters it. In a real permutation test the GLM matters too — and the
permuted GLM here never forms the permuted data, which removes what is otherwise the dominant term
(see [theory](theory.md#the-permuted-glm)). In the MATLAB toolbox that change alone took the GLM from
55% of the loop to 23%.

The thing actually worth optimising is **how many permutations you run**, and that is what the
[tail approximations](theory.md#fewer-permutations) address — roughly a 5× cut for the same accuracy.
A 72× faster transform on 10,000 permutations is slower than a 14× faster transform on 1,000.
