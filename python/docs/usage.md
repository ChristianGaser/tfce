# Usage

Everything is **arrays in, arrays out**. No image objects, no file I/O, no design parsing - those
belong in a layer above, so the core can be dropped into anything.

```python
import tfce
```

---

## Volumes

A single map, shape `(nx, ny, nz)`:

```python
t = tfce.tfce(stat_map, E=0.5, H=2.0)          # -> (nx, ny, nz), float32
```

A block of maps - permutations, say - shape `(nx, ny, nz, B)`:

```python
t = tfce.tfce(perms, E=0.5, H=2.0, n_jobs=-1)  # -> (nx, ny, nz, B)
```

The output has the shape the input had. Values come back as `float32`: the map is a permuted
statistic, and single precision is all it is ever worth.

### Connectivity

```python
tfce.tfce(x, connectivity=26)   # default: shares a face, an edge OR a corner
tfce.tfce(x, connectivity=18)   # face or edge
tfce.tfce(x, connectivity=6)    # face only
```

| | who uses it |
| --- | --- |
| **26** | the MATLAB TFCE toolbox, `fslmaths` - and the default here |
| **18** | - |
| **6** | nilearn (`generate_binary_structure(3, 1)`) |

This matters more than it looks. A looser neighbourhood merges clusters, so extents grow and TFCE
values grow with them. **Two analyses with different connectivity are not comparable.** If you are
reproducing someone's result, match their neighbourhood.

---

## Surfaces

Surface data is a 1-D array of vertices plus the mesh, given as an adjacency:

```python
adj = tfce.adjacency_from_faces(faces, n_vertices)
t   = tfce.tfce(surf_map, adjacency=adj, E=1.0, H=2.0)
```

Faces are **1-based**, which is how GIFTI stores them and therefore what `nibabel` and
`SPM.xVol.G.faces` hand over. An index that names no vertex is refused, rather than read past the end
of the arrays:

```python
tfce.adjacency_from_faces(bad_faces, n_vertices)
# TfceError: faces name a vertex that is not in the data
```

If your faces are 0-based (some tools are), add one.

With `nibabel`:

```python
import nibabel as nib

g      = nib.load("lh.pial.surf.gii")
faces  = g.agg_data("triangle") + 1          # nibabel gives 0-based
n_vert = g.agg_data("pointset").shape[0]

adj = tfce.adjacency_from_faces(faces, n_vert)
```

`adjacency_from_faces` returns `(indptr, indices)` - CSR, 0-based - which you can hand straight back
to `tfce()`, or wrap in a scipy sparse matrix. A `scipy.sparse` matrix is accepted too:

```python
from scipy.sparse import csr_matrix
A = csr_matrix((np.ones(indices.size), indices, indptr), shape=(n_vert, n_vert))
t = tfce.tfce(surf_map, adjacency=A, E=1.0)
```

Batches work the same way, shape `(n_vertices, B)`.

Note the **connectivity argument is ignored** for surfaces: the mesh already says who the neighbours
are.

---

## `E` and `H`

`E` weights **extent**, `H` weights **height**. They set how focal and broad effects trade off.

| data | `E` | `H` |
| --- | --- | --- |
| 3D volume | **0.5** | 2 |
| surface | **1.0** | 2 |
| TBSS (effectively 2D) | **1.0** | 2 |

`E=0.5, H=2` are the values Smith and Nichols established empirically for volumes. The surface and
TBSS values follow the same reasoning for data of lower effective dimension.

**Do not tune these.** They are exposed because different data types need different values, not
because they are knobs. Changing `E` changes what "significant" means, in a way that is not
comparable across studies. (The MATLAB toolbox does not expose them at all, for exactly this reason.)

---

## Two-sided

```python
tfce.tfce(x, two_sided=True)    # default
```

Enhances the negative part of the map as well, giving it a **negative** TFCE value. A map with no
negative values - an F-statistic, say - is completely unaffected by this, so leaving it on costs a
little time and changes nothing.

```python
tfce.tfce(x, two_sided=False)   # positive clusters only
```

---

## Threads

```python
tfce.tfce(perms, n_jobs=-1)     # one thread per map
tfce.tfce(perms, n_jobs=8)      # eight
tfce.tfce(perms, n_jobs=1)      # serial (default)
```

Only the batched forms use this. Permutations are independent of one another, so that is where the
parallelism belongs - every thread stays on its own map and never locks.

Two things to know:

- **The GIL is released** for the whole call, so this is real parallelism and composes with a thread
  pool above it.
- **Do not expect it to scale with your core count.** The max-tree is bound by memory rather than by
  arithmetic. On an 8-core machine it plateaus around **2×** from six threads onwards, and gets
  *worse* beyond eight. `n_jobs=-1` sits on that plateau.

---

## Common patterns

**Only the maximum, for an FWE null:**

```python
null_max = np.abs(tfce.tfce(perms, n_jobs=-1)).reshape(-1, n_perm).max(axis=0)
```

**Masked data.** TFCE needs the spatial structure, so pass the full volume and mask afterwards.
Elements at zero are not part of any cluster:

```python
stat_map[~mask] = 0
t = tfce.tfce(stat_map)
```

**Memory.** A batch costs `n_elements × B` floats twice over (in and out). For a 91×109×91 volume
that is about 3.6 MB per map, so a block of 16 is ~115 MB. If you are permuting thousands of times,
do it in blocks rather than all at once.
