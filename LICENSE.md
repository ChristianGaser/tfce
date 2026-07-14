# Licensing

This repository holds two things with two different licences. The split is not
arbitrary: it follows exactly what each part is derived from.

| | Licence | Why |
| --- | --- | --- |
| [`c/`](c/) — the TFCE core | **BSD-3-Clause** ([LICENSE](c/LICENSE)) | Original work. Contains no SPM, SnPM or PALM code, and includes no header from any of them — only libc and pthreads. Nothing obliges it to be GPL. |
| [`python/`](python/) — the `tfce` package | **BSD-3-Clause** ([LICENSE](python/LICENSE)) | Original work on top of numpy and scipy. It compiles in `c/`, and nothing else. |
| [`matlab/`](matlab/) — the SPM toolbox | **GPL-2.0-or-later** ([COPYING](COPYING)) | Genuinely derived from GPL code: `snpm_P_FDR.m` is from SnPM, `cat_spm_results_ui.m` from SPM's `spm_results_ui`, the `palm_*` subfunctions from PALM, and every `.m` file calls SPM. |

## What this means in practice

**If you want to use the TFCE core, or the Python package: it is BSD-3-Clause.**
Depend on it, vendor it, ship it in a closed product — the three clauses are the
whole of it. This is deliberate. `nilearn` is BSD-3-Clause and `nipreps` is
Apache-2.0, and permissive projects rightly refuse to take on a GPL dependency; a
GPL licence here would have shut both of them out of a piece of code that owes
them nothing.

**If you use the MATLAB/SPM toolbox: it is GPL-2.0-or-later**, as it always has
been. Nothing about it changes. BSD-3-Clause is GPLv2-compatible, so the same C
core can be, and is, compiled into the mex-files that the GPL toolbox ships.

## Why the C is not GPL, even though the toolbox is

The MATLAB toolbox is GPL because it is *built on* SPM — it calls `spm_*`
throughout and carries code taken from SnPM and PALM. None of that is true of the
C core. It is a max-tree and a union-find, written from scratch, and it depends on
nothing but the C standard library. It is compiled into the MATLAB mex-files and
into the Python extension alike, and being permissive lets it serve both without
forcing the GPL onto a Python ecosystem that neither wants it nor is derived from
anything that requires it.

Copyright (c) 2020-2026, Christian Gaser
Structural Brain Mapping Group, Departments of Psychiatry and Neurology,
Jena University Hospital.
