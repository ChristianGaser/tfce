"""
Build the exact TFCE max-tree as a Python extension.

The C core is shared with the MATLAB toolbox and lives at the repository root.
There is exactly one implementation of TFCE in this repository and both bindings
sit on it -- but a source distribution cannot reach outside its own directory, so
the core is *vendored* into src/tfce/_c/ at build time.

That gives two cases, and both work:

  * building from the git checkout -- the core is copied down from the repository
    root, so an edit to the C is picked up by the next build with nothing to
    remember;
  * building from an sdist -- there is no repository root, and the copy already
    made when the sdist was created is used as it stands.

src/tfce/_c/ is therefore generated, and is git-ignored.
"""

import os
import shutil

import numpy as np
from setuptools import Extension, setup
from Cython.Build import cythonize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
VENDOR = os.path.join("src", "tfce", "_c")

CORE_FILES = [
    "tfce_capi.c",
    "tfce_capi.h",
    "tfce_maxtree.h",
    "tfce_batch.h",
    "tfce_threads.h",
]


def vendor_core():
    """Copy the C core into the package, so the sdist is self-contained."""
    dest = os.path.join(HERE, VENDOR)
    os.makedirs(dest, exist_ok=True)

    for name in CORE_FILES:
        src = os.path.join(ROOT, name)
        dst = os.path.join(dest, name)

        if os.path.exists(src):
            # a git checkout: the root is the source of truth, always
            shutil.copyfile(src, dst)
        elif not os.path.exists(dst):
            raise RuntimeError(
                f"{name} is neither at the repository root nor vendored in "
                f"{VENDOR}. A source distribution should carry it; a checkout "
                f"should have it one directory up."
            )


vendor_core()

extra_compile_args = ["-O3"]
extra_link_args = []
if os.name == "nt":
    extra_compile_args = ["/O2"]
else:
    extra_compile_args += ["-std=c99", "-pthread"]
    extra_link_args += ["-pthread"]

ext = Extension(
    "tfce._maxtree",
    sources=[
        os.path.join("src", "tfce", "_maxtree.pyx"),
        os.path.join(VENDOR, "tfce_capi.c"),
    ],
    include_dirs=[np.get_include(), VENDOR],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
)

setup(ext_modules=cythonize([ext], language_level=3))
