/* SPDX-License-Identifier: BSD-3-Clause */
/* Copyright (c) 2020-2026, Christian Gaser. See LICENSE. */
/*
 * Keep a Linux binary loadable on an older glibc than the one it was built on.
 *
 * An ELF shared object records, for every symbol it imports, the glibc symbol
 * version it was linked against, and the dynamic loader treats that as a hard
 * minimum. A mex-file built on Ubuntu 22.04 imports pthread_create@GLIBC_2.34,
 * so MATLAB on any older system refuses it with
 *
 *   Invalid MEX-file '.../tfceMex_maxtree.mexa64':
 *   /lib64/libc.so.6: version `GLIBC_2.34' not found
 *
 * although nothing in the toolbox needs anything new. Three symbols carry the
 * whole incompatibility, and neither of them by our doing:
 *
 *   pthread_create, pthread_join   glibc 2.34 merged libpthread into libc and
 *                                  gave every pthread symbol a fresh 2.34
 *                                  version; the functions did not change
 *   pow                            glibc 2.29 added a more accurate
 *                                  implementation under a new version
 *
 * glibc keeps the previous implementations as compat symbols indefinitely, so
 * requesting them by name costs nothing: the 2.2.5 pow is the one every binary
 * built before 2019 used, and the 2.2.5 pthread_create is the same function
 * that 2.34 merely renamed. .symver rewrites the *reference*, hence this must
 * be included by every translation unit that goes into a mex-file.
 *
 * The floor drops from glibc 2.34 (Ubuntu 22.04, 2022) to glibc 2.14 (2011) --
 * what remains is memcpy@GLIBC_2.14, well below the glibc any MATLAB release
 * requires of its own accord. compile_mex.yml enforces that floor on every
 * build, because nothing on the build machine notices when it slips again.
 *
 * The failure mode is benign: if a version named here did not exist, the link
 * would fail loudly on the build machine. It cannot quietly emit a binary that
 * is broken somewhere else.
 *
 * Not needed by the Python binding, whose wheels are built in a manylinux_2_28
 * container and therefore link against glibc 2.28 to begin with.
 *
 * Define TFCE_NO_GLIBC_COMPAT to switch this off.
 *
 * Christian Gaser
 */

#ifndef TFCE_GLIBC_COMPAT_H
#define TFCE_GLIBC_COMPAT_H

#if defined(__linux__) && !defined(TFCE_NO_GLIBC_COMPAT)

/* __GLIBC__ is only defined once a libc header has been seen. */
#include <limits.h>

#if defined(__GLIBC__) && defined(__ELF__)

/* The baseline version node is the one the ABI of that architecture started
   at, so it exists in every glibc that can run the binary at all. It is not
   the same everywhere: x86-64 entered at 2.2.5, aarch64 only at 2.17. On an
   architecture that is not listed the header does nothing, which is the safe
   direction -- an unknown baseline would be a link error, not a fallback. */
#if defined(__x86_64__) && !defined(__ILP32__)
#define TFCE_GLIBC_BASE "GLIBC_2.2.5"
#elif defined(__aarch64__)
#define TFCE_GLIBC_BASE "GLIBC_2.17"
#endif

#ifdef TFCE_GLIBC_BASE
__asm__(".symver pthread_create,pthread_create@" TFCE_GLIBC_BASE);
__asm__(".symver pthread_join,pthread_join@" TFCE_GLIBC_BASE);
__asm__(".symver pow,pow@" TFCE_GLIBC_BASE);
#endif

#endif /* __GLIBC__ && __ELF__ */

#endif /* __linux__ && !TFCE_NO_GLIBC_COMPAT */

#endif /* TFCE_GLIBC_COMPAT_H */
