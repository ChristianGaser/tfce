/* SPDX-License-Identifier: BSD-3-Clause */
/* Copyright (c) 2020-2026, Christian Gaser. See LICENSE. */
/*
 * Threaded batch driver for the TFCE max-tree: many maps at once, one map per
 * thread.
 *
 * This used to live inside tfceMex_maxtree_batch.c, which meant it could only be
 * reached from MATLAB. It is plain C and has nothing to do with MATLAB, so it
 * lives here instead and is shared by the mex-file and by the C API that the
 * language bindings (Python/Cython, and anything else) sit on.
 *
 * Parallelising across permutations rather than within a single map keeps every
 * thread on its own data, so there is no locking on the hot path -- the mutex is
 * touched once per map, only to hand out the next index.
 *
 * Christian Gaser
 */

#ifndef TFCE_BATCH_H
#define TFCE_BATCH_H

#include "tfce_threads.h"
#include "tfce_maxtree.h"

#define TFCE_MAX_THREADS 256

typedef struct {
  const tfce_val *T;         /* N x B, one map per column   */
  tfce_val *out;             /* N x B, written by this call */
  int N, B;
  const Neigh *nb;
  double E, H;
  int calc_neg;
  int *next;                 /* next map to claim, shared   */
  tfce_mutex_t *mtx;
  int ok;
} TfceBatchArgs;

static void *tfce_batch_thread(void *pa)
{
  TfceBatchArgs *A = (TfceBatchArgs *)pa;
  Workspace ws;

  /* one workspace per thread, reused across every map that thread claims */
  if (!ws_alloc(&ws, A->N)) { A->ok = 0; return NULL; }

  for (;;) {
    int b;

    tfce_mutex_lock(A->mtx);
    b = (*A->next)++;
    tfce_mutex_unlock(A->mtx);

    if (b >= A->B) break;

    maxtree_map(&ws, A->T + (size_t)b * A->N, A->out + (size_t)b * A->N,
                A->nb, A->E, A->H, A->calc_neg);
  }

  ws_free(&ws);
  return NULL;
}

/*
 * TFCE of B maps. Returns 1 on success, 0 if a worker could not allocate.
 *
 * n_threads <= 0 means one per map. The count is capped at B, because a thread
 * with no map to work on is only overhead.
 */
static int tfce_batch_run(const tfce_val *T, tfce_val *out, int N, int B,
                          const Neigh *nb, double E, double H, int calc_neg,
                          int n_threads)
{
  TfceBatchArgs args[TFCE_MAX_THREADS];
  tfce_thread_t th[TFCE_MAX_THREADS];
  tfce_mutex_t mtx;
  int i, next = 0, started = 0, ok = 1;

  if (n_threads <= 0) n_threads = B;
  if (n_threads > B)  n_threads = B;
  if (n_threads > TFCE_MAX_THREADS) n_threads = TFCE_MAX_THREADS;
  if (n_threads < 1)  n_threads = 1;

  if (tfce_mutex_init(&mtx) != 0) return 0;

  for (i = 0; i < n_threads; i++) {
    args[i].T = T; args[i].out = out;
    args[i].N = N; args[i].B = B;
    args[i].nb = nb;
    args[i].E = E; args[i].H = H;
    args[i].calc_neg = calc_neg;
    args[i].next = &next;
    args[i].mtx  = &mtx;
    args[i].ok   = 1;
  }

  for (i = 0; i < n_threads; i++)
    if (tfce_thread_create(&th[i], tfce_batch_thread, &args[i]) == 0) started++;
    else break;

  /* if a thread could not be created, this one picks up what is left */
  if (started < n_threads) tfce_batch_thread(&args[0]);

  for (i = 0; i < started; i++) tfce_thread_join(th[i]);

  tfce_mutex_destroy(&mtx);

  for (i = 0; i < n_threads; i++) if (!args[i].ok) ok = 0;
  return ok;
}

#endif /* TFCE_BATCH_H */
