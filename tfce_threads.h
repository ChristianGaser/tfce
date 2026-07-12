/*
 * Minimal portable threading for the TFCE mex-files.
 *
 * POSIX uses pthreads. Windows has no pthread.h when MATLAB compiles with MSVC,
 * so the same handful of calls is mapped onto Win32 threads and critical
 * sections. Only what the mex-files actually need is wrapped: create a thread
 * that runs a void *(*)(void *), join it, and a mutex.
 *
 * Christian Gaser
 */

#ifndef TFCE_THREADS_H
#define TFCE_THREADS_H

#include <stdlib.h>

#ifdef _WIN32

#include <windows.h>

typedef HANDLE             tfce_thread_t;
typedef CRITICAL_SECTION   tfce_mutex_t;

/* Win32 thread functions have a different signature, so the worker and its
   argument are passed through a small trampoline that is freed on entry. */
typedef struct {
  void *(*fn)(void *);
  void *arg;
} tfce_trampoline;

static DWORD WINAPI tfce_thread_entry(LPVOID p)
{
  tfce_trampoline tr = *(tfce_trampoline *)p;
  free(p);
  tr.fn(tr.arg);
  return 0;
}

static int tfce_thread_create(tfce_thread_t *t, void *(*fn)(void *), void *arg)
{
  tfce_trampoline *tr = (tfce_trampoline *) malloc(sizeof(tfce_trampoline));
  if (!tr) return -1;
  tr->fn  = fn;
  tr->arg = arg;
  *t = CreateThread(NULL, 0, tfce_thread_entry, tr, 0, NULL);
  if (*t == NULL) { free(tr); return -1; }
  return 0;
}

static void tfce_thread_join(tfce_thread_t t)
{
  WaitForSingleObject(t, INFINITE);
  CloseHandle(t);
}

static int  tfce_mutex_init(tfce_mutex_t *m)    { InitializeCriticalSection(m); return 0; }
static void tfce_mutex_lock(tfce_mutex_t *m)    { EnterCriticalSection(m); }
static void tfce_mutex_unlock(tfce_mutex_t *m)  { LeaveCriticalSection(m); }
static void tfce_mutex_destroy(tfce_mutex_t *m) { DeleteCriticalSection(m); }

#else

#include <pthread.h>

typedef pthread_t        tfce_thread_t;
typedef pthread_mutex_t  tfce_mutex_t;

static int tfce_thread_create(tfce_thread_t *t, void *(*fn)(void *), void *arg)
{
  return pthread_create(t, NULL, fn, arg);
}

static void tfce_thread_join(tfce_thread_t t)   { pthread_join(t, NULL); }

static int  tfce_mutex_init(tfce_mutex_t *m)    { return pthread_mutex_init(m, NULL); }
static void tfce_mutex_lock(tfce_mutex_t *m)    { pthread_mutex_lock(m); }
static void tfce_mutex_unlock(tfce_mutex_t *m)  { pthread_mutex_unlock(m); }
static void tfce_mutex_destroy(tfce_mutex_t *m) { pthread_mutex_destroy(m); }

#endif

#endif /* TFCE_THREADS_H */
