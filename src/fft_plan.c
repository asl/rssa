/*
 *   R package for Singular Spectrum Analysis
 *   Copyright (c) 2013 Anton Korobeynikov <asl@math.spbu.ru>
 *   Copyright (c) 2013 Alexander Shlemov <shlemovalex@gmail.com>
 *
 *   This program is free software; you can redistribute it
 *   and/or modify it under the terms of the GNU General Public
 *   License as published by the Free Software Foundation;
 *   either version 2 of the License, or (at your option)
 *   any later version.
 *
 *   This program is distributed in the hope that it will be
 *   useful, but WITHOUT ANY WARRANTY; without even the implied
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *   PURPOSE.  See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public
 *   License along with this program; if not, write to the
 *   Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
 *   MA 02139, USA.
 */

#include <R.h>
#include <Rinternals.h>

#include "fft_plan.h"

unsigned valid_plan(const fft_plan *f, R_len_t d, const R_len_t *dim) {
  if (d != f->d)
    return 0;

  int i;
  for (i = 0; i < d; ++i) {
    if (dim[i] != f->dim[i])
      return 0;
  }

  return 1;
}

#if HAVE_FFTW3_H
static void free_plan(fft_plan *f) {
  fftw_destroy_plan(f->r2c_plan);
  fftw_destroy_plan(f->c2r_plan);

  Free(f->dim);
}

static void initialize_plan(fft_plan *f, R_len_t N) {
  fftw_complex *ocirc;
  double *circ;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Estimate the best plans for given input length */
  f->r2c_plan = fftw_plan_dft_r2c_1d(N, circ, ocirc, FFTW_ESTIMATE);
  f->c2r_plan = fftw_plan_dft_c2r_1d(N, ocirc, circ, FFTW_ESTIMATE);

  f->d = 1;
  f->dim = Calloc(1, R_len_t);
  f->dim[0] = N;

  fftw_free(circ);
  fftw_free(ocirc);
}
#else
static void initialize_plan(fft_plan *f, R_len_t N) {
  f->d = 1;
  f->dim = Calloc(1, R_len_t);
  f->dim[0] = N;
}

static void free_plan(fft_plan *f) {
  Free(f->dim);
}
#endif

static void fft_plan_finalizer(SEXP ptr) {
  fft_plan *f;

  if (TYPEOF(ptr) != EXTPTRSXP)
    return;

  f = R_ExternalPtrAddr(ptr);
  if (!f)
    return;

  free_plan(f);
  Free(f);

  R_ClearExternalPtr(ptr);
}

SEXP initialize_fft_plan(SEXP rN) {
  R_len_t N;
  fft_plan *f;
  SEXP res;

  N = INTEGER(rN)[0];

  /* Allocate memory */
  f = Calloc(1, fft_plan);

  /* Do actual plan initialization */
  initialize_plan(f, N);

  /* Make an external pointer envelope */
  PROTECT(res = R_MakeExternalPtr(f, install("fft plan"), R_NilValue));
  R_RegisterCFinalizer(res, fft_plan_finalizer);
  UNPROTECT(1);

  return res;
}

SEXP is_fft_plan(SEXP ptr) {
  SEXP ans;
  fft_plan *f;

  PROTECT(ans = allocVector(LGLSXP, 1));
  LOGICAL(ans)[0] = 1;

  /* object is an external pointer */
  if (TYPEOF(ptr) != EXTPTRSXP)
    LOGICAL(ans)[0] = 0;

  /* tag should be 'fft plan' */
  if (LOGICAL(ans)[0] &&
      R_ExternalPtrTag(ptr) != install("fft plan"))
    LOGICAL(ans)[0] = 0;

  /* pointer itself should not be null */
  if (LOGICAL(ans)[0]) {
    f = R_ExternalPtrAddr(ptr);
    if (!f)
      LOGICAL(ans)[0] = 0;
  }

  UNPROTECT(1);

  return ans;
}


SEXP dim_fft_plan(SEXP fftplan) {
  /* Perform a type checking */
  if (!LOGICAL(is_fft_plan(fftplan))[0]) {
    error("pointer provided is not a fft plan");
    return NILSXP;
  }

  SEXP ans;
  fft_plan *f = R_ExternalPtrAddr(fftplan);

  PROTECT(ans = allocVector(INTSXP, f->d));
  memcpy(INTEGER(ans), f->dim, f->d * sizeof(R_len_t));

  UNPROTECT(1);
  return ans;
}
