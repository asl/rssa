/*
 *   R package for Singular Spectrum Analysis
 *   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
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

#include <complex.h>
#include <fftw3.h>
#include "extmat.h"
#include "hankel.h"

typedef struct {
  fftw_complex * circ_freq;
  fftw_plan r2c_plan;
  fftw_plan c2r_plan;
  R_len_t window;
  R_len_t length;
} hankel_matrix;

unsigned _hankel_nrow(const void *matrix) {
  const hankel_matrix *h = matrix;
  return h->window;
}

unsigned _hankel_ncol(const void *matrix) {
  const hankel_matrix *h = matrix;
  return h->length - h->window + 1;
}

void _free_circulant(hankel_matrix *h) {
  fftw_free(h->circ_freq);
  fftw_destroy_plan(h->r2c_plan);
  fftw_destroy_plan(h->c2r_plan);
}

void _initialize_circulant(hankel_matrix *h,
                           const double *F, R_len_t N, R_len_t L) {
  R_len_t K = N - L + 1, i;
  fftw_complex *ocirc;
  fftw_plan p1, p2;
  double *circ;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Estimate the best plans for given input length */
  p1 = fftw_plan_dft_r2c_1d(N, circ, ocirc, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_c2r_1d(N, ocirc, circ, FFTW_ESTIMATE);

  /* Fill input buffer */
  for (i = K-1; i < N; ++i)
    circ[i - K + 1] = F[i];

  for (i = 0; i < K-1; ++i) {
    circ[L + i] = F[i];
  }

  /* Run the plan on input data */
  fftw_execute(p1);

  /* Cleanup and return */
  fftw_free(circ);

  h->circ_freq = ocirc;
  h->r2c_plan = p1;
  h->c2r_plan = p2;
  h->window = L;
  h->length = N;
}

void _hankel_matmul(double* out,
                    const double* v,
                    const void* matrix) {
  const hankel_matrix *h = matrix;
  R_len_t N = h->length, L = h->window;
  R_len_t K = N - L + 1, i;
  double *circ;
  fftw_complex *ocirc;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Fill the arrays */
  for (i = 0; i < K; ++i)
    circ[i] = v[K - i - 1];
  memset(circ + K, 0, (L - 1)*sizeof(double));

  /* Compute the FFT of the reversed vector v */
  fftw_execute_dft_r2c(h->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < (N/2 + 1); ++i)
    ocirc[i] = ocirc[i] * h->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(h->c2r_plan, ocirc, circ);

  /* Cleanup and return */
  for (i = 0; i < L; ++i)
    out[i] = circ[i] / N;

  fftw_free(circ);
  fftw_free(ocirc);
}

void _hankel_tmatmul(double* out,
                     const double* v,
                     const void* matrix) {
  const hankel_matrix *h = matrix;
  R_len_t N = h->length, L = h->window;
  R_len_t K = N - L + 1, i;
  double *circ;
  fftw_complex *ocirc;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Fill the arrays */
  memset(circ, 0, (K - 1)*sizeof(double));
  for (i = 0; i < L; ++i)
    circ[i + K - 1] = v[L - i - 1];

  /* Compute the FFT of the reversed vector v */
  fftw_execute_dft_r2c(h->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < (N/2 + 1); ++i)
    ocirc[i] = ocirc[i] * h->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(h->c2r_plan, ocirc, circ);

  /* Cleanup and return */
  for (i = 0; i < K; ++i)
    out[i] = circ[i + L - 1] / N;

  fftw_free(circ);
  fftw_free(ocirc);
}

static void hmatFinalizer(SEXP ptr) {
  ext_matrix *e;
  hankel_matrix *h;

  if (TYPEOF(ptr) != EXTPTRSXP)
    return;

  e = R_ExternalPtrAddr(ptr);
  if (!e)
    return;

  if (!strcmp(e->type, "hankel matrix"))
    return;

  h = e->matrix;

  _free_circulant(h);
  Free(h);

  Free(e);
  R_ClearExternalPtr(ptr);
}

SEXP initialize_hmat(SEXP F, SEXP window) {
  R_len_t N, L;
  hankel_matrix *h;
  ext_matrix *e;
  SEXP hmat;

  N = length(F);
  L = INTEGER(window)[0];

  /* Allocate memory */
  e = Calloc(1, ext_matrix);
  e->type = "hankel matrix";
  e->mulfn = _hankel_matmul;
  e->tmulfn = _hankel_tmatmul;
  e->ncol = _hankel_ncol;
  e->nrow = _hankel_nrow;

  /* Build toeplitz circulants for hankel matrix */
  h = Calloc(1, hankel_matrix);
  _initialize_circulant(h, REAL(F), N, L);
  e->matrix = h;

  /* Make an external pointer envolope */
  hmat = R_MakeExternalPtr(e, install("external matrix"), R_NilValue);
  R_RegisterCFinalizer(hmat, hmatFinalizer);

  return hmat;
}

SEXP is_hmat(SEXP ptr) {
  SEXP ans;
  ext_matrix *e = NULL;

  PROTECT(ans = allocVector(LGLSXP, 1));
  LOGICAL(ans)[0] = 1;

  /* hmat object is an external pointer */
  if (TYPEOF(ptr) != EXTPTRSXP)
    LOGICAL(ans)[0] = 0;

  /* tag should be 'external matrix' */
  if (LOGICAL(ans)[0] &&
      R_ExternalPtrTag(ptr) != install("external matrix"))
    LOGICAL(ans)[0] = 0;

  /* pointer itself should not be null */
  if (LOGICAL(ans)[0]) {
    e = R_ExternalPtrAddr(ptr);
    if (!e)
      LOGICAL(ans)[0] = 0;
  }

  /* finally, type should be `hankel matrix' */
  if (LOGICAL(ans)[0] && e &&
      strcmp(e->type, "hankel matrix") != 0)
    LOGICAL(ans)[0] = 0;

  UNPROTECT(1);

  return ans;
}

SEXP hankel_rows(SEXP ptr) {
  SEXP tchk;
  SEXP ans = NILSXP;

  /* Perform a type checking */
  PROTECT(tchk = is_hmat(ptr));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e = R_ExternalPtrAddr(ptr);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = _hankel_nrow(e->matrix);
    UNPROTECT(1);
  } else
    error("pointer provided is not a hankel matrix");

  UNPROTECT(1);

  return ans;
}

SEXP hankel_cols(SEXP ptr) {
  SEXP tchk;
  SEXP ans = NILSXP;

  /* Perform a type checking */
  PROTECT(tchk = is_hmat(ptr));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e = R_ExternalPtrAddr(ptr);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = _hankel_ncol(e->matrix);
    UNPROTECT(1);
  } else
    error("pointer provided is not a hankel matrix");

  UNPROTECT(1);

  return ans;
}

SEXP hmatmul(SEXP hmat, SEXP v, SEXP transposed) {
  SEXP Y = NILSXP, tchk;

  /* Perform a type checking */
  PROTECT(tchk = is_hmat(hmat));

  if (LOGICAL(tchk)[0]) {
    R_len_t K, L;
    ext_matrix *e;
    hankel_matrix *h;

    /* Grab needed data */
    e = R_ExternalPtrAddr(hmat);
    h = e->matrix;

    L = (LOGICAL(transposed)[0] ? h->length - h->window + 1 : h->window);

    /* Check agains absurd values of inputs */
    K = length(v);
    if (K + L - 1 != h->length)
      error("invalid length of input vector 'v'");

    /* Allocate output buffer */
    PROTECT(Y = allocVector(REALSXP, L));

    /* Calculate the product */
    if (LOGICAL(transposed)[0])
      _hankel_tmatmul(REAL(Y), REAL(v), h);
    else
      _hankel_matmul(REAL(Y), REAL(v), h);

    UNPROTECT(1);
  } else
    error("pointer provided is not a hankel matrix");

  UNPROTECT(1);

  return Y;
}
