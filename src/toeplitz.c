/*
 *   R package for Singular Spectrum Analysis
 *   Copyright (c) 2009-2010 Anton Korobeynikov <asl@math.spbu.ru>
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

#include "extmat.h"
#include "config.h"
#if HAVE_FFTW3_H
#include <fftw3.h>
#else
#include <R_ext/Applic.h>
#endif

typedef struct {
#if HAVE_FFTW3_H
  fftw_complex * circ_freq;
  fftw_plan r2c_plan;
  fftw_plan c2r_plan;
#else
  complex double * circ_freq;
#endif
  R_len_t window;
  R_len_t length;
} toeplitz_matrix;

static unsigned toeplitz_nrow(const void *matrix) {
  const toeplitz_matrix *t = matrix;
  return t->window;
}

static unsigned toeplitz_ncol(const void *matrix) {
  const toeplitz_matrix *t = matrix;
  return t->length - t->window + 1;
}

#if HAVE_FFTW3_H
static void free_circulant(toeplitz_matrix *t) {
  fftw_free(t->circ_freq);
  fftw_destroy_plan(t->r2c_plan);
  fftw_destroy_plan(t->c2r_plan);
}

static void initialize_circulant(toeplitz_matrix *t,
                                 const double *R, R_len_t L) {
  R_len_t N = 2*L - 1, i;
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
  for (i = 0; i < L; ++i)
    circ[i] = R[i];

  for (i = 0; i < L-1; ++i)
    circ[L + i] = R[L-i-1];

  /* Run the plan on input data */
  fftw_execute(p1);

  /* Cleanup and return */
  fftw_free(circ);

  t->circ_freq = ocirc;
  t->r2c_plan = p1;
  t->c2r_plan = p2;
  t->window = L;
  t->length = N;
}

static void toeplitz_matmul(double* out,
                            const double* v,
                            const void* matrix) {
  const toeplitz_matrix *t = matrix;
  R_len_t N = t->length, L = t->window;
  R_len_t K = N - L + 1, i;
  double *circ;
  fftw_complex *ocirc;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Fill the arrays */
  for (i = 0; i < K; ++i)
    circ[i] = v[i];
  memset(circ + K, 0, (L - 1)*sizeof(double));

  /* Compute the FFT of the vector v */
  fftw_execute_dft_r2c(t->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < (N/2 + 1); ++i)
    ocirc[i] = ocirc[i] * t->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(t->c2r_plan, ocirc, circ);

  /* Cleanup and return */
  for (i = 0; i < L; ++i)
    out[i] = circ[i] / N;

  fftw_free(circ);
  fftw_free(ocirc);
}

static void toeplitz_tmatmul(double* out,
                             const double* v,
                             const void* matrix) {
  const toeplitz_matrix *t = matrix;
  R_len_t N = t->length, L = t->window;
  R_len_t K = N - L + 1, i;
  double *circ;
  fftw_complex *ocirc;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Fill the arrays */
  memset(circ, 0, (K - 1)*sizeof(double));
  for (i = 0; i < L; ++i)
    circ[i + K - 1] = v[i];

  /* Compute the FFT of the reversed vector v */
  fftw_execute_dft_r2c(t->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < (N/2 + 1); ++i)
    ocirc[i] = ocirc[i] * t->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(t->c2r_plan, ocirc, circ);

  /* Cleanup and return */
  for (i = 0; i < K; ++i)
    out[i] = circ[i + L - 1] / N;

  fftw_free(circ);
  fftw_free(ocirc);
}

static void calc_Lcor(double *R, const double *F, R_len_t N, R_len_t L) {
  R_len_t No, i;
  double *circ;
  fftw_complex *ocirc;
  fftw_plan p1, p2;

  /* Allocate needed memory */
  No = N + L - 1;
  circ = (double*) fftw_malloc(No * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((No/2 + 1) * sizeof(fftw_complex));

  /* Estimate the best plans for given input length */
  p1 = fftw_plan_dft_r2c_1d(No, circ, ocirc, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_c2r_1d(No, ocirc, circ, FFTW_ESTIMATE);

  memcpy(circ, F, N*sizeof(double));
  memset(circ + N, 0, (L - 1)*sizeof(double));

  /* Run the plan on input data */
  fftw_execute(p1);

  /* Auto dot-product */
  for (i = 0; i < No/2 + 1; ++i)
    ocirc[i] = ocirc[i] * conj(ocirc[i]);

  fftw_execute(p2);

  /* Return */
  for (i = 0; i < L; ++i)
    R[i] = circ[i] / (N - i) / No;

  /* Cleanup */
  fftw_free(circ);
  fftw_free(ocirc);
  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
}

#else

static void free_circulant(toeplitz_matrix *t) {
  Free(t->circ_freq);
}

static void initialize_circulant(toeplitz_matrix *t,
                                 const double *R, R_len_t L) {
  R_len_t N = 2*L - 1, i;
  int *iwork, maxf, maxp;
  double *work;
  complex double *circ;

  /* Allocate needed memory */
  circ = Calloc(N, complex double);

  /* Estimate the best plans for given input length */
  fft_factor(N, &maxf, &maxp);
  if (maxf == 0)
    error("fft factorization error");

  work = Calloc(4 * maxf, double);
  iwork = Calloc(maxp, int);

  /* Fill input buffer */
  for (i = 0; i < L; ++i)
    circ[i] = R[i];

  for (i = 0; i < L-1; ++i)
    circ[L + i] = R[L-i-1];

  /* Run the plan on input data */
  fft_work((double*)circ, ((double*)circ)+1, 1, N, 1, -2, work, iwork);

  /* Cleanup and return */
  Free(work);
  Free(iwork);

  t->circ_freq = circ;
  t->window = L;
  t->length = N;
}

static void toeplitz_matmul(double* out,
                            const double* v,
                            const void* matrix) {
  const toeplitz_matrix *t = matrix;
  R_len_t N = t->length, L = t->window;
  R_len_t K = N - L + 1, i;
  double *work;
  complex double *circ;
  int *iwork, maxf, maxp;

  /* Estimate the best plans for given input length */
  fft_factor(N, &maxf, &maxp);
  if (maxf == 0)
    error("fft factorization error");

  /* Allocate needed memory */
  circ = Calloc(N, complex double);
  work = Calloc(4 * maxf, double);
  iwork = Calloc(maxp, int);

  /* Fill the arrays */
  for (i = 0; i < K; ++i)
    circ[i] = v[i];

  /* Compute the FFT of the vector v */
  fft_work((double*)circ, ((double*)circ)+1, 1, N, 1, -2, work, iwork);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < N; ++i)
    circ[i] = circ[i] * t->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fft_factor(N, &maxf, &maxp);
  fft_work((double*)circ, ((double*)circ)+1, 1, N, 1, +2, work, iwork);

  /* Cleanup and return */
  for (i = 0; i < L; ++i)
    out[i] = creal(circ[i]) / N;

  Free(circ);
  Free(work);
  Free(iwork);
}

static void toeplitz_tmatmul(double* out,
                             const double* v,
                             const void* matrix) {
  const toeplitz_matrix *t = matrix;
  R_len_t N = t->length, L = t->window;
  R_len_t K = N - L + 1, i;
  double *work;
  complex double *circ;
  int *iwork, maxf, maxp;

  /* Estimate the best plans for given input length */
  fft_factor(N, &maxf, &maxp);
  if (maxf == 0)
    error("fft factorization error");

  /* Allocate needed memory */
  circ = Calloc(N, complex double);
  work = Calloc(4 * maxf, double);
  iwork = Calloc(maxp, int);

  /* Fill the arrays */
  for (i = 0; i < L; ++i)
    circ[i + K - 1] = v[i];

  /* Compute the FFT of the reversed vector v */
  fft_work((double*)circ, ((double*)circ)+1, 1, N, 1, -2, work, iwork);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < N; ++i)
    circ[i] = circ[i] * t->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fft_factor(N, &maxf, &maxp);
  fft_work((double*)circ, ((double*)circ)+1, 1, N, 1, +2, work, iwork);

  /* Cleanup and return */
  for (i = 0; i < K; ++i)
    out[i] = creal(circ[i + L - 1]) / N;

  Free(circ);
  Free(work);
  Free(iwork);
}

static void calc_Lcor(double *R, const double *F, R_len_t N, R_len_t L) {
  R_len_t No = N + L - 1, i;
  complex double *circ;
  double *work;
  int *iwork, maxf, maxp;

  /* Estimate the best plans for given input length */
  fft_factor(No, &maxf, &maxp);
  if (maxf == 0)
    error("fft factorization error");

  /* Allocate needed memory */
  circ = Calloc(No, complex double);
  work = Calloc(4 * maxf, double);
  iwork = Calloc(maxp, int);

  for (i = 0; i < N; ++i)
    circ[i] = F[i];

  /* Compute the FFT of the padded input vector */
  fft_work((double*)circ, ((double*)circ)+1, 1, No, 1, -2, work, iwork);

  /* Auto dot-product */
  for (i = 0; i < No; ++i)
    circ[i] = circ[i] * conj(circ[i]);

  /* Compute the reverse transform to obtain result */
  fft_factor(No, &maxf, &maxp);
  fft_work((double*)circ, ((double*)circ)+1, 1, No, 1, +2, work, iwork);

  /* Return */
  for (i = 0; i < L; ++i)
    R[i] = creal(circ[i]) / (N - i) / No;

  /* Cleanup */
  Free(circ);
  Free(work);
  Free(iwork);
}
#endif

static void tmat_finalizer(SEXP ptr) {
  ext_matrix *e;
  toeplitz_matrix *h;

  if (TYPEOF(ptr) != EXTPTRSXP)
    return;

  e = R_ExternalPtrAddr(ptr);
  if (!e)
    return;

  if (strcmp(e->type, "toeplitz matrix"))
    return;

  h = e->matrix;

  free_circulant(h);
  Free(h);

  Free(e);
  R_ClearExternalPtr(ptr);
}

SEXP initialize_tmat(SEXP R) {
  R_len_t L;
  toeplitz_matrix *t;
  ext_matrix *e;
  SEXP tmat;

  L = length(R);

  /* Allocate memory */
  e = Calloc(1, ext_matrix);
  e->type = "toeplitz matrix";
  e->mulfn = toeplitz_matmul;
  e->tmulfn = toeplitz_tmatmul;
  e->ncol = toeplitz_ncol;
  e->nrow = toeplitz_nrow;

  /* Build toeplitz circulants for toeplitz matrix */
  t = Calloc(1, toeplitz_matrix);
  initialize_circulant(t, REAL(R), L);
  e->matrix = t;

  /* Make an external pointer envelope */
  tmat = R_MakeExternalPtr(e, install("external matrix"), R_NilValue);
  R_RegisterCFinalizer(tmat, tmat_finalizer);

  return tmat;
}

SEXP is_tmat(SEXP ptr) {
  SEXP ans = NILSXP, tchk;
  ext_matrix *e = NULL;

  PROTECT(ans = allocVector(LGLSXP, 1));
  LOGICAL(ans)[0] = 1;

  /* Object should be external matrix */
  PROTECT(tchk = is_extmat(ptr));

  /* pointer itself should not be null */
  if (LOGICAL(tchk)[0]) {
    e = R_ExternalPtrAddr(ptr);
    if (!e)
      LOGICAL(ans)[0] = 0;
  } else
    LOGICAL(ans)[0] = 0;

  /* finally, type should be `toeplitz matrix' */
  if (LOGICAL(ans)[0] && e &&
      strcmp(e->type, "toeplitz matrix") != 0)
    LOGICAL(ans)[0] = 0;

  UNPROTECT(2);

  return ans;
}

SEXP toeplitz_rows(SEXP ptr) {
  SEXP tchk;
  SEXP ans = NILSXP;

  /* Perform a type checking */
  PROTECT(tchk = is_tmat(ptr));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e = R_ExternalPtrAddr(ptr);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = toeplitz_nrow(e->matrix);
    UNPROTECT(1);
  } else
    error("pointer provided is not a toeplitz matrix");

  UNPROTECT(1);

  return ans;
}

SEXP toeplitz_cols(SEXP ptr) {
  SEXP tchk;
  SEXP ans = NILSXP;

  /* Perform a type checking */
  PROTECT(tchk = is_tmat(ptr));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e = R_ExternalPtrAddr(ptr);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = toeplitz_ncol(e->matrix);
    UNPROTECT(1);
  } else
    error("pointer provided is not a toeplitz matrix");

  UNPROTECT(1);

  return ans;
}

SEXP tmatmul(SEXP tmat, SEXP v, SEXP transposed) {
  SEXP Y = NILSXP, tchk;

  /* Perform a type checking */
  PROTECT(tchk = is_tmat(tmat));

  if (LOGICAL(tchk)[0]) {
    R_len_t K, L;
    ext_matrix *e;
    toeplitz_matrix *t;

    /* Grab needed data */
    e = R_ExternalPtrAddr(tmat);
    t = e->matrix;

    L = (LOGICAL(transposed)[0] ? toeplitz_ncol(t) : toeplitz_nrow(t));
    K = (LOGICAL(transposed)[0] ? toeplitz_nrow(t) : toeplitz_ncol(t));

    /* Check agains absurd values of inputs */
    if (K != length(v))
      error("invalid length of input vector 'v'");

    /* Allocate output buffer */
    PROTECT(Y = allocVector(REALSXP, L));

    /* Calculate the product */
    if (LOGICAL(transposed)[0])
      toeplitz_tmatmul(REAL(Y), REAL(v), t);
    else
      toeplitz_matmul(REAL(Y), REAL(v), t);

    UNPROTECT(1);
  } else
    error("pointer provided is not a toeplitz matrix");

  UNPROTECT(1);

  return Y;
}

SEXP Lcor(SEXP F, SEXP L) {
  SEXP R = NILSXP;

  R_len_t N = length(F), intL = INTEGER(L)[0];

  if (intL == 0 || intL > N - 1)
    error("invalid length of inpur vector 'F'");

  /* Allocate output buffer */
  PROTECT(R = allocVector(REALSXP, intL));

  calc_Lcor(REAL(R), REAL(F), length(F), intL);

  UNPROTECT(1);

  return R;
}
