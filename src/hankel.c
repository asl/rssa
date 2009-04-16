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

static unsigned hankel_nrow(const void *matrix) {
  const hankel_matrix *h = matrix;
  return h->window;
}

static unsigned hankel_ncol(const void *matrix) {
  const hankel_matrix *h = matrix;
  return h->length - h->window + 1;
}

static void free_circulant(hankel_matrix *h) {
  fftw_free(h->circ_freq);
  fftw_destroy_plan(h->r2c_plan);
  fftw_destroy_plan(h->c2r_plan);
}

static void initialize_circulant(hankel_matrix *h,
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

static void hankel_matmul(double* out,
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

static void hankel_tmatmul(double* out,
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

/* This is just direct R-to-C translation  */
static R_INLINE void hankelize(double *F,
                               double *U, double *V,
                               R_len_t L, R_len_t K) {
  R_len_t i, N = K + L - 1;
  R_len_t leftu, rightu, leftv, rightv, l, j;

  for (i = 0; i < N; ++i) {
    double s = 0;

    if (i < L) {
      leftu = i;
      leftv = 0;
    } else {
      leftu = L - 1;
      leftv = i - L + 1;
    }
    if (i < K) {
      rightu = 0;
      rightv = i;
    } else {
      rightu = i - K + 1;
      rightv = K - 1;
    }

    l = leftu - rightu + 1;

    for (j = 0; j < l; ++j) {
      s += U[leftu - j] * V[leftv + j];
    }

    F[i] = s / (double) l;
  }
}

static R_INLINE void hankelize_fft(double *F,
                                   double *U, double *V,
                                   size_t L, size_t K) {
  size_t i, N = K + L - 1;
  double *iU, *iV;
  fftw_complex *cU, *cV;
  fftw_plan p1, p2;

  /* Allocate needed memory */
  iU = (double*) fftw_malloc(N * sizeof(double));
  iV = (double*) fftw_malloc(N * sizeof(double));
  cU = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));
  cV = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Estimate the best plans for given input length */
  p1 = fftw_plan_dft_r2c_1d(N, iU, cU, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_c2r_1d(N, cU, iU, FFTW_ESTIMATE);

  /* Fill in buffers */
  memcpy(iU, U, L*sizeof(double));
  memset(iU+L, 0, (K - 1)*sizeof(double));

  memcpy(iV, V, K*sizeof(double));
  memset(iV+K, 0, (L - 1)*sizeof(double));

  /* Compute the FFTs */
  fftw_execute_dft_r2c(p1, iU, cU);
  fftw_execute_dft_r2c(p1, iV, cV);

  /* Dot-multiply */
  for (i = 0; i < N/2 + 1; ++i)
    cU[i] = cU[i] * cV[i];

  /* Compute the inverse FFT */
  fftw_execute_dft_c2r(p2, cU, iU);

  /* Form the result */
  for (i = 0; i < N; ++i) {
    size_t leftu, rightu, l;

    if (i < L) leftu = i; else leftu = L - 1;
    if (i < K) rightu = 0; else  rightu = i - K + 1;

    l = N*(leftu - rightu + 1);

    F[i] = iU[i] / l;
  }

  fftw_free(iU);
  fftw_free(iV);
  fftw_free(cU);
  fftw_free(cV);
  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
}

static void hmat_finalizer(SEXP ptr) {
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

  free_circulant(h);
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
  e->mulfn = hankel_matmul;
  e->tmulfn = hankel_tmatmul;
  e->ncol = hankel_ncol;
  e->nrow = hankel_nrow;

  /* Build toeplitz circulants for hankel matrix */
  h = Calloc(1, hankel_matrix);
  initialize_circulant(h, REAL(F), N, L);
  e->matrix = h;

  /* Make an external pointer envelope */
  hmat = R_MakeExternalPtr(e, install("external matrix"), R_NilValue);
  R_RegisterCFinalizer(hmat, hmat_finalizer);

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
    INTEGER(ans)[0] = hankel_nrow(e->matrix);
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
    INTEGER(ans)[0] = hankel_ncol(e->matrix);
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
      hankel_tmatmul(REAL(Y), REAL(v), h);
    else
      hankel_matmul(REAL(Y), REAL(v), h);

    UNPROTECT(1);
  } else
    error("pointer provided is not a hankel matrix");

  UNPROTECT(1);

  return Y;
}

SEXP hankelize_one(SEXP U, SEXP V) {
  double *rU = REAL(U), *rV = REAL(V), *rF;
  R_len_t L, K;
  SEXP F;

  /* Calculate length of inputs */
  L = length(U); K = length(V);

  /* Allocate buffer for output */
  PROTECT(F = allocVector(REALSXP, K + L - 1));
  rF = REAL(F);

  /* Perform the actual hankelization */
  hankelize(rF, rU, rV, L, K);

  UNPROTECT(1);
  return F;
}

SEXP hankelize_one_fft(SEXP U, SEXP V) {
  double *rU = REAL(U), *rV = REAL(V), *rF;
  R_len_t L, K;
  SEXP F;

  /* Calculate length of inputs */
  L = length(U); K = length(V);

  /* Allocate buffer for output */
  PROTECT(F = allocVector(REALSXP, K + L - 1));
  rF = REAL(F);

  /* Perform the actual hankelization */
  hankelize_fft(rF, rU, rV, L, K);

  UNPROTECT(1);
  return F;
}

SEXP hankelize_multi(SEXP U, SEXP V) {
  double *rU = REAL(U), *rV = REAL(V), *rF;
  R_len_t L, K, N, i, count;
  SEXP F;
  int *dimu, *dimv;

  /* Calculate length of inputs and output */
  dimu = INTEGER(getAttrib(U, R_DimSymbol));
  dimv = INTEGER(getAttrib(V, R_DimSymbol));
  L = dimu[0]; K = dimv[0];
  if ((count = dimu[1]) != dimv[1])
    error("Both 'U' and 'V' should have equal number of columns");
  N = K + L - 1;

  /* Allocate buffer for output */
  PROTECT(F = allocMatrix(REALSXP, N, count));
  rF = REAL(F);

  /* Perform the actual hankelization */
  for (i = 0; i < count; ++i) {
    R_CheckUserInterrupt();
    /* TODO: nice target for OpenMP stuff */
    hankelize(rF+i*N, rU+i*L, rV+i*K, L, K);
  }

  UNPROTECT(1);
  return F;
}
