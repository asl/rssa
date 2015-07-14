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

#include "fft_plan.h"

typedef struct {
#if HAVE_FFTW3_H
  fftw_complex * circ_freq;
#else
  SEXP circ_freq;
#endif
  R_len_t window;
  R_len_t length;
  fft_plan *fft_plan;
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
}

static void initialize_circulant(toeplitz_matrix *t, fft_plan *f,
                                 const double *R, R_len_t L) {
  R_len_t N = 2*L - 1, i;
  fftw_complex *ocirc;
  double *circ;

  if (!valid_plan(f, N))
    error("invalid FFT plan for given FFT length");

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Fill input buffer */
  for (i = 0; i < L; ++i)
    circ[i] = R[i];

  for (i = 0; i < L-1; ++i)
    circ[L + i] = R[L-i-1];

  /* Run the plan on input data */
  fftw_execute_dft_r2c(f->r2c_plan, circ, ocirc);

  /* Cleanup and return */
  fftw_free(circ);

  t->circ_freq = ocirc;
  t->fft_plan = f;
  t->window = L;
  t->length = N;
}

static void toeplitz_matmul(double* out,
                            const double* v,
                            const void* matrix) {
  const toeplitz_matrix *t = matrix;
  const fft_plan *f = t->fft_plan;
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
  fftw_execute_dft_r2c(f->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < (N/2 + 1); ++i)
    ocirc[i] = ocirc[i] * t->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(f->c2r_plan, ocirc, circ);

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
  const fft_plan *f = t->fft_plan;
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
  fftw_execute_dft_r2c(f->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < (N/2 + 1); ++i)
    ocirc[i] = ocirc[i] * t->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(f->c2r_plan, ocirc, circ);

  /* Cleanup and return */
  for (i = 0; i < K; ++i)
    out[i] = circ[i + L - 1] / N;

  fftw_free(circ);
  fftw_free(ocirc);
}

static void calc_Lcor(double *R, const double *F,
                      R_len_t N, R_len_t L,
                      int circular) {
  R_len_t No, i;
  double *circ;
  fftw_complex *ocirc;
  fftw_plan p1, p2;

  /* Allocate needed memory */
  No = circular ? N : N + L - 1;
  circ = (double*) fftw_malloc(No * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((No/2 + 1) * sizeof(fftw_complex));

  /* Estimate the best plans for given input length */
  p1 = fftw_plan_dft_r2c_1d(No, circ, ocirc, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_c2r_1d(No, ocirc, circ, FFTW_ESTIMATE);

  memcpy(circ, F, N*sizeof(double));
  memset(circ + N, 0, (No - N)*sizeof(double));

  /* Run the plan on input data */
  fftw_execute(p1);

  /* Auto dot-product */
  for (i = 0; i < No/2 + 1; ++i)
    ocirc[i] = ocirc[i] * conj(ocirc[i]);

  fftw_execute(p2);

  /* Return */
  for (i = 0; i < L; ++i)
    R[i] = circ[i] / (circular ? N : N - i) / No;

  /* Cleanup */
  fftw_free(circ);
  fftw_free(ocirc);
  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
}

#else

static void free_circulant(toeplitz_matrix *t) {
  R_ReleaseObject(t->circ_freq);
}

static void initialize_circulant(toeplitz_matrix *t, fft_plan *f,
                                 const double *R, R_len_t L) {
  R_len_t N = 2*L - 1, i;
  Rcomplex *circ;
  SEXP rcirc;

  if (!valid_plan(f, N))
    error("invalid FFT plan for given FFT length");

  /* Allocate needed memory */
  PROTECT(rcirc = allocVector(CPLXSXP, N));
  circ = COMPLEX(rcirc);

  /* Fill input buffer */
  for (i = 0; i < L; ++i) {
    circ[i].r = R[i];
    circ[i].i = 0;
  }

  for (i = 0; i < L - 1; ++i) {
    circ[L + i].r = R[L - i - 1];
    circ[L + i].i = 0;
  }

  /* Run the plan on input data */
  R_PreserveObject(t->circ_freq = eval(lang2(install("fft"), rcirc), R_GlobalEnv));

  /* Cleanup and return */
  UNPROTECT(1);

  t->fft_plan = f;
  t->window = L;
  t->length = N;
}

static void toeplitz_matmul(double* out,
                            const double* v,
                            const void* matrix) {
  const toeplitz_matrix *t = matrix;
  R_len_t N = t->length, L = t->window;
  R_len_t K = N - L + 1, i;
  SEXP rcirc, res, rV1, rTrue;

  /* Allocate needed memory */
  PROTECT(rcirc = allocVector(CPLXSXP, N));
  memset(COMPLEX(rcirc), 0, N * sizeof(Rcomplex));

  /* Fill the arrays */
  for (i = 0; i < K; ++i)
    COMPLEX(rcirc)[i].r = v[i];

  /* Compute the FFT of the vector v */
  PROTECT(rV1 = eval(lang2(install("fft"), rcirc), R_GlobalEnv));

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < N; ++i) {
    Rcomplex x1 = COMPLEX(rV1)[i], x2 = COMPLEX(t->circ_freq)[i];
    COMPLEX(rV1)[i].r = x1.r * x2.r - x1.i * x2.i;
    COMPLEX(rV1)[i].i = x1.r * x2.i + x1.i * x2.r;
  }

  /* Compute the reverse transform to obtain result */
  PROTECT(rTrue = allocVector(LGLSXP, 1));
  LOGICAL(rTrue)[0] = 1;
  PROTECT(res = eval(lang3(install("fft"), rV1, rTrue), R_GlobalEnv));

  /* Cleanup and return */
  for (i = 0; i < L; ++i)
    out[i] = COMPLEX(res)[i].r / N;

  UNPROTECT(4);
}

static void toeplitz_tmatmul(double* out,
                             const double* v,
                             const void* matrix) {
  const toeplitz_matrix *t = matrix;
  R_len_t N = t->length, L = t->window;
  R_len_t K = N - L + 1, i;
  SEXP rcirc, res, rV1, rTrue;

  /* Allocate needed memory */
  PROTECT(rcirc = allocVector(CPLXSXP, N));
  memset(COMPLEX(rcirc), 0, N * sizeof(Rcomplex));

  /* Fill the arrays */
  for (i = 0; i < L; ++i)
    COMPLEX(rcirc)[i + K - 1].r = v[i];

  /* Compute the FFT of the vector v */
  PROTECT(rV1 = eval(lang2(install("fft"), rcirc), R_GlobalEnv));

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < N; ++i) {
    Rcomplex x1 = COMPLEX(rV1)[i], x2 = COMPLEX(t->circ_freq)[i];
    COMPLEX(rV1)[i].r = x1.r * x2.r - x1.i * x2.i;
    COMPLEX(rV1)[i].i = x1.r * x2.i + x1.i * x2.r;
  }

  /* Compute the reverse transform to obtain result */
  PROTECT(rTrue = allocVector(LGLSXP, 1));
  LOGICAL(rTrue)[0] = 1;
  PROTECT(res = eval(lang3(install("fft"), rV1, rTrue), R_GlobalEnv));

  /* Cleanup and return */
  for (i = 0; i < K; ++i)
    out[i] = COMPLEX(res)[i + L - 1].r / N;

  UNPROTECT(4);
}

static void calc_Lcor(double *R, const double *F,
                      R_len_t N, R_len_t L,
                      int circular) {
  R_len_t No = circular ? N : N + L - 1, i;
  SEXP rcirc, res, rV1, rTrue;

  /* Allocate needed memory */
  PROTECT(rcirc = allocVector(CPLXSXP, No));
  memset(COMPLEX(rcirc), 0, No * sizeof(Rcomplex));

  for (i = 0; i < N; ++i)
    COMPLEX(rcirc)[i].r = F[i];

  /* Compute the FFT of the padded input vector */
  PROTECT(rV1 = eval(lang2(install("fft"), rcirc), R_GlobalEnv));

  /* Auto dot-product */
  for (i = 0; i < No; ++i) {
    Rcomplex x = COMPLEX(rV1)[i];
    COMPLEX(rV1)[i].r = x.r * x.r + x.i * x.i;
    COMPLEX(rV1)[i].i = 0;
  }

  /* Compute the reverse transform to obtain result */
  PROTECT(rTrue = allocVector(LGLSXP, 1));
  LOGICAL(rTrue)[0] = 1;
  PROTECT(res = eval(lang3(install("fft"), rV1, rTrue), R_GlobalEnv));

  /* Return */
  for (i = 0; i < L; ++i)
    R[i] = COMPLEX(res)[i].r / (circular ? N : N - i) / No;

  /* Cleanup */
  UNPROTECT(4);
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

SEXP initialize_tmat(SEXP R, SEXP fft_plan) {
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
  initialize_circulant(t, R_ExternalPtrAddr(fft_plan), REAL(R), L);
  e->matrix = t;

  /* Make an external pointer envelope */
  tmat = R_MakeExternalPtr(e, install("external matrix"), fft_plan);
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

SEXP Lcor(SEXP F, SEXP L, SEXP circular) {
  SEXP R = NILSXP;

  R_len_t N = length(F), intL = INTEGER(L)[0];

  if (intL <= 0 || intL > N)
    error("invalid length of input vector 'F'");

  /* Allocate output buffer */
  PROTECT(R = allocVector(REALSXP, intL));

  calc_Lcor(REAL(R), REAL(F), length(F), intL, LOGICAL(circular)[0]);

  UNPROTECT(1);

  return R;
}
