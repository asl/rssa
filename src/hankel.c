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
#include <R_ext/Complex.h>
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
} hankel_matrix;

static unsigned hankel_nrow(const void *matrix) {
  const hankel_matrix *h = matrix;
  return h->window;
}

static unsigned hankel_ncol(const void *matrix) {
  const hankel_matrix *h = matrix;
  return h->length - h->window + 1;
}

#if HAVE_FFTW3_H
static void free_plan(fft_plan *f) {
  fftw_destroy_plan(f->r2c_plan);
  fftw_destroy_plan(f->c2r_plan);
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

  f->N = N;

  fftw_free(circ);
  fftw_free(ocirc);
}

static void free_circulant(hankel_matrix *h) {
  fftw_free(h->circ_freq);
}

static void initialize_circulant(hankel_matrix *h, fft_plan *f,
                                 const double *F, R_len_t N, R_len_t L) {
  R_len_t K = N - L + 1, i;
  fftw_complex *ocirc;
  double *circ;

  if (!valid_plan(f, N))
    error("invalid FFT plan for given FFT length");

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Fill input buffer */
  for (i = K-1; i < N; ++i)
    circ[i - K + 1] = F[i];

  for (i = 0; i < K-1; ++i) {
    circ[L + i] = F[i];
  }

  /* Run the plan on input data */
  fftw_execute_dft_r2c(f->r2c_plan, circ, ocirc);

  /* Cleanup and return */
  fftw_free(circ);

  h->circ_freq = ocirc;
  h->window = L;
  h->length = N;
  h->fft_plan = f;
}

static void hankel_matmul(double* out,
                          const double* v,
                          const void* matrix) {
  const hankel_matrix *h = matrix;
  const fft_plan *f = h->fft_plan;
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
  fftw_execute_dft_r2c(f->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < (N/2 + 1); ++i)
    ocirc[i] = ocirc[i] * h->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(f->c2r_plan, ocirc, circ);

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
  const fft_plan *f = h->fft_plan;
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
  fftw_execute_dft_r2c(f->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < (N/2 + 1); ++i)
    ocirc[i] = ocirc[i] * h->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(f->c2r_plan, ocirc, circ);

  /* Cleanup and return */
  for (i = 0; i < K; ++i)
    out[i] = circ[i + L - 1] / N;

  fftw_free(circ);
  fftw_free(ocirc);
}

static R_INLINE void hankelize_fft(double *F,
                                   const double *U, const double *V,
                                   R_len_t L, R_len_t K,
                                   const fft_plan *f) {
  R_len_t N = K + L - 1;
  R_len_t i;
  double *iU, *iV;
  fftw_complex *cU, *cV;

  if (!valid_plan(f, N))
    error("invalid FFT plan for given FFT length");

  /* Allocate needed memory */
  iU = (double*) fftw_malloc(N * sizeof(double));
  iV = (double*) fftw_malloc(N * sizeof(double));
  cU = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));
  cV = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Fill in buffers */
  memcpy(iU, U, L*sizeof(double));
  memset(iU + L, 0, (K - 1)*sizeof(double));

  memcpy(iV, V, K*sizeof(double));
  memset(iV + K, 0, (L - 1)*sizeof(double));

  /* Compute the FFTs */
  fftw_execute_dft_r2c(f->r2c_plan, iU, cU);
  fftw_execute_dft_r2c(f->r2c_plan, iV, cV);

  /* Dot-multiply */
  for (i = 0; i < N/2 + 1; ++i)
    cU[i] = cU[i] * cV[i];

  /* Compute the inverse FFT */
  fftw_execute_dft_c2r(f->c2r_plan, cU, iU);

  /* Form the result */
  for (i = 0; i < N; ++i) {
    R_len_t leftu, rightu, l;

    if (i < L) leftu = i; else leftu = L - 1;
    if (i < K) rightu = 0; else  rightu = i - K + 1;

    l = (leftu - rightu + 1);

    F[i] = iU[i] / l / N;
  }

  fftw_free(iU);
  fftw_free(iV);
  fftw_free(cU);
  fftw_free(cV);
}

static void compute_L_covariation_matrix_first_row(const double *F, R_len_t N, R_len_t L,
                                                   double *R,
                                                   const fft_plan *f) {
  double *oF;
  fftw_complex *iF, *iFc;
  R_len_t i, K;

  /* Length's check */
  if ((L >= N) || (L < 1))
    error("must be 'N' > 'L' >= 1");

  K = N - L + 1;

  /* Allocate needed memory */
  oF = (double*) fftw_malloc(N * sizeof(double));
  iF = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));
  iFc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  memcpy(oF, F, N*sizeof(double));

  fftw_execute_dft_r2c(f->r2c_plan, oF, iF);

  memcpy(oF, F, K*sizeof(double));
  memset(oF + K, 0, (N - K)*sizeof(double));

  fftw_execute_dft_r2c(f->r2c_plan, oF, iFc);

  /* Dot-product */
  for (i = 0; i < N/2 + 1; ++i)
    iF[i] = iF[i] * conj(iFc[i]);

  fftw_execute_dft_c2r(f->c2r_plan, iF, oF);

  /* Return */
  for (i = 0; i < L; ++i)
    R[i] = oF[i] / N;

  /* Cleanup */
  fftw_free(oF);
  fftw_free(iF);
  fftw_free(iFc);
}
#else
static void free_plan(fft_plan *f) {
  (void)f;
}

static void initialize_plan(fft_plan *f, R_len_t N) {
  f->N = N;
}

static void free_circulant(hankel_matrix *h) {
  R_ReleaseObject(h->circ_freq);
}

static void initialize_circulant(hankel_matrix *h, fft_plan *f,
                                 const double *F, R_len_t N, R_len_t L) {
  R_len_t K = N - L + 1, i;
  Rcomplex* circ;
  SEXP rcirc;

  if (!valid_plan(f, N))
    error("invalid FFT plan for given FFT length");

  /* Allocate needed memory */
  PROTECT(rcirc = allocVector(CPLXSXP, N));
  circ = COMPLEX(rcirc);

  /* Fill input buffer */
  for (i = K - 1; i < N; ++i) {
    circ[i - K + 1].r = F[i];
    circ[i - K + 1].i = 0;
  }

  for (i = 0; i < K - 1; ++i) {
    circ[L + i].r = F[i];
    circ[L + i].i = 0;
  }

  /* Run the FFT on input data */
  R_PreserveObject(h->circ_freq = eval(lang2(install("fft"), rcirc), R_GlobalEnv));

  UNPROTECT(1);

  h->window = L;
  h->length = N;
  h->fft_plan = f;
}

static void hankel_matmul(double* out,
                          const double* v,
                          const void* matrix) {
  const hankel_matrix *h = matrix;
  R_len_t N = h->length, L = h->window;
  R_len_t K = N - L + 1, i;
  SEXP rcirc, rV1, res, rTrue;

  /* Allocate needed memory */
  PROTECT(rcirc = allocVector(CPLXSXP, N));
  memset(COMPLEX(rcirc), 0, N * sizeof(Rcomplex));

  /* Fill the arrays */
  for (i = 0; i < K; ++i)
    COMPLEX(rcirc)[i].r = v[K - i - 1];

  /* Compute the FFT of the reversed vector v */
  PROTECT(rV1 = eval(lang2(install("fft"), rcirc), R_GlobalEnv));

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < N; ++i) {
    Rcomplex x1 = COMPLEX(rV1)[i], x2 = COMPLEX(h->circ_freq)[i];
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

static void hankel_tmatmul(double* out,
                           const double* v,
                           const void* matrix) {
  const hankel_matrix *h = matrix;
  R_len_t N = h->length, L = h->window;
  R_len_t K = N - L + 1, i;
  SEXP rcirc, rV1, res, rTrue;

  /* Allocate needed memory */
  PROTECT(rcirc = allocVector(CPLXSXP, N));
  memset(COMPLEX(rcirc), 0, N * sizeof(Rcomplex));

  /* Fill the arrays */
  for (i = 0; i < L; ++i) {
    COMPLEX(rcirc)[i + K - 1].r = v[L - i - 1];
    COMPLEX(rcirc)[i + K - 1].i = 0;
  }

  /* Compute the FFT of the reversed vector v */
  PROTECT(rV1 = eval(lang2(install("fft"), rcirc), R_GlobalEnv));

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  /* FIXME! => rV1 */
  for (i = 0; i < N; ++i) {
    Rcomplex x1 = COMPLEX(rV1)[i], x2 = COMPLEX(h->circ_freq)[i];
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

static R_INLINE void hankelize_fft(double *F,
                                   const double *U, const double *V,
                                   R_len_t L, R_len_t K,
                                   const fft_plan *f) {
  R_len_t N = K + L - 1;
  R_len_t i;
  SEXP rU, rU1, rV, rV1, res, rTrue;

  if (!valid_plan(f, N))
    error("invalid FFT plan for given FFT length");

  /* Allocate needed memory */
  PROTECT(rU = allocVector(CPLXSXP, N));
  memset(COMPLEX(rU), 0, N * sizeof(Rcomplex));
  PROTECT(rV = allocVector(CPLXSXP, N));
  memset(COMPLEX(rV), 0, N * sizeof(Rcomplex));

  /* Fill in buffers */
  for (i = 0; i < L; ++i)
    COMPLEX(rU)[i].r = U[i];

  for (i = 0; i < K; ++i)
    COMPLEX(rV)[i].r = V[i];

  /* Compute the FFTs */
  PROTECT(rU1 = eval(lang2(install("fft"), rU), R_GlobalEnv));
  PROTECT(rV1 = eval(lang2(install("fft"), rV), R_GlobalEnv));

  /* Dot-multiply */
  for (i = 0; i < N; ++i) {
    Rcomplex x1 = COMPLEX(rU1)[i], x2 = COMPLEX(rV1)[i];
    COMPLEX(rU1)[i].r = x1.r * x2.r - x1.i * x2.i;
    COMPLEX(rU1)[i].i = x1.r * x2.i + x1.i * x2.r;
  }

  /* Compute the inverse FFT */
  PROTECT(rTrue = allocVector(LGLSXP, 1));
  LOGICAL(rTrue)[0] = 1;
  PROTECT(res = eval(lang3(install("fft"), rU1, rTrue), R_GlobalEnv));

  /* Form the result */
  for (i = 0; i < N; ++i) {
    R_len_t leftu, rightu, l;

    if (i < L) leftu = i; else leftu = L - 1;
    if (i < K) rightu = 0; else  rightu = i - K + 1;

    l = (leftu - rightu + 1);

    F[i] = COMPLEX(res)[i].r / l / N;
  }

  UNPROTECT(6);
}

static void compute_L_covariation_matrix_first_row(const double *F, R_len_t N, R_len_t L,
                                                   double *R,
                                                   const fft_plan *f) {
  R_len_t K = N - L + 1;
  R_len_t i;
  SEXP rF, rFc, rF1, rFc1, rTrue, res;

  if (!valid_plan(f, N))
    error("invalid FFT plan for given FFT length");

  /* Length check */
  if ((L >= N) || (L < 1))
    error("must be 'N' > 'L' >= 1");

  /* Allocate needed memory */
  PROTECT(rF = allocVector(CPLXSXP, N));
  memset(COMPLEX(rF), 0, N * sizeof(Rcomplex));
  PROTECT(rFc = allocVector(CPLXSXP, N));
  memset(COMPLEX(rFc), 0, N * sizeof(Rcomplex));

  /* Fill in buffers */
  for (i = 0; i < N; ++i)
    COMPLEX(rF)[i].r = F[i];

  for (i = 0; i < K; ++i)
    COMPLEX(rFc)[i].r = F[i];

  /* Compute the FFTs */
  PROTECT(rF1 = eval(lang2(install("fft"), rF), R_GlobalEnv));
  PROTECT(rFc1 = eval(lang2(install("fft"), rFc), R_GlobalEnv));

  /* Dot-product */
  for (i = 0; i < N; ++i) {
    Rcomplex x1 = COMPLEX(rF1)[i], x2 = COMPLEX(rFc1)[i];
    COMPLEX(rF1)[i].r = x1.r * x2.r + x1.i * x2.i;
    COMPLEX(rF1)[i].i = -x1.r * x2.i + x1.i * x2.r;
  }

  /* Compute the inverse FFT */
  PROTECT(rTrue = allocVector(LGLSXP, 1));
  LOGICAL(rTrue)[0] = 1;
  PROTECT(res = eval(lang3(install("fft"), rF1, rTrue), R_GlobalEnv));

  /* Form the result */
  for (i = 0; i < L; ++i)
    R[i] = COMPLEX(res)[i].r / N;

  UNPROTECT(6);
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

static void hmat_finalizer(SEXP ptr) {
  ext_matrix *e;
  hankel_matrix *h;

  if (TYPEOF(ptr) != EXTPTRSXP)
    return;

  e = R_ExternalPtrAddr(ptr);
  if (!e)
    return;

  if (strcmp(e->type, "hankel matrix"))
    return;

  h = e->matrix;

  free_circulant(h);
  Free(h);

  Free(e);
  R_ClearExternalPtr(ptr);
}

SEXP initialize_hmat(SEXP F, SEXP window, SEXP fft_plan) {
  /* Perform a type checking */
  if (!LOGICAL(is_fft_plan(fft_plan))[0]) {
    error("pointer provided is not a fft plan");
    return NILSXP;
  }

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

  initialize_circulant(h, R_ExternalPtrAddr(fft_plan), REAL(F), N, L);
  e->matrix = h;

  /* Make an external pointer envelope */
  hmat = R_MakeExternalPtr(e, install("external matrix"), fft_plan);
  R_RegisterCFinalizer(hmat, hmat_finalizer);

  return hmat;
}

SEXP is_hmat(SEXP ptr) {
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

  /* finally, type should be `hankel matrix' */
  if (LOGICAL(ans)[0] && e &&
      strcmp(e->type, "hankel matrix") != 0)
    LOGICAL(ans)[0] = 0;

  UNPROTECT(2);

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

    L = (LOGICAL(transposed)[0] ? hankel_ncol(h) : hankel_nrow(h));
    K = (LOGICAL(transposed)[0] ? hankel_nrow(h) : hankel_ncol(h));

    /* Check agains absurd values of inputs */
    if (K != length(v))
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

SEXP hankelize_one_fft(SEXP U, SEXP V, SEXP fftplan) {
  /* Perform a type checking */
  if (!LOGICAL(is_fft_plan(fftplan))[0]) {
    error("pointer provided is not a fft plan");
    return NILSXP;
  }

  SEXP F = NILSXP;
  fft_plan *plan;
  double *rU = REAL(U), *rV = REAL(V), *rF;
  R_len_t L, K, N;

  L = length(U); K = length(V);
  N = K + L - 1;

  /* Grab needed data */
  plan = R_ExternalPtrAddr(fftplan);

  /* Allocate buffer for output */
  PROTECT(F = allocVector(REALSXP, N));
  rF = REAL(F);

  /* Perform the actual hankelization */
  hankelize_fft(rF, rU, rV, L, K, plan);

  UNPROTECT(1);

  return F;
}

SEXP hankelize_multi_fft(SEXP U, SEXP V, SEXP fftplan) {
  /* Perform a type checking */
  if (!LOGICAL(is_fft_plan(fftplan))[0]) {
    error("pointer provided is not a fft plan");
    return NILSXP;
  }

  double *rU = REAL(U), *rV = REAL(V), *rF;
  R_len_t L, K, N, i, count;
  SEXP F = NILSXP;
  fft_plan *plan;
  int *dimu, *dimv;

  /* Calculate length of inputs and output */
  dimu = INTEGER(getAttrib(U, R_DimSymbol));
  dimv = INTEGER(getAttrib(V, R_DimSymbol));
  L = dimu[0]; K = dimv[0];
  if ((count = dimu[1]) != dimv[1])
    error("Both 'U' and 'V' should have equal number of columns");
  N = K + L - 1;

  /* Grab needed data */
  plan = R_ExternalPtrAddr(fftplan);

  /* Allocate buffer for output */
  PROTECT(F = allocMatrix(REALSXP, N, count));
  rF = REAL(F);

  /* Perform the actual hankelization */
  for (i = 0; i < count; ++i) {
    R_CheckUserInterrupt();
    /* TODO: nice target for OpenMP stuff */
    hankelize_fft(rF+i*N, rU+i*L, rV+i*K, L, K, plan);
  }

  UNPROTECT(1);
  return F;
}

SEXP Lcov_matrix(SEXP F, SEXP L, SEXP fft_plan) {
  /* Perform a type checking */
  if (!LOGICAL(is_fft_plan(fft_plan))[0]) {
    error("pointer provided is not a fft plan");
    return NILSXP;
  }

  R_len_t intL = INTEGER(L)[0];
  R_len_t i, j;
  R_len_t K = length(F) - intL + 1;
  SEXP ans;

  PROTECT(ans = allocMatrix(REALSXP, intL, intL));
  double *rans = REAL(ans);
  double *pF = REAL(F);
  compute_L_covariation_matrix_first_row(REAL(F), length(F), intL, rans, R_ExternalPtrAddr(fft_plan));

  for (j = 1; j < intL; ++j)
    rans[intL*j] = rans[j];

  for (j = 1; j < intL; ++j)
    for (i = j; i < intL; ++i)
      rans[i + intL * j] = rans[j + intL * i] = rans[(i - 1) + (j - 1) * intL] - pF[i - 1] * pF[j - 1] + pF[i + K - 1] * pF[j + K - 1];

  UNPROTECT(1);
  return ans;
}
