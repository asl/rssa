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
  R_len_t factor;
  R_len_t length;
  fft_plan *fft_plan;
} hankel_matrix;

static unsigned hankel_nrow(const void *matrix) {
  const hankel_matrix *h = matrix;
  return h->fft_plan->col_ind != NULL ? h->fft_plan->col_ind->num : h->window;
}

static unsigned hankel_ncol(const void *matrix) {
  const hankel_matrix *h = matrix;
  return h->fft_plan->row_ind != NULL ? h->fft_plan->row_ind->num : h->factor;
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
                                 const double *F, R_len_t N, R_len_t L,
                                 const int *circular) {
  fftw_complex *ocirc;
  double *circ;

  if (!valid_plan(f, N))
    error("invalid FFT plan for given FFT length");

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Fill input buffer */
  memcpy(circ, F, N * sizeof(double));

  /* Run the plan on input data */
  fftw_execute_dft_r2c(f->r2c_plan, circ, ocirc);

  /* Cleanup and return */
  fftw_free(circ);

  h->circ_freq = ocirc;
  h->window = L;
  h->length = N;
  h->factor = circular[0] ? N : N - L + 1;
  h->fft_plan = f;
}

static void hankel_matmul(double* out,
                          const double* v,
                          const void* matrix) {
  const hankel_matrix *h = matrix;
  const fft_plan *f = h->fft_plan;
  R_len_t N = h->length, L = h->window;
  R_len_t K = h->factor, i;
  double *circ;
  fftw_complex *ocirc;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Fill the arrays */
  memset(circ, 0, N * sizeof(double));
  if (f->row_ind == NULL) {
    memcpy(circ, v, K * sizeof(double));
  } else {
    for (i = 0; i < f->row_ind->num; ++i) {
      circ[f->row_ind->ind[i]] = v[i];
    }
  }

  /* Compute the FFT of the reversed vector v */
  fftw_execute_dft_r2c(f->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < (N/2 + 1); ++i)
    ocirc[i] = conj(ocirc[i]) * h->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(f->c2r_plan, ocirc, circ);

  /* Cleanup and return */
  if (f->col_ind == NULL) {
    for (i = 0; i < L; ++i) {
      out[i] = circ[i] / N;
    }
  } else {
    for (i = 0; i < f->col_ind->num; ++i) {
      out[i] = circ[f->col_ind->ind[i]] / N;
    }
  }

  fftw_free(circ);
  fftw_free(ocirc);
}

static void hankel_tmatmul(double* out,
                           const double* v,
                           const void* matrix) {
  const hankel_matrix *h = matrix;
  const fft_plan *f = h->fft_plan;
  R_len_t N = h->length, L = h->window;
  R_len_t K = h->factor, i;
  double *circ;
  fftw_complex *ocirc;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(N * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));

  /* Fill the arrays */
  memset(circ, 0, N * sizeof(double));
  if (f->col_ind == NULL) {
    memcpy(circ, v, L * sizeof(double));
  } else {
    for (i = 0; i < f->col_ind->num; ++i) {
      circ[f->col_ind->ind[i]] = v[i];
    }
  }

  /* Compute the FFT of the reversed vector v */
  fftw_execute_dft_r2c(f->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < (N/2 + 1); ++i)
    ocirc[i] = conj(ocirc[i]) * h->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(f->c2r_plan, ocirc, circ);

  /* Cleanup and return */
  if (f->row_ind == NULL) {
    for (i = 0; i < K; ++i) {
      out[i] = circ[i] / N;
    }
  } else {
    for (i = 0; i < f->row_ind->num; ++i) {
      out[i] = circ[f->row_ind->ind[i]] / N;
    }
  }

  fftw_free(circ);
  fftw_free(ocirc);
}

static R_INLINE void hankelize_fft(double *F,
                                   const double *U, const double *V,
                                   R_len_t L, R_len_t K,
                                   const fft_plan *f) {
  R_len_t N = f->N, i;
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
  memset(iU, 0, N * sizeof(double));
  if (f->col_ind == NULL) {
    memcpy(iU, U, L * sizeof(double));
  } else {
    for (i = 0; i < f->col_ind->num; ++i) {
      iU[f->col_ind->ind[i]] = U[i];
    }
  }

  memset(iV, 0, N * sizeof(double));
  if (f->row_ind == NULL) {
    memcpy(iV, V, K * sizeof(double));
  } else {
    for (i = 0; i < f->row_ind->num; ++i) {
      iV[f->row_ind->ind[i]] = V[i];
    }
  }

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
    if (f->weights[i]) {
      F[i] = iU[i] / N / f->weights[i];
    } else {
      F[i] = NA_REAL;
    }
  }

  fftw_free(iU);
  fftw_free(iV);
  fftw_free(cU);
  fftw_free(cV);
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
                                 const double *F, R_len_t N, R_len_t L,
                                 const int *circular) {
  R_len_t i;
  Rcomplex* circ;
  SEXP rcirc;

  if (!valid_plan(f, N))
    error("invalid FFT plan for given FFT length");

  /* Allocate needed memory */
  PROTECT(rcirc = allocVector(CPLXSXP, N));
  circ = COMPLEX(rcirc);

  /* Fill input buffer */
  for (i = 0; i < N; ++i) {
    circ[i].r = F[i];
    circ[i].i = 0;
  }

  /* Run the FFT on input data */
  R_PreserveObject(h->circ_freq = eval(lang2(install("fft"), rcirc), R_GlobalEnv));

  UNPROTECT(1);

  h->window = L;
  h->length = N;
  h->factor = circular[0] ? N : N - L + 1;
  h->fft_plan = f;
}

static void hankel_matmul(double* out,
                          const double* v,
                          const void* matrix) {
  const hankel_matrix *h = matrix;
  const fft_plan *f = h->fft_plan;
  R_len_t N = h->length, L = h->window;
  R_len_t K = h->factor, i;
  SEXP rcirc, rV1, res, rTrue;

  /* Allocate needed memory */
  PROTECT(rcirc = allocVector(CPLXSXP, N));
  memset(COMPLEX(rcirc), 0, N * sizeof(Rcomplex));

  /* Fill the arrays */
  if (f->row_ind == NULL) {
    for (i = 0; i < K; ++i)
      COMPLEX(rcirc)[i].r = v[i];
  } else {
    for (i = 0; i < f->row_ind->num; ++i) {
      COMPLEX(rcirc)[f->row_ind->ind[i]].r = v[i];
    }
  }

  /* Compute the FFT of the reversed vector v */
  PROTECT(rV1 = eval(lang2(install("fft"), rcirc), R_GlobalEnv));

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < N; ++i) {
    Rcomplex x1 = COMPLEX(rV1)[i], x2 = COMPLEX(h->circ_freq)[i];
    COMPLEX(rV1)[i].r = x1.r * x2.r + x1.i * x2.i;
    COMPLEX(rV1)[i].i = x1.r * x2.i - x1.i * x2.r;
  }

  /* Compute the reverse transform to obtain result */
  PROTECT(rTrue = allocVector(LGLSXP, 1));
  LOGICAL(rTrue)[0] = 1;
  PROTECT(res = eval(lang3(install("fft"), rV1, rTrue), R_GlobalEnv));

  /* Cleanup and return */
  if (f->col_ind == NULL) {
    for (i = 0; i < L; ++i)
      out[i] = COMPLEX(res)[i].r / N;
  } else {
    for (i = 0; i < f->col_ind->num; ++i) {
      out[i] = COMPLEX(res)[f->col_ind->ind[i]].r / N;
    }
  }

  UNPROTECT(4);
}

static void hankel_tmatmul(double* out,
                           const double* v,
                           const void* matrix) {
  const hankel_matrix *h = matrix;
  const fft_plan *f = h->fft_plan;
  R_len_t N = h->length, L = h->window;
  R_len_t K = h->factor, i;
  SEXP rcirc, rV1, res, rTrue;

  /* Allocate needed memory */
  PROTECT(rcirc = allocVector(CPLXSXP, N));
  memset(COMPLEX(rcirc), 0, N * sizeof(Rcomplex));

  /* Fill the arrays */
  if (f->col_ind == NULL) {
    for (i = 0; i < L; ++i)
      COMPLEX(rcirc)[i].r = v[i];
  } else {
    for (i = 0; i < f->col_ind->num; ++i) {
      COMPLEX(rcirc)[f->col_ind->ind[i]].r = v[i];
    }
  }

  /* Compute the FFT of the reversed vector v */
  PROTECT(rV1 = eval(lang2(install("fft"), rcirc), R_GlobalEnv));

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  /* FIXME! => rV1 */
  for (i = 0; i < N; ++i) {
    Rcomplex x1 = COMPLEX(rV1)[i], x2 = COMPLEX(h->circ_freq)[i];
    COMPLEX(rV1)[i].r = x1.r * x2.r + x1.i * x2.i;
    COMPLEX(rV1)[i].i = x1.r * x2.i - x1.i * x2.r;
  }

  /* Compute the reverse transform to obtain result */
  PROTECT(rTrue = allocVector(LGLSXP, 1));
  LOGICAL(rTrue)[0] = 1;
  PROTECT(res = eval(lang3(install("fft"), rV1, rTrue), R_GlobalEnv));

  /* Cleanup and return */
  if (f->row_ind == NULL) {
    for (i = 0; i < K; ++i)
      out[i] = COMPLEX(res)[i].r / N;
  } else {
    for (i = 0; i < f->row_ind->num; ++i) {
      out[i] = COMPLEX(res)[f->row_ind->ind[i]].r / N;
    }
  }

  UNPROTECT(4);
}

static R_INLINE void hankelize_fft(double *F,
                                   const double *U, const double *V,
                                   R_len_t L, R_len_t K,
                                   const fft_plan *f) {
  R_len_t N = f->N;
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
  if (f->col_ind == NULL) {
    for (i = 0; i < L; ++i)
      COMPLEX(rU)[i].r = U[i];
  } else {
    for (i = 0; i < f->col_ind->num; ++i) {
      COMPLEX(rU)[f->col_ind->ind[i]].r = U[i];
    }
  }

  if (f->row_ind == NULL) {
    for (i = 0; i < K; ++i)
      COMPLEX(rV)[i].r = V[i];
  } else {
    for (i = 0; i < f->row_ind->num; ++i) {
      COMPLEX(rV)[f->row_ind->ind[i]].r = V[i];
    }
  }

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
    if (f->weights[i]) {
      F[i] = COMPLEX(res)[i].r / N / f->weights[i];
    } else {
      F[i] = NA_REAL;
    }
  }

  UNPROTECT(6);
}

#endif

static area_indices *alloc_area1d(SEXP mask) {
  if (mask == R_NilValue) {
    return NULL;
  }
  area_indices *area = Calloc(1, area_indices);
  int *maskValues = LOGICAL(mask);
  R_len_t max_ind = length(mask);

  /* Count the number of nonzero elements and allocate the arrays */
  R_len_t ind;
  area->num = 0;
  for (ind = 0; ind < max_ind; ++ind) {
    area->num += maskValues[ind];
  }

  area->ind = Calloc(area->num, R_len_t);

  /* Fill in the arrays of indices */
  R_len_t k;
  for (ind = 0, k = 0; ind < max_ind; ++ind) {
    if (maskValues[ind]) {
      area->ind[k] = ind;
      ++k;
    }
  }

  return area;
}

static void fft_plan_finalizer(SEXP ptr) {
  fft_plan *f;

  if (TYPEOF(ptr) != EXTPTRSXP)
    return;

  f = R_ExternalPtrAddr(ptr);
  if (!f)
    return;

  free_plan(f);
  free_area(f->col_ind);
  free_area(f->row_ind);
  Free(f->weights);

  Free(f);

  R_ClearExternalPtr(ptr);
}

SEXP initialize_fft_plan(SEXP rN, SEXP wmask, SEXP fmask, SEXP weights) {
  R_len_t N;
  fft_plan *f;
  SEXP res;

  N = INTEGER(rN)[0];

  /* Allocate memory */
  f = Calloc(1, fft_plan);

  /* Do actual plan initialization */
  initialize_plan(f, N);

  /* TODO: add a check for correct window sizes */
  f->col_ind = alloc_area1d(wmask);
  f->row_ind = alloc_area1d(fmask);
  f->weights = alloc_weights(weights);

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

SEXP initialize_hmat(SEXP F, SEXP window, SEXP circular, SEXP fft_plan) {
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

  initialize_circulant(h, R_ExternalPtrAddr(fft_plan), REAL(F), N, L, LOGICAL(circular));
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

  /* Grab needed data */
  plan = R_ExternalPtrAddr(fftplan);
  N = plan->N;

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

  /* Grab needed data */
  plan = R_ExternalPtrAddr(fftplan);
  N = plan->N;

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
