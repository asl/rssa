/*
 *   R package for Singular Spectrum Analysis
 *   Copyright (c) 2009 Konstantin Usevich <usevich.k.d@gmail.com>
 *   Copyright (c) 2009-2010 Anton Korobeynikov <asl@math.spbu.ru>
 *   Copyright (c) 2014-2015 Alex Shlemov <shlemovalex@gmail.com>
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
#include "masks.h"
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
  Rcomplex *circ_freq;
#endif
  R_len_t rank;
  R_len_t *window;
  R_len_t *factor;
  R_len_t *length;
  area_indices *col_ind;
  area_indices *row_ind;
  unsigned *weights;
} hbhankel_matrix;

static inline R_len_t prod(R_len_t rank, const R_len_t *N) {
  R_len_t r, res = 1;

  for (r = 0; r < rank; ++r)
    res *= N[r];

  return res;
}

static inline R_len_t hprod(R_len_t rank, const R_len_t *N) {
  R_len_t r, res = N[0]/2 + 1;

  for (r = 1; r < rank; ++r)
    res *= N[r];

  return res;
}

static unsigned hbhankel_nrow(const void *matrix) {
  const hbhankel_matrix *h = matrix;
  return h->col_ind != NULL ? h->col_ind->num : prod(h->rank, h->window);
}

static unsigned hbhankel_ncol(const void *matrix) {
  const hbhankel_matrix *h = matrix;
  return h->row_ind != NULL ? h->row_ind->num : prod(h->rank, h->factor);
}

/* TODO Get rid off this ugly slow complicated function */
static inline void fill_subarray(double *big,
                                 double *small,
                                 R_len_t rank,
                                 const R_len_t *Nbig,
                                 const R_len_t *Nsmall,
                                 int direction) {
  R_len_t pNsmall = prod(rank, Nsmall);
  R_len_t *mul, i, j, ii, r;

  mul = Calloc(rank, R_len_t);
  mul[0] = 1;
  for (r = 1; r < rank; ++r)
    mul[r] = mul[r-1] * Nbig[r-1];

  for (i = 0; i < pNsmall; ++i) {
    ii = i;
    j = 0;
    for (r = 0; r < rank; ++r) {
      j += ii % Nsmall[r] * mul[r];
      ii /= Nsmall[r];
    }

    if (direction)
      big[j] = small[i];
    else
      small[i] = big[j];
  }

  Free(mul);
}

#if HAVE_FFTW3_H
static void free_circulant(hbhankel_matrix *h) {
  Free(h->window);
  Free(h->factor);
  Free(h->length);

  fftw_free(h->circ_freq);
  fftw_destroy_plan(h->r2c_plan);
  fftw_destroy_plan(h->c2r_plan);
}

static void initialize_circulant(hbhankel_matrix *h,
                                 const double *F,
                                 R_len_t rank,
                                 const R_len_t *N,
                                 const R_len_t *L,
                                 const int *circular) {
  fftw_complex *ocirc;
  fftw_plan p1, p2;
  double *circ;
  R_len_t *revN, r;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(prod(rank, N) * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc(hprod(rank, N) * sizeof(fftw_complex));

  /* Estimate the best plans for given input length, note, that input data is
     stored in column-major mode, that's why we're passing dimensions in
     *reverse* order */
  revN = Calloc(rank, R_len_t);
  for (r = 0; r < rank; ++r) revN[r] = N[rank - 1 - r];
  p1 = fftw_plan_dft_r2c(rank, revN, circ, ocirc, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_c2r(rank, revN, ocirc, circ, FFTW_ESTIMATE);
  Free(revN);

  /* Fill input buffer */
  memcpy(circ, F, prod(rank, N) * sizeof(double));

  /* Run the plan on input data */
  fftw_execute(p1);

  /* Cleanup and return */
  fftw_free(circ);

  h->circ_freq = ocirc;
  h->r2c_plan = p1;
  h->c2r_plan = p2;

  h->rank = rank;

  h->window = Calloc(rank, R_len_t);
  memcpy(h->window, L, rank * sizeof(R_len_t));

  h->length = Calloc(rank, R_len_t);
  memcpy(h->length, N, rank * sizeof(R_len_t));

  h->factor = Calloc(rank, R_len_t);
  for (r = 0; r < rank; ++r) h->factor[r] = circular[r] ? N[r] : N[r] - L[r] + 1;
}

static void convolveNd_half(const fftw_complex *ox,
                            double *y,
                            R_len_t rank,
                            const R_len_t *N,
                            int conjugate,
                            fftw_plan r2c_plan,
                            fftw_plan c2r_plan) {
  R_len_t i;
  fftw_complex *oy;
  R_len_t pN = prod(rank, N), phN = hprod(rank, N);

  /* Allocate needed memory */
  oy = (fftw_complex*) fftw_malloc(phN * sizeof(fftw_complex));

  /* Compute the Nd-FFT of the matrix y */
  fftw_execute_dft_r2c(r2c_plan, y, oy);

  /* Compute conjugation if needed */
  if (conjugate)
    for (i = 0; i < phN; ++i)
      oy[i] = conj(oy[i]);

  /* Dot-multiply ox and oy, and divide by Nx*...*Nz*/
  for (i = 0; i < phN; ++i)
    oy[i] *= ox[i] / pN;

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(c2r_plan, oy, y);

  /* Cleanup */
  fftw_free(oy);
}

static void convolveNd(double *x,
                       double *y,
                       R_len_t rank,
                       const R_len_t *N,
                       int conjugate,
                       fftw_plan r2c_plan,
                       fftw_plan c2r_plan) {
  fftw_complex *ox;
  R_len_t phN = hprod(rank, N);

  /* Allocate needed memory */
  ox = (fftw_complex*) fftw_malloc(phN * sizeof(fftw_complex));

  /* Compute the NdFFT of the arrays x and y */
  fftw_execute_dft_r2c(r2c_plan, x, ox);

  convolveNd_half(ox, y, rank, N, conjugate, r2c_plan, c2r_plan);

  /* Cleanup */
  fftw_free(ox);
}

static void matmul(double* out,
                   const double* v,
                   const void* matrix,
                   int transposed) {
  const hbhankel_matrix *h = matrix;
  R_len_t rank = h->rank;
  R_len_t *N = h->length;
  R_len_t *L, *K;
  R_len_t pN = prod(rank, N);
  R_len_t i;
  area_indices *col_ind;
  area_indices *row_ind;
  double *circ;

  if (transposed) {
    K = h->window;
    L = h->factor;
    col_ind = h->row_ind;
    row_ind = h->col_ind;
  } else {
    L = h->window;
    K = h->factor;
    row_ind = h->row_ind;
    col_ind = h->col_ind;
  }

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(pN * sizeof(double));

  /* Fill the arrays */
  memset(circ, 0, pN * sizeof(double));
  if (row_ind == NULL) {
    fill_subarray(circ, (double*)v, rank, N, K, 1);
  } else {
    for (i = 0; i < row_ind->num; ++i) {
      circ[row_ind->ind[i]] = v[i];
    }
  }

  convolveNd_half(h->circ_freq, circ,
                  rank, h->length,
                  1,
                  h->r2c_plan, h->c2r_plan);

  /* Cleanup and return */
  if (col_ind == NULL) {
    fill_subarray(circ, out, rank, N, L, 0);
  } else {
    for (i = 0; i < col_ind->num; ++i) {
      out[i] = circ[col_ind->ind[i]];
    }
  }

  fftw_free(circ);
}

static R_INLINE void hbhankelize_fft(double *F,
                                     const double *U, const double *V,
                                     const hbhankel_matrix* h) {
  R_len_t *N = h->length;
  R_len_t *L = h->window;
  R_len_t *K = h->factor;
  R_len_t rank = h->rank;
  R_len_t pN = prod(rank, N);
  R_len_t i;

  double *iU, *iV;

  /* Allocate needed memory */
  iU = (double*) fftw_malloc(pN * sizeof(double));
  iV = (double*) fftw_malloc(pN * sizeof(double));

  /* Fill the arrays */
  memset(iU, 0, pN * sizeof(double));
  if (h->col_ind == NULL) {
    fill_subarray(iU, (double*)U, rank, N, L, 1);
  } else {
    for (i = 0; i < h->col_ind->num; ++i) {
      iU[h->col_ind->ind[i]] = U[i];
    }
  }

  memset(iV, 0, pN * sizeof(double));
  if (h->row_ind == NULL) {
    fill_subarray(iV, (double*)V, rank, N, K, 1);
  } else {
    for (i = 0; i < h->row_ind->num; ++i) {
      iV[h->row_ind->ind[i]] = V[i];
    }
  }

  /* Compute convolution */
  convolveNd(iV, iU, rank, N, 0, h->r2c_plan, h->c2r_plan);

  /* Form the result */
  for (i = 0; i < pN; ++i) {
    if (h->weights[i]) {
      F[i] = iU[i] / h->weights[i];
    } else {
      F[i] = NA_REAL;
    }
  }

  fftw_free(iU);
  fftw_free(iV);
}

SEXP convolveN(SEXP x, SEXP y,
               SEXP input_dim, SEXP output_dim,
               SEXP Conj) {
  SEXP x_dim = getAttrib(x, R_DimSymbol);
  SEXP y_dim = getAttrib(y, R_DimSymbol);
  R_len_t rank = length(input_dim);
  R_len_t *N = INTEGER(input_dim);
  R_len_t pN = prod(rank, N), phN = hprod(rank, N);
  int conjugate = LOGICAL(Conj)[0];

  fftw_complex *ox, *oy;
  fftw_plan r2c_plan, c2r_plan;
  double *circ;
  R_len_t *revN, r, i;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(pN * sizeof(double));
  ox = (fftw_complex*) fftw_malloc(phN * sizeof(fftw_complex));
  oy = (fftw_complex*) fftw_malloc(phN * sizeof(fftw_complex));

  /* Estimate the best plans for given input length, note, that input data is
     stored in column-major mode, that's why we're passing dimensions in
     *reverse* order */
  revN = Calloc(rank, R_len_t);
  for (r = 0; r < rank; ++r) revN[r] = N[rank - 1 - r];
  r2c_plan = fftw_plan_dft_r2c(rank, revN, circ, ox, FFTW_ESTIMATE);
  c2r_plan = fftw_plan_dft_c2r(rank, revN, ox, circ, FFTW_ESTIMATE);
  Free(revN);

  /* Fill input buffer by X values*/
  memset(circ, 0, pN * sizeof(double));
  fill_subarray(circ, REAL(x), rank, N, INTEGER(x_dim), 1);

  /* Run the plan on X-input data */
  fftw_execute_dft_r2c(r2c_plan, circ, ox);

  /* Fill input buffer by Y values*/
  memset(circ, 0, pN * sizeof(double));
  fill_subarray(circ, REAL(y), rank, N, INTEGER(y_dim), 1);

  /* Run the plan on Y-input data */
  fftw_execute_dft_r2c(r2c_plan, circ, oy);

  /* Compute conjugation if needed */
  if (conjugate)
    for (i = 0; i < phN; ++i)
      oy[i] = conj(oy[i]);

  /* Dot-multiply ox and oy, and divide by Nx*...*Nz*/
  for (i = 0; i < phN; ++i)
    oy[i] *= ox[i] / pN;

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(c2r_plan, oy, circ);

  SEXP res;
  PROTECT(res = allocVector(REALSXP, prod(rank, INTEGER(output_dim))));
  fill_subarray(circ, REAL(res), rank, N, INTEGER(output_dim), 0);
  /* setAttrib(output_dim, R_NamesSymbol, R_NilValue); */
  setAttrib(res, R_DimSymbol, output_dim);
  /* setAttrib(res, R_DimNamesSymbol, R_NilValue); */

  /* Cleanup */
  fftw_free(ox);
  fftw_free(oy);
  fftw_free(circ);

  /* Return */
  UNPROTECT(1);
  return res;
}
#else
static void fftn_r2c(const double *z, R_len_t rank, const R_len_t *N,
                     Rcomplex *res) {
  SEXP rA, dim, Res;
  R_len_t n = prod(rank, N);

  rA = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(rA), z, sizeof(double) * n);

  dim = PROTECT(allocVector(INTSXP, rank));
  memcpy(INTEGER(dim), N, sizeof(R_len_t) * rank);
  setAttrib(rA, R_DimSymbol, dim);

  Res = PROTECT(eval(lang2(install("fft"), rA), R_GlobalEnv));

  /* Return result */
  memcpy(res, COMPLEX(Res), sizeof(Rcomplex) * n);

  /* Unprotect all */
  UNPROTECT(3);
}

static void fftn_c2r(const Rcomplex *z, R_len_t rank, const R_len_t *N,
                     double *res) {
  SEXP rTrue, cA, dim, Res;
  R_len_t n = prod(rank, N), i;

  rTrue = PROTECT(allocVector(LGLSXP, 1));
  LOGICAL(rTrue)[0] = 1;

  cA = PROTECT(allocVector(CPLXSXP, n));
  memcpy(COMPLEX(cA), z, sizeof(Rcomplex) * n);

  dim = PROTECT(allocVector(INTSXP, rank));
  memcpy(INTEGER(dim), N, sizeof(R_len_t) * rank);
  setAttrib(cA, R_DimSymbol, dim);

  Res = PROTECT(eval(lang3(install("fft"), cA, rTrue), R_GlobalEnv));

  /* Return result */
  for (i = 0; i < n; ++i)
    res[i] = COMPLEX(Res)[i].r;

  /* Unprotect all */
  UNPROTECT(4);
}

static void free_circulant(hbhankel_matrix *h) {
  Free(h->window);
  Free(h->factor);
  Free(h->length);

  Free(h->circ_freq);
}

static void initialize_circulant(hbhankel_matrix *h,
                                 const double *F,
                                 R_len_t rank,
                                 const R_len_t *N,
                                 const R_len_t *L,
                                 const int *circular) {
  Rcomplex *ocirc;
  R_len_t r;

  /* Allocate needed memory */
  ocirc = Calloc(prod(rank, N), Rcomplex);

  /* Perform FFTn on input data */
  /* We don't need input buffer here */
  fftn_r2c(F, rank, N, ocirc);

  /* Return */
  h->circ_freq = ocirc;

  h->rank = rank;

  h->window = Calloc(rank, R_len_t);
  memcpy(h->window, L, rank * sizeof(R_len_t));

  h->length = Calloc(rank, R_len_t);
  memcpy(h->length, N, rank * sizeof(R_len_t));

  h->factor = Calloc(rank, R_len_t);
  for (r = 0; r < rank; ++r) h->factor[r] = circular[r] ? N[r] : N[r] - L[r] + 1;
}

static void convolveNd_half(const Rcomplex *ox,
                            double *y,
                            R_len_t rank,
                            const R_len_t *N,
                            int conjugate) {
  R_len_t i;
  Rcomplex *oy;
  R_len_t pN = prod(rank, N);

  /* Allocate needed memory */
  oy = Calloc(pN, Rcomplex);

  /* Compute the Nd-FFT of the matrix y */
  fftn_r2c(y, rank, N, oy);

  /* Compute conjugation if needed */
  if (conjugate)
    for (i = 0; i < pN; ++i)
      oy[i].i = -oy[i].i;

  /* Dot-multiply ox and oy, and divide by Nx*...*Nz*/
  for (i = 0; i < pN; ++i) {
    Rcomplex x1 = oy[i], x2 = ox[i];
    oy[i].r = (x1.r * x2.r - x1.i * x2.i) / pN;
    oy[i].i = (x1.r * x2.i + x1.i * x2.r) / pN;
  }

  /* Compute the reverse transform to obtain result */
  fftn_c2r(oy, rank, N, y);

  /* Cleanup */
  Free(oy);
}

static void convolveNd(double *x,
                       double *y,
                       R_len_t rank,
                       const R_len_t *N,
                       int conjugate) {
  Rcomplex *ox;
  R_len_t pN = prod(rank, N);

  /* Allocate needed memory */
  ox = Calloc(pN, Rcomplex);

  /* Compute the NdFFT of the arrays x and y */
  fftn_r2c(x, rank, N, ox);

  convolveNd_half(ox, y, rank, N, conjugate);

  /* Cleanup */
  Free(ox);
}

static void matmul(double* out,
                   const double* v,
                   const void* matrix,
                   int transposed) {
  const hbhankel_matrix *h = matrix;
  R_len_t rank = h->rank;
  R_len_t *N = h->length;
  R_len_t *L, *K;
  R_len_t pN = prod(rank, N);
  R_len_t i;
  area_indices *col_ind;
  area_indices *row_ind;
  double *circ;

  if (transposed) {
    K = h->window;
    L = h->factor;
    col_ind = h->row_ind;
    row_ind = h->col_ind;
  } else {
    L = h->window;
    K = h->factor;
    row_ind = h->row_ind;
    col_ind = h->col_ind;
  }

  /* Allocate needed memory */
  circ = Calloc(pN, double);

  /* Fill the arrays */
  memset(circ, 0, pN * sizeof(double));
  if (row_ind == NULL) {
    fill_subarray(circ, (double*)v, rank, N, K, 1);
  } else {
    for (i = 0; i < row_ind->num; ++i) {
      circ[row_ind->ind[i]] = v[i];
    }
  }

  convolveNd_half(h->circ_freq, circ,
                  rank, h->length,
                  1);

  /* Cleanup and return */
  if (col_ind == NULL) {
    fill_subarray(circ, out, rank, N, L, 0);
  } else {
    for (i = 0; i < col_ind->num; ++i) {
      out[i] = circ[col_ind->ind[i]];
    }
  }

  Free(circ);
}


static R_INLINE void hbhankelize_fft(double *F,
                                     const double *U, const double *V,
                                     const hbhankel_matrix* h) {
  R_len_t *N = h->length;
  R_len_t *L = h->window;
  R_len_t *K = h->factor;
  R_len_t rank = h->rank;
  R_len_t pN = prod(rank, N);
  R_len_t i;

  double *iU, *iV;

  /* Allocate needed memory */
  iU = Calloc(pN, double);
  iV = Calloc(pN, double);

  /* Fill the arrays */
  memset(iU, 0, pN * sizeof(double));
  if (h->col_ind == NULL) {
    fill_subarray(iU, (double*)U, rank, N, L, 1);
  } else {
    for (i = 0; i < h->col_ind->num; ++i) {
      iU[h->col_ind->ind[i]] = U[i];
    }
  }

  memset(iV, 0, pN * sizeof(double));
  if (h->row_ind == NULL) {
    fill_subarray(iV, (double*)V, rank, N, K, 1);
  } else {
    for (i = 0; i < h->row_ind->num; ++i) {
      iV[h->row_ind->ind[i]] = V[i];
    }
  }

  /* Compute convolution */
  convolveNd(iV, iU, rank, N, 0);

  /* Form the result */
  for (i = 0; i < pN; ++i) {
    if (h->weights[i]) {
      F[i] = iU[i] / h->weights[i];
    } else {
      F[i] = NA_REAL;
    }
  }

  Free(iU);
  Free(iV);
}

SEXP convolveN(SEXP x, SEXP y,
               SEXP input_dim, SEXP output_dim,
               SEXP Conj) {
  SEXP x_dim = getAttrib(x, R_DimSymbol);
  SEXP y_dim = getAttrib(y, R_DimSymbol);
  R_len_t rank = length(input_dim);
  R_len_t *N = INTEGER(input_dim);
  R_len_t pN = prod(rank, N);
  int conjugate = LOGICAL(Conj)[0];

  Rcomplex *ox, *oy;
  double *circ;
  R_len_t i;

  /* Allocate needed memory */
  circ = Calloc(pN, double);
  ox = Calloc(pN, Rcomplex);
  oy = Calloc(pN, Rcomplex);

  /* Fill input buffer by X values*/
  memset(circ, 0, pN * sizeof(double));
  fill_subarray(circ, REAL(x), rank, N, INTEGER(x_dim), 1);

  /* Run the plan on X-input data */
  fftn_r2c(circ, rank, N, ox);

  /* Fill input buffer by Y values*/
  memset(circ, 0, pN * sizeof(double));
  fill_subarray(circ, REAL(y), rank, N, INTEGER(y_dim), 1);

  /* Run the plan on Y-input data */
  fftn_r2c(circ, rank, N, oy);

  /* Compute conjugation if needed */
  if (conjugate)
    for (i = 0; i < pN; ++i)
      oy[i].i = -oy[i].i;

  /* Dot-multiply ox and oy, and divide by Nx*...*Nz*/
  for (i = 0; i < pN; ++i) {
    Rcomplex x1 = oy[i], x2 = ox[i];
    oy[i].r = (x1.r * x2.r - x1.i * x2.i) / pN;
    oy[i].i = (x1.r * x2.i + x1.i * x2.r) / pN;
  }

  /* Compute the reverse transform to obtain result */
  fftn_c2r(oy, rank, N, circ);

  SEXP res;
  PROTECT(res = allocVector(REALSXP, prod(rank, INTEGER(output_dim))));
  fill_subarray(circ, REAL(res), rank, N, INTEGER(output_dim), 0);
  setAttrib(res, R_DimSymbol, output_dim);

  /* Cleanup */
  Free(ox);
  Free(oy);
  Free(circ);

  /* Return */
  UNPROTECT(1);
  return res;
}
#endif

static void hbhankel_matmul(double* out,
                            const double* v,
                            const void* matrix) {
  matmul(out, v, matrix, 0);
}

static void hbhankel_tmatmul(double* out,
                             const double* v,
                             const void* matrix) {
  matmul(out, v, matrix, 1);
}

static area_indices *alloc_area2d(SEXP mask, SEXP N) {
  R_len_t r, j, ii;
  if (mask == R_NilValue) {
    return NULL;
  }
  area_indices *area = Calloc(1, area_indices);
  int *maskValues = LOGICAL(mask);
  SEXP DIM = getAttrib(mask, R_DimSymbol);
  R_len_t *dimMask = INTEGER(DIM);
  R_len_t rank = length(DIM);
  R_len_t max_ind = prod(rank, dimMask);

  /* Count the number of nonzero elements and allocate the arrays */
  R_len_t ind;
  area->num = 0;
  for (ind = 0; ind < max_ind; ++ind) {
    area->num += maskValues[ind];
  }

  area->ind = Calloc(area->num, R_len_t);

  R_len_t *mul = Calloc(rank, R_len_t);
  mul[0] = 1;
  for (r = 1; r < rank; ++r)
    mul[r] = mul[r-1] * INTEGER(N)[r-1];

  /* Fill in the arrays of indices (not optimal) */
  R_len_t k;
  for (ind = 0, k = 0; ind < max_ind; ++ind) {
    if (maskValues[ind]) {
      j = 0; ii = ind;
      for (r = 0; r < rank; ++r) {
        j += ii % dimMask[r] * mul[r];
        ii /= dimMask[r];
      }
      area->ind[k] = j;
      ++k;
    }
  }

  Free(mul);

  return area;
}

static void hbhmat_finalizer(SEXP ptr) {
  ext_matrix *e;
  hbhankel_matrix *h;

  if (TYPEOF(ptr) != EXTPTRSXP)
    return;

  e = R_ExternalPtrAddr(ptr);
  if (!e)
    return;

  if (strcmp(e->type, "hbhankel matrix"))
    return;

  h = e->matrix;

  free_area(h->col_ind);
  free_area(h->row_ind);
  Free(h->weights);

  free_circulant(h);
  Free(h);

  Free(e);
  R_ClearExternalPtr(ptr);
}

SEXP initialize_hbhmat(SEXP F, SEXP window,
                       SEXP wmask, SEXP fmask, SEXP weights,
                       SEXP circular) {
  hbhankel_matrix *h;
  ext_matrix *e;
  SEXP hbhmat;

  SEXP N = getAttrib(F, R_DimSymbol);

  /* Allocate memory */
  e = Calloc(1, ext_matrix);
  e->type = "hbhankel matrix";
  e->mulfn = hbhankel_matmul;
  e->tmulfn = hbhankel_tmatmul;
  e->ncol = hbhankel_ncol;
  e->nrow = hbhankel_nrow;

  /* Build toeplitz circulants for hankel matrix */
  h = Calloc(1, hbhankel_matrix);
  initialize_circulant(h, REAL(F), length(N), INTEGER(N), INTEGER(window), LOGICAL(circular));
  /* TODO: add a check for correct window sizes */
  h->col_ind = alloc_area2d(wmask, N);
  h->row_ind = alloc_area2d(fmask, N);
  h->weights = alloc_weights(weights);
  e->matrix = h;

  /* Make an external pointer envelope */
  hbhmat = R_MakeExternalPtr(e, install("external matrix"), R_NilValue);
  R_RegisterCFinalizer(hbhmat, hbhmat_finalizer);

  return hbhmat;
}

SEXP is_hbhmat(SEXP ptr) {
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
      strcmp(e->type, "hbhankel matrix") != 0)
    LOGICAL(ans)[0] = 0;

  UNPROTECT(2);

  return ans;
}

SEXP hbhankelize_one_fft(SEXP U, SEXP V, SEXP hmat) {
  SEXP F = NILSXP, tchk;

  /* Perform a type checking */
  PROTECT(tchk = is_hbhmat(hmat));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e;
    hbhankel_matrix *h;
    double *rU = REAL(U), *rV = REAL(V), *rF;

    /* Grab needed data */
    e = R_ExternalPtrAddr(hmat);
    h = e->matrix;

    /* Allocate buffer for output */
    PROTECT(F = allocVector(REALSXP, prod(h->rank, h->length)));
    rF = REAL(F);

    /* Perform the actual hankelization */
    hbhankelize_fft(rF, rU, rV, h);

    UNPROTECT(1);
  } else
    error("pointer provided is not a hankel block-hankel matrix");

  UNPROTECT(1);

  return F;
}
