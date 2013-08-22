/*
 *   R package for Singular Spectrum Analysis
 *   Copyright (c) 2009 Konstantin Usevich <usevich.k.d@gmail.com>
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
  R_len_t num;
  R_len_t *ind; /* Indices in an Nx x Ny array */
} area2d_indices;

typedef struct {
#if HAVE_FFTW3_H
  fftw_complex * circ_freq;
  fftw_plan r2c_plan;
  fftw_plan c2r_plan;
#else
#endif
  struct {R_len_t x; R_len_t y;} window;
  struct {R_len_t x; R_len_t y;} length;
  area2d_indices *row_ind;
  area2d_indices *col_ind;
  unsigned *weights;
} hbhankel_matrix;

static unsigned hbhankel_nrow(const void *matrix) {
  const hbhankel_matrix *h = matrix;
  return h->row_ind != NULL ? h->row_ind->num : h->window.x * h->window.y;
}

static unsigned hbhankel_ncol(const void *matrix) {
  const hbhankel_matrix *h = matrix;
  return h->col_ind != NULL ? h->col_ind->num :
         (h->length.x - h->window.x + 1) * (h->length.y - h->window.y + 1);
}

#if HAVE_FFTW3_H
static void free_circulant(hbhankel_matrix *h) {
  fftw_free(h->circ_freq);
  fftw_destroy_plan(h->r2c_plan);
  fftw_destroy_plan(h->c2r_plan);
}

static void initialize_circulant(hbhankel_matrix *h,
                                 const double *F,
                                 R_len_t Nx, R_len_t Ny,
                                 R_len_t Lx, R_len_t Ly) {
  R_len_t Kx = Nx - Lx + 1, Ky = Ny - Ly + 1, i, j;
  fftw_complex *ocirc;
  fftw_plan p1, p2;
  double *circ;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(Nx * Ny * sizeof(double));
  ocirc = (fftw_complex*) fftw_malloc(Ny * (Nx/2 + 1) * sizeof(fftw_complex));

  /* Estimate the best plans for given input length, note, that input data is
     stored in column-major mode, that's why we're passing dimensions in
     *reverse* order */
  p1 = fftw_plan_dft_r2c_2d(Ny, Nx, circ, ocirc, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_c2r_2d(Ny, Nx, ocirc, circ, FFTW_ESTIMATE);

  /* Fill input buffer */
  for (j = 0; j < Ny; ++j)
    for (i = 0; i < Nx; ++i)
      /* This is pretty ad-hoc solution and needs to be fixed in the future */
      circ[i + Nx*j] = F[(i + Kx - 1) % Nx + Nx*((j + Ky - 1) % Ny)];

  /* Run the plan on input data */
  fftw_execute(p1);

  /* Cleanup and return */
  fftw_free(circ);

  h->circ_freq = ocirc;
  h->r2c_plan = p1;
  h->c2r_plan = p2;
  h->window.x = Lx; h->window.y = Ly;
  h->length.x = Nx; h->length.y = Ny;
}

static void convolve(double *circ,
                     const hbhankel_matrix *h) {
  R_len_t Nx = h->length.x, Ny = h->length.y, i;
  fftw_complex *ocirc;

  /* Allocate needed memory */
  ocirc = (fftw_complex*) fftw_malloc(Ny*(Nx / 2 + 1) * sizeof(fftw_complex));

  /* Compute the FFT of the reversed vector v */
  fftw_execute_dft_r2c(h->r2c_plan, circ, ocirc);

  /* Dot-multiply with pre-computed FFT of toeplitz circulant */
  for (i = 0; i < Ny * (Nx/2 + 1); ++i)
    ocirc[i] = ocirc[i] * h->circ_freq[i];

  /* Compute the reverse transform to obtain result */
  fftw_execute_dft_c2r(h->c2r_plan, ocirc, circ);

  /* Cleanup */
  fftw_free(ocirc);
}

static void hbhankel_matmul(double* out,
                            const double* v,
                            const void* matrix) {
  const hbhankel_matrix *h = matrix;
  R_len_t Nx = h->length.x, Ny = h->length.y;
  R_len_t Lx = h->window.x, Ly = h->window.y;
  R_len_t Kx = Nx - Lx + 1, Ky = Ny - Ly + 1, i, j;
  double *circ;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(Nx * Ny * sizeof(double));

  /* Fill the arrays */
  memset(circ, 0, Nx * Ny * sizeof(double));
  if (h->col_ind == NULL) {
    for (j = 0; j < Ky; ++j)
      for (i = 0; i < Kx; ++i)
        circ[i + j*Nx] = v[Kx*Ky - i - j*Kx - 1];
  } else {
    for (i = 0; i < h->col_ind->num; ++i) {
      circ[(Kx - 1) + (Ky - 1)*Nx - h->col_ind->ind[i]] = v[i];
    }
  }

  convolve(circ, h);

  /* Cleanup and return */
  if (h->row_ind == NULL) {
    for (j = 0; j < Ly; ++j)
      for (i = 0; i < Lx; ++i)
        out[i + j*Lx] = circ[i + j*Nx] / (Nx * Ny);
  } else {
    for (i = 0; i < h->row_ind->num; ++i) {
      out[i] = circ[h->row_ind->ind[i]] / (Nx * Ny);
    }
  }

  fftw_free(circ);
}

static void hbhankel_tmatmul(double* out,
                             const double* v,
                             const void* matrix) {
  const hbhankel_matrix *h = matrix;
  R_len_t Nx = h->length.x, Ny = h->length.y;
  R_len_t Lx = h->window.x, Ly = h->window.y;
  R_len_t Kx = Nx - Lx + 1, Ky = Ny - Ly + 1, i, j;
  double *circ;

  /* Allocate needed memory */
  circ = (double*) fftw_malloc(Nx * Ny * sizeof(double));

  /* Fill the arrays */
  memset(circ, 0, Nx * Ny * sizeof(double));
  if (h->row_ind == NULL) {
    for (j = 0; j < Ly; ++j)
      for (i = 0; i < Lx; ++i)
        circ[(i + Kx - 1) + (j + Ky - 1)*Nx] = v[Lx*Ly - i - j*Lx - 1];
  } else {
    for (i = 0; i < h->row_ind->num; ++i) {
      circ[Nx*Ny - 1 - h->row_ind->ind[i]] = v[i];
    }
  }

  convolve(circ, h);

  /* Cleanup and return */
  if (h->col_ind == NULL) {
    for (j = 0; j < Ky; ++j)
      for (i = 0; i < Kx; ++i)
        out[i + j * Kx] = circ[(i + Lx - 1) + (j + Ly - 1)*Nx] / (Nx * Ny);
  } else {
    for (i = 0; i < h->col_ind->num; ++i) {
      out[i] =  circ[h->col_ind->ind[i] + (Lx-1) + (Ly-1)*Nx] / (Nx * Ny);
    }
  }

  fftw_free(circ);
}

static R_INLINE void hbhankelize_fft(double *F,
                                     const double *U, const double *V,
                                     const hbhankel_matrix* h) {
  R_len_t Nx = h->length.x, Ny = h->length.y;
  R_len_t Lx = h->window.x, Ly = h->window.y;
  R_len_t Kx = Nx - Lx + 1, Ky = Ny - Ly + 1;
  R_len_t i, j;

  double *iU, *iV;
  fftw_complex *cU, *cV;

  /* Allocate needed memory */
  iU = (double*) fftw_malloc(Nx * Ny * sizeof(double));
  iV = (double*) fftw_malloc(Nx * Ny * sizeof(double));
  cU = (fftw_complex*) fftw_malloc(Ny*(Nx / 2 + 1) * sizeof(fftw_complex));
  cV = (fftw_complex*) fftw_malloc(Ny*(Nx / 2 + 1) * sizeof(fftw_complex));

  /* Fill the arrays */
  memset(iU, 0, Nx * Ny * sizeof(double));
  if (h->row_ind == NULL) {
    for (j = 0; j < Ly; ++j)
      for (i = 0; i < Lx; ++i)
        iU[i + j*Nx] = U[i + j*Lx];
  } else {
    for (i = 0; i < h->row_ind->num; ++i) {
      iU[h->row_ind->ind[i]] = U[i];
    }
  }

  memset(iV, 0, Nx * Ny * sizeof(double));
  if (h->col_ind == NULL) {
    for (j = 0; j < Ky; ++j)
      for (i = 0; i < Kx; ++i)
        iV[i + j*Nx] = V[i + j*Kx];
  } else {
    for (i = 0; i < h->col_ind->num; ++i) {
      iV[h->col_ind->ind[i]] = V[i];
    }
  }

  /* Compute the FFTs */
  fftw_execute_dft_r2c(h->r2c_plan, iU, cU);
  fftw_execute_dft_r2c(h->r2c_plan, iV, cV);

   /* Dot-multiply */
  for (i = 0; i < Ny * (Nx/2 + 1); ++i)
    cU[i] = cU[i] * cV[i];

  /* Compute the inverse FFT */
  fftw_execute_dft_c2r(h->c2r_plan, cU, iU);

  /* Form the result */
  for (i = 0; i < Nx * Ny; ++i) {
    if (h->weights[i]) {
      F[i] = iU[i] / h->weights[i] / Nx / Ny;
    }
  }

  fftw_free(iU);
  fftw_free(iV);
  fftw_free(cU);
  fftw_free(cV);
}
#else
static void free_circulant(hbhankel_matrix *h) {
  error("FFTW-less version of 2D-SSA is not implemented yet!");
}

static void initialize_circulant(hbhankel_matrix *h,
                                 const double *F,
                                 R_len_t Nx, R_len_t Ny,
                                 R_len_t Lx, R_len_t Ly) {
  error("FFTW-less version of 2D-SSA is not implemented yet!");
}

static void hbhankel_matmul(double* out,
                            const double* v,
                            const void* matrix) {
  error("FFTW-less version of 2D-SSA is not implemented yet!");
}

static void hbhankel_tmatmul(double* out,
                             const double* v,
                             const void* matrix) {
  error("FFTW-less version of 2D-SSA is not implemented yet!");
}

static R_INLINE void hbhankelize_fft(double *F,
                                     const double *U, const double *V,
                                     const hbhankel_matrix* h) {
  error("FFTW-less version of 2D-SSA is not implemented yet!");
}
#endif

static area2d_indices *alloc_area2d(SEXP mask, R_len_t Nx) {
  if (mask == R_NilValue) {
    return NULL;
  }
  area2d_indices *area = Calloc(1, area2d_indices);
  int *maskValues = LOGICAL(mask);
  R_len_t *dimMask = INTEGER(getAttrib(mask, R_DimSymbol));
  R_len_t max_ind = dimMask[0] * dimMask[1];

  /* Count the number of nonzero elements and allocate the arrays */
  R_len_t ind;
  area->num = 0;
  for (ind = 0; ind < max_ind; ++ind) {
    area->num += maskValues[ind];
  }

  area->ind = Calloc(area->num, R_len_t);

  /* Fill in the arrays of indices (not optimal) */
  R_len_t k;
  for (ind = 0, k = 0; ind < max_ind; ++ind) {
    if (maskValues[ind]) {
      area->ind[k] = ind % dimMask[0] + (ind / dimMask[0]) * Nx;
      ++k;
    }
  }

  return area;
}

static void free_area2d(area2d_indices *area) {
  if (area == NULL) {
    return;
  }
  Free(area->ind);
  Free(area);
}

static unsigned *alloc_weights(SEXP weights) {
  if (weights == R_NilValue) {
    error("the weights should be precomputed.");
  }
  unsigned *wcopy = Calloc(length(weights), unsigned);
  memcpy(wcopy, INTEGER(weights), sizeof(unsigned) * length(weights));
  return wcopy;
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

  free_area2d(h->row_ind);
  free_area2d(h->col_ind);
  Free(h->weights);

  free_circulant(h);
  Free(h);

  Free(e);
  R_ClearExternalPtr(ptr);
}

SEXP initialize_hbhmat(SEXP F, SEXP windowx, SEXP windowy,
                       SEXP wmask, SEXP fmask, SEXP weights) {
  R_len_t Nx, Ny, Lx, Ly;
  hbhankel_matrix *h;
  ext_matrix *e;
  SEXP hbhmat;

  int *dimF = INTEGER(getAttrib(F, R_DimSymbol));
  Nx = dimF[0]; Ny = dimF[1];
  Lx = INTEGER(windowx)[0]; Ly = INTEGER(windowy)[0];

  /* Allocate memory */
  e = Calloc(1, ext_matrix);
  e->type = "hbhankel matrix";
  e->mulfn = hbhankel_matmul;
  e->tmulfn = hbhankel_tmatmul;
  e->ncol = hbhankel_ncol;
  e->nrow = hbhankel_nrow;

  /* Build toeplitz circulants for hankel matrix */
  h = Calloc(1, hbhankel_matrix);
  initialize_circulant(h, REAL(F), Nx, Ny, Lx, Ly);
  /* TODO: add a check for correct window sizes */
  h->row_ind = alloc_area2d(wmask, Nx);
  h->col_ind = alloc_area2d(fmask, Nx);
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

SEXP hbhankel_rows(SEXP ptr) {
  SEXP tchk;
  SEXP ans = NILSXP;

  /* Perform a type checking */
  PROTECT(tchk = is_hbhmat(ptr));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e = R_ExternalPtrAddr(ptr);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = hbhankel_nrow(e->matrix);
    UNPROTECT(1);
  } else
    error("pointer provided is not a hankel block-hankel matrix");

  UNPROTECT(1);

  return ans;
}

SEXP hbhankel_cols(SEXP ptr) {
  SEXP tchk;
  SEXP ans = NILSXP;

  /* Perform a type checking */
  PROTECT(tchk = is_hbhmat(ptr));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e = R_ExternalPtrAddr(ptr);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = hbhankel_ncol(e->matrix);
    UNPROTECT(1);
  } else
    error("pointer provided is not a hankel block hankel matrix");

  UNPROTECT(1);

  return ans;
}

SEXP hbhmatmul(SEXP hmat, SEXP v, SEXP transposed) {
  SEXP Y = NILSXP, tchk;

  /* Perform a type checking */
  PROTECT(tchk = is_hbhmat(hmat));

  if (LOGICAL(tchk)[0]) {
    R_len_t K, L;
    ext_matrix *e;
    hbhankel_matrix *h;

    /* Grab needed data */
    e = R_ExternalPtrAddr(hmat);
    h = e->matrix;

    L = (LOGICAL(transposed)[0] ? hbhankel_ncol(h) : hbhankel_nrow(h));
    K = (LOGICAL(transposed)[0] ? hbhankel_nrow(h) : hbhankel_ncol(h));

    /* Check agains absurd values of inputs */
    if (K != length(v))
      error("invalid length of input vector 'v'");

    /* Allocate output buffer */
    PROTECT(Y = allocVector(REALSXP, L));

    /* Calculate the product */
    if (LOGICAL(transposed)[0])
      hbhankel_tmatmul(REAL(Y), REAL(v), h);
    else
      hbhankel_matmul(REAL(Y), REAL(v), h);

    UNPROTECT(1);
  } else
    error("pointer provided is not a hankel block-hankel matrix");

  UNPROTECT(1);

  return Y;
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
    PROTECT(F = allocVector(REALSXP, h->length.x * h->length.y));
    rF = REAL(F);

    /* Perform the actual hankelization */
    hbhankelize_fft(rF, rU, rV, h);

    UNPROTECT(1);
  } else
    error("pointer provided is not a hankel block-hankel matrix");

  UNPROTECT(1);

  return F;
}
