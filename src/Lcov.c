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



#if HAVE_FFTW3_H
static void compute_L_covariation_matrix_first_row(const double *F, R_len_t N, R_len_t L,
                                 double *R) {
  double *oF;
  fftw_complex *iF, *iFc;
  fftw_plan p1, p2;
  R_len_t i, K;

  /* Length's check */
  if ((L >= N) || (L < 1))
    error("must be 'N' > 'L' >= 1");
  
  K = N - L + 1;
  
  /* Allocate needed memory */
   
  oF = (double*) fftw_malloc(N * sizeof(double));
  iF = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));
  iFc = (fftw_complex*) fftw_malloc((N/2 + 1) * sizeof(fftw_complex));
  
  /* Estimate the best plans for given input length */
  p1 = fftw_plan_dft_r2c_1d(N, oF, iF, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_c2r_1d(N, iF, oF, FFTW_ESTIMATE);

  memcpy(oF, F, N*sizeof(double));
  
  fftw_execute_dft_r2c(p1, oF, iF);
  
  
  memcpy(oF, F, K*sizeof(double));
  memset(oF + K, 0, (N - K)*sizeof(double));  
   
  fftw_execute_dft_r2c(p1, oF, iFc);
 
  /* Dot-product */
  for (i = 0; i < N/2 + 1; ++i)
    iF[i] = iF[i] * conj(iFc[i]);
    
  fftw_execute_dft_c2r(p2, iF, oF);

  /* Return */
  for (i = 0; i < L; ++i)
    R[i] = oF[i] / N;  
  
  /* Cleanup */
  fftw_free(oF);
  fftw_free(iF);
  fftw_free(iFc);
}
#else
static void compute_L_covariation_matrix_first_row(const double *F, R_len_t N, R_len_t L,
                                 double *R) {
  error("Sorry, but computation L-covariations without FFTW is not implemented yet =(");
}
#endif

/* Useless function. Will be deleted soon*/
// SEXP Lcov(SEXP F, SEXP L) {
  // SEXP R;
  
  // R_len_t intL = INTEGER(L)[0];
  
  // /* Allocate output buffer */
  // PROTECT(R = allocVector(REALSXP, intL));

  // compute_L_covariation_matrix_first_row(REAL(F), length(F), intL, REAL(R));

  // UNPROTECT(1);
  // return R;
// }


SEXP Lcov_matrix(SEXP F, SEXP L) {
  R_len_t intL = INTEGER(L)[0];
  R_len_t i, j;
  R_len_t K = length(F) - intL + 1;
  SEXP ans;
  
  PROTECT(ans = allocMatrix(REALSXP, intL, intL));
  double *rans = REAL(ans);
  double *pF = REAL(F);
  compute_L_covariation_matrix_first_row(REAL(F), length(F), intL, rans);
  
  for(j = 1; j < intL; ++j)
    rans[intL*j] = rans[j];
  
  for(j = 1; j < intL; ++j)
    for(i = j; i < intL; ++i)
      rans[i + intL*j] = rans[j + intL*i] = rans[(i-1) + (j-1)*intL] - pF[i-1]*pF[j-1] + pF[i+K-1]*pF[j+K-1];
      
  UNPROTECT(1);
  return ans;
}