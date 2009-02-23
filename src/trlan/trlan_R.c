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
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include <R_ext/BLAS.h>

#include "hankel.h"
#include "trlan.h"

typedef struct {
  void *matrix;
  double *tmp;
  R_len_t m, n;
} op_param;

#define UNUSED(x) (void)(x)

void hankel_op(int *pnrow, int *pncol,
               double *xin, int *pldx,
               double *yout, int *pldy,
               void *lparam) {
  op_param *param = lparam;
  hankel_matrix *hmat = param->matrix;
  double *tmp = param->tmp;

  int ncol = *pncol, ldx  = *pldx, ldy  = *pldy, i;
  UNUSED(pnrow);

  for (i = 0; i < ncol; ++i) {
    _hmatmul2(tmp, xin+i*ldx, hmat, 1);
    _hmatmul2(yout+i*ldy, tmp, hmat, 0);
  }
}

void dense_op(int *pnrow, int *pncol,
              double *xin, int *pldx,
              double *yout, int *pldy,
              void *lparam) {
  op_param *param = lparam;
  double *A   = param->matrix;
  double *tmp = param->tmp;
  double one = 1.0, zero = 0.0; int i1 = 1;
  int ncol = *pncol, ldx  = *pldx, ldy  = *pldy, i;
  int m = param->m, n = param->n;
  char transt = 'T', transn = 'N';

  UNUSED(pnrow);

  for (i = 0; i < ncol; ++i) {
    F77_CALL(dgemv)(&transt, &m, &n, &one, A, &m,
                    xin+i*ldx, &i1, &zero, tmp, &i1);
    F77_CALL(dgemv)(&transn, &m, &n, &one, A, &m,
                    tmp, &i1, &zero, yout+i*ldy, &i1);
  }
}

/* Get the list element named str, or return NULL */
static SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = GET_NAMES(list);
  int i;

  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

#define getScalarListElement(trg, list, str, coerce, def)         \
  do {                                                            \
    SEXP __tmp = getListElement(list, str);                       \
    trg = (__tmp != R_NilValue ? coerce(__tmp) : (def));          \
  } while(0)

/* Main driver routine for TRLAN */
SEXP trlan_svd(SEXP A, SEXP ne, SEXP opts) {
  R_len_t m = 0, n = 0, kmax, lwrk;
  int neig = *INTEGER(ne), maxiter, i, verbose;
  double *wrk, *eval, *evec, tol, *rF, *rU;
  op_param param;
  trl_info info;
  trl_matprod opfn = NULL;
  SEXP F, U, res;

  /* Check source time and grab dimensions */
  if (isMatrix(A)) {
    /* Ordinary matrix case */
    int *dimA = INTEGER(getAttrib(A, R_DimSymbol));
    m = dimA[0]; n = dimA[1];
    param.matrix = REAL(A);
    opfn = dense_op;
  } else if (TYPEOF(A) == EXTPTRSXP &&
             R_ExternalPtrTag(A) == install("hankel matrix")) {
    /* Hankel matrix case */
    hankel_matrix *h = R_ExternalPtrAddr(A);
    m = _hankel_rows(h); n = _hankel_cols(h);
    param.matrix = h;
    opfn = hankel_op;
  } else
    error("unsupported input matrix 'A' type");

  param.tmp = (double*)R_alloc(n, sizeof(double));
  param.m = m; param.n = n;

  /* Compute needed options */

  /* Fix number of requested eigentriples */
  if (neig > m) neig = m;
  if (neig > n) neig = n;

  /* Maximum number of iterations */
  getScalarListElement(kmax, opts, "kmax", asInteger, 5*neig);
  kmax = imin2(kmax, n+1);
  kmax = imin2(kmax, m+1);

  /* Tolerance */
  getScalarListElement(tol, opts, "tol", asReal, sqrt(DOUBLE_EPS));

  /* Maximum number of matrix-vector products */
  getScalarListElement(maxiter, opts, "maxiter", asInteger, neig*m);

  /* Verboseness */
  getScalarListElement(verbose, opts, "verbose", asInteger, 0);

  lwrk = kmax*(kmax+10);
  wrk  = Calloc(lwrk, double);
  eval = Calloc(kmax, double);
  evec = Calloc(kmax*m, double);

  trl_init_info(&info, m, kmax, +1, neig, tol, 7, maxiter, -1);
  info.verbose = verbose;

  /* The Lanczos recurrence is set to start with [1,1,...,1]^T */
  trl_set_iguess(&info, 0, -1, 0, NULL);

  trlan(opfn, &info, m, kmax, eval, evec, m, lwrk, wrk, &param);

  trl_print_info(&info);

  /* Cleanup */
  Free(wrk);

  if (info.stat == 0) {
    if (info.nec < neig)
      warning("%d singular triplets did not converge within %d iterations.",
              neig, maxiter);
  } else
    error("nu-TRLan returned error code %d", info.stat);

  /* Form the result */
  neig = info.nec;

  PROTECT(F = allocVector(REALSXP, neig)); rF = REAL(F);
  PROTECT(U = allocMatrix(REALSXP, m, neig)); rU = REAL(U);

  for (i = 0; i < neig; ++i) {
    R_len_t idx = neig - i - 1;
    rF[i] = sqrt(eval[idx]);
    Memcpy((rU+m*i), (evec+m*idx), m);
  }

  PROTECT(res = list2(F, U));
  SET_TAG(res, install("d"));
  SET_TAG(CDR(res), install("u"));

  /* Cleanup */
  Free(eval); Free(evec);

  UNPROTECT(3);
  return res;
}

