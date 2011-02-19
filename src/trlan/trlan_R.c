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

#include "extmat.h"
#include "trlan.h"

typedef struct {
  void *matrix;
  double *tmp;
  R_len_t m, n;
} op_param;

#define UNUSED(x) (void)(x)

void extmat_op(int *pnrow, int *pncol,
               double *xin, int *pldx,
               double *yout, int *pldy,
               void *lparam) {
  op_param *param = lparam;
  const ext_matrix *e = param->matrix;
  double *tmp = param->tmp;

  int ncol = *pncol, ldx  = *pldx, ldy  = *pldy, i;
  UNUSED(pnrow);

  for (i = 0; i < ncol; ++i) {
    e->tmulfn(tmp, xin+i*ldx, e->matrix);
    e->mulfn(yout+i*ldy, tmp, e->matrix);
  }
}

void extmat_op_eigen(int *pnrow, int *pncol,
                     double *xin, int *pldx,
                     double *yout, int *pldy,
                     void *lparam) {
  op_param *param = lparam;
  const ext_matrix *e = param->matrix;

  int ncol = *pncol, ldx  = *pldx, ldy  = *pldy, i;
  UNUSED(pnrow);

  for (i = 0; i < ncol; ++i)
    e->mulfn(yout+i*ldy, xin+i*ldx, e->matrix);
}

void extmat_op2(int *pnrow, int *pncol,
                double *xin, int *pldx,
                double *yout, int *pldy,
                void *lparam) {
  op_param *param = lparam;
  const ext_matrix *e = param->matrix;
  R_len_t m = param->m;

  int ncol = *pncol, ldx  = *pldx, ldy  = *pldy, i;
  UNUSED(pnrow);

  for (i = 0; i < ncol; ++i) {
    e->mulfn(yout+i*ldy, xin+i*ldx+m, e->matrix);
    e->tmulfn(yout+i*ldy+m, xin+i*ldx, e->matrix);
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

void dense_op_eigen(int *pnrow, int *pncol,
                    double *xin, int *pldx,
                    double *yout, int *pldy,
                    void *lparam) {
  op_param *param = lparam;
  double *A   = param->matrix;
  double one = 1.0, zero = 0.0; int i1 = 1;
  int ncol = *pncol, ldx  = *pldx, ldy  = *pldy, i;
  int m = param->m, n = param->n;
  char transn = 'N';

  UNUSED(pnrow);

  for (i = 0; i < ncol; ++i)
    F77_CALL(dgemv)(&transn, &m, &n, &one, A, &m,
                    xin+i*ldx, &i1, &zero, yout+i*ldy, &i1);
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
SEXP trlan_svd(SEXP A, SEXP ne, SEXP opts,
               SEXP ilambda, SEXP iU) {
  R_len_t m = 0, n = 0, kmax, lwrk, computed = 0;
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
             R_ExternalPtrTag(A) == install("external matrix")) {
    /* External matrix case */
    ext_matrix *e = R_ExternalPtrAddr(A);
    m = e->nrow(e->matrix); n = e->ncol(e->matrix);
    param.matrix = e;
#ifdef CYCLIC
    opfn = extmat_op2;
#else
    opfn = extmat_op;
#endif
  } else
    error("unsupported input matrix 'A' type");

  /* Compute needed options */
#ifdef CYCLIC
  param.m = m; param.n = n;
  m = n = m + n;
#else
  param.m = m; param.n = n;
#endif
  param.tmp = (double*)R_alloc(n, sizeof(double));

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

  /* Check, whether we have something to start from */
  if (!isNull(ilambda) && !isNull(iU)) {
    int *dimU;

    /* Check, whether we have vector and matrix */
    if (!isVector(ilambda))
      error("lambda provided should be vector");
    if (!isVector(iU))
      error("U provided should be matrix");

    /* Check dimensions */
    dimU = INTEGER(getAttrib(iU, R_DimSymbol));
    if (dimU[0] != m)
      error("invalid row dimension of provided matrix U, expected %d", m);
    if (dimU[1] != length(ilambda))
      warning("column dimension of matrix U differs from the length of lambda");

    /* Determine the safe upper bound of number of already computed vectors */
    computed = length(ilambda);
    computed = imin2(computed, dimU[1]);
    computed = imin2(computed, kmax);
    computed = imin2(computed, 3*neig/4);

    for (i = 0; i < computed; ++i) {
      double lambda = REAL(ilambda)[i];
      eval[i] = lambda * lambda;
    }
    Memcpy(evec, REAL(iU), computed*m);
  }

  /* The Lanczos recurrence is set to start with [1,1,...,1]^T */
  trl_set_iguess(&info, computed, -1, 0, NULL);

  trlan(opfn, &info, m, kmax, eval, evec, m, lwrk, wrk, &param);

  /* Cleanup */
  Free(wrk);

  if (info.stat == 0) {
    if (info.nec < neig) {
      warning("%d singular triplets did not converge within %d iterations.",
              neig, maxiter);
      neig = info.nec;
    }
  } else
    error("nu-TRLan returned error code %d", info.stat);

  /* Form the result */
  PROTECT(F = allocVector(REALSXP, neig)); rF = REAL(F);
  PROTECT(U = allocMatrix(REALSXP, m, neig)); rU = REAL(U);

  for (i = 0; i < neig; ++i) {
    R_len_t idx = info.nec - i - 1;
#ifdef CYCLIC
    rF[i] = eval[idx];
#else
    rF[i] = sqrt(eval[idx]);
#endif
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

SEXP trlan_eigen(SEXP A, SEXP ne, SEXP opts,
                 SEXP ilambda, SEXP iU) {
  R_len_t m = 0, n = 0, kmax, lwrk, computed = 0;
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
    opfn = dense_op_eigen;
  } else if (TYPEOF(A) == EXTPTRSXP &&
             R_ExternalPtrTag(A) == install("external matrix")) {
    /* External matrix case */
    ext_matrix *e = R_ExternalPtrAddr(A);
    m = e->nrow(e->matrix); n = e->ncol(e->matrix);
    param.matrix = e;
    opfn = extmat_op_eigen;
  } else
    error("unsupported input matrix 'A' type");

  /* Compute needed options */
  param.m = m; param.n = n;
  param.tmp = NULL;

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

  /* Check, whether we have something to start from */
  if (!isNull(ilambda) && !isNull(iU)) {
    int *dimU;

    /* Check, whether we have vector and matrix */
    if (!isVector(ilambda))
      error("lambda provided should be vector");
    if (!isVector(iU))
      error("U provided should be matrix");

    /* Check dimensions */
    dimU = INTEGER(getAttrib(iU, R_DimSymbol));
    if (dimU[0] != m)
      error("invalid row dimension of provided matrix U, expected %d", m);
    if (dimU[1] != length(ilambda))
      warning("column dimension of matrix U differs from the length of lambda");

    /* Determine the safe upper bound of number of already computed vectors */
    computed = length(ilambda);
    computed = imin2(computed, dimU[1]);
    computed = imin2(computed, kmax);
    computed = imin2(computed, 3*neig/4);

    for (i = 0; i < computed; ++i)
      eval[i] = REAL(ilambda)[i];

    Memcpy(evec, REAL(iU), computed*m);
  }

  /* The Lanczos recurrence is set to start with [1,1,...,1]^T */
  trl_set_iguess(&info, computed, -1, 0, NULL);

  trlan(opfn, &info, m, kmax, eval, evec, m, lwrk, wrk, &param);

  trl_print_info(&info);

  /* Cleanup */
  Free(wrk);

  if (info.stat == 0) {
    if (info.nec < neig) {
      warning("%d singular triplets did not converge within %d iterations.",
              neig, maxiter);
      neig = info.nec;
    }
  } else
    error("nu-TRLan returned error code %d", info.stat);

  /* Form the result */
  PROTECT(F = allocVector(REALSXP, neig)); rF = REAL(F);
  PROTECT(U = allocMatrix(REALSXP, m, neig)); rU = REAL(U);

  for (i = 0; i < neig; ++i) {
    R_len_t idx = info.nec - i - 1;
    rF[i] = eval[idx];
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
