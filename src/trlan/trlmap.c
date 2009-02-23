/*
  ZTRLan routine (version 1.0)
  Lawrence Berkeley National Lab.
*/

#include <R_ext/BLAS.h>

double trl_ddot(int n, const double *dx, int incx,
                const double *dy, int incy) {
  return F77_CALL(ddot)(&n, dx, &incx, dy, &incy);
}

void trl_dcopy(int n, double *dx, int incx, double *dy, int incy) {

  F77_CALL(dcopy)(&n, dx, &incx, dy, &incy);
}

void trl_dgemv(char *trans, int m, int n, double alpha, double *a, int lda,
               double *x, int incx, double beta, double *y, int incy) {
  F77_CALL(dgemv)(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

void trl_daxpy(int n, double da, double *dx, int incx, double *dy,
               int incy) {
  F77_CALL(daxpy)(&n, &da, dx, &incx, dy, &incy);
}

void trl_dgemm(char *transa, char *transb, int m, int n, int k,
               double alpha, double *a, int lda, double *b, int ldb,
               double beta, double *c, int ldc) {
  F77_CALL(dgemm)(transa, transb,
                  &m, &n, &k,
                  &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

void trl_dscal(int n, double da, double *dx, int incx) {
  F77_CALL(dscal)(&n, &da, dx, &incx);
}
