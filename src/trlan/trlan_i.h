
#ifndef __TRLANI_H
#define __TRLANI_H

/* The maximum number of string allowed for a information titile. */
#define TRLAN_STRING_LEN 132


void trl_clear_counter(trl_info * info, int nrow, int mev, int lde);
/*
// Purpose:
// ========
// Clears the counters inside info and performs a minimal check on the input parameters.
//
// Arguments:
// ==========
// info    (input/output) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//          On exit, the information is cleared.
//
// nrow    (input) integer
//          On entry, specifies the number of rows that is on this proccesor.
//
// mev     (input) integer
//          On entry, specifies the number of Ritz pairs, that can be stored in
//          eval and evec.
//
// lde     (input) integer
//          On entry, specifies the leading dimension of the array evec, i.e.,
//          (lde >= nrow).
//
////
*/
void trl_print_setup(trl_info * info, int lbas, int lmis, int lwrk);
/*
// Purpose:
// ========
// Print the definition of the eigenvalue problme.
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// lbas    (input) integer
//          On entry, specifies the size of workspace required to store Lanczos basis, i.e.,
//          nrow*(maxlan-mev).
//
// lmis    (input) integer
//          On entry, specifies the size of miscellenious workspace requred to solve the
//          eigen problem.
//
// lwrk    (input) integer
//          On entry, specifies the size of workspace provided by the user.
//
*/

void trl_shuffle_eig(int nd, int mnd, double *lambda, double *res,
                      trl_info * info, int *kept, int locked);
//
// Purpose
// =======
// The subroutine trl_shuffle_eig accomplishes two tasks:
//  (1) sort the Ritz values so that those to be saved after
//      restart are ordered in the front of the array,
//  (2) determine the number of Ritz values to be saved.
// On return from this subroutine, the Ritz values are in ascending
// order so that DSTEIN can be used to compute the eigenvectors
//
// Arguments
// =========
// nd         (input) INTEGER
//             On entry, specifies the current size of the Lanczos basis (not locked).
//
// mnd        (input) INTEGER
//             On entry, specifies the maximum size of the Lanczos basis.
//             (lanczos basis size used)-(Ritz values locked).
//
// lambda     (input/output) DOUBLE PRECISION ARRAY of DIMENSION (nd)
//             On entry, contains the computed Ritz values.
//
// res        (input/output) DOUBLE PRECISION ARRAY of DIMENSION (nd)
//             On entry, contains the residual norm of the computed Ritz values.
//
// info       (input/output) POINTER to TRLINFO structure
//             On entry, points to the current TRLINFO structure.
//
// kept       (input/output) POINTER to INTEGER
//             On entry, specifies the number of Ritz values saved.
//
////

double trl_ddot(int, const double *, int, const double *, int);
void trl_daxpy(int n, double da, double *dx, int incx, double *dy,
               int incy);
void trl_dcopy(int n, double *dx, int incx, double *dy, int incy);
void trl_dgemm(char *transa, char *transb, int m, int n, int k,
               double alpha, double *a, int lda, double *b, int ldb,
               double beta, double *c, int ldc);
void trl_dgemv(char *trans, int m, int n, double alpha, double *a, int lda,
               double *x, int incx, double beta, double *y, int incy);
void trl_dscal(int n, double da, double *dx, int incx);

#endif
