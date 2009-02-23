
#ifndef __DSTQRB
#define __DSTQRB

int dstqrb(int n, double * d__, double * e,
           double * z__, double * work, int * info);
/*

    Purpose
    =======
    DSTQRB computes all eigenvalues and the last component of the eigenvectors of a
    symmetric tridiagonal matrix using the implicit QL or QR method.
    This is mainly a modification of the CLAPACK subroutine dsteqr.c

    Arguments
    =========

    N       (input) INT
            The order of the matrix.  N >= 0.

    D       (input/output) DOUBLE PRECISION array, dimension (N)
            On entry, the diagonal elements of the tridiagonal matrix.
            On exit, if INFO = 0, the eigenvalues in ascending order.

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)
            On entry, the (n-1) subdiagonal elements of the tridiagonal
            matrix.
            On exit, E has been destroyed.

    Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
            On entry, if  COMPZ = 'V', then Z contains the orthogonal
            matrix used in the reduction to tridiagonal form.
            On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
            orthonormal eigenvectors of the original symmetric matrix,
            and if COMPZ = 'I', Z contains the orthonormal eigenvectors
            of the symmetric tridiagonal matrix.
            If COMPZ = 'N', then Z is not referenced.

    WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
            If COMPZ = 'N', then WORK is not referenced.

    INFO    (output) integer
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
            > 0:  the algorithm has failed to find all the eigenvalues in
                  a total of 30*N iterations; if INFO = i, then i
                  elements of E have not converged to zero; on exit, D
                  and E contain the elements of a symmetric tridiagonal
                  matrix which is orthogonally similar to the original
                  matrix.
    =====================================================================
*/

#endif

