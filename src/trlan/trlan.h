
#ifndef __TRLAN_H
#define __TRLAN_H

#include <time.h>

/*
  The data structure to store the current information about
  the eigenvalue problem and the progress of TRLAN.
*/
typedef struct strct_trl_info {
    int stat;			/* status  (error code) of TRLAN */

    /* specification of eigenpairs wanted */
    int lohi;			/* which end of spectrum to compute                  */
    /* lohi < 0 --> the smallest eigenvalues             */
    /* lohi = 0 --> whichever converge first             */
    /* lohi > 0 --> the largest eigenvalues              */
    int ned;			/* number of eigenpairs wanted                       */
    int nec;			/* number of eigenpairs converged                    */
    /* if the user has nec correct eigenvectors, then    */
    /* they are expected to be stored at the beginning   */
    /* of the eigenvector array                          */
    double tol;			/* an eigenpair is declared converged if its         */
    /* residual norm is less than tol*||OP||             */

    /* specification of resource allowed to use by TRLAN */
    int mpicom;			/* the MPI communicator                              */
    int maxlan;			/* maximum basis size to be used                     */
    int klan;			/* the actual basis size, this value may be smaller  */
    /* than maxlan.  It is set when restarting.          */
    int maxmv;			/* maximum number of MATVEC allowed                  */
    /* one MATVEC == one application of the operator on  */
    /* one vector                                        */

    int restart;		/* index of restarting schemes                       */
    int locked;			/* number of eigenvalue locked                       */
    int guess;			/* initial guess                                     */
    /* <= 0, user did not provide initial guess, use     */
    /*       [1,1,..,1]                                  */
    /* =  1, user has supplied initial guess, will only  */
    /*       use the first one                           */
    /* >  1, restart with previous check-point file      */

    /* some information about the progress and resouce comsumption    */
    int matvec;			/* number of MATVEC used by TRLAN                    */
    int nloop;			/* number of restart of the Lanczos iterations       */
    int north;			/* number of full orthogonalization invoked          */
    int nrand;			/* number of times random elements are introduced.   */
    /* Random elements are introduced when an invariant  */
    /* subspace is found but not all wanted eigenvalues  */
    /* are computed.                                     */

    /* variables to store timing results */
    clock_t clk_rate;		/* system clock rate (SYSTEM_CLOCK)                 */
    clock_t clk_max;		/* maximum counter value                            */
    clock_t clk_tot;		/* total time spent in TRLAN (in clock ticks)       */
    clock_t clk_op;		/* time in applying the operator (MATVEC)           */
    clock_t clk_orth;		/* time in re-orthogonalization                     */
    clock_t clk_res;		/* time in restarting the Lanczos iterations        */
    double tick_t;		/* the sum of clk_tot and tick_t is the actual time */
    double tick_o;
    double tick_h;
    double tick_r;
    int clk_in;			/* time spent in reading input data file           */
    int wrds_in;		/* number of real(8) words read                    */
    int clk_out;		/* time spent in writing output data file          */
    int wrds_out;		/* number of real(8) words written to file         */

    double anrm;		/* norm of the operator used                       */

    int my_pe;			/* the PE number of current processor              */
    /* (start with 0)                                  */
    int npes;			/* number of PEs in the group                      */
    int nloc;			/* local problem size                              */
    int ntot;			/* global problem size                             */

    /* how much inforation to output during the execution of TRLAN    */
    int verbose;		/* default only print information related to       */
    /* fatal errors                                    */
    /* if verbose > 0, more diagnostic messages        */
    /* are printed to the following files              */
    /* variables needed to measure convergence factor (crat)          */
    /* convergence rate of the restarted Lanczos algorithm is measure */
    /* by the reduction of residual norm per MATVEC.  The residual    */
    /* norm of the target is used.                                    */
    double crat;
    double trgt;		/* the Ritz value that might convege next          */
    double tres;		/*  residual norm of the target                    */
    int tmv;			/* MATVEC used when target and tres were recorded  */
    double avgm;

    /* Stores some convergence history for restart scheme 9           */
    double cfac;
    double ptres;
    double predicted_crate;
    double mgap;
    double mgamma, gamma0;
    double old_target;
    int target_id;
    int old_locked;
    int k1, k2, k;
    double rfact;

    /* Store "shift" */
    double ref;

    int log_io;			/* Fortran I/O unit number to be used for          */
    /* debug output. Used if verbose > 0.              */
    FILE *log_fp;
    char log_file[128];
    /* base of the file names used by TRLAN to store debug info       */
    /* if verbose > 0, the filenames are computed by appending 'PE#'  */
    /* to this basis.                                                 */

    /* check-pointing parameters
       when cpflag is greater than 0, the basis vectors will be written
       out roughly 'cpflag' times.  For simplicitly, each PE writes its
       own portion of the basis vectors to a file with cpfile followed by
       the processor number.  The file is written as unformatted fortran
       files with the following content:
       nrow, kb(basis size)
       alpha(1:kb)
       beta(1:kb)
       1st basis vector, 2nd basis vector, ..., last basis vector
       the residual vector
    */
    int cpflag, cpio;
    FILE *cpt_fp;
    char cpfile[128], oldcpf[128];
    int dummy;			/* a dummy variable to fill space */
/*
// .. end of trl_info_ ..
*/
} trl_info;

typedef void (*trl_matprod) (int *, int *, double *, int *, double *, int *, void *);

void trl_init_info(trl_info * info, int nrow, int mxlan, int lohi,
                   int ned, double tol, int restart, int maxmv,
                   int mpicom);
/*
//
// Purpose:
// ========
// Initializes a TRL_INFO variable. This function must be called before calling
// any other user level routine in TRLAN package.
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN.
//          On exit, points to the initialized data structure.
//
// nrow    (input) integer
//          On entry, specifies the local dimension of the problem.
//
// mxlan   (input) integer
//          On entry, specifies the maximum number of basis vectors to be used.
//
// lohi    (input)  integer
//          On entry, specifies, which end of the spectrum to compute:
//           lohi < 0, then lower end, the smallest eigenvalues
//           lohi > 0, then high end, the largest eigenvalues
//           lohi = 0, then either lower and or high end, whoever converges first
//          are computed.
//
// ned      (input) integer
//           On entry, specifies the number of wanted eigenvalues and eigenvectors.
//
// tol      (optional) double precision
//           If provided, specifies the tolerance on residual norm. By default,
//           tol is set to be sqrt(epsilon).
//
// trestart (optional) integer
//           If provided, specifies the thick-restarting scheme, 1-4. Default is 1.
//
// mxmv     (optionial) integer
//           If provided, specifies the maximum number of matrix-vector multiplication
//           allowed. By default, mxmv is set to be info%ntot*info%ned.
//
// mpicom   (optional) integer
//           If provided, specifites the MPI communicator. By default, it is a duplicate
//           of MPI_COMM_WORLD. In sequential case, this is set to 0.
*/
//
void trl_set_restart(trl_info * info, double rfact);
/*
//
// Purpose
// =======
// Set the (minimum) basis size for the restart 7 and 8 schemes, i.e., the (minimum) basis
// size if rfact * (number of Ritz vectors kept)
//
// Arguments:
// ==========
// info    (input/output) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// rfact   (input) double precision
//          On entry, specify the (minimum) basis size.
//
*/
////
void trl_set_debug(trl_info * info, int msglvl, char *filename);
/*
//
// Purpose:
// ========
// Set information related to debugging, the initialization routine trl_init_info sets the
// parameters so that  no debug information is printed.
//
// Arguments:
// ==========
// info    (input/output) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// msglvl  (input) integer
//          On entry, specifies the new message/verbose level:
//           msglvl <  0         : nothing is printed
//           msglvl = 1, .., 10  : the higher the level, the more debug information is
//                                 printed.
//
// file    (input) character string
//          On entry, specifies the leading part of the debug file name. Each processor will
//          generate its own debug file with the file name formed by appending the processor
//          number to the string FILE. The actual name is composed by the routine
//          TRL_PE_FILENAME in trlaux
//
////
*/
void trl_set_checkpoint(trl_info * info, int cpflag, char *file);
/*
//
// Purpose:
// ========
// Set up the information related to check-pointing
//
// Arguments:
// ==========
// info    (input/output) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// cpflag  (input) integer
//          On entry, spcifies roughly how many titmes checkpoints are written.
//
// file    (input) character string
//          On entry, specifies the leading part of the checkpoint file name. Each processor
//          will generate its own debug file with the file name formed by appending the
//          processor number to the string FILE. The actual name is composed by the routine
//          TRL_PE_FILENAME in trlaux
//
////
*/
void trl_set_iguess(trl_info * info, int nec, int iguess, int ncps,
                    char *cpf );
/*
//
// Purpose:
// ========
// Set up parameters related to initial guesses of the Lanczos iterations, i.e., the number of
// eigenvector already converged (initially assumed to be zero) and whether the user has
// provided initial guess vector (initially assumed no).  It can also tell TRLan to read check
// point file that is different from the default name.
//
// Arguments:
// ==========
// info    (input/output) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// nec     (input) integer
//          On entry, specifies the number of eigenvalues, that have been converged.
//
// iguess  (input) integer
//          On entry, specifies one of the following options:
//           iguess <= 0, user did not provide initial guess, use [1,1,..,1].
//           iguess =  1, user has supplied initial guess, will only use the first one.
//           iguess >  1, restart with previous check-point file.
//
// ncps    (input)
//          If provided, it specifies the name of checkpoints file.
//
////
*/
void trl_print_info(trl_info * info);
/*
//
// Purpose:
// ========
// Provides an uniform way of printing information stored in TRL_INFO_T.  It needs to be
// called by all PEs. In parallel environment, when writing to standard outputd device, only
// PE0 will actually write its local summary information. Note that this function must be
// called on all PEs at the same time ***
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
////
*/
void trl_terse_info(trl_info * info, FILE * iou);
/*
//
// Purpose:
// ========
// It is a more compact version of trl_print_info, i.e., this is a local routine, indivadual
// PE can call it without regard of whether other PEs do the same and the output may be
// written to a different I/O unit number than the log_io
//
// Arguments:
// ==========
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the current information about
//          the eigenvalue problem and the progress of TRLAN.
//
// iou     (input) pointer to a file
//          On entry, points to the file, here the information is outputed.
//
////
*/
void
trl_ritz_projection_(trl_matprod op,
                     trl_info * info, int lde, int mev, double *evec,
                     double *eres, int lwrk, double *wrk, double *base,
                     void *lparam);
/*
//
// Purpose
// =======
// A separate Rayleigh-Ritz projection routine
// Given a set of approximately orthonormal vectors (V), this routine
// performs the following operations
//  (1) V'*V ==> G
//  (2) R'*R :=  G
//  (3) V'*A*V => H1, inv(R')*H1*inv(R) => H
//  (4) Y*D*Y' := H
//  (5) V*inv(R)*Y => V, diag(D) => lambda,
//      r(i) = ||A*V(:,i)-lambda(i)*V(:,i)||
//
// Arguments:
// ==========
// op       (input) function pointer
//           On entry, points to the matrix-vector multiplication routine.
//
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN
//
// evec    (output) double precision vector of lenvth (nrow*mev)
//          On exit, stores the eigenvectors.
//
// lde     (input) integer
//          On entry, specifies the leading dimension of the array evec, i.e.,
//
// eres    (output) double precision vector of length (2*nev)
//          the array to store new Ritz values and residual norms
//
// base    (optional) double precision vector
//          Workspace to store the result of matrix-vector operation.
//
// wrk     (optional) double precision vector of length (lwrk)
//          Workspace to store projection matrix, etc.
//
////
*/
void
trl_rayleigh_quotients(trl_matprod op,
                       trl_info * info, int ncol, double *evec,
                       double *eres, double *base, void *lparam);
/*
//
// Purpose:
// ========
// Compute Rayleigh quotients, when it is given a set of Ritz vectors and Ritz values,
// normalize the Ritz vectors, and compute their Rayleigh quotients to replace the Ritz values.
//
// Arguments:
// ==========
// op       (input) function pointer
//           On entry, points to the matrix-vector multiplication routine.
//
// info     (input) pointer to the structure trl_info_
//           On entry, points to the data structure to store the current information about
//           the eigenvalue problem and the progress of TRLAN.
//
// evec     (input) double precision array of dimension (nloc,ncol)
//           On entry, stores the portion of eigenvectors on this PE.
//
// base     (workspace)
//           The workspace used to store results of MATVEC
//
// eres     (output) double precision array of dimension (ncol)
//           On exist, store new Ritz values and new residual norms, i.e., if there are NEV
//           Ritz pairs, eres(1:NEV) stores the new Rayleigh quotient and eres(nev+1:nev+nev)
//           stores the new residual norms.
//
// base     (optional) double precision array od dimension (nloc)
//           If provided, double precision workspace.
*/

////
void trlan(trl_matprod op,
           trl_info * info, int nrow, int mev, double *eval,
           double *evec, int lde, int lwrk, double *wrk, void *lparam);
//
// Purpose: Top (user) level routines
// ========
// A thick-restart Lanczos routine for computing eigenvalues and
// eigenvectors of a real symmetric operator/matrix (A).
// -- only accept one input vector, the input vector is expected
//    to be stored in the (nec+1)st column of EVEC.
// -- it extends the Lanczos basis one vector at a time.
// -- orthogonality among the Lanczos vectors are maintained using
//    full re-orthogonalization when necessary.
//
// Requirements:
// =============
// 1) User supplies OP with the specified interface.
// 2) If (info%nec>0), evec(1:nrow, 1:info%nec) must contain valid
//    eigenvectors and eval(1:nec) must be the corresponding eigenvalues.
//    These eigenpairs are assumed to have zero residual norms inside
//    TRLAN.
// 3) lde >= nrow.
// 4) The arrays evec and eval must be large enough to store the
//    solutions, i.e., mev >= info%ned and mev >= info%nec.
// 5) The array wrk may be of arbitrary size.  Internally, the workspace
//    size is
//        nrow*max(0,info%ned-size(evec,2))+maxlan*(maxlan+10)
//
//    If wrk is smaller than this, trlan routine will allocate additional
//    workspace to accommodate.
//
// Arguments:
// ==========
// op      (input) function pointer
//         On entry, points to a function that comptues op(X) == A*X,
//         when given a set of vectors X.
//         The operator that defines the eigenvalue problem is expected to
//         have the following interface
//          void op(nrow, ncol, xin, ldx, yout, ldy)
//            nrow  (input) integer
//                   On entry, specifies the number of rows in xin and yout.
//            ncol  (input) integer
//                   On entry, specifies, the number of columns in Xin and
//                   yout.
//            xin   (input) double precision vector of length (ldx*ncol)
//                   On entry, contatins the input vector to be multiplied.
//                   It consists of ncol column vectors with each column
//                   stored in consecutive order.
//            ldx   (input) integer
//                   On entry, specifies the leading dimension of the array
//                   xin, i.e., the i-th column vector starts with element
//                   (i-1)*ldx+1 and ends with element (i-1)*ldx+nrow in xin.
//            yout  (output) double precision vector of length (ldy*ncol)
//                   On exit, contains the result array, i.e., it stores the
//                   result of matrix-vector multiplications.
//            ldy   (input) integer
//                   On entry, specifies the leading dimension of the array
//                   yout, i.e., the i-th column vector starts with element
//                   (i-1)*ldy+1 in Yout and ends with element (i-1)*ldy+nrow.
//
// info    (input) pointer to the structure trl_info_
//          On entry, points to the data structure to store the information
//          about the eigenvalue problem and the progress of TRLAN
//
// nrow    (input) integer
//          On entry, specifies the number of rows that is on this processor.
//
// mev     (input) integer
//          On entry, specifies the number of Ritz pairs, that can be stored in
//          eval and evec.
//
// eval    (output) double precision vector of length (mev)
//          On exist, stores the eigenvalues.
//
// evec    (output) double precision vector of lenvth (nrow*mev)
//          On exit, stores the eigenvectors.
//
// lde     (input) integer
//          On entry, specifies the leading dimension of the array evec, i.e.,
//
// lwrk    (optional) integer
//          On entry, specifies, the size of WRK.  When both WRK and LWRK are
//          present, then LWRK should correctly indicate the size of WRK. If WRK
//          is present by not LWRK, the size of WRK is assumed to be MEV which is
//          only enough to store the residual norms on exit.  If WRK is not
//          present, LWRK is not used even if it is present.
//          (lde >= nrow).
//
// wrk     (optional) workspace
//          If it is provided and there is enough space, the residual norm of
//          the converged eigenpairs will be stored at wrk(1:info%nec) on exit.
//
////
void trl_check_ritz(trl_matprod op,
                    trl_info * info, int nrow, int ncol, double *rvec,
                    int ldrvec, double *alpha, int *check, double *beta,
                    double *eval, int lwrk, double *wrk, void* lparam);
//
// Purpose:
// ========
// Performs a standard check on the computed Ritz pairs.
//
// Arguments:
// ==========
// op       (input) function pointer
//           On entry, points to the matrix-vector multiplication routine.
//
// info     (input) pointer to the structure trl_info_
//           On entry, points to the data structure to store the information
//           about the eigenvalue problem and the progress of TRLAN
//
// nrow     (input) integer
//           On entry, specifies the problem size.
//
// ncol     (input) integer
//           On entry, specifies the number of Ritz values computed.
//
// rvec     (input) double precision array of dimension (nrow,ncol)
//           On entry, specifies the array storing the Ritz vectors.
//
// alpha    (input) double precision array of dimension (ncol)
//           On entry, contains the Ritz values computed.
//
// beta     (optional) double precision array of dimension (ncol)
//           If provided, contaions the residual norms returned from a Lanczos routine.
//
// eval     (optional) double precision array of dimension (ncol)
//           If provided, contains the actual eigenvalues to compute the error in
//           Ritz values.
//
// lwrk     (optional) integer
//           If provided, specifies the size of workspace provided.
//
// wrk      (optional) double precision array of size(lwrk)
//           If provided, double precidion workspace.

#endif

