#ifndef __TRLANCORE_H
#define __TRLANCORE_H

void add_clock_ticks(trl_info * info, clock_t *time, double *rtime,
                     clock_t clk1);
/*
  Purpose
  =======
  add clock ticks to a integer variable, if there is potential for
  overflow convert it to floating-point numbers

  Arguments
  =========
  info   (input) Pointer to structure trl_info.
  On entry, points the current TRL_INFO. Used only to pass the maximum counter
  value for the clock.

  time   (input/output) Pointer to integer
  On entry, contains the processor time spent in one of following operations:
  clk_op    (in matrix-vector operation)
  clk_orth  (in re-orthogonalization)
  clk_res   (in restarting the Lanczos iterations)
  In the case of overflow, the time is saved in rtime, i.e.,
  total time = time + rtime.

  rtime  (input/output) Pointer to double
  On entry, contains the processor time spent in one of following operations:
  tick_o   (in matrix-vector operation)
  tick_h   (in re-orthogonalization)
  tick_r   (in restarting the Lanczos iterations)

  clk1   (input) clock_t
  On entry, contains the processor time at the start of the time interval
  we are interested in.

*/
void print_alpha_beta(trl_info * info, char *title, int i,
                      double *alpha, double *beta);
/*
  Purpose
  =======
  Print the Ith alpha and beta value to the log file. Function trl_print_real_ is used.

  Arguments
  =========
  info    (input) Pointer to structure trl_info_
  On entry, points the current TRL_INFO. The information is printed out to the
  log file specified in trl_info.

  title   (workspace) String of length (STRING_LEN)
  On entry, provides the space to store the title to print out, i.e., "alpha(jnd) ="
  and "beta(jnd) =".

  i       (input) Integer
  On entry, specifies the index of alpha and beta to print out.

  alpha   (input) Double array of dimension (info->maxlan)
  On entry, contains the alpha values.

  beta    (input) Double array of dimension (info->maxlan)
  On entry, contains the beta values.


*/
void print_all_alpha_beta(trl_info * info, char *title, int jnd,
                          double *alfrot, double *betrot);
/*
  Purpose
  =======
  Print all computed alpha and beta to the log file. Function trl_print_real_ is used.

  Arguments
  =========
  info     (input) Pointer to structure trl_info_
  On entry, points to the current TRL_INFO. The information is printed out to
  the log file specified in info.

  title    (workspace) String of length (12+digit of jnd)
  On entry, provides the space to store the title to print out, i.e., "alfrot(1:jnd)..",
  and "beta(1:jnd).."

  jnd      (input) Integer
  On entry, specifies the number of alpha and beta computed so far.

  alfrot   (input) Double precision array of dimension (info->maxlan)
  On entry, contains the alpha computed so far.

  betrot   (input) Double precision array of dimension (info->maxlan)
  On entry, contains the beta computed so far.


*/
void print_lambda_res(trl_info * info, int jnd, double *lambda,
                      double *res);
/*
  Purpose
  =======
  Print the lambda and its residual norm computed so far. Function trl_print_real_ is used.

  Arguments
  =========
  info      (input) Pointer to sructure trl_info
  On entry, points to the current TRL_INFO. The information is printed out to
  the log file specified in info.

  jnd       (input) Integer
  On entry, specifies the number of lambda computed so far.

  lambda    (input) Double precision array of dimension (info->maxlan)
  On entry, contains the lambda computed so far.

  res       (input) Double precision array of dimension (info->maxlen)
  On entry, contains the residual norm of lambda computed so far.

*/
void
trlanczos(trl_matprod op,
          trl_info * info, int nrow, int mev, double *eval,
          double *evec, int lde, double *base, int ldb, int nbas,
          double *wrk, int lwrk, void *lparam);
/*
  Purpose
  =======
  The actual work routine of restarted Lanczos program for real
  symmetric eigenvalue problems

  user may directly invoke this sunroutine but she/he is responsible
  for input correct arguments as needed

  1) info needs to be initialized
  2) if info%nec>0, the first nec elements of eval and first nec
  columns of evec must contain valid eigenpairs
  3) workspace of currect size
  eval(mev)
  evec(lde, mev) (lde >= nrow, mev >= ned)
  base(ldb, info%maxlan-mev+1) (ldb>=nrow, not used if mev>maxlan)
  wrk(lwrk) minimum amount of memory required by TRLANCZOS is
  maxlan*(maxlan+10)
  4) if log files are to be written, the user need to open files on IO
  unit log_io so that the log gile may be written correctly.
  5) the user must set the timing variable info%clk_tot and
  info%clk_max using system_clock function call in order for this
  subroutine to track timing results correctly

  Algorithm
  =========
  0. initialize input vector
  1. do while (more eigenvalues to compute .and. more MATVEC allowed)
  2.    first step
  o   alpha(k+1) = dot_product(q_{k+1}, Aq_{k+1})
  o   rr = A*q_{k+1}-alpha(k+1)*q_{k+1}-\sum_{i=1}^{k} beta(i)q_i
  o   (re-orthogonalization)
  3.    do j = k+2, m
  o     rr = Aq_j
  o     alpha(j) = dot_product(q_j, rr)
  o     rr = rr - alpha(j)*q_j - beta(j-1)*q_{j-1}
  o     (re-orthogonalization)
  end do j = k+2, m
  4.    restarting
  o   call dstqrb to decompose the tridiagonal matrix
  o   perform convergence test
  o   determine what and how many Ritz pairs to save
  o   compute the Ritz pairs to be saved
  end do while

  The re-orthogonalization procedure is implemented in trl_orth.  it
  produces a normalized vector rr that is guaranteed to be orthogonal
  to the current basis.  An error will be reported if it can not
  achieve its goal.

  Arguments
  =========
  ops     (input) Functin pointer.
  On entry, points to the function that performs the matrix-vector operation.
  The operator that defines the eigenvalue problem is expected to have
  the following interface
  void op(nrow, ncol, xin, ldx, yout, ldy)
  nrow  (input) Integer
  On entry, specifies the number of rows in xin and xout.
  ncol  (input) Integer
  On entry, specifies the number of columns in xin and xout.
  xin   (input) double precision array of dimension (ldx,ncol)
  On entry, specifies the vector/vectors for which the matrix-vector
  is performed.
  ldx   (input) Integer
  On entry, specifies the leading dimension of xin
  yout  (output) Double precision array of diimension (ldy,ncol)
  On exit, specifies the resulting vector/vectors.
  ldy   (input) Integer
  On entry, specifies the leading dimension of yout.

  info    (input) Pointer to structure trl_info_
  On entry, points to the current TRL_INFO.

  nrow    (input) Integer
  On entry, specifies the number of rows in eigenvectors.

  mev     (input) Integer
  On entry, specifies the number of columns allocated to store eigenvectors.

  eval    (output) Double array of dimension (mev)
  On successful exit, contains the eigenvalues.

  evec    (output) Double array of dimension (lde,mev)
  On successful exit, contains the eigenvectors.

  lde     (input) Integer
  On entry, specifies the leading dimension of evec.

  base    (workspace) Double precision array of dimension (ldb,nbas)
  Used to hold the lanczos vectors not fit in evec, i.e., nbas=info->maxlan-mev+1.

  ldb     (input) Integer
  On entry, specifies the leading dimension of base.

  nbas    (input) Integer
  On entry, specifies the number of columns in base.

  wrk     (workspace) Double precision array of dimension (lwrk)
  Workspace for lanczos iterations.

  lwrk    (input) Integer
  On entry, specifies the size of workspace provided.


*/
void trl_orth(int nrow, double *v1, int ld1, int m1, double *v2, int ld2,
              int m2, double *rr, int kept, double *alpha, double *beta,
              double *wrk, int lwrk, trl_info * info);
/*
  Purpose
  =======
  Applies full re-orthogonalization;
  1. if (global re-orthogonalization is needed)
  call trl_cgs
  else
  perform extended local re-reorthogonalization
  endif
  2. perform normalization

  Arguments:
  ==========
  nrow   (input) Integer
  On entry, specifies the number of rows in eigenvectors.

  v1     (input) double precision array (ld1,m1)
  On entry, contains the first part of Lanczos basis computed.

  ld1    (input) Integer
  On entry, specifies the leading dimention of v1.

  m1     (input) Integer
  On entry, specifies the number of Lanczos basis in v1.

  v2     (input) double precision array (ld2,m2)
  On entry, contains the second part of Lanczos basis computed.

  ld2    (input) Integer
  On entry, specifies the leading dimention of v2.

  m2     (input) Integer
  On entry, specifies the number of Lanczos basis in v2.

  rr     (input/output) double precision array (nrow)
  On entry, contains the new Lanczos basis computed.
  On exit, contains the next Lanczos basis computed after the orthogonalization.

  kept   (input) Integer
  On etnry, specifies the number of Ritz vectors kept.

  alpha  (input/output) double precision array (m1+m2)
  On entry, contains the alpha values, on exit, they are updated.

  beta   (input/output) double precision array (m1+m2)
  On entry, contains the beta values, on exit, they are updated if necessary,
  (full orthogonalization).

  wrk    (workspace) complex array (lwrk)

  lwrk   (input) Integer
  Specifies the size of workspace.

  info   (input) Pointer to structure trl_info_
  On entry, points to the current TRL_INFO.


*/
void trl_initial_guess(int nrow, double *evec, int lde, int mev,
                       double *base, int ldb, int nbas, double *alpha,
                       double *beta, int *j1, int *j2, trl_info * info,
                       double *wrk, int lwrk);
/*
  Purpose
  =======
  check to make sure the initial guess vector contains valid nonzero numbers if not fill with
  random numbers this routine will also read the checkpoint files to restore the previous state
  of the Lancozs iterations

  Arguments
  =========
  nrow   (input) Integer
  On entry, specifies the number of rows in eigenvectors.

  evec   (input/output) Double array of dimension (lde,mev)
  On entry, the (nec+1)st column contains the initial guess.

  lde    (input) Integer
  On entry, specifies the leading dimention of evec.

  mev    (input) Integer
  On entry, specifies the number of Ritz vectors that can be stored in evec.

  base   (input/output) Double array of dimension (ldb,nbas)
  Stores the Ritz vectors, that cannot be stored in evec.

  ldb    (input) Integer
  On entry, specifies the leading dimention of base.

  nbas   (input) Integer
  On entry, specifies the number of Ritz vectors that can be stored in base

  alpha  (input/output) Double array of dimension (mev+nbas-1)
  On exit, stores alpha values if checkpoint files are provided.

  beta   (input/output) Double array of dimension (mev+nbas-1)
  On exit, stores beta values if checkpoint files are provided.

  j1     (output) Pointer to integer
  On exit, stores j1 (number of Ritz vectors in evec) if checkpoint files are
  provided.

  j2     (output) Pointer to integer
  On exit, stores j1 (number of Ritz vectors in base) if checkpoint files are
  provided.

  info   (input/output) Pointer to trl_info structure
  On entry, points to the data structure to store the current information about
  the eigenvalue problem and the progress of TRLAN.

  wrk    (workspace) Double array of dimension (lwrk)

  lwrk   (input) Integer
  Specifies the size of the workspace.


*/
void trl_tridiag(int nd, double *alpha, double *beta, double *rot,
                 double *alfrot, double *betrot, double *wrk, int lwrk,
                 int *ierr);
/*
  Purpose
  =======
  transforms an real symmetric arrow matrix into a
  symmetric tridiagonal matrix

  Arguments
  =========
  nd       (input) integer
  On entry, specifies the dimention of the arrow matrix.

  alpha    (input) double precision array (nd)
  On entry, contains the alpha values.

  beta     (input) double precision array (nd)
  On entry, contains the beta values

  rot      (workspace) double precision array (nd*nd)
  Used to store the arrow matrix.

  alfrot   (output) double precision array (nd)
  On exit, contains alpha values after rotation.

  betrot   (output) double precision array (nd)
  On exit, contains beta values after rotation.

  wrk      (workspace) double precision array (lwrk)

  lwrk     (input) integer
  Specifies the size of workspace.

  ierr     (output) integer
  Returns the error from LAPACK calls.
*/
void trl_sort_eig(int nd, int lohi, int nec, double ref, double *lambda,
                  double *res);
/*
  Purpose
  =======
  sort the eigenvalues so that the wanted eigenvalues are ouputed to the user in
  front of the arrays. the final Ritz values are in ascending order so that DSTEIN
  can be used to compute the eigenvectors

  Arguments;
  ==========
  nd       (input) integer
  On entry, specifies the size of lambda.

  lohi     (input) integer
  On entry, specifies which eigenvalues are desired.

  nec      (input) integer
  On entry, specifies how many Ritz values have been converged.

  lambda   (input) double precision array (nd)
  On entry, contains the Ritz values.

  res      (input) double precision array (nd)
  On entry, contains the residual norm of the Ritz values.

*/
void trl_get_tvec(int nd, double *alpha, double *beta, int irot, int nrot,
                  double *rot, int nlam, double *lambda, double *yy,
                  int *iwrk, double *wrk, int lwrk, int *ierr);
/*
  Purpose:
  ========
  generating eigenvectors of the projected eigenvalue problem acorrding to the given
  Ritz values using LAPACK routine DSTEIN (inverse iterations).

  Arguments;
  ==========
  nd        (input) integer
  On entry, specifies the size of alpha and beta.

  alpha     (input) doubel precision array (nd)
  On entry, contains the alpha values.

  beta      (input) double precision array (nd)
  On entry, contains the beta values.

  irot      (input) integer
  On entry, specifies the starting column index of yy to apply the rotation.

  nrot      (input) integer
  On entry, specifies the ending column index of yy to apply the rotation.

  rot       (input) double precision array (nrot, nrot)
  On entry, contains the rotation matrix.

  nlam      (input) integer
  On entry, specifies the size of lambda.

  lambda    (input) double precision array (nlam)
  On entry, contains the Ritz values.

  yy        (output) double precision array (nd,nlam)
  On exit, contains the eigenvector of the tri-diagonal matrix.

  iwrk      (workspace) integer array (4nd)
  wrk       (workspace) double precision (lwrk>=5nd)
  lwrk      (input) integer
  specifies the size of workspace.

  ierr      (output) integer
  Error from Lapack subroutines.

*/
void trl_get_tvec_a(int nd, int kept, double *alpha, double *beta,
                    int nlam, double *lambda, double *yy, double *wrk,
                    int lwrk, int *iwrk, int *ierr);
/*
  Purpose
  =======
  compute all eigenvalues and eigenvectors of the projected matrix, using LAPACK routine DSYEV
  The eigenvectors corresponding to lambda(1:nlam) are placed at the first nlam*nd locations
  of yy on exit.

  Arguments;
  ==========
  nd        (input) integer
  On entry, specifies the size of alpha and beta.

  kept      (input) integer
  On entry, specifies the number of Ritz values kept.

  alpha     (input) doubel precision array (nd)
  On entry, contains the alpha values.

  beta      (input) double precision array (nd)
  On entry, contains the beta values.

  nlam      (input) integer
  On entry, specifies the size of lambda.

  lambda    (input) double precision array (nlam)
  On entry, contains the Ritz values.

  yy        (output) double precision array (nd,nlam)
  On exit, contains the eigenvector of the arrow-head matrix.

  iwrk      (workspace) integer array (nd)
  wrk       (workspace) double precision (lwrk)
  lwrk      (input) integer
  specifies the size of workspace.

  ierr      (output) integer
  Error from Lapack subroutines.
*/
void trl_get_eval(int nd, int locked, double *alpha, double *beta,
                  double *lambda, double *res, double *wrk, int lwrk,
                  int *ierr);
/*
  Purpose
  =======
  Evaluates the eigenvalues and their corresponding residual norms of a
  real symmetric tridiagonal matrix

  it returns eigenvalues in two sections
  1) the first section is the locked eigenvalues, their residual norms are zero
  2) the second section contains the new Ritz values in their ascending order.
  res will contain corresponding residual norms

  Arguments;
  ==========
  nd         (input) integer
  On entry specifies, the size of alpha and beta.

  locked     (input) integer
  On entry, specifies the number of Ritz values locked.

  alpha      (input) double preicsion array (nd)
  On entry, contains the alpha values.

  beta       (input) double precision array (nd)
  On entry, contains the beta values.

  lambea     (output) double precision array (nd)
  On exit, contains the Ritz values.

  res        (output) double precision array (nd)
  On exit, contains the residual norm of the Ritz values.

  wrk        (workspace) double precision array (lwrk)

  lwrk       (input) integer
  Specifies the size of workspace.


*/
void trl_set_locking(int jnd, int nlam, double *lambda, double *res,
                     double *yy, int anrm, int *locked);
/*
  Purpose
  =======
  Move the Ritz pairs with extremely small residual norms to the front of the arrays so that
  locking can be performed cleanly.

  Arguments
  =========
  double: lambda(nlam), res(nlam), yy(jnd*nlam)


*/
void trl_ritz_vectors(int nrow, int lck, int ny, double *yy, int ldy,
                      double *vec1, int ld1, int m1, double *vec2,
                      int ld2, int m2, double *wrk, int lwrk);
/*
  Purpose
  =======
  compute the Ritz vectors from the basis vectors and the eigenvectors of the projected system
  the basis vectors may be stored in two separete arrays the result need to be stored back in them
  lwrk should be no less than ny (lwrk>=ny) ***NOT checked inside***

  Arguments
  =========
  nrow   (input) Integer
  On entry, specifies the number of rows in eigenvectors.

  lck    (input) Integer
  On entry, specifies the number of Ritz values locked.

  ny     (input) Integer
  On entry, specifies the number of columns in yy.

  yy     (input) double precision array (ldy,ny)
  On entry, contains the eigenvector of the "tri-diagonal" matrix.

  ldy    (input) Integer
  On entry. specify the leading dimention of yy.

  vec1   (input) double precision array (ld1,m1)
  On entry, contains the first part of Lanczos basis.

  m1     (input) Integer
  On entry, specifies the number of Lanczos basis stored in vec1.

  ld1    (input) Integer
  On entry, specifies the leading dimention of vec1.

  vec2   (input) double precision array (ld2,m2)
  On entry, contains the second part of Lanczos basis.

  m2     (input) Integer
  On entry, specifies the number of Lanczos basis stored in vec2.

  ld2    (input) Integer
  On entry, specifies the leading dimention of vec2.

  wrk    (workspace) double precision array (lwrk)
  yy2    (workspace) double precision array (ldy,ny)
  lwrk   (input)
  Specifies the size of the workspace.


*/
int trl_cgs(trl_info * info, int nrow, double *v1, int ld1, int m1,
            double *v2, int ld2, int m2, double *rr, double *rnrm,
            double *alpha, int *north, double *wrk);
/*
  Purpose
  =======
  Perform full Gram-Schmidt routine -- orthogonalize a new vector against all existing vectors.

  Arguments
  =========
  info   (input) Pointer to structure trl_info_
  On entry, points to the current TRL_INFO.

  nrow   (input) Integer
  On entry, specifies the number of rows in eigenvectors.

  v1     (input) double precision array (ld1,m1)
  On entry, contains the first part of Lanczos basis computed.

  ld1    (input) Integer
  On entry, specifies the leading dimention of v1.

  m1     (input) Integer
  On entry, specifies the number of Lanczos basis in v1.

  v2     (input) double precision array (ld2,m2)
  On entry, contains the second part of Lanczos basis computed.

  ld2    (input) Integer
  On entry, specifies the leading dimention of v2.

  m2     (input) Integer
  On entry, specifies the number of Lanczos basis in v2.

  rr     (input/output) double precision array (nrow)
  On entry, contains the new Lanczos basis computed.
  On exit, contains the next Lanczos basis computed after the orthogonalization.

  rnrm   (output) double precision
  On entry, specifies the norm of the current Lanczos basis.
  On exit, specifies the norm of the new Lanczos basis.

  alpha  (input/output) double precision array (m1+m2)
  On entry, contains the alpha values, on exit, they are updated.

  north  (output)
  On exit, specifies the number of times the full-orthogonalization is applied.

  wrk    (workspace) complex array (m1+m2)


*/
void trl_convergence_test(int nd, double *lambda, double *res,
                          trl_info * info, double *wrk);
/*
  Purpose
  =======
  count the numer of wanted eigenvalues that have small residual norms
  eigenvalues are assumed to be order from small to large

  Arguments;
  ==========
  nd       (input) integer
  On entry, specifies the size of lambda.

  lambda   (input) double precision array (nd)
  On entry, contains the Ritz values.

  res      (input) double precision array (nd)
  On entry, contains the residual norm of the Ritz values.

  info     (input) Pointer to trl_info structure
  On entry, points to the data structure to store the current information about
  the eigenvalue problem and the progress of TRLAN.

  wrk      (workspace) double precision array (2*nd)

*/
int trl_check_dgen(trl_info * info, int jnd, double *lambda,
                   double *res);

#endif
