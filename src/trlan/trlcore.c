/*
  ZTRLan routine (version 1.0)
  Lawrence Berkeley National Lab.
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "dsort2_i.h"
#include "dstqrb_i.h"
#include "trlan.h"
#include "trlan_i.h"
#include "trlaux_i.h"
#include "trlcore_i.h"
#include "trl_comm_i.h"

/* Following are the internal subroutines for printing etc used in function
   lanczos_ */

void add_clock_ticks(trl_info * info, clock_t *time, double *rtime,
                     clock_t clk1)
{
  clock_t clk2, clk3;

  clk2 = clock();
#ifdef __CLK_RESTART
  if (clk2 <= clk1) {
    clk3  = (info->clk_max - clk1);
    clk3 += clk2;
  } else {
    clk3 = clk2 - clk1;
  }
#else
  if (clk2 < clk1) {
    clk3  = (info->clk_max - clk1);
    clk3 += (info->clk_max + clk2);
  } else if (clk1 < 0 && clk2 >= 0) {
    clk3  = -clk1;
    clk3 +=  clk2;
  } else {
    clk3 = clk2 - clk1;
  }
#endif
  if (clk3 + (*time) >= (*time)) {
    *time = clk3 + *time;
  } else {
    *rtime = (*rtime) + ((*time) + clk3) / (double) (info->clk_rate);
    *time = 0;
  }
}

void print_alpha_beta(trl_info * info, char *title, int i,
                      double *alpha, double *beta)
{
  sprintf(title, " alpha(%d) =", i);
  trl_print_real(info, title, 1, &alpha[i - 1], 1);
  sprintf(title, "  beta(%d) =", i);
  trl_print_real(info, title, 1, &beta[i - 1], 1);
}

void print_all_alpha_beta(trl_info * info, char *title, int jnd,
                          double *alfrot, double *betrot)
{
  sprintf(title, "alfrot(1:%d)..", jnd);
  trl_print_real(info, title, jnd, alfrot, 1);
  sprintf(title, "betrot(1:%d)..", jnd);
  trl_print_real(info, title, jnd, betrot, 1);
}

void print_lambda_res(trl_info * info, int jnd, double *lambda,
                      double *res)
{
  trl_print_real(info, "Current eigenvalues..", jnd, lambda, 1);
  trl_print_real(info, "Current residual norms..", jnd, res, 1);
}

void print_restart_state(trl_info * info, char *title, int nrow,
                         int mev, double *alpha, double *beta,
                         double *betrot, double *evec, double *yy,
                         int kept, int locked, int *iwrk, double *wrk2,
                         int i2, int jml)
{
  long c__1 = 1;
  int i, j1, j2;

  iwrk[0] = kept + locked;
  iwrk[1] = locked + i2;
  strcpy(title, "Number of saved and locked Ritz pairs ..");
  trl_print_int(info, title, 2, iwrk, 1);
  if (info->verbose > 2) {
    if (iwrk[1] == 0) {
      strcpy(title, "Ritz values saved (ascending ordered) ..");
    } else {
      strcpy(title, "Ritz values saved (may not be ordered) ..");
    }
    trl_print_real(info, title, kept + locked, alpha, 1);
    strcpy(title, "Residual norms of the saved Ritz pairs ..");
    for (i = 0; i < (kept + locked); i++) {
      betrot[i] = fabs(beta[i]);
    }
    trl_print_real(info, title, kept + locked, betrot, 1);
  }
  if (info->verbose > 7) {
    for (j1 = 0; j1 < imin2(kept, info->verbose); j1++) {
      for (j2 = 0; j2 <= j1; j2++) {
        wrk2[j2] =
          trl_ddot(jml, &yy[j2 * jml], c__1, &yy[j1 * jml],
                   c__1);
      }
      wrk2[j1] = wrk2[j1] - 1;
      sprintf(title, "Orthogonality level of y(%d) ..", j1 + 1);
      trl_print_real(info, title, j1 + 1, wrk2, 1);
    }
  }
  if (info->verbose > 10) {
    for (j1 = 0; j1 < imin2(kept, info->verbose); j1++) {
      sprintf(title, "eigenvector %d of Q'AQ ..", j1);
      trl_print_real(info, title, jml, &yy[(j1 - 1) * jml], 1);
    }
  }
  if (info->verbose > 10) {
    int j1n = imin2(nrow, info->verbose);
    for (j1 = 0; j1 < imin2(kept + locked, mev); j1++) {
      sprintf(title, "Ritz vector %d (1:%d) ..", j1, j1n);
      trl_print_real(info, title, j1n, &evec[j1 * nrow], 1);
    }
  }
}

void print_final_state(trl_info * info, char *title, int nrow, int mev,
                       double *eval, double *beta, double *evec,
                       double *yy, int kept, int jml)
{
  int j1;

  strcpy(title, "Final eigenvalues  (in ascending order)..");
  trl_print_real(info, title, kept, eval, 1);
  if (info->verbose > 4) {
    strcpy(title, "Final residual norms..");
    trl_print_real(info, title, kept, beta, 1);
  }
  if (info->verbose > 8) {
    for (j1 = 0; j1 < imin2(kept, info->verbose); j1++) {
      sprintf(title, "Eigenvector %d of Q''AQ ..", j1);
      trl_print_real(info, title, jml, &yy[j1 * jml], 1);
    }
  }
  if (info->verbose > 10) {
    int j1n = imin2(nrow, info->verbose);
    for (j1 = 0; j1 < imin2(kept, mev); j1++) {
      sprintf(title, "Ritz vector %d (1:%d) ..", j1, j1n);
      trl_print_real(info, title, j1n, &evec[j1 * nrow], 1);
    }
  }
}

void write_checkpoint(trl_info * info, char *title, int nrow,
                      double *alpha, double *beta, double *evec,
                      double *base, int lde, int j1n, int jnd, int ldb,
                      int j2n)
{
  int ii, c1, c2;

  trl_pe_filename(138, title, info->cpfile, info->my_pe, info->npes);
  c1 = clock();
  ii = trl_write_checkpoint(title, nrow, alpha, beta, evec, lde, j1n,
                            base, ldb, j2n);
  c2 = clock();
  if (c2 > c1) {
    info->clk_out = info->clk_out + (c2 - c1);
  } else {
    info->clk_out = info->clk_out + ((info->clk_max - c1) + c2);
  }
  info->wrds_out = info->wrds_out + jnd * (nrow + nrow + 2) + nrow + 2;
  info->stat = trl_sync_flag(info->mpicom, ii);
}

void log_error_state(trl_info * info, int kept, int j1, int j2, int jnd,
                     int nrow, int mev, double *eval, double *alpha,
                     double *alfrot, double *beta, double *betrot,
                     double *evec, double *base, double *qa, double *qb,
                     double *rr, char *title, int *iwrk)
{
  FILE *fp = info->log_fp;

  trl_time_stamp(fp);
  strcpy(title, "Dumping the content of the variables on error..");
  iwrk[0] = info->stat;
  trl_print_int(info, title, 1, iwrk, 1);
  trl_terse_info(info, fp);
  fprintf(fp, "This Lanczos iteration started with %d vectors.\n", kept);
  fprintf(fp, "There are %d (%d, %d) Lanczos vectors currently.\n", jnd,
          j1, j2);
  if (jnd != j1 + j2)
    jnd = j1 + j2;
  if (jnd < 0 || jnd > info->klan)
    jnd = 0;
  strcpy(title, "Content of eval ..");
  trl_print_real(info, title, mev, eval, 1);
  if (jnd > 0) {
    sprintf(title, "Alpha(1:%d) .. ", jnd);
    trl_print_real(info, title, jnd, alpha, 1);
    sprintf(title, " Beta(1:%d) .. ", jnd);
    trl_print_real(info, title, jnd, beta, 1);
    sprintf(title, "Alfrot(1:%d) .. ", jnd);
    trl_print_real(info, title, jnd, alfrot, 1);
    sprintf(title, "betrot(1:%d) .. ", jnd);
    trl_print_real(info, title, jnd, betrot, 1);
  }
  if (j1 > 0) {
    strcpy(title, "the First row of evec ..");
    trl_print_real(info, title, j1, evec, nrow);
    sprintf(title, "row %d of evec ..", nrow);
    trl_print_real(info, title, j1, &evec[nrow - 1], nrow);
  }
  if (j2 > 0) {
    strcpy(title, "the First row of base ..");
    trl_print_real(info, title, j2, base, nrow);
    sprintf(title, "row %d of base ..", nrow);
    trl_print_real(info, title, j2, &base[nrow - 1], nrow);
  }
  if (qb != NULL) {
    sprintf(title, "Content of qb (q_%d) ..", jnd - 1);
    trl_print_real(info, title, nrow, qb, 1);
  }
  if (qa != NULL) {
    sprintf(title, "Content of qa (q_%d) ..", jnd);
    trl_print_real(info, title, nrow, qa, 1);
  }
  if (rr != NULL) {
    sprintf(title, "Content of rr (residual == q_%d) ..", jnd + 1);
    trl_print_real(info, title, nrow, rr, 1);
  }
  if (info->my_pe == 0 && info->log_fp != stdout) {
    warning("TRLanczos returned with error\n""Contents of most variables are dumped to log file %s.\n",
            info->log_file);
  }
}

void
trlanczos(trl_matprod op,
          trl_info * info, int nrow, int mev, double *eval,
          double *evec, int lde, double *base, int ldb, int nbas,
          double *wrk, int lwrk, void *lparam) {
  char notrans = 'N';
  int c__1 = 1;
  int i__1 = 1;
  double one = 1.0;

  char title[TRLAN_STRING_LEN];
  int i, i1, i2, j1, j2, jnd, jml, j1n, j2n, kept, prek, ldqa, ldrr,
    count;
  int next_test, lwrk2, chkpnt, locked, degen;
  clock_t clk1;
  int *iwrk;
  double d__1;
  double *alpha, *beta, *rot, *alfrot, *betrot, *lambda, *res, *yy, *wrk2;
  double *rr = NULL, *qa = NULL, *qb = NULL;

  strcpy(title, "");
  clk1 = 0;
  /*
    alpha, beta, alfrot and betrot have fixed locations in wrk, i.e.,
    alpha: wrk(1:maxlan), beta: wrk(maxlan+1:2*maxlan),
    alfrot: wrk(2*maxlan+1:3*maxlan), betrot: wrk(3*maxlan+1:4*maxlan)
  */
  alpha = &wrk[0];
  i1 = info->maxlan + 1;
  i2 = info->maxlan + info->maxlan;
  beta = &wrk[i1 - 1];
  i1 = i2 + 1;
  i2 = i2 + info->maxlan;
  alfrot = &wrk[i1 - 1];
  i1 = i2 + 1;
  i2 += info->maxlan;
  betrot = &wrk[i1 - 1];
  /* allocate an integer workspace. iwrk holds... */
  iwrk = Calloc(4 * info->maxlan, int);
  memset(iwrk, 0, (4 * info->maxlan) * sizeof(int));
  /* chkpnt specifies how often the check point should be written. */
  if (info->cpflag <= 0) {
    /* check point is not written */
    chkpnt = info->maxmv + info->maxlan;
  } else {
    /* check point is written at every chpnt matrix-vector operations. */
    chkpnt = info->maxmv / info->cpflag;
  }
  /* locked specifies the number of eigenvalues converged. */
  locked = info->nec;
  /* assign values to alpha, beta uses the assumption that the content of
     eval(1:nec) are eigenvalues and their residual norms are zero */
  if (locked > 0) {
    memcpy(alpha, eval, locked * sizeof(double));
    memset(beta, 0, locked * sizeof(double));
  }
  /*
    get valid initial guess for the Lanczos iterations
    wrk2 points to the end of available wrk
    (first 4*maxlan hold alpha, beta, alfrot, and betrot)
    to the end of wrk.
    *** j1 and j2 are the number of colums, and not indices ***
    */
  wrk2 = &wrk[i2];
  lwrk2 = lwrk - i2;
  trl_initial_guess(nrow, evec, lde, mev, base, ldb, nbas, alpha,
                    beta, &j1, &j2, info, wrk2, lwrk2);
  /* On return from trl_initial_guess_, j1 is the last column index of evec
     used, and j2 is the last column index of base used, i.e., jnd specifies
     the size of the current lanczos basis. */
  jnd = j1 + j2;
  kept = jnd;
  if (info->stat != 0) {
    if (info->stat < 0 && (info->verbose > 0 || info->my_pe == 0)) {
      qa = NULL;
      qb = NULL;
      rr = NULL;
      kept = 0;
      j1 = 0;
      j2 = 0;
      jnd = 0;
      log_error_state(info, kept, j1, j2, jnd, nrow, mev, eval,
                      alpha, alfrot, beta, betrot, evec, base, qa,
                      qb, rr, title, iwrk);
    }
    Free(iwrk);
    return;
  }
  /* we will perform the first convergence test after next_testmatrix-vector
     multiplications */
  i1 = info->ned - jnd;
  next_test = i1 + imin2(i1, imin2(6, info->ned / 2));
/************************************************
//            -- the TRLan outer loop --           //
// restart if                                      //
//    1. there is more eigenvalues to compute, and //
//    2. more matrix-vector operations are allowed //
// *********************************************** //
*/
  //tick1 = clock();
  count = 0;
  degen = -1;
  while (info->matvec < info->maxmv &&
         (degen != -1 || info->nec < info->ned)) {
    /*
    // jnd is the size of the current lanczos basis, so increment it.
    */
    //tick2 = clock();
    //tick1 = tick2;
    jnd++;
    /*
    // specify the workspace to hold the rotation matrix, that transform the
    // tridiagonal matrix to a diagonal matrix, i.e., the size of the matrix is
    // (jnd-locked) or the size of the current lanczos basis minus the number
    // of the eigenvalues converged so far. The workspace of the rotation matrix
    // is located at the end of wrk.
    */
    i2 = lwrk - (jnd - locked) * (jnd - locked);
    rot = &wrk[i2];
    /*
    // spcify the workspace required for the orthogonalization procedure.
    // the workspace is after the space storing alpha, beta, alfrot, and betrot,
    // but before the space storing the rotation matrx.
    */
    i1 = 4 * info->maxlan + 1;
    wrk2 = &wrk[i1 - 1];
    lwrk2 = i2 - i1 + 1;
    /*
    // check if there are enough workspace for the orthogonalization procedure.
    */
    i1 = imax2(5 * info->maxlan - 3 * locked, 4 * info->maxlan);
    if (lwrk2 < i1) {
      info->stat = -11;
      return;
    }
    /*
    // the first iteration of TRLan
    // qa points to the last lanczos base vector.
    */
    if (j1 < mev) {
      /* there is still enough space in evec, so use it. */
      j1++;
      qa = &evec[(j1 - 1) * lde];
      ldqa = lde;
    } else {
      /* no more space in evec, so use base. */
      j2++;
      qa = &base[(j2 - 1) * ldb];
      ldqa = nrow;
    }
    /*
    // j1n and j2n specify the location of the next lanczos basis, and
    // rr points the next lanczos base vector.
    */
    if (j1 < mev) {
      /* there is still enough space in evec, so use it. */
      j1n = j1 + 1;
      j2n = 0;
      rr = &evec[(j1n - 1) * lde];
      ldrr = lde;
    } else {
      /* no more space in evec, so use base. */
      j1n = mev;
      j2n = j2 + 1;
      rr = &base[(j2n - 1) * ldb];
      ldrr = ldb;
    }
    /*
    // perform matrix-vector multiplication, i.e., rr = A*qa
    // record the total time and the time in MATVEC
    */
    clk1 = clock();
#ifdef __CLK_RESTART
    if (clk1 <= info->clk_tot) {
      info->tick_t +=
        (info->clk_max -
         info->clk_tot) / (double) (info->clk_rate);
      info->tick_t += clk1 / (double) (info->clk_rate);
      info->clk_tot = clk1;
    }
#else
    if (clk1 < info->clk_tot) {
      info->tick_t +=
        (info->clk_max -
         info->clk_tot) / (double) (info->clk_rate);
      info->tick_t +=
        (info->clk_max + clk1) / (double) (info->clk_rate);
      info->clk_tot = clk1;
    } else if (info->clk_tot < 0 && clk1 >= 0) {
      info->tick_t -= info->clk_tot / (double) (info->clk_rate);
      info->tick_t += clk1 / (double) (info->clk_rate);
      info->clk_tot = clk1;
    }
#endif
    op(&nrow, &i__1, qa, &ldqa, rr, &ldrr, lparam);
    add_clock_ticks(info, &(info->clk_op), &(info->tick_o), clk1);
    (info->matvec)++;
    /*
    // computed the next alpha = qa' * A * qa
    */
    alpha[jnd - 1] = trl_ddot(nrow, qa, c__1, rr, c__1);
    trl_g_sum(info->mpicom, 1, &alpha[jnd - 1], wrk2);
    /*
    // Perform the Lanczos orthogonalization.
    // rr = rr - sum_{i=1,...j1}
    //             beta(i)*evec(:,i) - sum_{1,...,j2} beta(j1+i)*base(:,i)
    // Just for a convenience beta[jnd-1]=alpha[jnd-1] just computed.
    */
    beta[jnd - 1] = alpha[jnd - 1];
    /*
    // orthogonalize with lanczos vectors stored in evec, first.
    */
    if (j1 > 2) {
      /*
      // compute rr = rr - [evec(:,1),...,evec(:,i1)]*[beta(1),...,beta(i1)]'
      */
      d__1 = -one;
      trl_dgemv(&notrans, nrow, j1, d__1, evec, lde, beta, c__1, one,
                rr, c__1);
    } else if (j1 == 1) {
      /*
      // there is no beta, so just compute
      //   rr = rr - alpha(1)*qa
      */
      d__1 = -alpha[0];
      trl_daxpy(nrow, d__1, qa, c__1, rr, c__1);
    } else if (j1 == 2) {
      /*
      // there is only one beta, so just do
      //   rr = rr - beta(1)*evec(1:nrow,1) - beta(2)*evec(1:nrow,2)
      */
      d__1 = -beta[0];
      trl_daxpy(nrow, d__1, evec, c__1, rr, c__1);
      d__1 = -beta[1];
      trl_daxpy(nrow, d__1, &evec[lde], c__1, rr, c__1);
    }
    /*
    // orthogonalize with lanczos vectors stored in base, now.
    */
    if (j2 > 2) {
      /*
      // compute rr = rr - [evec(:,1),...,evec(:,i1)]*[beta(j1+1),...,beta(j1+j2)]'
      */
      d__1 = -one;
      trl_dgemv(&notrans, nrow, j2, d__1, base, ldb, &beta[j1], c__1,
                one, rr, c__1);
    } else if (j2 == 1) {
      /*
      // there is no beta, so just compute
      //   rr = rr - beta(jnd)*qa
      */
      d__1 = -beta[jnd - 1];
      trl_daxpy(nrow, d__1, qa, c__1, rr, c__1);
    } else if (j2 == 2) {
      /*
      // there is only one beta, so just do
      // rr = rr - beta(j1+1)*base(1:nrow,1) - beta(jnd)*base(1:nrow,2)
      */
      d__1 = -beta[j1];
      trl_daxpy(nrow, d__1, base, c__1, rr, c__1);
      d__1 = -beta[jnd - 1];
      trl_daxpy(nrow, d__1, &base[ldb], c__1, rr, c__1);
    }
    /*
    // perform re-orthogonalization (full-orthogonalization)
    */
    clk1 = clock();
    trl_orth(nrow, evec, lde, j1, base, ldb, j2, rr, kept, alpha,
             beta, wrk2, lwrk2, info);
    if (info->verbose > 8) {
      /* check orthogonality after the initilization step */
      trl_check_orth(info, nrow, evec, lde, j1n, base, ldb, j2n,
                     wrk2, lwrk2);
    }
    add_clock_ticks(info, &(info->clk_orth), &(info->tick_h), clk1);
    if (info->stat != 0)
      goto end;
    if (info->verbose > 5) {
      print_alpha_beta(info, title, jnd, alpha, beta);
    }
    /*
    // transform the matrix formed by alpha and beta into a
    // tridiagonal matrix, rot stores the transformation matrix
    */
    /* the already-converged part is just diagonal. */
    memcpy(alfrot, alpha, locked * sizeof(double));
    memset(betrot, 0, locked * sizeof(double));
    /*
    // now, diagonalize the rest of matrix.
    */
    i1 = jnd - locked;
    trl_tridiag(i1, &alpha[locked], &beta[locked], rot,
                &alfrot[locked], &betrot[locked], wrk2, lwrk2,
                &(info->stat));
    if (info->stat != 0)
      goto end;
    betrot[jnd - 1] = beta[jnd - 1];
    /*
    // **************************************************** //
    // regular iterations of Lanczos algorithm (inner loop) //
    // loop if                                              //
    //   1. there is space to store lanczos basis           //
    //   2. there is more eigenvalues to compute, and       //
    // **************************************************** //
    */
    while (jnd < info->klan && (degen != -1 || info->nec < info->ned)) {
      /*
      // compute the kth lanczos vector.
      //  qb is (k-2)nd lanczos vector, and qa is (k-1)st lanczos vector.
      */
      qb = qa;
      qa = rr;
      /*
      // increment j1, j2, and jnd.
      */
      j1 = j1n;
      j2 = j2n;
      jnd++;
      /*
      // find the next available space for the kth lanczos vector.
      */
      if (j1n < mev) {
        /* there is still a space in evec. */
        j1n++;
        rr = &evec[(j1n - 1) * lde];
      } else {
        /* no more space in evec, so use base. */
        j2n++;
        if (j2n <= nbas) {
          rr = &base[(j2n - 1) * ldb];
        } else {
          info->stat = -1111;
          goto end;
        }
      }
      /*
      // perform the matrix-vector operation.
      */
      clk1 = clock();
      op(&nrow, &i__1, qa, &ldqa, rr, &ldrr, lparam);
      add_clock_ticks(info, &(info->clk_op), &(info->tick_o), clk1);
      info->matvec = info->matvec + 1;
      //
      /* compute alpha(jnd) = qa' * A * qa */
      alpha[jnd - 1] = trl_ddot(nrow, qa, c__1, rr, c__1);
      trl_g_sum(info->mpicom, 1, &alpha[jnd - 1], wrk2);
      /*
      // the Lanczos orthogonalization (three-term recurrence).
      //   rr = rr - alpha(jnd)*qa - beta(jnd-1)*qb
      */
      d__1 = -alpha[jnd - 1];
      trl_daxpy(nrow, d__1, qa, c__1, rr, c__1);
      d__1 = -beta[jnd - 2];
      trl_daxpy(nrow, d__1, qb, c__1, rr, c__1);
      /*
      // re-orthogonalization, and compute beta(jnd)
      */
      clk1 = clock();
      trl_orth(nrow, evec, lde, j1, base, ldb, j2, rr, kept, alpha,
               beta, wrk2, lwrk2, info);
      add_clock_ticks(info, &(info->clk_orth), &(info->tick_h),
                      clk1);
      /*
      // copy alpha and beta into alfrot and betrot
      */
      alfrot[jnd - 1] = alpha[jnd - 1];
      betrot[jnd - 1] = beta[jnd - 1];
      if (info->stat != 0)
        goto end;
      if (info->verbose > 4) {
        print_alpha_beta(info, title, jnd, alpha, beta);
      }
      /*
      // perform convergence test once in a while
      */
      if (info->matvec >= next_test) {
        if (info->verbose > 5) {
          print_all_alpha_beta(info, title, jnd, alfrot,
                               betrot);
        }
        lambda = wrk2;
        res = &wrk2[jnd];
        /*
        // At return of get_eval lambda are order in the ascending order
        */
        trl_get_eval(jnd, locked, alfrot, betrot, lambda, res,
                     &wrk2[jnd + jnd + 1], lwrk2 - jnd - jnd,
                     &(info->stat));
        if (info->stat != 0)
          goto end;
        if (info->verbose > 2) {
          print_lambda_res(info, jnd, lambda, res);
        }
        i1 = imin2(mev, jnd);
        memcpy(eval, wrk2, i1 * sizeof(double));
        /*
        // At return from convergence_test, lambda are order in the ascending
        // order of the distance from ref if lohi < -1
        */
        trl_convergence_test(jnd, lambda, res, info,
                             &wrk2[jnd + jnd]);
        /*
        // decide when to perform the next test,
        */
        degen = trl_check_dgen(info, jnd, lambda, res);
        if ((degen != -1 || info->nec < info->ned)
            && info->nec > 0) {
          /*
          // assuming a same number of matrix-vector product is required for
          // each eigenvalues to converge.
          */
          next_test =
            (double) (info->ned * info->matvec) /
            (double) (info->nec);
        } else if (info->nec == 0) {
          next_test = next_test + next_test;
          if (info->maxlan == info->ntot) {
            next_test =
              (int) ceil(0.5 *
                         (info->maxlan + info->matvec));
          }
        }
        if (info->verbose > 0)
          trl_print_progress(info);
      }
    }
    /*
    // ************************************************************* //
    // end of inner (regular Lanczos three-term recurrence) loop     //
    // ************************************************************* //
    */
    /*
    // error checking for debugging use
    */
    //tick1 = clock();
    lambda = wrk2;
    res = &wrk2[jnd];
    if (info->verbose > 6) {
      wrk2 = &wrk2[jnd + jnd];
      i2 = lwrk2 - jnd - jnd;
      trl_check_orth(info, nrow, evec, lde, j1n, base, ldb, j2n,
                     wrk2, i2);
      if (info->verbose > 7) {
        trl_check_recurrence(op, info, nrow, mev, evec, lde, j1n,
                             base, ldb, j2n, kept, alpha, beta,
                             wrk2, i2, lparam);
      }
    }
    /*
    // convert the integer counters to floating-point counters
    */
    i2 = info->clk_max / 4;
    if (info->clk_op > i2) {
      info->tick_o =
        info->tick_o + info->clk_op / (double) (info->clk_rate);
      info->clk_op = 0;
    }
    if (info->clk_orth > i2) {
      info->tick_h =
        info->tick_h + info->clk_orth / (double) (info->clk_rate);
      info->clk_orth = 0;
    }
    if (info->clk_res > i2) {
      info->tick_r =
        info->tick_r + info->clk_res / (double) (info->clk_rate);
      info->clk_res = 0;
    }
    /*
    // *** Determine whether to restart ***
    // compute the Ritz values and Ritz vectors if they are not up to
    // date
    */
    clk1 = clock();
    prek = kept;
    jml = jnd - locked;
    i2 = kept - locked + 1;
    if (degen != -1 || info->nec < info->ned) {
      /* need to compute the updated Ritz values and residual norms */
      wrk2 = &wrk[4 * info->maxlan + 2 * jnd];
      lwrk2 = lwrk - i2 * i2 - 4 * info->maxlan - 2 * jnd;
      if (lwrk2 < 3 * jnd) {
        info->stat = -12;
        goto end;
      }
      if (info->verbose > 5) {
        print_all_alpha_beta(info, title, jnd, alfrot, betrot);
      }
      /*
      // Given tridiagonal matrix
      // (diagonals stored in alfrot, and off-diagonals stored in betrot),
      // computes Ritz value (approximate eigenvalues), using dstqrb, and
      // retrned in lambda. the last components of the eigenvectors are also
      // computed.
      // At return, lambda are stored in the ascending order.
      */
      trl_get_eval(jnd, locked, alfrot, betrot, lambda, res, wrk2,
                   lwrk2, &(info->stat));
      if (info->stat != 0)
        goto end;
      if (info->verbose > 2) {
        print_lambda_res(info, jnd, lambda, res);
      }
      /*
      // At return, lambda are stored in the ascending order of the distance
      // from ref if lohi < -1
      // otherwise, they are sorted in the ascending order of lambda.
      */
      trl_convergence_test(jnd, lambda, res, info, wrk2);
      if (info->verbose > 0) {
        trl_print_progress(info);
      }
      degen = trl_check_dgen(info, jnd, lambda, res);
    }
    /*
    // Given the tridiagonal matrix and Ritz values, compute the Ritz vectors
    // (rotational matrix, used for tridiagonalization, is also applied).
    // Also, decide how many vectors to save if restart
    */
    //tick2 = clock();
    //time4 += (tick2-tick1);
    if ((degen != -1 || info->nec < info->ned)
        && info->matvec < info->maxmv) {
      /*
      // prepare to restart, reorder the eigenvales based on the input parameter.
      // At return, lambda kept are ordered in the ascending order.
      */
      trl_shuffle_eig(jml, info->klan - locked, &lambda[locked],
                      &res[locked], info, &kept, locked);
      /*
      // compute eigenvectors using dstein (inverse interations)
      */
      if (kept * 3 < jml) {
        i1 = 4 * info->maxlan + jnd + kept * jml;
        yy = &wrk[4 * info->maxlan + jnd];
        wrk2 = &wrk[i1];
        lwrk2 = lwrk - i1 - i2 * i2;
        trl_get_tvec(jml, &alfrot[locked], &betrot[locked], 0, i2,
                     rot, kept, &lambda[locked], yy, iwrk, wrk2,
                     lwrk2, &(info->stat));
        if (info->stat == 0 && (locked + kept) > 0) {
          memcpy(alpha, lambda,
                 (locked + kept) * sizeof(double));
        }
      }
      /*
      // compute eigenvectors using dsyev (QR)
      */
      if (kept * 3 >= jml || info->stat != 0) {
        if ((locked + kept) > 0)
          memcpy(alfrot, lambda,
                 (locked + kept) * sizeof(double));
        i1 = 4 * info->maxlan + jml * jml;
        yy = &wrk[4 * info->maxlan];
        wrk2 = &wrk[i1];
        lwrk2 = lwrk - i1;

        trl_get_tvec_a(jml, prek - locked, &alpha[locked],
                       &beta[locked], kept, &alfrot[locked], yy,
                       wrk2, lwrk2, iwrk, &(info->stat));
      }
      if (info->stat != 0)
        goto end;
      for (i = 0; i < kept; i++) {
        beta[locked + i] = yy[(i + 1) * jml - 1] * betrot[jnd - 1];
      }
      if (jml > info->ned + (info->ned / 5 + 6)) {
        trl_set_locking(jml, kept, &alpha[locked], &beta[locked],
                        yy, info->anrm, &i2);
      } else {
        i2 = 0;
      }
      /*
      // generate Ritz vectos, reclaim the space pointed by ROT
      */
      i1 = 4 * info->maxlan + kept * jml + jnd;
      wrk2 = &wrk[i1];
      lwrk2 = lwrk - i1;
      trl_ritz_vectors(nrow, locked, kept, yy, jml, evec, lde, j1,
                       base, ldb, j2, wrk2, lwrk2);
      if (info->verbose > 0) {
        print_restart_state(info, title, nrow, mev, alpha, beta,
                            betrot, evec, yy, kept, locked, iwrk,
                            wrk2, i2, jml);
      }
      /*
      // reset the counters and indices to the correct values for restarting
      */
      kept += locked;
      locked += i2;
      info->locked = locked;
      jnd = kept;
      if (jnd <= mev) {
        j1 = jnd;
        j2 = 0;
      } else {
        j1 = mev;
        j2 = jnd - mev;
        if (j2 >= (nbas - 1)) {
          info->stat = -1111;
          goto end;
        }
      }
      if (info->nec > 0) {
        next_test =
          (int) (((double) (info->matvec * info->ned)) /
                 ((double) info->nec));
      } else {
        next_test = next_test + info->maxlan;
      }
      i1 = imin2(mev, jnd);
      if (i1 > 0)
        memcpy(eval, lambda, i1 * sizeof(double));
      /* copying the last Lanczos vector at the end of kept Ritz vectors */
      if (jnd < mev) {
        j1n = j1 + 1;
        j2n = 0;
        memcpy(&evec[(j1n - 1) * lde], rr, nrow * sizeof(double));
      } else {
        j1n = mev;
        j2n = j2 + 1;
        memcpy(&base[(j2n - 1) * ldb], rr, nrow * sizeof(double));
      }
      /* write checkpoint files */
      if (info->matvec >= chkpnt) {
        write_checkpoint(info, title, nrow, alpha, beta, evec,
                         base, lde, j1n, jnd, ldb, j2n);
        chkpnt = chkpnt + info->maxmv / info->cpflag;
      }
    } else {
      /*
        all wanted eigenpairs converged or maximum MATVEC used sort the
        eigenvalues in final output order
      */
      kept = imin2(info->nec, imax2(info->ned, mev - 1));
      info->nec = kept;
      if (kept == 0)
        kept = imin2(mev - 1, info->ned);
      trl_sort_eig(jnd, info->lohi, kept, info->ref, lambda, res);
      memcpy(eval, lambda, kept * sizeof(double));
      if (kept * 3 < jnd) {
        /*
        // eigenvectors of the projection matrix (try inverse
        // interations)
        */
        i1 = kept * jnd + 4 * info->maxlan;
        yy = &wrk[4 * info->maxlan];
        wrk2 = &wrk[i1];
        lwrk2 = lwrk - i1 - i2 * i2;
        trl_get_tvec(jnd, alfrot, betrot, locked, i2, rot,
                     kept, eval, yy, iwrk, wrk2, lwrk2,
                     &(info->stat));
      }
      if (kept * 3 >= jnd || info->stat != 0) {
        /*
        // too many eigenvectors or inverse iterations have failed,
        // try QR
        */
        i1 = 4 * info->maxlan + jnd * jnd;
        yy = &wrk[4 * info->maxlan];
        wrk2 = &wrk[i1];
        lwrk2 = lwrk - i1;
        trl_get_tvec_a(jnd, prek, alpha, beta, kept, eval, yy,
                       wrk2, lwrk2, iwrk, &(info->stat));
        if (info->stat != 0)
          goto end;
      }
      if (kept > 0)
        memcpy(alpha, eval, kept * sizeof(double));
      for (i = 0; i < kept; i++) {
        beta[i] = betrot[jnd - 1] * yy[(1 + i) * jnd - 1];
      }
      /*
      // generate eigenvectos, reclaim the space pointed by ROT
      */
      i1 = kept * jnd + 4 * info->maxlan;
      wrk2 = &wrk[i1];
      lwrk2 = lwrk - i1;

      trl_ritz_vectors(nrow, 0, kept, yy, jnd, evec, lde, j1, base,
                       ldb, j2, wrk2, lwrk2);

      if (info->verbose > 1) {
        print_final_state(info, title, nrow, mev, eval, beta,
                          evec, yy, kept, jml);
      }
      /*
      // reset the counters and indices to be used by check_orth and
      // check_recurrence
      */
      jnd = kept;
      j1 = kept;
      j2 = 0;
      if (j1 < mev) {
        j1n = j1 + 1;
        j2n = 0;
        memcpy(&evec[(j1n - 1) * lde], rr, nrow * sizeof(double));
      } else {
        j1n = mev;
        j2n = 1;
        memcpy(base, rr, nrow * sizeof(double));
      }
      /*
      // write checkpoint files
      */
      if (info->cpflag > 0) {
        write_checkpoint(info, title, nrow, alpha, beta, evec,
                         base, lde, j1n, jnd, ldb, j2n);
      }
    }
    /*
    // check the orthogonality of the basis vectors before restarting
    */
    if (info->verbose > 6) {
      trl_check_orth(info, nrow, evec, lde, j1n, base, ldb, j2n,
                     wrk2, lwrk2);
      if (info->verbose > 7) {
        trl_check_recurrence(op, info, nrow, mev, evec, lde, j1n,
                             base, ldb, j2n, kept, alpha, beta,
                             wrk2, lwrk2, lparam);
      }
    }
    add_clock_ticks(info, &(info->clk_res), &(info->tick_r), clk1);
    info->nloop = info->nloop + 1;
  }
  /*
  // ******************* //
  // end of restart_loop //
  // ******************* //
  */
  /* write the estimated residual norms to the beginning of WRK */
  for (i = 0; i < j1; i++) {
    wrk[i] = fabs(beta[i]);
  }
end:
  if (info->stat < 0 && (info->verbose > 0 || info->my_pe == 0)) {
    log_error_state(info, kept, j1, j2, jnd, nrow, mev, eval, alpha,
                    alfrot, beta, betrot, evec, base, qa, qb, rr,
                    title, iwrk);
  }
  Free(iwrk);
  return;
}

/*
  Purpose
  =======
  Smooth out a vector, i.e.,
  rr = rr + rr + Cshift(rr, 1) + Cshift(rr, -1) in Fortran.
  Used in trl_initial_guess_.

  Arguments
  =========
  n   (input) Integer
  On entry, specifies the dimension of rr.

  rr  (input/output) Double precision array of dimension (n)
  On entry, contains the initial state of the vector. On exit, contain the
  vector after the smothing is applied.
*/

static void smoothrr(int n, double *rr)
{
  int i;

  if (n <= 0)
    return;
  double rr1 = rr[0], rr2;
  rr2 = rr1;
  rr[0] = 2 * rr[0] + rr[2] + rr[n - 1];
  for (i = 1; i < n - 1; i++) {
    double tmp = rr[i];
    rr[i] = 2 * rr[i] + rr[i + 1] + rr2;
    rr2 = tmp;
  }
  rr[n - 1] = 2 * rr[n - 1] + rr[1] + rr2;
}

void trl_initial_guess(int nrow, double *evec, int lde, int mev,
                       double *base, int ldb, int nbas, double *alpha,
                       double *beta, int *j1, int *j2, trl_info * info,
                       double *wrk, int lwrk)
{
  long c__1 = 1;
  int i, j, k, nran, north;
  double tmp, rnrm;
  clock_t ii, jj;
  char file[TRLAN_STRING_LEN];

  // generate random seeds based on current clock ticks
  ii = clock();
  if (info->my_pe > 0) {
    ii = ii - (int) (info->my_pe * sqrt((double) ii));
  }

  j = info->nec;
  if (info->guess > 1) {
    // retrieve a check-point file
    i = info->cpio;
    if (info->oldcpf != 0 && strlen(info->oldcpf) > 0) {
      trl_pe_filename(TRLAN_STRING_LEN, file, info->oldcpf, info->my_pe,
                      info->npes);
    } else {
      trl_pe_filename(TRLAN_STRING_LEN, file, info->cpfile, info->my_pe,
                      info->npes);
    }

    ii = clock();
    i = trl_read_checkpoint(file, nrow, &evec[j * lde], lde,
                            mev - info->nec, j1, base, ldb, nbas, j2,
                            (mev + nbas - 1 - j), &alpha[j],
                            (mev + nbas - 1 - j), &beta[j]);
    info->stat = trl_sync_flag(info->mpicom, i);
    jj = clock();
    if (jj > ii) {
      info->clk_in = jj - ii;
    } else {
      info->clk_in = (info->clk_max - ii) + jj;
    }
    info->wrds_in = (*j1 + *j2) * (nrow + nrow + 2) + nrow + 2;
    *j1 = *j1 + info->nec;
    if (info->stat != 0)
      return;
  } else {
    if (info->guess <= 0) {
      // generate an arbitrary initial starting vector
      // if (info->guess == 0), use the vector [1, 1, ...]^T
      // else perturb some random locations of the above vector
      for (k = 0; k < nrow; k++) {
        evec[j * lde + k] = 1.0;
      }
      nran = imin2(1 - info->guess, lwrk);
      nran = 2 * (nran / 2);
      GetRNGstate();
      if (nran > 0 && nran < nrow) {
        for (k = 0; k < nran; k++) {
          wrk[k] = unif_rand();
        }
        for (i = 0; i < nran - 1; i += 2) {
          ii = (int) (nrow * wrk[i]);
          evec[j * lde + ii] =
            evec[j * lde + ii] + wrk[i + 1] - 0.5;
        }
      } else if (nran >= nrow) {
        for (i = 0; i < nrow; i++) {
          evec[j * lde + i] = unif_rand();
        }
        smoothrr(nrow, &evec[(info->nec) * lde]);
        info->nrand++;
      }
      PutRNGstate();
    }
    *j1 = info->nec;
    *j2 = 0;
  }
  tmp = 0.0;
  // make sure the norm of the next vector can be computed
  wrk[0] = trl_ddot(nrow, &evec[j * lde], c__1, &evec[j * lde], c__1);
  trl_g_sum(info->mpicom, 1, &wrk[0], &wrk[1]);
  if (wrk[0] >= DBL_MIN && wrk[0] <= DBL_MAX) {
    // set rnrm to let trl_CGS normalize evec(1:nrow, j)
    rnrm = sqrt(wrk[0]);
  } else {
    GetRNGstate();
    for (i = 0; i < nrow; i++) {
      evec[j * lde + i] = unif_rand();
    }
    PutRNGstate();
    smoothrr(nrow, &evec[(info->nec) * lde]);
    info->nrand++;
  }
  //
  // orthogonalize initial guess against all existing vectors
  //
  i = 0;
  tmp = 1.0;
  nran = info->nrand;
  north = info->north;
  if (*j1 < mev) {
    info->stat = trl_cgs(info, nrow, evec, lde, *j1, base, ldb, 0,
                         &evec[(*j1) * lde], &rnrm, &tmp, &i, wrk);
  } else if (*j2 <= 0) {
    info->stat = trl_cgs(info, nrow, evec, lde, *j1, evec, lde, 0,
                         base, &rnrm, &tmp, &i, wrk);

  } else {
    info->stat = trl_cgs(info, nrow, evec, lde, *j1, base, ldb, *j2,
                         &base[(*j2) * ldb], &rnrm, &tmp, &i, wrk);
  }

  if (info->verbose > 6) {
    if (*j1 < mev) {
      i = *j1 + 1;
      ii = *j2;
    } else {
      i = *j1;
      ii = *j2 + 1;
    }
    trl_check_orth(info, nrow, evec, lde, *j1, base, ldb, ii, wrk,
                   lwrk);
  }
  return;
}

void trl_orth(int nrow, double *v1, int ld1, int m1, double *v2, int ld2,
              int m2, double *rr, int kept, double *alpha, double *beta,
              double *wrk, int lwrk, trl_info * info)
{
  double zero = 0.0, one = 1.0;
  long c__1 = 1;
  double d__1;
  int i, usecgs, jnd, jm1, no, nr;
  double tmp;
  double *qa, *qb;

  // check for workspace size
  jnd = m1 + m2;
  jm1 = jnd - 1;
  tmp = zero;
  if (ld1 >= nrow && ld2 >= nrow && lwrk >= imax2(4, jnd + jnd)) {
    info->stat = 0;
  } else {
    info->stat = -101;
    return;
  }
//
// compute the norm of the vector RR
//
  wrk[0] = trl_ddot(nrow, rr, c__1, rr, c__1);
  trl_g_sum(info->mpicom, 1, &wrk[0], &wrk[1]);
  if (!(wrk[0] >= zero) || !(wrk[0] <= DBL_MAX)) {
    info->stat = -102;
    return;
  }
  beta[jnd - 1] = sqrt(wrk[0]);
  tmp = alpha[jnd - 1] * alpha[jnd - 1];
  if (jm1 > kept) {
    tmp += (beta[jm1 - 1] * beta[jm1 - 1]);
  } else if (kept > 0) {
    tmp += trl_ddot(jm1, beta, c__1, beta, c__1);
  }

  if (jm1 == kept) {
    usecgs = 1;
  } else if (jnd >= info->ntot) {
    usecgs = 0;
  } else if (DBL_EPSILON * wrk[0] >= tmp) {
    double anorm = 0.0;
    for (i = 0; i < jnd; ++i) {
      d__1 = fabs(alpha[i]) + fabs(beta[i]);
      if (d__1 > anorm)
        anorm = d__1;
    }
    usecgs = (beta[jm1] < DBL_EPSILON * anorm * info->ntot);
  } else {
    usecgs = 1;
  }
//
// whether to perform full re-orthogonalization or extended local
// re-orthogonalization
  if (usecgs != 0) {
    // perform global re-orthogonalization
    nr = info->nrand;
    no = info->north;
    info->stat = trl_cgs(info, nrow, v1, ld1, m1, v2, ld2, m2, rr,
                         &beta[jnd - 1], &alpha[jnd - 1],
                         &(info->north), wrk);
  } else if (jnd > 1) {
    // perform local re-orthogonalization against two previous vectors
    if (m2 > 1) {
      qa = &v2[(m2 - 1) * ld2];
      qb = &v2[(m2 - 2) * ld2];
    } else if (m2 == 1) {
      qa = v2;
      qb = &v1[(m1 - 1) * ld1];
    } else {
      qa = &v1[(m1 - 1) * ld1];
      qb = &v1[(jm1 - 1) * ld1];
    }
    wrk[0] = zero;
    wrk[1] = zero;
    for (i = 0; i < nrow; i++) {
      wrk[0] = wrk[0] + qa[i] * rr[i];
      wrk[1] = wrk[1] + qb[i] * rr[i];
    }
    trl_g_sum(info->mpicom, 2, &wrk[0], &wrk[2]);
    alpha[jnd - 1] = alpha[jnd - 1] + wrk[0];
    d__1 = -wrk[0];
    trl_daxpy(nrow, d__1, qa, c__1, rr, c__1);
    d__1 = -wrk[1];
    trl_daxpy(nrow, d__1, qb, c__1, rr, c__1);
    tmp = one / beta[jnd - 1];
    trl_dscal(nrow, tmp, rr, c__1);
  } else {
    // perform local re-orthogonalization against the only vector
    if (m1 == 1) {
      qa = v1;
    } else {
      qa = v2;
    }
    wrk[0] = trl_ddot(nrow, qa, c__1, rr, c__1);
    trl_g_sum(info->mpicom, 1, &wrk[0], &wrk[1]);
    alpha[jnd - 1] = alpha[jnd - 1] + wrk[0];
    d__1 = -wrk[0];
    trl_daxpy(nrow, d__1, qa, c__1, rr, c__1);
    tmp = one / beta[jnd - 1];
    trl_dscal(nrow, tmp, rr, c__1);
  }
  // when beta(jnd) is exceedingly small, it should be treated as zero
  if (info->stat == 0) {
    if (beta[jnd - 1] <= DBL_EPSILON * fabs(alpha[jnd - 1])) {
      beta[jnd - 1] = zero;
    } else if (jnd >= info->ntot) {
      beta[jnd - 1] = zero;
    }
  }
}

void trl_tridiag(int nd, double *alpha, double *beta, double *rot,
                 double *alfrot, double *betrot, double *wrk, int lwrk,
                 int *ierr) {
  char upper = 'U';
  int i, lwrk2;
  /* special case, nd == 1; */
  if (nd == 0) {
    return;
  } else if (nd <= 1) {
    rot[0] = 1.0;
    alfrot[0] = alpha[0];
    betrot[0] = beta[0];
    *ierr = 0;
    return;
  }
  if (lwrk < nd + nd) {
    *ierr = -11;
    return;
  } else {
    *ierr = 0;
  }

  /* first form the array matrix as a full matrix in rot
     alpha on diagonal, and beta on the last off-diagonal column and row */
  memset(rot, 0, (nd * nd) * sizeof(double));
  for (i = 0; i < nd; i++) {
    rot[i * nd + i] = alpha[i];
  }
  for (i = 0; i < nd - 1; i++) {
    rot[(nd - 1) * nd + i] = beta[i];
    rot[(i + 1) * nd - 1] = beta[i];
  }
  lwrk2 = lwrk - nd;

  /* call LAPACK routines to reduce the matrix into tridiagonal form and
     generate the rotation matrix */
  F77_CALL(dsytrd)(&upper, &nd,
                   rot, &nd, alfrot, betrot,
                   wrk, &wrk[nd], &lwrk2, ierr);
  if (*ierr != 0) {
    *ierr = -112;
    return;
  }
  betrot[nd - 1] = beta[nd - 1];
  F77_CALL(dorgtr)(&upper, &nd, rot, &nd, wrk, &wrk[nd], &lwrk2, ierr);
  if (*ierr != 0) {
    *ierr = -113;
    return;
  }
}

void trl_get_eval(int nd, int locked, double *alpha, double *beta,
                  double *lambda, double *res, double *wrk, int lwrk,
                  int *ierr) {
  int i;
  int d__1, d__2;

  if (lwrk > 3 * nd) {
    *ierr = 0;
  } else {
    *ierr = -121;
    return;
  }
  memcpy(lambda, alpha, nd * sizeof(double));
  memcpy(wrk, &beta[locked], (nd - locked) * sizeof(double));
  d__1 = (long) (nd - locked);
  dstqrb(d__1, &lambda[locked], wrk, &res[locked], &wrk[nd], &d__2);
  *ierr = (int) d__2;
  if (*ierr == 0) {
    memset(res, 0, locked * sizeof(double));
    for (i = locked; i < nd; i++) {
      res[i] = beta[nd - 1] * fabs(res[i]);
    }
  } else {
    *ierr = -122;
  }
}

void trl_convergence_test(int nd, double *lambda, double *res,
                          trl_info * info, double *wrk)
{
  double bnd;
  int i, j, ncl, ncr;
// copy lambda and res to wrk, sort them in ascending order of lambda
  if (info->lohi < -1) {
    dsort2s(nd, info->ref, lambda, res);
  }
  memcpy(&wrk[nd], lambda, nd * sizeof(double));
  for (i = 0; i < nd; i++) {
    //wrk[nd+i] = lambda[i];
    wrk[i] = fabs(res[i]);
  }
  // sort in the ascending order of wrk[nd:2*nd]=lambda
  if (info->lohi == -2) {
    // around ref
    dsort2s(nd, info->ref, &wrk[nd], wrk);
  } else if (info->lohi == -3) {
    // larger than ref
    dsort2su(nd, info->ref, &wrk[nd], wrk);
  } else if (info->lohi == -4) {
    // smaller than ref
    dsort2sd(nd, info->ref, &wrk[nd], wrk);
  } else {
    dsort2(nd, &wrk[nd], wrk);
  }
  //
  // determine the convergence rate of the previous target
  if (info->tmv > 0 && info->matvec > info->tmv) {
    j = 0;
    bnd = fabs(lambda[j] - info->trgt);
    for (i = 0; i < nd; i++) {
      if (fabs(lambda[i] - info->trgt) < bnd) {
        bnd = fabs(lambda[i] - info->trgt);
        j = i;
      }
    }
    if (info->tres > res[j]) {
      bnd = res[j] / info->tres;
      if (bnd > 0.0) {
        info->crat =
          exp(log(bnd) / (double) (info->matvec - info->tmv));
        info->cfac = bnd;
      } else {
        info->crat = 1.0;
        info->cfac = bnd;
      }
    } else {
      info->crat = 1.0;
      info->cfac = bnd;
    }
  }
  // find out who has converged at the lower end of the spectrum
  info->anrm =
    fmax2(info->anrm, fmax2(fabs(wrk[nd + 1]), fabs(wrk[nd + nd - 1])));
  bnd = DBL_MIN + info->tol * info->anrm;
  ncl = 0;
  ncr = nd;
  if (info->lohi <= 0) {
    ncl = nd - 1;
    i = 0;
    while (i < nd) {
      if (wrk[i] < bnd) {
        if (info->lohi == -3 && wrk[i + nd] < info->ref) {
          ncl = i - 1;
          i = nd;
        } else if (info->lohi == -4 && wrk[i + nd] > info->ref) {
          ncl = i - 1;
          i = nd;
        } else {
          i++;
        }
      } else {
        ncl = i - 1;
        i = nd;
      }
    }
  }
  // find out who has converged at the high end of the spectrum
  if (info->lohi >= 0) {
    ncr = 0;
    i = nd - 1;
    while (i >= 0) {
      if (wrk[i] < bnd) {
        i--;
      } else {
        ncr = i + 1;
        i = -1;
      }
    }
  }

  /* determine the number of wanted eigenvalues that have converged
     compute the next target
     ncl = index of wrk corresponding to the smallest eig converged.
     ncr = index of wrk corresponding to the laragest eig converged. */
  info->tmv = info->matvec;
  info->ptres = info->trgt;
  if (info->lohi < 0) {
    info->nec = ncl + 1;
    info->trgt = wrk[nd + imin2(nd - 1, ncl + 1)];
    info->tres = wrk[imin2(nd - 1, ncl + 1)];
  } else if (info->lohi > 0) {
    info->nec = nd - ncr;
    info->trgt = wrk[nd + imax2(0, ncr - 1)];
    info->tres = wrk[imax2(0, ncr - 1)];
  } else {
    if (ncr <= ncl) {
      ncl = nd / 2;
      ncr = ncl + 1;
      info->trgt = wrk[nd + (nd + 1) / 2 - 1];
      info->tres = wrk[(nd + 1) / 2 - 1];
    } else if (wrk[ncl + 1] <= wrk[ncr - 1]) {
      info->trgt = wrk[nd + ncl + 1];
      info->tres = wrk[ncl + 1];
    } else {
      info->trgt = wrk[nd + ncr - 1];
      info->tres = wrk[ncr - 1];
    }
    info->nec = ncl + nd - ncr + 1;
    //for( i=ncl; i<ncr-1; i++ ) {
    for (i = ncl + 1; i < ncr; i++) {
      if (wrk[i] < bnd)
        info->nec = info->nec + 1;
    }
  }
}

void trl_sort_eig(int nd, int lohi, int nec, double ref, double *lambda,
                  double *res)
{
  int i, j;

  if (lohi == 0) {
    // sort the eigenvalues according to their absolute residual values
    // to get those converged first
    dsort2a(nd, res, lambda);
    // sort the first nec eigenvalue in the order of lambda
    dsort2(nec, lambda, res);
  } else {
    // sort the eigenvalues and residual norms in ascending order of the
    // eigenvalues
    if (lohi == -2) {
      // around ref
      dsort2s(nd, ref, lambda, res);
      dsort2(nec, lambda, res);
    } else if (lohi == -3) {
      // larger than ref
      dsort2su(nd, ref, lambda, res);
      dsort2(nec, lambda, res);
    } else if (lohi == -4) {
      // smaller than ref
      dsort2sd(nd, ref, lambda, res);
      dsort2(nec, lambda, res);
    } else {
      dsort2(nd, lambda, res);

      if (lohi > 0) {
        // move the largest ones to the front (still ascending order)
        j = nd - nec;
        for (i = 0; i < nec; i++) {
          res[i] = res[j];
          lambda[i] = lambda[j];
          j++;
        }
      }
    }
  }
}

void trl_get_tvec(int nd, double *alpha, double *beta, int irot, int nrot,
                  double *rot, int nlam, double *lambda, double *yy,
                  int *iwrk, double *wrk, int lwrk, int *ierr)
{
  char notrans = 'N';
  int c__1 = 1;
  double zero = 0.0, one = 1.0;
  int i, j, k, ncol, ii, ioff;

  if (nlam <= 0) {
    *ierr = 0;
    return;
  }
  if (lwrk > 5 * nd) {
    *ierr = 0;
  } else {
    *ierr = -131;
    return;
  }

  /* set up IBLOCK and ISPLIT for calling dstein */
  for (i = 0; i < nd; i++) {
    iwrk[i] = 1;
    iwrk[nd + i] = nd;
  }
  F77_CALL(dstein)(&nd, alpha, beta,
                   &nlam, lambda,
                   iwrk, &iwrk[nd], yy, &nd, wrk,
                   &iwrk[2 * nd], &iwrk[3 * nd], ierr);

  if (*ierr != 0)
    error("TRL_GET_TVEC: dstein failed with error code %d\n", *ierr);

  /* apply the rotations to the IROT+1:IROT+NROT rows of YY
     generates results 'NCOL' columns at a time */
  if (nrot > 1) {
    ncol = lwrk / nrot;
    for (i = 1; i <= nlam; i += ncol) {
      j = imin2(nlam, i + ncol - 1);
      k = j - i + 1;
      if (k > 1) {
        trl_dgemm(&notrans, &notrans, nrot, k, nrot, one, rot,
                  nrot, &yy[(i - 1) * nd + irot], nd, zero, wrk,
                  nrot);
        for (ii = i - 1; ii < j; ii++) {
          ioff = (ii - i + 1) * nrot;
          memcpy(&yy[ii * nd + irot], &wrk[ioff],
                 nrot * sizeof(double));
        }
      } else {
        trl_dgemv(&notrans, nrot, nrot, one, rot, nrot,
                  &yy[(i - 1) * nd + irot], c__1, zero, wrk, c__1);
        memcpy(&yy[(i - 1) * nd + irot], wrk,
               nrot * sizeof(double));
      }
    }
  }
}

void trl_get_tvec_a(int nd, int kept, double *alpha, double *beta,
                    int nlam, double *lambda, double *yy, double *wrk,
                    int lwrk, int *iwrk, int *ierr) {
  char job = 'V';
  char upl = 'U';

  int i, j, i2, j2, ii;
  double tmp;

  /* fill yy with the projection matrix, then call DSYEV */
  if (nlam <= 0) {
    *ierr = 0;
    return;
  }
  if (lwrk >= nd + nd + nd) {
    *ierr = 0;
  } else {
    *ierr = -141;
    return;
  }
  memset(yy, 0, (nd * nd) * sizeof(double));
  j = 0;
  for (i = 0; i < nd; i++) {
    yy[j] = alpha[i];
    j += (nd + 1);
  }
  if (kept > 0) {
    memcpy(&yy[kept * nd], beta, kept * sizeof(double));
  }
  for (i = kept; i < nd - 1; i++) {
    yy[(i + 1) * (nd + 1) - 1] = beta[i];
  }

  F77_CALL(dsyev)(&job, &upl, &nd, yy, &nd, alpha, wrk, &lwrk, ierr);

  if (*ierr != 0)
    error("Error from dsyev: %d.\n", *ierr);

  if (nlam >= nd)
    return;

  /* reorder the eigenvectors
     both lambda(1:kept) and alpha are in ascending order */
  tmp = fmax2(alpha[nd - 1] - alpha[0],
              fmax2(fabs(alpha[nd - 1]), fabs(alpha[0])));
  tmp = DBL_EPSILON * tmp * nd;
  j = 0;
  i = 0;
  while (i < nlam) {
    // move j so that alpha(j) is within tmp distance away
    ii = j;
    j = nd - 1;
    while (ii < nd) {
      if (alpha[ii] < lambda[i] - tmp) {
        ii++;
      } else {
        j = ii;
        ii = nd;
      }
    }
    if (alpha[j] > lambda[i] + tmp) {
      *ierr = -143;
      return;
    }
    // identify the group size in lambda
    ii = i + 1;
    i2 = nlam - 1;
    while (ii < nlam) {
      if (lambda[ii] <= lambda[i] + tmp) {
        ii++;
      } else {
        i2 = ii - 1;
        ii = nd;
      }
    }
    // identify the group size in alpha
    ii = j + 1;
    j2 = nd - 1;
    while (ii < nd) {
      if (alpha[ii] <= lambda[i] + tmp) {
        ii++;
      } else {
        j2 = ii - 1;
        ii = nd;
      }
    }
    // assign the index values
    if (j2 == j && i2 == i) {
      iwrk[i] = j;
    } else if (j2 - j == i2 - i) {
      //iwrk(i:i2) = (/(ii, ii=j, j2)/)
      for (ii = i; ii < i2; ii++) {
        iwrk[ii] = j + ii - i;
      }
    } else if (j2 - j > i2 - i) {
      //j2 = j + i2 - i;
      //iwrk(i:i2) = (/(ii, ii=j, j2)/)
      for (ii = i; ii < i2; ii++) {
        iwrk[ii] = j + ii - i;
      }
    } else if (j2 < nd) {
      i2 = i + j2 - j;
      //iwrk(i:i2) = (/(ii, ii=j, j2)/)
      for (ii = i; ii < i2; ii++) {
        iwrk[ii] = j + ii - i;
      }
    } else {
      *ierr = -144;
      return;
    }
    i = i2 + 1;
    j = j2 + 1;
  }
  // perform the actual copying
  for (i = 0; i < nlam; i++) {
    // for( i=1; i<=nlam; i++ ) {
    j = iwrk[i];
    if (j > i) {
      alpha[i] = alpha[j];
      memcpy(&yy[i * nd], &yy[j * nd], nd * sizeof(double));
    }
  }
//
// .. end of trl_get_tvec_a_ ..
//
}

#define small(tmp,eps) (fabs(tmp) >= (eps) ? (eps)*fabs(tmp) : (eps)*(eps)*(anrm))
void trl_set_locking(int jnd, int nlam, double *lambda, double *res,
                     double *yy, int anrm, int *locked)
{
  double zero = 0.0;
  int i, j, ii, ioff;
  double tmp, eps;
  unsigned ti, tj;

  eps = DBL_EPSILON;
  i = 0;
  j = nlam - 1;
  ti = (fabs(1) <= small(2, 1));
  ti = (fabs(res[i]) <= small(lambda[i], eps));
  tj = (fabs(res[j]) <= small(lambda[j], eps));
  while (i < j) {
    if (ti) {
      // res[i] is very small, so lock ith lambda
      res[i] = zero;
      i = i + 1;
      if (i <= j) {
        ti = (fabs(res[i]) <= small(lambda[i], eps));
      } else {
        ti = 0;
      }
    } else {
      if (tj) {
        // res[i] (small ones) is still large,
        // but res[j] (big ones) is very small,
        // so swap res[j] with res[i], and lock res[i].
        // swap the eigenvectors accordingly.
        tmp = lambda[i];
        lambda[i] = lambda[j];
        lambda[j] = tmp;
        res[j] = res[i];
        res[i] = zero;
        ioff = (j - i) * jnd;
        for (ii = (i + 1) * jnd - jnd; ii < (i + 1) * jnd; ii++) {
          tmp = yy[ii];
          yy[ii] = yy[ii + ioff];
          yy[ii + ioff] = tmp;
        }
        i++;
        if (i <= j) {
          ti = (fabs(res[i]) <= small(lambda[i], eps));
        } else {
          ti = 0;
        }
      }
      j--;
      if (j > i) {
        tj = (fabs(res[j]) <= small(lambda[j], eps));
      } else {
        tj = 0;
      }
    }
  }
  if (ti) {
    *locked = i + 1;
  } else {
    *locked = i;
  }
}

void trl_ritz_vectors(int nrow, int lck, int ny, double *yy, int ldy,
                      double *vec1, int ld1, int m1, double *vec2,
                      int ld2, int m2, double *wrk, int lwrk)
{
  char notrans = 'N';
  double zero = 0.0, one = 1.0;
  int c__1 = 1;
  int i, j, k, stride, ii, jl1, jl2, il1, il2, kv1;

  // vec1*yy and vec2*yy where vec1 and vec2 are kept-locked
  // m1 number of kept in vec1
  // m2 number of kept in vec2
  if (lck <= m1) {
    // all locked are in vec1
    il1 = lck + 1;
    jl1 = m1 - lck;
    il2 = 1;
    jl2 = m2;
  } else {
    // all kept in vec1 are locked
    il1 = m1 + 1;
    jl1 = 0;
    il2 = lck - m1 + 1;
    jl2 = m1 + m2 - lck;
  }

  if (jl1 == 0 && jl2 == 0)
    return;

  kv1 = imin2(m1 - il1 + 1, ny);
  memset(wrk, 0, lwrk * sizeof(double));
  if (ny > 1) {
    stride = lwrk / ny;
    for (i = 0; i < nrow; i += stride) {
      j = imin2(nrow - 1, i + stride - 1);
      k = j - i + 1;

      if (jl1 > 0) {
        // compute wrk = vec1(i:j,:)*yy
        // (Note the leading dimension of vec1 is ld1. This effectively shift
        // the vec1(i:j) to the top of the matrix.
        trl_dgemm(&notrans, &notrans, k, ny, jl1, one,
                  &vec1[(il1 - 1) * ld1 + i], ld1, yy, ldy, zero,
                  wrk, k);
      } else {
        memset(wrk, 0, lwrk * sizeof(double));
      }

      if (jl2 > 0) {
        trl_dgemm(&notrans, &notrans, k, ny, jl2, one,
                  &vec2[(il2 - 1) * ld2 + i], ld2, &yy[jl1], ldy,
                  one, wrk, k);
      }

      for (ii = 0; ii <= kv1 - 1; ii++) {
        memcpy(&vec1[(ii + il1 - 1) * ld1 + i], &wrk[ii * k],
               k * sizeof(double));
      }

      for (ii = 0; ii <= (ny - kv1 - 1); ii++) {
        memcpy(&vec2[(ii + il2 - 1) * ld2 + i],
               &wrk[(kv1 + ii) * k], k * sizeof(double));
      }

    }
  } else if (ny == 1) {
    stride = lwrk;
    for (i = 0; i < nrow; i += stride) {
      j = imin2(nrow - 1, i + stride - 1);
      k = j - i + 1;
      if (jl1 > 0) {
        trl_dgemv(&notrans, k, jl1, one,
                  &vec1[(il1 - 1) * ld1 + i], ld1, yy, c__1, zero,
                  wrk, c__1);
        if (jl2 > 0) {
          trl_dgemv(&notrans, k, jl2, one,
                    &vec2[(il2 - 1) * ld2 + i], ld2, &yy[jl1],
                    c__1, one, wrk, c__1);
        }
      } else {
        trl_dgemv(&notrans, k, jl2, one,
                  &vec2[(il2 - 1) * ld2 + i], ld2, &yy[jl1], c__1,
                  zero, wrk, c__1);
      }
      if (kv1 > 0) {
        memcpy(&vec1[(il1 - 1) * ld1 + i], wrk,
               k * sizeof(double));
      } else {
        memcpy(&vec2[(il2 - 1) * ld2 + i], wrk,
               k * sizeof(double));
      }
    }
  }
}

int trl_cgs(trl_info * info, int nrow, double *v1, int ld1, int m1,
            double *v2, int ld2, int m2, double *rr, double *rnrm,
            double *alpha, int *north, double *wrk)
{
  long c__1 = 1;
  char notrans = 'N';
  double one = 1.0, zero = 0.0;
  const int maxorth = 3;
  double d__1;
  int mpicom, myid, i, j, k, nold, irnd, cnt, ierr = 0;
  double tmp, old_rnrm;

  mpicom = info->mpicom;
  myid = info->my_pe;
  nold = m1 + m2;
  if (ld1 < nrow || (ld2 < nrow && m2 > 0)) {
    return -201;
  }
  irnd = 0;
  ierr = 0;
  if (nold > 0) {
    cnt = 0;
    while (cnt <= maxorth) {
      // compute [v1 v2]'*rr=wrk
      trl_g_dot(mpicom, nrow, v1, ld1, m1, v2, ld2, m2, rr, wrk);
      if (m1 > 1) {
        d__1 = -one;
        trl_dgemv(&notrans, nrow, m1, d__1, v1, ld1, wrk, c__1,
                  one, rr, c__1);
      } else if (m1 == 1) {
        d__1 = -wrk[0];
        trl_daxpy(nrow, d__1, v1, c__1, rr, c__1);
      }
      if (m2 > 1) {
        d__1 = -one;
        trl_dgemv(&notrans, nrow, m2, d__1, v2, ld2, &wrk[m1],
                  c__1, one, rr, c__1);
      } else if (m2 == 1) {
        d__1 = -wrk[nold - 1];
        trl_daxpy(nrow, d__1, v2, c__1, rr, c__1);
      }

      (*north)++;
      cnt = cnt + 1;
      tmp = trl_ddot(nold, wrk, c__1, wrk, c__1);
      wrk[0] = trl_ddot(nrow, rr, c__1, rr, c__1);
      trl_g_sum(mpicom, 1, wrk, &wrk[1]);
      *rnrm = sqrt(wrk[0]);
      old_rnrm = sqrt(wrk[0] + tmp);

      /* decisions about whether to re-orthogonalize is based on
         relative size of tmp and wrk(1) (R norm square) */
      if (DBL_EPSILON * wrk[0] > tmp) {
        /* no need for more orthogonalization */
        cnt = maxorth + 1;
      } else if (((cnt > 1 && !
                   (tmp <=
                    info->ntot * DBL_EPSILON * DBL_EPSILON *
                    (wrk[0] + tmp))) || !(wrk[0] > DBL_MIN))
                 && irnd < maxorth) {
        // the input vector is so nearly linear dependent on the
        // existing vectors, we need to perturb it in order to
        // generate a new vector that is orthogonal to the existing
        // ones
        // the perturbation is done in two stages:
        // -- perturbing only one number
        // -- call random_number to generate a whole random vector
        cnt = 0;
        irnd++;
        info->nrand++;
        GetRNGstate();
        if (irnd == 1 && *rnrm > 0.0 &&
            *rnrm > DBL_EPSILON * old_rnrm) {
          /* old version:  modify only one element
             new version:  modify (nrow * epsilon * rnrm / old_rnrm ) elements */
          tmp = unif_rand();
          i = (int) (nrow * tmp);
          k = i +  (int) (fmax2(1.0,
                                (nrow * DBL_EPSILON * old_rnrm / *rnrm)));
          for (j = i; j < k; j++) {
            tmp = unif_rand();
            while (fabs(tmp - 0.5) <= DBL_EPSILON) {
              tmp = unif_rand();
            }
            rr[j] += (*rnrm) * (tmp - 0.5);
          }
        } else {
          /* fill with random numbers produced by intrinsic function */
          for (i = 0; i <= myid; i++) {
            tmp = unif_rand();
          }
          i = (int) (nrow * tmp);
          tmp = unif_rand();
          j = (int) (nrow * tmp);
          if (i < j) {
            for (k = i; k <= j; k++) {
              rr[k] = unif_rand();
            }
          } else if (i > j) {
            for (k = j; k <= i; k++) {
              rr[k] = unif_rand();
            }
          } else {
            for (k = 0; k < nrow; k++) {
              rr[k] = unif_rand();
            }
          }
        }
        PutRNGstate();
        /* rr = rr + rr + Cshift(rr, 1) + Cshift(rr, -1) */
        smoothrr(nrow, rr);
      }
    }
    /* failed to reduce the dot-products between the new vector
       and the old vectors to small number -- usually an indication of
       problem with orthogonality of the old vectors */
    if (!(wrk[0] >= tmp))
      ierr = -203;
  }

  /* normalization */
  if (ierr == 0) {
    if (*rnrm > DBL_MIN) {
      tmp = one / *rnrm;
      trl_dscal(nrow, tmp, rr, c__1);
    } else {
      return -204;
    }
  }
  if (irnd > 0)
    *rnrm = zero;
  return ierr;
}

int trl_check_dgen(trl_info * info, int jnd, double *lambda,
                   double *res) {
  /*
    if( info->nec >= info->ned ) {
    if( (info->lohi == -1 || info->lohi == 0 ) &&
    (lambda[info->nec]  +res[info->nec]   > lambda[info->nec+1] &&
    lambda[info->nec+1]-res[info->nec+1] < lambda[info->nec]) ) {
    return 1;
    }
    if( (info->lohi == 1 || info->lohi == 0 ) &&
    (lambda[jnd-info->nec+1]-res[jnd-info->nec+1] < lambda[jnd-info->nec] &&
    lambda[jnd-info->nec]  +res[jnd-info->nec]   > lambda[jnd-info->nec+1]) ) {
    return 1;
    }
    }
  */
  return -1;
}
