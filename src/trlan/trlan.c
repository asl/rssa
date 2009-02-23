/*
  ZTRLan routine (version 1.0)
  Lawrence Berkeley National Lab.
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "trlan.h"
#include "trlan_i.h"
#include "trlcore_i.h"
#include "trlaux_i.h"
#include "trl_comm_i.h"

void trl_clear_counter(trl_info * info, int nrow, int mev, int lde) {
  int ntmp;

  info->stat = 0;
  if (nrow != info->nloc || nrow > info->ntot)
    error("TRLAN: info not setup for this problem.\n       Please reset info before calling TRLAN.\n");

  if (info->nec < 0)
    info->nec = 0;

  if (lde < nrow)
    error("TRLAN: leading dimension of EVEC to small.\n");

  if (info->tol >= 1.0) {
    info->tol = sqrt(DBL_EPSILON);
  } else if (info->tol <= DBL_MIN) {
    info->tol = DBL_EPSILON;
  }

  if (info->ned + info->ned >= info->ntot) {
    warning("TRLAN: info->ned (%d) is large relative to the matrix dimension (%d)\n",
            info->ned, info->ntot);
    warning(" **    It is more appropriate to use LAPACK dsyev/ssyev.\n");
    if (info->ned > info->ntot) {
      info->ned = imin2(info->ntot - 1, info->maxlan - 3);
      warning("TRLAN: ** reduced ned to %d **\n", info->ned);
    }
  }
  if (mev < info->ned)
    error("TRLAN: array EVAL and EVEC can not hold wanted number of eigenpairs.\n");

  if (info->ntot < 10)
    error("TRLAN is not designed to work with such a small matrix(%dx%d).  Please use LAPACK or EISPACK instead.\n",
          info->ntot, info->ntot);

  info->nrand = info->stat;
  info->stat = trl_sync_flag(info->mpicom, info->nrand);

  /* decide what is a good maximum basis size to use */
  if (info->maxlan < info->ned + 3) {
    info->maxlan = info->ned + imin2(info->ned, 20) +
      (int)(log((double)info->ntot));
    info->maxlan = imin2(info->maxlan, info->ntot);
    warning("TRLAN: ** reset maxlan to %d! **\n", info->maxlan);
  }

  if (info->maxlan < mev) {
    ntmp = imin2(info->ntot, imax2(100 + info->ned, 10 * info->ned));
    if (mev < ntmp) {
      info->maxlan = mev;
    } else {
      info->maxlan = ntmp;
    }
  }
  if (info->maxlan < 5)
    error("TRLAN must have at least 5 basis vectors, it is currently %d.\n",
          info->maxlan);

  /* clear regular counters */
  info->tmv = -1;
  info->klan = imin2(info->maxlan, info->ntot);
  if (info->restart >= 7) {
    info->klan = imin2(info->maxlan,
                       imax2(100,
                             imin2(info->klan, 2 * (info->ned))));
  }
  info->locked = info->nec;
  info->matvec = 0;
  info->nloop = 0;
  info->north = 0;
  info->nrand = 0;
  info->tick_t = 0.0;
  info->clk_op = 0;
  info->tick_o = 0.0;
  info->clk_orth = 0;
  info->tick_h = 0.0;
  info->clk_res = 0;
  info->tick_r = 0.0;
  info->clk_in = 0;
  info->clk_out = 0;
  info->wrds_in = 0;
  info->wrds_out = 0;
  info->avgm = 0.0;
  return;
}

void trl_print_setup(trl_info * info, int lbas, int lmis, int lwrk)
{
  if (info->lohi > 0) {
    fprintf(info->log_fp,
            "TRLAN is to compute %6d largest eigenpair(s).\n",
            info->ned);
  } else if (info->lohi < 0) {
    fprintf(info->log_fp,
            "TRLAN is to compute %6d smallest eigenpair(s).\n",
            info->ned);
  } else {
    fprintf(info->log_fp,
            "TRLAN is to compute %6d first converged eigenpair(s).\n",
            info->ned);
  }
  fprintf(info->log_fp,
          "Problem dimension: %9d (PE:%4d) %12d (Global)\n", info->nloc,
          info->my_pe, info->ntot);
  fprintf(info->log_fp, "Maximum basis size:                   %10d\n",
          info->maxlan);
  fprintf(info->log_fp, "Dynamic restarting scheme:            %10d\n",
          info->restart);
  fprintf(info->log_fp, "Maximum applications of the operator: %10d\n",
          info->maxmv);
  fprintf(info->log_fp, "Relative convergence tolerance: %10e\n",
          info->tol);
  /* initial guess */
  if (info->guess == 1) {
    fprintf(info->log_fp, "User provided the starting vector.\n");
  } else if (info->guess == 0) {
    fprintf(info->log_fp, "TRLAN uses [1,1,...] as starting vctor.\n");
  } else if (info->guess < 0) {
    fprintf(info->log_fp,
            "TRLAN generates a random starting vector.\n");
  } else if (info->oldcpf == 0 || strlen(info->oldcpf) == 0) {
    fprintf(info->log_fp,
            "Restarting with existing checkpoint files %s ####\n",
            info->oldcpf);
  } else {
    fprintf(info->log_fp,
            "Restarting with existing checkpoint files %s ####\n",
            info->cpfile);
  }
  if (info->cpflag > 0) {
    fprintf(info->log_fp,
            "TLRAN will write about %d sets of checkpointing files %s ####.\n",
            info->cpflag, info->cpfile);
  }
  /* print the workspace size parameters */
  fprintf(info->log_fp, "(required) array BASE size is %d\n", lbas);
  fprintf(info->log_fp, "(required) array MISC size is %d\n", lmis);
  if (lwrk > 0) {
    fprintf(info->log_fp,
            "Caller has supplied a work array with %d elements.\n",
            lwrk);
  } else {
    fprintf(info->log_fp, "Caller did not supply work array.\n");
  }
}

void
trl_ritz_projection(trl_matprod op,
                     trl_info * info, int lde, int mev, double *evec,
                     double *eres, int lwrk, double *wrk, double *base,
                     void *lparam)
{
  char trans = 'T', notrans = 'N', upper = 'U', job = 'V';
  double one = 1.0, zero = 0.0;
  int i__1 = 1;
  int i, j, ierr, nev, nsqr, nrow, iuau, irvv, lwrk2;
  double d__1;
  double *rvv, *uau, *wrk2, *avec;

  nrow = info->nloc;
  if (info->nec > 0) {
    nev = info->nec + 1;
  } else {
    nev = imin2(info->ned, mev - 1);
    if (info->lohi != 0)
      nev++;
  }
  nsqr = nev * nev;
  if (lwrk < 0) {
    lwrk = 0;
  }
  if (base != NULL) {
    avec = base;
  } else if (mev > nev) {
    avec = &evec[(mev - 1) * nrow];
  } else {
    avec = Calloc(nrow, double);
  }
  if (info->verbose >= 0) {
    if (info->log_fp == NULL) {
      trl_reopen_logfile(info);
    }
    fprintf(info->log_fp,
            "TRLAN performing a separate Rayleigh-Ritz project for %d vectors.",
            nev);
  }
  /* memory allocation -- need 3*nev*nev elements, will allocate them     */
  /* in two consecutive blocks, uau(nev*nev), rvv(2*nev*nev)              */
  /* in actual use, rvv is further split in two until the last operation  */
  iuau = nsqr;
  irvv = nsqr + nsqr;
  if (lwrk >= iuau + irvv) {
    uau = wrk;
    rvv = &wrk[nsqr];
    wrk2 = &wrk[nsqr + nsqr];
    lwrk2 = lwrk - nsqr - nsqr;
  } else if (lwrk >= irvv) {
    rvv = wrk;
    wrk2 = &wrk[nsqr];
    lwrk2 = lwrk - nsqr;
    uau = Calloc(nsqr, double);
  } else if (lwrk >= iuau) {
    uau = wrk;
    rvv = Calloc(nsqr + nsqr, double);
    wrk2 = &rvv[nsqr];
    lwrk2 = nsqr;
  } else {
    uau = Calloc(nsqr, double);
    rvv = Calloc(nsqr + nsqr, double);
    wrk2 = &rvv[nsqr];
    lwrk2 = nsqr;
  }
  /* step (1) : V'*V ==> G */

  trl_dgemm(&trans, &notrans, nev, nev, nrow, one, evec, lde, evec, lde,
            zero, rvv, nev);
  trl_g_sum(info->mpicom, nsqr, rvv, wrk2);

  /* step (2) : Choleskey factorization of G */
  F77_CALL(dpotrf)(&upper, &nev, rvv, &nev, &ierr);
  if (ierr != 0) {
    info->stat = -234;
    goto end;
  }
  /* step (3) : compute H_1 = V'*A*V                              */
  /* use the first nrow elements of avec to store the results of  */
  /* matrix-vector multiplication                                 */
  memset(wrk2, 0, lwrk2 * sizeof(double));
  for (i = 1; i <= nev; i++) {
    op(&nrow, &i__1, &evec[(i - 1) * nrow], &lde, avec, &nrow, lparam);
    trl_dgemv(&trans, nrow, i, one, evec, lde, avec, i__1, zero,
              &wrk2[(i - 1) * nev], i__1);
  }
  trl_g_sum(info->mpicom, nsqr, wrk2, uau);
  for (i = 1; i < nev; i++) {
    for (j = 0; j < i; j++) {
      wrk2[i + j * nev] = wrk2[(i - 1) * nev + j];
    }
  }
  /* compute solution of R^T H_2 = H_1 */
  F77_CALL(dtrtrs)(&upper, &trans, &notrans,
                   &nev, &nev, rvv, &nev, wrk2, &nev,
                   &ierr);
  if (ierr != 0) {
    info->stat = -235;
    goto end;
  }
  /* compute solution of R^T H = H_2^T */
  for (i = 1; i < nev; i++) {
    for (j = 0; j < nev; j++) {
      uau[i + j * nev] = wrk2[(i - 1) * nev + j];
    }
  }
  F77_CALL(dtrtrs)(&upper, &trans, &notrans,
                   &nev, &nev, rvv, &nev, uau, &nev,
                   &ierr);
  if (ierr != 0) {
    info->stat = -236;
    goto end;
  }
  /* solve the small symmetric eigenvalue problem */
  F77_CALL(dsyev)(&job, &upper, &nev, uau, &nev, eres, wrk2, &nsqr, &ierr);
  if (ierr != 0) {
    info->stat = -237;
    goto end;
  }
  /* solve R Y = Y to prepare for multiplying with V */
  F77_CALL(dtrtrs)(&upper, &notrans, &notrans,
                   &nev, &nev, rvv, &nev, uau, &nev,
                   &ierr);
  if (ierr != 0) {
    info->stat = -238;
    goto end;
  }
  /* call trl_ritz_vector to do the final multiplication */
  if (lwrk >= 3 * nsqr) {
    wrk2 = &wrk[nsqr];
  } else if (lwrk >= nsqr + nsqr) {
    wrk2 = wrk;
  } else {
    wrk2 = rvv;
  }
  i = lwrk2;
  trl_ritz_vectors(nrow, 0, nev, uau, nev, evec, lde, nev, avec, nrow,
                    0, wrk2, i);
  /* compute the residual norms */
  for (i = 0; i < nev; i++) {
    op(&nrow, &i__1, &evec[i * nrow], &lde, avec, &nrow, lparam);
    d__1 = eres[i];
    trl_daxpy(nrow, d__1, &evec[i * nrow], i__1, avec, i__1);
    eres[nev + i] = trl_ddot(nrow, avec, i__1, avec, i__1);
  }
  trl_g_sum(info->mpicom, nev, &eres[nev], avec);
  for (i = nev; i < nev + nev; i++) {
    if (eres[i] > 0.0) {
      eres[i] = sqrt(eres[i]);
    } else {
      eres[i] = -DBL_MAX;
    }
  }
  if (info->lohi < 0) {
    for (i = nev - 1; i < nev + nev - 2; i++) {
      eres[i] = eres[i + 1];
    }
  } else if (info->lohi > 0) {
    for (i = 0; i < nev - 1; i++) {
      eres[i] = eres[i + 1];
      memcpy(&evec[i * nrow], &evec[(i + 1) * nrow], nrow);
    }
    for (i = nev - 1; i < nev + nev - 2; i++) {
      eres[i] = eres[i + 2];
    }
  }
end:
  if (lwrk < iuau) {
    Free(uau);
    Free(rvv);
  } else if (lwrk < irvv) {
    Free(rvv);
  } else if (lwrk < iuau + irvv) {
    Free(uau);
  }
}


void trlan(trl_matprod op,
           trl_info * info, int nrow, int mev, double *eval,
           double *evec, int lde, int lwrk, double *wrk, void *lparam)
{
  clock_t clk1;
  int ii, nbas, nmis, ibas, imis, ldb, lwrk0;
  double *base, *misc;

  imis = -1; /* if this routine allocated misc, imis will be 0 */
  ibas = -1; /* if this routine allocated base, ibas will be 0 */
  clk1 = clock();
  info->clk_tot = clk1;
  if (info->ned > mev) {
    warning("info->ned (%d) is larger than mev (%d) reducing info->ned to %d\n",
            info->ned, mev, mev);
    info->ned = mev;
  }
  /* there is nothing to do if there is no more eigenvalue to compute */
  if (info->ned <= info->nec || info->ned <= 0)
    goto end;

  lwrk0 = lwrk;
  info->stat = 0;
  ldb = ((nrow + 3) / 4) * 4;
  if ((ldb % 4096) == 0)
    ldb = ldb + 4;
  trl_clear_counter(info, nrow, mev, lde);
  if (info->stat != 0)
    goto end;
  /*
    Internally, the workspace is broken into two parts
    one to store (maxlan+1) Lanczos vectors, and the other to
    store all others (size maxlan*(maxlan+ned+14))
    The next If-block decides how the two arrays are mapped.
  */
  nbas = imax2(1, info->maxlan - mev + 1);
  ii = nbas * ldb;
  nmis = info->maxlan * (info->maxlan + 10);
  if (lwrk0 >= imin2(ii, nmis)) {
    /* use wrk either as base or misc or both depending its size    */
    if (lwrk0 >= ii + nmis) {
      /* WRK is large enough for both arrays */
      base = wrk;
      misc = &wrk[ii];
      nmis = lwrk0 - ii;
      /* printf( "\n\n ******** Large enough workspace ********* \n\n" ); */
    } else if (lwrk0 >= imax2(ii, nmis)) {
      /* associate the larger one of base and misc to WRK */
      if (ii >= nmis) {
        base = wrk;
        misc = Calloc(nmis, double);
        imis = 0;
      } else {
        misc = wrk;
        nmis = lwrk0;
        base = Calloc(ii, double);
        ibas = 0;
      }
    } else if (ii <= nmis) {
      /* base is smaller, associate base with WRK */
      base = wrk;
      misc = Calloc(nmis, double);
      imis = 0;
    } else {
      /* misc is smaller, associate misc with WRK */
      misc = wrk;
      nmis = lwrk0;
      base = Calloc(ii, double);
      ibas = 0;
    }
  } else {
    /* have to allocate both base and misc */
    base = Calloc(ii, double);
    ibas = 0;
    misc = Calloc(nmis, double);
    imis = 0;
  }
  memset(base, 0, ii * sizeof(double));
  memset(misc, 0, nmis * sizeof(double));
  /* make sure every process is successful so far */
  ii = trl_sync_flag(info->mpicom, info->stat);
  info->stat = ii;
  if (ii != 0)
    goto end;
  /* open log and checkpoint files */
  trl_open_logfile(info);
  /* trl_open_cptfile(info); */
  if (info->verbose > 0) {
    trl_time_stamp(info->log_fp);
    trl_print_setup(info, nbas * ldb, nmis, lwrk0);
  }

  /* call trlanczos to do the real work  */
  /* printf( "calling trlanczso (%d)\n",info->cpflag ); */
  trlanczos(op, info, nrow, mev, eval, evec, lde, base, ldb, nbas, misc,
             nmis, lparam);
  /* printf( " ** out of trlanczos (locked=%d) **\n",info->locked ); */

  /* close log and checkpoint files */
  trl_close_logfile(info);
  if (lwrk0 >= mev) {
    memmove(wrk, misc, mev * sizeof(double));
  } else {
    memmove(wrk, misc, lwrk0 * sizeof(double));
  }

  /* DONE, reclaim the space allocated */
end:
  if (imis == 0)
    Free(misc);
  if (ibas == 0)
    Free(base);
  clk1 = clock();
  if (clk1 < info->clk_tot) {
    info->tick_t +=
      (info->clk_max -
       info->clk_tot) / (double) (info->clk_rate);
    info->tick_t +=
      (info->clk_max + clk1) / (double) (info->clk_rate);
    info->clk_tot = 0;
  } else if (info->clk_tot < 0 && clk1 >= 0) {
    info->tick_t -= info->clk_tot / (double) (info->clk_rate);
    info->tick_t += clk1 / (double) (info->clk_rate);
    info->clk_tot = 0;
  } else {
    info->tick_t  += (clk1 - info->clk_tot) / (double) (info->clk_rate);
    info->clk_tot  = 0;
  }

  return;
}

void trl_set_restart(trl_info * info, double rfact)
{
  info->rfact = rfact;
}

void trl_set_debug(trl_info * info, int msglvl, char *filename)
{
  info->verbose = msglvl;
  if (filename != NULL) {
    strcpy(info->log_file, filename);
    if (msglvl > 0) {
      printf
        ("TRLAN will write diagnostic messages to files with prefix %s.\n",
         info->log_file);
    }
  }
}

void trl_set_checkpoint(trl_info * info, int cpflag, char *file)
{
  info->cpflag = cpflag;
  if (file != NULL) {
    strcpy(info->cpfile, file);
  }
}

void trl_set_iguess(trl_info * info, int nec, int iguess, int nopts,
                    char *cpf)
{
  /* assign nec and iguess flags to info */
  info->nec = nec;
  info->guess = iguess;
  if (strlen(info->oldcpf) > 0 && info->guess > 1) {
    /* check to make sure the files exist */
    trl_pe_filename(TRLAN_STRING_LEN, cpf, info->oldcpf, info->my_pe,
                     info->npes);
    if ((info->cpt_fp = fopen(cpf, "r")) != NULL) {
      if (fclose(info->cpt_fp) != 0) {
        info->stat = -9;
      }
    } else {
      info->stat = -8;
    }
    info->stat = trl_sync_flag(info->mpicom, info->stat);
  } else {
    info->stat = 0;
  }
}

void trl_print_info(trl_info * info) {
  double tmp1[12], tmp2[12];
  int i;
  double t_tot, t_op, t_orth, t_res, t_in, t_out, rinv, r_in, r_out;

  if (info->clk_rate > 0) {
    rinv = 1.0 / (double) (info->clk_rate);
  } else {
    /* get clock rate */
    rinv = 1.0 / CLOCKS_PER_SEC;
  }
  t_op   = info->tick_o + info->clk_op   * rinv;
  t_tot  = info->tick_t + info->clk_tot  * rinv;
  t_res  = info->tick_r + info->clk_res  * rinv;
  t_orth = info->tick_h + info->clk_orth * rinv;
  t_in   = info->clk_in * rinv;
  t_out  = info->clk_out * rinv;
  /* printf( "tick_t=%e clk_tot=%d\n",info->tick_t,info->clk_tot ); */
  if (info->clk_in > 0) {
    r_in = 8.0 * info->wrds_in;
  } else {
    r_in = 0;
  }
  if (info->clk_out > 0) {
    r_out = 8.0 * info->wrds_out;
  } else {
    r_out = 0;
  }
  tmp2[0] = t_tot;
  tmp2[1] = t_op;
  tmp2[2] = t_orth;
  tmp2[3] = t_res;
  tmp2[4] = t_in;
  tmp2[5] = t_out;
  tmp2[10] = r_in;
  tmp2[11] = r_out;
  /* printf("print info\n" ); */
  /* printf("calling g_sum\n" ); */
  trl_g_sum(info->mpicom, 12, tmp2, tmp1);
  if (info->log_fp == NULL) {
    trl_reopen_logfile(info);
  }
  trl_time_stamp(info->log_fp);
  /* printf("printing\n"); */
  if (info->npes > 1) {
    fprintf(info->log_fp,
            "TRLAN execution summary (exit status = %d) on PE %d\n",
            info->stat, info->my_pe);
  } else {
    fprintf(info->log_fp,
            "TRLAN execution summary (exit status =%d)\n", info->stat);
  }
  if (info->lohi > 0) {
    fprintf(info->log_fp,
            "Number of LARGEST eigenpairs      %10d (computed) %11d (wanted)\n",
            info->nec, info->ned);
  } else if (info->lohi < 0) {
    fprintf(info->log_fp,
            "Number of SMALLEST eigenpairs    %10d (computed) %11d (wanted)\n",
            info->nec, info->ned);
  } else {
    fprintf(info->log_fp,
            "Number of EXTREME eigenpairs     %10d (computed) %11d (wanted)\n",
            info->nec, info->ned);
  }
  fprintf(info->log_fp,
          "Times the operator is applied:   %10d (MAX: %16d )\n",
          info->matvec, info->maxmv);
  fprintf(info->log_fp,
          "Problem size:                    %10d (PE: %4d) %11d (Global)\n",
          info->nloc, info->my_pe, info->ntot);
  fprintf(info->log_fp,
          "Convergence tolerance:           %10.3e (rel) %16.3e (abs)\n",
          info->tol, info->tol * info->anrm);
  fprintf(info->log_fp, "Maximum basis size:              %10d\n",
          info->maxlan);
  fprintf(info->log_fp, "Restarting scheme:               %10d\n",
          info->restart);
  fprintf(info->log_fp, "Number of re-orthogonalizations: %10d\n",
          info->north);
  fprintf(info->log_fp, "Number of (re)start loops:       %10d\n",
          info->nloop);
  if (info->nrand > 0) {
    fprintf(info->log_fp, "Number of random vectors used:   %10d\n",
            info->nrand);
  }
  if (info->npes > 1) {
    fprintf(info->log_fp, "Number of MPI processes:         %10d\n",
            info->npes);
  }
  fprintf(info->log_fp, "Number of eigenpairs locked:     %10d\n",
          info->locked);
  fprintf(info->log_fp, "time in OP:            %12.4e sec\n", t_op);
  fprintf(info->log_fp, "time in orth:          %12.4e sec\n", t_orth);
  fprintf(info->log_fp, "time in restarting:    %12.4e sec\n", t_res);
  fprintf(info->log_fp, "total time in TRLAN:   %12.4e sec\n", t_tot);

  /*
    if (info->verbose > 0 && info->log_fp != info->log_fp) {
    fprintf( info->log_fp, "Debug infomation written to files %s ####\n",info->log_file );
    }
  */
  if (info->guess > 1 && info->wrds_in > 0) {
    if (strlen(info->oldcpf) <= 0) {
      fprintf(info->log_fp,
              "TRLAN restarted with checkpoint files %s ####\n",
              info->oldcpf);
    } else {
      fprintf(info->log_fp,
              "TRLAN restarted with checkpoint files %s ####\n",
              info->cpfile);
    }
    fprintf(info->log_fp,
            "Bytes read   %12.5e, Time(sec): %12.5e, Rate(B/s): %12.5e\n",
            r_in, t_in, r_in / t_in);
  }
  if (info->clk_out > 0 && info->wrds_out > 0) {
    fprintf(info->log_fp, "Checkpoint files are %s ####\n",
            info->cpfile);
    fprintf(info->log_fp,
            "Bytes read   %12.5e, Time(sec): %12.5e, Rate(B/s): %12.5e\n",
            r_out, t_out, r_out / t_out);
  }
  if (info->npes > 1) {
    /* write global performance information */
    rinv = 1.0 / info->npes;
    for (i = 0; i < 12; i++) {
      tmp1[i] = tmp1[i] * rinv;
    }
    for (i = 0; i < 6; i++) {
      if (tmp1[i] > 0) {
        tmp1[i + 6] = tmp1[i + 6] / tmp1[i];
      } else {
        tmp1[i + 6] = 0.0;
      }
    }
    if (tmp1[4] == tmp1[5] && tmp1[4] == 0) {
      fprintf(info->log_fp,
              " -- Global summary -- \n" );
      fprintf(info->log_fp,
              "                       Overall,\t\t  MATVEC,\t  Re-orth,\t  Restart,\n");
      fprintf(info->log_fp,
              "Time(ave)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
              tmp1[0], tmp1[1], tmp1[2], tmp1[3]);
      fprintf(info->log_fp,
              "Rate(tot)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
              tmp1[6], tmp1[7], tmp1[8], tmp1[9]);
    } else {
      fprintf(info->log_fp,
              " -- Global summary -- \n" );
      fprintf(info->log_fp,
              "                       Overall,\t\t  MATVEC,\t  Re-orth,\t  Restart,\t  Read,\t  Write\n");
      fprintf(info->log_fp,
              "Time(ave)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
              tmp1[0], tmp1[1], tmp1[2], tmp1[3], tmp1[4], tmp1[5]);
      fprintf(info->log_fp,
              "Rate(tot)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
              tmp1[6], tmp1[7], tmp1[8], tmp1[9], tmp1[10],
              tmp1[11]);
    }
  }
  trl_close_logfile(info);
  return;
}

void trl_terse_info(trl_info * info, FILE * iou) {
  int rate;
  double t_tot, t_op, t_orth, t_res;

  if (iou == NULL) {
    if (info->log_fp == NULL) {
      iou = stdout;
    } else {
      iou = info->log_fp;
    }
  }
  if (info->clk_rate > 0) {
    t_op = info->tick_o + info->clk_op / (double) (info->clk_rate);
    t_tot = info->tick_t + info->clk_tot / (double) (info->clk_rate);
    t_res = info->tick_r + info->clk_res / (double) (info->clk_rate);
    t_orth = info->tick_h + info->clk_orth / (double) (info->clk_rate);
  } else {
    /* get clock rate */
    rate = CLOCKS_PER_SEC;
    t_op = info->tick_o + info->clk_op / (double) (rate);
    t_tot = info->tick_t + info->clk_tot / (double) (rate);
    t_res = info->tick_r + info->clk_res / (double) (rate);
    t_orth = info->tick_h + info->clk_orth / (double) (rate);
  }
  if (info->lohi > 0) {
    fprintf(iou,
            "MAXLAN:%10d, Restart:%10d,   NED: + %7d,      NEC:%10d\n",
            info->maxlan, info->restart, info->ned, info->nec);
  } else if (info->lohi < 0) {
    fprintf(iou,
            "MAXLAN:%10d, Restart:%10d,   NED: - %7d,      NEC:%10d\n",
            info->maxlan, info->restart, info->ned, info->nec);
  } else {
    fprintf(iou,
            "MAXLAN:%10d, Restart:%10d,   NED: 0 %7d,      NEC:%10d\n",
            info->maxlan, info->restart, info->ned, info->nec);
  }
  fprintf(iou,
          "MATVEC:%10d,  Reorth:%10d, Nloop:   %7d,  Nlocked:%10d\n",
          info->matvec, info->north, info->nloop, info->locked);
  if (t_tot > 0.001 && fmax2(t_tot, fmax2(t_op, fmax2(t_res, t_orth))) < 1000) {
    fprintf(iou,
            "Ttotal:%10.6f,    T_op:%10.6f, Torth:%10.6f,   Tstart:%10.6f\n",
            t_tot, t_op, t_orth, t_res);
  } else {
    fprintf(iou,
            "Ttotal:%10.3e,    T_op:%10.3e, Torth:%10.3e,   Tstart:%10.3e\n",
            t_tot, t_op, t_orth, t_res);
  }
}

void
trl_check_ritz(trl_matprod op,
               trl_info * info, int nrow, int ncol, double *rvec,
               int ldrvec, double *alpha, int *check, double *beta,
               double *eval, int lwrk, double *wrk, void *lparam)
{
  long c__1 = 1;
  int i__1 = 1;
  double d__1;
/*
  aq -- store the result of matrix-vector multiplication, size nrow
  rq -- store the Rayleigh-Quotient and the residual norms
  gsumwrk -- workspace left over for trl_g_sum to use dimension of the input arrays
*/
  double *aq, *rq, *gsumwrk, *res, *err;
  int i, aqi, rqi, gsumwrki, icheck;
  double gapl, gapr;

  if (ncol <= 0)
    return; /* there is nothing to do */

  /* figure out whether it is necessary to allocate new workspace */
  *check = 0;
  aqi = 0;
  rqi = 0;
  gsumwrki = 0;
  if (lwrk > nrow + (4 * ncol)) {
    aq = &wrk[0];
    rq = &wrk[nrow];
    gsumwrk = &wrk[nrow + (3 * ncol)];
  } else if (lwrk >= (nrow + ncol)) {
    aq = &wrk[0];
    gsumwrk = &wrk[nrow];
    rq = Calloc(3 * ncol, double);
    rqi = 1;
  } else if (lwrk >= (4 * ncol)) {
    rq = &wrk[0];
    gsumwrk = &wrk[3 * ncol];
    aq = Calloc(nrow, double);
    aqi = 1;
  } else if (lwrk >= ncol) {
    gsumwrk = wrk;
    aq = Calloc(nrow, double);
    aqi = 1;
    rq = Calloc(3 * ncol, double);
    rqi = 1;
  } else {
    /* WRK not provided -- allocate space for AQ and RQ,  */
    /* gsumwrk points to the last third of RQ             */
    aq = Calloc(nrow, double);
    aqi = 1;
    rq = Calloc(3 * ncol, double);
    rqi = 1;
    gsumwrk = Calloc(ncol, double);
    gsumwrki = 1;
  }
  memset(aq, 0, nrow * sizeof(double));
  memset(rq, 0, 2 * ncol * sizeof(double));
  memset(gsumwrk, 0, ncol * sizeof(double));
  /* go through each Ritz pair one at a time, compute Rayleigh  */
  /* quotient and the corresponding residual norm               */
  res = &rq[ncol];
  for (i = 0; i < ncol; i++) {
    op(&nrow, &i__1, &rvec[i * ldrvec], &nrow, aq, &nrow, lparam);
    /* Rayleigh quotient -- assuming rvec(:,i) has unit norm */
    rq[i] = trl_ddot(nrow, &rvec[i * ldrvec], c__1, aq, c__1);
    trl_g_sum(info->mpicom, 1, &rq[i], gsumwrk);
    d__1 = -rq[i]; /* indent separated =- into = - */
    trl_daxpy(nrow, d__1, &rvec[i * ldrvec], c__1, aq, c__1);
    res[i] = trl_ddot(nrow, aq, c__1, aq, c__1);
  }
  trl_g_sum(info->mpicom, ncol, res, gsumwrk);
  for (i = 0; i < ncol; i++) {
    res[i] = sqrt(res[i]);
  }
  /* compute the error estimate based on computed residual norms */
  /*  and the Ritz values                                        */
  err = &rq[2 * ncol];
  gapl = alpha[ncol - 1] - alpha[0];
  for (i = 0; i < ncol - 1; i++) {
    gapr = alpha[i + 1] - alpha[i];
    gapl = fmin2(gapl, gapr);
    if (res[i] >= gapl) {
      err[i] = res[i];
    } else {
      err[i] = res[i] * res[i] / gapl;
    }
    gapl = gapr;
  }
  if (res[ncol - 1] >= gapl) {
    err[ncol - 1] = res[ncol - 1];
  } else {
    err[ncol - 1] = res[ncol - 1] * res[ncol - 1] / gapl;
  }

  /* if writing to stdout, only PE 0 does it */
  if (info->log_fp == NULL)
    trl_reopen_logfile(info);

  if (info->log_fp != stdout || info->my_pe <= 0) {
    if (info->stat != 0) {
      *check = -4;
    }
    /* print out the information */
    fprintf(info->log_fp, "TRL_CHECK_RITZ: \n");
    fprintf(info->log_fp,
            "           Ritz value       res norm   res diff  est error  diff w rq  act. error\n");
    if (beta != NULL && eval != NULL) {
      for (i = 0; i < ncol; i++) {
        icheck = 0;
        fprintf(info->log_fp,
                "%21.14f    %11.3e%11.3e%11.3e%11.3e %11.3e%11.3e\n",
                alpha[i], res[i], beta[i] - res[i], err[i],
                rq[i] - alpha[i], eval[i] - alpha[i], eval[i]);
        /* check the accuracy of results.. */
        if (fabs(beta[i] - res[i]) > 0.00001) {
          *check = *check - 1;
          icheck++;
        }

        if (fabs(rq[i] - alpha[i]) > nrow * nrow * info->tol) {
          *check = *check - 1;
          icheck++;
        }

        if ((fabs(eval[i] - alpha[i]) >
             10 * nrow * nrow * info->tol) ||
            (fabs(eval[i] - alpha[i]) > 10 * err[i])) {
          *check = *check - 1;
          icheck++;
        }
      }

    } else if (beta != NULL) {
      for (i = 0; i < ncol; i++) {
        fprintf(info->log_fp, "%21.14f    %11.3e%11.3e%11.3e%11.3e\n",
                alpha[i], res[i], beta[i] - res[i], err[i],
                rq[i] - alpha[i]);
        /* check the accuracy of results.. */
        if (fabs(beta[i] - res[i]) > 0.00001) {
          *check = *check - 1;
          icheck++;
        }

        if (fabs(rq[i] - alpha[i]) > nrow * nrow * info->tol) {
          *check = *check - 1;
          icheck++;
        }
      }
    } else if (eval != NULL) {
      for (i = 0; i < ncol; i++) {
        fprintf(info->log_fp,
                "%21.14f     %11.3e           %11.3e%11.3e%11.3e%11.3e\n",
                alpha[i], res[i], err[i], rq[i] - alpha[i],
                eval[i] - alpha[i], eval[i]);
      }
    } else {
      for (i = 0; i < ncol; i++) {
        fprintf(info->log_fp, "%21.14f    %11.5e           %11.3e%11.3e\n",
                alpha[i], res[i], err[i], rq[i] - alpha[i]);
      }
    }
  }
  if (info->nec < info->ned)
    *check = 1;

  if (rqi > 0)
    Free(rq);

  if (aqi > 0)
    Free(aq);

  if (gsumwrki > 0)
    Free(gsumwrk);

  trl_close_logfile(info);
}

void
trl_rayleigh_quotients(trl_matprod op,
                       trl_info * info, int ncol, double *evec,
                       double *eres, double *base, void *lparam) {
  int c__1 = 1;
  int i__1 = 1;
  double d__1;
  int i, nrow;
  double wrk[4], *avec;

  nrow = info->nloc;
  if (ncol <= 0)
    return;
  if (base != NULL) {
    avec = base;
  } else {
    avec = Calloc(nrow, double);
  }
  memset(avec, 0, nrow * sizeof(double));
  if (info->verbose >= 0) {
    FILE *fp = info->log_fp;
    if (fp == NULL) {
      trl_reopen_logfile(info);
    }
    fp = info->log_fp;
    fprintf(fp,
            "TRLAN computing Rayleigh Quotients for %d Ritz pairs\n",
            ncol);
  }
  /* loop through each vector to normalize the vector, compute Rayleigh  */
  /* quotient and compute residual norm of the new Ritz pairs            */
  for (i = 0; i < ncol; i++) {
    wrk[0] =
      trl_ddot(nrow, &evec[i * nrow], c__1, &evec[i * nrow], c__1);
    op(&nrow, &i__1, &evec[i * nrow], &nrow, avec, &nrow, lparam);
    wrk[1] = trl_ddot(nrow, &evec[i * nrow], c__1, avec, c__1);
    trl_g_sum(info->mpicom, 2, wrk, &wrk[2]);
    info->matvec = info->matvec + 1;
    if (wrk[0] > 0.0) {
      eres[i] = wrk[1] / wrk[0];
      d__1 = -eres[i];
      trl_daxpy(nrow, d__1, &evec[i * nrow], c__1, avec, c__1);
      wrk[1] = trl_ddot(nrow, avec, c__1, avec, c__1);
      trl_g_sum(info->mpicom, 1, &wrk[1], &wrk[2]);
      wrk[0] = 1.0 / sqrt(wrk[0]);
      eres[ncol + i] = wrk[0] * sqrt(wrk[1]);
      d__1 = wrk[0];
      trl_dscal(nrow, d__1, &evec[i * nrow], c__1);
    } else {
      eres[i] = -DBL_MAX;
      eres[ncol + i] = -DBL_MAX;
    }
  }
  if (base == NULL)
    Free(avec);

  trl_close_logfile(info);
}
