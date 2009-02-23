/*
  ZTRLan routine (version 1.0)
  Lawrence Berkeley National Lab.
*/

/*
  All communication routines used by TRLAN are in this file. This file contains
  routines to be used on sequential or shared memory environments.  No actual
  data exchange is performed.
*/

#include <R.h>
#include <Rmath.h>

#include "trlan.h"
#include "trlan_i.h"
#include "trl_comm_i.h"

void trl_init_info(trl_info * info, int nrow, int mxlan, int lohi,
                   int ned, double tol, int restart, int maxmv,
                   int mpicom) {
  if (tol > 0) {
    info->tol = tol;
    if (info->tol <= DBL_MIN) {
      info->tol = DBL_EPSILON;
    } else if (info->tol > 1.0) {
      info->tol = fmin2(0.1, 1.0 / (info->tol));
    }
  } else {
    info->tol = sqrt(DBL_EPSILON);
  }
  if (restart > 0) {
    info->restart = restart;
  } else {
    info->restart = 0;
  }
  if (maxmv > 0) {
    info->maxmv = maxmv;
  } else {
    info->maxmv = imin2(imax2(info->ntot, 1000), 1000 * info->ned);
  }
  info->mpicom = -INT_MAX;

  info->maxlan = mxlan;
  if (mxlan <= ned) {
    info->maxlan = ned + imax2(ned, 6);
  }
  info->lohi = lohi;
  info->ned = ned;
  info->nloc = nrow;
  info->ntot = nrow;
  info->guess = 0;
  info->nec = 0;
  info->locked = info->nec;
  info->matvec = 0;
  info->nloop = 0;
  info->north = 0;
  info->nrand = 0;
  info->clk_rate = CLOCKS_PER_SEC;
#ifdef __64INT
  info->clk_max = 9223372036854775807LL;
#else
  info->clk_max = (clock_t) (pow(2.0, (8.0 * sizeof(clock_t) - 1.0))-1.0);
  if (info->clk_max < 0) {
    if (sizeof(clock_t) == 8) {
      info->clk_max = 9223372036854775807LL;
    } else {
      error("error initializing clock.");
    }
  }
#endif
  if ((double)(info->clk_max) <= 0)
    error( "??\n" );
  info->clk_tot = 0;
  info->clk_op = 0;
  info->clk_orth = 0;
  info->clk_res = 0;
  info->tick_t = 0;
  info->tick_o = 0;
  info->tick_h = 0;
  info->tick_r = 0;
  info->clk_in = 0;
  info->clk_out = 0;
  info->wrds_in = 0;
  info->wrds_out = 0;
  info->verbose = 0;
  info->stat = 0;
  info->anrm = 0;
  info->tmv = -1;
  info->trgt = -DBL_MAX;
  info->tres = -1.0;
  info->crat = 1.0;

  info->predicted_crate = 0.0;
  info->old_target = 0.0;
  info->target_id = 0;
  info->ref = 0.0;
  info->avgm = 0.0;
  info->k1 = 0;
  info->k2 = 0;
  info->k = 0;

  info->my_pe = 0;
  info->npes = 1;
  info->cpflag = 0;
  strcpy(info->oldcpf, "");

  // log file pointer
  info->log_io = 99;
  strcpy(info->log_file, "");

  // checkpoint file pointer
  info->cpio = 98;
  strcpy(info->cpfile, "");
}

void trl_g_sum(int mpicom, int nelm, double *x, double *y)
{
}

int trl_sync_flag(int mpicom, int inflag) {
  return inflag;
}

void trl_g_dot(int mpicom, int nrow, double *v1, int ld1, int m1,
               double *v2, int ld2, int m2, double *rr, double *wrk) {
  char trans = 'T';
  double one = 1.0, zero = 0.0;
  int c__1 = 1;
  int i, nd;

  nd = m1 + m2;
  // nothing to do if both m1 and m2 are zero
  if (nd <= 0)
    return;
  // make sure the array sizes are correct
  if (ld1 < nrow || ld2 < nrow) {
    error("trl_g_dot: incorrect array sizes");
  }
  if (m1 > 2) {
    trl_dgemv(&trans, nrow, m1, one, v1, ld1, rr, c__1, zero, wrk,
              c__1);
  } else if (m1 == 2) {
    wrk[0] = zero;
    wrk[1] = zero;
    for (i = 0; i < nrow; i++) {
      wrk[0] += v1[i] * rr[i];
      wrk[1] += v1[ld1 + i] * rr[i];
    }
  } else if (m1 == 1) {
    wrk[0] = trl_ddot(nrow, v1, c__1, rr, c__1);
  }
  if (m2 > 2) {
    trl_dgemv(&trans, nrow, m2, one, v2, ld2, rr, c__1, zero, &wrk[m1],
              c__1);
  } else if (m2 == 2) {
    wrk[m1] = zero;
    wrk[nd - 1] = zero;
    for (i = 0; i < nrow; i++) {
      wrk[m1]     += v2[i]        * rr[i];
      wrk[nd - 1] += v2[ld2 + i] * rr[i];
    }
  } else if (m2 == 1) {
    wrk[m1] = trl_ddot(nrow, v2, c__1, rr, c__1);
  }
}
