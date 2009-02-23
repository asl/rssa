/*
  ZTRLan routine (version 1.0)
  Lawrence Berkeley National Lab.
*/

#include <R.h>
#include <Rmath.h>

#include "dsort2_i.h"
#include "trlan.h"
#include "trlan_i.h"

#define nint(a) (ceil(a)-a > 0.5  ? (int)(a) : ceil(a))

/*
  TRLAN Low level utility routines
  This file contains a number of routines related to decisions of how
  many Ritz pairs to save when restarting the Thick-Restart Lanczos
  method. The subroutine trl_shuffle_eig is the main access point.
*/

static void trl_restart_max_gap_cost_ratio(int n, int tind, trl_info * info,
                                           double *lambda, double *res, int *kl,
                                           int *kr);
static void trl_restart_max_gap_cost_ratio_static(int n, int tind,
                                                  trl_info * info,
                                                  double *lambda, double *res,
                                                  int *kl, int *kr);
static void trl_restart_fixed(int nd, int mnd, int tind, double *lambda,
                              double *res, trl_info * info, int *kl, int *kr);
static void trl_restart_scan(int nd, double *res, trl_info * info, int kept,
                             int *kl, int *kr);
static void trl_restart_small_res(int nd, int mnd, int tind, double *lambda,
                                  double *res, trl_info * info, int *kl,
                                  int *kr);
static void trl_restart_max_gap_ratio(int nd, int tind, int kept, double *lambda,
                                      double *res, trl_info * info, int *kl,
                                      int *kr);
static void trl_restart_max_progress(int nd, int tind, int kept, double *lambda,
                                     double *res, trl_info * info, int *kl,
                                     int *kr);
static void trl_restart_max_reduction(int nd, int tind, int kept, double *lambda,
                                      double *res, trl_info * info, int *kl,
                                      int *kr);
static void trl_restart_search_range(int nd, double *lambda, double *res,
                                     trl_info * info, int ncl, int ncr,
                                     int *lohi, int tind, int *klm, int *krm);
static double trl_min_gap_ratio(trl_info * info, int nd, int tind,
                                double *res);

void trl_shuffle_eig(int nd, int mnd, double *lambda, double *res,
                      trl_info * info, int *kept, int locked) {
  int i, ncl, ncr, kl, kr, tind, minsep;
  double bnd;

  // very small basis -- save the half with the smallest residual norms
  if (nd <= 5) {
    dsort2(nd, res, lambda);
    if (nd > 3) {
      *kept = 2;
    } else if (nd > 0) {
      *kept = 1;
    } else {
      *kept = 0;
    }
    if (*kept > 1)
      dsort2(*kept, lambda, res);
    return;
  }

  /*
    preparation for normal case, first determine what are converged.
    ncl are the index (base zero) of res converged from the left
    ncr are the index (base zero) of res converged from the right */
  bnd = fmin2(info->tol, DBL_EPSILON * info->anrm);
  ncr = 0;
  ncl = nd - 1;
  i = nd - 1;

  /* determine how many has converged from the right */
  while (i >= 0) {
    if (res[i] <= bnd) {
      i--;
    } else {
      ncr = i + 1;
      i = -1;
    }
  }
  i = 0;

  /* determine how many has converged from the left */
  while (i < nd) {
    if (res[i] <= bnd) {
      i++;
    } else {
      ncl = i - 1;
      i = nd;
    }
  }
  kl = ncl;
  kr = ncr;
  if (ncr > ncl) {
    /* find the one that is closest to info->trgt */
    tind = (kl + kr) / 2;
    while (lambda[tind] != info->trgt && kr > kl) {
      if (lambda[tind] < info->trgt) {
        kl = tind + 1;
        tind = (kl + kr) / 2;
      } else if (lambda[tind] > info->trgt) {
        kr = tind - 1;
        tind = (kl + kr) / 2;
      } else {
        kl = tind;
        kr = tind;
      }
    }

    /* assign kl to the largest index of lambda that is smaller than
       info->trgt */
    if (lambda[tind] == info->trgt) {
      kl = tind - 1;
      while (kl >= 0 && lambda[kl] == info->trgt) {
        kl--;
      }
      /* assign kr to the smallest index of lambda that is greater than
         info->trgt */
      kr = tind + 1;
      while (kr < nd && lambda[kr] == info->trgt) {
        kr++;
      }
    } else {
      kl = tind - 1;
      kr = tind + 1;
    }
    /* initial assignment of kl and kr */
    if (info->lohi > 0) {
      /* large eigenvalues. */
      kr = kl;
      kl = imin2(ncl, imax2(0, nd - info->ned) - 1);
    } else if (info->lohi < 0) {
      /* small eigenvalues. */
      kl = kr;
      kr = imax2(ncr, imin2(nd - info->nec, info->ned + 1) - 1);
    } else if (ncr - tind > tind - ncl) {
      /* tind is closer to smallest lambda converged from left */
      kl = kr;
      kr = imax2(ncr, imin2(nd - info->nec, info->ned + 1) - 1);
    } else {
      /* tind is closer to largest lambda converged from right */
      kr = kl;
      kl = imin2(ncl, imax2(0, nd - info->ned) - 1);
    }
  } else {
    /* all have converged, keep all -- should not happen */
    *kept = nd;
    return;
  }

  /* We keep at least one from each end. */
  #if 0
  if (kl < kr) {
    kl ++;
  }
  if (kr > kl) {
    kr --;
  }
  #endif

  /* We also save Ritz values with the smallest residual */
  #if 0
  if (kl < kr) {
    if (res[kl+1] < res[kr-1]) {
      kl ++;
    } else {
      kr --;
    }
  }
  #endif

  /* normal cases, call subroutines to complete the tasks
     [1 .. kl] and [kr .. nd] are saved for later the initial values of kl and
     kr are simply ncl and ncr they are further analyzed according to the
     restarting strategy requested
  */
  switch (info->restart) {
   case 1:
    /* fixed number beyond the currently converged ones */
    trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
    break;
   case 2:
    /* add the ones with smallest residual norms */
    trl_restart_small_res(nd, mnd, tind, lambda, res, info, &kl, &kr);
    break;
   case 3:
    if (info->nloop > 0) {
      /* maximize the gap ratio */
      trl_restart_max_gap_ratio(nd, tind, *kept, lambda, res, info,
                                 &kl, &kr);
    } else {
      /* this is the first restart -- use trl_restart_fixed instead */
      trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
    }
    break;
   case 4:
    if (info->nloop > 0) {
      /* maximize [gap-ratio * (m-k)] */
      trl_restart_max_progress(nd, tind, *kept, lambda, res, info,
                                &kl, &kr);
    } else {
      /* this is the first restart -- use trl_restart_fixed instead */
      trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
    }
    break;
   case 5:
    if (info->nloop > 0) {
      /* maximize [sqrt(gap tatio) * (m-k)] */
      trl_restart_max_reduction(nd, tind, *kept, lambda, res, info,
                                 &kl, &kr);
    } else {
      /* this is the first restart -- use trl_restart_fixed instead */
      trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
    }
    break;
   case 6:
    /* progressively vary the thickness */
    trl_restart_scan(nd, res, info, *kept, &kl, &kr);
    break;
   case 7:
    if (info->nloop > 0) {
      trl_restart_max_gap_cost_ratio(nd, tind, info, lambda, res,
                                      &kl, &kr);
    } else {
      trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
    }
    break;
   case 8:
    if (info->nloop > 0) {
      trl_restart_max_gap_cost_ratio_static(nd, tind, info, lambda,
                                             res, &kl, &kr);
    } else {
      trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
    }
    break;
   default:
    if (info->restart <= -info->ned) {
      if (info->lohi >= 0) {
        kl = -1;
        kr = imax2(2, nd + info->restart);
      } else if (info->lohi < 0) {
        kl = imin2(-info->restart, nd - 3) - 1;
        kr = nd;
      } else {
        kl = imin2(nd - 3, -info->restart) / 2 - 1;
        kr = nd - kl - 2;
      }
    } else {
      trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
    }
    break;
  }

  /* make sure kr > kl+minsep */
  minsep = imax2(3, imax2(nd / 6, nd - 6 * info->ned));
  if (kr <= kl + minsep || kl + nd - kr + minsep > mnd) {
    if (ncl < kl && kl < kr && kr < ncr) {
      kl--;
      kr++;
    } else if (info->lohi > 0) {
      kr = imax2(minsep, imin2(nd / 3, ncr)) - 1;
      kl = -1;
    } else if (info->lohi < 0) {
      kl = imin2(nd - minsep, imax2((nd + nd) / 3, ncl + 2)) - 1;
      kr = nd;
    } else {
      kl = (nd - minsep) / 2 - 2;
      kr = (nd - minsep + 1) / 2;
    }
  }

  /* copy the (kr:nd) elements to (kl+1:kl+nd-kr+2) */
  /* kr is temporarily decreased by 1 to make indices easier to compute */
  kr--;
  for (i = 1; i < nd - kr; i++) {
    lambda[kl + i] = lambda[kr + i];
    res[kl + i] = res[kr + i];
  }
  *kept = kl + imax2(1, nd - kr);
  if (info->lohi < -1) {
    dsort2(*kept, lambda, res);
  }
  return;
}

/* Internal functio used in trl_restart_fixed_ */
static inline double gap_ratio(int i, int j, int tind, double *lambda) {
  return (lambda[i] - lambda[tind]) / (lambda[j] - lambda[tind]);
}

//
// Purpose
// =======
// Determine the size of the maximum Lanczos basis size and the number of Ritz
// pairs to keep. The decision is base on the following criteria
//   argmax( m, k ) exp( - gamma * (m-k) )/ (m-k)(m-k-1)
// where gamma = (lambda(kl)-lambda(1))/(lambda(kr)-lambda(1)
// the numerator of the criteria approximates the expected improvement in the
// residual norm of the Ritz pairs at the next restart, and the denominator
// approximates the cost of the reorthogonalization till the next restart.
// The next basis size m is chosen between info-rfact * k and maxlan, where
// k is the number of Ritz vectors kept.
//
// Arguments
// =========
// n         (input/output) INTEGER
//            On entry, specifies the number of Lanczos basis considered.
//            On exit, specifies the number of maximum Lanczos basis size
//            for the next iterations.
//
// tind      (input) INTEGER
//            On entry, specifies the index of the next target Ritz value.
//
//
// lambda    (input) DOUBLE PRECISION ARRAY
//            On entry, contains the computed Ritz values.
//
// res       (input) DOUBLE PRECISION ARRAY
//            On entry, contains the residueal norm.
//
// kl        (input/output) Pointer to INTEGER
//            On entry, specifies the largest Ritz value from left, that has
//            been converged.
//            On exit, specifies the largest Ritz value from left, that are kept.
//
// kr        (input/output) Pointer to INTEGER
//            On entry, specifies the smallest Ritz value from right, that has
//            been converged.
//            On exit, specifies the smallest Ritz value from left, that are kept.
//
////
void trl_restart_max_gap_cost_ratio(int n, int tind, trl_info * info,
                                     double *lambda, double *res, int *kl,
                                     int *kr) {
  int i, j, k, nd, l1, l2, mn, mn1, mn2, t, k1, k2,
    min_k1, min_k2, min_m, min_l;
  double gamma0, gamma, val, min_val, min_gamma, min_gamma0, min_ratio,
    tmp, tmp2, tmp3, up, dw, def1, def2, z1, z2;

  nd = info->ned;
  mn = info->maxlan;
  trl_restart_search_range(n, lambda, res, info, (*kl), (*kr),
                            &(info->lohi), tind, &k1, &k2);
  // ** Static approach to decide the minimum gap ratio **
  // minimum gap beteween kl and kr for ok-conditioned systems,
  // i.e, diag(i) and diag(i^2)
  //t = nint(4.0*abs(k1-k2)/5.0);
  //t = nint((info->mgap)*abs(k1-k2));
  //t = nint(2.0*abs(k1-k2)/5.0);
  /*
    if( t > n-nd ) t=n-nd;
    if( t < 2 ) t = 2;
    if( info->lohi == -1 ) {
    mcnv0 = info->klan - (*kr);
    } else if( info->lohi == 1 ) {
    mcnv0 = (*kl)+1;
    }
    if( t > 2 && t+kept > nd && info->crat > 0.0 ) {
    min_val = trl_min_gap_ratio(info, nd, tind, res);
    if( min_val > info->crat ) t = max(2, nd-kept-1);
    }
  */
  //
  // ** Dynamic approach to decide the minimum gap ratio **
  def1 = 1.0;
  def2 = 0.7;
  up = 1.0;
  dw = 0.7;
  info->avgm =
    (info->avgm * ((double) info->nloop - 1) +
     (double) (info->klan - info->k + 1)) / ((double) info->nloop);
  z1 = (info->ptres) / (info->tol * info->anrm);
  z2 = 1.0 / info->cfac;
  if (z2 < 1.0) {
    tmp = def2;
    t = abs(k2 - k1) * def2;
  } else if (z1 < 1.0) {
    tmp = def1;
    t = abs(k2 - k1) * def1;
  } else {
    tmp3 =
      log(z1 +
          sqrt(z1 - 1.0) * sqrt(z1 + 1.0)) / (2.0 * (info->avgm));
    tmp2 =
      log(z2 + sqrt(z2 - 1.0) * sqrt(z2 + 1.0)) / (info->klan -
                                                   info->k + 1);
    tmp = tmp2 / tmp3;
    tmp = pow(tmp, 2.0);
    tmp = atan(tmp) * M_2_PI;
    tmp = dw + (up - dw) * tmp;

    t = abs(k2 - k1) * tmp;
  }
  //
  // ** Static approach to decided the minimum gap ratio **
  //t = max (min (n - info->ned, nint ((k2 - k1) * info->rfact)), 2);
  //
  mn1 = nint(2.0 * (info->klan) / 5.0);
  min_val = 0.0;
  min_k1 = k1;
  min_k2 = k2;
  min_m = info->klan;
  //printf( "k1=%d k2=%d\n",k1,k2 );
  for (i = k1; i <= k2 - t; i++) {
    // considering how many Ritz vectors to keep from left.
    for (j = i + t; j <= k2; j++) {
      // considering how many Ritz vectors to keep from right.
      //
      /* note the two option, i.e., square-rooted, are hard-coded here !!! */
      gamma0 =
        sqrt((lambda[i + 1] - lambda[tind]) / (lambda[j - 1] -
                                               lambda[i + 1]));
      /* original object function */
      //gamma0 = sqrt((lambda[i+1]-lambda[tind])/(lambda[j-1]-lambda[tind]));
      /* Without square */
      //gamma0 = (lambda[i+1]-lambda[tind])/(lambda[j-1]-lambda[tind]);
      l1 = (i + 1) + (n - j);
      l2 = l1 + info->locked;

      //mn2 = mn1;
      //if( mn1 < (info->rfact*l2) ) mn2 = info->rfact*l2;
      //if( mn1 < (1.5*l2) ) mn2 = 1.5*l2;
      //if( mn1 < (l2+1) ) mn2 = l2+1;
      mn2 = l2 + 1;
      if (mn2 > mn)
        mn2 = mn;
      if (mn2 < nd)
        mn2 = nd;

      for (k = mn2; k <= mn; k++) {
        // considering k for the next Lanczos sizes
        gamma = ((k - l2) * gamma0);
        //* Approximate cost function *//
        //val = ( k * k - l2 * l2 + 2 * k * l1 ) / ( gamma );
        //* More exact cost function *//
        val = ((k - l2) * (k + l2 - 1) + k * l1) / (gamma);
        //* Including cosh function *//
        //val = ( (k-l2) * (k+l2-1) + k * l1 ) / ( cosh( 2.0 * gamma ) );
        //printf( "%d %e %e\n",( (k-l2) * (k+l2-1) + k * l1 ), gamma, cosh(2.0*gamma) );

        if (min_val == 0.0 || val < min_val) {
          min_k1 = i;
          min_k2 = j;
          min_m = k;
          min_val = val;
          min_gamma = gamma;
          min_gamma0 = gamma0;
          min_ratio =
            (lambda[i + 1] - lambda[tind]) / (lambda[j - 1] -
                                              lambda[tind]);
          min_l = l2;
        }
      }
    }
  }
  *kl = min_k1;
  *kr = min_k2;
  info->klan = min_m;
  info->k1 = min_k1;
  info->k2 = min_k2;
  info->k = min_l;
  //printf( "%e %d %d %d\n",tmp,min_k1,min_k2,min_m );
  info->predicted_crate = exp(-min_gamma);
  info->mgamma = min_gamma;
  info->gamma0 = min_gamma0;
  info->old_target = lambda[tind];
  info->target_id = tind;
  info->old_locked = info->locked;
}

//
// Purpose
// =======
// Determine the size of the maximum Lanczos basis size and the number of Ritz
// pairs to keep. The decision is base on the following criteria
//   argmax( m, k ) exp( - gamma * (m-k) )/ (m-k)(m-k-1)
// where gamma = (lambda(kl)-lambda(1))/(lambda(kr)-lambda(1)
// the numerator of the criteria approximates the expected improvement in the
// residual norm of the Ritz pairs at the next restart, and the denominator
// approximates the cost of the reorthogonalization till the next restart.
// The next basis size m is fixed at info-rfact * k and maxlan, where k is the number
// of Ritz vectors kept.
//
// Arguments
// =========
// n         (input/output) INTEGER
//            On entry, specifies the number of Lanczos basis considered.
//            On exit, specifies the number of maximum Lanczos basis size
//            for the next iterations.
//
// tind      (input) INTEGER
//            On entry, specifies the index of the next target Ritz value.
//
//
// lambda    (input) DOUBLE PRECISION ARRAY
//            On entry, contains the computed Ritz values.
//
// res       (input) DOUBLE PRECISION ARRAY
//            On entry, contains the residueal norm.
//
// kl        (input/output) Pointer to INTEGER
//            On entry, specifies the largest Ritz value from left, that has
//            been converged.
//            On exit, specifies the largest Ritz value from left, that are kept.
//
// kr        (input/output) Pointer to INTEGER
//            On entry, specifies the smallest Ritz value from right, that has
//            been converged.
//            On exit, specifies the smallest Ritz value from left, that are kept.
//
////
void trl_restart_max_gap_cost_ratio_static(int n, int tind,
                                            trl_info * info,
                                            double *lambda, double *res,
                                            int *kl, int *kr)
{
  int i, j, nd, l1, l2, mn, mn1, mn2, t, k1, k2,
    min_k1, min_k2, min_m, min_l;
  double gamma0, gamma, val, min_val, min_gamma, min_gamma0, min_ratio,
    tmp2, tmp3, tmp;

  nd = info->ned;
  mn = info->maxlan;
  trl_restart_search_range(n, lambda, res, info, (*kl), (*kr),
                            &(info->lohi), tind, &k1, &k2);
  //
  // ** Static approach to decide the minimum gap ratio **
  //t = nint(4.0*abs(k1-k2)/5.0);
  /*
    t = nint((info->mgap)*abs(k1-k2));
    if( t > n-nd ) t=n-nd;
    if( t < 2 ) t = 2;

    if( t > 2 && t+kept > nd && info->crat > 0.0 ) {
    min_val = trl_min_gap_ratio(info, nd, tind, res);
    if( min_val > info->crat ) t = max(2, nd-kept-1);
    }
  */
  //
  // ** Dynamic approach to decide the minimum gap ratio **
  if (info->crat > 0.0) {
    tmp2 = -log(info->crat);
    tmp3 = -(tmp2 * (info->maxlan)) / log(info->tol * info->anrm);
    if (tmp3 >= 0.5 || info->klan < info->maxlan) {
      tmp = 0.8;
      t = abs(k2 - k1) * 0.8;
    } else {
      tmp = trl_min_gap_ratio(info, nd, tind, res);
      tmp = tmp2 / tmp;
      tmp = pow(tmp, 0.25);
      tmp = atan(tmp) * (2.0 * M_2_PI);
      t = abs(k2 - k1) * tmp;
    }
  } else {
    tmp = 0.8;
    t = abs(k2 - k1) * 0.8;
  }
  mn1 = nint(2.0 * (info->klan) / 5.0);
  min_val = 0.0;
  min_k1 = k1;
  min_k2 = k2;
  min_m = info->klan;
  for (i = k1; i <= k2 - t; i++) {
    // considering how many Ritz vectors to keep from left.
    for (j = i + t; j <= k2; j++) {
      // considering how many Ritz vectors to keep from right.
      gamma0 =
        sqrt((lambda[i + 1] - lambda[tind]) / (lambda[j - 1] -
                                               lambda[i + 1]));
      //gamma0 = sqrt((lambda[i+1]-lambda[tind])/(lambda[j-1]-lambda[tind]));
      //gamma0 = (lambda[i+1]-lambda[tind])/(lambda[j-1]-lambda[tind]);
      l1 = (i + 1) + (n - j);
      l2 = l1 + info->locked;
      mn2 = mn1;
      if (mn1 < ((info->rfact) * l2))
        mn2 = (info->rfact) * l2;
      if (mn2 > mn)
        mn2 = mn;
      if (mn2 < nd)
        mn2 = nd;
      gamma = (mn2 - l2 + 1) * gamma0;
      val = (mn2 * mn2 - l2 * l2 + 2 * mn2 * l1) / (gamma);
      if (min_val == 0.0 || val < min_val) {
        min_k1 = i;
        min_k2 = j;
        min_m = mn2;
        min_val = val;
        min_gamma = gamma;
        min_gamma0 = gamma0;
        min_ratio =
          (lambda[i + 1] - lambda[tind]) / (lambda[j - 1] -
                                            lambda[tind]);
        min_l = l2;
      }
    }
  }
  *kl = min_k1;
  *kr = min_k2;
  info->klan = min_m;
}

//
// Purpose
// =======
// Save fixed number of extra Ritz pairs.
//
// Arguments
// =========
// nd         (input) INTEGER
//             On entry, specifies the number of Lanczos basis considered.
//
// mnd        (input) INTEGER
//             On entry, specifies the number of maximum Lanczos basis.
//             (lanczos basis size used)-(Ritz values locked).
//
// tind       (input) INTEGER
//             On entry, specifies the index of lambda, that is closest to target.
//
// lambda     (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the computed Ritz values.
//
// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the residual of the computed Ritz values.
//
// info       (input) POINTER to TRLANINFO structure
//             On entry, points to the current TRLANINFO structure.
//
// kl         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from left.
//
// kr         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from right.
//
////
void trl_restart_fixed(int nd, int mnd, int tind, double *lambda,
                        double *res, trl_info * info, int *kl, int *kr) {
  int extra, i, kl0, kr0, minsep, kli, kri;
  double gmin;

  // the number of extra Ritz pairs to be saved
  kli = *kl;
  kri = *kr;
  kl0 = kli;
  kr0 = kri;
  extra =
    nint((mnd - info->nec) * (0.4 + 0.1 * info->ned / (double) mnd));
  if (extra > info->ned + info->ned && extra > 5) {
    gmin = ((double) mnd) / ((double) info->ned);
    extra =
      nint((extra + (log(gmin) * info->ned) * gmin) / (1.0 + gmin));
  }
  minsep = imax2(3, imax2(nd / 5, nd - 4 * info->ned));
  gmin = trl_min_gap_ratio(info, nd, tind, res);
  if (info->lohi > 0) {
    kri = imin2(tind - 1, kri) - extra;
    kli = -1;
  } else if (info->lohi < 0) {
    kli = imax2(tind + 1, kli) + extra;
    kri = nd;
  } else {
    extra++;
    kli = kli + extra / 2;
    kri = kri - extra / 2;
    i = 0;
    while (kli > kl0 && kri < kr0 && i == 0) {
      if (10.0 * res[kli] < res[kri]) {
        // lambda converged much more from left, so shift to right
        if (res[kli + 1] < res[kri + 1]) {
          kli++;
          kri++;
        } else {
          i = -1;
        }
      } else if (10.0 * res[kri] < res[kli]) {
        // lambda converged much more from right, so shift to left
        if (res[kri - 1] < res[kli - 1]) {
          kri--;
          kli--;
        } else {
          i = -1;
        }
      } else {
        i = -1;
      }
    }
  }
  // adjust kl and kr until the minimal gap ratio is satisfied
  while (kli + minsep < kri &&
         gap_ratio(imax2(0, kli), imin2(kri, nd - 1), tind,
                    lambda) < gmin) {
    if (info->lohi < 0) {
      kli++;
    } else if (info->lohi > 0) {
      kri--;
    } else if (res[kli] < res[kri]) {
      kli++;
    } else {
      kri++;
    }
  }
  // make sure nearly degenerate/duplicated Ritz pairs are included
  // lambda(kl)+r(kl) > lambda(j) and
  //                lambda(kl)-r(kl) > lambda(j)-r(j)
  if (info->lohi > 0) {
    i = kri - 1;
    // (lambda[i],lambda[i]+res[i]) is
    //   in (lambda[kri]-res[kri],lambda[kri]+res[kri]).
    while (i > kli + minsep && lambda[kri] - res[kri] < lambda[i] &&
           lambda[kri] + res[kri] < lambda[i] + res[i]) {
      i--;
    }
    kri = i + 1;
  } else {
    kl0 = kli;
    i = kli + 1;
    while (i < kri - minsep && lambda[kli] + res[kli] > lambda[i] &&
           lambda[kli] - res[kli] > lambda[i] - res[i]) {
      i++;
    }
    kli = i - 1;
  }
  *kl = kli;
  *kr = kri;
}

// Purpose
// =======
// This subroutine determines the number Ritz vectors to be saved by
// adding some quantity to the current thickness.  If the thickness is
// larger than nd-2, it is reset to something smaller.
// The thickness is varied progressively to scan all possible values.
//
// Arguments
// =========
// nd         (input) INTEGER
//             On entry, specifies the number of Lanczos basis considered.
//
// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the residual of the computed Ritz values.
//
// info       (input) POINTER to TRLANINFO structure
//             On entry, points to the current TRLANINFO structure.
//
// kept       (input) INTEGER
//             On entry, specifies the number of Ritz values kept.
//
// kl         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from left.
//
// kr         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from right.
//
////
void trl_restart_scan(int nd, double *res, trl_info * info, int kept,
                       int *kl, int *kr)
{
  int extra, i, kl0, kr0, kli, kri;

  kli = *kl;
  kri = *kr;
  // three cases have to be dealt with separately
  if (info->lohi < 0) {
    // smallest eigenvalues
    kri = nd;
    kli = kept + imin2(imax2(info->nec, 1), (nd - kept) / 2) - 1;
    if (kli <= 0) {
      if (nd > 6) {
        kli = nd / 2 - 1;
      } else if (nd > 2) {
        kli = 1;
      }
    } else if (kli + 3 > nd) {
      kli =
        info->nec + imin2(info->ned,
                          imin2(10, (nd - info->ned) / 2)) - 1;
    }
  } else if (info->lohi > 0) {
    kli = -1;
    kri = kept + imin2(imax2(info->nec, 1), (nd - kept) / 2) - 1;
    if (kri <= 0) {
      if (nd > 6) {
        kri = nd / 2 - 1;
      } else if (nd > 2) {
        kri = 1;
      }
    } else if (kri + 4 > nd) {
      kri =
        info->nec + imin2(info->ned, imin2(10, (nd - info->ned) / 2));
    }
    kri = nd - kri - 1;
  } else {
    kl0 = kli;
    kr0 = kri;
    extra = kept + imin2(info->nec, (nd - kept) / 2) + 1;
    if (extra <= 0) {
      if (nd > 6) {
        extra = nd / 2;
      } else if (nd > 2) {
        extra = 2;
      }
    } else if (extra + 3 > nd) {
      extra =
        info->nec + imin2(info->ned, imin2(10, (nd - info->ned) / 2));
    }
    kli = imax2(kli, (extra / 2) - 1);
    kri = imin2(kri, (nd - extra / 2) - 1);
    i = 0;
    while (kli > kl0 && kri < kr0 && i == 0) {
      if (10.0 * res[kli] < res[kri]) {
        if (res[kli + 1] < res[kri + 1]) {
          kli++;
          kri++;
        } else {
          i = -1;
        }
      } else if (10.0 * res[kri] < res[kli]) {
        if (res[kri - 1] < res[kli - 1]) {
          kri--;
          kli--;
        } else {
          i = -1;
        }
      } else {
        i = -1;
      }
    }
  }
  *kl = kli;
  *kr = kri;
}

/*
  Purpose
  =======
  Returns maximum value in a

  Arguments
  =========
  n   (input) INTEGER
  On entry, specifies the size of the array.

  a   (input) DOUBLE PREICISION ARRAY of LENGTH n
  On entry, contains the values to search for the maximum value.
*/
static inline double maxval(int n, double *a)
{
  int i;
  double val;
  if (n <= 0) {
    return 0.0;
  }
  val = a[0];
  for (i = 1; i < n; i++) {
    val = fmax2(val, a[i]);
  }
  return val;
}

//
// Arguments
// =========
// nd         (input) INTEGER
//             On entry, specifies the number of Lanczos basis considered.
//
// mnd        (input) INTEGER
//             On entry, specifies the number of maximum Lanczos basis.
//             (lanczos basis size used)-(Ritz values locked).
//
// tind       (input) INTEGER
//             On entry, specifies the index of lambda, that is closest to target.
//
// lambda     (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the computed Ritz values.
//
// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the residual of the computed Ritz values.
//
// info       (input) POINTER to TRLANINFO structure
//             On entry, points to the current TRLANINFO structure.
//
// kl         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from left.
//
// kr         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from right.
//
//
//
////
void trl_restart_small_res(int nd, int mnd, int tind, double *lambda,
                            double *res, trl_info * info, int *kl,
                            int *kr)
{
  int extra, i, j, ii, kl0, kr0, kli, kri, minsep, done;
  double acpt, resmax, gmin;

  // the number of extra Ritz pairs to be saved
  minsep = imax2(3, imax2(nd / 5, nd - 4 * info->ned));
  extra =
    nint((mnd - info->nec) * (0.4 + 0.1 * info->ned / ((double) mnd)));
  if (extra > info->ned + info->ned && extra > 5) {
    gmin = ((double) mnd) / ((double) info->ned);
    extra =
      nint((extra + (log(gmin) * info->ned) * gmin) / (1.0 + gmin));
  }
  kli = *kl;
  kri = *kr;
  kl0 = kli;
  kr0 = kri;
  resmax = maxval(nd, res);
  acpt = resmax / res[tind];
//
// determine the number of Ritz pairs that has to be saved
  if (info->lohi > 0) {
    if (acpt < 0.999 && acpt >= 0.0) {
      ii = tind - 1;
      acpt = fmax2(sqrt(acpt) * res[tind], res[ii] + res[ii]);
      acpt = fmin2(acpt, resmax);
      kri = ii - 1;
      while (res[kri] < acpt && kri > kli + 3) {
        kri--;
      }
    } else {
      kri = kr0 - extra;
    }
    kri = imax2(2, kri);
    kli = imin2(kli, kri - 2);
  } else if (info->lohi < 0) {
    if (acpt < 0.999 && acpt >= 0.0) {
      ii = tind + 1;
      acpt = fmax2(sqrt(acpt) * res[tind], res[ii] + res[ii]);
      acpt = fmin2(acpt, resmax);
      kli = ii + 1;
      while (res[kli] < acpt && kli < kri - 3) {
        kli++;
      }
    } else {
      kli = kl0 + extra;
    }
    kli = imin2(nd - 4, kli);
    kri = imax2(kri, kli + 2);
  } else {
    // save whoever has smaller residual norms
    i = kli + 1;
    j = kri - 1;
    for (ii = 1; ii <= extra; ii++) {
      if (res[i] < res[j]) {
        kli = i;
        i++;
      } else if (res[i] > res[j]) {
        kri = j;
        j--;
      } else if (i < nd - j) {
        kli = i;
        i++;
      } else {
        kri = j;
        j--;
      }
    }
  }
  // adjust kl and kr until the minimal gap ratio is satisfied
  kl0 = kli;
  kr0 = kri;
  gmin = trl_min_gap_ratio(info, nd, tind, res);
  done = 0;
  while (kli + minsep < kri &&
         gap_ratio(imax2(0, kli), imin2(nd - 1, kri), tind, lambda) < gmin
         && (done == 0)) {
    if (info->lohi < 0) {
      kli++;
    } else if (info->lohi > 0) {
      kri--;
    } else if (res[kli] < res[kri]) {
      kli++;
    } else if (kri < nd - 1) {
      kri++;
    } else {
      done = 1;
    }
  }
  // make sure nearly degenerate Ritz pairs are included
  // lambda(kl)+r(kl) > lambda(j) and
  //                lambda(kl)-r(kl) > lambda(j)-r(j)
  if (info->lohi > 0) {
    i = kr0 - 1;
    while (i > kli + minsep && lambda[kri] - res[kri] < lambda[i] &&
           lambda[kri] + res[kri] < lambda[i] + res[i]) {
      i--;
    }
    kri = imin2(kri, i + 1);
  } else {
    i = kl0 + 1;
    while (i < kri - minsep && lambda[kli] + res[kli] > lambda[i] &&
           lambda[kli] - res[kli] > lambda[i] - res[i]) {
      i++;
    }
    kli = imax2(kli, i - 1);
  }
  *kl = kli;
  *kr = kri;
}

/* statement function for computing gap ratio */
#define gratio(i,j)                                     \
  (lambda[(j)] == info->trgt ?                          \
   DBL_MAX :                                            \
   (lambda[(i)]-info->trgt)/(lambda[(j)]-info->trgt) )

//
// Purpose
// =======
// Search throught all pairs of (kl, kr) for the one with the maximum
// gap ratio for the next Ritz pair(target).
// This is an optimized version of the original version.  It only search
// through nd number once. (Only single loop!)
//
// Arguments
// =========
// nd         (input) INTEGER
//             On entry, specifies the number of Lanczos basis considered.
//
// mnd        (input) INTEGER
//             On entry, specifies the number of maximum Lanczos basis.
//             (lanczos basis size used)-(Ritz values locked).
//
// tind       (input) INTEGER
//             On entry, specifies the index of lambda, that is closest to target.
//
// kept       (input) INTEGER
//             On entry, specifies the number of Ritz values kept.
//
// lambda     (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the computed Ritz values.
//
// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the residual of the computed Ritz values.
//
// info       (input) POINTER to TRLANINFO structure
//             On entry, points to the current TRLANINFO structure.
//
// kl         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from left.
//
// kr         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from right.
//
//
////
void trl_restart_max_gap_ratio(int nd, int tind, int kept, double *lambda,
                                double *res, trl_info * info, int *kl,
                                int *kr)
{
  int i, j, lohi, klm, krm, kli, kri, igap;
  double bnd, tmp;

  // determine the search range
  kli = *kl;
  kri = *kr;
  trl_restart_search_range(nd, lambda, res, info, kli, kri, &lohi, tind,
                            &klm, &krm);
  kli = klm;
  kri = krm;
  igap = imax2(imin2(nd - info->ned, nint((krm - klm) * 0.4)), 2);
  if (igap > 2 && info->matvec > info->maxlan) {
    if (info->clk_op + info->tick_o >
        10.0 * (info->clk_orth + info->tick_h + info->clk_res +
                info->tick_r)) {
      igap = imax2(2, nd - kept - 1);
    } else {
      bnd = trl_min_gap_ratio(info, nd, tind, res);
      if (info->crat < bnd)
        igap = imax2(2, nd - kept - 1);
    }
  }
  if (kli + igap > kri) {
    *kl = kli;
    *kr = kri;
    return;
  }
  // two cases depending on lohi
  if (lohi > 0) {
    // target is at the high end of spectrum
    bnd = gratio(kri, kli);
    for (i = klm; i <= krm - igap; i++) {
      j = i + igap;
      tmp = gratio(j, i);
      if (tmp > bnd) {
        kli = i;
        kri = j;
        bnd = tmp;
      }
    }
  } else {
    bnd = gratio(kli, kri);
    for (i = klm; i <= krm - igap; i++) {
      j = i + igap;
      tmp = gratio(i, j);
      if (tmp > bnd) {
        kli = i;
        kri = j;
        bnd = tmp;
      }
    }
  }
  *kl = kli;
  *kr = kri;
}

/* merit measure the factor of residual norm reduction */
#define merit(i,j)                                                  \
  ((lambda[(i)]-info->trgt) * abs(j-i) / (lambda[(j)]-info->trgt))

//
// Purpose
// =======
// Search for a pair (kl, kr) such that the reduction in residual norm
// of the target (info%trgt) will be the largest before next restart
// The merit function is [gap-ratio * (m-k)]
//
// Arguments
// =========
// nd         (input) INTEGER
//             On entry, specifies the number of Lanczos basis considered.
//
// tind       (input) INTEGER
//             On entry, specifies the index of lambda, that is closest to target.
//
// kept       (input) INTEGER
//             On entry, specifies the number of Ritz values kept.
//
// lambda     (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the computed Ritz values.
//
// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the residual of the computed Ritz values.
//
// info       (input) POINTER to TRLANINFO structure
//             On entry, points to the current TRLANINFO structure.
//
// kl         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from left.
//
// kr         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from right.
//
////
void trl_restart_max_progress(int nd, int tind, int kept, double *lambda,
                               double *res, trl_info * info, int *kl,
                               int *kr) {
  int i, j, lohi, kli, kri, klm, krm, igap;
  double tmp, ss;

  /* determine the search range */
  trl_restart_search_range(nd, lambda, res, info, *kl, *kr, &lohi, tind,
                            &klm, &krm);

  /* perform the brute-force search */
  kli = klm;
  kri = krm;
  igap = imax2(imin2(nd - info->ned, nint((kri - kli) * 0.4)), 2);
  if (igap > 2 && igap + kept > nd && info->crat > 0.0) {
    ss = trl_min_gap_ratio(info, nd, tind, res);
    if (ss > info->crat)
      igap = imax2(2, nd - kept - 1);
  }
  if (lohi > 0) {
    ss = merit(kri, kli);
    for (i = klm; i <= krm - igap; i++) {
      for (j = i + igap; j <= krm; j++) {
        tmp = merit(j, i);
        if (tmp > ss) {
          ss = tmp;
          kli = i;
          kri = j;
        }
      }
    }
  } else {
    ss = merit(kli, kri);
    for (i = klm; i <= krm - igap; i++) {
      for (j = i + igap; j <= krm; j++) {
        tmp = merit(i, j);
        if (tmp > ss) {
          ss = tmp;
          kli = i;
          kri = j;
        }
      }
    }
  }
  *kl = kli;
  *kr = kri;
}

/* merit measure the factor of residual norm reduction */
/*#define merit_maxred(i,j) ( sqrt( (lambda[(i)]-info->trgt)/(lambda[(j)]-info->trgt)) * abs(j-i) );*/
#define merit_maxred(i,j)                                               \
  (sqrt((lambda[(i)]-info->trgt)/(lambda[(j)]-lambda[(i)])) * abs(j-i));

//
// Purpose
// =======
// Search for a pair (kl, kr) such that the reduction in residual norm
// of the target (info%trgt) will be the largest before next restart
// the merit function is [sqrt(gap ratio) * (m-k)]
//
// Arguments
// =========
// nd         (input) INTEGER
//             On entry, specifies the number of Lanczos basis considered.
//
// tind       (input) INTEGER
//             On entry, specifies the index of lambda, that is closest to target.
//
// kept       (input) INTEGER
//             On entry, specifies the number of Ritz values kept.
//
// lambda     (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the computed Ritz values.
//
// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
//             On entry, contains the residual of the computed Ritz values.
//
// info       (input) POINTER to TRLANINFO structure
//             On entry, points to the current TRLANINFO structure.
//
// kl         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from left.
//
// kr         (input/output) POINTER to INTEGER
//             On exit, points to the starting index of lambda to keep from right.
////
void trl_restart_max_reduction(int nd, int tind, int kept, double *lambda,
                               double *res, trl_info * info, int *kl,
                               int *kr) {
  int i, j, lohi, kli, kri, klm, krm, t, igap;
  double tmp, tmp2, tmp3, ss, up, dw, def1, def2, z1, z2;

  // determine the search range
  kli = *kl;
  kri = *kr;
  trl_restart_search_range(nd, lambda, res, info, kli, kri, &lohi, tind,
                            &klm, &krm);
  // perform the brute-force search
  kli = klm;
  kri = krm;
  //
  // ** Static approach to decide the minimum gap ratio **
  //igap = max( min(nd-info->ned, nint((krm-klm)*info->rfact)), 2);
  //
  // ** Dynamic approach to decide the minimum gap ratio **
  def1 = 1.0;
  def2 = 0.7;
  up = 1.0;
  dw = 0.7;
  info->avgm =
    (info->avgm * ((double) info->nloop - 1) +
     (double) (info->klan - info->k + 1)) / ((double) info->nloop);
  z1 = (info->ptres) / (info->tol * info->anrm);
  z2 = 1.0 / info->cfac;
  if (z2 < 1.0) {
    tmp = def2;
    t = abs(krm - klm) * def2;
  } else if (z1 < 1.0) {
    tmp = def1;
    t = abs(krm - klm) * def1;
  } else {
    tmp3 =
      log(z1 +
          sqrt(z1 - 1.0) * sqrt(z1 + 1.0)) / (2.0 * (info->avgm));
    tmp2 =
      log(z2 + sqrt(z2 - 1.0) * sqrt(z2 + 1.0)) / (info->klan -
                                                   info->k + 1);

    tmp = tmp2 / tmp3;
    tmp = pow(tmp, 2.0);
    tmp = atan(tmp) * M_2_PI;
    tmp = dw + (up - dw) * tmp;

    t = abs(krm - klm) * tmp;
  }
  igap = t;
  if (igap > 2 && igap + kept > nd && info->crat > 0.0) {
    ss = trl_min_gap_ratio(info, nd, tind, res);
    if (ss > info->crat)
      igap = imax2(2, nd - kept - 1);
  }
  if (lohi > 0) {
    ss = merit_maxred(kri, kli);
    for (i = klm; i <= krm - igap; i++) {
      for (j = i + igap; j <= krm; j++) {
        tmp = merit_maxred(j, i);
        if (tmp > ss) {
          ss = tmp;
          kli = i;
          kri = j;
        }
      }
    }
  } else {
    ss = merit_maxred(kli, kri);
    for (i = klm; i <= krm - igap; i++) {
      for (j = i + igap; j <= krm; j++) {
        tmp = merit_maxred(i, j);
        if (tmp > ss) {
          ss = tmp;
          kli = i;
          kri = j;
        }
      }
    }
  }
  *kl = kli;
  *kr = kri;
}

/*
  Purpose
  =======
  Determine the search range --
  used by the schemes that performs brute-force search.
*/
void trl_restart_search_range(int nd, double *lambda, double *res,
                               trl_info * info, int ncl, int ncr,
                               int *lohi, int tind, int *klm, int *krm) {
  int j, klmi, krmi;
  double bnd;

  klmi = imax2(ncl, 0);
  krmi = imin2(ncr, nd - 1);
  bnd = info->tol * info->anrm;
  *lohi = info->lohi;
  // make sure there is some distance between the boundary and the
  // target Ritz value
  if (info->lohi > 0) {
    // save high end
    krmi = imin2(imax2(info->maxlan - info->ned,
                       (info->maxlan + info->nec) / 2) - 1,
                 imin2(krmi, tind - 1));
    while (krmi + krmi > ncl + ncr && res[krmi] < bnd) {
      krmi--;
    }
  } else if (info->lohi < 0) {
    // save lower end
    klmi = imax2(imin2(info->ned, (info->maxlan + info->nec) / 2) - 1,
                 imax2(tind + 1, klmi));
    while (klmi + klmi < ncl + ncr && res[klmi] < bnd) {
      klmi++;
    }
  } else {
    // save both ends
    if (tind - klmi < krmi - tind) {
      // tind is closer to klmi
      *lohi = -1;
      klmi = tind + 1;
    } else {
      // tind is closer to krmi
      *lohi = 1;
      krmi = tind - 1;
    }
    j = info->locked + klmi + nd - krmi + 1;
    if (j > 0) {
      j = j / 2;
      // should be bounded by 0 and nd-1
      klmi = imax2(0, klmi - j);
      krmi = imin2(krmi + j, nd - 1);
    }
  }
  *klm = klmi;
  *krm = krmi;
}

/*
  Purpose
  =======
  Try to determine the minimal gap ratio need to compute all wanted
  eigenvalues
*/
double trl_min_gap_ratio(trl_info * info, int nd, int tind, double *res) {
  double gamma;

  gamma = info->maxmv * (info->nec + 1.0) / info->ned - info->matvec;
  if (gamma < info->klan) {
    gamma =
      fmax2(2.0,
            (info->maxmv - info->matvec) / ((double) (info->ned - info->nec)));
  }
  return fmin2(log(res[tind] / (info->tol * info->anrm)) / gamma, 0.5);
}
