/*
  ZTRLan routine (version 1.0)
  Lawrence Berkeley National Lab.
*/

#include <R.h>
#include <Rmath.h>

#include "trlan.h"
#include "trlan_i.h"
#include "trlaux_i.h"
#include "trl_comm_i.h"

/* This file contains most auxillary routines which are not extensively
   used in during the normal operations of the TRLan, e.g., printing,
   error checking, etc. */

int close_file(FILE * fp, int err1, int err2) {
  if (fclose(fp) != 0) {
    return err2;
  }
  return err1;
}

void trl_open_logfile(trl_info * info) {
  char filename[TRLAN_STRING_LEN];

  if (info->log_file != 0 && strlen(info->log_file) > 0) {
    trl_pe_filename(TRLAN_STRING_LEN, filename, info->log_file, info->my_pe,
                    info->npes);
    info->log_fp = fopen(filename, "w");
  } else {
    info->log_fp = stdout;
  }
}

void trl_reopen_logfile(trl_info * info)
{
  char filename[TRLAN_STRING_LEN];

  if (info->log_file != 0 && strlen(info->log_file) > 0) {
    trl_pe_filename(TRLAN_STRING_LEN, filename, info->log_file, info->my_pe,
                     info->npes);
    info->log_fp = fopen(filename, "a");
  } else {
    info->log_fp = stdout;
  }
}

void trl_close_logfile(trl_info * info)
{
  if (info->log_fp != NULL && info->log_fp != stdout) {
    fclose(info->log_fp);
  }
  info->log_fp = NULL;
}

void trl_open_cptfile(trl_info * info)
{
  char filename[TRLAN_STRING_LEN];

  if (info->cpfile != 0 && strlen(info->cpfile) > 0) {
    trl_pe_filename(TRLAN_STRING_LEN, filename, info->cpfile, info->my_pe,
                    info->npes);
    info->cpt_fp = fopen(filename, "w");
  } else {
    info->cpt_fp = stdout;
  }
}

void trl_close_cptfile(trl_info * info)
{
  if (info->cpt_fp != stdout) {
    fclose(info->cpt_fp);
  }
  info->cpt_fp = NULL;
}

void trl_print_int(trl_info * info, char *title, int size_array,
                   int *array, int inc)
{
  int i;

  fprintf(info->log_fp, "PE %d : %s", info->my_pe, title);
  if (size_array > 2) {
    fprintf(info->log_fp, "\n");
  }
  for (i = 0; i < size_array; i += inc) {
    fprintf(info->log_fp, "%10d", array[i]);
    if ((i % 8) == 7)
      fprintf(info->log_fp, "\n");
  }
  if (((size_array - 1) % 8) != 7)
    fprintf(info->log_fp, "\n");
}

void trl_print_real(trl_info * info, char *title, int size_array,
                    double *array, int inc)
{
  int i;

  fprintf(info->log_fp, "PE %d : %s", info->my_pe, title);
  if (size_array > 1) {
    fprintf(info->log_fp, "\n");
  }
  for (i = 0; i < size_array; i += inc) {
    fprintf(info->log_fp, " %10.7e", array[i]);
    if ((i % 8) == 7)
      fprintf(info->log_fp, "\n");
  }
  if (((size_array - 1) % 8) != 7)
    fprintf(info->log_fp, "\n");
}

void trl_print_progress(trl_info * info)
{
  fprintf(info->log_fp, "MATVEC: %10d,    Nloop: %10d,     Nec: %10d\n",
          info->matvec, info->nloop, info->nec);
  fprintf(info->log_fp, "Reorth: %10d,    Nrand: %10d,    Ierr: %10d\n",
          info->north, info->nrand, info->stat);
  fprintf(info->log_fp,
          "Target: %10.3e,   ResNrm: %10.3e,    CFact: %10.3e\n",
          info->trgt, info->tres, info->crat);
}

void trl_check_orth(trl_info * info, int nrow, double *v1, int ld1,
                    int j1, double *v2, int ld2, int j2, double *wrk,
                    int lwrk)
{
  double one = 1.0, zero = 0.0;
  long c__1 = 1;
  int i, j, k, jnd;
  double nrmfro, nrminf;

  jnd = j1 + j2;
  nrmfro = zero;
  nrminf = zero;

  if (jnd <= 0)
    return;

  if (lwrk < (jnd + jnd))
    error("TRL_CHECK_ORTH: workspace too small.\n");

  fprintf(info->log_fp,
          "TRL_CHECK_ORTH: check orthogonality of %d basis vectors.\n",
          jnd);
  /* check orthognality of the basis vectors */
  for (i = 0; i < j1; i++) {
    trl_g_dot(info->mpicom, nrow, v1, ld1, i + 1, v2, ld2, 0,
              &v1[i * nrow], wrk);
    wrk[i] = wrk[i] - one;
    if (info->verbose > 7) {
      fprintf(info->log_fp, "Orthogonality level of v(%d) ..\n",
              i + 1);
      for (j = 0; j <= i; j++) {
        fprintf(info->log_fp, " %10.3e", wrk[j]);
        if ((j % 8) == 7)
          fprintf(info->log_fp, "\n");
      }
      if ((i % 8) != 7)
        fprintf(info->log_fp, "\n");
    }
    nrmfro =
      nrmfro + 2 * trl_ddot(i, wrk, c__1, wrk,
                            c__1) + wrk[i] * wrk[i];
    if (i == 0) {
      wrk[i + 1] = fabs(wrk[i]);
    } else {
      wrk[i + 1] = fmax2(wrk[i], wrk[i - 1]);
    }
    nrminf = fmax2(nrminf, wrk[i + 1]);
  }
  for (i = 0; i < j2; i++) {
    j = j1 + i;
    trl_g_dot(info->mpicom, nrow, v1, ld1, j1, v2, ld2, i + 1,
              &v2[i * nrow], wrk);
    wrk[j] = wrk[j] - one;
    if (info->verbose > 7) {
      fprintf(info->log_fp, "Orthogonality level of v(%d) ..\n",
              j + 1);
      for (k = 0; k <= j; k++) {
        fprintf(info->log_fp, " %10.3e", wrk[k]);
        if ((k % 8) == 7)
          fprintf(info->log_fp, "\n");
      }
      if ((j % 8) != 7)
        fprintf(info->log_fp, "\n");
    }
    nrmfro =
      nrmfro + 2 * trl_ddot(j, wrk, c__1, wrk,
                            c__1) + wrk[j] * wrk[j];
    nrminf = fmax2(nrminf, fabs(wrk[j]));
  }
  fprintf(info->log_fp,
          "Frobenius norm of orthogonality level %10i %4i  %14.5e\n",
          info->matvec, jnd, sqrt(nrmfro));
  fprintf(info->log_fp,
          "Maximum abs. value of orthogonality level is  %14.5e\n",
          nrminf);
}

void
trl_check_recurrence(trl_matprod op,
                     trl_info * info, int nrow, int ncol, double *v1,
                     int ld1, int m1, double *v2, int ld2, int m2,
                     int kept, double *alpha, double *beta, double *wrk,
                     int lwrk, void *lparam)
{
  long c__1 = 1;
  int i__1 = 1;
  double zero = 0.0, one = 1.0;
  int i, ii, j, j1, j2, jnd, mv1;
  char title[TRLAN_STRING_LEN];
  double d__1;
  double *aq = NULL, *qkp1, *cs, *alf, *bet;
  mv1 = m1;
  if (m2 > 0) {
    j2 = m2 - 1;
    j1 = m1;
  } else {
    j2 = 0;
    j1 = m1 - 1;
  }
  jnd = j1 + j2;
  if (lwrk < jnd * 4 + imax2(jnd * 4, nrow))
    error("TRL_CHECK_RECURRENCE: not enough workspace.\n");

  if (lwrk >= jnd * 4 + nrow)
    aq = &wrk[lwrk - nrow];
  else if (lwrk >= jnd * 4)
    aq = Calloc(nrow, double);

  memset(wrk, 0, 4 * jnd * sizeof(double));
  cs = &wrk[jnd];
  alf = &wrk[2 * jnd];
  bet = &wrk[3 * jnd];

  /*
    first type of relation
    A q_i = Alpha_i q_i + Beta_i q_{k+1}
  */
  if (kept < ncol)
    qkp1 = &v1[kept * nrow];
  else
    qkp1 = &v2[(kept - j1) * nrow];

  for (i = 0; i < imin2(j1, kept); i++) {
    op(&nrow, &i__1, &v1[i * nrow], &nrow, aq, &nrow, lparam);
    for (ii = 0; ii < nrow; ii++) {
      alf[i] += aq[ii] * v1[i * nrow + ii];
      aq[ii] -= alpha[i] * v1[i * nrow + ii];
      bet[i] += aq[ii] * aq[ii];
      cs[i] += aq[ii] * qkp1[ii];
      aq[ii] -= beta[i] * qkp1[ii];
      wrk[i] += aq[ii] * aq[ii];
    }
  }
  for (i = 0; i < (kept - j1); i++) {
    j = i + j1;
    op(&nrow, &i__1, &v2[i * nrow], &nrow, aq, &nrow, lparam);
    for (ii = 0; ii < nrow; ii++) {
      alf[j] += aq[ii] * v2[i * nrow + ii];
      aq[ii] -= alpha[j] * v2[i * nrow + ii];
      bet[j] += aq[ii] * aq[ii];
      cs[j] += aq[ii] * qkp1[ii];
      aq[ii] -= beta[j] * qkp1[ii];
      wrk[j] += aq[ii] * aq[ii];
    }
  }
  /*
    the (k+1)st base vector need to orthogonalize against all previous
    vectors
  */
  if (jnd > kept) {
    op(&nrow, &i__1, qkp1, &nrow, aq, &nrow, lparam);
    alf[kept] = trl_ddot(nrow, aq, c__1, qkp1, c__1);
    d__1 = -alpha[kept];
    trl_daxpy(nrow, d__1, qkp1, c__1, aq, c__1);
    for (i = 0; i < imin2(j1, kept); i++) {
      d__1 = -beta[i];
      trl_daxpy(nrow, d__1, &v1[i * nrow], c__1, aq, c__1);
    }
    for (i = 0; i < kept - j1; i++) {
      j = j1 + i;
      d__1 = -beta[j];
      trl_daxpy(nrow, d__1, &v2[i * nrow], c__1, aq, c__1);
    }
    bet[kept] = trl_ddot(nrow, aq, c__1, aq, c__1);
    if (kept + 2 <= j1) {
      cs[kept] =
        trl_ddot(nrow, aq, c__1, &v1[(kept + 1) * nrow], c__1);
      d__1 = -beta[kept];
      trl_daxpy(nrow, d__1, &v1[(kept + 1) * nrow], c__1, aq, c__1);
    } else {
      cs[kept] =
        trl_ddot(nrow, aq, c__1, &v2[(kept + 1 - j1) * nrow],
                 c__1);
      d__1 = -beta[kept];
      trl_daxpy(nrow, d__1, &v2[(kept + 1 - j1) * nrow], c__1, aq,
                c__1);
    }
    wrk[kept] = trl_ddot(nrow, aq, c__1, aq, c__1);
  }
  /*
    the third kind of relation -- normal three term recurrence
    depending the fact that if the lower-bound of loop is less than
    upper bound, the look should not be executed
  */
  for (i = kept + 1; i < j1; i++) {
    op(&nrow, &i__1, &v1[nrow * i], &nrow, aq, &nrow, lparam);
    if (i < (mv1 - 1)) {
      for (ii = 0; ii < nrow; ii++) {
        alf[i] += aq[ii] * v1[i * nrow + ii];
        aq[ii] -=
          (alpha[i] * v1[i * nrow + ii] +
           beta[i - 1] * v1[(i - 1) * nrow + ii]);
        bet[i] += aq[ii] * aq[ii];
        cs[i] += aq[ii] * v1[(i + 1) * nrow + ii];
        aq[ii] -= beta[i] * v1[(i + 1) * nrow + ii];
        wrk[i] += aq[ii] * aq[ii];
      }
    } else {
      for (ii = 0; ii < nrow; ii++) {
        alf[i] += aq[ii] * v1[i * nrow + ii];
        aq[ii] -=
          (alpha[i] * v1[i * nrow + ii] +
           beta[i - 1] * v1[(i - 1) * nrow + ii]);
        bet[i] += aq[ii] * aq[ii];
        cs[i] += aq[ii] * v2[ii];
        aq[ii] -= beta[i] * v2[ii];
        wrk[i] += aq[ii] * aq[ii];
      }
    }
  }
  for (i = imax2(0, kept - j1 + 1); i < j2; i++) {
    j = i + j1;
    op(&nrow, &i__1, &v2[i * nrow], &nrow, aq, &nrow, lparam);
    if (i > 0) {
      for (ii = 0; ii < nrow; ii++) {
        alf[j] += aq[ii] * v2[i * nrow + ii];
        aq[ii] -=
          (beta[j - 1] * v2[(i - 1) * nrow + ii] +
           alpha[j] * v2[i * nrow + ii]);
        bet[j] += aq[ii] * aq[ii];
        cs[j] += aq[ii] * v2[(i + 1) * nrow + ii];
        aq[ii] -= beta[j] * v2[(i + 1) * nrow + ii];
        wrk[j] += aq[ii] * aq[ii];
      }
    } else {
      for (ii = 0; ii < nrow; ii++) {
        alf[j] += aq[ii] * v2[ii];
        aq[ii] -=
          (beta[j - 1] * v1[(j1 - 1) * nrow + ii] +
           alpha[j] * v2[ii]);
        bet[j] += aq[ii] * aq[ii];
        cs[j] += aq[ii] * v2[nrow + ii];
        aq[ii] -= beta[j] * v2[nrow + ii];
        wrk[j] += aq[ii] * aq[ii];
      }
    }
  }

  trl_g_sum(info->mpicom, jnd * 4, wrk, &wrk[jnd * 4]);
  aq[0] = zero;
  for (ii = 0; ii < jnd; ii++) {
    aq[0] += wrk[ii];
  }
  aq[0] = sqrt(aq[0]);
  for (ii = 0; ii < jnd; ii++) {
    wrk[ii] = sqrt(wrk[ii]);
  }
  for (ii = 0; ii < jnd; ii++) {
    if (bet[ii] > zero) {
      if (beta[ii] < zero) {
        bet[ii] = -sqrt(bet[ii]);
      } else {
        bet[ii] = sqrt(bet[ii]);
      }
      cs[ii] = cs[ii] / bet[ii];
    } else {
      bet[ii] = zero;
    }
  }
  strcpy(title, "Alpha computed by TRLAN ..");
  trl_print_real(info, title, jnd, alpha, 1);
  strcpy(title, "Alpha computed explicitly in TRL_CHECK_RECURRENCE ..");
  trl_print_real(info, title, jnd, alf, 1);
  strcpy(title, "Differences in alpha ..");
  d__1 = -one;
  trl_daxpy(jnd, d__1, alpha, c__1, alf, c__1);
  trl_print_real(info, title, jnd, alf, 1);
  strcpy(title, "Beta computed by TRLAN ..");
  trl_print_real(info, title, jnd, beta, 1);
  strcpy(title, "Beta computed explicitly in TRL_CHECK_RECURRENCE ..");
  trl_print_real(info, title, jnd, bet, 1);
  strcpy(title, "Differences in beta ..");
  d__1 = -one;
  trl_daxpy(jnd, d__1, beta, c__1, bet, c__1);
  trl_print_real(info, title, jnd, bet, 1);
  strcpy(title, "Error in Lanczos recurrence (overall) =");
  trl_print_real(info, title, 1, aq, 1);
  if (info->verbose > 7) {
    strcpy(title,
           "|| A q_i - alpha_i q_i - beta_{i-1} q_{i-1} - beta_i q_{i+1} ||..");
    trl_print_real(info, title, jnd, wrk, 1);
    strcpy(title,
           "(A q_i - alpha_i q_i - beta_{i-1} q_{i-1})*q_{i+1}/beta_i ..");
    trl_print_real(info, title, jnd, cs, 1);
    strcpy(title, "Sine of the angles ..");
    for (ii = 0; ii < jnd; ii++) {
      cs[ii] = cs[ii] * cs[ii];
      if (cs[ii] < one) {
        cs[ii] = sqrt(one - cs[ii]);
      } else {
        cs[ii] = -one;
      }
    }
    trl_print_real(info, title, jnd, cs, 1);
  }
  if (lwrk < jnd * 4 + nrow)
    Free(aq);
}

int trl_write_checkpoint(char *filename, int nrow, double *alpha,
                         double *beta, double *evec, int lde, int me,
                         double *base, int ldb, int nb)
{
  int jnd, i, j;
  FILE *io_fp;

  jnd = me + nb - 1;
  io_fp = fopen(filename, "w");
  if (io_fp == NULL)
    error("TRL_WRITE_CHECKPOINT: failed to open file: %s.\n", filename);

  if (fwrite(&nrow, sizeof(nrow), 1, io_fp) < 1) {
    return close_file(io_fp, -223, -222);
  }
  if (fwrite(&jnd, sizeof(jnd), 1, io_fp) < 1) {
    return close_file(io_fp, -223, -222);
  }

  for (i = 0; i < jnd; i++) {
    if (fwrite(&alpha[i], sizeof(alpha[i]), 1, io_fp) < 1) {
      return close_file(io_fp, -223, -222);
    }
  }
  for (i = 0; i < jnd; i++) {
    if (fwrite(&beta[i], sizeof(beta[i]), 1, io_fp) < 1) {
      return close_file(io_fp, -223, -222);
    }
  }
  for (i = 0; i < me; i++) {
    for (j = 0; j < nrow; j++) {
      if (fwrite
          (&evec[i * nrow + j], sizeof(evec[i * nrow + j]), 1,
           io_fp) < 1) {
        return close_file(io_fp, -223, -222);
      }
    }
  }
  for (i = 0; i < nb; i++) {
    for (j = 0; j < nrow; j++) {
      if (fwrite
          (&base[i * nrow + j], sizeof(base[i * nrow + j]), 1,
           io_fp) < 1) {
        return close_file(io_fp, -223, -222);
      }
    }
  }
  return close_file(io_fp, 0, -223);
}

int trl_read_checkpoint(char *filename, int nrow, double *evec, int lde,
                        int mev, int *j1, double *base, int ldb, int nbas,
                        int *j2, int nalpha, double *alpha, int nbeta,
                        double *beta)
{
  int i, j;
  FILE *io_fp;

  if (lde < nrow || ldb < nrow)
    error("TRL_READ_CHECKPOINT: leading dimensions too small.\n");

  /* open file */
  io_fp = fopen(filename, "r");
  if (io_fp == NULL)
    error("TRL_READ_CHECKPOINT: failed to open check-point file %s.\n",
          filename);

  /* read size information */
  if (fread(j1, sizeof(*j1), 1, io_fp) <= 0) {
    return close_file(io_fp, -215, -216);
  }
  if (fread(j2, sizeof(*j2), 1, io_fp) <= 0) {
    return close_file(io_fp, -215, -216);
  }
  if (*j1 != nrow)
    error("TRL_READ_CHECKPOINT: Nrow mismatch.\n");

  if (*j2 > imin2(nalpha, imin2(nbeta, mev + nbas - 1)))
    error("TRL_READ_CHECKPOINT: MAXLAN too small.");

  /* can continue read all data */
  for (i = 0; i < *j2; i++) {
    if (fread(&alpha[i], sizeof(alpha[i]), 1, io_fp) <= 0) {
      return close_file(io_fp, -215, -216);
    }
  }
  for (i = 0; i < *j2; i++) {
    if (fread(&beta[i], sizeof(beta[i]), 1, io_fp) <= 0) {
      return close_file(io_fp, -215, -216);
    }
  }
  *j1 = imin2(mev, *j2);
  *j2 = *j2 - *j1;
  if (*j1 < mev) {
    for (i = 0; i <= *j1; i++) {
      for (j = 0; j < nrow; j++) {
        if (fread
            (&evec[i * nrow + j], sizeof(evec[i * nrow + j]), 1,
             io_fp) <= 0) {
          return close_file(io_fp, -215, -216);
        }
      }
    }
  } else {
    for (i = 0; i < *j1; i++) {
      for (j = 0; j < nrow; j++) {
        if (fread
            (&evec[i * nrow + j], sizeof(evec[i * nrow + j]), 1,
             io_fp) <= 0) {
          return close_file(io_fp, -215, -216);
        }
      }
    }
    for (i = 0; i < *j2; i++) {
      for (j = 0; j < nrow; j++) {
        if (fread
            (&base[i * nrow + j], sizeof(base[i * nrow + j]), 1,
             io_fp) <= 0) {
          return close_file(io_fp, -215, -216);
        }
      }
    }
  }
  return close_file(io_fp, 0, -215);
}

int indchar(char *a, char b) {
  char *t = strchr(a, b);

  if (t != NULL) {
    return (1 + (t - a));
  } else {
    return strlen(a) + 1;
  }
}

void trl_pe_filename(int nlen, char *filename, char *base, int my_rank,
                     int npe) {
  int lead, ndig, len, off;
  char *format;

  ndig = 1;
  lead = npe;
  while (lead > 9) {
    lead /= 10;
    ++ndig;
  }
  len = indchar(base, ' ') - 1;

  if (nlen < len + ndig + 1)
    error("error: not enough space for filename (%d+%d chars).\n",
          len, ndig);

  memset(filename, 0, nlen * sizeof(char));
  strncpy(filename, base, len);
  off = 1 + ndig % 10;
  off = 5 + 2 * off;
  format = Calloc(off, char);
  sprintf(format, "%%s%%0%d.%dd", ndig, ndig);
  sprintf(filename, format, filename, my_rank);
  Free(format);
  return;
}

void trl_time_stamp(FILE * fp) {
  time_t clck;

  clck = time(NULL);
  fprintf(fp, "                                                  %s",
          asctime(localtime(&clck)));
}
