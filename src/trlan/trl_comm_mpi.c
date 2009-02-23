/*
// ZTRLan routine (version 1.0)
// Lawrence Berkeley National Lab.
//
*/
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "mpi.h"

#include "trlan.h"
#include "trlan_i.h"
#include "trlaux_i.h"
#include "trl_comm_i.h"


////
//void trl_init_info(trl_info *info, int nrow, int mxlan, int lohi, int ned, 
//                    int nopts, ... ) {
void trl_init_info(trl_info * info, int nrow, int mxlan, int lohi,
		   int ned, double tol, int restart, int maxmv,
		   int mpicom)
{
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
//           If provided, specifies the thick-restarting scheme, 1-8. Default is 1.
//
// mxmv     (optionial) integer
//           If provided, specifies the maximum number of matrix-vector 
//           multiplication allowed. By default, mxmv is set to be 
//           (info->ntot)*(info->ned).
//
// mpicom   (optional) integer
//           If provided, specifites the MPI communicator. By default, it is a 
//           duplicate of MPI_COMM_WORLD. In sequential case, this is set to 0.
//
// ..
// .. executable statements ..
    int id, ie;
    //va_list argptr;
    //va_start( argptr,nopts );
    //if( nopts > 0 ) {
    if (tol > 0) {
	//info->tol = va_arg( argptr, double );
	info->tol = tol;
	if (info->tol <= DBL_MIN) {
	    info->tol = DBL_EPSILON;
	} else if (info->tol > 1.0) {
	    info->tol = min(0.1, 1.0 / (info->tol));
	}
    } else {
	info->tol = sqrt(DBL_EPSILON);
    }
    //if( nopts > 1 ) {
    if (restart > 0) {
	//info->restart = va_arg( argptr,int );
	info->restart = restart;
	//printf( "info->restart=%d\n",info->restart );
    } else {
	info->restart = 0;
    }
    //if( nopts > 2 ) {
    if (maxmv > 0) {
	//info->maxmv = va_arg( argptr,int );
	info->maxmv = maxmv;
	//printf( "info->maxmv=%d\n",info->maxmv );
    } else {
	info->maxmv = min(max(info->ntot, 1000), 1000 * info->ned);
    }
    //if( nopts > 3 ) {
    if (mpicom > 0) {
	//info->mpicom = va_arg( argptr,int );
	info->mpicom = mpicom;
    } else {
	info->stat =
	    MPI_Comm_dup(MPI_COMM_WORLD, (MPI_Comm *) (&(info->mpicom)));
	if (info->stat == MPI_SUCCESS) {
	    info->stat = 0;
	} else {
	    printf
		("TRL_INIT_INFO: unable to duplication MPI_COMM_WORLD.\n");
	    printf("This probably indicate MPI_INIT is not called.\n");
	    return;
	}
    }
    //printf( "info->mpicom=%d\n",info->mpicom );
    //va_end( argptr );
    info->maxlan = mxlan;
    if (mxlan <= ned) {
	info->maxlan = ned + max(ned, 6);
    }
    info->lohi = lohi;
    info->ned = ned;
    info->nloc = nrow;
    info->stat = MPI_Allreduce(&nrow, &(info->ntot), 1, MPI_INT, MPI_SUM,
			       (MPI_Comm) (info->mpicom));
    if (info->stat == MPI_SUCCESS) {
	info->stat = 0;
    } else {
	printf("TRL_INIT_INFO: failed to perform MPI_ALLREDUCE, ierr=%d\n",
	       info->stat);
	info->stat = -1;
	return;
    }
    info->guess = 0;
    info->nec = 0;
    info->locked = info->nec;
    info->matvec = 0;
    info->nloop = 0;
    info->north = 0;
    info->nrand = 0;
    info->flop = 0;
    info->rflp = 0;
    info->flop_h = 0;
    info->rflp_h = 0;
    info->flop_r = 0;
    info->rflp_r = 0;
    info->clk_rate = CLOCKS_PER_SEC;
    info->clk_max = (int) (pow(2.0, (8.0 * sizeof(clock_t) - 1.0)));;
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

    info->stat = MPI_Comm_rank((MPI_Comm) (info->mpicom), &(info->my_pe));
    if (info->stat != MPI_SUCCESS) {
	printf("TRL_INIT_INFO: failed to perform MPI_COMM_RANK, ierr=%d\n",
	       info->stat);
	info->stat = -2;
	return;
    }
    info->stat = MPI_Comm_size((MPI_Comm) (info->mpicom), &(info->npes));
    if (info->stat != MPI_SUCCESS) {
	printf("TRL_INIT_INFO: failed to perform MPI_COMM_SIZE, ierr=%d\n",
	       info->stat);
	info->stat = -3;
	return;
    }
    ie = MPI_Allreduce(&(info->stat), &(id), 1, MPI_INT, MPI_MIN,
		       (MPI_Comm) (info->mpicom));
    if (ie == MPI_SUCCESS) {
	info->stat = id;
    } else {
	printf("TRL_INIT_INFO: failed to perform MPI_ALLREDUCE, ierr=%d\n",
	       info->stat);
	info->stat = -1;
    }
    info->cpflag = 0;
    strcpy(info->oldcpf, "");
    info->log_io = 99;
    strcpy(info->log_file, "");
    info->cpio = 98;
    strcpy(info->cpfile, "");
    return;
}

////
void trl_g_sum_(int mpicom, int nelm, double *x, double *y)
{
//
// Purpose
// =======
// trl_g_sum performs global sum in the parallel environment by calling MPI_ALLREDUCE 
// to perform the actual task. This wrapper is needed to avoid actually linking with 
// MPI library when MPI is not used.
//
// Arguments
// =========
// mpicom   (input) integer
//           On entry, specifites the MPI communicator.
//
// nelm     (input) integer
//           On entry, specify the number of element in x.
//           
// x        (input) double
//           On entry, stores the variable to be summed.
//           
// y        (output) double
//           On exit, stores the results.
//           
// .. 
// .. Parameters ..
    static int c__1 = 1;
//
// ..
// .. Local scalars ..
    int ierr;
//
// ..
// .. Executable statements ..
    ierr =
	MPI_Allreduce(x, y, nelm, MPI_DOUBLE, MPI_SUM, (MPI_Comm) mpicom);
    trl_dcopy(nelm, y, c__1, x, c__1);
    if (ierr != MPI_SUCCESS) {
	printf("TRL_G_SUM: MPI_ALLREDUCE failed with error code %d.\n",
	       ierr);
	MPI_Abort((MPI_Comm) mpicom, ierr);
    }
    return;
}

////
int trl_sync_flag_(int mpicom, int inflag)
{
//
// Purpose
// =======
// Given an integer value, this function returns the minimum value of
// all the PEs.
//
// Arguments
// =========
// mpicom   (input) integer
//           On entry, specifites the MPI communicator.
// 
// inflag   (input) integer
//           On entry, contains the value to be compared with. 
// ..
// .. Local scalars ..
    int outflag, ierr;
//
// ..
// .. Executable statements ..
    ierr =
	MPI_Allreduce(&inflag, &outflag, 1, MPI_INT, MPI_MIN,
		      (MPI_Comm) mpicom);
    if (ierr != MPI_SUCCESS) {
	printf("TRL_SYNC_FLAG: MPI_ALLREDUCE failed with error code %d.",
	       ierr);
	MPI_Abort((MPI_Comm) mpicom, ierr);
    }
    return outflag;
}

////
void trl_g_dot_(int mpicom, int nrow, double *v1, int ld1, int m1,
		double *v2, int ld2, int m2, double *rr, double *wrk)
{
//
// Purpose
// =======
// Implements a distributed version of BLAS routine dgemv, which is used 
// to compute dot-products by TRLAN. This function performs wrk = [V1, V2]'*rr.
//
// Arguments
// =========
// mpicom     (input) integer
//             On entry, specifies MPI communicator.
//
// nrow       (input) integer
//             On entry, specifies the number of rows on the local processor.
//             
// v1         (input) double precision array of dimenstion (ld1,m1)
//             On entry, stores the first part of the matrix.
//             
// ld1        (input) integer
//             On entry, specifies the leading dimension of the array v1.
//             
// m1         (input) integer
//             On entry, specifies the number of columns in v1.
//             
// v2         (input) double precision array of dimension (ld2,m2)
//             On entry, stores the second part of the matrix.
//             
// ld2        (input) integer
//             On entry, specifies the leading dimension of v2.
//             
// m2         (input) integer
//             On entry, specifies the number of columns in v2.
//             
// rr         (input) double precision array of length (m1+m2)
//             On entry, stores the vector to be multiplied.
//             
// wrk        (output) double precision array of length nrow.
//             On exit, stores the results of this operation, partial summs are 
//             stored at the end of this array.
//
// ..
// .. Parametesrs ..
    static char trans = 'T';
    static int c__1 = 1;
    static double zero = 0.0, one = 1.0;
//
// ..
// .. Local scalars ..
    int i, j, npe, nd, m1p1;
//
// ..
// .. Executable statements ..
    nd = m1 + m2;
    if (nd <= 0)
	return;
    if (ld1 < nrow || ld2 < nrow) {
	printf("trl_g_dot incorrect array sizes\n");
	return;
    }
    m1p1 = m1 + 1;
    if (m1 > 2) {
	trl_dgemv(&trans, nrow, m1, one, v1, ld1, rr, c__1, zero, &wrk[nd],
		  c__1);
    } else if (m1 == 2) {
	wrk[nd] = zero;
	wrk[nd + 1] = zero;
	for (i = 0; i < nrow; i++) {
	    wrk[nd] += v1[i] * rr[i];
	    wrk[nd + 1] += v1[i + ld1] * rr[i];
	}
    } else if (m1 == 1) {
	wrk[nd] = trl_ddot(nrow, v1, c__1, rr, c__1);
    }
    if (m2 > 2) {
	trl_dgemv(&trans, nrow, m2, one, v2, ld2, rr, c__1, zero,
		  &wrk[nd + m1p1 - 1], c__1);
    } else if (m2 == 2) {
	wrk[nd + m1p1 - 1] = zero;
	wrk[nd + m1p1] = zero;
	for (i = 0; i < nrow; i++) {
	    wrk[nd + m1p1 - 1] += v2[i] * rr[i];
	    wrk[nd + m1p1] += v2[i + ld2] * rr[i];
	}
    } else if (m2 == 1) {
	wrk[nd + m1p1 - 1] = trl_ddot(nrow, v2, c__1, rr, c__1);
    }
    if ((MPI_Comm) mpicom == MPI_COMM_SELF) {
	npe = 1;
    } else {
	i = MPI_Comm_size((MPI_Comm) mpicom, &npe);
    }
    if (npe > 1) {
	i = MPI_Allreduce(&wrk[nd], wrk, nd, MPI_DOUBLE, MPI_SUM,
			  (MPI_Comm) mpicom);
	if (i != MPI_SUCCESS) {
	    printf("TRL_G_DOT: MPI_ALLREDUCE failed with error code %d.\n",
		   i);
	    MPI_Abort((MPI_Comm) mpicom, i);
	}
    } else {
	trl_dcopy(nd, &wrk[nd], c__1, wrk, c__1);
    }
}
