#ifndef __TRLANAUX_H
#define __TRLANAUX_H

/*
  Purpose
  =======
  Closes the file handler.

  Arguments
  =========
  fp      (input) pointer to file
          On entry, points to the file to be closed.

   err1   (input) integer
          On entry, specifies the return value when closing file failed.

   err2   (input) integer
          On entry, specifies the return value when closing file succeeded.
*/
int close_file(FILE * fp, int err1, int err2);

/*
   Purpose:
   ========
   Opens log file.

   Arguments:
   ==========
   info    (input) pointer to the structure trl_info_
            On entry, points to the data structure to store the information
            about the eigenvalue problem and the progress of TRLAN.
*/
void trl_open_logfile(trl_info * info);

/*
   Purpose:
   ========
   Opens log file for append.

   Arguments:
   ==========
   info    (input) pointer to the structure trl_info_
            On entry, points to the data structure to store the information
            about the eigenvalue problem and the progress of TRLAN.
*/

void trl_reopen_logfile(trl_info * info);

/*
   Purpose:
   ========
   Closes log file.

   Arguments:
   ==========
   info    (input) pointer to the structure trl_info_
            On entry, points to the data structure to store the information
            about the eigenvalue problem and the progress of TRLAN.
*/

void trl_close_logfile(trl_info * info);

/*
   Purpose:
   ========
   Opens check points file.

   Arguments:
   ==========
   info    (input) pointer to the structure trl_info_
            On entry, points to the data structure to store the information
            about the eigenvalue problem and the progress of TRLAN.
*/
void trl_open_cptfile(trl_info * info);

/*
   Purpose:
   ========
   Closes check point file.

   Arguments:
   ==========
   info    (input) pointer to the structure trl_info_
            On entry, points to the data structure to store the information
            about the eigenvalue problem and the progress of TRLAN.


*/
void trl_close_cptfile(trl_info * info);

/*
   Purpose:
   ========
   Print an integer array for debugging.

   Arguments:
   ==========
   info        (input) pointer to the structure trl_info_
                On entry, points to the data structure to store the information
                about the eigenvalue problem and the progress of TRLAN

   title       (input) character string
                On entry, specifies the title of the information to be printed.

   size_array  (input) integer
                On entry specifies, the number of integers to be printed.

   array       (input) integer array of length ((size_array-1)*inc+1)
                On entry, contains the integer to be printed.

   inc         (input) integer
                On entry, specifies how the index to array should be incremented.
*/
void trl_print_int(trl_info * info, char *title, int size_array,
                   int *array, int inc);

/*
   Purpose
   =======
   Print a double precision array for debugging.

   Arguments:
   ==========
   info        (input) pointer to the structure trl_info_
                On entry, points to the data structure to store the information
                about the eigenvalue problem and the progress of TRLAN

   title       (input) character string
               On entry, specifies the title of the information to be printed.

   size_array  (input) integer
               On entry specifies, the number of doubles to be printed.

   array       (input) double array of length ((size_array-1)*inc+1)
               On entry, contains the doubles to be printed.

   inc         (input) integer
               On entry, specifies how the index to array should be incremented.
*/
void trl_print_real(trl_info * info, char *title, int size_array,
                    double *array, int inc);

/*
   Purpose
   =======
   Print the current progress of eigenvalue solution

   Arguments:
   ==========
   info        (input) pointer to the structure trl_info_
                On entry, points to the data structure to store the information
                about the eigenvalue problem and the progress of TRLAN
*/
void trl_print_progress(trl_info * info);

/*
   Purpose:
   ========
   Check orthogonality of the basis.

   Arguments:
   ==========
   info     (input) pointer to the structure trl_info_
            On entry, points to the data structure to store the information
            about the eigenvalue problem and the progress of TRLAN.

   nrow     (input) integer
            On entry, specifies the problem size, i.e., the number of rows in
            v1 and v2.

   v1       (input) double precision array of diimension (ld1,j1)
            On entry, contains the first part of the basis.

   ld1      (input) integer
            On entry, specifies the leading diimension of v1.

   j1       (input) integer
            On entry, specifies the last column index of v1, containing the
            basis.

   v2       (input) double precision array of dimension (ld2,j2)
            On entry, contains the second part of the basis.

   ld2      (input) integer
            On entry, specifies the leading dimension of v2.

   j2       (input) integer
            On entry, specifies the last column index of v2, containing the
            basis.

   wrk      (workspace) double precision array of length (lwrk)

   lwrk     (input) integer
            On entry, specifies the size of workspace.*/
void trl_check_orth(trl_info * info, int nrow, double *v1, int ld1,
                    int j1, double *v2, int ld2, int j2, double *wrk,
                    int lwrk);


/*
   Purpose
   =======
   Check Lanczos recurrence relation for debug purpose.

   Arguments:
   ==========
   op       (input) function pointer
            On entry, points to the matrix-vector multiplication routine.

   info     (input) pointer to the structure trl_info_
            On entry, points to the data structure to store the information
            about the eigenvalue problem and the progress of TRLAN.

   nrow     (input) integer
            On entry, specifies the problem size, i.e., the number of ros in v1
            and v2.

   ncol     (input) integer
            On entry, specifies the maximum number of eigenvalues that can be
            stored.

   v1       (input) double precision array of dimension (ld1,m1)
            On entry, contains the first part of basis.

   ld1      (input) integer
            On entry, specifies the leading dimension of v1.

   m1       (input) integer
            On entry, specifies the last column index of v1 that contains a base
            vector.

   v2       (input) double precision array of dimension (ld2,m2)
            On entry, contains the second part of basis.

   ld2      (input) integer
            On entry, specifies the leading dimension of v2.

   m2       (input) integer
            On entry, specifies the last column index of v2 that contains a base
            vector.

   kept     (input) integer
             On entry, specifies the number of basis kept at the last restart.

   alpha    (input) integer
             On entry, contains the alpha values computed so far.

   beta     (input) integer
             On entry, contains the beta values computed so far.

   wrk      (workspace) double precision vector of length (lwrk)

   lwrk     (input) integer
             On entry, specifies the size of workspace.


*/
void trl_check_recurrence(trl_matprod op,
                          trl_info * info, int nrow, int ncol, double *v1,
                          int ld1, int m1, double *v2, int ld2, int m2,
                          int kept, double *alpha, double *beta, double *wrk,
                          int lwrk, void* lparam);

/*
   Purpose
   =======
   Write a check-point file.

   Arguments
   =========
   filename   (input) character string
               On entry, specifies the name of checkpoint file.

   nrow       (input) integer
               On entry, specifies the problem size.

   alpha      (input) double precision array of length (me+nb-1)
               On entry, contains the alpha values computed so far.

   beta       (input) double precision array of length (me+ne-1)
               On entry, contains the beta values computed so far.

   evec       (input) double precision array of dimensioni (lde,me)
               On entry, contains the first part of basis vectors.

   lde        (input) integer
               On entry, specifies the leading dimension of evec.

   me         (input) integer
               On entry, specifies the last column index of evec, that contains
               a base vector.

   base       (input) double precision array of dimension (ldb,nb)
               On entry, contains the second part of basis.

   ldb        (input) integer
               On entry, specifies the leading dimension of base.

   nb         (input) integer
              On entry, specifies the last column index of base, that contains
              a base vector.
*/
int trl_write_checkpoint(char *filename, int nrow, double *alpha,
                         double *beta, double *evec, int lde, int me,
                         double *base, int ldb, int nb);


/*
   Purpose
   =======
   Read check-point file

   Arguments:
   ==========
   filename   (input) character string
               On entry, specifies the name of checkpoint file.

   nrow       (input) integer
               On entry, specifies the problem size.

   evec       (output) double precision array of dimensioni (lde,j1)
               On exit, contains the first part of basis vectors stored in checkpoint.

   lde         (input) integer
               On entry, specifies the leading dimension of evec.

   mev        (input) integer
               On entry, specifies the number of eigenvalues converged.

   j1         (input) integer
               On entry, specifies the last column index of evec, that contains a base vector.

   base       (output) double precision array of dimension (ldb,nbas)
               On exit, contains the second part of basis stored in the checkpoint.

   ldb        (input) integer
               On entry, specifies the leading dimension of base.

   nbas       (input) integer
               On entry, specifies the number of columns in base.

   j2         (input) integer
               On entry, specifies the last column index of base, that contains a base vector.

   nalpha     (input) integer
               On entry, specifies the size of alpha

   alpha      (output) double precision array of length (nalpha)
               On exit, contains the alpha values stored in checkpoint.

   nbeta      (input) integer
               On entry, specifies the size of beta.

   beta       (output) double precision array of length (nbeta)
               On exit, contains the beta values stored in checkpoint.
*/
int trl_read_checkpoint(char *filename, int nrow, double *evec, int lde,
                        int mev, int *j1, double *base, int ldb, int nbas,
                        int *j2, int nalpha, double *alpha, int nbeta,
                        double *beta);

/*
   Purpose
   =======
   Generates file name from a base name and the PE number.

   Arguments:
   ==========
   nlen         (input) integer
                 On entry, specifies the size of filiename.

   filename     (output) character string of length <= nlen.
                 On exit, specifies the file name.

   base         (input) character string
                 On entry, specifies the leading part of the file name.

   my_rank      (input) integer
                 On entry, specifies the PE number.

   npe          (input) integer
                 On entry, specifies the number of processors.
*/
void trl_pe_filename(int nlen, char *filename, char *base, int my_rank,
                     int npe);

/*
   Purpose
   =======
   Print the current date and time.

   Arguments:
   ==========
   fp         (input) pointer to a file
               On entry, points to the file to output.
*/
void trl_time_stamp(FILE * iou);

#endif
