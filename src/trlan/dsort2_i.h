#ifndef __DSORT2_H
#define __DSORT2_H

////
void dsort2(int N, double *ARRAY1, double *ARRAY2);
//
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing
// order of ARRAY1.
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//

void dsort2i(int N, double *ARRAY1, double *ARRAY2);
//
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the decreasing
// order of ARRAY1.
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//


void dsort2a(int N, double *ARRAY1, double *ARRAY2);
//
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing
// order of abs(ARRAY1).
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
////
void dsort2s(int N, double s, double *ARRAY1, double *ARRAY2);
//
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing
// order of abs(ARRAY1-s).
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// s       (input) DOUBLE PRECISION
//          On entry, specifies the reference value s.
//
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
//
////
void dsort2su(int N, double s, double *ARRAY1, double *ARRAY2);
//
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing order
// of ARRAY1-s if ARRAY1-s is non-negative. Negative ARRAY1-s are ordered after
// those with non-negative ARRAY1-s and in the increasing order of ARRAY1.
//
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
////
void dsort2sd(int N, double s, double *ARRAY1, double *ARRAY2);
//
// Purpose
// =======
// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing order
// of abs(ARRAY1-s) if ARRAY1-s is non-positive. Positive ARRAY1-s are ordered
// after those with non-positive ARRAY1-s and in the increasing order of -ARRAY1.
//
// Arguments
// =========
// N       (input) INTEGER
//          On entry, specifies the size of the arrays.
//
// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.
//
// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
//          On entry, contains the array to be sorted.
//          On exit, contains the sorted array.

#endif
