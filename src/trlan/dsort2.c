/*
  ZTRLan routine (version 1.0)
  Lawrence Berkeley National Lab.
*/

#include <math.h>
#include "dsort2_i.h"

void dsort2(int N, double *ARRAY1, double *ARRAY2)
{
  int IGAP, I, J;
  double TEMP;

  IGAP = N / 2;
  while (IGAP > 0) {
    for (I = IGAP; I < N; I++) {
      J = I - IGAP;
      while (J >= 0) {
        if (ARRAY1[J] > ARRAY1[J + IGAP]) {
          TEMP = ARRAY1[J];
          ARRAY1[J] = ARRAY1[J + IGAP];
          ARRAY1[J + IGAP] = TEMP;
          TEMP = ARRAY2[J];
          ARRAY2[J] = ARRAY2[J + IGAP];
          ARRAY2[J + IGAP] = TEMP;
          J = J - IGAP;
        } else {
          break;
        }
      }
    }
    IGAP = IGAP / 2;
  }
}

void dsort2i(int N, double *ARRAY1, double *ARRAY2)
{
  int IGAP, I, J;
  double TEMP;

  IGAP = N / 2;
  while (IGAP > 0) {
    for (I = IGAP; I < N; I++) {
      J = I - IGAP;
      while (J >= 0) {
        if (ARRAY1[J] < ARRAY1[J + IGAP]) {
          TEMP = ARRAY1[J];
          ARRAY1[J] = ARRAY1[J + IGAP];
          ARRAY1[J + IGAP] = TEMP;
          TEMP = ARRAY2[J];
          ARRAY2[J] = ARRAY2[J + IGAP];
          ARRAY2[J + IGAP] = TEMP;
          J = J - IGAP;
        } else {
          break;
        }
      }
    }
    IGAP = IGAP / 2;
  }
}

void dsort2a(int N, double *ARRAY1, double *ARRAY2)
{
  int IGAP, I, J;
  double TEMP;

  IGAP = N / 2;
  while (IGAP > 0) {
    for (I = IGAP; I < N; I++) {
      J = I - IGAP;
      while (J >= 0) {
        if (fabs(ARRAY1[J]) > fabs(ARRAY1[J + IGAP])) {
          TEMP = ARRAY1[J];
          ARRAY1[J] = ARRAY1[J + IGAP];
          ARRAY1[J + IGAP] = TEMP;
          TEMP = ARRAY2[J];
          ARRAY2[J] = ARRAY2[J + IGAP];
          ARRAY2[J + IGAP] = TEMP;
          J = J - IGAP;
        } else {
          break;
        }
      }
    }
    IGAP = IGAP / 2;
  }
}

void dsort2s(int N, double s, double *ARRAY1, double *ARRAY2)
{
  int IGAP, I, J;
  double TEMP;

  IGAP = N / 2;
  while (IGAP > 0) {
    for (I = IGAP; I < N; I++) {
      J = I - IGAP;
      while (J >= 0) {
        if (fabs(ARRAY1[J] - s) > fabs(ARRAY1[J + IGAP] - s)) {
          TEMP = ARRAY1[J];
          ARRAY1[J] = ARRAY1[J + IGAP];
          ARRAY1[J + IGAP] = TEMP;
          TEMP = ARRAY2[J];
          ARRAY2[J] = ARRAY2[J + IGAP];
          ARRAY2[J + IGAP] = TEMP;
          J = J - IGAP;
        } else {
          break;
        }
      }
    }
    IGAP = IGAP / 2;
  }
}

void dsort2su(int N, double s, double *ARRAY1, double *ARRAY2)
{
  int IGAP, I, J;
  double TEMP, v1, v2, d1, d2, maxd;

  IGAP = N / 2;
  maxd = fabs(ARRAY1[0]);
  for (I = 1; I < N; I++) {
    if (maxd < fabs(ARRAY1[I])) {
      maxd = fabs(ARRAY1[I]);
    }
  }
  while (IGAP > 0) {
    for (I = IGAP; I < N; I++) {
      J = I - IGAP;
      while (J >= 0) {
        v1 = fabs(ARRAY1[J]);
        v2 = fabs(ARRAY1[J + IGAP]);
        d1 = v1 - s;
        d2 = v2 - s;
        if (d1 < 0.0) {
          d1 = maxd + v1;
        }
        if (d2 < 0.0) {
          d2 = maxd + v2;
        }
        if (d1 > d2) {
          TEMP = ARRAY1[J];
          ARRAY1[J] = ARRAY1[J + IGAP];
          ARRAY1[J + IGAP] = TEMP;
          TEMP = ARRAY2[J];
          ARRAY2[J] = ARRAY2[J + IGAP];
          ARRAY2[J + IGAP] = TEMP;
          J = J - IGAP;
        } else {
          break;
        }
      }
    }
    IGAP = IGAP / 2;
  }
/*
// .. end of dsort2a_c_
*/
}

void dsort2sd(int N, double s, double *ARRAY1, double *ARRAY2)
{
  int IGAP, I, J;
  double TEMP, v1, v2, d1, d2, maxd;

  IGAP = N / 2;
  maxd = fabs(ARRAY1[0]);
  for (I = 1; I < N; I++) {
    if (maxd < fabs(ARRAY1[I])) {
      maxd = fabs(ARRAY1[I]);
    }
  }
  maxd = maxd + 1.0;
  while (IGAP > 0) {
    for (I = IGAP; I < N; I++) {
      J = I - IGAP;
      while (J >= 0) {
        v1 = fabs(ARRAY1[J]);
        v2 = fabs(ARRAY1[J + IGAP]);
        d1 = s - v1;
        d2 = s - v2;
        if (d1 < 0.0) {
          d1 = s + (maxd - v1);
        }
        if (d2 < 0.0) {
          d2 = s + (maxd - v2);
        }
        if (d1 > d2) {
          TEMP = ARRAY1[J];
          ARRAY1[J] = ARRAY1[J + IGAP];
          ARRAY1[J + IGAP] = TEMP;
          TEMP = ARRAY2[J];
          ARRAY2[J] = ARRAY2[J + IGAP];
          ARRAY2[J + IGAP] = TEMP;
          J = J - IGAP;
        } else {
          break;
        }
      }
    }
    IGAP = IGAP / 2;
  }
}
