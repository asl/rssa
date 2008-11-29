/*
 *   R package for Singular Spectrum Analysis
 *   Copyright (c) 2008 Anton Korobeynikov <asl@math.spbu.ru>
 *
 *   This program is free software; you can redistribute it
 *   and/or modify it under the terms of the GNU General Public
 *   License as published by the Free Software Foundation;
 *   either version 2 of the License, or (at your option)
 *   any later version.
 *
 *   This program is distributed in the hope that it will be
 *   useful, but WITHOUT ANY WARRANTY; without even the implied
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *   PURPOSE.  See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public
 *   License along with this program; if not, write to the
 *   Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
 *   MA 02139, USA.
 */

/* This is just direct R-to-C translation and will need to be rethought in the
 * future */
void hankelize_one(double *F, double *U, double *V, int *l, int *k) {
  int L = *l, K = *k;
  int N = K + L - 1;
  int i;

  for (i = 0; i < N; ++i) {
    int leftu, rightu, leftv, rightv, l, j;
    double s = 0;

    if (i < L) {
      leftu = i;
      leftv = 0;
    } else {
      leftu = L - 1;
      leftv = i - L + 1;
    }
    if (i < K) {
      rightu = 0;
      rightv = i;
    } else {
      rightu = i - K + 1;
      rightv = K - 1;
    }

    l = leftu - rightu + 1;

    for (j = 0; j < l; ++j) {
      s += U[leftu - j] * V[leftv + j];
    }

    F[i] = s / (double) l;
  }
}
