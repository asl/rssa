/*
 *   R package for Singular Spectrum Analysis
 *   Copyright (c) 2009-2010 Anton Korobeynikov <asl@math.spbu.ru>
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

#include <R.h>
#include <Rinternals.h>

#include "extmat.h"

SEXP is_extmat(SEXP ptr) {
  SEXP ans;
  ext_matrix *e = NULL;

  PROTECT(ans = allocVector(LGLSXP, 1));
  LOGICAL(ans)[0] = 1;

  /* object is an external pointer */
  if (TYPEOF(ptr) != EXTPTRSXP)
    LOGICAL(ans)[0] = 0;

  /* tag should be 'external matrix' */
  if (LOGICAL(ans)[0] &&
      R_ExternalPtrTag(ptr) != install("external matrix"))
    LOGICAL(ans)[0] = 0;

  /* pointer itself should not be null */
  if (LOGICAL(ans)[0]) {
    e = R_ExternalPtrAddr(ptr);
    if (!e)
      LOGICAL(ans)[0] = 0;
  }

  /* finally, type should be nonnull */
  if (LOGICAL(ans)[0] && e && e->type == NULL)
    LOGICAL(ans)[0] = 0;

  UNPROTECT(1);

  return ans;
}

SEXP extmat_rows(SEXP ptr) {
  SEXP tchk;
  SEXP ans = NILSXP;

  /* Perform a type checking */
  PROTECT(tchk = is_extmat(ptr));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e = R_ExternalPtrAddr(ptr);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = e->nrow(e->matrix);
    UNPROTECT(1);
  } else
    error("pointer provided is not an external matrix");

  UNPROTECT(1);

  return ans;
}

SEXP extmat_cols(SEXP ptr) {
  SEXP tchk;
  SEXP ans = NILSXP;

  /* Perform a type checking */
  PROTECT(tchk = is_extmat(ptr));

  if (LOGICAL(tchk)[0]) {
    ext_matrix *e = R_ExternalPtrAddr(ptr);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = e->ncol(e->matrix);
    UNPROTECT(1);
  } else
    error("pointer provided is not an external matrix");

  UNPROTECT(1);

  return ans;
}
