#   R package for Singular Spectrum Analysis
#   Copyright (c) 2013 Anton Korobeynikov <asl@math.spbu.ru>
#
#   This program is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public
#   License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option)
#   any later version.
#
#   This program is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied
#   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#   PURPOSE.  See the GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public
#   License along with this program; if not, write to the
#   Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
#   MA 02139, USA.

igapfill <- function(x, L,
                     groups,
                     kind = c("1d-ssa"),
                     fill = mean(x, na.rm = TRUE), eps = 1e-6, numiter = 0,
                     ..., cache = TRUE) {
  kind <- match.arg(kind)
  if (!identical(kind, "1d-ssa"))
    stop("gapfilling is supported for 1D SSA only")

  N <- length(x)
  # Determine the indices of missing values
  na.idx <- which(is.na(x))

  # Obtain the initial approximation
  F <- x
  F[na.idx] <- (if (length(fill) > 1) fill[na.idx] else fill)
  s <- ssa(F, kind = kind, L = L, ...)
  groups <- unlist(groups)
  r <- reconstruct(s, groups = list(groups), cache = cache)
  stopifnot(length(r) == 1)
  F[na.idx] <- r[[1]][na.idx]

  # Do the actual iterations until the convergence (or stoppping due to number
  # of iterations)
  it <- 0
  repeat {
    is <- clone(s, copy.cache = FALSE, copy.storage = FALSE)
    .set(s, "F", F)
    .set(s, "Fattr", attributes(F))
    r <- reconstruct(s, groups = list(groups), ..., cache = cache)
    stopifnot(length(r) == 1)
    rF <- x
    rF[na.idx] <- r[[1]][na.idx]

    it <- it + 1
    if ((numiter > 0 && it >= numiter) || max((F-rF)^2) < eps)
      break
    F <- rF
  }

  F
}
