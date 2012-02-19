#   R package for Singular Spectrum Analysis
#   Copyright (c) 2012 Alexandr Shlemov <shlemovalex@gmail.com>
#   Copyright (c) 2012 Anton Korobeynikov <asl@math.spbu.ru>
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

basis2lrf <- function(U) {
  N <- nrow(U);
  lpf <- U %*% t(U[N, , drop = FALSE]);

  (lpf[-N]) / (1 - lpf[N])
}

apply.lrf <- function(F, lrf, len = 1) {
  N <- length(F)
  r <- length(lrf)

  # Sanity check of inputs
  if (r > N)
    stop("Wrong length of LRF")

  # Run the actual LRF
  F <- c(F, rep(NA, len))
  for (i in 1:len)
    F[N+i] <- sum(F[(N+i-r) : (N+i-1)]*lrf)

  F
}

"rforecast.1d-ssa" <- function(this, groups, len = 1,
                               base = c("reconstructed", "original"),
                               ..., cache = TRUE) {
  base <- match.arg(base)
  if (missing(groups))
    groups <- as.list(1:min(nlambda(this), nu(this)))

  # Determine the upper bound of desired eigentriples
  desired <- max(unlist(groups));

  # Continue decomposition, if necessary
  if (desired > min(nlambda(this), nu(this)))
    decompose(this, ..., neig = desired);

  # Grab the reconstructed series if we're basing on them
  if (identical(base, "reconstructed"))
    r <- reconstruct(this, groups = groups, ..., cache = cache)

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]

    # Calculate the LRF corresponding to group
    U <- .get(this, "U")[, group, drop = FALSE]
    lrf <- basis2lrf(U)

    # Calculate the forecasted values
    out[[i]] <- apply.lrf(if (identical(base, "reconstructed")) r[[i]] else .get(this, "F"),
                          lrf, len)
    # FIXME: try to fixup the attributes
  }

  names(out) <- paste(sep = "", "F", 1:length(groups))

  out
}

rforecast.ssa <- function(x, groups, len = 1,
                          base = c("reconstructed", "original"),
                          ..., cache = TRUE) {
  stop("generic recurrent forecast not implemented yet!")
}
