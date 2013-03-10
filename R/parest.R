#   R package for Singular Spectrum Analysis
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

parestimate.pairs <- function(U) {
  # Sanity check
  stopifnot(ncol(U) == 2)

  # Now calculate the cosines between the consecutive segments
  U1 <- apply(U[-1, ], 2, diff)
  U2 <- apply(U[-nrow(U), ], 2, diff)
  scos <- rowSums(U1*U2) / sqrt(rowSums(U1*U1)) / sqrt(rowSums(U2*U2))

  # Some ad-hoc test for checking the sanity of the results
  mres <- mad(2*pi/acos(scos))
  if (mres > 1)
    warning("too big deviation of estimates, period estimates might be unreliable")

  list(periods=2*pi/acos(median(scos)))
}

parestimate.esprit <- function(U) {
  Z <- qr.solve(U[-nrow(U),], U[-1, ])
  r <- eigen(Z, only.values = TRUE)$values
  list(periods=2*pi/Arg(r), moduli = Mod(r))
}

parestimate.1d.ssa <- function(x, groups, method = c("pairs", "esprit-ls"),
                               ...,
                               drop = TRUE) {
  method <- match.arg(method)

  if (missing(groups))
    groups <- 1:min(nlambda(x), nu(x))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    if (identical(method, "pairs")) {
      if (length(group) != 2)
        stop("can estimate for pair of eigenvectors only using `pairs' method")
      res <- parestimate.pairs(x$U[, group])
    } else if (identical(method, "esprit-ls")) {
      res <- parestimate.esprit(x$U[, group])
    }
    out[[i]] <- res
  }

  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

parestimate.toeplitz.ssa <- `parestimate.1d.ssa`

parestimate <- function(x, ...)
  UseMethod("parestimate")
