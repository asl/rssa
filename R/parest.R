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

roots2pars <- function(roots) {
  out <- list(roots = roots,
              periods = 2*pi / Arg(roots),
              frequencies = Arg(roots) / (2*pi),
              moduli = Mod(roots),
              rates = log(Mod(roots)))

  # Fix erroneous BIG period in case of real root
  out$periods[abs(Arg(roots)) < .Machine$double.eps^.5] <- Inf

  class(out) <- "fdimpars.1d"
  out
}

print.fdimpars.1d <- function(x, ...) {
  cat("   period     rate   |    Mod     Arg  |     Re        Im\n")
  for (i in seq_along(x$roots)) {
    cat(sprintf("% 9.3f  % 8.6f | % 7.5f  % 3.2f | % 7.5f  % 7.5f\n",
                x$periods[i],
                x$rates[i],
                x$moduli[i],
                x$frequencies[i] * 2*pi,
                Re(x$roots[i]),
                Im(x$roots[i])))
  }
}

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

  r <- exp(1i * acos(median(scos)))
  roots2pars(r)
}

tls.solve <- function(A, B) {
  stopifnot(ncol(A) == ncol(B))
  r <- ncol(A)
  V <- svd(cbind(A, B))$v[, 1:r, drop = FALSE]
  qr.solve(V[1:r,, drop = FALSE], V[-(1:r),, drop = FALSE])
}

parestimate.esprit <- function(U, method = c("esprit-ls", "esprit-tls")) {
  method <- match.arg(method)
  solver <- if (identical(method, "esprit-ls")) {
        qr.solve
      } else if (identical(method, "esprit-tls")) {
        tls.solve
      }

  Z <- solver(U[-nrow(U), ], U[-1, ])
  r <- eigen(Z, only.values = TRUE)$values
  roots2pars(r)
}

parestimate.1d.ssa <- function(x, groups, method = c("pairs", "esprit-ls", "esprit-tls"),
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
    } else if (identical(method, "esprit-ls") || identical(method, "esprit-tls")) {
      res <- parestimate.esprit(x$U[, group], method = method)
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
