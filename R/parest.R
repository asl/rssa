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
              rate = log(Mod(roots)))

  class(out) <- "fdimpars.1d"
  out
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
      res <- parestimate.esprit(x$U[, group, drop = FALSE], method = method)
    }
    out[[i]] <- res
  }

  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

parestimate.toeplitz.ssa <- `parestimate.1d.ssa`

shift.matrices.2d <- function(U, umask, solve.method = c("ls", "tls")) {
  solve.method <- match.arg(solve.method)
  solver <- if (identical(solve.method, "ls")) {
        qr.solve
      } else if (identical(solve.method, "tls")) {
        tls.solve
      }

  lm.mask <- umask[-1, , drop = FALSE] & umask[-nrow(umask),, drop = FALSE]
  lm1.mask <- as.vector(rbind(lm.mask, FALSE)[umask])
  lm2.mask <- as.vector(rbind(FALSE, lm.mask)[umask])

  mu.mask <- umask[, -1, drop = FALSE] & umask[, -ncol(umask), drop = FALSE]
  mu1.mask <- as.vector(cbind(mu.mask, FALSE)[umask])
  mu2.mask <- as.vector(cbind(FALSE, mu.mask)[umask])


  lmA <- U[lm1.mask,, drop = FALSE]
  lmB <- U[lm2.mask,, drop = FALSE]


  muA <- U[mu1.mask,, drop = FALSE]
  muB <- U[mu2.mask,, drop = FALSE]


  Zx = solver(lmA, lmB)
  Zy = solver(muA, muB)

  list(Zx = Zx, Zy = Zy)
}

est_exp_2desprit <- function(Zs, beta = 8) {
  Z <- (1-beta) * Zs$Zx + beta * Zs$Zy
  Ze <- eigen(Z, symmetric = FALSE)
  Tinv <- Ze$vectors

  list(diag(qr.solve(Tinv, Zs$Zx %*% Tinv)), diag(qr.solve(Tinv, Zs$Zy %*% Tinv)))
}

est_exp_memp_new <- function(Zs, beta = 8) {
  Z <- (1-beta) * Zs$Zx + beta * Zs$Zy
  Ze <- eigen(Z)
  Zxe <- eigen(Zs$Zx)
  Zye <- eigen(Zs$Zy)
  Px <- max.col(t(abs(qr.solve(Ze$vectors, Zxe$vectors))))
  Py <- max.col(t(abs(qr.solve(Ze$vectors, Zxe$vectors))))

  list(Zxe$values[Px], Zye$values[Py])
}

parestimate.esprit2d <- function(U, umask,
                                 method = c("esprit-diag-ls", "esprit-diag-tls",
                                            "esprit-memp-ls", "esprit-memp-tls"),
                                 beta = 8) {
  method <- match.arg(method)
  solve.method <- if (identical(method, "esprit-diag-ls") || identical(method, "esprit-memp-ls")) {
        "ls"
      } else if (identical(method, "esprit-diag-tls") || identical(method, "esprit-memp-tls")) {
        "tls"
      }

  Zs <- shift.matrices.2d(U, umask, solve.method = solve.method)

  r <- if (identical(method, "esprit-diag-ls") || identical(method, "esprit-diag-tls")) {
        est_exp_2desprit(Zs, beta = beta)
      } else if (identical(method = "esprit-memp-ls") || identical(method = "esprit-memp-tls")) {
        est_exp_memp_new(Zs, beta = beta)
      }

  out <- lapply(r, roots2pars)
  class(out) <- "fdimpars.2d"
  out
}

parestimate.2d.ssa <- function(x, groups,
                               method = c("esprit-diag-ls", "esprit-diag-tls",
                                          "esprit-memp-ls", "esprit-memp-tls"),
                               ...,
                               beta = 8,
                               drop = TRUE) {
  method <- match.arg(method)
  if (missing(groups))
    groups <- 1:min(nlambda(x), nu(x))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]

    umask <- matrix(TRUE, x$window[1], x$window[2])
    out[[i]] <- parestimate.esprit2d(x$U[, group, drop = FALSE],
                                     umask,
                                     method = method,
                                     beta = beta)
  }

  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

parestimate <- function(x, ...)
  UseMethod("parestimate")
