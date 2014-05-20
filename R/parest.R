#   R package for Singular Spectrum Analysis
#   Copyright (c) 2012 Anton Korobeynikov <asl@math.spbu.ru>
#   Copyright (c) 2013 Alex Shlemov <shlemovalex@gmail.com>
#   Copyright (c) 2013 Konstantin Usevich <konstantin.usevich@statmod.ru>
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

parestimate.pairs <- function(U, normalize = FALSE) {
  # Sanity check
  stopifnot(ncol(U) == 2)

  # Now calculate the cosines between the consecutive segments
  U1 <- apply(U[-1, ], 2, diff)
  U2 <- apply(U[-nrow(U), ], 2, diff)
  scos <- rowSums(U1*U2) / sqrt(rowSums(U1*U1)) / sqrt(rowSums(U2*U2))

  # Some ad-hoc test for checking the sanity of the results
  mres <- mad(2*pi/acos(scos)) / median(2*pi/acos(scos))
  if (mres > 0.3)
    warning("too big deviation of estimates, period estimates might be unreliable")

  r <- exp(1i * acos(median(scos)))

  if (normalize) r <- r / abs(r)

  roots2pars(r)
}

tls.solve <- function(A, B) {
  stopifnot(ncol(A) == ncol(B))
  r <- ncol(A)
  V <- svd(cbind(A, B))$v[, 1:r, drop = FALSE]
  qr.solve(t(V[1:r,, drop = FALSE]), t(V[-(1:r),, drop = FALSE]))
}

shift.matrix <- function(U,
                         wmask = NULL,
                         topology = Inf,
                         solve.method = c("ls", "tls")) {
  solve.method <- match.arg(solve.method)
  solver <- switch(solve.method,
                   ls = qr.solve,
                   tls = tls.solve)

  if (is.null(wmask))
    wmask <- rep(TRUE, nrow(U))

  lm.mask <- wmask[-1] & wmask[-length(wmask)]
  lm1.mask <- c(lm.mask, FALSE)[wmask]
  lm2.mask <- c(FALSE, lm.mask)[wmask]

  lmA <- U[lm1.mask,, drop = FALSE]
  lmB <- U[lm2.mask,, drop = FALSE]
  if (topology == length(wmask) && wmask[1] && wmask[length(wmask)]) {
    lmA <- rbind(lmA, U[length(wmask),, drop = FALSE])
    lmB <- rbind(lmB, U[1,, drop = FALSE])
  }

  Conj(solver(lmA, lmB))
}

parestimate.esprit <- function(U,
                               wmask = NULL,
                               topology = Inf,
                               normalize = FALSE,
                               method = c("esprit-ls", "esprit-tls")) {
  method <- match.arg(method)
  Z <- shift.matrix(U,
                    wmask = wmask,
                    topology = topology,
                    solve.method = switch(method,
                                          `esprit-ls` = "ls",
                                          `esprit-tls` = "tls"))

  r <- eigen(Z, only.values = TRUE)$values

  if (normalize) r <- r / abs(r)

  roots2pars(r)
}

parestimate.1d.ssa <- function(x, groups, method = c("pairs", "esprit-ls", "esprit-tls"),
                               subspace = c("column", "row"),
                               normalize.roots = NULL,
                               ...,
                               drop = TRUE) {
  method <- match.arg(method)

  if (missing(groups))
    groups <- 1:min(nsigma(x), nu(x))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  subspace <- match.arg(subspace)
  if (identical(subspace, "column")) {
    span <- .colspan
    wmask <- x$wmask
  } else if (identical(subspace, "row")) {
    span <- .rowspan
    wmask <- x$fmask
  }

  if (is.null(normalize.roots))
    normalize.roots <- x$circular || inherits(x, "toeplitz.ssa")

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    if (identical(method, "pairs")) {
      if (is.shaped(x))
        stop("`pairs' parameter estimation method is not implemented for shaped SSA case yet")
      if (length(group) != 2)
        stop("can estimate for pair of eigenvectors only using `pairs' method")
      res <- parestimate.pairs(span(x, group), normalize = normalize.roots)
    } else if (identical(method, "esprit-ls") || identical(method, "esprit-tls")) {
      res <- parestimate.esprit(span(x, group),
                                wmask = wmask,
                                topology = ifelse(x$circular, x$length, Inf),
                                normalize = normalize.roots,
                                method = method)
    }
    out[[i]] <- res
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

parestimate.toeplitz.ssa <- `parestimate.1d.ssa`
parestimate.mssa <- function(x, groups, method = c("pairs", "esprit-ls", "esprit-tls"),
                             subspace = c("column", "row"),
                             normalize.roots = NULL,
                             ...,
                             drop = TRUE) {
  subspace <- match.arg(subspace)

  if (identical(subspace, "row"))
    stop("Row space parameter estimation is not implemented for MSSA yet")
  parestimate.1d.ssa(x = x, groups = groups, method = method,
                     subspace = subspace,
                     normalize.roots = normalize.roots,
                     ...,
                     drop = drop)
}

shift.matrices.2d <- function(U, L,
                              wmask = NULL,
                              topology = c(Inf, Inf),
                              solve.method = c("ls", "tls")) {
  solve.method <- match.arg(solve.method)
  solver <- switch(solve.method,
                   ls = qr.solve,
                   tls = tls.solve)

  if (is.null(wmask)) {
    wmask <- matrix(TRUE, L[1], L[2])
  }

  lm.mask <- wmask[-1, , drop = FALSE] & wmask[-nrow(wmask),, drop = FALSE]
  lm1.mask <- as.vector(rbind(lm.mask, FALSE)[wmask])
  lm2.mask <- as.vector(rbind(FALSE, lm.mask)[wmask])

  mu.mask <- wmask[, -1, drop = FALSE] & wmask[, -ncol(wmask), drop = FALSE]
  mu1.mask <- as.vector(cbind(mu.mask, FALSE)[wmask])
  mu2.mask <- as.vector(cbind(FALSE, mu.mask)[wmask])

  lmA <- U[lm1.mask,, drop = FALSE]
  lmB <- U[lm2.mask,, drop = FALSE]
  if (topology[1] == L[1]) {
    lmc.mask <- wmask[1,, drop = FALSE] & wmask[L[1],, drop = FALSE]
    lmc1.mask <- as.vector(rbind(matrix(FALSE, L[1] - 1, L[2]), lmc.mask)[wmask])
    lmc2.mask <- as.vector(rbind(lmc.mask, matrix(FALSE, L[1] - 1, L[2]))[wmask])
    lmA <- rbind(lmA, U[lmc1.mask,, drop = FALSE])
    lmB <- rbind(lmB, U[lmc2.mask,, drop = FALSE])
  }

  muA <- U[mu1.mask,, drop = FALSE]
  muB <- U[mu2.mask,, drop = FALSE]
  if (topology[2] == L[2]) {
    muc.mask <- wmask[, 1, drop = FALSE] & wmask[, L[2], drop = FALSE]
    muc1.mask <- as.vector(cbind(matrix(FALSE, L[1], L[2] - 1), muc.mask)[wmask])
    muc2.mask <- as.vector(cbind(muc.mask, matrix(FALSE, L[1], L[2] - 1))[wmask])
    muA <- rbind(muA, U[muc1.mask,, drop = FALSE])
    muB <- rbind(muB, U[muc2.mask,, drop = FALSE])
  }

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

parestimate.esprit2d <- function(U, L,
                                 wmask = NULL,
                                 topology = c(Inf, Inf),
                                 normalize = c(FALSE, FALSE),
                                 method = c("esprit-diag-ls", "esprit-diag-tls",
                                            "esprit-memp-ls", "esprit-memp-tls"),
                                 beta = 8) {
  method <- match.arg(method)
  solve.method <- switch(method,
                         `esprit-diag-ls`  =, `esprit-memp-ls`  = "ls",
                         `esprit-diag-tls` =, `esprit-memp-tls` = "tls")

  Zs <- shift.matrices.2d(U,
                          L = L,
                          wmask = wmask,
                          topology = topology,
                          solve.method = solve.method)

  r <- switch(method,
              `esprit-diag-ls` =, `esprit-diag-tls` = est_exp_2desprit(Zs, beta = beta),
              `esprit-memp-ls` =, `esprit-memp-tls` = est_exp_memp_new(Zs, beta = beta))

  for (k in 1:2)
    if (normalize[k]) r[[k]] <- r[[k]] / abs(r[[k]])

  out <- lapply(r, roots2pars)
  class(out) <- "fdimpars.2d"
  out
}

parestimate.2d.ssa <- function(x, groups,
                               method = c("esprit-diag-ls", "esprit-diag-tls",
                                          "esprit-memp-ls", "esprit-memp-tls"),
                               subspace = c("column", "row"),
                               normalize.roots = NULL,
                               ...,
                               beta = 8,
                               drop = TRUE) {
  method <- match.arg(method)
  if (missing(groups))
    groups <- 1:min(nsigma(x), nu(x))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  subspace <- match.arg(subspace)
  if (identical(subspace, "column")) {
    span <- .colspan
    wmask <- x$wmask
    window <- x$window
  } else if (identical(subspace, "row")) {
    span <- .rowspan
    wmask <- x$fmask
    window <- ifelse(x$circular, x$length - x$window + 1, x$window)
  }

  if (is.null(normalize.roots))
    normalize.roots <- x$circular | inherits(x, "toeplitz.ssa")
  if (length(normalize.roots) > 2)
    warning("Incorrect argument length: length(normalize.roots) > 2, two leading values will be used")
  normalize.roots <- rep(normalize.roots, 2)[1:2]

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]

    out[[i]] <- parestimate.esprit2d(span(x, group),
                                     L = window,
                                     wmask = wmask,
                                     topology = ifelse(x$circular, x$length, Inf),
                                     normalize = normalize.roots,
                                     method = method,
                                     beta = beta)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

print.fdimpars.2d <- function(x, ...) {
  cat("x: period     rate   | y: period     rate\n")
  for (i in seq_along(x[[1]]$roots)) {
    cat(sprintf("% 9.3f  % 8.6f | % 9.3f  % 8.6f \n",
                x[[1]]$periods[i],
                x[[1]]$rate[i],
                x[[2]]$periods[i],
                x[[2]]$rate[i]))
  }
}

parestimate <- function(x, ...)
  UseMethod("parestimate")
