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

.cycle.permutation <- function(v, k = 0) {
  n <- length(v)
  k <- k %% n
  if (k) {
    v <- c(v[(k + 1) : n], v[1 : k])

  }

  v
}

.mdim.cycle.permutation <- function(v, ndim, k = 0) {
  v <- as.array(v)
  d <- dim(v)
  idx <- lapply(d, seq_len)
  idx[[ndim]] <- .cycle.permutation(idx[[ndim]], k = k)
  do.call("[", c(list(v), idx, list(drop = FALSE)))
}

.annulate.row <- function(v, ndim, i = 1, value = 0) {
  v <- as.array(v)
  d <- dim(v)
  idx <- lapply(d, seq_len)
  idx[[ndim]] <- i
  do.call("[<-", c(list(v), idx, list(value = value)))
}

.shifted.matrix.masks <- function(wmask,
                                  ndim,
                                  circular = FALSE) {
  wmask <- as.array(wmask)
  d <- dim(wmask)
  mask <- wmask & .mdim.cycle.permutation(wmask, ndim, 1)

  if (!circular) {
    mask <- .annulate.row(mask, ndim, i = d[ndim], FALSE)
  }

  lind <- array(NA_integer_, dim = d)
  lind[wmask] <- seq_len(sum(wmask))
  rind <- .mdim.cycle.permutation(lind, ndim, 1)

  list(left.mask = lind[mask], right.mask = rind[mask])
}

shift.matrix <- function(U, ...) {
  wmask <- rep(TRUE, nrow(U))
  Conj(.shift.matrix(U, wmask, ndim = 1, ...))
}

.shift.matrix <- function(U, wmask,
                          ndim,
                          circular = FALSE,
                          solve.method = c("ls", "tls")) {
  solve.method <- match.arg(solve.method)
  solver <- switch(solve.method,
                   ls = qr.solve,
                   tls = tls.solve)

  smxs <- .shifted.matrix.masks(wmask, ndim, circular = circular)

  lm.left <- U[smxs$left.mask,, drop = FALSE]
  lm.right <- U[smxs$right.mask,, drop = FALSE]
  solver(lm.left, lm.right)
}

parestimate.esprit <- function(U,
                               wmask = NULL,
                               circular = FALSE,
                               normalize = FALSE,
                               method = c("esprit-ls", "esprit-tls")) {
  method <- match.arg(method)

  if (is.null(wmask))
    wmask <- rep(TRUE, nrow(U))

  solve.method <- switch(method,
                         `esprit-ls` = "ls",
                         `esprit-tls` = "tls")

  Z <- .shift.matrix(U,
                     wmask = wmask,
                     ndim = 1,
                     circular = circular,
                     solve.method = solve.method)

  r <- eigen(Z, only.values = TRUE)$values

  if (normalize) r <- r / abs(r)

  roots2pars(r)
}

parestimate.1d.ssa <- function(x, groups, method = c("esprit-ls", "esprit-tls", "pairs"),
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
                                circular = x$circular,
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
parestimate.mssa <- function(x, groups, method = c("esprit-ls", "esprit-tls", "pairs"),
                             subspace = c("column", "row"),
                             normalize.roots = NULL,
                             ...,
                             drop = TRUE) {
  method <- match.arg(method)
  subspace <- match.arg(subspace)

  if (identical(subspace, "row"))
    stop("Row space parameter estimation is not implemented for MSSA yet")
  parestimate.1d.ssa(x = x, groups = groups, method = method,
                     subspace = subspace,
                     normalize.roots = normalize.roots,
                     ...,
                     drop = drop)
}

.matrix.linear.combination <- function(Zs, beta = 8) {
  if (length(beta) == 1) {
    beta <- beta ^ seq_len(length(Zs) - 1)
  }

  if (length(beta) == length(Zs) - 1) {
    beta <- c(1 - sum(beta), beta)
  }

  Z <- matrix(0., ncol = ncol(Zs[[1]]), nrow = nrow(Zs[[1]]))
  for (i in seq_along(Zs)) {
    Z <- Z + beta[i] * Zs[[i]]
  }

  Z
}

.est.exp.2desprit <- function(Zs, beta = 8) {
  Z <- .matrix.linear.combination(Zs, beta)
  Ze <- eigen(Z, symmetric = FALSE)
  Tinv <- Ze$vectors

  lapply(seq_along(Zs),
         function(i) diag(qr.solve(Tinv, Zs[[i]] %*% Tinv)))
}

.est.exp.memp.new <- function(Zs, beta = 8) {
  Z <- .matrix.linear.combination(Zs, beta)
  Ze <- eigen(Z)
  Zse <- lapply(Zs, eigen, symmetric = FALSE)
  Ps <- lapply(seq_along(Zs),
               function(i) max.col(t(abs(qr.solve(Ze$vectors, Zse[[i]]$vectors)))))  # TODO Use assignment problem here

  # for (P in Ps) {
  #   stopifnot(length(P) == length(unique(P)))
  # }

  lapply(seq_along(Zs),
         function(i) Zse[[i]]$values[Ps[[i]]])
}

parestimate.esprit2d <- function(U, L,
                                 wmask = NULL,
                                 circular = c(FALSE, FALSE),
                                 normalize = c(FALSE, FALSE),
                                 solve.method = c("ls", "tls"),
                                 pairing.method = c("diag", "memp"),
                                 beta = 8) {
  if (is.null(wmask))
    wmask <- array(TRUE, dim = L)

  wmask <- as.array(wmask)
  d <- dim(wmask)

  solve.method <- match.arg(solve.method)
  pairing.method <- match.arg(pairing.method)

  Zs <- lapply(seq_along(d),
               function(ndim) {
                 .shift.matrix(U,
                               wmask = wmask,
                               ndim = ndim,
                               circular = circular[ndim],
                               solve.method = solve.method)
               })


  pairer <- switch(pairing.method, diag = .est.exp.2desprit, memp = .est.exp.memp.new)
  r <- pairer(Zs, beta = beta)

  for (k in seq_along(d))
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
  solve.method <- switch(method,
                         `esprit-diag-ls`  =, `esprit-memp-ls`  = "ls",
                         `esprit-diag-tls` =, `esprit-memp-tls` = "tls")
  pairing.method <- switch(method,
                    `esprit-diag-ls` =, `esprit-diag-tls` = "diag",
                    `esprit-memp-ls` =, `esprit-memp-tls` = "memp")

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
                                     circular = x$circular,
                                     normalize = normalize.roots,
                                     solve.method = solve.method,
                                     pairing.method = pairing.method,
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
