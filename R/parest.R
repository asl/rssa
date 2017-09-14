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
  Conj(qr.solve(t(V[1:r,, drop = FALSE]), t(V[-(1:r),, drop = FALSE])))
}

.cycle.permutation <- function(v, k = 0) {
  n <- length(v)
  k <- k %% n
  if (k) {
    v <- c(v[(k + 1):n], v[1:k])

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

.pairs <- function(x, groups,
                   subspace = c("column", "row"),
                   normalize.roots = NULL,
                   ...,
                   drop) {
  if (missing(groups))
    groups <- 1:min(nsigma(x), nu(x))

  subspace <- match.arg(subspace)

  if (is.null(normalize.roots))
    normalize.roots <- x$circular || inherits(x, "toeplitz.ssa")

  if (is.shaped(x)) {
    stop("`pairs' parameter estimation method is not implemented for shaped SSA case yet")
  }

  if (inherits(x, "cssa")) {
    stop("`pairs' parameter estimation method is not implemented for Complex SSA case yet")
  }

  if (identical(subspace, "column")) {
    span <- .colspan
  } else if (identical(subspace, "row")) {
    if (inherits(x, "mssa")) {
      stop("row space `pairs' parameter estimation method is not implemented for MSSA yet")
    }
    span <- .rowspan
  }

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    if (length(group) != 2)
      stop("can estimate for pair of eigenvectors only using `pairs' method")
    out[[i]] <- parestimate.pairs(span(x, group), normalize = normalize.roots)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

parestimate.1d.ssa <- function(x, groups,
                               method = c("esprit", "pairs"),
                               subspace = c("column", "row"),
                               normalize.roots = NULL,
                               dimensions = NULL,
                               solve.method = c("ls", "tls"),
                               ...,
                               drop = TRUE) {
  method <- match.arg(method)
  solve.method <- match.arg(solve.method)

  if (missing(groups))
    groups <- 1:min(nsigma(x), nu(x))

  subspace <- match.arg(subspace)

  if (is.null(normalize.roots))
    normalize.roots <- x$circular || inherits(x, "toeplitz.ssa")

  if (identical(method, "pairs")) {
    .pairs(x, groups = groups,
           subspace = subspace,
           normalize.roots = normalize.roots,
           ...,
           drop = drop)
  } else if (identical(method, "esprit")) {
    parestimate.nd.ssa(x, groups = groups,
                       subspace = subspace,
                       normalize.roots = normalize.roots,
                       dimensions = c(x = 1),
                       ...,
                       solve.method = solve.method,
                       drop = drop)
  }
}

parestimate.toeplitz.ssa <- parestimate.1d.ssa
parestimate.mssa <- parestimate.1d.ssa
parestimate.cssa <- parestimate.1d.ssa

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

# TODO Use a solution of assignment problem
.simple.assignment <- function(mx) {
  mx <- as.matrix(mx)
  stopifnot(ncol(mx) == nrow(mx))
  d <- nrow(mx)
  stopifnot(all(mx > -Inf))
  res <- rep(0, d)
  for (k in seq_len(d)) {
    maxij <- which(mx  == max(mx), arr.ind = TRUE)[1, ]
    res[maxij[1]] <- maxij[2]
    mx[maxij[1], ] <- -Inf
    mx[, maxij[2]] <- -Inf
  }

  stopifnot(isTRUE(all.equal(sort(res), seq_len(ncol(mx)))))  # Ensure res is permutation

  res
}

.est.exp.memp.new <- function(Zs, beta = 8) {
  Z <- .matrix.linear.combination(Zs, beta)
  Ze <- eigen(Z, symmetric = FALSE)
  Zse <- lapply(Zs, eigen, symmetric = FALSE)
  Ps <- lapply(seq_along(Zs),
               function(i) .simple.assignment(t(abs(qr.solve(Ze$vectors, Zse[[i]]$vectors)))))

  for (P in Ps) {
    stopifnot(length(P) == length(unique(P)))
  }

  lapply(seq_along(Zs),
         function(i) Zse[[i]]$values[Ps[[i]]])
}

.esprit <- function(U,
                    wmask,
                    circular,
                    normalize,
                    dimensions = NULL,
                    solve.method = c("ls", "tls"),
                    pairing.method = c("diag", "memp"),
                    beta = 8) {
  wmask <- as.array(wmask)
  d <- dim(wmask)

  solve.method <- match.arg(solve.method)
  pairing.method <- match.arg(pairing.method)

  if (is.null(dimensions)) {
    dimensions <- seq_along(d)
  }

  if (max(dimensions) > length(d)) {
    stop(sprintf("some of input dimension indices exceed the actual number of object dimensions (%d)",
                 length(d)))
  }

  Zs <- lapply(dimensions,
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

  names(out) <- names(dimensions)
  if (length(names(out)) == 0 || any(names(out) == "")) {
    default.names <- c("x", "y", "z", "t", "u", "s",
                       paste("x",
                             seq_len(max(dimensions)),
                             sep = "_"))
    names(out) <- default.names[dimensions]
  }

  if (length(out) == 1) {
    out <- out[[1]]
  } else {
    class(out) <- "fdimpars.nd"
  }

  out
}

parestimate.nd.ssa <- function(x, groups,
                               method = c("esprit"),
                               subspace = c("column", "row"),
                               normalize.roots = NULL,
                               dimensions = NULL,
                               solve.method = c("ls", "tls"),
                               pairing.method = c("diag", "memp"),
                               beta = 8,
                               ...,
                               drop = TRUE) {
  method <- match.arg(method)
  stopifnot(identical(method, "esprit"))
  solve.method <- match.arg(solve.method)
  pairing.method <- match.arg(pairing.method)

  if (missing(groups))
    groups <- seq_len(min(nsigma(x), nu(x)))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  subspace <- match.arg(subspace)
  if (identical(subspace, "column")) {
    span <- .colspan
    wmask <- .wmask(x)
  } else if (identical(subspace, "row")) {
    span <- function(...) Conj(.rowspan(...))  # TODO Mb include it into the .rowspan method?
    wmask <- .fmask(x)
  }

  if (is.null(dimensions)) {
    dimensions <- seq_len(.dim(x))
  }

  if (is.null(normalize.roots))
    normalize.roots <- x$circular | inherits(x, "toeplitz.ssa")
  if (length(normalize.roots) > .dim(x))
    warning("incorrect argument length: length(normalize.roots) > .dim(x), only leading values will be used")
  normalize.roots <- rep(normalize.roots, .dim(x))[seq_len(.dim(x))]

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    out[[i]] <- .esprit(span(x, group),
                        wmask = wmask,
                        circular = x$circular,
                        normalize = normalize.roots,
                        solve.method = solve.method,
                        pairing.method = pairing.method,
                        beta = beta,
                        dimensions = dimensions)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

print.fdimpars.nd <- function(x, ...) {
  if (length(names(x)) == 0 || any(names(x) == "")) {
    names(x) <- paste("x", seq_along(x), sep = "_")
  }

  header <- paste(sapply(names(x),
                         function(name) sprintf("%8s: period       rate", name)),
                  sep = "", collapse = " | ")
  cat(header)
  cat("\n")

  for (i in seq_along(x[[1]]$roots)) {
    row <- paste(sapply(seq_along(x),
                        function(k) sprintf("% 16.3f  % 8.6f",
                                            x[[k]]$periods[i],
                                            x[[k]]$rate[i])),
                 sep = "", collapse = " | ")
    cat(row)
    cat("\n")
  }
}

parestimate <- function(x, ...)
  UseMethod("parestimate")
