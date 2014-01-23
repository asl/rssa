#   R package for Singular Spectrum Analysis
#   Copyright (c) 2014 Alex Shlemov <shlemovalex@gmail.com>
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


orthopoly <- function(d, L) {
  if (is.character(d)) {
    d <- match.arg(d,
                   choices = c("none", "constant", "linear", "quadratic", "qubic"))
    d <- switch(d,
                none = 0, constant = 1, linear = 2, quadratic = 3, qubic = 4)
  }

  stopifnot(is.numeric(d) && length(d) == 1)

  if (d == 0) {
    # Return matrix with zero columns
    matrix(NA_real_, L, 0)
  } else if (d == 1) {
    matrix(1 / sqrt(L), L, 1)
  } else {
    cbind(1 / sqrt(L), poly(seq_len(L), d - 1))
  }
}

.phmat <- function(x) {
  hmat <- .get.or.create.hmat(x)
  column.projector <- .get(x, "column.projector")
  row.projector <- .get(x, "row.projector")

  matmul <- function(v) {
    v <- v - row.projector %*% crossprod(row.projector, v)
    v <- hmatmul(hmat, v, transposed = FALSE)
    v <- v - column.projector %*% crossprod(column.projector, v)

    v
  }

  tmatmul <- function(v) {
    v <- v - column.projector %*% crossprod(column.projector, v)
    v <- hmatmul(hmat, v, transposed = TRUE)
    v <- v - row.projector %*% crossprod(row.projector, v)

    v
  }

  extmat(matmul, tmatmul, hrows(hmat), hcols(hmat))
}

.get.or.create.phmat <- function(x)
  .get.or.create(x, "phmat", .phmat(x))

hmatmul.wrap <- function(h, mx, transposed = TRUE) {
  if (ncol(mx) == 0)
    return(matrix(NA_real_, if (transposed) hcols(h) else hrows(h), 0))

  apply(mx, 2, hmatmul, hmat = h, transposed = transposed)
}

.init.pssa <- function(x) {
  hmat <- .get.or.create.hmat(x)
  LU <- .get(x, "column.projector")
  RV <- .get(x, "row.projector")

  RU <- hmatmul.wrap(hmat, RV, transposed = FALSE)
  LV <- hmatmul.wrap(hmat, LU, transposed = TRUE) - RV %*% crossprod(RU, LU)
  Rsigma <- sqrt(colSums(RU^2))
  RU <- RU / rep(Rsigma, each = nrow(RU))
  Lsigma <- sqrt(colSums(LV^2))
  LV <- LV / rep(Lsigma, each = nrow(LV))

  .set.decomposition(x,
                     nPR = ncol(RV), nPL = ncol(LU),
                     sigma = c(Rsigma, Lsigma), U = cbind(RU, LU), V = cbind(RV, LV))
}

decompose.pssa <- function(x,
                           neig = min(50, L, K),
                           ...,
                           force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1
  stop("Unsupported SVD method for SSA with projection!")
}

nspecial.pssa <- function(x) {
  sum(unlist(.decomposition(x, c("nPL", "nPR"))))
}

decompose.pssa.svd <- function(x,
                               neig = min(50, L, K),
                               ...,
                               force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1
  nspecial <- nspecial(x)

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > nspecial)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  h <- hankel(.F(x), L = L)

  # Subtract special components
  sigma <- .sigma(x)[seq_len(nspecial)]
  U <- .U(x)[, seq_len(nspecial), drop = FALSE]
  V <- .V(x)[, seq_len(nspecial), drop = FALSE]
  h <- h - U %*% (sigma * t(V))

  # Do decomposition
  S <- svd(h, nu = neig, nv = neig)

  # Save results
  .set.decomposition(x,
                     sigma = c(sigma, S$d), U = cbind(U, S$u), V = cbind(V, S$v))

  x
}

decompose.pssa.eigen <- function(x,
                                 neig = min(50, L, K),
                                 ...,
                                 force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1
  nspecial <- nspecial(x)

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > nspecial)
    stop("Continuation of decomposition is not supported for this method.")

  sigma <- .sigma(x)[seq_len(nspecial)]
  U <- .U(x)[, seq_len(nspecial), drop = FALSE]
  V <- .V(x)[, seq_len(nspecial), drop = FALSE]

  # We will compute (X - P_X) %*% t(X - P_X) = X %*% t(X) - P_X %*% t(X) - X %*% t(P_X) + P_X %*% t(P_X)
  # Get common Lcov matrix, i.e. X %*% t(X)
  Lcov <- Lcov.matrix(.F(x), L = L, fft.plan = .get.or.create.fft.plan(x))
  # Get hankel circulant
  h <- .get.or.create.hmat(x)
  # Compute X %*% t(P_X)
  XtPX <- hmatmul.wrap(h, V, transposed = FALSE) %*% (sigma * t(U))
  # Compute P_X %*% t(P_X)
  PXtPX <- U %*% crossprod(V * rep(sigma, each = nrow(V))) %*% t(U)

  C <- Lcov - XtPX - t(XtPX) + PXtPX

  # Do decomposition
  S <- eigen(C, symmetric = TRUE)

  # Fix small negative values
  S$values[S$values < 0] <- 0

  # Save results
  .set.decomposition(x,
                     sigma = c(sigma, sqrt(S$values[seq_len(neig)])),
                     U = cbind(U, S$vectors[, seq_len(neig), drop = FALSE]),
                     V = V)

  x
}

decompose.pssa.propack <- function(x,
                                   neig = min(50, L, K),
                                   ...,
                                   force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1
  nspecial <- nspecial(x)

  # We will use special (first nspecial) entries below
  sigma <- .sigma(x)
  U <- .U(x)
  V <- .V(x)

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > nspecial)
    stop("Continuation of decompostion is not yet implemented for this method.")

  ph <- .get.or.create.phmat(x)
  S <- propack.svd(ph, neig = neig, ...)

  # Form results
  sigma <- c(sigma[seq_len(nspecial)], S$d)
  U <- cbind(U[, seq_len(nspecial)], S$u)
  V <- cbind(V[, seq_len(nspecial)], S$v)

  # Save results
  .set.decomposition(x,
                     sigma = sigma, U = U, V = V)

  x
}

decompose.pssa.nutrlan <- function(x,
                                   neig = min(50, L, K),
                                   ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  nspecial <- nspecial(x)
  # We will use special (first nspecial) entries below
  sigma <- .sigma(x)
  U <- .U(x)
  V <- .V(x)

  ph <- .get.or.create.phmat(x)
  S <- trlan.svd(ph, neig = neig, ...,
                 lambda = sigma[-seq_len(nspecial)], U = U[, -seq_len(nspecial), drop = FALSE])

  # Form results
  sigma <- c(sigma[seq_len(nspecial)], S$d)
  U <- cbind(U[, seq_len(nspecial)], S$u)
  V <- cbind(V[, seq_len(nspecial)], S$v)

  # Save results
  .set.decomposition(x,
                     sigma = sigma, U = U, V = V)

  x
}

calc.v.pssa <- function(x, idx, ...) {
  N <- x$length; L <- x$window; K <- N - L + 1
  nV <- nv(x)

  V <- matrix(NA_real_, K, length(idx))
  idx.old <- idx[idx <= nV]
  idx.new <- idx[idx > nV]

  if (length(idx.old) > 0) {
    V[, idx <= nV] <- .V(x)[, idx.old]
  }

  if (length(idx.new) > 0) {
    sigma <- .sigma(x)[idx.new]
    U <- .U(x)[, idx.new, drop = FALSE]
    ph <- .get.or.create.phmat(x)

    V[, idx > nV] <- sapply(seq_along(idx.new),
                            function(i) ematmul(ph, U[, i], transposed = TRUE) / sigma[i])
  }

  invisible(V)
}

.colspan.pssa <- function(x, idx) {
  decomposition <- .decomposition(x, c("U", "nPR"))

  span <- decomposition$U[, idx, drop = FALSE]

  # Perform orthogonalization only if it is really needed
  if (any(idx <= decomposition$nPR)) {
    # TODO It can be implemented little more efficiently
    span <- qr.Q(qr(span))
  }

  span
}

.rowspan.pssa <- function(x, idx) {
  decomposition <- .decomposition(x, c("nPR", "nPR"))

  span <- calc.v(x, idx)

  # Perform orthogonalization only if it is really needed
  if (any((idx > decomposition$nPR) & (idx <= decomposition$nPR + decomposition$nPL))) {
    # TODO It can be implemented little more efficiently
    span <- qr.Q(qr(span))
  }

  span
}

enlarge.basis <- function(B, len, solve.method = c("ls", "tls")) {
  solve.method <- match.arg(solve.method)
  solver <- switch(solve.method,
                   ls = qr.solve,
                   tls = tls.solve)
  N <- nrow(B)

  if (ncol(B) == 0) {
    return(matrix(NA_real_, nrow(B) + len, 0))
  }

  B.head <- B[-N,, drop = FALSE]
  B.tail <- B[-1,, drop = FALSE]
  P <- solver(B.head, B.tail)

  B <- rbind(B, matrix(NA, len, ncol(B)))
  for (i in seq_len(len)) {
    B[N + i, ] <- B[N + i - 1, ] %*% P
  }

  B
}

rforecast.pssa <- function(x, groups, len = 1,
                           base = c("reconstructed", "original"),
                           only.new = TRUE,
                           ...,
                           drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  L <- x$window
  K <- x$length - L + 1

  base <- match.arg(base)
  if (missing(groups))
    groups <- as.list(seq_len(min(nsigma(x), nu(x))))

  # Continue decomposition, if necessary
  desired <- .maybe.continue(x, groups = groups, ...)

  # Grab the reconstructed series if we're basing on them
  if (identical(base, "reconstructed"))
    r <- reconstruct(x, groups = groups, ..., cache = cache)

  right.special.triples <- seq_len(.decomposition(x, "nPR"))
  right.groups <- lapply(groups, function(group) intersect(group, right.special.triples))
  nonright.groups <- lapply(groups, function(group) setdiff(group, right.special.triples))

  # Calculate the LRR corresponding to groups
  lf <- lrr(x, groups = nonright.groups, drop = FALSE)
  stopifnot(length(lf) == length(groups))

  all.sigma <- .sigma(x)[seq_len(desired)]
  all.U <- .U(x)[, seq_len(desired), drop = FALSE]
  all.V <- calc.v(x, seq_len(desired))

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    right.group <- if (identical(base, "reconstructed")) right.groups[[i]] else right.special.triples

    # Calculate drifts
    sigma <- all.sigma[right.group]
    U <- all.U[, right.group, drop = FALSE]
    V <- all.V[, right.group, drop = FALSE]
    V <- enlarge.basis(V, len)
    drift <- ((c(-lf[[i]], 1) %*% U) * sigma) %*% t(V[K + seq_len(len),, drop = FALSE])

    # Calculate the forecasted values
    out[[i]] <- apply.lrr(if (identical(base, "reconstructed")) r[[i]] else .get(x, "F"),
                          lf[[i]], len, only.new = only.new, drift = drift)
    out[[i]] <- .apply.attributes(x, out[[i]],
                                  fixup = TRUE,
                                  only.new = only.new, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  # Forecasted series can be pretty huge...
  invisible(out)
}

vforecast.pssa <- function(x, groups, len = 1,
                           only.new = TRUE,
                           ...,
                           drop = TRUE, drop.attributes = FALSE) {
  L <- x$window
  N <- x$length
  K <- x$length - L + 1

  if (missing(groups))
    groups <- as.list(seq_len(min(nsigma(x), nu(x))))

  # Continue decomposition, if necessary
  desired <- .maybe.continue(x, groups = groups, ...)

  right.special.triples <- seq_len(.decomposition(x, "nPR"))
  right.groups <- lapply(groups, function(group) intersect(group, right.special.triples))
  nonright.groups <- lapply(groups, function(group) setdiff(group, right.special.triples))

  # Calculate the LRR corresponding to groups
  lf <- lrr(x, groups = nonright.groups, drop = FALSE)
  stopifnot(length(lf) == length(groups))

  all.sigma <- .sigma(x)[seq_len(desired)]
  all.U <- .U(x)[, seq_len(desired), drop = FALSE]
  all.V <- calc.v(x, seq_len(desired))

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    right.group <- right.groups[[i]]
    nonright.group <- nonright.groups[[i]]

    U <- all.U[, right.group, drop = FALSE]
    V <- all.V[, right.group, drop = FALSE]
    V <- enlarge.basis(V, len + L - 1)
    Z <- V * rep(all.sigma[right.group], each = nrow(V))
    C <- U %*% t(Z)

    U.head.nr <- all.U[-L, nonright.group, drop = FALSE]
    P.head <- tcrossprod(qr.Q(qr(U.head.nr)))

    X <- all.U[, group, drop = FALSE] %*% (all.sigma[group] * t(all.V[, group]))
    X <- cbind(X, matrix(NA, L, len + L - 1))
    for (j in (K + 1) : ncol(X)) {
      tt <- X[-1, j - 1] - C[-L, j]
      tt <- P.head %*% tt
      tt <- c(tt, sum(lf[[i]] * tt))
      X[, j] <- tt + C[, j]
    }
    ser <- hankel(X)

    # Calculate the forecasted values
    out[[i]] <- if (only.new) ser[N + seq_len(len)] else ser[seq_len(N + len)]
    out[[i]] <- .apply.attributes(x, out[[i]],
                                  fixup = TRUE,
                                  only.new = only.new, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  # Forecasted series can be pretty huge...
  invisible(out)
}
