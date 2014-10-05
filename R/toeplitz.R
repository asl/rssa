#   R package for Singular Spectrum Analysis
#   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
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

#   Routines for toeplitz SSA

Lcor <- function(F, L, circular = FALSE) {
  storage.mode(F) <- "double"
  storage.mode(L) <- "integer"
  storage.mode(circular) <- "logical"
  .Call("Lcor", F, L, circular)
}

new.tmat <- function(F, L = (N + 1) %/% 2,
                     circular = FALSE,
                     fft.plan = NULL) {
  N <- length(F)
  R <- Lcor(F, L, circular = circular)

  storage.mode(R) <- "double"
  t <- .Call("initialize_tmat", R, if (is.null(fft.plan)) fft.plan.1d(2*L - 1, L = L) else fft.plan)
}

tcols <- function(t) {
  .Call("toeplitz_cols", t)
}

trows <- function(t) {
  .Call("toeplitz_rows", t)
}

is.tmat <- function(t) {
  .Call("is_tmat", t)
}

tmatmul <- function(tmat, v, transposed = FALSE) {
  storage.mode(v) <- "double"
  storage.mode(transposed) <- "logical"
  .Call("tmatmul", tmat, v, transposed)
}

.hankelize.one.toeplitz.ssa <- .hankelize.one.1d.ssa

.get.or.create.tmat <- function(x) {
  .get.or.create(x, "tmat", new.tmat(F = x$F, L = x$window,
                                     circular = x$circular))
}

.traj.dim.toeplitz.ssa <- .traj.dim.1d.ssa

decompose.toeplitz.ssa.nutrlan <- function(x,
                                           neig = min(50, L, K),
                                           ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  F <- .F(x)
  h <- .get.or.create.hmat(x)

  lambda <- .decomposition(x)$lambda
  U <- .U(x)

  T <- .get.or.create.tmat(x)

  S <- trlan.eigen(T, neig = neig, ...,
                   lambda = lambda, U = U)

  # Save results
  num <- length(S$d)
  sigma <- numeric(num)
  V <- matrix(nrow = .traj.dim(x)[2], ncol = num)
  for (i in 1:num) {
    Z <- hmatmul(h, S$u[, i], transposed = TRUE)
    sigma[i] <- sqrt(sum(Z^2))
    V[, i] <- Z / sigma[i]
  }
  o <- order(sigma, decreasing = TRUE)
  sigma <- sigma[o]
  U <- S$u[, o, drop = FALSE]
  V <- V[, o, drop = FALSE]
  lambda <- S$d[o]

  # Save results
  .set.decomposition(x,
                     sigma = sigma, U = U, V = V, lambda = lambda,
                     kind = "toeplitz.decomposition")

  x
}

decompose.toeplitz.ssa.eigen <- function(x,
                                         neig = min(50, L, K),
                                         ...,
                                         force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decompostion is not supported for this method.")

  # Build hankel matrix
  F <- .F(x)
  h <- .get.or.create.hmat(x)

  # Do decomposition
  C <- toeplitz(Lcor(F, L, circular = x$circular))
  S <- eigen(C, symmetric = TRUE)

  sigma <- numeric(L)
  V <- matrix(nrow = .traj.dim(x)[2], ncol = L)
  for (i in 1:L) {
    Z <- hmatmul(h, S$vectors[,i], transposed = TRUE)
    sigma[i] <- sqrt(sum(Z^2))
    V[, i] <- Z / sigma[i]
  }

  o <- order(sigma[seq_len(neig)], decreasing = TRUE)
  sigma <- sigma[o]
  U <- S$vectors[, o, drop = FALSE]
  V <- V[, o, drop = FALSE]

  # Save results
  .set.decomposition(x,
                     sigma = sigma, U = U, V = V,
                     kind = "toeplitz.decomposition")

  x
}

decompose.toeplitz.ssa.svd <- function(x,
                                       neig = min(50, L, K),
                                       ...,
                                       force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decompostion is not supported for this method.")

  # Build hankel matrix
  F <- .F(x)
  h <- .get.or.create.hmat(x)

  # Do decomposition
  C <- toeplitz(Lcor(F, L, circular = x$circular))
  S <- svd(C, nu = neig, nv = neig)

  sigma <- numeric(neig)
  V <- matrix(nrow = .traj.dim(x)[2], ncol = neig)
  for (i in 1:neig) {
    Z <- hmatmul(h, S$u[,i], transposed = TRUE)
    sigma[i] <- sqrt(sum(Z^2))
    V[, i] <- Z / sigma[i]
  }

  o <- order(sigma, decreasing = TRUE)
  sigma <- sigma[o]
  U <- S$u[, o, drop = FALSE]
  V <- V[, o, drop = FALSE]

  # Save results
  .set.decomposition(x,
                     sigma = sigma, U = U, V = V,
                     kind = "toeplitz.decomposition")

  x
}

decompose.toeplitz.ssa.propack <- function(x,
                                           neig = min(50, L, K),
                                           ...,
                                           force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.");

  S <- propack.svd(.get.or.create.tmat(x), neig = neig, ...)

  h <- .get.or.create.hmat(x)
  num <- length(S$d)
  sigma <- numeric(num)
  V <- matrix(nrow = .traj.dim(x)[2], ncol = num)
  for (i in 1:num) {
    Z <- hmatmul(h, S$u[, i], transposed = TRUE)
    sigma[i] <- sqrt(sum(Z^2))
    V[, i] <- Z / sigma[i]
  }

  o <- order(sigma, decreasing = TRUE)
  sigma <- sigma[o]
  U <- S$u[, o, drop = FALSE]
  V <- V[, o, drop = FALSE]

  # Save results
  .set.decomposition(x,
                     sigma = sigma, U = U, V = V,
                     kind = "toeplitz.decomposition")

  x
}

decompose.toeplitz.ssa <- function(x,
                                   neig = min(50, L, K),
                                   ...,
                                   force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1
  stop("Unsupported SVD method for Toeplitz SSA!")
}

calc.v.toeplitz.ssa <- calc.v.1d.ssa

.rowspan.toeplitz.ssa <- function(x, idx) {
  qr.Q(qr(.V(x)[, idx, drop = FALSE]))
}
