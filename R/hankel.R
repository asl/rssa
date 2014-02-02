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

#   Routines for normal hankel SSA

hankel <- function(X, L) {
  if (is.matrix(X) && nargs() == 1) {
     L <- nrow(X); K <- ncol(X); N <- K + L - 1
     left  <- c(1:L, L*(2:K))
     right <- c(1+L*(0:(K-1)), ((K-1)*L+2):(K*L))
     v <- sapply(1:N, function(i) mean(X[seq.int(left[i], right[i], by = L-1)]))
     return (v)
  }

  # Coerce output to vector, if necessary
  if (!is.vector(X))
    X <- as.vector(X)
  N <- length(X)
  if (missing(L))
    L <- (N - 1) %/% 2
  K <- N - L + 1
  outer(1:L, 1:K, function(x,y) X[x+y-1])
}

.get.or.create.fft.plan <- function(x) {
  .get.or.create(x, "fft.plan", fft.plan.1d(x$length))
}

.get.or.create.hmat <- function(x) {
  .get.or.create(x, "hmat",
                 new.hmat(.F(x), L = x$window,
                          fft.plan = .get.or.create.fft.plan(x)))
}

.hankelize.one.default <- function(U, V, fft.plan = NULL) {
  L <- length(U); K <- length(V); N = K + L - 1
  fft.plan <- (if (is.null(fft.plan)) fft.plan.1d(N) else fft.plan)
  storage.mode(U) <- storage.mode(V) <- "double"
  .Call("hankelize_one_fft", U, V, fft.plan)
}

.hankelize.one.1d.ssa <- function(x, U, V, fft.plan = NULL) {
  fft.plan <- (if (is.null(fft.plan)) .get.or.create.fft.plan(x) else fft.plan)
  storage.mode(U) <- storage.mode(V) <- "double"
  .Call("hankelize_one_fft", U, V, fft.plan)
}

.hankelize.multi.default <- function(U, V, fft.plan) {
  stopifnot(is.numeric(V))
  storage.mode(U) <- storage.mode(V) <- "double"
  .Call("hankelize_multi_fft", U, V, fft.plan)
}

fft.plan.1d <- function(N) {
  storage.mode(N) <- "integer"
  .Call("initialize_fft_plan", N)
}

is.fft.plan <- function(fft.plan) {
  .Call("is_fft_plan", fft.plan)
}

new.hmat <- function(F,
                     L = (N + 1) %/% 2,
                     fft.plan = NULL) {
  N <- length(F)
  storage.mode(F) <- "double"
  storage.mode(L) <- "integer"
  h <- .Call("initialize_hmat", F, L, if (is.null(fft.plan)) fft.plan.1d(N) else fft.plan)
}

hcols <- function(h) {
  .Call("hankel_cols", h)
}

hrows <- function(h) {
  .Call("hankel_rows", h)
}

is.hmat <- function(h) {
  .Call("is_hmat", h)
}

hmatmul <- function(hmat, v, transposed = FALSE) {
  storage.mode(v) <- "double";
  storage.mode(transposed) <- "logical";
  .Call("hmatmul", hmat, v, transposed);
}

.traj.dim.default <- function(x) {
  c(x$window, x$length - x$window + 1)
}

decompose.1d.ssa <- function(x,
                             neig = min(50, L, K),
                             ...,
                             force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1
  stop("Unsupported SVD method for 1D SSA!")
}

decompose.1d.ssa.svd <- function(x,
                                 neig = min(L, K),
                                 ...,
                                 force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  h <- hankel(.F(x), L = L)

  # Do decomposition
  S <- svd(h, nu = neig, nv = neig)

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)

  x
}

Lcov.matrix <- function(F,
                        L = (N + 1) %/% 2,
                        fft.plan = NULL) {
  N <- length(F)
  storage.mode(F) <- "double"
  storage.mode(L) <- "integer"
  .Call("Lcov_matrix", F, L, if (is.null(fft.plan)) fft.plan.1d(N) else fft.plan)
}

decompose.1d.ssa.eigen <- function(x,
                                   neig = min(50, L, K),
                                   ...,
                                   force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  fft.plan <- .get.or.create.fft.plan(x)

  # Do decomposition
  S <- eigen(Lcov.matrix(.F(x), L = L, fft.plan = fft.plan), symmetric = TRUE)

  # Fix small negative values
  S$values[S$values < 0] <- 0

  # Save results
  .set.decomposition(x,
                     sigma = sqrt(S$values[1:neig]),
                     U = S$vectors[, 1:neig, drop = FALSE])

  x
}

decompose.1d.ssa.propack <- function(x,
                                     neig = min(50, L, K),
                                     ...,
                                     force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.")

  h <- .get.or.create.hmat(x)
  S <- propack.svd(h, neig = neig, ...)

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)

  x
}

decompose.1d.ssa.nutrlan <- function(x,
                                     neig = min(50, L, K),
                                     ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  h <- .get.or.create.hmat(x)

  S <- trlan.svd(h, neig = neig, ...,
                 lambda = .sigma(x), U = .U(x))

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u);

  x
}

calc.v.1d.ssa <- function(x, idx, ...) {
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

    if (any(sigma <= .Machine$double.eps)) {
      sigma[sigma <= .Machine$double.eps] <- Inf
      warning("some sigmas are equal to zero. The corresponding vectors will be zeroed")
    }

    U <- .U(x)[, idx.new, drop = FALSE]
    h <- .get.or.create.hmat(x)
    V[, idx > nV] <- sapply(seq_along(idx.new),
                            function(i) hmatmul(h, U[, i], transposed = TRUE) / sigma[i])
  }

  invisible(V)
}
