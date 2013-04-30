#   R package for Singular Spectrum Analysis
#   Copyright (c) 2013 Anton Korobeynikov <anton at korobeynikov dot info>
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

.get.or.create.mfft.plan <- function(x) {
  .get.or.create(x, "fft.plan", lapply(x$length, fft.plan.1d))
}

.traj.dim.mssa <- function(x) {
  c(x$window, sum(x$length - x$window + 1))
}

.hmat.striped <- function(x, fft.plan) {
  # FIXME: think about NA's in the end
  N <- x$length; L <- x$window; K <- N - L + 1

  h <- lapply(1:length(N),
              function(idx) new.hmat(x$F[, idx], L = L,
                                     fft.plan = fft.plan[[idx]]))
  b <- c(0, cumsum(K))
  matmul <- function(v) rowSums(sapply(1:length(h),
                                       function(idx) hmatmul(h[[idx]], v[(b[idx]+1):b[idx+1]], transposed = FALSE)))
  tmatmul <- function(v) unlist(lapply(h, hmatmul, v = v, transposed = TRUE))

  extmat(matmul, tmatmul, nrow = L, ncol = sum(K))
}

.get.or.create.mhmat <- function(x) {
  fft.plan <- .get.or.create.mfft.plan(x)
  .get.or.create(x, "hmat",
                 .hmat.striped(x, fft.plan = fft.plan))
}

decompose.mssa.svd <- function(x,
                               neig = min(L, K),
                               ...,
                               force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  F <- .get(x, "F")
  h <- do.call(cbind, apply(x$F, 2, hankel, L = L))

  # Do decomposition
  S <- svd(h, nu = neig, nv = neig)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)
  if (!is.null(S$v))
    .set(x, "V", S$v)

  x
}

decompose.mssa.eigen <- function(x, ...,
                                 force.continue = FALSE) {
  N <- x$length; L <- x$window

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  F <- .get(x, "F")
  fft.plan <- .get.or.create.mfft.plan(x)
  h <- do.call(cbind, apply(x$F, 2, hankel, L = L))

  # Do decomposition
  # FIXME: Build the L-covariance matrix properly
  S <- eigen(tcrossprod(h, h), symmetric = TRUE)

  # Fix small negative values
  S$values[S$values < 0] <- 0

  # Save results
  .set(x, "lambda", sqrt(S$values))
  .set(x, "U", S$vectors)

  x
}

decompose.mssa.propack <- function(x,
                                   neig = min(50, L, prod(K)),
                                   ...,
                                   force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.")

  h <- .get.or.create.mhmat(x)
  S <- propack.svd(h, neig = neig, ...)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)
  if (!is.null(S$v))
    .set(x, "V", S$v)

  x
}

decompose.mssa.nutrlan <- function(x,
                                   neig = min(50, L, prod(K)),
                                   ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  h <- .get.or.create.mhmat(x)

  lambda <- .get(x, "lambda", allow.null = TRUE)
  U <- .get(x, "U", allow.null = TRUE)

  S <- trlan.svd(h, neig = neig, ...,
                 lambda = lambda, U = U)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)

  x
}

.init.mssa <- function(x, ...) {
  # Initialize FFT plan
  .get.or.create.mfft.plan(x)

  .get.or.create.mhmat(x)

  x
}

calc.v.mssa<- function(x, idx, env = .GlobalEnv, ...) {
  lambda <- .get(x, "lambda")[idx]
  U <- .get(x, "U")[, idx, drop = FALSE]
  h <- .get.or.create.mhmat(x)

  invisible(sapply(1:length(idx),
                   function(i) ematmul(h, U[, i], transposed = TRUE) / lambda[i]))
}

.hankelize.one.mssa <- function(x, U, V) {
  N <- x$length; L <- x$window; K <- N - L + 1

  b <- c(0, cumsum(K))
  fft.plan <- .get.or.create.mfft.plan(x)

  # FIXME: All these apply's are really ugly. Switch to C version...
  unlist(lapply(1:length(K),
                function(idx) .hankelize.one.1d.ssa(x, U, V[(b[idx]+1):b[idx+1]], fft.plan[[idx]])))
}

.slength.mssa <- function(x)
  sum(x$length)
