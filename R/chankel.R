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

.get.or.create.cfft.plan <- function(x) {
  .get.or.create(x, "fft.plan", fft.plan.1d(x$length, L = x$window, circular = x$circular,
                                            wmask = x$wmask, fmask = x$fmask, weights = x$weights))
}

.traj.dim.cssa <- function(x) {
  c(x$window, sum(x$length - x$window + 1))
}

.chmat <- function(x, fft.plan) {
  N <- x$length; L <- x$window; K <- N - L + 1
  F <- .F(x)

  R <- new.hmat(Re(F), L = L, fft.plan = fft.plan)
  I <- new.hmat(Im(F), L = L, fft.plan = fft.plan)

  matmul <- function(v)
    c(hmatmul(R, v[1:K], transposed = FALSE) - hmatmul(I, v[(K+1):(2*K)], transposed = FALSE),
      hmatmul(I, v[1:K], transposed = FALSE) + hmatmul(R, v[(K+1):(2*K)], transposed = FALSE))
  tmatmul <- function(v)
    c( hmatmul(R, v[1:L], transposed = TRUE) + hmatmul(I, v[(L+1):(2*L)], transposed = TRUE),
      -hmatmul(I, v[1:L], transposed = TRUE) + hmatmul(R, v[(L+1):(2*L)], transposed = TRUE))

  extmat(matmul, tmatmul, nrow = 2*L, ncol = 2*K)
}

.get.or.create.chmat <- function(x) {
  .get.or.create(x, "hmat",
                 .chmat(x, fft.plan = .get.or.create.cfft.plan(x)))
}

decompose.cssa <- function(x,
                           neig = min(50, L, K),
                           ...,
                           force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1
  stop("Unsupported SVD method for Complex SSA!")
}

decompose.cssa.svd <- function(x,
                               neig = min(L, K),
                               ...,
                               force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  F <- .F(x)
  h <- hankel(F, L = L)

  # Do decomposition
  S <- svd(h)

  # Save results
  .set.decomposition(x,
                     sigma = S$d[seq_len(neig)],
                     U = if (!is.null(S$u)) S$u[, seq_len(neig), drop = FALSE] else NULL,
                     V = if (!is.null(S$v)) S$v[, seq_len(neig), drop = FALSE] else NULL)

  x
}

.traj.dim.cssa.svd <- function(x) {
  c(x$window, x$length - x$window + 1)
}

cssa.to.complex <- function(values, vectors) {
  # First, make sure values come into the pairs
  d1 <- values[c(TRUE, FALSE)]
  d2 <- values[c(FALSE, TRUE)]
  if (any((d1 - d2) / d2 > 1e-3 & d2 > 0))
    warning("Too big difference between consecutive eigenvalues. CSSA might not converge")

  # And vectors
  Y1 <- vectors[1:(nrow(vectors)/2),    c(TRUE,  FALSE)]
  Z1 <- vectors[-(1:(nrow(vectors)/2)), c(TRUE,  FALSE)]
  Y2 <- vectors[1:(nrow(vectors)/2),    c(FALSE, TRUE)]
  Z2 <- vectors[-(1:(nrow(vectors)/2)), c(FALSE, TRUE)]

  V1 <- Y1 + 1i*Z1
  V2 <- Y2 + 1i*Z2

  # Sanity check
  if (any((Mod(V1 - 1i*V2) > 1e-6) & (Mod(V1 + 1i*V2) > 1e-6)))
    warning("Too big difference between consecutive eigenvectors. CSSA might not converge")

  list(d = d2, u = V2
       # , vectors2 = V2,
       # dd = (d1 - d2) / d2 > 1e-3 & d2 > 0,
       # vv = (Mod(V1 - 1i*V2) > 1e-6) & (Mod(V1 + 1i*V2) > 1e-6)
       )
}

decompose.cssa.eigen <- function(x, ...,
                                 neig = min(L, K),
                                 force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  F <- .F(x)

  R <- hankel(Re(F), L = L)
  I <- hankel(Im(F), L = L)
  h <- cbind(rbind(R, I), rbind(-I, R))

  # Do decomposition
  # FIXME: Build the L-covariance matrix properly
  S <- eigen(tcrossprod(h, h), symmetric = TRUE)

  # Fix small negative values
  S$values[S$values < 0] <- 0

  S <- cssa.to.complex(sqrt(S$values), S$vectors)

  # Save results
  .set.decomposition(x,
                     sigma = S$d[1:neig], U = S$u[, 1:neig, drop = FALSE])

  x
}

decompose.cssa.propack <- function(x,
                                   neig = min(50, L, K),
                                   ...,
                                   force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.")

  h <- .get.or.create.chmat(x)
  S <- propack.svd(h, neig = 2*neig, ...)

  S <- cssa.to.complex(S$d, S$u)

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)

  x
}

decompose.cssa.nutrlan <- function(x,
                                   neig = min(50, L, K),
                                   ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  h <- .get.or.create.chmat(x)

  sigma <- .sigma(x)
  U <- .U(x)

  S <- trlan.svd(h, neig = 2*neig, ...,
                 lambda = sigma, U = U)

  S <- cssa.to.complex(S$d, S$u)

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u)

  x
}

.traj.dim.cssa <- function(x) {
  c(2*x$window, 2*(x$length - x$window + 1))
}

calc.v.cssa<- function(x, idx, env = .GlobalEnv, ...) {
  sigma <- .sigma(x)[idx]

  if (any(sigma <= .Machine$double.eps)) {
    sigma[sigma <= .Machine$double.eps] <- Inf
    warning("some sigmas are equal to zero. The corresponding vectors will be zeroed")
  }

  U <- .U(x)[, idx, drop = FALSE]
  h <- .get.or.create.chmat(x)

  invisible(sapply(1:length(idx),
                   function(i) {
                     v <- ematmul(h, c(Re(U[, i]), Im(U[, i])), transposed = TRUE) / sigma[i]
                     v[1:(length(v) / 2)] + 1i*v[-(1:(length(v) / 2))]
                   }))
}

.hankelize.one.cssa <- function(x, U, V, fft.plan = .get.or.create.cfft.plan(x)) {
  R1 <- .hankelize.one.default(Re(U), Re(V), fft.plan = fft.plan)
  R2 <- .hankelize.one.default(Im(U), Im(V), fft.plan = fft.plan)
  I1 <- .hankelize.one.default(Re(U), Im(V), fft.plan = fft.plan)
  I2 <- .hankelize.one.default(Im(U), Re(V), fft.plan = fft.plan)

  (R1 + R2) + 1i*(-I1 + I2)
}

.hankelize.multi.complex <- function(U, V, fft.plan) {
  ReU <- Re(U); ReV <- Re(V); ImU <- Im(U); ImV <- Im(V)
  storage.mode(ReU) <- storage.mode(ReV) <- storage.mode(ImU) <- storage.mode(ImV) <- "double"
  R1 <- .Call("hankelize_multi_fft", ReU, ReV, fft.plan)
  R2 <- .Call("hankelize_multi_fft", ImU, ImV, fft.plan)
  I1 <- .Call("hankelize_multi_fft", ReU, ImV, fft.plan)
  I2 <- .Call("hankelize_multi_fft", ImU, ReV, fft.plan)

  (R1 + R2) + 1i*(-I1 + I2)
}

plot.cssa.reconstruction <- function(x,
                                     slice = list(),
                                     ...,
                                     type = c("raw", "cumsum"),
                                     plot.method = c("native", "matplot"),
                                     na.pad = c("left", "right"),
                                     base.series = NULL,
                                     add.original = TRUE,
                                     add.residuals = TRUE) {
  # Adopt CSSA to MSSA case - construct new reconstruction object
  original <- attr(x, "series")
  res <- attr(x, "residuals")

  x <- lapply(x, function(el) list(Re = Re(el), Im = Im(el)))
  attr(x, "series") <- list(Re = Re(original), Im = Im(original))
  attr(x, "residuals") <- list(Re = Re(res), Im = Im(res))

  class(x) <- paste(c("mssa", "ssa"), "reconstruction", sep = ".")

  # Do the call with the same set of arguments
  mplot <- match.call(expand.dots = TRUE)
  mplot[[1L]] <- as.name("plot")
  mplot[[2L]] <- x
  eval(mplot, parent.frame())
}
