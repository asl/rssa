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

.get.or.create.trajmat.cssa <- .get.or.create.chmat

decompose.cssa <- function(x,
                           neig = NULL,
                           ...,
                           force.continue = FALSE) {
  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  if (is.null(neig))
    neig <- .default.neig(x, ...)


  if (identical(x$svd.method, "svd")) {
    S <- svd(hankel(.F(x), L = x$window), nu = neig, nv = neig)
    .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else if (identical(x$svd.method, "eigen")) {
    h <- hankel(.F(x), L = x$window)

    ## FIXME: Build the complex L-covariance matrix properly
    S <- eigen(tcrossprod(h, Conj(h)), symmetric = TRUE)

    ## Fix small negative values
    S$values[S$values < 0] <- 0

    ## Save results
    .set.decomposition(x,
                       sigma = S$values[seq_len(neig)],
                       U = S$vectors[, seq_len(neig), drop = FALSE])
  } else
    stop("unsupported SVD method")

  x
}

.traj.dim.cssa <- function(x) {
  c(x$window, x$length - x$window + 1)
}

calc.v.cssa <- function(x, idx, ...) {
  nV <- nv(x)

  V <- matrix(NA_complex_, .traj.dim(x)[2], length(idx))
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

    h <- .get.or.create.chmat(x)
    rV <- crossprod(h, rbind(Re(U), Im(U))) / rep(sigma, each = 2*nrow(V))
    V[, idx > nV] <- rV[seq_len(nrow(V)),, drop = FALSE] + 1i*rV[-seq_len(nrow(V)),, drop = FALSE]
  }

  invisible(V)
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

.init.fragment.cssa <- function(this)
  expression({
    if (any(circular))
      stop("Circular variant of complex SSA isn't implemented yet")

    # Sanity check - the input series should be complex
    if (!is.complex(x))
      stop("complex SSA should be performed on complex time series")
    N <- length(x)

    wmask <- fmask <- weights <- NULL

    column.projector <- row.projector <- NULL
  })
