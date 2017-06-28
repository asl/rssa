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
  .Call("Lcor_", F, L, circular)
}

new.tmat <- function(F, L = (N + 1) %/% 2,
                     circular = FALSE,
                     fft.plan = NULL) {
  N <- length(F)
  R <- Lcor(F, L, circular = circular)

  storage.mode(R) <- "double"

  new("extmat",
      .Call("initialize_tmat", R, if (is.null(fft.plan)) fft.plan.1d(2*L - 1, L = L) else fft.plan))
}

tcols <- function(t) {
  ncol(t)
}

trows <- function(t) {
  nrow(t)
}

is.tmat <- function(t) {
  is.extmat(t) && .Call("is_tmat", t@.xData)
}

tmatmul <- function(tmat, v, transposed = FALSE) {
  ematmul(tmat, v, transposed = transposed)
}

.hankelize.one.toeplitz.ssa <- .hankelize.one.1d.ssa

.get.or.create.tmat <- function(x) {
  .get.or.create(x, "tmat", new.tmat(F = x$F, L = x$window,
                                     circular = x$circular))
}

.traj.dim.toeplitz.ssa <- .traj.dim.1d.ssa

decompose.toeplitz.ssa <- function(x,
                                   neig = NULL,
                                   ...,
                                   force.continue = FALSE) {
  ## Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0 &&
       capable(x, "decompose.continue"))
    stop(paste0("Continuation of decomposition is not yet implemented for this method: ", x$svd.method))

  if (is.null(neig))
    neig <- .default.neig(x, ...)

  if (identical(x$svd.method, "svd")) {
    S <- svd(toeplitz(Lcor(.F(x), x$window, circular = x$circular)), nu = neig, nv = neig)
    U <- S$u
    lambda <- NULL
  } else if (identical(x$svd.method, "eigen")) {
    S <- eigen(toeplitz(Lcor(.F(x), x$window, circular = x$circular)), symmetric = TRUE)
    U <- S$vectors
    lambda <- NULL
  } else if (identical(x$svd.method, "nutrlan")) {
    S <- trlan.eigen(.get.or.create.tmat(x), neig = neig, ...,
                     lambda = .decomposition(x)$lambda, U = .U(x))
    U <- S$u
    lambda <- S$d
  } else if (identical(x$svd.method, "propack")) {
    S <- propack.svd(.get.or.create.tmat(x), neig = neig, ...)
    U <- S$u
    lambda <- NULL
  } else
    stop("unsupported SVD method")

  Z <- crossprod(.get.or.create.hmat(x), U)
  sigma <- apply(Z, 2, function(x) sqrt(sum(x^2)))
  V <- sweep(Z, 2, sigma, FUN = "/")

  neig <- min(neig, length(sigma))

  o <- order(sigma[seq_len(neig)], decreasing = TRUE)
  sigma <- sigma[o]
  U <- U[, o, drop = FALSE]
  V <- V[, o, drop = FALSE]
  if (!is.null(lambda))
    lambda <- lambda[o]

  .set.decomposition(x,
                     sigma = sigma, U = U, V = V, lambda = lambda,
                     kind = "toeplitz.decomposition")

  x
}

.rowspan.toeplitz.ssa <- function(x, idx) {
  qr.Q(qr(.V(x)[, idx, drop = FALSE]))
}

.init.fragment.toeplitz.ssa <- function(this)
  expression({
    eval(.init.fragment.1d.ssa(this))
    ## Disallow shaped
    if (!all(wmask) || !all(fmask))
      stop("gaps are not allowed in Toeplitz SSA")
  })
