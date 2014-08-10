#   R package for Singular Spectrum Analysis
#   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
#   Copyright (c) 2009 Konstantin Usevich <usevich.k.d@gmail.com>
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

#   Routines for hankel-block hankel (aka 2d) SSA


fft2 <- function(X, inverse = FALSE) {
  # TODO Use FTTW here
  t(mvfft(t(mvfft(X, inverse = inverse)), inverse = inverse))
}

.convolution.dims <- function(x.dim, y.dim, type = "circular") {
  type <- sapply(type, match.arg, choices = c("circular", "open", "filter"))

  stopifnot(length(x.dim) == length(y.dim))
  rank <- length(x.dim)

  if (length(type) != rank) {
    # Use recycling

    if (rank %% length(type) != 0)
      warning("longer object length is not a multiple of shorter object length")

    type <- type[(seq_len(rank) - 1) %% length(type) + 1]
  }

  stopifnot(x.dim[type == "filter"] >= y.dim[type == "filter"])

  output.dim <- ifelse(type == "circular",
                       x.dim,
                       ifelse(type == "open",
                              x.dim + y.dim - 1,
                              x.dim - y.dim + 1))

  input.dim <- ifelse(type == "open",
                      x.dim + y.dim - 1,
                      x.dim)

  names(input.dim) <- names(output.dim) <- names(x.dim)

  list(input.dim = input.dim,
       output.dim = output.dim)
}

.convolven <- function(x, y, conj = TRUE, type = "circular") {
  if (is.null(dim(x))) dim(x) <- length(x)
  if (is.null(dim(y))) dim(y) <- length(y)

  io.dims <- .convolution.dims(dim(x), dim(y), type)
  input.dim <- io.dims$input.dim
  output.dim <- io.dims$output.dim

  storage.mode(x) <- storage.mode(y) <- "double"
  storage.mode(conj) <- "logical"
  storage.mode(input.dim) <- storage.mode(output.dim) <- "integer"

  res <- .Call("convolveN", x, y, input.dim, output.dim, conj)
  if (length(output.dim) == 1)
    dim(res) <- NULL
  else
    dim(res) <- output.dim

  res
}

.factor.mask.2d <- function(field.mask, window.mask, circular = FALSE) {
  field.mask[] <- as.numeric(field.mask)
  window.mask[] <- as.numeric(window.mask)
  tmp <- .convolven(field.mask, window.mask, conj = TRUE,
                   type = ifelse(circular, "circular", "filter"))

  abs(tmp - sum(window.mask)) < 0.5 # ==0, but not exact in case of numeric error
}

.field.weights.2d <- function(window.mask, factor.mask, circular = FALSE) {
  window.mask[] <- as.numeric(window.mask)
  factor.mask[] <- as.numeric(factor.mask)
  res <- .convolven(factor.mask, window.mask, conj = FALSE,
                   type = ifelse(circular, "circular", "open"))
  res[] <- as.integer(round(res))

  res
}

# TODO generalize to ndim case
circle.mask <- function(R) {
  I <- matrix(seq_len(2*R - 1), 2*R - 1, 2*R - 1)
  J <- t(I)

  (I - R)^2 + (J - R)^2 < R^2
}

triangle.mask <- function(side) {
  I <- matrix(seq_len(side), side, side)
  J <- t(I)

  I + J <= side + 1
}

new.hbhmat <- function(F, L = (N + 1) %/% 2,
                       wmask = NULL, fmask = NULL, weights = NULL,
                       circular = FALSE) {
  rank <- length(dim(F))
  if (length(circular) > rank)
    warning("Incorrect argument length: length(circular) > rank, two leading values will be used")
  if (length(circular) != rank)
    circular <- circular[(seq_len(rank) - 1) %% length(circular) + 1]

  N <- dim(F)

  storage.mode(F) <- "double"
  storage.mode(L) <- "integer"
  storage.mode(circular) <- "logical"

  if (!is.null(wmask)) {
    storage.mode(wmask) <- "logical"
  }

  if (!is.null(fmask)) {
    storage.mode(fmask) <- "logical"
  }

  if (!is.null(weights)) {
    mask <- weights > 0
    F[!mask] <- mean(F[mask]) # Improve FFT stability & remove NAs
  } else {
    weights <- .hweightsn(N, L)
  }
  storage.mode(weights) <- "integer"

  h <- .Call("initialize_hbhmat", F, L, wmask, fmask, weights, circular)
}

hbhcols <- function(h) {
  .Call("hbhankel_cols", h)
}

hbhrows <- function(h) {
  .Call("hbhankel_rows", h)
}

is.hbhmat <- function(h) {
  .Call("is_hbhmat", h)
}

hbhmatmul <- function(hmat, v, transposed = FALSE) {
  storage.mode(v) <- "double"
  storage.mode(transposed) <- "logical"
  .Call("hbhmatmul", hmat, v, transposed)
}

.get.or.create.hbhmat <- function(x) {
  .get.or.create(x, "hmat",
                 new.hbhmat(.F(x), L = x$window, wmask = x$wmask, fmask = x$fmask,
                            weights = x$weights, circular = x$circular))
}

as.matrix.hbhmat <- function(x) {
  apply(diag(hbhcols(x)), 2, hbhmatmul, hmat = x)
}

tcrossprod.hbhmat <- function(x) {
  apply(diag(hbhrows(x)), 2, function(u) hbhmatmul(hmat = x, hbhmatmul(hmat = x, u, transposed = TRUE)))
}

.traj.dim.nd.ssa <- function(x) {
  Ldim <- sum(x$wmask)
  if (Ldim == 0)
    Ldim <- prod(x$window)

  Kdim <- sum(x$fmask)
  if (Kdim == 0)
    Kdim <- prod(x$length - ifelse(x$circular, 0, x$window - 1))

  c(Ldim, Kdim)
}

decompose.nd.ssa <- function(x,
                             neig = min(50, prod(L), prod(K)),
                             ...) {
  N <- x$length; L <- x$window; K <- N - L + 1
  stop("Unsupported SVD method for 2D.SSA!")
}

decompose.nd.ssa.svd <- function(x,
                                 neig = min(50, prod(L), prod(K)),
                                 ...,
                                 force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not yet implemented for this method.")

  # Create circulant and convert it to ordinary matrix
  h <- as.matrix.hbhmat(.get.or.create.hbhmat(x))

  # Do decompostion
  S <- svd(h, nu = neig, nv = neig)

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)

  x
}

decompose.nd.ssa.eigen <- function(x,
                                   neig = min(50, prod(L), prod(K)),
                                   ...,
                                   force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not yet implemented for this method.")

  # Create circulant and compute XX^T in form of ordinary matrix
  C <- tcrossprod.hbhmat(.get.or.create.hbhmat(x))

  # Do decompostion
  S <- eigen(C, symmetric = TRUE)

  # Fix small negative values
  S$values[S$values < 0] <- 0

  # Save results
  .set.decomposition(x,
                     sigma = sqrt(S$values[1:neig]),
                     U = S$vectors[, 1:neig, drop = FALSE])

  x
}

decompose.nd.ssa.nutrlan <- function(x,
                                     neig = min(50, prod(L), prod(K)),
                                     ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  h <- .get.or.create.hbhmat(x)

  sigma <- .sigma(x)
  U <- .U(x)

  S <- trlan.svd(h, neig = neig, ...,
                 lambda = sigma, U = U)

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u)

  x
}

decompose.nd.ssa.propack <- function(x,
                                     neig = min(50, prod(L), prod(K)),
                                     ...,
                                     force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not yet implemented for this method.")

  h <- .get.or.create.hbhmat(x)

  S <- propack.svd(h, neig = neig, ...)

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)

  x
}

calc.v.nd.ssa <- function(x, idx, ...) {
  sigma <- .sigma(x)[idx]

  if (any(sigma <= .Machine$double.eps)) {
    sigma[sigma <= .Machine$double.eps] <- Inf
    warning("some sigmas are equal to zero. The corresponding vectors will be zeroed")
  }

  U <- .U(x)[, idx, drop = FALSE]
  h <- .get.or.create.hbhmat(x)

  invisible(sapply(1:length(idx),
                   function(i) hbhmatmul(h, U[, i], transposed = TRUE) / sigma[i]))
}

.hankelize.one.nd.ssa <- function(x, U, V) {
  h <- .get.or.create.hbhmat(x)
  storage.mode(U) <- storage.mode(V) <- "double"
  .Call("hbhankelize_one_fft", U, V, h)
}
