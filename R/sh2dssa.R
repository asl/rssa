#   R package for Singular Spectrum Analysis
#   Copyright (c) 2013 Alexander Shlemov <shlemovalex@gmail.com>
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

fft2 <- function(X, inverse = FALSE) {
  # TODO Use FTTW here
  t(mvfft(t(mvfft(X, inverse = inverse)), inverse = inverse))
}

convolve2.open <- function(X, Y, conj = FALSE) {
  new.dim <- dim(X) + dim(Y) - 1

  x <- matrix(0, new.dim[1], new.dim[2])
  x[1:nrow(X), 1:ncol(X)] <- X

  y <- matrix(0, new.dim[1], new.dim[2])
  y[1:nrow(Y), 1:ncol(Y)] <- Y

  Re(fft2(fft2(x) * if (conj) Conj(fft2(y)) else fft2(y), inverse = TRUE) / prod(new.dim))
}

convolve2.filter <- function(X, Y, conj = TRUE) {
  new.dim <- dim(X) - dim(Y) + 1

  x <- X
  y <- matrix(0, nrow(X), ncol(X))
  y[1:nrow(Y), 1:ncol(Y)] <- Y

  tmp <- Re(fft2(fft2(x) * if (conj) Conj(fft2(y)) else fft2(y), inverse = TRUE) / prod(dim(X)))
  tmp[seq_len(new.dim[1]), seq_len(new.dim[2])]
}

factor.mask <- function(field.mask, window.mask) {
  field.mask[] <- as.numeric(field.mask)
  window.mask[] <- as.numeric(window.mask)
  tmp <- convolve2.filter(field.mask, window.mask)

  abs(tmp - sum(window.mask)) < 0.5 # ==0, but not exact in case of numeric error
}

field.weights <- function(window.mask, factor.mask) {
  window.mask[] <- as.numeric(window.mask)
  factor.mask[] <- as.numeric(factor.mask)
  res <- convolve2.open(factor.mask, window.mask)
  res[] <- as.integer(round(res))

  res
}

circle.mask <- function(R) {
  I <- matrix(seq_len(2*R - 1), 2*R - 1, 2*R - 1)
  J <- t(I)

  (I - R)^2 + (J - R)^2 < R^2
}

fmatmul <- function(field, X, umask, vmask, transposed = FALSE) {
  if (transposed) {
    tmp <- umask; umask <- vmask; vmask <- tmp
  }

  stopifnot(length(X) == sum(vmask))

  x <- vmask
  x[vmask] <- X
  y <- convolve2.filter(field, x)
  Y <- y[umask]

  Y
}

.sh2hbhmat <- function(x) {
  umask <- .get(x, "umask")
  vmask <- .get(x, "vmask")
  weights <- .get(x, "weights")

  F <- .get(x, "F")
  mask <- weights > 0
  F[!mask] <- mean(F[mask]) # Improve FFT stability & remove NAs

  matmul <- function(v) fmatmul(F, v, umask, vmask)
  tmatmul <- function(u) fmatmul(F, u, umask, vmask, transposed = TRUE)

  extmat(matmul, tmatmul, sum(umask), sum(vmask))
}

.get.or.create.sh2hbhmat <- function(x) {
  .get.or.create(x, "hmat", .sh2hbhmat(x))
}

decompose.shaped2d.ssa <- function(x,
                                   neig = min(50, prod(L), prod(K)),
                                   ...) {
  N <- x$length; L <- x$window; K <- N - L + 1
  stop("Unsupported SVD method for shaped 2D.SSA!")
}

decompose.shaped2d.ssa.propack <- function(x,
                                       neig = min(50, prod(L), prod(K)),
                                       ...,
                                       force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.")

  sh <- .get.or.create.sh2hbhmat(x)
  S <- propack.svd(sh, neig = neig, ...)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)
  if (!is.null(S$v))
    .set(x, "V", S$v)

  x
}

decompose.shaped2d.ssa.nutrlan <- function(x,
                                       neig = min(50, prod(L), prod(K)),
                                       ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  sh <- .get.or.create.sh2hbhmat(x)

  lambda <- .get(x, "lambda", allow.null = TRUE)
  U <- .get(x, "U", allow.null = TRUE)

  S <- trlan.svd(sh, neig = neig, ...,
                 lambda = lambda, U = U)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)

  x
}

calc.v.shaped2d.ssa <- function(x, idx, ...) {
  lambda <- .get(x, "lambda")[idx]
  U <- .get(x, "U")[, idx, drop = FALSE]
  sh <- .get.or.create.sh2hbhmat(x)

  invisible(sapply(1:length(idx),
                   function(i) ematmul(sh, U[, i], transposed = TRUE) / lambda[i]))
}

.hankelize.one.shaped2d.ssa <- function(x, U, V) {
  umask <- .get(x, "umask")
  vmask <- .get(x, "vmask")
  w <- .get(x, "weights")

  x <- umask
  x[umask] <- U

  y <- vmask
  y[vmask] <- V

  res <- convolve2.open(x, y)
  res <- res / w
  res[w == 0] <- NA

  res
}

