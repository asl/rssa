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

new.hbhmat <- function(F,
                       L = (N + 1) %/% 2) {
  N <- dim(F)
  storage.mode(F) <- "double"
  storage.mode(L) <- "integer"
  h <- .Call("initialize_hbhmat", F, L[1], L[2])
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
                 new.hbhmat(x$F, L = x$window))
}

as.matrix.hbhmat <- function(x) {
  apply(diag(hbhcols(x)), 2, hbhmatmul, hmat = x)
}

tcrossprod.hbhmat <- function(x) {
  apply(diag(hbhrows(x)), 2, function(u) hbhmatmul(hmat = x, hbhmatmul(hmat = x, u, transposed = TRUE)))
}

decompose.2d.ssa <- function(x,
                             neig = min(50, prod(L), prod(K)),
                             ...) {
  N <- x$length; L <- x$window; K <- N - L + 1
  stop("Unsupported SVD method for 2D.SSA!")
}

decompose.2d.ssa.svd <- function(x,
                                 neig = min(50, prod(L), prod(K)),
                                 ...,
                                 force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decomposition is not yet implemented for this method.")

  # Create circulant and convert it to ordinary matrix
  h <- as.matrix.hbhmat(.get.or.create.hbhmat(x))

  # Do decompostion
  S <- svd(h, nu = neig, nv = neig)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)
  if (!is.null(S$v))
    .set(x, "V", S$v)

  x
}

decompose.2d.ssa.eigen <- function(x,
                                   neig = min(50, prod(L), prod(K)),
                                   ...,
                                   force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decomposition is not yet implemented for this method.")

  # Create circulant and compute XX^T in form of ordinary matrix
  C <- tcrossprod.hbhmat(.get.or.create.hbhmat(x))

  # Do decompostion
  S <- eigen(C, symmetric = TRUE)

  # Fix small negative values
  S$values[S$values < 0] <- 0

  # Save results
  .set(x, "lambda", sqrt(S$values[1:neig]))
  .set(x, "U", S$vectors[, 1:neig, drop = FALSE])

  x
}

decompose.2d.ssa.nutrlan <- function(x,
                                     neig = min(50, prod(L), prod(K)),
                                     ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  h <- .get.or.create.hbhmat(x)

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

decompose.2d.ssa.propack <- function(x,
                                     neig = min(50, prod(L), prod(K)),
                                     ...,
                                     force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decomposition is not yet implemented for this method.")

  h <- .get.or.create.hbhmat(x)

  S <- propack.svd(h, neig = neig, ...)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)
  if (!is.null(S$v))
    .set(x, "V", S$v)

  x
}

calc.v.2d.ssa <- function(x, idx, ...) {
  lambda <- .get(x, "lambda")[idx]
  U <- .get(x, "U")[, idx, drop = FALSE]
  h <- .get.or.create.hbhmat(x)

  invisible(sapply(1:length(idx),
                   function(i) hbhmatmul(h, U[, i], transposed = TRUE) / lambda[i]))
}

.hankelize.one.2d.ssa <- function(x, U, V) {
  h <- .get.or.create.hbhmat(x)
  storage.mode(U) <- storage.mode(V) <- "double"
  .Call("hbhankelize_one_fft", U, V, h)
}
