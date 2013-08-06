#   R package for Singular Spectrum Analysis
#   Copyright (c) 2013 Alex Shlemov <shlemovalex@gmail.com>
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

orthopoly <- function(degree, L) {
  if (is.character(degree)) {
    degree <- match.arg(degree,
                        choices = c("none", "constant", "linear", "quadratic", "qubic"))
    degree <- switch(degree,
                     none = -1, constant = 0, linear = 1, quadratic = 2, qubic = 3)
  }

  stopifnot(is.numeric(degree) && length(degree) == 1)

  if (degree == -1) {
    # Return matrix with zero columns
    matrix(NA_real_, L, 0)
  } else if (degree == 0) {
    matrix(1 / sqrt(L), L, 1)
  } else {
    cbind(1 / sqrt(L), poly(1:L, degree))
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

.get.or.create.hmat <- function(x) {
  .get.or.create(x, "hmat",
                 new.hmat(x$F, L = x$window,
                          fft.plan = .get.or.create.fft.plan(x)))
}

.compute.proj.triples <- function(x) {
  hmat <- .get.or.create.hmat(x)
  column.projector <- .get(x, "column.projector")
  row.projector <- .get(x, "row.projector")

  RF <- hmatmul.wrap(hmat, row.projector, transposed = FALSE)
  LF <- hmatmul.wrap(hmat, column.projector, transposed = TRUE) - row.projector %*% crossprod(RF, column.projector)
  R.lambda <- sqrt(colSums(RF^2))
  RF <- RF / rep(R.lambda, each = nrow(RF))
  L.lambda <- sqrt(colSums(LF^2))
  LF <- LF / rep(L.lambda, each = nrow(LF))

  row.proj.triples <- list(lambda = R.lambda, U = RF, V = row.projector)
  column.proj.triples <- list(lambda = L.lambda, U = column.projector, V = LF)

  # Store special components
  .set(x, "row.proj.triples", row.proj.triples)
  .set(x, "column.proj.triples", column.proj.triples)
}

decompose.pssa <- function(x,
                           neig = min(50, L, K),
                           ...,
                           force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1
  stop("Unsupported SVD method for SSA with projection!")
}

decompose.pssa.propack <- function(x,
                                   neig = min(50, L, K),
                                   ...,
                                   force.continue = FALSE) {
  .compute.proj.triples(x)
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.")

  ph <- .get.or.create.phmat(x)
  S <- propack.svd(ph, neig = neig, ...)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)
  if (!is.null(S$v))
    .set(x, "V", S$v)

  x
}

decompose.pssa.nutrlan <- function(x,
                                   neig = min(50, L, K),
                                   ...) {
  .compute.proj.triples(x)
  N <- x$length; L <- x$window; K <- N - L + 1

  ph <- .get.or.create.phmat(x)

  lambda <- .get(x, "lambda", allow.null = TRUE)
  U <- .get(x, "U", allow.null = TRUE)

  S <- trlan.svd(ph, neig = neig, ...,
                 lambda = lambda, U = U)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)

  x
}

calc.v.pssa <- function(x, idx, ...) {
  lambda <- .get(x, "lambda")[idx]
  U <- .get(x, "U")[, idx, drop = FALSE]
  ph <- .get.or.create.phmat(x)

  invisible(sapply(1:length(idx),
                   function(i) ematmul(ph, U[, i], transposed = TRUE) / lambda[i]))
}

.hankelize.one.pssa <- .hankelize.one.1d.ssa
