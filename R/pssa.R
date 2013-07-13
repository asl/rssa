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
    return(matrix(NA_real_, L, 0))
  }

  mx <- sapply(0:degree, function(deg) (1:L) ^ deg)

  # G-S orthogonalization and return
  qr.Q(qr(mx))
}

.phmat <- function(x) {
  hmat <- .get.or.create.hmat(x)
  left.projector <- .get(x, "left.projector")
  right.projector <- .get(x, "right.projector")

  matmul <- function(v) {
    v <- v - right.projector %*% crossprod(right.projector, v)
    v <- hmatmul(hmat, v, transposed = FALSE)
    v <- v - left.projector %*% crossprod(left.projector, v)

    v
  }

  tmatmul <- function(v) {
    v <- v - left.projector %*% crossprod(left.projector, v)
    v <- hmatmul(hmat, v, transposed = TRUE)
    v <- v - right.projector %*% crossprod(right.projector, v)

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
  left.projector <- .get(x, "left.projector")
  right.projector <- .get(x, "right.projector")

  RF <- hmatmul.wrap(hmat, right.projector, transposed = FALSE)
  LF <- hmatmul.wrap(hmat, left.projector, transposed = TRUE) - right.projector %*% crossprod(RF, left.projector)
  R.lambda <- sqrt(colSums(RF^2))
  RF <- RF / rep(R.lambda, each = nrow(RF))
  L.lambda <- sqrt(colSums(LF^2))
  LF <- LF / rep(L.lambda, each = nrow(LF))

  right.proj.triples <- list(lambda = R.lambda, U = RF, V = right.projector)
  left.proj.triples <- list(lambda = L.lambda, U = left.projector, V = LF)

  # Store special components
  .set(x, "right.proj.triples", right.proj.triples)
  .set(x, "left.proj.triples", left.proj.triples)
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
