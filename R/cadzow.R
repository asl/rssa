#   R package for Singular Spectrum Analysis
#   Copyright (c) 2013 Anton Korobeynikov <asl@math.spbu.ru>
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

.extend.series <- function(x, alpha) {
  if (is.list(x)) lapply(x, sys.function(), alpha = alpha) else alpha * x
}

.series.dist <- function(F1, F2, norm, mask = TRUE) {
  mask <- as.logical(mask)
  F1 <- as.vector(unlist(F1))[mask]
  F2 <- as.vector(unlist(F2))[mask]

  norm(F1 - F2)
}

.series.winnerprod <- function(F1, F2, weights = 1) {
  mask <- weights > 0

  weights <- weights[mask]
  F1 <- as.vector(unlist(F1))[mask]
  F2 <- as.vector(unlist(F2))[mask]

  sum(weights * F1 * F2)
}

.inner.fmt.conversion <- function(x, ...)
  UseMethod(".inner.fmt.conversion")

.inner.fmt.conversion.ssa <- function(x, ...)
  identical

.inner.fmt.conversion.1d.ssa <- .inner.fmt.conversion.toeplitz.ssa <- function(x, ...)
  as.numeric

.inner.fmt.conversion.cssa <- function(x, ...)
  as.complex

.inner.fmt.conversion.2d.ssa <- function(x, ...)
  as.matrix

.inner.fmt.conversion.nd.ssa <- function(x, ...)
  as.array

.inner.fmt.conversion.mssa <- function(x, ...) {
  template <- x$F

  # Prevent storing huge ssa-object in closure
  x <- NULL

  function(x) .to.series.list(x, template = template)
}

cadzow.ssa <- function(x, rank,
                       correct = TRUE,
                       tol = 1e-6, maxiter = 0,
                       norm = function(x) sqrt(max(x^2)),
                       trace = FALSE,
                       ..., cache = TRUE) {
  # Get conversion
  conversion <- .inner.fmt.conversion(x)

  # Get weights and mask
  weights <- .hweights(x)
  mask <- weights > 0

  # Obtain the initial reconstruction of rank r
  r <- reconstruct(x, groups = list(1:rank), ..., cache = cache)
  stopifnot(length(r) == 1)
  F <- r[[1]]

  # Do the actual iterations until the convergence (or stoppping due to number
  # of iterations)
  it <- 0
  repeat {
    s <- clone(x, copy.cache = FALSE, copy.storage = FALSE)
    .set(s, "F", conversion(F))
    r <- reconstruct(s, groups = list(1:rank), ..., cache = FALSE)
    stopifnot(length(r) == 1)
    rF <- r[[1]]

    it <- it + 1
    if ((maxiter > 0 && it >= maxiter) || (sqd <- .series.dist(F, rF, norm, mask)) < tol)
      break
    if (trace)
      cat(sprintf("Iteration: %d, distance: %s\n", it, format(sqd)))
    F <- rF
  }

  if (correct) {
    alpha <- .series.winnerprod(.F(x), F, weights) / .series.winnerprod(F, F, weights)
    F <- .extend.series(F, alpha)
  }

  F
}

cadzow <- function(x, ...)
  UseMethod("cadzow")
