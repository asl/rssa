#   R package for Singular Spectrum Analysis
#   Copyright (c) 2008, 2009 Anton Korobeynikov <asl@math.spbu.ru>
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

wcor.default <- function(x, L = (N + 1) %/% 2, ..., weights = NULL) {
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")

  if (is.null(weights)) {
    # Compute weights
    N <- nrow(x)
    weights <- .hweights(x, L)
  }

  # Compute w-covariation
  cov <- crossprod(weights * x, x)

  # Convert to correlations
  cor <- cov2cor(cov)

  # Fix possible numeric error
  cor[cor > 1] <- 1; cor[cor < -1] <- -1

  # Add class
  class(cor) <- "wcor.matrix"

  # Return
  cor
}

wcor.toeplitz.ssa <- wcor.1d.ssa <- function(x, groups, ..., cache = TRUE) {
  N <- prod(x$length)
  if (missing(groups))
    groups <- as.list(1:nlambda(x))

  # Compute reconstruction.
  F <- reconstruct(x, groups, ..., cache = cache)
  mx <- matrix(unlist(F), nrow = N, ncol = length(groups))
  colnames(mx) <- names(F)

  # Finally, compute w-correlations and return
  wcor.default(mx, weights = .hweights(x))
}

wcor.2d.ssa <- function(x, groups, ..., cache = TRUE) {
  N <- prod(x$length)
  if (missing(groups))
    groups <- as.list(1:nlambda(x))

  # Compute reconstruction.
  F <- reconstruct(x, groups, ..., cache = cache)
  mx <- matrix(unlist(F), nrow = N, ncol = length(groups))
  colnames(mx) <- names(F)

  # Get weights
  w <- .hweights(x)

  if (any(w == 0)) {
    # Omit uncovered elements
    mx <- mx[as.vector(w > 0),, drop = FALSE]
    w <- as.vector(w[w > 0])
  }

  # Finally, compute w-correlations and return
  wcor.default(mx, weights = w)
}

wcor.mssa <- function(x, groups, ..., cache = TRUE) {
  N <- sum(x$length)
  if (missing(groups))
    groups <- as.list(1:nlambda(x))

  # Compute reconstruction.
  F <- lapply(reconstruct(x, groups, ..., cache = cache), .to.series.list)
  mx <- matrix(unlist(F), nrow = N, ncol = length(groups))
  colnames(mx) <- names(F)

  # Finally, compute w-correlations and return
  wcor.default(mx, weights = .hweights(x))
}

wcor.ssa <- function(x, groups, ..., cache = TRUE)
  stop("Unsupported SVD method for SSA!")

wcor <- function(x, ...) {
  UseMethod("wcor")
}

clusterify.wcor.matrix <- function(x,
                                   nclust = N,
                                   ...,
                                   dist = function(X) (1 - X) / 2) {
  N <- nrow(x)
  h <- cutree(hclust(as.dist(dist(x)), ...), k = nclust)
  split(1:N, h)
}

.hweights <- function(x, ...) {
  UseMethod(".hweights")
}

.hweights.default <- function(x, L = (N + 1) %/% 2, ...) {
  N <- if (length(x) == 1) x else length(x)
  K <- N - L + 1
  Ls <- min(L, K); Ks <- max(L, K)

  # Compute and return weights
  if (Ls > 1)
    c(1:(Ls-1), rep(Ls, Ks-Ls+1), (Ls-1):1)
  else
    rep(1, N)
}

.hweights.matrix <- function(x, L = (N + 1) %/% 2, ...) {
  N <- nrow(x)

  .hweights.default(N, L)
}

.hweights.1d.ssa <- .hweights.toeplitz.ssa <- function(x, ...) {
  .hweights.default(x$length, x$window)
}

.hweights.2d.ssa <- function(x, ...) {
  w <- .get(x, "weights")

  if (!is.null(w)) {
    # Just return stored weights
    w
  } else {
    N <- x$length; L <- x$window

    as.vector(tcrossprod(.hweights.default(N[1], L[1]),
                         .hweights.default(N[2], L[2])))
  }
}

.hweights.mssa <- function(x, ...) {
  w <- .get(x, "weights")

  if (!is.null(w)) {
    # Return positive weights
    w[w > 0]
  } else {
    N <- x$length; L <- x$window

    rep(.hweights.default(N[1], L), length(N))
  }
}

wnorm.default <- function(x, L = (N + 1) %/% 2, ...) {
  N <- length(x)

  # Compute weights
  w <- .hweights.default(x, L)

  # Compute wnorm
  sqrt(sum(w * x^2))
}

wnorm.1d.ssa <- wnorm.toeplitz.ssa <- function(x, ...) {
  # Compute weights
  w <- .hweights(x)

  # Compute wnorm
  sqrt(sum(w * as.vector(x$F)^2))
}

wnorm.2d.ssa <- function(x, ...) {
  # Get F
  F <- .get(x, "F")

  # Compute weights
  w <- .hweights(x)

  if (any(w == 0)) {
    # Omit uncovered elements
    F <- as.vector(F[w > 0])
    w <- as.vector(w[w > 0])
  }

  # Compute wnorm
  sqrt(sum(w * F^2))
}

wnorm.mssa <- function(x, ...) {
  # Compute weights
  w <- .hweights(x)

  # Get series
  F <- .get(x, "F")

  # Compute wnorm
  sqrt(sum(w * unlist(F)^2))
}

#N = 399;
#a = 1.005;
#T = 200;
#F1 <- (1/a)^(1:N);
#F2 <- (a^(1:N))*cos(2*pi*(1:N)/T);
#.wcor(cbind(F1, F2));
