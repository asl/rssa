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
  cov <- crossprod(weights * Conj(x), x)

  # Convert to correlations
  Is <- sqrt(1 / abs(diag(cov)))
  cor <- cov
  cor[] <- Is * cov * rep(Is, each = nrow(cov))
  cor[cbind(seq_len(nrow(cov)), seq_len(ncol(cov)))] <- 1

  # Fix possible numeric error
  cor[cor > 1] <- 1; cor[cor < -1] <- -1

  # Add class
  class(cor) <- "wcor.matrix"

  # Return
  cor
}

wcor.ssa <- function(x, groups, Fs, ..., cache = TRUE) {
  # Get conversion
  conversion <- .inner.fmt.conversion(x)

  if (!missing(Fs) && !missing(groups)) {
    stop("Only one of `groups' and `Fs' shoud be passed")
  }

  if (missing(Fs)) {
    if (missing(groups))
      groups <- as.list(1:nsigma(x))

    # Compute reconstruction.
    Fs <- reconstruct(x, groups = groups, ..., cache = cache)
  }

  Fs <- lapply(Fs, function(x) as.vector(unlist(conversion(x))))
  mx <- do.call(cbind, Fs)

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

wcor <- function(x, ...) {
  UseMethod("wcor")
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

.hweightsn <- function(N, L) {
  stopifnot(length(N) == length(L))
  ws <- lapply(seq_along(N),
               function(r) .hweights.default(N[r], L[r]))

  if (length(N) > 1)
    for (r in 2:length(N))
      ws[[1]] <- as.vector(tcrossprod(ws[[1]], ws[[r]]))

  ws[[1]]
}

.hweights.matrix <- function(x, L = (N + 1) %/% 2, ...) {
  N <- nrow(x)

  .hweights.default(N, L)
}

.hweights.1d.ssa <- .hweights.toeplitz.ssa <- .hweights.cssa <- .hweights.nd.ssa <- function(x, ...) {
  w <- .get(x, "weights")

  if (!is.null(w)) {
    # Just return stored weights
    w
  } else {
    .hweightsn(x$length, x$window)
  }
}

.hweights.mssa <- function(x, ...) {
  w <- .get(x, "weights")

  if (!is.null(w)) {
    ## Return meaningfull weights
    N <- x$length; mN <- max(N)
    cidx <- unlist(lapply(seq_along(N), function(idx) seq_len(N[idx]) + mN * (idx - 1)))
    w[cidx]
  } else {
    N <- x$length; L <- x$window

    rep(.hweights.default(N[1], L), length(N))
  }
}

wnorm.default <- wnorm.complex <- function(x, L = (N + 1) %/% 2, ...) {
  N <- length(x)

  # Compute weights
  w <- .hweights.default(x, L)

  # Compute wnorm
  sqrt(sum(w * abs(x)^2))
}

wnorm.nd.ssa <- wnorm.1d.ssa <- wnorm.toeplitz.ssa <- wnorm.cssa <- wnorm.mssa <- function(x, ...) {
  # Compute weights
  w <- .hweights(x)

  # Get series
  F <- unlist(.F(x))

  if (any(w == 0)) {
    # Omit uncovered elements
    F <- as.vector(F[w > 0])
    w <- as.vector(w[w > 0])
  }

  # Compute wnorm
  sqrt(sum(w * abs(F)^2))
}

frobenius.cor <- function(x, groups, ...) {
  # MB, add class check here too? We can just return identical matrix for non-ossa decompositions

  if (missing(groups))
    groups <- as.list(seq_len(nsigma(x)))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  idx <- unique(unlist(groups))
  sigma <- .sigma(x)[idx]
  U <- .U(x)[, idx, drop = FALSE]
  V <- calc.v(x, idx)

  # Compute frobenius covariation of elementary matrices
  cov <- crossprod(U) * crossprod(V) * tcrossprod(sigma)

  # Summing by groups (by columns and then by rows)
  cov <- sapply(groups, function(group) rowSums(cov[, match(group, idx), drop = FALSE]))
  cov <- t(cov)
  cov <- as.matrix(sapply(groups, function(group) rowSums(cov[, match(group, idx), drop = FALSE])))

  # Convert to correlations
  cor <- cov2cor(cov)

  # Fix possible numeric error
  cor[cor > 1] <- 1; cor[cor < -1] <- -1

  # Add class
  class(cor) <- "wcor.matrix"

  # Set names
  colnames(cor) <- rownames(cor) <- .group.names(groups)

  # Return
  cor
}

.is.frobenius.orthogonal <- function(x, groups, eps = sqrt(.Machine$double.eps), ...) {
  if (!inherits(x, "ossa"))
    return(TRUE)

  if (missing(groups))
    groups <- as.list(seq_len(nsigma(x)))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  fcor <- frobenius.cor(x, groups, ...)
  diag(fcor) <- 0

  max.fcor <- fcor[which.max(abs(fcor))]

  if (abs(max.fcor) < eps)
    TRUE
  else
    max.fcor
}

wcor.ossa <- function(x, groups, Fs, ..., cache = TRUE) {
  if (!missing(groups) && missing(Fs)) {
    isfcor <- .is.frobenius.orthogonal(x, groups, ...)
    if (!isTRUE(isfcor))
      warning(sprintf("Component matrices are not F-orthogonal (max F-cor is %s). W-cor matrix can be irrelevant",
                      format(isfcor, digits = 3)))
  }

  NextMethod()
}

#N = 399;
#a = 1.005;
#T = 200;
#F1 <- (1/a)^(1:N);
#F2 <- (a^(1:N))*cos(2*pi*(1:N)/T);
#.wcor(cbind(F1, F2));
