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

wcor.default <- function(x, L = (N + 1) %/% 2, ...) {
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")

  N <- nrow(x);

  K <- N - L + 1;
  Ls <- min(L, K); Ks <- max(L, K);

  # Compute weights
  w <- c(1:(Ls-1), rep(Ls, Ks-Ls+1), seq(from = Ls-1, to = 1, by = -1));

  # Compute w-covariation
  cov <- crossprod(w * x, x)

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
  L <- x$window; N <- x$length;
  if (missing(groups))
    groups <- as.list(1:nlambda(x));

  # Compute reconstruction.
  F <- reconstruct(x, groups, ..., cache = cache);
  x <- matrix(unlist(F), nrow = N, ncol = length(groups));
  colnames(x) <- names(F);

  # Finally, compute w-correlations and return
  wcor.default(x, L = L)
}

wcor.ssa <- function(x, groups, ..., cache = TRUE)
  stop("Unsupported SVD method for SSA!");

wcor <- function(x, ...) {
  UseMethod("wcor");
}

clusterify.wcor.matrix <- function(x,
                                   nclust = N,
                                   ...,
                                   dist = function(X) (1 - X) / 2) {
  N <- nrow(x);
  h <- cutree(hclust(as.dist(dist(x)), ...), k = nclust);
  split(1:N, h);
}

wnorm.default <- function(x, L = (N + 1) %/% 2, ...) {
  # TODO Implement wnorm for 2dSSA
  N <- length(x)

  K <- N - L + 1
  Ls <- min(L, K); Ks <- max(L, K)

  # Compute weights
  w <- c(1:(Ls-1), rep(Ls, Ks-Ls+1), (Ls-1):1)

  # Compute wnorm
  sqrt(sum(w * x^2))
}

wnorm.1d.ssa <- wnorm.toeplitz.ssa <- function(x, ...)
  wnorm.default(x$F, x$window)

#N = 399;
#a = 1.005;
#T = 200;
#F1 <- (1/a)^(1:N);
#F2 <- (a^(1:N))*cos(2*pi*(1:N)/T);
#.wcor(cbind(F1, F2));
