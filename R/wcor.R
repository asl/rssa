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

wcor.default <- function(X, L = (N + 1) %/% 2, ...) {
  if (is.data.frame(X))
    X <- as.matrix(X)
  else if (!is.matrix(X))
    stop("'X' must be a matrix or a data frame")
  if (!all(is.finite(X)))
    stop("'X' must contain finite values only")

  N <- nrow(X);

  K <- N - L + 1;
  Ls <- min(L, K); Ks <- max(L, K);

  # Compute weights
  w <- c(1:(Ls-1), rep(Ls, Ks-Ls+1), seq(from = Ls-1, to = 1, by = -1));

  # Compute w-covariation
  cov <- crossprod(sqrt(w) * X);

  # Convert to correlations
  Is <- 1/sqrt(diag(cov));
  R <- Is * cov * rep(Is, each = nrow(cov));
  class(R) <- "wcor.matrix";

  return (R);
}

"wcor.toeplitz-ssa" <- "wcor.1d-ssa" <- function(X, groups, ..., cache = TRUE) {
  L <- X$window; N <- X$length;
  if (missing(groups))
    groups <- as.list(1:nlambda(X));

  # Compute reconstruction.
  F <- reconstruct(X, groups, ..., cache = cache);
  X <- matrix(unlist(F), nrow = N, ncol = length(groups));
  colnames(X) <- names(F);

  # Finally, compute w-correlations and return
  wcor.default(X, L = L)
}

wcor.ssa <- function(X, groups, ..., cache = TRUE)
  stop("Unsupported SVD method for SSA!");

wcor <- function(X, ...) {
  UseMethod("wcor");
}

# FIXME: Add legend
plot.wcor.matrix <- function(x, col = rev(gray(seq(0, 1, len = 20))),
                             xlab = "", ylab = "",
                             main = "W-correlation Matrix",
                             ...) {
  image(1:ncol(x), 1:nrow(x), abs(x), col = col,
        axes = FALSE, xlab = xlab, ylab = ylab, main = main,
        ...);
  axis(1, at = 1:ncol(x), labels = colnames(x));
  axis(2, at = 1:nrow(x), labels = rownames(x), las = 2);
  box();
}

clusterify.wcor.matrix <- function(this,
                                   nclust = N,
                                   ...,
                                   dist = function(X) (1 - X) / 2) {
  N <- nrow(this);
  h <- cutree(hclust(as.dist(dist(this)), ...), k = nclust);
  split(1:N, h);
}

#N = 399;
#a = 1.005;
#T = 200;
#F1 <- (1/a)^(1:N);
#F2 <- (a^(1:N))*cos(2*pi*(1:N)/T);
#.wcor(cbind(F1, F2));
