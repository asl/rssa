#   R package for Static Singular Analysis
#   Copyright (c) 2008 Anton Korobeynikov <asl@math.spbu.ru>
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

# FIXME: more checks
wcor.default <- function(X, L) {
  if (is.data.frame(X))
    X <- as.matrix(X)
  else if (!is.matrix(X))
    stop("'X' must be a matrix or a data frame")
  if (!all(is.finite(X)))
    stop("'X' must contain finite values only")

  if (missing(L))
    L <- ncol(X);
  N <- nrow(X);

  K <- N - L + 1;
  Ls <- min(L, K); Ks <- max(L, K);

  # Compute weights
  w <- c(1:Ls, rep(Ls, Ks-Ls), seq(from = N-Ks, to = 1, by = -1));

  # Compute w-covariation
  cov <- crossprod(sqrt(w) * X);

  # Convert to correlations
  Is <- 1/sqrt(diag(cov));
  R <- Is * cov * rep(Is, each = nrow(cov));
  class(R) <- "wcor.matrix";

  return (R);
}

wcor.ssa <- function(S, groups) {
  L <- S$window;
  if (missing(groups))
    groups <- as.list(1:L);

  # Compute reconstruction.
  # FIXME: Modify ssa.reconstruct to return matrix
  F <- ssa.reconstruct(S, groups);
  X <- matrix(nrow = S$length, ncol = 0);
  for (i in seq_along(F)) {
    X <- cbind(X, F[[i]]);
  }
  colnames(X) <- names(F);

  # Finally, compute w-correlations and return
  wcor.default(X);
}

wcor <- function(X, ...) {
  UseMethod("wcor");
}

#N = 399;
#a = 1.005;
#T = 200;
#F1 <- (1/a)^(1:N);
#F2 <- (a^(1:N))*cos(2*pi*(1:N)/T);
#.wcor(cbind(F1, F2));
