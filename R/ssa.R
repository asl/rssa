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

# FIXME: add version, which accepts first column and last row
hankel <- function(X, L) {
  if (is.matrix(X) && nargs() == 1) {
     L <- dim(X)[1]; K <- dim(X)[2]; N <- K + L - 1;
     left  <- c(1:L, L*(2:K));
     right <- c(1+L*(0:(K-1)), ((K-1)*L+2):(K*L));
     v <- sapply(1:N, function(i) mean(X[seq(left[i], right[i], by = L-1)]));
     return (v);
  }

  # Coerce output to vector, if necessary
  if (!is.vector(X))
    X <- as.vector(X);
  N <- length(X);
  if (missing(L))
    L <- (N + 1) %/% 2;
  K <- N - L + 1;
  outer(1:L, 1:K, function(x,y) X[x+y-1]);
}

ssa.decompose <- function(x,
                          L = (N + 1) %/% 2,
                          method = c("hankel"),
                          centering = c("none", "row", "both")) {
  method <- match.arg(method);
  centering <- match.arg(centering);
  N <- length(x); K <- N - L + 1;
  
  if (identical(method, "hankel")) {   
    X <- hankel(x, L = L);

    # Perform extra centering, if needed
    rmeans <- cmeans <- cl <- cU <- cV <- NULL;
    .norm <- function(X) sqrt(sum(X^2));
    if (identical(centering, "row") || identical(centering, "both")) {
      rmeans <- rowMeans(X);
      X <- sweep(X, 1, rmeans);
    }
    if (identical(centering, "both")) {
      cmeans <- colMeans(X);
      X <- sweep(X, 2, cmeans);
    }

    # Compute eigentriples for first and second means
    if (!is.null(rmeans)) {
      mnorm <- .norm(rmeans);
      cl <- sqrt(mnorm * sqrt(K));
      cU <- rmeans / mnorm;
      cV <- rep(1, K) * mnorm / cl;
    }
    if (!is.null(cmeans)) {
      mnorm <- .norm(cmeans);
      cl <- c(cl, sqrt(mnorm * sqrt(L)));
      cU <- cbind(cU, rep(1, L) / sqrt(L));
      cV <- cbind(cV, sqrt(L) * cmeans / cl[2]);
    }

    # FIXME: Use special SVD for hankel matrixes
    S <- svd(X);
    names(S) <- c("lambda", "U", "V");
  } else {
    stop("Unknown method in SSA")
  }


  # Add centering eigentriples, if any
  if (!is.null(cl)) {
    S$lambda <- c(cl, S$lambda);
    S$U <- cbind(cU, S$U);
    S$V <- cbind(cV, S$V);
  }

  S$call   <- match.call();
  S$method <- method;
  S$centering <- centering;
  S$series <- deparse(substitute(x));
  S$window <- L;
  S$length <- length(x);
  class(S) <- "ssa";

  return (S);
}

ssa.reconstruct <- function(S, groups) {
  out <- list();

  if (missing(groups))
    groups <- as.list(1:length(S$lambda));

  for (i in seq_along(groups)) {
    group <- groups[[i]];
    X <- S$U[, group] %*% diag(S$lambda[group], nrow = length(group)) %*% t(S$V[,group])
    out[[i]] <- hankel(X);
  }
  names(out) <- paste("F", 1:length(groups), sep="");
  return (out);
}

#.F <- function(x) exp(-.01 * x)*cos(x/100);
#F <- .F(1:100);
#a = 1; b = 0; T = 10;
#.F <- function(x) a*x + b + 5*sin(2*pi*x/T);
#F <- .F(1:19);
#Rprof();
#SS <- ssa.decompose(F, centering = "both");
#L <- ssa.reconstruct(SS, list(1:2));
#Rprof(NULL);
#summaryRprof();
