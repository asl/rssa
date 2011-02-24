#   R package for Singular Spectrum Analysis
#   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
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

#   Routines for toeplitz SSA

tcirc.old <- function(F, L = (N - 1) %/% 2) {
  N <- length(F);

  R <- as.vector(acf(F, lag.max = L - 1, type = "covariance", plot = FALSE, demean = FALSE)$acf);

  .res <- list();
  .res$C <- as.vector(fft(c(R, rev(R[-1]))));
  .res$L <- L;
  return (.res);
}

tmatmul.old <- function(C, v) {
  v <- as.vector(fft(C$C * fft(c(v, rep(0, C$L-1))), inverse = TRUE));
  Re((v/length(C$C))[1:C$L]);
}

new.tmat <- function(F,
                     L = (N - 1) %/% 2) {
  N <- length(F);

  # FIXME: Perform estimation via FFT (to drop complexity from O(N^2) to O(N log N)
  R <- as.vector(acf(F, lag.max = L - 1, type = "covariance", plot = FALSE, demean = FALSE)$acf);

  storage.mode(R) <- "double";
  t <- .Call("initialize_tmat", R);
}

tcols <- function(t) {
  .Call("toeplitz_cols", t)
}

trows <- function(t) {
  .Call("toeplitz_rows", t)
}

is.tmat <- function(t) {
  .Call("is_tmat", t)
}

tmatmul <- function(tmat, v, transposed = FALSE) {
  storage.mode(v) <- "double";
  storage.mode(transposed) <- "logical";
  .Call("tmatmul", tmat, v, transposed);
}

"decompose.toeplitz-ssa.nutrlan" <- function(x,
                                             neig = min(50, L, K),
                                             ...) {
  N <- x$length; L <- x$window; K <- N - L + 1;

  h <- .get(x, "hmat", allow.null = TRUE);
  if (is.null(h)) {
    F <- .get(x, "F");
    h <- new.hmat(F, L = L);
  }

  olambda <- .get(x, "olambda", allow.null = TRUE);
  U <- .get(x, "U", allow.null = TRUE);

  T <- new.tmat(F, L = L);

  S <- trlan.eigen(T, neig = neig, ...,
                   lambda = olambda, U = U);

  # Fix small negative values
  S$values[S$values < 0] <- 0;

  # Save results
  .set(x, "hmat", h);
  .set(x, "olambda", S$d);
  if (!is.null(S$u))
    .set(x, "U", S$u);

  num <- length(S$d);
  lambda <- numeric(num);
  V <- matrix(nrow = K, ncol = num);
  for (i in 1:num) {
    Z <- hmatmul(h, S$u[, i], transposed = TRUE);
    lambda[i] <- sum(Z^2);
    V[, i] <- Z / lambda[i];
  }

  # Save results
  .set(x, "lambda", lambda);
  .set(x, "V", V);

  x;
}

"decompose.toeplitz-ssa.eigen" <- function(x, ...,
                                           force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1;

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decompostion is not supported for this method.")

  # Build hankel matrix
  F <- .get(x, "F");
  h <- new.hmat(F, L = L);

  # Do decomposition
  if ("neig" %in% names(list(...)))
    warning("'neig' option ignored for SSA method 'eigen', computing EVERYTHING",
            immediate. = TRUE)

  R <- acf(F, lag.max = L - 1, type = "covariance", plot = FALSE, demean = FALSE);
  # FIXME: find a better way to construct toeplitz matrix
  C <- toeplitz(as.vector(R$acf));
  S <- eigen(C, symmetric = TRUE);

  # Fix small negative values
  S$values[S$values < 0] <- 0;
  
  .set(x, "U", S$vectors);

  lambda <- numeric(L);
  V <- matrix(nrow = K, ncol = L);
  for (i in 1:L) {
    Z <- hmatmul(h, S$vectors[,i], transposed = TRUE);
    lambda[i] <- sum(Z^2);
    V[, i] <- Z / lambda[i]; 
  }
  
  # Save results
  .set(x, "lambda", lambda);
  .set(x, "V", V);
  
  x;
}

"decompose.toeplitz-ssa.svd" <- function(x, ...) {
  stop("'SVD' method is not applicable to toeplitz SSA");
}

"decompose.toeplitz-ssa.propack" <- function(x, ...) {
  stop("'PROPACK' method is not applicable to toeplitz SSA");
}