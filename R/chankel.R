#   R package for Singular Spectrum Analysis
#   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
#   Copyright (c) 2012 Alexander Shlemov <shlemovalex@gmail.com>
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

#   Routines for SSA with centering

chankel <- function(X, L, lcv, rcv) {
  N <- length(X);
  K <- N - L + 1;
  h <- outer(1:L, 1:K, function(x, y) X[x + y - 1]);

  h - outer(lcv, rep(1, K)) - outer(rep(1, L), rcv);
}

new.chmat <- function(F, L, lcv, rcv) {
  storage.mode(F) <- "double";
  storage.mode(L) <- "integer";
  storage.mode(lcv) <- "double";
  storage.mode(rcv) <- "double";
  .Call("initialize_chmat", F, L, lcv, rcv);
}

chcols <- function(ch) {
  .Call("chankel_cols", ch)
}

chrows <- function(ch) {
  .Call("chankel_rows", ch)
}

is.chmat <- function(ch) {
  .Call("is_chmat", ch)
}

chmatmul <- function(chmat, v, transposed = FALSE) {
  storage.mode(v) <- "double";
  storage.mode(transposed) <- "logical";
  .Call("chmatmul", chmat, v, transposed);
}

.hankelize.one.chankel <- function(U, V, ch) {
  storage.mode(U) <- storage.mode(V) <- "double";
  .Call("hankelize_one_fft_chankel", U, V, ch);
}

.hankelize.one.both.centering.ssa.nutrlan <-
.hankelize.one.both.centering.ssa.propack <- function(x, U, V) {
  ch <- .get(x, "chmat");
  .hankelize.one.chankel(U, V, ch);
}

ssums <- function(F, L){
  diff(c(0, cumsum(F)), lag = L);
}

decompose.both.centering.ssa <- function(x, ...) {
  N <- x$length; L <- x$window; K <- N - L + 1;

  F <- .get(x, "F");

  L.s <- min(L, K);
  dv <- c(1:(L.s-1), rep(L.s, N-2*L.s+2), (L.s-1):1);

  M <- sum(F * dv) / (L * K);
  ###########################
  m <- ssums(F - M, K) / K;

  m0 <- c(rep(0, K - 1), m, rep(0, K - 1));
  sm0 <- ssums(m0, K);
  F.t1 <- sm0 / dv + M;

  lam1 <- sqrt(sum(m^2) * K + M^2 * L * K);

  mrow <- m + M;
  u1 <- mrow / sqrt(sum(mrow^2));
  v1 <- rep(1/sqrt(K), K);

  mcol <- ssums(F - M, L) / L;

  m0 <- c(rep(0, L - 1), mcol, rep(0, L - 1));
  sm0 <- ssums(m0, L);

  F.t2 <- sm0 / dv;

  lam2 <- sqrt(sum(mcol^2) * L);
  v2 <- mcol / sqrt(sum(mcol^2));
  u2 <- rep(1/sqrt(L), L);

  .set(x, "F.trow", F.t1);
  .set(x, "F.tcol", F.t2);
  .set(x, "lam1", lam1);
  .set(x, "lam2", lam2);
  .set(x, "u1", u1);
  .set(x, "u2", u2);
  .set(x, "v1", v1);
  .set(x, "v2", v2);

  .decompose.centering(x, ...);
}

.decompose.centering.svd <- function(x,
                                 neig = min(L, K),
                                 ...,
                                 force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1;

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  F <- .get(x, "F");
  ch <- chankel(F, L, x$u1 * x$lam1 / sqrt(K), x$v2 * x$lam2 / sqrt(L));

  # Do decomposition
  S <- svd(ch, nu = neig, nv = neig);

  # Save results
  .set(x, "lambda", S$d);
  if (!is.null(S$u))
    .set(x, "U", S$u);
  if (!is.null(S$v))
    .set(x, "V", S$v);

  x;
}

.decompose.centering.both.centering.ssa.svd <- .decompose.centering.svd

.decompose.centering.both.centering.ssa.propack <- function(x,
                                         neig = min(50, L, K),
                                         ...,
                                         force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1;

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.")

  F <- .get(x, "F");
  ch <- new.chmat(F, L, x$u1 * x$lam1 / sqrt(K), x$v2 * x$lam2 / sqrt(L));

  S <- propack.svd(ch, neig = neig, ...);

  # Save results
  .set(x, "chmat", ch);
  .set(x, "lambda", S$d);
  if (!is.null(S$u))
    .set(x, "U", S$u);
  if (!is.null(S$v))
    .set(x, "V", S$v);

  x;
}

.decompose.centering.both.centering.ssa.nutrlan <- function(x,
    neig = min(50, L, K),
    ...) {
  N <- x$length; L <- x$window; K <- N - L + 1;

  ch <- .get(x, "chmat", allow.null = TRUE);
  if (is.null(ch)) {
    F <- .get(x, "F");
    ch <- new.chmat(F, L, x$u1 * x$lam1 / sqrt(K), x$v2 * x$lam2 / sqrt(L));
  }

  lambda <- .get(x, "lambda", allow.null = TRUE);
  U <- .get(x, "U", allow.null = TRUE);

  S <- trlan.svd(ch, neig = neig, ...,
      lambda = lambda, U = U);

  # Save results
  .set(x, "chmat", ch);
  .set(x, "lambda", S$d);
  if (!is.null(S$u))
    .set(x, "U", S$u);

  x;
}

.calc.v.centering.chankel <- function(x, idx) {
  lambda <- .get(x, "lambda")[idx];
  U <- .get(x, "U")[, idx, drop = FALSE];
  ch <- .get(x, "chmat");

  invisible(sapply(1:length(idx),
                   function(i) chmatmul(ch, U[, i], transposed = TRUE) / lambda[i]));
}

.calc.v.centering.svd <- function(x, idx, env) {
  # Check, if there is garbage-collected storage to hold some pre-calculated
  # stuff.
  if (identical(env, .GlobalEnv) ||
      !exists(".ssa.temporary.storage", envir = env, inherits = FALSE)) {
    F <- .get(x, "F");

    # Build chankel matrix.
    N <- x$length; L <- x$window; K <- N - L + 1;
    X <- chankel(F, L, x$u1 * x$lam1 / sqrt(K), x$v2 * x$lam2 / sqrt(L));

    # Save to later use, if possible.
    if (!identical(env, .GlobalEnv)) {
      assign(".ssa.temporary.storage", X, envir = env, inherits = FALSE);
    }
  } else {
    X <- get(".ssa.temporary.storage", envir = env, inherits = FALSE);
  }

  lambda <- .get(x, "lambda")[idx];
  U <- .get(x, "U")[, idx, drop = FALSE];

  invisible(sapply(1:length(idx),
                   function(i) crossprod(X, U[, i]) / lambda[i]));
}

calc.v.both.centering.ssa.svd <- function(x, idx, env = .GlobalEnv, ...) .calc.v.centering.svd(x, idx, env)
calc.v.both.centering.ssa.eigen <- function(x, idx, env = .GlobalEnv, ...) .calc.v.centering.chankel(x, idx)
calc.v.both.centering.ssa.propack <- function(x, idx, env = .GlobalEnv, ...) .calc.v.centering.chankel(x, idx)
calc.v.both.centering.ssa.nutrlan <- function(x, idx, env = .GlobalEnv, ...) .calc.v.centering.chankel(x, idx)
