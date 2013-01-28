# R package for Singular Spectrum Analysis
# Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
# Copyright (c) 2013 Alexander Shlemov <shlemovalex@gmail.com>
#
# This program is free software; you can redistribute it
# and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation;
# either version 2 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the
# Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
# MA 02139, USA.

# Routines for SSA with centering

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

  .Call("initialize_hmat", F, L, lcv, rcv);
}

.hankelize.one.both.centering.ssa.eigen <-
.hankelize.one.both.centering.ssa.nutrlan <-
.hankelize.one.both.centering.ssa.propack <- function(x, U, V) {
  ch <- .get(x, "chmat");
  .hankelize.one.hankel(U, V, ch);
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

decompose.row.centering.ssa <- function(x, ...) {
  ssums <- function(F, L){
    diff(c(0, cumsum(F)), lag = L);
  }

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

  F.t2 <- rep(0, N);

  lam2 <- 0;
  v2 <- rep(0, K);
  u2 <- rep(0, L);

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

.decompose.centering <- function(x, ...)
  UseMethod(".decompose.centering")

.decompose.centering.both.centering.ssa.svd <-
.decompose.centering.row.centering.ssa.svd <- function(x,
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

.decompose.centering.both.centering.ssa.eigen <-
.decompose.centering.row.centering.ssa.eigen <- function(x, ...,
        force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1;

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decompostion is not supported for this method.")

  # Build hankel matrix (this can be done more efficiently!)
  F <- .get(x, "F");
  ch <- chankel(F, L, x$u1 * x$lam1 / sqrt(K), x$v2 * x$lam2 / sqrt(L));

  # Do decomposition
  S <- eigen(tcrossprod(ch));

  # Fix small negative values
  S$values[S$values < 0] <- 0;

  # Save results
  .set(x, "lambda", sqrt(S$values));
  .set(x, "U", S$vectors);

  # Save chmat object for efficient hankelization
  .set(x, "chmat",  new.chmat(F, L, x$u1 * x$lam1 / sqrt(K), x$v2 * x$lam2 / sqrt(L)));

  x;
}

.decompose.centering.both.centering.ssa.propack <-
.decompose.centering.row.centering.ssa.propack <- function(x,
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

.decompose.centering.both.centering.ssa.nutrlan <-
.decompose.centering.row.centering.ssa.nutrlan <- function(x,
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

  invisible(sapply(seq_along(idx),
          function(i) hmatmul(ch, U[, i], transposed = TRUE) / lambda[i]));
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

  invisible(sapply(seq_along(idx),
          function(i) crossprod(X, U[, i]) / lambda[i]));
}

calc.v.both.centering.ssa.svd <-
calc.v.row.centering.ssa.svd <- function(x, idx, env = .GlobalEnv, ...) .calc.v.centering.svd(x, idx, env)

calc.v.both.centering.ssa.eigen <-
calc.v.row.centering.ssa.eigen <- function(x, idx, env = .GlobalEnv, ...) .calc.v.centering.chankel(x, idx)

calc.v.both.centering.ssa.propack <-
calc.v.row.centering.ssa.propack <- function(x, idx, env = .GlobalEnv, ...) .calc.v.centering.chankel(x, idx)

calc.v.both.centering.ssa.nutrlan <-
calc.v.row.centering.ssa.nutrlan <- function(x, idx, env = .GlobalEnv, ...) .calc.v.centering.chankel(x, idx)


divide.groups <- function(x, groups) {
  all.special.components.idx <- if ("row.centering.ssa" %in% class(x)) {
        -1;
      } else {
        c(-1, -2);
      }

  special.groups <- list();
  SVD.groups <- list();

  for (i in seq_along(groups)) {
    group <- groups[[i]];

    # Expand zeroes
    if (0 %in% group) {
      group <- group[group != 0];
      group <- c(group, all.special.components.idx);
    }
    group <- unique(group);

    # Extract special componemts
    special.groups[[i]] <- group[group < 0];

    # Exclude special components
    SVD.groups[[i]] <- group[group > 0];
  }

  list(special.groups = special.groups, SVD.groups = SVD.groups);
}

reconstruct.both.centering.ssa <-
reconstruct.row.centering.ssa <- function(x, groups, ..., drop = FALSE, cache = TRUE) {
  if (missing(groups))
    groups <- as.list(0:min(nlambda(x), nu(x)));

  divided.groups <- divide.groups(x, groups);
  SVD.groups <- divided.groups$SVD.groups;
  special.groups <- divided.groups$special.groups;

  # Call trival reconstruct
  out <- reconstruct.ssa(x, SVD.groups, ..., drop = drop, cache = cache);

  # Add centering components
  F.trow <- .get(x, "F.trow", allow.null = TRUE);
  F.tcol <- .get(x, "F.tcol", allow.null = TRUE);
  for (i in seq_along(groups)) {
    special.group <- special.groups[[i]];

    if (-1 %in% special.group) {
      out[[i]] <- out[[i]] + F.trow;
    }

    if (-2 %in% special.group) {
      out[[i]] <- out[[i]] + F.tcol;
    }
  }

  # Fix residuals
  rspecial.groups <- unique(unlist(special.groups));
  residuals <- attr(out, "residuals");
  if (-1 %in% rspecial.groups) {
    residuals <- residuals - F.trow;
  }
  if (-2 %in% rspecial.groups) {
    residuals <- residuals - F.tcol;
  }

  attr(out, "residuals") <- residuals;

  invisible(out);
}
