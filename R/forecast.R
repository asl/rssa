#   R package for Singular Spectrum Analysis
#   Copyright (c) 2012 Alexander Shlemov <shlemovalex@gmail.com>
#   Copyright (c) 2012 Anton Korobeynikov <asl@math.spbu.ru>
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

check.for.groups <- function(use.group = TRUE) {
  # R magic to extract the call to parent function
  call <- as.list(match.call(definition = sys.function(-1),
                             call = sys.call(sys.parent())))
  if (use.group) {
    if (!is.null(call[["groups"]]))
      stop("`groups' found in the arguments, however function expect `group'")
  } else {
    if (!is.null(call[["group"]]))
      warning("`group' found in the arguments, however function expect `groups'")
  }
}

basis2lrr <- function(U, eps = sqrt(.Machine$double.eps)) {
  N <- nrow(U);
  lpf <- U %*% t(U[N, , drop = FALSE]);

  divider <- 1 - lpf[N]
  if (divider < eps)
    stop("Verticality coefficient equals to 1");

  lpf[-N] / divider
}

"lrr.1d-ssa" <- function(this, group, ...) {
  if (missing(group))
    group <- 1:min(nlambda(this), nu(this))

  check.for.groups(use.group = TRUE)

  # Determine the upper bound of desired eigentriples
  desired <- max(group)

  # Continue decomposition, if necessary
  if (desired > nu(this))
    decompose(this, ..., neig = desired)

  U <- .get(this, "U")[, group, drop = FALSE]

  res <- basis2lrr(U)
  class(res) <- "lrr"

  res
}

companion.matrix.lrr <- function(x) {
  n <- length(x)
  res <- matrix(0, n, n)
  res[, n] <- x
  res[seq(from = 2, by = n + 1, length.out = n - 1)] <- 1
  res
}

roots.lrr <- function(x, ..., method = c("companion", "polyroot")) {
  method <- match.arg(method)

  res <-
    if (identical(method, "polyroot")) {
      polyroot(c(-x, 1))
    } else {
      eigen(companion.matrix.lrr(x), only.values = TRUE)$values
    }

  res[order(abs(res), decreasing = TRUE)]
}

plot.lrr <- function(x, ..., raw = FALSE) {
  r <- roots(x)
  if (raw) {
    plot(r, ...)
  } else {
    xlim <- range(c(Re(r), +1, -1))
    ylim <- range(c(Im(r), +1, -1))

    plot(r, ...,
         xlim = xlim, ylim = ylim,
         main = "Roots of Linear Recurrence Formula",
         xlab = "Real Part",
         ylab = "Imaginary Part",
         asp = 1)
    symbols(0, 0, circles = 1, add = TRUE, inches = FALSE)
  }
}

apply.lrr <- function(F, lrr, len = 1, only.new = FALSE) {
  N <- length(F)
  r <- length(lrr)

  # Sanity check of inputs
  if (r > N)
    stop("Wrong length of LRR")

  # Run the actual LRR
  F <- c(F, rep(NA, len))
  for (i in 1:len)
    F[N+i] <- sum(F[(N+i-r) : (N+i-1)]*lrr)

  if (only.new) F[(N+1):(N+len)] else F
}

"rforecast.1d-ssa" <- function(this, groups, len = 1,
                               base = c("reconstructed", "original"),
                               only.new = TRUE,
                               ..., cache = TRUE) {
  base <- match.arg(base)
  if (missing(groups))
    groups <- as.list(1:min(nlambda(this), nu(this)))

  check.for.groups(use.group = FALSE)

  # Determine the upper bound of desired eigentriples
  desired <- max(unlist(groups));

  # Continue decomposition, if necessary
  if (desired > min(nlambda(this), nu(this)))
    decompose(this, ..., neig = desired);

  # Grab the reconstructed series if we're basing on them
  if (identical(base, "reconstructed"))
    r <- reconstruct(this, groups = groups, ..., cache = cache)

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]

    # Calculate the LRR corresponding to group
    lf <- lrr(this, group)

    # Calculate the forecasted values
    out[[i]] <- apply.lrr(if (identical(base, "reconstructed")) r[[i]] else .get(this, "F"),
                          lf, len, only.new = only.new)
    # FIXME: try to fixup the attributes
  }

  names(out) <- paste(sep = "", "F", 1:length(groups))

  # Forecasted series can be pretty huge...
  invisible(out)
}

"vforecast.1d-ssa" <- function(this, groups, len = 1,
                               only.new = TRUE,
                               ...) {
  L <- this$window
  K <- this$length - L + 1
  N <- K + L - 1 + len + L - 1
  N.res <- K + L - 1 + len

  if (missing(groups))
    groups <- as.list(1:min(nlambda(this), nu(this)))

  check.for.groups(use.group = FALSE)

  # Determine the upper bound of desired eigentriples
  desired <- max(unlist(groups))

  # Continue decomposition, if necessary
  if (desired > min(nlambda(this), nu(this)))
    decompose(this, ..., neig = desired)

  lambda <- .get(this, "lambda")
  U <- .get(this, "U")

  V <- if (nv(this) >= desired) .get(this, "V") else NULL

  # Make hankel matrix for fast hankelization (we use it for plan)
  h <- new.hmat(double(N), L)

  out <- list()
  for (i in seq_along(groups)) {
    group <- unique(groups[[i]])

    Uet <- U[, group, drop = FALSE]
    Vet <- if (is.null(V)) calc.v(this, idx = group) else V[, group, drop = FALSE]

    Z <- rbind(t(lambda[group] * t(Vet)), matrix(NA, len + L - 1, length(group)))

    U.head <- Uet[-L, , drop = FALSE]
    U.tail <- Uet[-1, , drop = FALSE]
    Pi <- Uet[L, ]
    tUhUt <- t(U.head) %*% U.tail
    P <- tUhUt + 1 / (1 - sum(Pi^2)) * Pi %*% (t(Pi) %*% tUhUt)

    for (j in (K + 1):(K + len + L - 1)) {
      Z[j, ] <- P %*% Z[j - 1, ]
    }

    res <- double(N)
    for (j in seq_along(group)) {
      res <- res + .hankelize.one.hankel(Uet[ , j], Z[ , j], h)
    }

    out[[i]] <- res[(if (only.new) (K+L):N.res else 1:N.res)]
    # FIXME: try to fixup the attributes
  }

  names(out) <- paste(sep = "", "F", 1:length(groups))

  # Forecasted series can be pretty huge...
  invisible(out)
}

"bforecast.1d-ssa" <- function(this, group,
                               len = 1, R = 100, level = 0.95,
                               type = c("recurrent", "vector"),
                               ...,
                               cache = TRUE) {
  type <- match.arg(type)
  check.for.groups(use.group = TRUE)
  dots <- list(...)

  # First, perform the reconstruction and calculate the residuals.
  r <- reconstruct(this, groups = list(group), ..., cache = cache)
  stopifnot(length(r) == 1)
  res <- residuals(r)

  forecast.fun <- if (identical(type, "recurrent")) rforecast else vforecast
  boot.forecast <- function(F, base) {
    s <- clone(base, copy.cache = FALSE, copy.storage = FALSE)
    .set(s, "F", F)
    .set(s, "Fattr", attributes(F))
    do.call(forecast.fun,
            c(list(s,
                   groups = list(group), len = len, only.new = TRUE),
              dots))[[1]]
  }

  # Do the actual bootstrap forecast
  bF <- matrix(nrow = len, ncol = R)
  bF[] <- replicate(R,
                    boot.forecast(r[[1]] + sample(res, replace = TRUE), this))

  # Finally, calculate the statistics of interest
  cf <- apply(bF, 1, quantile, probs = c((1-level) / 2, (1 + level) / 2))
  cbind(Value = rowMeans(bF), t(cf))
}

"sforecast.1d-ssa" <- function(this, group,
                               len = 1,
                               ...,
                               cache = TRUE) {
  check.for.groups(use.group = TRUE)

  # First, perform the reconstruction.
  r <- reconstruct(this, groups = list(group), ..., cache = cache)
  stopifnot(length(r) == 1)
  r <- r[[1]]

  # Calculate the LRR
  lf <- lrr(this, group)

  # Now calculate the recurrent forecast on sliding part of the series
  N <- this$length
  L <- this$window
  K <- N - L + 1
  K.l <- K - len

  if (N < L + len)
    stop("too large `len' for sliding forecast")

  res <- numeric(K.l)
  for (i in 1:K.l)
    res[i] <- apply.lrr(r[i:(i+L)], lf, len = len, only.new = TRUE)[len]

  res
}

"lrr.toeplitz-ssa" <- `lrr.1d-ssa`;
"vforecast.toeplitz-ssa" <- `vforecast.1d-ssa`;
"rforecast.toeplitz-ssa" <- `rforecast.1d-ssa`;
"bforecast.toeplitz-ssa" <- `bforecast.1d-ssa`;
"sforecast.toeplitz-ssa" <- `sforecast.1d-ssa`;

rforecast.ssa <- function(x, groups, len = 1,
                          base = c("reconstructed", "original"),
                          only.new = TRUE,
                          ..., cache = TRUE) {
  stop("generic recurrent forecast not implemented yet!")
}

bforecast.ssa <- function(x, group,
                          len = 1, R = 100, level = 0.95,
                          type = c("recurrent", "vector"),
                          ...,
                          cache = TRUE) {
  stop("generic bootstrapped forecast not implemented yet!")
}

lrr.ssa <- function(x, group) {
  stop("generic LRR calculation not implemented yet!")
}

vforecast.ssa <- function(x, groups, len = 1,
                          only.new = TRUE,
                          ...) {
  stop("generic vector forecast not implemented yet!")
}

sforecast.ssa <- function(x, group, len = 1,
                          ...,
                          cache = TRUE) {
  stop("generic sliding forecast not implemented yet!")
}


lrr <- function(this, ...)
  UseMethod("lrr")
roots <- function(x, ...)
  UseMethod("roots")
rforecast <- function(this, ...)
  UseMethod("rforecast")
vforecast <- function(this, ...)
  UseMethod("vforecast")
bforecast <- function(this, ...)
  UseMethod("bforecast")
sforecast <- function(this, ...)
  UseMethod("sforecast")
