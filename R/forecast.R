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

lrr.default <- function(U, eps = sqrt(.Machine$double.eps), ...) {
  N <- nrow(U);
  lpf <- U %*% t(U[N, , drop = FALSE]);

  divider <- 1 - lpf[N]
  if (divider < eps)
    stop("Verticality coefficient equals to 1");

  lpf[-N] / divider
}

lrr.1d.ssa <- function(x, groups, ..., drop = TRUE) {
  if (missing(groups))
    groups <- 1:min(nlambda(x), nu(x))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    U <- .get(x, "U")[, group, drop = FALSE]

    res <- lrr.default(U, ...)
    class(res) <- "lrr"

    out[[i]] <- res
  }

  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
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


maybe.fixup.attributes <- function(x, v,
                                   only.new = TRUE, drop = FALSE) {
  # Grab the initial set of attributes
  res <- attributes(v)
  # Reconstruct the original series
  F <- .get(x, "F")
  if (!drop)
    attributes(F) <- .get(x, "Fattr")

  # Try to guess the indices of known time series classes
  if (is.ts(F)) {
    return (ts(v,
               start = if (only.new) tsp(F)[2] + 1/frequency(F) else start(F),
               frequency = frequency(F)))
  }

  return(v)
}

rforecast.1d.ssa <- function(x, groups, len = 1,
                             base = c("reconstructed", "original"),
                             only.new = TRUE,
                             ...,
                             drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  L <- x$window
  K <- x$length - L + 1

  base <- match.arg(base)
  if (missing(groups))
    groups <- as.list(1:min(nlambda(x), nu(x)))

  # Grab the reconstructed series if we're basing on them
  if (identical(base, "reconstructed"))
    r <- reconstruct(x, groups = groups, ..., cache = cache)

  # Calculate the LRR corresponding to groups
  lf <- lrr(x, groups = groups, drop = FALSE)
  stopifnot(length(lf) == length(groups))

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]

    # Calculate the forecasted values
    out[[i]] <- apply.lrr(if (identical(base, "reconstructed")) r[[i]] else .get(x, "F"),
                          lf[[i]], len, only.new = only.new)
    out[[i]] <- maybe.fixup.attributes(x, out[[i]], only.new = only.new, drop = drop.attributes)
  }

  names(out) <- paste(sep = "", "F", 1:length(groups))
  if (length(out) == 1 && drop)
    out <- out[[1]]

  # Forecasted series can be pretty huge...
  invisible(out)
}

vforecast.1d.ssa <- function(x, groups, len = 1,
                             only.new = TRUE,
                             ...,
                             drop = TRUE, drop.attributes = FALSE) {
  L <- x$window
  K <- x$length - L + 1
  N <- K + L - 1 + len + L - 1
  N.res <- K + L - 1 + len

  if (missing(groups))
    groups <- as.list(1:min(nlambda(x), nu(x)))

  # Continue decomposition, if necessary
  desired <- .maybe.continue(x, groups = groups, ...)

  lambda <- .get(x, "lambda")
  U <- .get(x, "U")

  V <- if (nv(x) >= desired) .get(x, "V") else NULL

  # Grab the FFT plan
  fft.plan <- fft.plan.1d(N)

  out <- list()
  for (i in seq_along(groups)) {
    group <- unique(groups[[i]])

    Uet <- U[, group, drop = FALSE]
    Vet <- if (is.null(V)) calc.v(x, idx = group) else V[, group, drop = FALSE]

    Z <- rbind(t(lambda[group] * t(Vet)), matrix(NA, len + L - 1, length(group)))

    U.head <- Uet[-L, , drop = FALSE]
    U.tail <- Uet[-1, , drop = FALSE]
    Pi <- Uet[L, ]
    tUhUt <- t(U.head) %*% U.tail
    P <- tUhUt + 1 / (1 - sum(Pi^2)) * Pi %*% (t(Pi) %*% tUhUt)

    for (j in (K + 1):(K + len + L - 1)) {
      Z[j, ] <- P %*% Z[j - 1, ]
    }

    res <- rowSums(.hankelize.multi.hankel(Uet, Z, fft.plan))

    out[[i]] <- res[(if (only.new) (K+L):N.res else 1:N.res)]
    out[[i]] <- maybe.fixup.attributes(x, out[[i]], only.new = only.new, drop = drop.attributes)
  }

  names(out) <- paste(sep = "", "F", 1:length(groups))

  if (length(out) == 1 && drop)
    out <- out[[1]]

  # Forecasted series can be pretty huge...
  invisible(out)
}

bforecast.1d.ssa <- function(x, groups,
                             len = 1, R = 100, level = 0.95,
                             type = c("recurrent", "vector"),
                             ...,
                             drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  type <- match.arg(type)
  dots <- list(...)
  if (missing(groups))
    groups <- list(1:min(nlambda(x), nu(x)))

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    # First, perform the reconstruction and calculate the residuals.
    r <- reconstruct(x, groups = list(group), ..., cache = cache)
    stopifnot(length(r) == 1)
    res <- residuals(r)

    forecast.fun <- if (identical(type, "recurrent")) rforecast else vforecast
    boot.forecast <- function(F, base) {
      s <- clone(base, copy.cache = FALSE, copy.storage = FALSE)
      .set(s, "F", F)
      .set(s, "Fattr", attributes(F))
      do.call(forecast.fun,
              c(list(s,
                     groups = list(group), len = len, drop = TRUE, only.new = TRUE),
                dots))
    }

    # Do the actual bootstrap forecast
    bF <- matrix(nrow = len, ncol = R)
    bF[] <- replicate(R,
                      boot.forecast(r[[1]] + sample(res, replace = TRUE), x))

    # Finally, calculate the statistics of interest
    cf <- apply(bF, 1, quantile, probs = c((1-level) / 2, (1 + level) / 2))
    out[[i]] <- maybe.fixup.attributes(x, cbind(Value = rowMeans(bF), t(cf)), drop = drop.attributes)
  }

  names(out) <- paste(sep = "", "F", 1:length(groups))
  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

predict.1d.ssa <- function(object,
                           groups, len = 1,
                           method = c("recurrent", "vector", "bootstrap-recurrent", "bootstrap-vector"),
                           ...,
                           drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  method <- match.arg(method)
  dots <- list(...)

  # Calculate a forecast
  switch (method,
          recurrent = do.call(rforecast, c(list(object, groups = groups, len = len, drop = drop, drop.attributes = drop.attributes, cache = cache), dots)),
          vector = do.call(vforecast, c(list(object, groups = groups, len = len, cache = cache), dots)),
          'bootstrap-recurrent' =  do.call(bforecast, c(list(object, type = "recurrent", groups = groups, len = len, drop = drop, drop.attributes = drop.attributes, cache = cache), dots)),
          'bootstrap-vector' =  do.call(bforecast, c(list(object, type = "vector", groups = groups, len = len, drop = drop, drop.attributes = drop.attributes, cache = cache), dots)))
}

forecast.1d.ssa <- function(object,
                            groups, len = 1,
                            method = c("recurrent", "vector", "bootstrap-recurrent", "bootstrap-vector"),
                            ...,
                            drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  method <- match.arg(method)
  dots <- list(...)

  # Perform the forecast
  f <- do.call(predict, c(list(object, groups = groups, len = len, method = method, drop = drop, drop.attributes = drop.attributes, cache = cache), dots))

  # Now perform a "cast" to forecast object
  require(forecast)
  F <- .get(object, "F")
  if (!drop.attributes)
    attributes(F) <- .get(object, "Fattr")
  out <- list()
  for (i in seq_along(groups)) {
    # Perform the reconstruction. We cannot do all-at-once, because we need proper residuals as well.
    r <- reconstruct(object, groups = groups[i], ..., drop.attributes = drop.attributes, cache = cache)
    stopifnot(length(r) == 1)

    res <- list(model = object,
                method = switch(method,
                                recurrent = "SSA (recurrent)",
                                recurrent = "SSA (vector)",
                                `bootstrap-recurrent` = "SSA (bootstrap recurrent)",
                                `bootstrap-vector` = "SSA (bootstrap vector)"),
                fitted = r[[1]],
                residuals = residuals(r),
                x = F)
    # Handle bootstrap forecast separately
    if (method %in% c( "bootstrap-recurrent", "bootstrap-vector")) {
      nbnd <- (ncol(f) - 1) / 2
      res$mean  <- f[, "Value"]
      res$lower <- f[, seq(from = 2, by = 1, length.out = nbnd)]
      res$upper <- f[, seq(from = 2 + nbnd, by = 1, length.out = nbnd)]
      # HACK! Need to change if bforecast defaults will be changed!
      if (is.element("level", names(dots))) res$level <- 100*dots$level else res$level <- 100*0.95
    } else {
      res$mean <- f
    }

    class(res) <- "forecast"
    out[[i]] <- res
  }

  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

"lrr.toeplitz.ssa" <- `lrr.1d.ssa`;
"vforecast.toeplitz.ssa" <- `vforecast.1d.ssa`;
"rforecast.toeplitz.ssa" <- `rforecast.1d.ssa`;
"bforecast.toeplitz.ssa" <- `bforecast.1d.ssa`;
"forecast.toeplitz.ssa" <- `forecast.1d.ssa`;
"predict.toeplitz.ssa" <- `predict.1d.ssa`;

lrr <- function(x, ...)
  UseMethod("lrr")
roots <- function(x, ...)
  UseMethod("roots")
rforecast <- function(x, ...)
  UseMethod("rforecast")
vforecast <- function(x, ...)
  UseMethod("vforecast")
bforecast <- function(x, ...)
  UseMethod("bforecast")
