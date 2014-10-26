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

lrr.default <- function(x, eps = sqrt(.Machine$double.eps),
                        reverse = FALSE,
                        ...,
                        orthonormalize = TRUE) {
  if (orthonormalize) {
    U <- qr.Q(qr(x))
  } else {
    U <- x
  }

  N <- nrow(U)

  # Return zero LRR coefficients for zero subspace
  if (ncol(U) == 0) return(rep(0, N - 1))

  idx <- if (!reverse) N else 1
  lpf <- Conj(U) %*% t(U[idx, , drop = FALSE])

  divider <- 1 - lpf[idx]
  if (Mod(divider) < eps)
    stop("Verticality coefficient equals to 1")

  lpf[-idx] / divider
}

lrr.1d.ssa <- function(x, groups,
                       reverse = FALSE,
                       ..., drop = TRUE) {
  if (is.shaped(x))
    stop("`LRR is not implemented for shaped SSA case yet")

  if (missing(groups))
    groups <- 1:min(nsigma(x), nu(x))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  out <- list()
  for (i in seq_along(groups)) {
      res <- lrr.default(.colspan(x, groups[[i]]), reverse = reverse,
                         ..., orthonormalize = FALSE)
    class(res) <- "lrr"

    out[[i]] <- res
  }

  names(out) <- .group.names(groups)
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

apply.lrr <- function(F, lrr, len = 1, only.new = FALSE,
                      drift = 0, reverse = FALSE) {
  # Recycle drifts if needed
  if (length(drift) != len) {
    drift <- rep(drift, len)[seq_len(len)]
  }

  N <- length(F)
  r <- length(lrr)

  # Sanity check of inputs
  if (r > N)
    stop("Wrong length of LRR")

  # Run the actual LRR
  if (!reverse) {
    F <- c(F, rep(NA, len))
    for (i in 1:len)
      F[N+i] <- sum(F[(N+i-r) : (N+i-1)]*lrr) + drift[i]

    if (only.new) F[(N+1):(N+len)] else F
  } else {
    F <- c(rep(NA, len), F)

    for (i in 1:len)
      F[len-i+1] <- sum(F[(len-i+1 + 1) : (len-i+1 + r)]*lrr) + drift[len-i+1]

    if (only.new) F[1:len] else F
  }
}

rforecast.1d.ssa <- function(x, groups, len = 1,
                             base = c("reconstructed", "original"),
                             only.new = TRUE,
                             reverse = FALSE,
                             ...,
                             drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  if (is.shaped(x))
    stop("`forecasting is not implemented for shaped SSA case yet")

  if (x$circular)
    stop("forecasting is not properly defined for circular SSA")

  L <- x$window

  base <- match.arg(base)
  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))

  # Grab the reconstructed series if we're basing on them
  if (identical(base, "reconstructed"))
    r <- reconstruct(x, groups = groups, ..., cache = cache)

  # Calculate the LRR corresponding to groups
  lf <- lrr(x, groups = groups, reverse = reverse, drop = FALSE)
  stopifnot(length(lf) == length(groups))

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]

    F <- if (identical(base, "reconstructed")) as.vector(r[[i]]) else .F(x)

    # Calculate the forecasted values
    out[[i]] <- apply.lrr(F, lf[[i]], len, only.new = only.new, reverse = reverse)
    out[[i]] <- .apply.attributes(x, out[[i]],
                                  fixup = TRUE, reverse = reverse,
                                  only.new = only.new, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  # Forecasted series can be pretty huge...
  invisible(out)
}

.na.bind <- function(original, new,
                     update.method = c("append", "replace")) {
  update.method <- match.arg(update.method)

  removed <- attr(original, "na.action")

  res <- switch(update.method,
                append = c(original, new),
                replace = new)

  if (!is.null(removed)) {
    all.old <- seq_len(length(removed) + length(original))
    full.old <- setdiff(all.old, removed)
    full.new <- c(full.old, max(full.old) + seq_len(length(res) - length(original)))
    res.na <- setdiff(all.old, full.new)
    attr(res, "na.action") <- if (length(res.na) > 0) res.na else NULL
  }

  res
}

rforecast.mssa <- function(x, groups, len = 1,
                           base = c("reconstructed", "original"),
                           direction = c("row", "column"),
                           only.new = TRUE,
                           ...,
                           drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  L <- x$window; N <- x$length; K <- N - L + 1

  cK <- cumsum(K)
  cKr <- cumsum(K - 1)
  cKl <- cKr - (K - 1) + 1

  base <- match.arg(base)
  direction <- match.arg(direction)
  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))

  # Grab the reconstructed series if we're basing on them
  if (identical(base, "reconstructed"))
    r <- reconstruct(x, groups = groups, ..., cache = cache)

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]

    F <- if (identical(base, "reconstructed")) .to.series.list(r[[i]]) else .F(x)

    # Calculate the forecasted values
    if (identical(direction, "column")) {
      # Calculate the LRR corresponding to group
      lf <- lrr(x, groups = list(group), drop = FALSE)
      stopifnot(length(lf) == 1)
      R <- matrix(NA, nrow = len, ncol = length(N))
      R[] <- sapply(F,
                    apply.lrr,
                    lrr = lf[[1]], len = len, only.new = TRUE)
    } else {
      V <- calc.v(x, idx = group)

      # Build W
      W <- V[cK,, drop = FALSE]
      # Build Q
      Q <- V[-cK,, drop = FALSE]

      # Calculate the forecasted values
      qIWWt <- qr(diag(length(N)) - tcrossprod(W))
      WtQ <- W %*% t(Q)

      R <- matrix(NA, nrow = len, ncol = length(N))
      # Build initial Z
      Z <- unlist(lapply(seq_along(F), function(idx) F[[idx]][seq(to = N[[idx]], length.out = K[[idx]] -1)]))
      for (idx in seq_len(len)) {
        # Calculate the projection
        cR <- t(qr.coef(qIWWt, WtQ %*% Z))
        R[idx, ] <- cR

        # Shift the Z vector
        Z[-cKr] <- Z[-cKl]
        Z[cKr] <- cR
      }
    }

    out[[i]] <- if (only.new) .to.series.list(R) else {
      res <- lapply(seq_along(F), function(idx) .na.bind(F[[idx]], R[, idx], update.method = "append"))
      class(res) <- "series.list"
      res
    }
    out[[i]] <- .apply.attributes(x, out[[i]],
                                  fixup = TRUE,
                                  only.new = only.new, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  # Forecasted series can be pretty huge...
  invisible(out)
}

vforecast.1d.ssa <- function(x, groups, len = 1,
                             only.new = TRUE,
                             ...,
                             drop = TRUE, drop.attributes = FALSE) {
  if (is.shaped(x))
    stop("`forecasting is not implemented for shaped SSA case yet")

  if (x$circular)
    stop("forecasting is not properly defined for circular SSA")

  L <- x$window
  K <- x$length - L + 1
  N <- K + L - 1 + len + L - 1
  N.res <- K + L - 1 + len

  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))

  # Continue decomposition, if necessary
  desired <- .maybe.continue(x, groups = groups, ...)

  sigma <- .sigma(x)
  U <- .U(x)
  V <- if (nv(x) >= desired) .V(x) else NULL

  # Grab the FFT plan
  fft.plan <- fft.plan.1d(N, L = L)

  out <- list()
  for (i in seq_along(groups)) {
    group <- unique(groups[[i]])

    Uet <- U[, group, drop = FALSE]
    Vet <- if (is.null(V)) calc.v(x, idx = group) else V[, group, drop = FALSE]
    Z <- rbind(t(sigma[group] * t(Vet)), matrix(NA, len + L - 1, length(group)))

    P <- shift.matrix(Uet)

    for (j in (K + 1):(K + len + L - 1)) {
      Z[j, ] <- P %*% Z[j - 1, ]
    }

    res <- rowSums(.hankelize.multi(Uet, Z, fft.plan))

    out[[i]] <- res[(if (only.new) (K+L):N.res else 1:N.res)]
    out[[i]] <- .apply.attributes(x, out[[i]],
                                  fixup = TRUE,
                                  only.new = only.new, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  # Forecasted series can be pretty huge...
  invisible(out)
}

vforecast.mssa <- function(x, groups, len = 1,
                           direction = c("row", "column"),
                           only.new = TRUE,
                           ...,
                           drop = TRUE, drop.attributes = FALSE) {
  direction <- match.arg(direction)
  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))

  # Continue decomposition, if necessary
  desired <- .maybe.continue(x, groups = groups, ...)

  F <- .F(x)

  dec <- .decomposition(x)
  sigma <- .sigma(dec)
  V <- if (nv(x) >= desired) .rowspan(dec) else NULL

  L <- x$window
  K <- x$length - L + 1
  N.res <- K + L - 1 + len
  N <- N.res + switch(direction, column = L, row = K) - 1

  # Grab the FFT plan
  fft.plan <- switch(direction,
                     column = lapply(N, fft.plan.1d, L = L),
                     row = mapply(fft.plan.1d, N = N, L = L + len + K - 1))

  cK <- cumsum(K)
  cKs <- cK - K + 1

  out <- list()

  for (i in seq_along(groups)) {
    group <- unique(groups[[i]])

    Uet <- .colspan(dec, group)
    Vet <- if (is.null(V)) calc.v(x, idx = group) else V[, group, drop = FALSE]

    if (identical(direction, "column")) {
      U.head <- Uet[-L, , drop = FALSE]
      U.tail <- Uet[-1, , drop = FALSE]
      Pi <- Uet[L, ]
      tUhUt <- crossprod(U.head, U.tail)
      P <- tUhUt + 1 / (1 - sum(Pi^2)) * Pi %*% (t(Pi) %*% tUhUt)

      R <- lapply(seq_along(N), function(idx) {
          Z <- rbind(t(sigma[group] * t(Vet[cKs[idx] : cK[idx], , drop = FALSE])), matrix(NA, len + L - 1, length(group)))

          for (j in (K[idx] + 1) : (K[idx] + len + L - 1)) {
            Z[j, ] <- P %*% Z[j - 1, ]
          }

          rowSums(.hankelize.multi(Uet,
                                   Z,
                                   fft.plan[[idx]]))
        })
    } else if (identical(direction, "row")) {
      V.head <- Vet[-cK, , drop = FALSE]
      V.tail <- Vet[-cKs, , drop = FALSE]
      Pi <- Vet[cK, , drop = FALSE]
      tVhVt <- crossprod(V.head, V.tail)
      P <- tVhVt + t(Pi) %*% (solve(diag(length(N)) - tcrossprod(Pi), Pi) %*% tVhVt)

      Z <- rbind(t(sigma[group] * t(Uet)), matrix(NA, len + max(K) - 1, length(group)))

      for (j in (L + 1) : (L + len + max(K) - 1)) {
        Z[j, ] <- P %*% Z[j - 1, ]
      }

      R <- lapply(seq_along(N), function(idx) {
          rowSums(.hankelize.multi(Z[1 : (L + len + K[idx] - 1), , drop = FALSE],
                                   Vet[cKs[idx] : cK[idx], , drop = FALSE],
                                   fft.plan[[idx]]))
        })
    }

    out[[i]] <- if (only.new) {
      .to.series.list(lapply(seq_along(N), function(idx) R[[idx]][seq(to = N.res[idx], length.out = len)]))
    } else {
      for (idx in seq_along(N)) {
        length(R[[idx]]) <- N.res[idx]
        R[[idx]] <- .na.bind(F[[idx]], R[[idx]], update.method = "replace")
      }

      class(R) <- "series.list"
      R
    }
    out[[i]] <- .apply.attributes(x, out[[i]],
                                  fixup = TRUE,
                                  only.new = only.new, drop = drop.attributes)

  }

  names(out) <- .group.names(groups)
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
    groups <- list(1:min(nsigma(x), nu(x)))

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
    out[[i]] <- .apply.attributes(x, cbind(Value = rowMeans(bF), t(cf)),
                                  fixup = TRUE, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
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

predict.mssa <- function(object,
                         groups, len = 1,
                         method = c("recurrent-column", "recurrent-row", "vector-column", "vector-row"),
                         ...,
                         drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  method <- match.arg(method)
  dots <- list(...)

  # Calculate a forecast
  switch (method,
          'recurrent-row' = do.call(rforecast, c(list(object, direction = "row", groups = groups, len = len, drop = drop, drop.attributes = drop.attributes, cache = cache), dots)),
          'recurrent-column' = do.call(rforecast, c(list(object, direction = "column", groups = groups, len = len, drop = drop, drop.attributes = drop.attributes, cache = cache), dots)),
          'vector-row' = do.call(vforecast, c(list(object, direction = "row", groups = groups, len = len, cache = cache), dots)),
          'vector-column' = do.call(vforecast, c(list(object, direction = "column", groups = groups, len = len, cache = cache), dots)))
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
                                vector = "SSA (vector)",
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

"lrr.mssa" <- `lrr.1d.ssa`

"lrr.cssa" <- `lrr.1d.ssa`
"rforecast.cssa" <- `rforecast.1d.ssa`;
"vforecast.cssa" <- `vforecast.1d.ssa`;

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
