#   R package for Singular Spectrum Analysis
#   Copyright (c) 2014 Alex Shlemov <shlemovalex@gmail.com>
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


orthopoly <- function(d, L) {
  if (is.character(d)) {
    d <- match.arg(d,
                   choices = c("none", "constant", "centering", "linear", "quadratic", "qubic"))
    d <- switch(d,
                none = 0, constant =, centering = 1, linear = 2, quadratic = 3, qubic = 4)
  }

  stopifnot(is.numeric(d) && length(d) == 1)

  # Check dimension
  if (length(L) > 1 && d > 1) {
    stop("Polynomial projection is not implemented for multidimensional SSA yet")
  }

  L <- prod(L)  # TODO Think about MSSA case

  if (d == 0) {
    # Return matrix with zero columns
    matrix(NA_real_, L, 0)
  } else if (d == 1) {
    matrix(1 / sqrt(L), L, 1)
  } else {
    cbind(1 / sqrt(L), poly(seq_len(L), d - 1))
  }
}

.phmat <- function(x) {
  hmat <- .get.or.create.trajmat(x)
  column.projector <- .get(x, "column.projector")
  row.projector <- .get(x, "row.projector")

  matmul <- function(v) {
    v <- v - row.projector %*% crossprod(row.projector, v)
    v <- hmatmul(hmat, v, transposed = FALSE)
    v <- v - column.projector %*% crossprod(column.projector, v)

    v
  }

  tmatmul <- function(v) {
    v <- v - column.projector %*% crossprod(column.projector, v)
    v <- hmatmul(hmat, v, transposed = TRUE)
    v <- v - row.projector %*% crossprod(row.projector, v)

    v
  }

  extmat(matmul, tmatmul, hrows(hmat), hcols(hmat))
}

.get.or.create.phmat <- function(x)
  .get.or.create(x, "phmat", .phmat(x))

.calc.projections <- function(x, force.update = FALSE) {
  if (!is.null(.decomposition(x)) && !force.update)
    return(x)

  hmat <- .get.or.create.trajmat(x)
  LU <- .get(x, "column.projector")
  RV <- .get(x, "row.projector")

  RU <- hmat %*% RV
  LV <- crossprod(hmat, LU) - RV %*% crossprod(RU, LU)
  Rsigma <- sqrt(colSums(RU^2))
  RU <- RU / rep(Rsigma, each = nrow(RU))
  Lsigma <- sqrt(colSums(LV^2))
  LV <- LV / rep(Lsigma, each = nrow(LV))

  .set.decomposition(x,
                     nPR = ncol(RV), nPL = ncol(LU),
                     sigma = c(Rsigma, Lsigma), U = cbind(RU, LU), V = cbind(RV, LV))

  x
}

nspecial.pssa <- function(x) {
  sum(unlist(.decomposition(x, c("nPL", "nPR"))))
}

decompose.pssa <- function(x,
                           neig = NULL,
                           ...,
                           force.continue = FALSE) {
  ## Check, whether continuation of decomposition is requested
  ## FIXME: Check the caps
  if (!force.continue && nsigma(x) > nspecial(x) &&
       !identical(x$svd.method, "nutrlan"))
    stop("Continuation of decomposition is not yet implemented for this method.")

  if (is.null(neig))
    neig <- .default.neig(x, ...)

  # Compute special eigentriples if needed
  .calc.projections(x)

  nspecial <- nspecial(x)

  # Extract special components
  ssigma <- .sigma(x)[seq_len(nspecial)]
  sU <- .U(x)[, seq_len(nspecial), drop = FALSE]
  sV <- .V(x)[, seq_len(nspecial), drop = FALSE]

  if (identical(x$svd.method, "svd")) {
    S <- svd(as.matrix(.get.or.create.phmat(x)), nu = neig, nv = neig)
    .set.decomposition(x,
                       sigma = c(ssigma, S$d), U = cbind(sU, S$u), V = cbind(sV, S$v))
  } else if (identical(x$svd.method, "eigen")) {
    S <- eigen(tcrossprod(.get.or.create.phmat(x)), symmetric = TRUE)

    ## Fix small negative values
    S$values[S$values < 0] <- 0

    .set.decomposition(x,
                       sigma = c(ssigma, sqrt(S$values[seq_len(neig)])),
                       U = cbind(sU, S$vectors[, seq_len(neig), drop = FALSE]),
                       V = sV)
  } else if (identical(x$svd.method, "propack")) {
    S <- propack.svd(.get.or.create.phmat(x), neig = neig, ...)
    ## Save results
    .set.decomposition(x,
                       sigma = c(ssigma, S$d), U = cbind(sU, S$u), V = cbind(sV, S$v))
  } else if (identical(x$svd.method, "nutrlan")) {
    S <- trlan.svd(.get.or.create.phmat(x), neig = neig, ...,
                   lambda = .sigma(x)[-seq_len(nspecial)], U = .U(x)[, -seq_len(nspecial), drop = FALSE])
    .set.decomposition(x,
                       sigma = c(ssigma, S$d), U = cbind(sU, S$u), V = cbind(sV, S$v))
  } else
    stop("unsupported SVD method")

  x
}

calc.v.pssa <- function(x, idx, ...) {
  nV <- nv(x)

  V <- matrix(NA_real_, .traj.dim(x)[2], length(idx))
  idx.old <- idx[idx <= nV]
  idx.new <- idx[idx > nV]

  if (length(idx.old) > 0) {
    V[, idx <= nV] <- .V(x)[, idx.old]
  }

  if (length(idx.new) > 0) {
    sigma <- .sigma(x)[idx.new]

    if (any(sigma <= .Machine$double.eps)) {
      sigma[sigma <= .Machine$double.eps] <- Inf
      warning("some sigmas are equal to zero. The corresponding vectors will be zeroed")
    }

    U <- .U(x)[, idx.new, drop = FALSE]

    h <- .get.or.create.phmat(x)
    V[, idx > nV] <- crossprod(h, U) / rep(sigma, each = nrow(V))
  }

  invisible(V)
}

.colspan.pssa <- function(x, idx) {
  decomposition <- .decomposition(x, c("U", "nPR"))

  span <- decomposition$U[, idx, drop = FALSE]

  # Perform orthogonalization only if it is really needed
  if (any(idx <= decomposition$nPR)) {
    # TODO It can be implemented little more efficiently
    span <- qr.Q(qr(span))
  }

  span
}

.rowspan.pssa <- function(x, idx) {
  decomposition <- .decomposition(x, c("nPR", "nPR"))

  span <- calc.v(x, idx)

  # Perform orthogonalization only if it is really needed
  if (any((idx > decomposition$nPR) & (idx <= decomposition$nPR + decomposition$nPL))) {
    # TODO It can be implemented little more efficiently
    span <- qr.Q(qr(span))
  }

  span
}

enlarge.basis <- function(B, len, ...) {
  N <- nrow(B)

  if (ncol(B) == 0) {
    return(matrix(NA_real_, N + len, 0))
  }

  P <- shift.matrix(B, ...)

  B <- rbind(B, matrix(NA, len, ncol(B)))
  for (i in seq_len(len)) {
    B[N + i, ] <- B[N + i - 1, ] %*% P
  }

  B
}

rforecast.pssa.1d.ssa <- function(x, groups, len = 1,
                                  base = c("reconstructed", "original"),
                                  only.new = TRUE,
                                  reverse = FALSE,
                                  ...,
                                  drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  if (is.shaped(x))
    stop("`forecasting is not implemented for shaped SSA case yet")

  L <- x$window
  K <- x$length - L + 1

  base <- match.arg(base)
  if (missing(groups))
    groups <- as.list(seq_len(min(nsigma(x), nu(x))))

  # Continue decomposition, if necessary
  desired <- .maybe.continue(x, groups = groups, ...)

  # Grab the reconstructed series if we're basing on them
  if (identical(base, "reconstructed"))
    r <- reconstruct(x, groups = groups, ..., cache = cache)

  right.special.triples <- seq_len(.decomposition(x, "nPR"))
  right.groups <- lapply(groups, function(group) intersect(group, right.special.triples))
  nonright.groups <- lapply(groups, function(group) setdiff(group, right.special.triples))

  # Calculate the LRR corresponding to groups
  lf <- lrr(x, groups = nonright.groups, reverse = reverse, drop = FALSE)
  stopifnot(length(lf) == length(groups))

  rlf <- lapply(right.groups,
                function(group) {
                  lrr.default(.rowspan(x, group),
                              reverse = reverse,
                              orthonormalize = FALSE)
                })
  stopifnot(length(rlf) == length(groups))

  rlf.all <- lrr.default(.rowspan(x, right.special.triples),
                         reverse = reverse,
                         orthonormalize = FALSE)

  sigma <- .sigma(x)
  U <- .U(x)
  V <- if (nv(x) >= desired) .V(x) else NULL

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    right.group <- if (identical(base, "reconstructed")) right.groups[[i]] else right.special.triples
    right.lrr <- if (identical(base, "reconstructed")) rlf[[i]] else rlf.all

    # Calculate drifts
    Uet <- U[, right.group, drop = FALSE]
    Vet <- if (is.null(V)) calc.v(x, idx = right.group) else V[, right.group, drop = FALSE]
    drift <- (((if (!reverse) c(-lf[[i]], 1) else c(1, -lf[[i]])) %*% Uet) *
              sigma[right.group]) %*% t(Vet)
    drift <- apply.lrr(drift, right.lrr,
                       reverse = reverse,
                       len = len, only.new = TRUE)

    # Calculate the forecasted values
    out[[i]] <- apply.lrr(if (identical(base, "reconstructed")) r[[i]] else .get(x, "F"),
                          lf[[i]],
                          reverse = reverse,
                          len, only.new = only.new, drift = drift)
    out[[i]] <- .apply.attributes(x, out[[i]],
                                  fixup = TRUE,
                                  reverse = reverse,
                                  only.new = only.new, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  # Forecasted series can be pretty huge...
  invisible(out)
}

vforecast.pssa.1d.ssa <- function(x, groups, len = 1,
                                  only.new = TRUE,
                                  ...,
                                  drop = TRUE, drop.attributes = FALSE) {
  if (is.shaped(x))
    stop("`forecasting is not implemented for shaped SSA case yet")

  L <- x$window
  K <- x$length - L + 1
  N <- K + L - 1 + len + L - 1
  N.res <- K + L - 1 + len

  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))

  groups <- lapply(groups, unique)

  # Continue decomposition, if necessary
  desired <- .maybe.continue(x, groups = groups, ...)

  right.special.triples <- seq_len(.decomposition(x, "nPR"))
  right.groups <- lapply(groups, function(group) intersect(group, right.special.triples))
  nonright.groups <- lapply(groups, function(group) setdiff(group, right.special.triples))

  sigma <- .sigma(x)
  U <- .U(x)
  V <- if (nv(x) >= desired) .V(x) else NULL

  # Grab the FFT plan
  fft.plan <- fft.plan.1d(N, L = L)

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    right.group <- right.groups[[i]]
    nonright.group <- nonright.groups[[i]]

    Uet.right <- U[, right.group, drop = FALSE]
    Uet.nonright <- U[, nonright.group, drop = FALSE]
    Vet.right <- if (is.null(V)) calc.v(x, idx = right.group) else V[, right.group, drop = FALSE]
    Vet.nonright <- if (is.null(V)) calc.v(x, idx = nonright.group) else V[, nonright.group, drop = FALSE]

    Z.right <- Vet.right * rep(sigma[right.group], each = nrow(Vet.right))
    Z.nonright <- Vet.nonright * rep(sigma[nonright.group], each = nrow(Vet.nonright))

    Pright <- t(shift.matrix(Z.right))
    Pnonright <- shift.matrix(Uet.nonright)

    PP <- qr.solve(Uet.nonright[-L,, drop = FALSE],
                   Uet.right[-1,, drop = FALSE] - Uet.right[-L,, drop = FALSE] %*% Pright)

    Uet <- cbind(Uet.right, Uet.nonright)
    Z <- rbind(cbind(Z.right, Z.nonright), matrix(NA, len + L - 1, length(group)))

    P <- rbind(cbind(Pright, matrix(0, length(right.group), length(nonright.group))),
               cbind(PP, Pnonright))

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

.default.neig.pssa <- function(x, ...) {
  nPR <- max(0, ncol(.get(x, "column.projector")))
  nPL <- max(0, ncol(.get(x, "row.projector")))

  tjdim <- .traj.dim(x)

  min(50, tjdim - max(nPR, nPL))
}

.init.fragment.pssa <- function(this)
  expression({
    ## First, initialize the main object
    ## We cannot use NextMethod here due to non-standard evaluation
    class.wo.pssa <- class(this)[!grepl("^pssa", class(this))]
    eval(getS3method(".init.fragment", class.wo.pssa)(this))
    ## eval(.init.fragment.1d.ssa(this))

    # unwind all dimentions except last
    .unwind.dim.except.last <- function(x) {
      d <- dim(x)
      if (is.null(d) || length(d) <= 2) {
        return(x)
      }

      dim(x) <- c(prod(d[-length(d)]), d[length(d)])

      x
    }

    column.projector <- .unwind.dim.except.last(column.projector)
    row.projector <- .unwind.dim.except.last(row.projector)

    ## Next, calculate the projectors
    column.projector <- if (length(column.projector) == 1) orthopoly(column.projector, L) else qr.Q(qr(column.projector))
    row.projector <- if (length(row.projector) == 1) orthopoly(row.projector, K) else qr.Q(qr(row.projector))

    ## Check projector dimensions
    ## TODO Think about MSSA case
    stopifnot(nrow(column.projector) == prod(L))
    stopifnot(nrow(row.projector) == prod(K))

    ## Shape projectors if needed
    if (!is.null(wmask)) {
      column.projector <- column.projector[as.vector(wmask),, drop = FALSE]
      column.projector <- qr.Q(qr(column.projector))
    }
    if (!is.null(fmask)) {
      row.projector <- row.projector[as.vector(fmask),, drop = FALSE]
      row.projector <- qr.Q(qr(row.projector))
    }
  })
