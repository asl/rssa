#   R package for Singular Spectrum Analysis
#   Copyright (c) 2013 Anton Korobeynikov <asl@math.spbu.ru>
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

is.complete <- function(idx, na.idx, L)
  # Do not use intersect here because it does unique and various other
  # stuff we can optimize out here
  length(na.idx[match(seq(from = idx, length.out = L), na.idx, 0)]) == 0

clplot <- function(x, ...,
                   main = "Ratio of complete lag vectors given window length",
                   xlab = "window length, L", ylab = "Ratio",
                   type = "l") {
  N <- length(x)
  na.idx <- which(is.na(x))
  # FIME: This is really, really ugly
  cr <- sapply(2:(N + 1) %/% 2,
               function(L) sum(sapply(1:(N-L+1), is.complete, na.idx = na.idx, L = L)) / (N - L + 1))
  plot(2:(N + 1) %/% 2, cr, ..., main = main, xlab = xlab, type = type)
}

# FIXME: This is ugly
complete.idx <- function(na.idx, N, L)
  which(sapply(1:(N-L+1), is.complete, na.idx = na.idx, L = L))

sproject <- function(U, X)
  U %*% (t(U) %*% X)

pi.project.complete <- function(U, v) {
  midx <- is.na(v)
  fidx <- !midx

  # Skip the processing in case of no complete values
  if (all(midx))
    return (numeric(0))

  W <- U[midx, , drop = FALSE]
  V <- U[fidx, , drop = FALSE]
  VtW <- V %*% t(W)
  # FIXME: Cache a for later use
  a <- qr(diag(nrow = sum(midx)) - W %*% t(W))
  nc <- ncol(a$qr)
  if (a$rank != nc)
    return (numeric(0))

  V %*% t(V) %*% v[fidx] + VtW %*% qr.coef(a, t(VtW) %*% v[fidx])
}

pi.project.missing <- function(U, v) {
  midx <- is.na(v)
  fidx <- !midx

  if (all(midx))
    return (numeric(0))

  W <- U[midx, , drop = FALSE]
  V <- U[fidx, , drop = FALSE]

  a <- qr(diag(nrow = sum(midx)) - W %*% t(W))
  nc <- ncol(a$qr)
  if (a$rank != nc)
    return (numeric(0))

  qr.coef(a, W %*% t(V) %*% v[fidx])
}

classify.gaps <- function(na.idx, L, N) {
  ## First, split the na.idx into separate clusters
  left  <- na.idx[which(diff(c(-L, na.idx)) > L)]
  right <- na.idx[which(diff(c(na.idx, N+L)) > L)]

  stopifnot(length(left) == length(right))
  stopifnot(all(left <= right))

  # Iterate over borders classifying the gaps
  res <- list()
  for (idx in seq_along(left)) {
      l <- left[idx]
      r <- right[idx]
      out <- list(left = l, right = r)
      out$pos <- if (l < L) "left" else if (r >= N-L+1) "right" else "internal"
      out$kind <- if (length(intersect(na.idx, seq(from = l, to = r))) == r - l + 1) "dense" else "sparse"
      res[[idx]] <- out
  }

  res
}

gapfill.1d.ssa <- function(x, groups,
                           base = c("reconstructed", "original"),
                           method = c("sequential", "simultaneous"),
                           alpha = 0.5,
                           ...,
                           drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  method <- match.arg(method)
  base <- match.arg(base)

  if (!is.shaped(x))
    stop("gapfilling should start from shaped SSA object")

  if (alpha < 0 || alpha > 1)
    stop("`alpha' should be between 0 and 1")

  L <- x$window; N <- x$length; K <- N - L + 1

  # Grab the reconstructed series if we're basing on them
  if (identical(base, "reconstructed"))
    r <- reconstruct(x, groups = groups, ..., cache = cache)

  out <- list()
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    F <- if (identical(base, "reconstructed")) r[[i]] else .F(x)
    na.idx <- which(is.na(F))

    Ug <- .colspan(x, group)

    if (identical(method, "simultaneous")) {
      Xr <- hankel(F, L)

      complete <- complete.idx(na.idx, N, L)
      incomplete <- setdiff(1:K, complete)

      ## Project the complete trajectory vectors onto subspace
      Xr[, complete] <- sproject(Ug, Xr[, complete, drop = FALSE])

      ## Now iterate over incomplete vectors building the projections for
      ## entries via so-called Pi-projection operator
      for (vidx in incomplete) {
        v <- Xr[, vidx]
        fidx <- !is.na(v)
        midx <- is.na(v)

        if (length(resf <- pi.project.complete(Ug, v)) == 0)
                next
        if (length(resm <- pi.project.missing(Ug, v)) == 0)
                next

        Xr[fidx, vidx] <- resf
        Xr[midx, vidx] <- resm
      }
      out[[i]] <- hankel(Xr)
    } else {
      res <- F
      gaps <- classify.gaps(na.idx, L, N)
      blrr <- lrr(Ug, direction = "backward")
      flrr <- lrr(Ug, direction = "forward")

      for (gap in gaps) {
        leftpos <- gap$left; rightpos <- gap$right; len <- rightpos - leftpos + 1
        to.fill <- seq.int(leftpos, rightpos)
        if (identical(gap$pos, "left")) {
          res[to.fill] <- apply.lrr(F[seq.int(from = rightpos + 1, length.out = L)], blrr,
                                    len = len, only.new = TRUE, direction = "backward")
        } else if (identical(gap$pos, "right")) {
          res[to.fill] <- apply.lrr(F[seq.int(to = leftpos - 1, length.out = L)], flrr,
                                    len = len, only.new = TRUE, direction = "forward")
        } else {
          stopifnot(identical(gap$pos, "internal"))

          res[to.fill] <-
                  (alpha *
                   apply.lrr(F[seq.int(from = rightpos + 1, length.out = L)], blrr,
                             len = len, only.new = TRUE, direction = "backward")
                   +
                   (1 - alpha) *
                   apply.lrr(F[seq.int(to = leftpos - 1, length.out = L)], flrr,
                             len = len, only.new = TRUE, direction = "forward"))
        }
      }
    }
    out[[i]] <- .apply.attributes(x, res, fixup = FALSE, only.new = FALSE, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

gapfill <- function(x, ...)
  UseMethod("gapfill")
