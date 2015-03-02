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

num.complete <- function(N, L, na.idx) {
  left  <- na.idx[which(diff(c(-L, na.idx)) > L)]
  right <- na.idx[which(diff(c(na.idx, N + L + 1)) > L)]

  cl.left <- c(1, right + 1)
  cl.right <- c(left - 1, N)

  num <- cl.right - cl.left + 1 - L + 1
  sum(ifelse(num > 0, num, 0))
}

clplot <- function(x, ...) {
  N <- length(x)
  na.idx <- which(is.na(x))
  cr <- sapply(2:((N + 1) %/% 2),
               function(L) num.complete(N, L = L, na.idx = na.idx) / (N - L + 1) * 100)

  dots <- list(...)
  dots <- modifyList(list(main = "Proportion of complete lag vectors",
                          xlab = "window length, L", ylab = "Percents",
                          type = "l"), dots)
  do.call(plot, c(list(2:((N + 1) %/% 2), cr), dots))
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
  right <- na.idx[which(diff(c(na.idx, N + L + 1)) > L)]

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

summarize.gaps.default <- function(x, L) {
  na.idx <- which(is.na(x))
  N <- length(x)
  res <- list()

  L.range <- if (missing(L) || is.null(L)) 1:((N + 1) %/% 2) else L
  for (L in L.range)
    res[[L]] <- classify.gaps(na.idx, L, N)

  res$N <- N
  res$L.range <- L.range
  res$na.idx <- na.idx
  res$call <- match.call()

  class(res) <- "ssa.gaps"
  res
}

summarize.gaps.1d.ssa <- summarize.gaps.toeplitz.ssa <- summarize.gaps.cssa <- function(x, L = NULL) {
  summarize.gaps.default(x$F, L = L)
}

summarize.gaps.ssa <- function(x, L) {
  stop("this function is not available for this SSA type")
}

summarize.gaps <- function(x, L)
  UseMethod("summarize.gaps")

print.ssa.gaps <- function(x, ...) {
  N <- x$N
  L.range <- x$L.range
  cat("\nCall:\n", deparse(x$call), "\n\n", sep="");
  cat("Gaps summary:\n")
  nogaps <- TRUE
  for (i in L.range) {
    gaps <- x[[i]]
    if (is.null(gaps) || length(gaps) == 0)
      next
    nogaps <- FALSE
    cat("L =", i, "\n")
    for (gap in gaps)
      cat("  [", gap$left, ", ", gap$right, "], ", gap$pos, ", ", gap$kind, "\n", sep = "")
  }
  if (nogaps)
    cat("  no gaps\n")
}

plot.ssa.gaps <- function(x, ...) {
  N <- x$N
  na.idx <- x$na.idx

  cr <- sapply(2:((N + 1) %/% 2),
               function(L) num.complete(N, L = L, na.idx = na.idx) / (N - L + 1) * 100)

  dots <- list(...)
  dots <- modifyList(list(main = "Proportion of complete lag vectors",
                          xlab = "window length, L", ylab = "Percents",
                          type = "l"), dots)
  do.call(plot, c(list(2:((N + 1) %/% 2), cr), dots))
}

.fill.in <- function(F, L, gap, flrr, blrr, alpha) {
  leftpos <- gap$left; rightpos <- gap$right; len <- rightpos - leftpos + 1
  to.fill <- seq.int(leftpos, rightpos)
  if (identical(gap$pos, "left")) {
     apply.lrr(F[seq.int(from = rightpos + 1, length.out = L)], blrr,
               len = len, only.new = TRUE, reverse = TRUE)
   } else if (identical(gap$pos, "right")) {
     apply.lrr(F[seq.int(to = leftpos - 1, length.out = L)], flrr,
               len = len, only.new = TRUE, reverse = FALSE)
   } else {
     stopifnot(identical(gap$pos, "internal"))

     (alpha *
      apply.lrr(F[seq.int(from = rightpos + 1, length.out = L)], blrr,
                len = len, only.new = TRUE, reverse = TRUE)
      +
      (1 - alpha) *
      apply.lrr(F[seq.int(to = leftpos - 1, length.out = L)], flrr,
                len = len, only.new = TRUE, reverse = FALSE))
   }
}

gapfill.1d.ssa <- function(x, groups,
                           base = c("original", "reconstructed"),
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
      res <- hankel(Xr)
    } else {
      res <- F
      gaps <- classify.gaps(na.idx, L, N)
      blrr <- lrr(Ug, reverse = TRUE)
      flrr <- lrr(Ug, reverse = FALSE)

      for (gap in gaps) {
        to.fill <- seq.int(gap$left, gap$right)
        res[to.fill] <- .fill.in(F, L, gap, flrr, blrr, alpha)
      }
    }
    out[[i]] <- .apply.attributes(x, res, fixup = FALSE, only.new = FALSE, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

gapfill.mssa <- function(x, groups,
                         base = c("original", "reconstructed"),
                         alpha = 0.5,
                         ...,
                         drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
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
    F <- if (identical(base, "reconstructed")) .to.series.list(r[[i]], na.rm = FALSE) else .F(x)

    Ug <- .colspan(x, group)
    blrr <- lrr(Ug, reverse = TRUE)
    flrr <- lrr(Ug, reverse = FALSE)

    res <- F
    for (idx in seq_along(N)) {
      na.idx <- which(is.na(F[[idx]]))
      gaps <- classify.gaps(na.idx, L, N[[idx]])

      for (gap in gaps) {
        to.fill <- seq.int(gap$left, gap$right)
        res[[idx]][to.fill] <- .fill.in(F[[idx]], L, gap, flrr, blrr, alpha)
      }
    }
    out[[i]] <- .apply.attributes(x, res, fixup = FALSE, only.new = FALSE, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

gapfill.cssa <- gapfill.toeplitz.ssa <- gapfill.1d.ssa

gapfill <- function(x, ...)
  UseMethod("gapfill")
