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
  if (all(midx)) {
#    print(vidx)
    return (numeric(0))
  }
  
  W <- U[midx, , drop = FALSE]
  V <- U[fidx, , drop = FALSE]
  VtW <- V %*% t(W)
  # FIXME: Cache a for later use
  a <- qr(diag(nrow = sum(midx)) - W %*% t(W))
  nc <- ncol(a$qr)
  if (a$rank != nc) {
#    print(vidx)
#    print(c(a$rank, nc))
    return (numeric(0))
  }

  V %*% t(V) %*% v[fidx] + VtW %*% qr.coef(a, t(VtW) %*% v[fidx])
}

pi.project.missing <- function(U, v) {
  midx <- is.na(v)
  fidx <- !midx

  if (all(midx)) {
    # print(vidx)
    return (numeric(0))
  }
#    stopifnot(any(fidx))
    
  W <- U[midx, , drop = FALSE]
  V <- U[fidx, , drop = FALSE]

  a <- qr(diag(nrow = sum(midx)) - W %*% t(W))
  nc <- ncol(a$qr)
  if (a$rank != nc) {
#    print(vidx)
#    print(c(a$rank, nc))
    return (numeric(0))
  }

  qr.coef(a, W %*% t(V) %*% v[fidx])
}

classify.gaps <- function(na.idx, L, N) {
  ## First, split the na.idx into separate clusters
  d <- which(diff(na.idx) > L)
  cl.borders <- if (length(na.idx) == 1) c(na.idx, na.idx) else sort(unique(c(min(na.idx), na.idx[d], na.idx[d+1], max(na.idx))))
  stopifnot(length(cl.borders) %% 2 == 0)
  left <- cl.borders[c(TRUE, FALSE)]
  right <- cl.borders[c(FALSE, TRUE)]
  
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

fill.in.left <- function(Xr, blrr) {
  leftpos <- 1
  rightpos <- ncol(Xr) - 1
  len <- rightpos - leftpos + 1
  L <- nrow(Xr)

  f <- apply.lrr(Xr[, rightpos + 1], blrr,
                 len = len,
                 direction = "backward", only.new = TRUE)

  for (vidx in seq(rightpos, leftpos)) {
    v <- Xr[, vidx]
    # First, fill in full entries
    if (rightpos - vidx < L - 1) {
      fidx <- seq.int(rightpos - vidx + 1 + 1, L)
      Xr[fidx, vidx] <- Xr[, vidx + 1][fidx - 1]
    }

    # Now, add missed entries from the LRR prediction
    midx <- seq.int(1, min(L, rightpos - vidx + 1))
    mfidx <- seq.int(from = vidx, length.out = length(midx))
    Xr[midx, vidx] <- f[mfidx]
  }

  Xr
}

fill.in.right <- function(Xr, flrr, L) {
  leftpos <- 2
  rightpos <- ncol(Xr)
  len <- rightpos - leftpos + 1
  L <- nrow(Xr)

  f <- apply.lrr(Xr[, leftpos - 1], flrr,
                 len = len,
                 direction = "forward", only.new = TRUE)
  for (vidx in seq(leftpos, rightpos)) {
    v <- Xr[, vidx]

    # First, fill in full entries
    if (vidx <= L) { 
      fidx <- seq.int(1, L - vidx + 1)
      Xr[fidx, vidx] <- Xr[, vidx - 1][fidx + 1]
    }

    # Now, add missed entries from the LRR prediction
    midx <- seq.int(to = L, length.out = min(L, vidx - 1))
    mfidx <- seq.int(to = vidx - 1, length.out = length(midx))
    Xr[midx, vidx] <- f[mfidx]
  }

  Xr
}

gapfill <- function(F, groups, L = (N + 1) %/% 2,
                    base = c("reconstructed", "original"),
                    method = c("stepwise", "simultaneous"),
                    ...,
                    drop.attributes = FALSE) {
  method <- match.arg(method)
  base <- match.arg(base)
  N <- length(F)
  K <- N - L + 1

  groups <- unlist(groups)

  s <- suppressWarnings(ssa(F, L = L, kind = "1d-ssa"))
 
  # Grab the reconstructed series if we're basing on them
  if (identical(base, "reconstructed")) {
    r <- reconstruct(s, groups = list(groups), ...)
    stopifnot(length(r) == 1)
    F <- r[[1]]
  }

  na.idx <- which(is.na(F))
  complete <- complete.idx(na.idx, N, L)
  incomplete <- setdiff(1:K, complete)
 
  #Ug <- s$U[, groups, drop = FALSE]
  Ug <- .colspan(s, groups)

  ## Project the complete trajectory vectors onto subspace
  Xr <- hankel(F, L)
  Xr[, complete] <- sproject(Ug, Xr[, complete, drop = FALSE])

  if (identical(method, "simultaneous")) {
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
  } else {
    gaps <- classify.gaps(na.idx, L, N)
    blrr <- lrr(Ug, direction = "backward")
    flrr <- lrr(Ug, direction = "forward")

    for (gap in gaps) {
      if (identical(gap$pos, "left")) {
        leftpos <- gap$left
        rightpos <- gap$right
        to.fill <- seq(leftpos, rightpos + 1)
        len <- rightpos - leftpos + 1
        Xr[, to.fill] <- fill.in.left(Xr[, to.fill], blrr)

        ## Additional handling for non-border left gaps - fill in and project
        if (leftpos > 1) {
          cidx <- seq(leftpos - 1, 1)
          for (vidx in cidx) {
            fidx <- seq(leftpos - vidx + 1, L)
            Xr[fidx, vidx] <- Xr[, vidx + 1][fidx - 1]
          }
          Xr[, cidx] <- sproject(Ug, Xr[, cidx])
        }
      } else if (identical(gap$pos, "right")) {
        leftpos <- gap$left - L + 1
        rightpos <- gap$right - L + 1
        to.fill <- seq(leftpos - 1, rightpos)
        len <- rightpos - leftpos + 1
        Xr[, to.fill] <- fill.in.right(Xr[, to.fill], flrr)

        ## Additional handling for non-border right gaps - fill in and project
        if (rightpos < K) {
          cidx <- seq(rightpos + 1, K)
          for (vidx in cidx) {
            fidx <- seq(1, L - (vidx - rightpos))
            Xr[fidx, vidx] <- Xr[, vidx - 1][fidx + 1]
          }
          Xr[, cidx] <- sproject(Ug, Xr[, cidx])
        }

      } else {
        stopifnot(identical(gap$pos, "internal"))

        ## Fill in in left direction
        leftpos <- gap$left - L + 1
        rightpos <- gap$right
        to.fill <- seq(leftpos, rightpos + 1)
        Xr[, to.fill] <- fill.in.left(Xr[, to.fill], blrr)

        ## Fill in in right direction
        leftpos <- gap$left - L + 1
        rightpos <- gap$right
        to.fill <- seq(leftpos - 1, rightpos)
        Xr[, to.fill] <- (Xr[, to.fill] + fill.in.right(Xr[, to.fill], flrr)) / 2
#        Xr[, to.fill] <-  fill.in.right(Xr[, to.fill], flrr)
      }
    }
  }

  .apply.attributes(s, hankel(Xr), fixup = FALSE, drop = drop.attributes)
}
