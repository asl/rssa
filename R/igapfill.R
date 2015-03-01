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

igapfill.1d.ssa <- function(x,
                            groups,
                            fill = NULL, tol = 1e-6, maxiter = 0,
                            norm = function(x) sqrt(max(x^2)),
                            base = c("original", "reconstructed"),
                            ...,
                            trace = FALSE,
                            drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  base <- match.arg(base)
  N <- x$length

  if (!is.shaped(x))
    stop("gapfilling should start from shaped SSA object")

  ## Obtain the initial approximation
  ugroups <- seq_len(max(unique(unlist(groups))))
  if (identical(base, "reconstructed")) {
    r <- reconstruct(x, groups = list(ugroups), ..., cache = cache)
    stopifnot(length(r) == 1)
    F <- r[[1]]
  } else {
    F <- .F(x)
  }

  ## Determine the indices of missing values
  na.idx <- which(is.na(F))

  ## Obtain the initial approximation
  if (is.null(fill)) fill <- mean(F, na.rm = TRUE)
  ## Check filler for sanity
  if (length(fill) > 1 &&
      (length(fill) != length(F) ||
       any(is.na(fill[na.idx]))))
      stop("filler should have the same shape as series and provide non-NA filled values")
  F[na.idx] <- if (length(fill) > 1) fill[na.idx] else fill

  # Do the actual iterations until the convergence (or stoppping due to number
  # of iterations)
  it <- 0
  scall <- x$ecall

  repeat {
    scall$x <- F
    s <- eval(scall)

    r <- reconstruct(s, groups = list(ugroups), ..., cache = cache)
    stopifnot(length(r) == 1)
    rF <- F
    rF[na.idx] <- r[[1]][na.idx]

    rss <- norm(F - rF)
    if (trace) cat(sprintf("RSS(%d): %s\n", it, paste0(rss, collapse = " ")))
    it <- it + 1
    if ((maxiter > 0 && it >= maxiter) || rss < tol)
      break
    F <- rF
  }

  scall$x <- F
  s <- eval(scall)
  r <- reconstruct(s, groups = groups, ..., drop.attributes = drop.attributes, cache = cache)

  out <- list()
  for (i in seq_along(r)) {
    if (identical(base, "reconstructed")) {
      out[[i]] <- r[[i]]
    } else {
      out[[i]] <- .F(x)
      out[[i]][na.idx] <- r[[i]][na.idx]
    }
    out[[i]] <- .apply.attributes(x, out[[i]], fixup = FALSE, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

igapfill.nd.ssa <- igapfill.cssa <- igapfill.toeplitz.ssa <- igapfill.1d.ssa

igapfill.mssa <- function(x,
                          groups,
                          fill = NULL, tol = 1e-6, maxiter = 0,
                          norm = function(x) sqrt(max(x^2)),
                          base = c("original", "reconstructed"),
                          ...,
                          trace = FALSE,
                          drop = TRUE, drop.attributes = FALSE, cache = TRUE) {
  base <- match.arg(base)
  N <- x$length

  if (!is.shaped(x))
    stop("gapfilling should start from shaped SSA object")

  ## Obtain the initial approximation
  ugroups <- seq_len(max(unique(unlist(groups))))
  if (identical(base, "reconstructed")) {
    r <- reconstruct(x, groups = list(ugroups), ..., cache = cache)
    stopifnot(length(r) == 1)
    F <- .to.series.list(r[[1]], na.rm = FALSE)
  } else {
    F <- .F(x)
  }

  ## Determine the indices of missing values
  na.idx <- lapply(F, function(series) which(is.na(series)))
  for (idx in seq_along(F)) {
    cna.idx <- na.idx[[idx]]

    ## Obtain the initial approximation
    if (is.null(fill)) fill <- mean(F[[idx]], na.rm = TRUE)
    ## Check filler for sanity
    if (is.list(fill) &&
        (length(fill[[idx]]) != length(F[[idx]]) ||
         any(is.na(fill[[idx]][cna.idx]))))
      stop("filler should have the same shape as series and provide non-NA filled values")

    F[[idx]][cna.idx] <- if (is.list(fill)) fill[[idx]][cna.idx] else fill
  }

  # Do the actual iterations until the convergence (or stoppping due to number
  # of iterations)
  it <- 0
  scall <- x$ecall

  repeat {
    scall$x <- .apply.attributes(s, F, fixup = FALSE, drop = FALSE)
    s <- eval(scall)

    r <- reconstruct(s, groups = list(ugroups), ..., cache = cache)
    stopifnot(length(r) == 1)
    rF <- F
    for (idx in seq_along(F)) {
      cna.idx <- na.idx[[idx]]
      rF[[idx]][cna.idx] <- r[[1]][cna.idx]
    }
    rss <- norm(unlist(F) - unlist(rF))
    if (trace) cat(sprintf("RSS(%d): %s\n", it, paste0(rss, collapse = " ")))
    it <- it + 1
    if ((maxiter > 0 && it >= maxiter) || rss < tol)
      break
    F <- rF
  }

  scall$x <- .apply.attributes(s, F, fixup = FALSE, drop = FALSE)
  s <- eval(scall)
  r <- reconstruct(s, groups = groups, ..., drop.attributes = drop.attributes, cache = cache)

  out <- list()
  for (i in seq_along(r)) {
    cr <- .to.series.list(r[[i]], na.rm = FALSE)
    if (identical(base, "reconstructed")) {
      out[[i]] <- cr
    } else {
      out[[i]] <- .F(x)
      for (idx in seq_along(F)) {
        cna.idx <- na.idx[[idx]]
        out[[i]][[idx]][cna.idx] <- cr[[idx]][cna.idx]
      }
    }
    out[[i]] <- .apply.attributes(x, out[[i]], fixup = FALSE, drop = drop.attributes)
  }

  names(out) <- .group.names(groups)
  if (length(out) == 1 && drop)
    out <- out[[1]]

  out
}

igapfill.ssa <- function(x, ...)
  stop("iterative gap filling is not available for this kind of SSA yet")

igapfill <- function(x, ...)
  UseMethod("igapfill")
