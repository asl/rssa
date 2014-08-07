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
                            fill = NULL, eps = 1e-6, numiter = 0,
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

    rss <- max((F-rF)^2)
    if (trace) cat(sprintf("RSS(%d): %s\n", it, paste0(rss, collapse = " ")))
    it <- it + 1
    if ((numiter > 0 && it >= numiter) || rss < eps)
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

igapfill.2d.ssa <- igapfill.1d.ssa

igapfill.ssa <- function(x, ...)
  stop("iterative gap filling is not available for this kind of SSA yet")

igapfill <- function(x, ...)
  UseMethod("igapfill")
