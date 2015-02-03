#   R package for Singular Spectrum Analysis
#   Copyright (c) 2015 Alex Shlemov <shlemovalex@gmail.com>
#   Copyright (c) 2015 Nina Golyandina <nina@gistatgroup.com>
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


pgram <- function(x) {
  X <- fft(x)
  n <- length(x)

  if (n %% 2 == 0) {
    N <- n %/% 2 + 1
    pgram <- abs(X[seq_len(N)])^2
    if (N > 2) pgram[2:(N-1)] <- 2 * pgram[2:(N-1)]
  } else {
    N <- (n + 1) %/% 2
    pgram <- abs(X[seq_len(N)])^2
    if (N >= 2) pgram[2:N] <- 2 * pgram[2:N]
  }

  pgram <- pgram / sum(pgram)
  freqs <- seq(0, 1, length.out = n + 1)[seq_len(N)]

  list(pgram = pgram, freqs = freqs)
}

grouping.auto <- function(x, ...)
  UseMethod("grouping.auto")

grouping.auto.ssa <- function(x, ...)
  stop("grouping.auto is not implemented for this kind of SSA yet")

grouping.auto.1d.ssa <- function(x, groups,
                                 base = c("series", "vectors"),
                                 vectors = c("eigen", "factor"),
                                 freq.bins = 2,
                                 threshold = 0,
                                 drop = TRUE,
                                 ...) {
  base <- match.arg(base)
  vectors <- match.arg(vectors)

  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  groups <- sort(unique(unlist(groups)))

  if (identical(base, "vectors")) {
    if (identical(vectors, "eigen")) {
      Fs <- lapply(groups, function(i) .U(x)[, i])
    } else if (identical(vectors, "factor")) {
      Fs <- lapply(groups, function(i) calc.v(x, i))
    }
  } else if (identical(base, "series")) {
    Fs <- reconstruct(x, groups = as.list(groups))
  }

  pgs <- lapply(Fs, pgram)

  if (!is.list(freq.bins)) {
    if (length(freq.bins) == 1 && freq.bins >= 2) {
      freq.bins <- seq(0, 0.5, length.out = freq.bins + 1)[-1]
    }

    freq.lower.bound <- c(-Inf, head(freq.bins, -1))
    freq.upper.bound <- freq.bins
  } else {
    freq.bins <- lapply(freq.bins, function(x) if (length(x) == 1) c(-Inf, x) else x)
    freq.lower.bound <- lapply(freq.bins, function(x) x[1])
    freq.upper.bound <- lapply(freq.bins, function(x) x[2])
  }

  nresult <- max(length(threshold),
                 length(freq.bins))

  if (length(freq.lower.bound) < nresult) {
    freq.lower.bound <- rep_len(freq.lower.bound, nresult)
  }

  if (length(freq.upper.bound) < nresult) {
    freq.upper.bound <- rep_len(freq.upper.bound, nresult)
  }

  if (length(threshold) < nresult) {
    threshold <- rep_len(threshold, nresult)
  }

  contributions <- matrix(NA, length(groups), nresult)
  for (i in seq_len(nresult)) {
    contributions[, i] <- sapply(pgs,
                                 function(pg) sum(pg$pgram[pg$freqs <= freq.upper.bound[i] & pg$freqs >= freq.lower.bound[i]]))
  }

  if (all(threshold <= 0)) {
    gi <- max.col(contributions)
    result <- lapply(seq_len(nresult), function(i) groups[gi == i])
  } else {
    result <- lapply(seq_len(nresult), function(i) groups[contributions[, i] >= threshold[i]])
  }

  names(result) <- if (!is.null(names(freq.bins))) .group.names(freq.bins) else .group.names(threshold)

  if (drop) {
    result <- result[sapply(result, length) > 0]
  }

  result
}

grouping.auto.toeplitz.ssa <- grouping.auto.1d.ssa
