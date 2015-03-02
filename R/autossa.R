#   R package for Singular Spectrum Analysis
#   Copyright (c) 2008-2015 Anton Korobeynikov <anton@korobeynikov.info>
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


grouping.auto <- function(x, ...,
                          grouping.method = c("pgram", "wcor")) {
  switch(match.arg(grouping.method),
         pgram = grouping.auto.pgram(x, ...),
         wcor  = grouping.auto.wcor(x, ...))
}

grouping.auto.wcor <- function(x, ...)
  UseMethod("grouping.auto.wcor")

grouping.auto.wcor.ssa <- function(x,
                                   groups, nclust = length(groups) / 2,
                                   ...) {
  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))

  groups <- sort(unique(unlist(groups)))
  
  w <- wcor(x, groups = groups, ...)
  h <- hclust(as.dist((1 - w) / 2), ...)
  res <- split(groups, cutree(h, k = nclust))

  attr(res, "hclust") <- h
  attr(res, "wcor") <- w
  class(res) <- "grouping.auto.wcor"

  res
}

plot.grouping.auto.wcor <- function(x, ...)
  plot(attr(x, "hclust"), ...)

grouping.auto.pgram <- function(x, ...)
  UseMethod("grouping.auto.pgram")

grouping.auto.pgram.ssa <- function(x, ...)
  stop("grouping.auto.pgram is not implemented for this kind of SSA yet")

pgram <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  stopifnot(all(is.finite(x)))

  X <- mvfft(x)
  n <- nrow(x)

  N <- n %/% 2 + 1
  spec <- abs(X[seq_len(N),, drop = FALSE])^2

  if (n %% 2 == 0) {
    if (N > 2) spec[2:(N-1), ] <- 2 * spec[2:(N-1), ]
  } else {
    if (N >= 2) spec[2:N, ] <- 2 * spec[2:N, ]
  }

  freq <- seq(0, 1, length.out = n + 1)[seq_len(N)]

  cumspecfuns <- lapply(seq_len(ncol(x)),
                        function(j)
                          approxfun(c(0, freq[-N] + 1/(2*n), 0.5),
                                    c(0, cumsum(spec[, j])),
                                    rule = 2))

  list(spec = spec, freq = freq, cumspecfuns = cumspecfuns)
}

grouping.auto.pgram.1d.ssa <- function(x, groups,
                                       base = c("series", "eigen", "factor"),
                                       freq.bins = 2,
                                       threshold = 0,
                                       method = c("constant", "linear"),
                                       ...,
                                       drop = TRUE) {
  if (is.shaped(x)) {
    stop("grouping.auto is not implemented for shaped SSA yet")
  }

  base <- match.arg(base)
  method <- match.arg(method)

  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  groups <- sort(unique(unlist(groups)))

  if (identical(base, "eigen")) {
      Fs <- .U(x)[, groups, drop = FALSE]
  } else if (identical(base, "factor")) {
      Fs <- calc.v(x, groups, ...)
  } else if (identical(base, "series")) {
    N <- x$length
    Fs <- matrix(unlist(reconstruct(x, groups = as.list(groups), ...)), nrow = N)
  }

  pgs <- pgram(Fs)

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

  norms <- colSums(pgs$spec)
  contributions <- matrix(NA, length(groups), nresult)
  for (i in seq_len(nresult)) {
    contributions[, i] <- if (identical(method, "constant"))
      colSums(pgs$spec[pgs$freq < freq.upper.bound[i] & pgs$freq >= freq.lower.bound[i],, drop = FALSE]) / norms
    else if (identical(method, "linear"))
      sapply(pgs$cumspecfuns, function(f) diff(f(c(freq.lower.bound[i], freq.upper.bound[i])))) / norms
  }

  type <- if (all(threshold <= 0)) "splitting" else "independent"

  if (identical(type, "splitting")) {
    gi <- max.col(contributions, ties.method = "first")
    result <- lapply(seq_len(nresult), function(i) groups[gi == i])
  } else if (identical(type, "independent")) {
    result <- lapply(seq_len(nresult), function(i) groups[contributions[, i] >= threshold[i]])
  } else {
    stop("Unknown type for pgrouping.auto.pgram")
  }

  names(result) <- if (!is.null(names(freq.bins))) .group.names(freq.bins) else .group.names(threshold)
  colnames(contributions) <- names(result)
  rownames(contributions) <- as.character(groups)

  if (drop) {
    result <- result[sapply(result, length) > 0]
  }

  attr(result, "contributions") <- contributions
  attr(result, "type") <- type
  attr(result, "threshold") <- threshold

  class(result) <- "grouping.auto.pgram"

  result
}

grouping.auto.pgram.toeplitz.ssa <- grouping.auto.pgram.1d.ssa

plot.grouping.auto.pgram <- function(x, superpose, order, ...) {
  type <- attr(x, "type")

  if (missing(order))
    order <- switch(type,
                    splitting = FALSE,
                    independent = TRUE)

  if (missing(superpose))
    superpose <- switch(type,
                        splitting = FALSE,
                        independent = TRUE)

  if (order && identical(type, "splitting"))
    warning("Ordering has no sense for splitting grouping")

  contributions <- as.data.frame(attr(x, "contributions"))
  components <- rownames(contributions)

  indices <- lapply(contributions, seq_along)
  labels <- lapply(indices, function(i) components[i])

  if (order) {
    for (i in seq_along(contributions)) {
      idx <- base::order(contributions[[i]], decreasing = TRUE)
      contributions[[i]] <- contributions[[i]][idx]
      labels[[i]] <- labels[[i]][idx]
    }
  }

  gind <- lapply(seq_len(ncol(contributions)), function(i) rep(i, nrow(contributions)))
  groups <- lapply(names(contributions), function(name) rep(name, nrow(contributions)))

  data <- data.frame(contribution = unlist(contributions),
                     index = unlist(indices),
                     gind = unlist(gind),
                     group = unlist(groups))

  if (!order) {
    xscale <- list(limits = labels[[1]])
  } else {
    if (superpose) {
      xscale <- list(draw = FALSE)
    } else {
      xscale <- list(limits = labels, relation = "free")
    }
  }

  dots <- list(...)
  dots <- .defaults(dots,
                    auto.key = superpose,
                    xlab = "Component",
                    ylab = "Relative contribution",
                    scales = list(x = xscale),
                    as.table = !superpose)

  if (superpose) {
    do.call(xyplot, c(list(contribution ~ index, groups = data$group, data = data), dots))
  } else {
    do.call(xyplot, c(list(contribution ~ index | group, data = data), dots))
  }
}
