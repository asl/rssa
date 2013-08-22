#   R package for Singular Spectrum Analysis
#   Copyright (c) 2013 Alex Shlemov <shlemovalex@gmail.com>
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


panel.reconstruction.2d.ssa <- function(x, y, z, recon, subscripts, at, ...,
                                        ref = FALSE,
                                        symmetric = FALSE,
                                        .cuts = 20,
                                        .useRaster = FALSE,
                                        region, contour,
                                        fill.uncovered = "void") {
  panel <- if (.useRaster) panel.levelplot.raster else panel.levelplot
  N <- dim(recon[[subscripts]])
  data <- expand.grid(y = rev(seq_len(N[1])), x = seq_len(N[2]))
  data$z <- as.vector(recon[[z[subscripts]]])

  if (is.character(fill.uncovered)) {
    fill.uncovered <- match.arg(fill.uncovered, choices = c("mean", "original", "void"))
    fill.uncovered <- switch(fill.uncovered,
                             mean = mean(data$z, na.rm = TRUE),
                             original = attr(recon, "series"),
                             void = NA)
  }

  if (!is.matrix(fill.uncovered)) {
    fill.uncovered <- matrix(fill.uncovered, N[1], N[2])
  }

  stopifnot(all(dim(fill.uncovered) == N))
  data$z[is.na(data$z)] <- fill.uncovered[is.na(data$z)]

  if (identical(at, "free")) {
    z.range <- range(if (symmetric) c(data$z, -data$z) else data$z, na.rm = TRUE)
    at <- seq(z.range[1], z.range[2], length.out = .cuts + 2)
  }

  panel(data$x, data$y, data$z, subscripts = seq_len(nrow(data)),
        at = at, contour = FALSE, region = TRUE, ...)

  if (ref) {
    # Draw zerolevel isoline
    par <- trellis.par.get("reference.line")
    panel.contourplot(data$x, data$y, data$z, subscripts = seq_len(nrow(data)),
                      at = 0, contour = TRUE, region = FALSE,
                      col = par$col, lty = par$lty, lwd = par$lwd)
  }
}

plot.2d.ssa.reconstruction <- function(x, ...,
                                       type = c("raw", "cumsum"),
                                       base.series = NULL,
                                       add.original = TRUE,
                                       add.residuals = TRUE,
                                       add.ranges,
                                       at) {
  dots <- list(...)
  type <- match.arg(type)

  if (missing(at))
    at <- if (identical(type, "cumsum")) "same" else "free"
  if (is.character(at))
    at <- match.arg(at, c("same", "free"))

  if (missing(add.ranges))
    add.ranges <- identical(at, "free")

  # Save `x' attributes
  xattr <- attributes(x)

  if (identical(type, "cumsum") && (length(x) > 1)) {
    for (i in 2:length(x))
      x[[i]] <- x[[i]] + x[[i - 1]]

    names(x)[2:length(x)] <- paste(names(x)[1], names(x)[2:length(x)], sep = ":")
  }

  original <- attr(x, "series")
  residuals <- attr(x, "residuals")

  # Handle base series, if any
  if (!is.null(base.series)) {
    stopifnot(inherits(base.series, "ssa.reconstruction"))
    original <- attr(base.series, "series")
    names(base.series) <- paste("base", names(base.series), sep = " ")
    x <- c(base.series, x)
  }

  if (add.original)
    x <- c(list(Original = original), x)
  if (add.residuals)
    x <- c(x, list(Residuals = residuals))

  # Restore `x' attributes
  attributes(x)[c("series", "residuals")] <- xattr[c("series", "residuals")]

  idx <- seq_along(x)
  d <- data.frame(row = idx, column = idx, z = idx)

  # Provide convenient defaults
  dots <- .defaults(dots,
                    xlab = "",
                    ylab = "",
                    main = "Reconstructions",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "same"),
                    aspect = "iso",
                    par.settings = list(regions = list(col = colorRampPalette(grey(c(0, 1))))),
                    cuts = 20,
                    colorkey = !identical(at, "free"),
                    symmetric = FALSE,
                    ref = FALSE,
                    useRaster = TRUE,
                    fill.uncovered = "void")

  # Disable colorkey if subplots are drawing in different scales
  if (identical(at, "free"))
    dots$colorkey <- FALSE

  if (identical(at, "same")) {
    all.values <- unlist(x)
    at <- pretty(if (dots$symmetric) c(all.values, -all.values) else all.values, n = dots$cuts)
  }

  # Rename args for transfer to panel function
  names(dots)[names(dots) == "cuts"] <- ".cuts"
  names(dots)[names(dots) == "useRaster"] <- ".useRaster"

  labels <- names(x)
  if (add.ranges) {
    ranges <- sapply(x, range, na.rm = TRUE)
    ranges <- signif(ranges, 2)
    labels <- paste(labels, " [", ranges[1, ],", ", ranges[2, ], "]", sep = "")
  }

  res <- do.call("levelplot",
                 c(list(x = z ~ row * column | factor(z, labels = labels),
                        data = d, recon = x,
                        at = at,
                        useRaster = dots$.useRaster,
                        panel = panel.reconstruction.2d.ssa,
                        prepanel = prepanel.reconstruction.2d.ssa),
                   dots))
  print(res)
}

prepanel.eigenvectors.2d.ssa <- function(x, y, subscripts, ssaobj, ...) {
  L <- ssaobj$window
  y <- c(seq_len(L[1]), rep(1, L[2]))
  x <- c(seq_len(L[2]), rep(1, L[1]))
  prepanel.default.levelplot(x, y, subscripts = seq_along(x))
}

panel.eigenvectors.2d.ssa <- function(x, y, z, ssaobj, subscripts, at, ...,
                                      ref = FALSE,
                                      symmetric = FALSE,
                                      .cuts = 20,
                                      .useRaster = FALSE,
                                      region, contour) {
  panel <- if (.useRaster) panel.levelplot.raster else panel.levelplot
  L <- ssaobj$window
  wmask <- .get(ssaobj, "wmask")
  if (is.null(wmask))
    wmask <- matrix(TRUE, L[1], L[2])

  data <- expand.grid(y = rev(seq_len(L[1])), x = seq_len(L[2]))[as.vector(wmask), ]
  data$z <- ssaobj$U[, z[subscripts]]

  if (identical(at, "free")) {
    z.range <- range(if (symmetric) c(data$z, -data$z) else data$z)
    at <- seq(z.range[1], z.range[2], length.out = .cuts + 2)
  }
  panel.levelplot(data$x, data$y, data$z, subscripts = seq_len(nrow(data)), at = at, ...)

  panel(data$x, data$y, data$z, subscripts = seq_len(nrow(data)),
        at = at, contour = FALSE, region = TRUE, ...)

  if (ref) {
    # Draw zerolevel isoline
    par <- trellis.par.get("reference.line")
    panel.contourplot(data$x, data$y, data$z, subscripts = seq_len(nrow(data)),
                      at = 0, contour = TRUE, region = FALSE,
                      col = par$col, lty = par$lty, lwd = par$lwd)
  }
}

.plot.ssa.vectors.2d.ssa <- function(x, ..., plot.contrib = FALSE, idx, at) {
  dots <- list(...)

  if (missing(at))
    at <- "free"
  if (is.character(at))
    at <- match.arg(at, c("free", "same"))

  # FIXME: check for proper lengths
  d <- data.frame(row = idx, column = idx, z = idx)

  if (plot.contrib) {
    total <- wnorm(x)^2
    lambda <- round(100*x$lambda[idx]^2 / total, digits = 2);
  }

  # Provide convenient defaults
  dots <- .defaults(dots,
                    xlab = "",
                    ylab = "",
                    main = "Eigenvectors",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "same"),
                    aspect = "iso",
                    par.settings = list(regions = list(col = colorRampPalette(grey(c(0, 1))))),
                    cuts = 20,
                    symmetric = FALSE,
                    ref = FALSE,
                    useRaster = TRUE)

  # Disable colorkey if subplots are drawed in different scales
  if (identical(at, "free"))
    dots$colorkey <- FALSE

  if (identical(at, "same")) {
    all.values <- x$U[, idx]
    at <- pretty(if (dots$symmetric) c(all.values, -all.values) else all.values, n = dots$cuts)
  }

  # Rename args for transfer to panel function
  names(dots)[names(dots) == "cuts"] <- ".cuts"
  names(dots)[names(dots) == "useRaster"] <- ".useRaster"

  do.call("levelplot",
          c(list(x = z ~ row * column | factor(z,
                                               labels = if (!plot.contrib) z else paste(z, " (", lambda, "%)", sep = "")),
                 data = d, ssaobj = x,
                 at = at,
                 useRaster = dots$.useRaster,
                 panel = panel.eigenvectors.2d.ssa,
                 prepanel = prepanel.eigenvectors.2d.ssa),
            dots));
}
