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


prepanel.reconstruction.2d.ssa <- function(z, subscripts, recon, ...) {
  N <- dim(recon[[z[subscripts]]])
  y <- c(seq_len(N[1]), rep(1, N[2]))
  x <- c(seq_len(N[2]), rep(1, N[1]))
  prepanel.default.levelplot(x, y, subscripts = seq_along(x))
}

panel.reconstruction.2d.ssa <- function(x, y, z, recon, subscripts, at, ...,
                                        ref = FALSE,
                                        symmetric = FALSE,
                                        .cuts = 20,
                                        .useRaster = FALSE,
                                        region, contour,
                                        fill.uncovered = "void",
                                        fill.color = NULL) {
  if (is.list(fill.uncovered))
    fill.uncovered <- fill.uncovered[[(subscripts - 1) %% length(fill.uncovered) + 1]]
  if (is.list(fill.color))
    fill.color <- fill.color[[(subscripts - 1) %% length(fill.color) + 1]]
  if (!is.null(fill.color))
    panel.fill(col = fill.color)

  if (is.list(at))
    at <- at[[(subscripts - 1) %% length(at) + 1]]

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

    if (diff(z.range) < .Machine$double.eps^.5)
      z.range <- z.range + c(-1, 1) * .Machine$double.eps^.25

    at <- seq(z.range[1], z.range[2], length.out = .cuts + 2)
  }

  # Cutoff outstanding values
  data$z[data$z < min(at)] <- min(at)
  data$z[data$z > max(at)] <- max(at)

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
                                       col = grey(c(0, 1)),
                                       zlim,
                                       at) {
  dots <- list(...)
  type <- match.arg(type)

  if (!missing(zlim))
    at <- "same"
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
                    cuts = 20,
                    colorkey = !(identical(at, "free") || (is.list(at) && length(at) > 1)),
                    symmetric = FALSE,
                    ref = FALSE,
                    useRaster = TRUE,
                    fill.uncovered = "void")
  dots <- modifyList(dots,
                     list(par.settings = list(regions = list(col = colorRampPalette(col)))))

  # Disable colorkey if subplots are drawing in different scales
  if (identical(at, "free") || (is.list(at) && length(at) > 1))
    dots$colorkey <- FALSE

  if (identical(at, "same")) {
    all.values <- if (missing(zlim)) unlist(x) else zlim
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

  do.call("levelplot",
          c(list(x = z ~ row * column | factor(z, labels = labels),
                 data = d, recon = x,
                 at = at,
                 panel = panel.reconstruction.2d.ssa,
                 prepanel = prepanel.reconstruction.2d.ssa),
            dots))
}

prepanel.eigenvectors.2d.ssa <- function(x, y, subscripts, ssaobj, ..., what = "eigen") {
  if (identical(what, "eigen")) {
    Ls <- ssaobj$window
  } else if (identical(what, "factor")) {
    Ls <- ifelse(ssaobj$circular, ssaobj$length, ssaobj$length - ssaobj$window + 1)
  }

  y <- c(seq_len(Ls[1]), rep(1, Ls[2]))
  x <- c(seq_len(Ls[2]), rep(1, Ls[1]))

  prepanel.default.levelplot(x, y, subscripts = seq_along(x))
}

panel.eigenvectors.2d.ssa <- function(x, y, z, ssaobj, subscripts, at, ...,
                                      what = "eigen",
                                      ref = FALSE,
                                      symmetric = FALSE,
                                      .cuts = 20,
                                      .useRaster = FALSE,
                                      region, contour,
                                      fill.color = NULL) {
  if (is.list(fill.color))
    fill.color <- fill.color[[(subscripts - 1) %% length(fill.color) + 1]]
  if (!is.null(fill.color))
    panel.fill(col = fill.color)

  if (is.list(at))
    at <- at[[(subscripts - 1) %% length(at) + 1]]

  panel <- if (.useRaster) panel.levelplot.raster else panel.levelplot

  idx <- z[subscripts]
  if (identical(what, "eigen")) {
    Ls <- ssaobj$window
    vmask <- .get(ssaobj, "wmask")
    vectors <- ssaobj$U[, idx]
  } else if (identical(what, "factor")) {
    Ls <- ifelse(ssaobj$circular, ssaobj$length, ssaobj$length - ssaobj$window + 1)
    vmask <- .get(ssaobj, "fmask")
    vectors <- if (nv(ssaobj) >= max(idx)) ssaobj$V[, idx] else calc.v(ssaobj, idx)
  }

  if (is.null(vmask))
    vmask <- matrix(TRUE, Ls[1], Ls[2])

  data <- expand.grid(y = rev(seq_len(Ls[1])), x = seq_len(Ls[2]))[as.vector(vmask), ]
  data$z <- vectors

  if (identical(at, "free")) {
    z.range <- range(if (symmetric) c(data$z, -data$z) else data$z)

    if (diff(z.range) < .Machine$double.eps^.5)
      z.range <- z.range + c(-1, 1) * .Machine$double.eps^.25

    at <- seq(z.range[1], z.range[2], length.out = .cuts + 2)
  }

  # Cutoff outstanding values
  data$z[data$z < min(at)] <- min(at)
  data$z[data$z > max(at)] <- max(at)

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

.plot.ssa.vectors.2d.ssa <- function(x, ...,
                                     what = c("eigen", "factor"),
                                     col = grey(c(0, 1)),
                                     zlim,
                                     at,
                                     plot.contrib = FALSE, idx) {
  dots <- list(...)
  what <- match.arg(what)

  if (max(idx) > nsigma(x))
    stop("Too few eigentriples computed for this decomposition")

  if (!missing(zlim))
    at <- "same"
  if (missing(at))
    at <- "free"
  if (is.character(at))
    at <- match.arg(at, c("free", "same"))

  # FIXME: check for proper lengths
  d <- data.frame(row = idx, column = idx, z = idx)

  # Provide convenient defaults
  dots <- .defaults(dots,
                    xlab = "",
                    ylab = "",
                    main = if (identical(what, "eigen")) "Eigenvectors" else "Factor vectors",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "same"),
                    aspect = "iso",
                    cuts = 20,
                    symmetric = FALSE,
                    ref = FALSE,
                    useRaster = TRUE)
  dots <- modifyList(dots,
                     list(par.settings = list(regions = list(col = colorRampPalette(col)))))

  # Disable colorkey if subplots are drawed in different scales
  if (identical(at, "free") || (is.list(at) && length(at) > 1))
    dots$colorkey <- FALSE

  if (identical(at, "same")) {
    if (missing(zlim)) {
      if (identical(what, "eigen")) {
        all.values <- range(.U(x)[, idx])
      } else if (identical(what, "factor")) {
        all.values <- range(calc.v(x, idx))
      }
    } else {
      all.values <- zlim
    }
    at <- pretty(if (dots$symmetric) c(all.values, -all.values) else all.values, n = dots$cuts)
  }

  # Rename args for transfer to panel function
  names(dots)[names(dots) == "cuts"] <- ".cuts"
  names(dots)[names(dots) == "useRaster"] <- ".useRaster"

  do.call("levelplot",
          c(list(x = z ~ row * column | factor(z,
                                               labels = if (!plot.contrib) z else paste(z, " (", .contribution(x, z, ...), "%)", sep = "")),
                 data = d, ssaobj = x,
                 what = what,
                 at = at,
                 panel = panel.eigenvectors.2d.ssa,
                 prepanel = prepanel.eigenvectors.2d.ssa),
            dots));
}
