#   R package for Singular Spectrum Analysis
#   Copyright (c) 2008 Anton Korobeynikov <asl@math.spbu.ru>
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


prepanel.eigenvectors <- function(x, y, ssaobj, symmetric = FALSE) {
  V <- ssaobj$U[, y]
  U <- if (identical(x, y)) 1:length(V)
       else ssaobj$U[, x]

  res <- prepanel.default.xyplot(U, V)
  if (symmetric) {
    res$ylim <- range(V, -V)
    if (!identical(x, y))
      res$xlim <- range(U, -U)
  }

  res
}

panel.eigenvectors <- function(x, y, ssaobj, ..., ref = FALSE) {
  V <- ssaobj$U[, y]
  U <- if (identical(x, y)) 1:length(V)
       else ssaobj$U[, x]

  if (ref) {
    panel.abline(h = 0, ..., reference = TRUE)
    if (!identical(x, y))
      panel.abline(v = 0, ..., reference = TRUE)
  }

  panel.xyplot(U, V, ...)
}

.defaults <- function(x, ...) {
  dots <- list(...)
  modifyList(dots, x)
}

.plot.ssa.values <- function(x, ..., numvalues, plot.type = "b") {
  dots <- list(...)

  # FIXME: check for proper lengths
  d <- data.frame(A = 1:numvalues, B = x$lambda[1:numvalues])

  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = plot.type,
                    xlab =  "Index",
                    ylab = "log of eigenvalue",
                    main = "Eigenvalues",
                    grid = TRUE,
                    scales = list(y = list(log = TRUE)),
                    par.settings = list(plot.symbol = list(pch = 20)))

  res <- do.call("xyplot",
                 c(list(x = B ~ A , data = d, ssaobj = x), dots))
  print(res)
}

.plot.ssa.vectors <- function(x, ...)
  UseMethod(".plot.ssa.vectors")

.plot.ssa.vectors.ssa <- function(x, ...)
  stop("`.plot.ssa.vectors' is not implemented for this kind of SSA")

.plot.ssa.vectors.1d.ssa <- function(x, ..., plot.contrib, idx, plot.type = "l") {
  dots <- list(...)

  # FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idx)

  if (plot.contrib) {
    total <- wnorm(x)^2
    lambda <- round(100*x$lambda[idx]^2 / total, digits = 2)
  }

  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = plot.type,
                    xlab = "",
                    ylab =  "",
                    main = "Eigenvectors",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "free"),
                    aspect = 1,
                    symmetric = TRUE,
                    ref = TRUE)

  res <- do.call("xyplot",
                 c(list(x = A ~ B | factor(A,
                                           labels = if (!plot.contrib) A else paste(A, " (", lambda, "%)", sep = "")),
                        data = d, ssaobj = x,
                        panel = panel.eigenvectors,
                        prepanel = prepanel.eigenvectors),
                   dots))
  print(res)
}

.plot.ssa.vectors.toeplitz.ssa <- `.plot.ssa.vectors.1d.ssa`

.plot.ssa.paired <- function(x, ..., plot.contrib, idx, idy, plot.type = "l") {
  dots <- list(...)

  # FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idy)

  if (plot.contrib) {
    total <- wnorm(x)^2
    lambdax <- round(100*x$lambda[idx]^2 / total, digits = 2)
    lambday <- round(100*x$lambda[idy]^2 / total, digits = 2)
  }

  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = plot.type,
                    xlab = "",
                    ylab = "",
                    main = "Pairs of eigenvectors",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "free"),
                    aspect = 1,
                    symmetric = TRUE,
                    ref = TRUE)

  res <- do.call("xyplot",
                 c(list(x = A ~ B | factor(A,
                                           labels = if (!plot.contrib) paste(A, "vs", B)
                                                    else paste(A, " (", lambdax, "%) vs ", B, " (", lambday, "%)", sep = "")),
                        data = d, ssaobj = x,
                        panel = panel.eigenvectors,
                        prepanel = prepanel.eigenvectors),
                   dots))
  print(res)
}

prepanel.series <- function(x, y, recon, ..., symmetric = FALSE) {
  Y <- recon[[paste("F", y, sep = "")]]
  X <- if (identical(x, y)) time(Y)
       else  recon[[paste("F", x, sep = "")]]

  res <- prepanel.default.xyplot(X, Y, ...)
  if (symmetric) {
    res$ylim <- range(Y, -Y)
    if (!identical(x, y))
      res$xlim <- range(X, -X)
  }

  res
}

panel.series <- function(x, y, recon, ..., ref = FALSE) {
  Y <- recon[[paste("F", y, sep = "")]]
  X <- if (identical(x, y)) time(Y)
       else  recon[[paste("F", x, sep = "")]]

  if (ref) {
    panel.abline(h = 0, ..., reference = TRUE)
    if (!identical(x, y))
      panel.abline(v = 0, ..., reference = TRUE)
  }

  panel.xyplot(X, Y, ...)
}

.plot.ssa.series <- function(x, ..., groups, plot.type = "l") {
  dots <- list(...)

  # FIXME: check for proper lengths
  idx <- seq_along(groups)
  d <- data.frame(A = idx, B = idx)

  r <- reconstruct(x, groups = groups, drop = FALSE)

  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = plot.type,
                    xlab = "",
                    ylab = "",
                    main = "Reconstructed series",
                    as.table = TRUE,
                    scales = list(relation = "free"))

  res <- do.call("xyplot",
                 c(list(x = A ~ B | factor(A, labels = paste(groups)),
                        data = d, recon = r,
                        panel = panel.series,
                        prepanel = prepanel.series),
                   dots))
  print(res)
}

panel.levelplot.wcor <- function(x, y, z, ..., grid, .useRaster = FALSE) {
  panel <- if (.useRaster) panel.levelplot.raster else panel.levelplot
  panel(x, y, z, ...)

  if (!is.list(grid))
    grid <- list(h = grid, v = grid)

  panel.abline(v = grid$v - 0.5, ..., reference = TRUE)
  panel.abline(h = grid$h - 0.5, ..., reference = TRUE)
}

plot.wcor.matrix <- function(x,
                             grid = c(),
                             ...,
                             cuts = 20,
                             zlim = range(abs(x), 0, 1)) {
  # Provide convenient defaults
  dots <- list(...)
  dots <- .defaults(dots,
                    xlab = "",
                    ylab = "",
                    colorkey = FALSE,
                    main = "W-correlation matrix",
                    aspect = "iso",
                    xlim = rownames(x),
                    ylim = colnames(x),
                    par.settings = list(regions = list(col = colorRampPalette(grey(c(1, 0))))),
                    useRaster = TRUE)

  data <- expand.grid(row = seq_len(nrow(x)), column = seq_len(ncol(x)))
  data$x <- as.vector(as.numeric(x))

  # Rename args for transfer to panel function
  names(dots)[names(dots) == "useRaster"] <- ".useRaster"

  res <- do.call("levelplot",
                 c(list(x = abs(x) ~ row * column,
                        data = data,
                        at = seq(zlim[1], zlim[2], length.out = cuts + 2),
                        panel = panel.levelplot.wcor,
                        grid = grid,
                        useRaster = dots$.useRaster),
                 dots))
  print(res)
}

plot.ssa <- function(x,
                     type = c("values", "vectors", "paired", "series", "wcor"),
                     ...,
                     plot.contrib = TRUE,
                     numvalues = nlambda(x),
                     numvectors = min(nlambda(x), 10),
                     idx = 1:numvectors,
                     idy,
                     groups) {
  type <- match.arg(type)

  if (identical(type, "values")) {
    .plot.ssa.values(x, ..., numvalues = numvalues)
  } else if (identical(type, "vectors")) {
    .plot.ssa.vectors(x, ..., plot.contrib = plot.contrib, idx = idx)
  } else if (identical(type, "paired")) {
    if (missing(idy))
      idy <- idx + 1

    .plot.ssa.paired(x, ..., plot.contrib = plot.contrib, idx = idx, idy = idy)
  } else if (identical(type, "series")) {
    if (missing(groups))
      groups <- as.list(1:min(nlambda(x), nu(x)))

    .plot.ssa.series(x, ..., groups = groups)
  } else if (identical(type, "wcor")) {
    if (missing(groups))
      groups <- as.list(1:min(nlambda(x), nu(x)))

    plot(wcor(x, groups = groups), ...)
  } else {
    stop("Unsupported type of SSA plot!")
  }
}

plot.1d.ssa.reconstruction <- function(x, ...,
                                       type = c("raw", "cumsum"),
                                       plot.method = c("native", "matplot"),
                                       base.series = NULL,
                                       add.original = TRUE,
                                       add.residuals = TRUE) {
  type <- match.arg(type)
  plot.method <- match.arg(plot.method)
  original <- attr(x, "series")
  res <- attr(x, "residuals")

  # Handle base series, if any
  if (!is.null(base.series)) {
    stopifnot(inherits(base.series, "ssa.reconstruction"))
    m0 <- matrix(unlist(base.series), ncol = length(base.series))
    original <- attr(base.series, "series")
  }

  # Nifty defaults
  dots <- list(...)
  dots <- .defaults(dots,
                    main = "Reconstructed Series",
                    type = "l",
                    ylab = "")

  # Prepare the matrix with all the data
  m <- matrix(unlist(x), ncol = length(x))
  if (!is.null(base.series)) m <- cbind(m0, m)

  # Transform the matrix, if necessary
  if (identical(type, "cumsum"))
    m <- t(apply(m, 1, cumsum))

  # Merge the attributes in
  attributes(m) <- append(attributes(m), attributes(x[[1]]))

  mnames <- paste("Reconstructed", 1:ncol(m))
  if (add.original) {
    m <- cbind(original, m)
    mnames <- c("Original", mnames)
  }
  if (add.residuals) {
    m <- cbind(m, res)
    mnames <- c(mnames, "Residuals")
  }
  colnames(m) <- mnames

  # Plot'em'all!
  if (identical(plot.method, "matplot") || !is.object(m))
    do.call(matplot, c(list(x = m), dots))
  else if (identical(plot.method, "native"))
    do.call(plot, c(list(m), dots))
  else
    stop("Unknown plot method")
}

plot.toeplitz.ssa.reconstruction <- `plot.1d.ssa.reconstruction`

prepanel.reconstruction.2d.ssa <- function(z, subscripts, recon, ...) {
  N <- dim(recon[[z[subscripts]]])
  x <- c(seq_len(N[1]), rep(1, N[2]))
  y <- c(seq_len(N[2]), rep(1, N[1]))
  prepanel.default.levelplot(x, y, subscripts = seq_along(x))
}

panel.reconstruction.2d.ssa <- function(x, y, z, recon, subscripts, at, ...,
                                        ref = FALSE,
                                        symmetric = FALSE,
                                        .cuts = 20,
                                        .useRaster = FALSE,
                                        region, contour) {
  panel <- if (.useRaster) panel.levelplot.raster else panel.levelplot
  N <- dim(recon[[subscripts]])
  data <- expand.grid(x = seq_len(N[1]), y = seq_len(N[2]))
  data$z <- as.vector(recon[[z[subscripts]]])

  if (identical(at, "free")) {
    z.range <- range(if (symmetric) c(data$z, -data$z) else data$z)
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


  idx <- seq_along(x)
  d <- data.frame(row = idx, column = idx, z = idx)

  # Provide convenient defaults
  dots <- .defaults(dots,
                    xlab = "i",
                    ylab =  "j",
                    main = "Reconstructions",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "same"),
                    aspect = "iso",
                    par.settings = list(regions = list(col = colorRampPalette(grey(c(0, 1))))),
                    cuts = 20,
                    colorkey = !identical(at, "free"),
                    symmetric = FALSE,
                    ref = FALSE,
                    useRaster = TRUE)

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
  x <- c(seq_len(L[1]), rep(1, L[2]))
  y <- c(seq_len(L[2]), rep(1, L[1]))
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

  data <- expand.grid(x = seq_len(L[1]), y = seq_len(L[2]))
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
                    xlab = "i",
                    ylab =  "j",
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

  res <- do.call("levelplot",
                 c(list(x = z ~ row * column | factor(z,
                                                      labels = if (!plot.contrib) z else paste(z, " (", lambda, "%)", sep = "")),
                        data = d, ssaobj = x,
                        at = at,
                        useRaster = dots$.useRaster,
                        panel = panel.eigenvectors.2d.ssa,
                        prepanel = prepanel.eigenvectors.2d.ssa),
                   dots));
  print(res)
}
