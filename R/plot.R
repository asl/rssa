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
  if (is.complex(V)) V <- c(Re(V), Im(V))
  U <- if (identical(x, y)) 1:length(V) else ssaobj$U[, x]
  if (is.complex(U)) U <- c(Re(U), Im(U))

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
  if (is.complex(V)) V <- c(Re(V), Im(V))
  U <- if (identical(x, y)) 1:length(V) else ssaobj$U[, x]
  if (is.complex(U)) U <- c(Re(U), Im(U))

  if (ref) {
    panel.abline(h = 0, ..., reference = TRUE)
    if (!identical(x, y))
      panel.abline(v = 0, ..., reference = TRUE)
  }

  groups <- if (inherits(ssaobj, "cssa")) 2 else 1
  panel.superpose(U, V,
                  panel.groups = panel.xyplot,
                  groups = gl(n = groups, length(V) / groups), ..., subscripts = 1:length(U))
}

prepanel.factorvectors <- function(x, y, ssaobj, symmetric = FALSE) {
  V <- if (y <= nv(ssaobj)) ssaobj$V[, y] else calc.v(ssaobj, idx = y)
  if (is.complex(V)) V <- c(Re(V), Im(V))
  U <- if (identical(x, y)) 1:length(V) else if (x <= nv(ssaobj)) ssaobj$V[, x] else calc.v(ssaobj, idx = x)
  if (is.complex(U)) V <- c(Re(U), Im(U))

  res <- prepanel.default.xyplot(U, V)
  if (symmetric) {
    res$ylim <- range(V, -V)
    if (!identical(x, y))
      res$xlim <- range(U, -U)
  }

  res
}

panel.factorvectors <- function(x, y, ssaobj, ..., ref = FALSE) {
  V <- if (y <= nv(ssaobj)) ssaobj$V[, y] else calc.v(ssaobj, idx = y)
  if (is.complex(V)) V <- c(Re(V), Im(V))
  U <- if (identical(x, y)) 1:length(V) else if (x <= nv(ssaobj)) ssaobj$V[, x] else calc.v(ssaobj, idx = x)
  if (is.complex(U)) V <- c(Re(U), Im(U))

  if (ref) {
    panel.abline(h = 0, ..., reference = TRUE)
    if (!identical(x, y))
      panel.abline(v = 0, ..., reference = TRUE)
  }

  groups <- if (inherits(ssaobj, "cssa")) 2 else if (inherits(ssaobj, "mssa")) length(ssaobj$length) else 1
  panel.superpose(U, V,
                  panel.groups = panel.xyplot,
                  groups = gl(n = groups, length(V) / groups), ..., subscripts = 1:length(U))
}

.defaults <- function(x, ...) {
  dots <- list(...)
  modifyList(dots, x)
}

.plot.ssa.values <- function(x, ..., numvalues, plot.type = "b") {
  dots <- list(...)

  # FIXME: check for proper lengths
  d <- data.frame(A = 1:numvalues, B = x$sigma[1:numvalues])

  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = plot.type,
                    xlab =  "Index",
                    ylab = "log of singular value",
                    main = "Singular Values",
                    grid = TRUE,
                    scales = list(y = list(log = TRUE)),
                    par.settings = list(plot.symbol = list(pch = 20)))

  do.call("xyplot",
          c(list(x = B ~ A , data = d, ssaobj = x), dots))
}

.contribution <- function(x, idx, ...) {
  ## Check for F-orthogonality
  isfcor <- .is.frobenius.orthogonal(x, idx, ...)
  if (!isTRUE(isfcor))
    warning(sprintf("Elementary matrices are not F-orthogonal (max F-cor is %s). Contributions can be irrelevant",
                    format(isfcor, digits = 3)))

  total <- wnorm(x)^2
  round(100*x$sigma[idx]^2 / total, digits = 2)
}

.plot.ssa.vectors <- function(x, ...)
  UseMethod(".plot.ssa.vectors")

.plot.ssa.vectors.ssa <- function(x, ...)
  stop("`.plot.ssa.vectors' is not implemented for this kind of SSA")

.plot.ssa.vectors.1d.ssa <- function(x, ...,
                                     what = c("eigen", "factor"),
                                     plot.contrib, idx, plot.type = "l") {
  dots <- list(...)
  what <- match.arg(what)

  ## FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idx)
  plot.formula <-
    A ~ B | factor(A,
                   labels = if (!plot.contrib) A else paste(A, " (", .contribution(x, idx, ...), "%)", sep = ""))

  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = plot.type,
                    xlab = "",
                    ylab =  "",
                    main = if (identical(what, "eigen")) "Eigenvectors" else "Factor vectors",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "free"),
                    aspect = 1,
                    symmetric = TRUE,
                    ref = TRUE)

  do.call("xyplot",
          c(list(x = plot.formula,
                 data = d, ssaobj = x,
                 panel = if (identical(what, "eigen")) panel.eigenvectors else panel.factorvectors,
                 prepanel = if (identical(what, "eigen")) prepanel.eigenvectors else prepanel.factorvectors),
            dots))
}

.plot.ssa.vectors.toeplitz.ssa <- `.plot.ssa.vectors.1d.ssa`
.plot.ssa.vectors.cssa <- `.plot.ssa.vectors.1d.ssa`

.plot.ssa.paired <- function(x, ...,
                             what = c("eigen", "factor"),
                             plot.contrib, idx, idy, plot.type = "l") {
  dots <- list(...)
  what <- match.arg(what)

  # FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idy)
  plot.formula <- A ~ B | factor(A,
                                 labels =
                                 if (!plot.contrib) paste(A, "vs", B)
                                 else paste(A, " (", .contribution(x, idx, ...), "%) vs ", B, " (", .contribution(x, idy, ...), "%)", sep = ""))

  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = plot.type,
                    xlab = "",
                    ylab = "",
                    main = if (identical(what, "eigen")) "Pairs of eigenvectors" else "Pairs of factor vectors",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "free"),
                    aspect = 1,
                    symmetric = TRUE,
                    ref = TRUE)

  do.call("xyplot",
          c(list(x = plot.formula,
                 data = d, ssaobj = x,
                 panel = if (identical(what, "eigen")) panel.eigenvectors else panel.factorvectors,
                 prepanel = if (identical(what, "eigen")) prepanel.eigenvectors else prepanel.factorvectors),
            dots))
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

.plot.ssa.series <- function(x, ...)
  UseMethod(".plot.ssa.series")

.plot.ssa.series.ssa <- function(x, ..., groups)
  plot(reconstruct(x, groups = groups, drop = FALSE), ...)

.plot.ssa.series.1d.ssa <- function(x, ..., groups, plot.type = "l") {
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

  do.call("xyplot",
          c(list(x = A ~ B | factor(A, labels = paste(groups)),
                 data = d, recon = r,
                 panel = panel.series,
                 prepanel = prepanel.series),
            dots))
}

.plot.ssa.series.toeplitz.ssa <- .plot.ssa.series.1d.ssa

panel.levelplot.wcor <- function(x, y, z, at, ..., grid, .useRaster = FALSE) {
  panel <- if (.useRaster) panel.levelplot.raster else panel.levelplot

  # Cutoff outstanding values
  z[z < min(at)] <- min(at)
  z[z > max(at)] <- max(at)

  panel(x, y, z, at = at, ...)

  if (!is.list(grid))
    grid <- list(h = grid, v = grid)

  panel.abline(v = grid$v - 0.5, ..., reference = TRUE)
  panel.abline(h = grid$h - 0.5, ..., reference = TRUE)
}

plot.wcor.matrix <- function(x,
                             grid = c(),
                             ...,
                             col = grey(c(1, 0)),
                             cuts = 20,
                             zlim = range(abs(x), 0, 1),
                             at) {
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
                    useRaster = TRUE)
  dots <- modifyList(dots,
                     list(par.settings = list(regions = list(col = colorRampPalette(col)))))

  data <- expand.grid(row = seq_len(nrow(x)), column = seq_len(ncol(x)))
  data$x <- as.vector(as.numeric(x))

  # Rename args for transfer to panel function
  names(dots)[names(dots) == "useRaster"] <- ".useRaster"

  if (missing(at))
    at <- pretty(zlim, n = cuts)

  do.call("levelplot",
          c(list(x = abs(x) ~ row * column,
                 data = data,
                 at = at,
                 panel = panel.levelplot.wcor,
                 grid = grid),
            dots))
}

plot.ssa <- function(x,
                     type = c("values", "vectors", "paired", "series", "wcor"),
                     ...,
                     vectors = c("eigen", "factor"),
                     plot.contrib = TRUE,
                     numvalues = nsigma(x),
                     numvectors = min(nsigma(x), 10),
                     idx = 1:numvectors,
                     idy,
                     groups) {
  type <- match.arg(type)
  vectors <- match.arg(vectors)

  if (identical(type, "values")) {
    .plot.ssa.values(x, ..., numvalues = numvalues)
  } else if (identical(type, "vectors")) {
    .plot.ssa.vectors(x, ..., what = vectors, plot.contrib = plot.contrib, idx = idx)
  } else if (identical(type, "paired")) {
    if (missing(idy))
      idy <- idx + 1

    .plot.ssa.paired(x, ..., what = vectors, plot.contrib = plot.contrib, idx = idx, idy = idy)
  } else if (identical(type, "series")) {
    if (missing(groups))
      groups <- as.list(1:min(nsigma(x), nu(x)))

    .plot.ssa.series(x, ..., groups = groups)
  } else if (identical(type, "wcor")) {
    if (missing(groups))
      groups <- as.list(1:min(nsigma(x), nu(x)))

    plot(wcor(x, groups = groups), ...)
  } else {
    stop("Unsupported type of SSA plot!")
  }
}

plot.1d.ssa.reconstruction <- function(x, ...,
                                       type = c("raw", "cumsum"),
                                       plot.method = c("native", "matplot", "xyplot"),
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

  cnames <- names(x)
  mnames <- ifelse(cnames == "", 1:ncol(m), cnames)
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
  if (identical(plot.method, "xyplot"))
    do.call(xyplot, c(list(m), dots))
  else if (identical(plot.method, "matplot") || !is.object(m))
    do.call(matplot, c(list(x = m), dots))
  else if (identical(plot.method, "native"))
    do.call(plot, c(list(m), dots))
  else
    stop("Unknown plot method")
}

plot.toeplitz.ssa.reconstruction <- `plot.1d.ssa.reconstruction`

prepanel.roots <- function(x, y, ...) {
  lim <- range(x, y, -x, -y, 1, -1)
  list(xlim = lim, ylim = lim)
}

panel.roots <- function(...) {
  tt <- seq(0, 2 * pi, length.out = 100)
  par <- trellis.par.get("reference.line")
  panel.lines(sin(tt), cos(tt), col = par$col, lty = par$lty, lwd = par$lwd) # TODO MB use `TeachingDemos' package here?
  panel.xyplot(...)
}

plot.fdimpars.1d <- function(x, ...) {
  dots <- list(...)

  # Provide convenient defaults
  dots <- .defaults(dots,
                    main = "Roots",
                    xlab = "Real part",
                    ylab = "Imaginary part",
                    aspect = "iso",
                    pch = 19)

  do.call("xyplot",
          c(list(Im(roots) ~ Re(roots),
                 data = x,
                 panel = panel.roots,
                 prepanel = prepanel.roots),
            dots))
}

plot.fdimpars.2d <- function(x, ...) {
  dots <- list(...)

  # Provide convenient defaults
  dots <- .defaults(dots,
                    main = "Roots",
                    xlab = "Real part",
                    ylab = "Imaginary part",
                    aspect = 1,
                    pch = 19)

  data <- list()
  data$root <- c(x[[1]]$roots, x[[2]]$roots)
  data$ind <- rep(c("lambda", "mu"), each = length(x[[1]]$roots))

  do.call("xyplot",
          c(list(Im(root) ~ Re(root) | ind,
                 data = data,
                 panel = panel.roots,
                 prepanel = prepanel.roots),
            dots))
}

plot.lrr <- function(x, ..., raw = FALSE) {
  if (!raw) {
    r <- roots(x)
    xyplot(Im(r) ~ Re(r), ...)
  }

  dots <- list(...)

  # Provide convenient defaults
  dots <- .defaults(dots,
                    main = "Roots of Linear Recurrence Relation")

  do.call("plot", c(list(roots2pars(roots(x))), dots))
}
