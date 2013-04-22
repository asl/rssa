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


prepanel.eigenvectors <- function(x, y, ssaobj) {
  V <- ssaobj$U[,y];
  U <- if (identical(x, y)) 1:length(V)
       else ssaobj$U[,x]

  prepanel.default.xyplot(U, V);
}

panel.eigenvectors <- function(x, y, ssaobj, ...) {
  V <- ssaobj$U[,y];
  U <- if (identical(x, y)) 1:length(V)
       else ssaobj$U[,x]

  panel.xyplot(U, V, ...);
}

.defaults <- function(v, key, value) {
  if (!(key %in% names(v))) v[[key]] <- value;
  v;
}

.plot.ssa.values <- function(x, ..., numvalues) {
  dots <- list(...);

  # FIXME: check for proper lengths
  d <- data.frame(A = 1:numvalues, B = x$lambda[1:numvalues]);

  # Provide convenient defaults
  dots <- .defaults(dots, "type", c("b", "g"));
  dots <- .defaults(dots, "xlab", "Index");
  dots <- .defaults(dots, "ylab", "log of eigenvalue");
  dots <- .defaults(dots, "main", "Eigenvalues");
  dots <- .defaults(dots, "scales", list(y = list(log = TRUE)));
  dots <- .defaults(dots, "pch", 20);

  res <- do.call("xyplot",
                 c(list(x = B ~ A , data = d, ssaobj = x), dots));
  print(res)
}

.plot.ssa.vectors <- function(x, ..., plot.contrib, idx) {
  dots <- list(...);

  # FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idx);

  if (plot.contrib) {
    total <- sum(x$lambda);
    lambda <- round(100*x$lambda[idx] / total, digits = 2);
  }

  # Provide convenient defaults
  dots <- .defaults(dots, "type", "l");
  dots <- .defaults(dots, "xlab", "");
  dots <- .defaults(dots, "ylab", "");
  dots <- .defaults(dots, "main", "Eigenvectors");
  dots <- .defaults(dots, "as.table", TRUE);
  dots <- .defaults(dots, "scales", list(relation = "free"));

  res <- do.call("xyplot",
                 c(list(x = A ~ B | factor(A,
                                           labels = if (!plot.contrib) A else paste(A, " (", lambda, "%)", sep = "")),
                        data = d, ssaobj = x,
                        panel = panel.eigenvectors,
                        prepanel = prepanel.eigenvectors),
                   dots));
  print(res)
}

.plot.ssa.paired <- function(x, ..., plot.contrib, idx, idy) {
  dots <- list(...);

  # FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idy);

  if (plot.contrib) {
    total <- sum(x$lambda);
    lambdax <- round(100*x$lambda[idx] / total, digits = 2);
    lambday <- round(100*x$lambda[idy] / total, digits = 2);
  }

  # Provide convenient defaults
  dots <- .defaults(dots, "type", "l");
  dots <- .defaults(dots, "xlab", "");
  dots <- .defaults(dots, "ylab", "");
  dots <- .defaults(dots, "main", "Pairs of eigenvectors");
  dots <- .defaults(dots, "as.table", TRUE);
  dots <- .defaults(dots, "scales", list(relation = "free"));
  dots <- .defaults(dots, "aspect", "x");

  res <- do.call("xyplot",
                 c(list(x = A ~ B | factor(A,
                                           labels = if (!plot.contrib) paste(A, "vs", B)
                                                    else paste(A, " (", lambdax, "%) vs ", B, " (", lambday, "%)", sep = "")),
                        data = d, ssaobj = x,
                        panel = panel.eigenvectors,
                        prepanel = prepanel.eigenvectors),
                   dots));
  print(res)
}

prepanel.series <- function(x, y, recon, ...) {
  Y <- recon[[paste("F", y, sep = "")]];
  X <- if (identical(x, y)) 1:length(Y)
       else  recon[[paste("F", x, sep = "")]];

  prepanel.default.xyplot(X, Y, ...);
}

panel.series <- function(x, y, recon, ...) {
  Y <- recon[[paste("F", y, sep = "")]];
  X <- if (identical(x, y)) 1:length(Y)
       else  recon[[paste("F", x, sep = "")]];

  panel.xyplot(X, Y, ...);
}

.plot.ssa.series <- function(x, ..., groups) {
  dots <- list(...);

  # FIXME: check for proper lengths
  idx <- seq_along(groups);
  d <- data.frame(A = idx, B = idx);

  r <- reconstruct(x, groups = groups, drop = FALSE);

  # Provide convenient defaults
  dots <- .defaults(dots, "type", "l");
  dots <- .defaults(dots, "xlab", "");
  dots <- .defaults(dots, "ylab", "");
  dots <- .defaults(dots, "main", "Reconstructed series");
  dots <- .defaults(dots, "as.table", TRUE);
  dots <- .defaults(dots, "scales", list(relation = "free"));

  res <- do.call("xyplot",
                 c(list(x = A ~ B | factor(A, labels = paste(groups)),
                        data = d, recon = r,
                        panel = panel.series,
                        prepanel = prepanel.series),
                   dots));
  print(res);
}

panel.levelplot.wcor <- function(x, y, z, ..., grid) {
  panel.levelplot(x, y, z, ...)

  if (!is.list(grid))
    grid <- list(h = grid, v = grid)

  panel.abline(v = grid$v - 0.5, ..., reference = TRUE)
  panel.abline(h = grid$h - 0.5, ..., reference = TRUE)
}

plot.wcor.matrix <- function(x,
                             grid = c(),
                             ...,
                             cuts = 20,
                             zlim = c(0, 1 + .Machine$double.eps^.5)) {
  # Provide convenient defaults
  dots <- list(...)
  dots <- .defaults(dots, "par.settings", list())
  dots$par.settings <- .defaults(dots$par.settings, "regions", list(col = colorRampPalette(grey(c(1, 0)))))

  dots <- .defaults(dots, "xlab", "")
  dots <- .defaults(dots, "ylab", "")
  dots <- .defaults(dots, "colorkey", FALSE)
  dots <- .defaults(dots, "main", "W-correlation matrix")
  dots <- .defaults(dots, "aspect", "iso")
  dots <- .defaults(dots, "xlim", rownames(x))
  dots <- .defaults(dots, "ylim", colnames(x))

  data <- expand.grid(row = seq_len(nrow(x)), column = seq_len(ncol(x)))
  data$x <- as.vector(as.numeric(x))

  res <- do.call("levelplot", c(list(abs(x) ~ row * column,
                                     data = data,
                                     at = seq(zlim[1], zlim[2], length.out = cuts),
                                     panel = panel.levelplot.wcor,
                                     grid = grid),
                                dots))
  print(res)
}

plot.ssa <- function(x,
                     type = c("values", "vectors", "paired", "series"),
                     ...,
                     plot.contrib = TRUE,
                     numvalues = nlambda(x),
                     numvectors = min(nlambda(x), 10),
                     idx = 1:numvectors,
                     idy,
                     groups) {
  type <- match.arg(type);

  if (identical(type, "values")) {
    .plot.ssa.values(x, ..., numvalues = numvalues);
  } else if (identical(type, "vectors")) {
    .plot.ssa.vectors(x, ..., plot.contrib = plot.contrib, idx = idx);
  } else if (identical(type, "paired")) {
    if (missing(idy))
      idy <- idx + 1;

    .plot.ssa.paired(x, ..., plot.contrib = plot.contrib, idx = idx, idy = idy);
  } else if (identical(type, "series")) {
    if (missing(groups))
      groups <- as.list(1:min(nlambda(x), nu(x)));

    .plot.ssa.series(x, ..., groups = groups);
  } else {
    stop("Unsupported type of SSA plot!");
  }
}

plot.1d.ssa.reconstruction <- function(x, ...,
                                       type = c("raw", "cumsum"),
                                       plot.method = c("native", "matplot"),
                                       base.series = NULL,
                                       add.original = TRUE,
                                       add.residuals = TRUE) {
  type <- match.arg(type);
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
  dots <- .defaults(dots, "main", "Reconstructed Series");
  dots <- .defaults(dots, "type", "l");
  dots <- .defaults(dots, "ylab", "");

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
