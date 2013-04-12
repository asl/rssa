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
  V <- ssaobj$U[,y];
  U <- if (identical(x, y)) 1:length(V)
       else ssaobj$U[,x]

  res <- prepanel.default.xyplot(U, V);
  if (symmetric) {
    res$ylim <- range(V, -V);
    if (!identical(x, y))
      res$xlim <- range(U, -U);
  }

  res;
}

panel.eigenvectors <- function(x, y, ssaobj, ..., ref = FALSE) {
  V <- ssaobj$U[,y];
  U <- if (identical(x, y)) 1:length(V)
       else ssaobj$U[,x]

  if (ref) {
    panel.abline(h = 0, ..., reference = TRUE)
    if (!identical(x, y))
      panel.abline(v = 0, ..., reference = TRUE)
  }

  panel.xyplot(U, V, ...);
}

.defaults <- function(v, key, value) {
  if (!(key %in% names(v))) v[[key]] <- value;
  v;
}

.plot.ssa.values <- function(x, ..., numvalues, plot.type = "b") {
  dots <- list(...);

  # FIXME: check for proper lengths
  d <- data.frame(A = 1:numvalues, B = x$lambda[1:numvalues]);

  # Provide convenient defaults
  dots <- .defaults(dots, "type", plot.type);
  dots <- .defaults(dots, "xlab", "Index");
  dots <- .defaults(dots, "ylab", "log of eigenvalue");
  dots <- .defaults(dots, "main", "Eigenvalues");
  dots <- .defaults(dots, "scales", list(y = list(log = TRUE)));
  dots <- .defaults(dots, "grid", TRUE)
  dots <- .defaults(dots, "par.settings", list())
  dots$par.settings <- .defaults(dots$par.settings, "plot.symbol", list(pch = 20))

  res <- do.call("xyplot",
                 c(list(x = B ~ A , data = d, ssaobj = x), dots));
  print(res)
}

.plot.ssa.vectors <- function(x, ..., plot.contrib, idx, plot.type = "l") {
  dots <- list(...);

  # FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idx);

  if (plot.contrib) {
    total <- sum(x$lambda^2);
    lambda <- round(100*x$lambda[idx]^2 / total, digits = 2);
  }

  # Provide convenient defaults
  dots <- .defaults(dots, "type", plot.type);
  dots <- .defaults(dots, "xlab", "");
  dots <- .defaults(dots, "ylab", "");
  dots <- .defaults(dots, "main", "Eigenvectors");
  dots <- .defaults(dots, "as.table", TRUE);
  dots <- .defaults(dots, "scales", list(draw = FALSE, relation = "free"));
  dots <- .defaults(dots, "symmetric", TRUE);
  dots <- .defaults(dots, "ref", TRUE);

  res <- do.call("xyplot",
                 c(list(x = A ~ B | factor(A,
                                           labels = if (!plot.contrib) A else paste(A, " (", lambda, "%)", sep = "")),
                        data = d, ssaobj = x,
                        panel = panel.eigenvectors,
                        prepanel = prepanel.eigenvectors),
                   dots));
  print(res)
}

.plot.ssa.paired <- function(x, ..., plot.contrib, idx, idy, plot.type = "l") {
  dots <- list(...);

  # FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idy);

  if (plot.contrib) {
    total <- sum(x$lambda^2);
    lambdax <- round(100*x$lambda[idx]^2 / total, digits = 2);
    lambday <- round(100*x$lambda[idy]^2 / total, digits = 2);
  }

  # Provide convenient defaults
  dots <- .defaults(dots, "type", plot.type);
  dots <- .defaults(dots, "xlab", "");
  dots <- .defaults(dots, "ylab", "");
  dots <- .defaults(dots, "main", "Pairs of eigenvectors");
  dots <- .defaults(dots, "as.table", TRUE);
  dots <- .defaults(dots, "scales", list(draw = FALSE, relation = "free"));
  dots <- .defaults(dots, "aspect", "x");
  dots <- .defaults(dots, "symmetric", TRUE);
  dots <- .defaults(dots, "ref", TRUE);

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

.plot.ssa.series <- function(x, ..., groups, plot.type = "l") {
  dots <- list(...);

  # FIXME: check for proper lengths
  idx <- seq_along(groups);
  d <- data.frame(A = idx, B = idx);

  r <- reconstruct(x, groups = groups, drop = FALSE);

  # Provide convenient defaults
  dots <- .defaults(dots, "type", plot.type);
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

  if (is.list(x)) {
    labels <- if (is.null(names(x))) rep("", length(x)) else names(x)
    labels[labels == ""] <- paste("W", seq_len(sum(labels == "")), sep = "")

    data <- lapply(seq_along(x), function(i) {
          mx <- x[[i]]
          grid <- expand.grid(row = seq_len(nrow(mx)), column = seq_len(ncol(mx)))
          grid$x <- as.vector(as.numeric(mx))
          grid$i <- i

          grid
        })
    data <- do.call("rbind", data)
    data$i <- factor(data$i, labels = labels)
    formula <- abs(x) ~ row * column | i
    dots <- .defaults(dots, "scales", list(relation = "free"))
    dots$scales <- .defaults(dots$scales, "relation", "free")
    dots <- .defaults(dots, "aspect", 1)
    dots <- .defaults(dots, "xlim", lapply(x, rownames))
    dots <- .defaults(dots, "ylim", lapply(x, colnames))
  } else {
    data <- expand.grid(row = seq_len(nrow(x)), column = seq_len(ncol(x)))
    data$x <- as.vector(as.numeric(x))
    formula <- abs(x) ~ row * column
    dots <- .defaults(dots, "aspect", "iso")
    dots <- .defaults(dots, "xlim", rownames(x))
    dots <- .defaults(dots, "ylim", colnames(x))
  }

  res <- do.call("levelplot", c(list(formula,
                                     data = data,
                                     at = seq(zlim[1], zlim[2], length.out = cuts),
                                     panel = panel.levelplot.wcor,
                                     grid = grid),
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
  } else if (identical(type, "wcor")) {
    if (missing(groups))
      groups <- as.list(1:min(nlambda(x), nu(x)));

    plot(wcor(x, groups = groups), ...)
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
