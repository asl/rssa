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

.plot.ssa.values <- function(this, ..., numvalues) {
  dots <- list(...);

  # FIXME: check for proper lengths
  d <- data.frame(A = 1:numvalues, B = this$lambda[1:numvalues]);

  # Provide convenient defaults
  dots <- .defaults(dots, "type", c("b", "g"));
  dots <- .defaults(dots, "xlab", "Index");
  dots <- .defaults(dots, "ylab", "log of eigenvalue");
  dots <- .defaults(dots, "main", "Eigenvalues");
  dots <- .defaults(dots, "scales", list(y = list(log = TRUE)));
  dots <- .defaults(dots, "pch", 20);

  do.call("xyplot",
          c(list(x = B ~ A , data = d, ssaobj = this), dots));
}

.plot.ssa.vectors <- function(this, ..., plot.contrib, idx) {
  dots <- list(...);

  # FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idx);

  if (plot.contrib) {
    total <- sum(this$lambda);
    lambda <- round(100*this$lambda[idx] / total, digits = 2);
  }

  # Provide convenient defaults
  dots <- .defaults(dots, "type", "l");
  dots <- .defaults(dots, "xlab", "");
  dots <- .defaults(dots, "ylab", "");
  dots <- .defaults(dots, "main", "Eigenvectors");
  dots <- .defaults(dots, "as.table", TRUE);
  dots <- .defaults(dots, "scales", list(relation = "free"));

  do.call("xyplot",
          c(list(x = A ~ B | factor(A,
                                    labels = if (!plot.contrib) A else paste(A, " (", lambda, "%)", sep = "")),
                 data = d, ssaobj = this,
                 panel = panel.eigenvectors,
                 prepanel = prepanel.eigenvectors),
            dots));
}

.plot.ssa.paired <- function(this, ..., plot.contrib, idx, idy) {
  dots <- list(...);

  # FIXME: check for proper lengths
  d <- data.frame(A = idx, B = idy);

  if (plot.contrib) {
    total <- sum(this$lambda);
    lambdax <- round(100*this$lambda[idx] / total, digits = 2);
    lambday <- round(100*this$lambda[idy] / total, digits = 2);
  }

  # Provide convenient defaults
  dots <- .defaults(dots, "type", "l");
  dots <- .defaults(dots, "xlab", "");
  dots <- .defaults(dots, "ylab", "");
  dots <- .defaults(dots, "main", "Pairs of eigenvectors");
  dots <- .defaults(dots, "as.table", TRUE);
  dots <- .defaults(dots, "scales", list(relation = "free"));

  do.call("xyplot",
          c(list(x = A ~ B | factor(A,
                                    labels = if (!plot.contrib) paste(A, "vs", B)
                                             else paste(A, " (", lambdax, "%) vs ", B, " (", lambday, "%)", sep = "")),
                 data = d, ssaobj = this,
                 panel = panel.eigenvectors,
                 prepanel = prepanel.eigenvectors),
            dots));
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

.plot.ssa.series <- function(this, ..., groups) {
  dots <- list(...);

  # FIXME: check for proper lengths
  idx <- seq_along(groups);
  d <- data.frame(A = idx, B = idx);

  r <- reconstruct(this, groups = groups);

  # Provide convenient defaults
  dots <- .defaults(dots, "type", "l");
  dots <- .defaults(dots, "xlab", "");
  dots <- .defaults(dots, "ylab", "");
  dots <- .defaults(dots, "main", "Reconstructed series");
  dots <- .defaults(dots, "as.table", TRUE);
  dots <- .defaults(dots, "scales", list(relation = "free"));

  do.call("xyplot",
          c(list(x = A ~ B | factor(A, labels = paste(groups)),
                 data = d, recon = r,
                 panel = panel.series,
                 prepanel = prepanel.series),
            dots));

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
  this <- x;
 
  if (identical(type, "values")) {
    .plot.ssa.values(this, ..., numvalues = numvalues);
  } else if (identical(type, "vectors")) {
    .plot.ssa.vectors(this, ..., plot.contrib = plot.contrib, idx = idx);
  } else if (identical(type, "paired")) {
    if (missing(idy))
      idy <- idx + 1;

    .plot.ssa.paired(this, ..., plot.contrib = plot.contrib, idx = idx, idy = idy);
  } else if (identical(type, "series")) {
    if (missing(groups))
      groups <- as.list(1:min(nlambda(this), nu(this)));

    .plot.ssa.series(this, ..., groups = groups);
  } else {
    stop("Unsupported type of SSA plot!");
  }
}
