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

eigenplot.ssa <- function(s,
                          kind = c("values", "vectors", "paired"),
                          ...,
                          plot.contrib = TRUE,
                          numvalues = nlambda(s),
                          numvectors = min(nlambda(s), 10),
                          idx = 1:numvectors,
                          idy) {
  kind <- match.arg(kind);
  dots <- list(...);

  .defaults <- function(v, key, value) {
    if (!(key %in% names(v))) v[[key]] <- value;
    v;
  }
  
  if (identical(kind, "values")) {
    # FIXME: check for proper lengths
    d <- data.frame(A = 1:numvalues, B = s$lambda[1:numvalues]);

    # Provide convenient defaults
    dots <- .defaults(dots, "type", c("b", "g"));
    dots <- .defaults(dots, "xlab", "Index");
    dots <- .defaults(dots, "ylab", "log of eigenvalue");
    dots <- .defaults(dots, "main", "Eigenvalues");
    dots <- .defaults(dots, "scales", list(y = list(log = TRUE)));
    dots <- .defaults(dots, "pch", 20);

    do.call("xyplot",
            c(list(x = B ~ A , data = d, ssaobj = s), dots));
  } else if (identical(kind, "vectors")) {
    # FIXME: check for proper lengths
    d <- data.frame(A = idx, B = idx);

    if (plot.contrib) {
      total <- sum(s$lambda);
      lambda <- round(100*s$lambda[idx] / total, digits = 2);
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
                   data = d, ssaobj = s,
                   panel = panel.eigenvectors,
                   prepanel = prepanel.eigenvectors),
              dots));
  } else if (identical(kind, "paired")) {
    if (missing(idy))
      idy <- idx + 1;

    # FIXME: check for proper lengths
    d <- data.frame(A = idx, B = idy);

    if (plot.contrib) {
      total <- sum(s$lambda);
      lambdax <- round(100*s$lambda[idx] / total, digits = 2);
      lambday <- round(100*s$lambda[idy] / total, digits = 2);
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
                   data = d, ssaobj = s,
                   panel = panel.eigenvectors,
                   prepanel = prepanel.eigenvectors),
              dots));
  } else {
    stop("Unsupported kind of eigenplot!");
  }
}
