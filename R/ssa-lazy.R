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

new.ssa <- function(x,
                    L = (N - 1) %/% 2,
                    ..., 
                    method = c("hankel", "toeplitz"),
                    centering = c("none", "row", "both"),
                    force.decompose = TRUE) {
  method <- match.arg(method);
  centering <- match.arg(centering);
  N <- length(x);

  # Create information body
  this <- list(length = N,
               window = L,
               call = match.call(),
               series = deparse(substitute(x)),
               method = method,
               centering = centering);

  # Create data storage
  attr(this, ".env") <- new.env();

  # Save series
  assign("F", x, envir = .storage(this), inherits = FALSE);
  
  # Make this S3 object
  class(this) <- "ssa";

  # Decompose, if necessary
  if (force.decompose) {
    decompose(this, ...);
  }

  this;
}

# FIXME: add version, which accepts first column and last row
hankel <- function(X, L) {
  if (is.matrix(X) && nargs() == 1) {
     L <- dim(X)[1]; K <- dim(X)[2]; N <- K + L - 1;
     left  <- c(1:L, L*(2:K));
     right <- c(1+L*(0:(K-1)), ((K-1)*L+2):(K*L));
     v <- sapply(1:N, function(i) mean(X[seq(left[i], right[i], by = L-1)]));
     return (v);
  }

  # Coerce output to vector, if necessary
  if (!is.vector(X))
    X <- as.vector(X);
  N <- length(X);
  if (missing(L))
    L <- (N + 1) %/% 2;
  K <- N - L + 1;
  outer(1:L, 1:K, function(x,y) X[x+y-1]);
}

.decompose.ssa.hankel <- function(this,
                                  nu = min(L, K), nv = min(L, K)) {
  N <- this$length; L <- this$window; K <- N - L + 1;
  F <- get("F", envir = .storage(this));

  X <- hankel(F, L = L);

  # FIXME: Use special SVD for hankel matrixes
  S <- svd(X, nu = nu, nv = nv);

  # Save results
  assign("lambda", S$d, envir = .storage(this), inherits = FALSE);
  if (!is.null(S$u)) {
    assign("U", S$u, envir = .storage(this), inherits = FALSE);
  }
  if (!is.null(S$v)) {
    assign("V", S$v, envir = .storage(this), inherits = FALSE);
  }
}

decompose.ssa <- function(this, ...) {
  method <- this$method;

  if (identical(method, "hankel")) {
    .decompose.ssa.hankel(this, ...);
  } else if (identical(method, "toeplitz")) {
    .decompose.ssa.toeplitz(this, ...);
  } else {
    stop("Unknown method in SSA")
  }
}

reconstruct.ssa <- function(this, groups, ...) {
  out <- list();
  nu <- nu(this); nv <- nv(this);

  if (missing(groups))
    groups <- as.list(1:nlambda(this));

  # We're supporting only 'full' data for now
  U <- get("U", envir = .storage(this));
  V <- get("V", envir = .storage(this));
  lambda <- get("lambda", envir = .storage(this));

  for (i in seq_along(groups)) {
    group <- groups[[i]];

    out[[i]] <- hankel(U[, group] %*%
                       diag(lambda[group], nrow = length(group)) %*%
                       t(V[,group]));
  }

  names(out) <- paste("F", 1:length(groups), sep="");

  # Reconstructed series can be pretty huge...
  invisible(out);
}

nu.ssa <- function(this, ...) {
  ifelse(exists("U", envir = .storage(this), inherits = FALSE),
         dim(get("U", envir = .storage(this)))[2],
         0);
}

nv.ssa <- function(this, ...) {
  ifelse(exists("V", envir = .storage(this), inherits = FALSE),
         dim(get("V", envir = .storage(this)))[2],
         0);
}

nlambda.ssa <- function(this, ...) {
  ifelse(exists("lambda", envir = .storage(this), inherits = FALSE),
         dim(get("lambda", envir = .storage(this)))[2],
         0);
}

clone.ssa <- function(this, ...) {
  # Copy the information body
  obj <- this;

  # Make new storage
  clone.env <- new.env();
  this.env <- .storage(this);
  attr(obj, ".env") <- clone.env;

  # Copy the contents of data storage
  for (field in ls(envir = this.env, all.names = TRUE)) {
    value <- get(field, envir = this.env, inherits = FALSE);
    assign(field, value, envir = clone.env, inherits = FALSE);
  }

  obj;
}

.storage <- function(this) {
  attr(this, ".env");
}

'$.ssa' <- function(this, name) {
  if (ind <- charmatch(name, names(this), nomatch = 0))
    return (this[[ind]]);

  if (exists(name, envir = .storage(this), inherits = FALSE))
    return (get(name, envir = .storage(this)));

  NULL;
}

clone <- function(this, ...) {
  UseMethod("clone");
}

decompose <- function(this, ...) {
  UseMethod("decompose");
}

reconstruct <- function(this, ...) {
  UseMethod("reconstruct");
}

nu <- function(this, ...) {
  UseMethod("nu");
}

nv <- function(this, ...) {
  UseMethod("nv");
}

nlambda <- function(this, ...) {
  UseMethod("nlambda");
}
