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
  .set(this, "F", x);
  
  # Make this S3 object
  class(this) <- "ssa";

  # Decompose, if necessary
  if (force.decompose)
    decompose(this, ...);

  this;
}

# FIXME: add version, which accepts first column and last row
hankel <- function(X, L) {
  if (is.matrix(X) && nargs() == 1) {
     L <- dim(X)[1]; K <- dim(X)[2]; N <- K + L - 1;
     left  <- c(1:L, L*(2:K));
     right <- c(1+L*(0:(K-1)), ((K-1)*L+2):(K*L));
     v <- sapply(1:N, function(i) mean(X[seq.int(left[i], right[i], by = L-1)]));
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

.hankelize.one <- function(U, V) {
  storage.mode(U) <- storage.mode(V) <- "double";
  .Call("hankelize_one", U, V);
}

.hankelize.multi <- function(U, V) {
  storage.mode(U) <- storage.mode(V) <- "double";
  .Call("hankelize_multi", U, V);
}

.decompose.ssa.hankel <- function(this,
                                  nu = min(L, K), nv = min(L, K)) {
  N <- this$length; L <- this$window; K <- N - L + 1;
  F <- .get(this, "F");

  X <- hankel(F, L = L);

  # FIXME: Use special SVD for hankel matrixes
  S <- svd(X, nu = nu, nv = nv);

  # Save results
  .set(this, "lambda", S$d);
  if (!is.null(S$u))
    .set(this, "U", S$u);
  if (!is.null(S$v))
    .set(this, "V", S$v);
}

.decompose.ssa.toeplitz <- function(this, ...) {
  stop("Unimplemented!")
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

precache.ssa <- function(this, n, ...) {
  if (missing(n)) {
    warning("Amount of sub-series missed, precaching EVERYTHING",
            immediate. = TRUE);
    n <- nlambda(this);
  }

  # Calculate numbers of sub-series to be calculated
  info <- .get.series.info(this);
  new <- setdiff(1:n, info);

  # We're supporting only 'full' data for now
  U <- .get(this, "U");
  V <- .get(this, "V");
  lambda <- .get(this, "lambda");

  F <- .hankelize.multi(U[,new], V[,new]);

  # Return numbers of sub-series cached
  invisible(sapply(new,
                   function(i) {
                     .set.series(this,
                                 lambda[i] * F[,i],
                                 i)}));
}

.get.series.info <- function(this) {
  if (.exists(this, "series:info"))
    return (.get(this, "series:info"));

  numeric(0);
}

.append.series.info <- function(this, index) {
  .set(this, "series:info",
           union(.get.series.info(this), index));
}

.set.series <- function(this, F, index) {
  name <- paste("series:", index, sep = "");
  .set(this, name, F);
  .append.series.info(this, index);
  index;
}

.get.series <- function(this, index) {
  F <- numeric(this$length);
  for (i in index) {
    name <- paste("series:", i, sep = "");
    F <- F + .get(this, name);
  }
  F;
}

cleanup.ssa <- function(this, ...) {
  .remove(this, ls(.storage(this), pattern = "series:"));
  gc();
}

reconstruct.ssa <- function(this, groups, ..., cache = TRUE) {
  out <- list();
  nu <- nu(this); nv <- nv(this);

  if (missing(groups))
    groups <- as.list(1:nlambda(this));

  # We're supporting only 'full' data for now
  U <- .get(this, "U");
  V <- .get(this, "V");
  lambda <- .get(this, "lambda");

  # Grab indices of pre-cached values
  info <- .get.series.info(this);

  # Do actual reconstruction
  for (i in seq_along(groups)) {
    group <- groups[[i]];
    new <- setdiff(group, info);
    cached <- intersect(group, info);

    if (length(new) == 0) {
      # Nothing to compute, just create zero output
      out[[i]] <- numeric(this$length);
    } else if (length(new) == 1) {
      # Special case for rank one reconstruction
      out[[i]] <- lambda[new] * .hankelize.one(U[, new], V[, new]);

      # Cache the reconstructed series, if this was requested
      if (cache)
        .set.series(this, out[[i]], new);
    } else {
      out[[i]] <- hankel(U[, new] %*%
                         diag(lambda[new], nrow = length(group)) %*%
                         t(V[, new]));
    }
    # Add pre-cached series
    out[[i]] <- out[[i]] + .get.series(this, cached);
  }

  names(out) <- paste("F", 1:length(groups), sep="");

  # Reconstructed series can be pretty huge...
  invisible(out);
}

nu.ssa <- function(this, ...) {
  ifelse(.exists(this, "U"), dim(.get(this, "U"))[2], 0);
}

nv.ssa <- function(this, ...) {
  ifelse(.exists(this, "V"), dim(.get(this, "V"))[2], 0);
}

nlambda.ssa <- function(this, ...) {
  ifelse(.exists(this, "lambda"), length(.get(this, "lambda")), 0);
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

.get <- function(this, name) {
  get(name, envir = .storage(this));
}

.set <- function(this, name, value) {
  assign(name, value, envir = .storage(this), inherits = FALSE);
}

.exists <- function(this, name) {
  exists(name, envir = .storage(this), inherits = FALSE);
}

.remove <- function(this, name) {
  rm(list = name, envir = .storage(this), inherits = FALSE);
}

'$.ssa' <- function(this, name) {
  if (ind <- charmatch(name, names(this), nomatch = 0))
    return (this[[ind]]);

  if (.exists(this, name))
    return (.get(this, name));

  NULL;
}

# Generics
clone <- function(this, ...)
  UseMethod("clone");
reconstruct <- function(this, ...)
  UseMethod("reconstruct");
nu <- function(this, ...)
  UseMethod("nu");
nv <- function(this, ...)
  UseMethod("nv");
nlambda <- function(this, ...)
  UseMethod("nlambda");
precache <- function(this, ...)
  UseMethod("precache");
cleanup <- function(this, ...)
  UseMethod("cleanup");

# There is decompose() call in stats package, we need to take control over it
decompose <- function(this, ...) UseMethod("decompose");
decompose.default <- stats::decompose;
formals(decompose.default) <- c(formals(decompose.default), alist(... = ));

#.F <- function(x) exp(-.01 * x)*cos(x/100) + 0.05*rnorm(length(x));
# F <- .F(1:5000);
