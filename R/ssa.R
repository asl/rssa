#   R package for Singular Spectrum Analysis
#   Copyright (c) 2008, 2009 Anton Korobeynikov <asl@math.spbu.ru>
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
                    method = c("nutrlan", "propack", "svd", "eigen"),
                    force.decompose = TRUE) {
  method <- match.arg(method);
  N <- length(x);

  # Create information body
  this <- list(length = N,
               window = L,
               call = match.call(),
               series = deparse(substitute(x)),
               method = method);

  # Create data storage
  this <- .create.storage(this);

  # Save series
  .set(this, "F", x);
  
  # Make this S3 object
  class(this) <- c(paste("ssa", method, sep = "."), "ssa");

  # Decompose, if necessary
  if (force.decompose)
    decompose(this, ...);

  this;
}

decompose.ssa.svd <- function(this,
                              neig = min(L, K),
                              ...,
                              force.continue = FALSE) {
  N <- this$length; L <- this$window; K <- N - L + 1;

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(this) > 0)
    stop("Continuation of decompostion is not supported for this method.")

  # Build hankel matrix
  F <- .get(this, "F");
  h <- hankel(F, L = L);

  # Do decomposition
  S <- svd(h, nu = neig, nv = neig);

  # Save results
  .set(this, "lambda", S$d);
  if (!is.null(S$u))
    .set(this, "U", S$u);
  if (!is.null(S$v))
    .set(this, "V", S$v);
}

decompose.ssa.eigen <- function(this, ...,
                                force.continue = FALSE) {
  N <- this$length; L <- this$window; K <- N - L + 1;

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(this) > 0)
    stop("Continuation of decompostion is not supported for this method.")

  # Build hankel matrix (this can be done more efficiently!)
  F <- .get(this, "F");
  h <- hankel(F, L = L);

  # Do decomposition
  if ("neig" %in% names(list(...)))
    warning("'neig' option ignored for SSA method 'eigen', computing EVERYTHING",
            immediate. = TRUE)
  
  S <- eigen(tcrossprod(h));

  # Fix small negative values
  S$values[S$values < 0] <- 0;
  
  # Save results
  .set(this, "lambda", sqrt(S$values));
  .set(this, "U", S$vectors);
}

decompose.ssa.propack <- function(this,
                                  neig = min(50, L, K),
                                  ...,
                                  force.continue = FALSE) {
  N <- this$length; L <- this$window; K <- N - L + 1;

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(this) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.")

  F <- .get(this, "F");
  h <- new.hmat(F, L = L);

  S <- propack_svd(h, neig = neig, ...);

  # Save results
  .set(this, "hmat", h);
  .set(this, "lambda", S$d);
  if (!is.null(S$u))
    .set(this, "U", S$u);
  if (!is.null(S$v))
    .set(this, "V", S$v);
}

decompose.ssa.nutrlan <- function(this,
                                  neig = min(50, L, K),
                                  ...) {
  N <- this$length; L <- this$window; K <- N - L + 1;

  h <- .get(this, "hmat", allow.null = TRUE);
  if (is.null(h)) {
    F <- .get(this, "F");
    h <- new.hmat(F, L = L);
  }

  lambda <- .get(this, "lambda", allow.null = TRUE);
  U <- .get(this, "U", allow.null = TRUE);

  S <- trlan_svd(h, neig = neig, ...,
                 lambda = lambda, U = U);

  # Save results
  .set(this, "hmat", h);
  .set(this, "lambda", S$d);
  if (!is.null(S$u))
    .set(this, "U", S$u);
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

  if (missing(groups))
    groups <- as.list(1:min(nlambda(this), nu(this)));

  # Determine the upper bound of desired eigentriples
  desired <- max(unlist(groups));

  # Continue decomposition, if necessary
  if (desired > min(nlambda(this), nu(this)))
    decompose(this, ..., neig = desired);

  # Grab indices of pre-cached values
  info <- .get.series.info(this);

  # Hack-hack-hack! Some routines will work much more efficiently if we'll
  # pass space to store some data which is needed to be calculated only once.
  e <- new.env();

  # Do actual reconstruction
  for (i in seq_along(groups)) {
    group <- groups[[i]];
    new <- setdiff(group, info);
    cached <- intersect(group, info);

    if (length(new) == 0) {
      # Nothing to compute, just create zero output
      out[[i]] <- numeric(this$length);
    } else {
      # Do actual reconstruction (depending on method, etc)
      out[[i]] <- .do.reconstruct(this, new, env = e);

      # Cache the reconstructed series, if this was requested
      if (cache && length(new) == 1)
        .set.series(this, out[[i]], new);
    }

    # Add pre-cached series
    out[[i]] <- out[[i]] + .get.series(this, cached);
  }

  # Cleanup
  rm(list = ls(envir = e, all.names = TRUE),
     envir = e, inherits = FALSE);
  gc();

  names(out) <- paste("F", 1:length(groups), sep="");

  # Reconstructed series can be pretty huge...
  invisible(out);
}

.do.reconstruct <- function(this, idx, env = .GlobalEnv) {
  if (max(idx) > nlambda(this))
    stop("Too few eigentriples computed for this decompostion")

  lambda <- .get(this, "lambda");
  U <- .get(this, "U");

  if (nv(this) > 0) {
    # Check, whether we have factor vectors for reconstruction
    V <- .get(this, "V");

    if (length(idx) == 1) {
      # Special case for rank one reconstruction
      res <- lambda[idx] * .hankelize.one(U[, idx], V[, idx]);
    } else {
      # This won't work for lengthy series. Consider fixing :)
      res <- hankel(U[, idx] %*%
                    diag(lambda[idx], nrow = length(idx)) %*%
                    t(V[, idx]));
    }
  } else {
    # No factor vectors available. Calculate them on-fly.
    # FIXME: Should we consider caching them? Per-request?
    res <- numeric(this$length);

    for (i in idx) {
      V <- calc.v(this, i, env = env);
      res <- res + lambda[i] * .hankelize.one(U[, i], V);
    }
  }

  res;
}

.calc.v.hankel <- function(this, idx) {
  lambda <- .get(this, "lambda")[idx];
  U <- .get(this, "U")[, idx, drop = FALSE];
  h <- .get(this, "hmat");

  invisible(sapply(1:length(idx),
                   function(i) hmatmul(h, U[, i], transposed = TRUE) / lambda[i]));
}

.calc.v.svd <- function(this, idx, env) {
  # Check, if there is garbage-collected storage to hold some pre-calculated
  # stuff.
  if (identical(env, .GlobalEnv) ||
      !exists(".ssa.temporary.storage", envir = env, inherits = FALSE)) {
    F <- .get(this, "F");

    # Build hankel matrix.
    X <- hankel(F, L = this$window);

    # Save to later use, if possible.
    if (!identical(env, .GlobalEnv)) {
      assign(".ssa.temporary.storage", X, envir = env, inherits = FALSE);
    }
  } else {
    X <- get(".ssa.temporary.storage", envir = env, inherits = FALSE);
  }

  lambda <- .get(this, "lambda")[idx];
  U <- .get(this, "U")[, idx, drop = FALSE];

  invisible(sapply(1:length(idx),
                   function(i) crossprod(X, U[, i]) / lambda[i]));
}

calc.v.ssa.nutrlan <- function(this, idx, env = .GlobalEnv) .calc.v.hankel(this, idx)
calc.v.ssa.propack <- function(this, idx, env = .GlobalEnv) .calc.v.hankel(this, idx)
calc.v.ssa.svd <- function(this, idx, env = .GlobalEnv) .calc.v.svd(this, idx, env)
calc.v.ssa.eigen <- function(this, idx, env = .GlobalEnv) .calc.v.svd(this, idx, env)

nu.ssa <- function(this, ...) {
  ifelse(.exists(this, "U"), ncol(.get(this, "U")), 0);
}

nv.ssa <- function(this, ...) {
  ifelse(.exists(this, "V"), ncol(.get(this, "V")), 0);
}

nlambda.ssa <- function(this, ...) {
  ifelse(.exists(this, "lambda"), length(.get(this, "lambda")), 0);
}

clone.ssa <- function(this, copy.cache = TRUE, ...) {
  obj <- .clone(this);
  if (copy.cache == FALSE)
    cleanup(obj);

  obj;
}

clusterify.ssa <- function(this, groups, nclust = length(groups) / 2,
                           ...,
                           type = c("wcor"), cache = TRUE) {
  type <- match.arg(type)

  if (missing(groups))
    groups <- as.list(1:nlambda(this));

  if (identical(type, "wcor")) {
    w <- wcor(this, groups = groups, ..., cache = cache);
    g <- clusterify(w, nclust = nclust, ...);
    out <- lapply(g, function(idx) unlist(groups[idx]));
  } else {
    stop("Unsupported clusterification method!");
  }
  out;
}

'$.ssa' <- function(this, name) {
  if (ind <- charmatch(name, names(this), nomatch = 0))
    return (this[[ind]]);

  if (.exists(this, name))
    return (.get(this, name));

  NULL;
}

print.ssa <- function(this, ...) {
  cat("\nCall:\n", deparse(this$call), "\n\n", sep="");
  cat("Series length:", this$length);
  cat(",\tWindow length:", this$window);
  cat("\n\nComputed:\n");
  cat("Eigenvalues:", nlambda(this));
  cat(",\tEigenvectors:", nu(this));
  cat(",\tFactor vectors:", nv(this));
  cat("\n\nPrecached:", length(.get.series.info(this)));
  cat(" subseries");
  cat("\n");
  invisible(this);
}

#.F <- function(x) exp(-.01 * x)*cos(x/100) + 0.05*rnorm(length(x));
# F <- .F(1:5000);
