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
                    kind = c("1d-ssa", "2d-ssa", "toeplitz-ssa"),
                    svd_method = c("nutrlan", "propack", "svd", "eigen"),
                    force.decompose = TRUE) {
  svd_method <- match.arg(svd_method);
  kind <- match.arg(kind);
  xattr <- attributes(x);

  if (identical(kind, "1d-ssa") || identical(kind, "toeplitz-ssa")) {
    # Coerce input to vector if necessary
    if (!is.vector(x))
      x <- as.vector(x);

    N <- length(x);

    # Fix svd method, if needed
    if ((identical(svd_method, "nutrlan") ||
         identical(svd_method, "propack")) &&
        L < 50)
      svd_method <- "eigen";
  } else if (identical(kind, "2d-ssa")) {
    # Coerce input to matrix if necessary
    if (!is.matrix(x))
      x <- as.matrix(x);

    N <- dim(x);
  }

  # Create information body
  this <- list(length = N,
               window = L,
               call = match.call(),
               kind = kind,
               series = deparse(substitute(x)),
               svd_method = svd_method);

  # Create data storage
  this <- .create.storage(this);

  # Save series
  .set(this, "F", x);

  # Save attributes
  .set(this, "Fattr", xattr);
  
  # Make this S3 object
  class(this) <- c(paste(kind, svd_method, sep = "."), kind, "ssa");

  # Decompose, if necessary
  if (force.decompose)
    this <- decompose(this, ...);

  this;
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

  # Hack-hack-hack! Some routines will work much more efficiently if we'll
  # pass space to store some data which is needed to be calculated only once.
  e <- new.env();

  for (idx in new) {
    # Do actual reconstruction (depending on method, etc)
    .set.series(this,
                .do.reconstruct(this, idx, env = e), idx);
  }

  # Cleanup
  rm(list = ls(envir = e, all.names = TRUE),
     envir = e, inherits = FALSE);
  invisible(gc(verbose = FALSE));
}

cleanup.ssa <- function(this, ...) {
  .remove(this, ls(.storage(this), pattern = "series:"));
  invisible(gc(verbose = FALSE));
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

    # Propagate attributes (e.g. dimension for 2d-SSA)
    attributes(out[[i]]) <- .get(this, "Fattr");
  }

  # Cleanup
  rm(list = ls(envir = e, all.names = TRUE),
     envir = e, inherits = FALSE);
  gc(verbose = FALSE);

  names(out) <- paste("F", 1:length(groups), sep="");

  # Reconstructed series can be pretty huge...
  invisible(out);
}

.do.reconstruct <- function(this, idx, env = .GlobalEnv) {
  if (max(idx) > nlambda(this))
    stop("Too few eigentriples computed for this decompostion")

  lambda <- .get(this, "lambda");
  U <- .get(this, "U");

  res <- numeric(this$length);

  for (i in idx) {
    if (nv(this) >= i) {
      # FIXME: Check, whether we have factor vectors for reconstruction
      # FIXME: Get rid of .get call
      V <- .get(this, "V")[, i];
    } else {
      # No factor vectors available. Calculate them on-fly.
      V <- calc.v(this, i, env = env);
    }

    res <- res + lambda[i] * .hankelize.one(this, U = U[, i], V = V);
  }

  res;
}

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

.object.size.ssa <- function(this, pat = NULL) {
  env <- .storage(this);
  if (is.null(pat)) {
    members <- ls(envir = env, all.names = TRUE);
  } else {
    members <- ls(envir = env, pattern = pat);
  }

  l <- sapply(members, function(x) object.size(.get(this, x)))

  ifelse(length(l), sum(l), 0);
}

print.ssa <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n", sep="");
  cat("Series length:", paste(x$length, collapse = " x "));
  cat(",\tWindow length:", paste(x$window, collapse = " x "));
  cat(",\tSVD method:", x$svd_method);
  cat("\n\nComputed:\n");
  cat("Eigenvalues:", nlambda(x));
  cat(",\tEigenvectors:", nu(x));
  cat(",\tFactor vectors:", nv(x));
  cat("\n\nPrecached:",
      length(.get.series.info(x)),
      "subseries (")
  cat(format(.object.size(x, pat = "series:") / 1024 / 1024, digits = digits),
      "MiB)");
  cat("\n\nOverall memory consumption (estimate):",
      format(.object.size(x) / 1024 / 1024, digits = digits),
      "MiB");
  cat("\n");
  invisible(x);
}

summary.ssa <- function(object, digits = max(3, getOption("digits") - 3), ...)
  print.ssa(x = object, digits = digits, ...)

calc.v.ssa <- function(this, idx, env = .GlobalEnv, ...)
  stop("Unsupported SVD method for SSA!");

decompose.ssa <- function(x, neig = min(50, L, K), ...) {
  L <- K <- 0;
  stop("Unsupported SVD method for SSA!");
}
