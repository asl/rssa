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

fix.svd.method <- function(svd.method, L, N, ...) {
  dots <- list(...)
  neig <- dots$neig
  truncated <- (identical(svd.method, "nutrlan") || identical(svd.method, "propack"))

  if (is.null(neig)) neig <- min(50, L, N - L + 1)
  if (truncated) {
    # It's not wise to call truncated methods for small matrices at all
    if (L < 50) {
      truncated <- FALSE
      svd.method <- "eigen"
    } else if (neig > L /2) {
    # Check, whether desired eigentriples amount is too huge
      if (L < 200) {
        svd.method <- "eigen"
        truncated <- FALSE
      } else {
        warning("too many eigentriples requested")
      }
    }
  }

  svd.method
}

new.ssa <- function(...) {
  warning("`new.ssa' method is deprecated, use `ssa' instead")
  ssa(...)
}

ssa <- function(x,
                L = (N + 1) %/% 2,
                ...,
                kind = c("1d-ssa", "2d-ssa", "toeplitz-ssa"),
                svd.method = c("nutrlan", "propack", "svd", "eigen"),
                force.decompose = TRUE) {
  svd.method <- match.arg(svd.method);
  kind <- match.arg(kind);
  xattr <- attributes(x);

  if (identical(kind, "1d-ssa") || identical(kind, "toeplitz-ssa")) {
    # Coerce input to vector if necessary
    if (!is.vector(x))
      x <- as.vector(x);

    N <- length(x);

    # Fix svd method, if needed
    svd.method <- fix.svd.method(svd.method, L, N, ...)
  } else if (identical(kind, "2d-ssa")) {
    # Coerce input to matrix if necessary
    if (!is.matrix(x))
      x <- as.matrix(x);

    N <- dim(x);
  }

  # Normalized the kind to be used
  kind <- sub("-", ".", kind, fixed = TRUE)

  # Create information body
  this <- list(length = N,
               window = L,
               call = match.call(),
               kind = kind,
               series = deparse(substitute(x)),
               svd.method = svd.method);

  # Create data storage
  this <- .create.storage(this);

  # Save series
  .set(this, "F", x);

  # Save attributes
  .set(this, "Fattr", xattr);

  # Make this S3 object
  class(this) <- c(paste(kind, svd.method, sep = "."), kind, "ssa");

  # Perform additional init steps, if necessary
  .init(this)

  # Decompose, if necessary
  if (force.decompose)
    this <- decompose(this, ...);

  this;
}

.init.default <- function(x, ...) {
  # Do nothing
  x
}

.maybe.continue <- function(x, groups, ...) {
  L <- x$window
  K <- x$length - x$window + 1

  # Determine the upper bound of desired eigentriples
  desired <- max(unlist(groups))

  # Sanity check
  if (desired > min(prod(L), prod(K)))
    stop("Cannot decompose that much, desired elementary series index is too huge")

  # Continue decomposition, if necessary
  if (desired > min(nlambda(x), nu(x)))
    decompose(x, ..., neig = min(desired + 1, prod(L), prod(K)))

  desired
}

precache <- function(x, n, ...) {
  if (missing(n)) {
    warning("Amount of sub-series missed, precaching EVERYTHING",
            immediate. = TRUE);
    n <- nlambda(x);
  }

  # Calculate numbers of sub-series to be calculated
  info <- .get.series.info(x);
  new <- setdiff(1:n, info);

  # Hack-hack-hack! Some routines will work much more efficiently if we'll
  # pass space to store some data which is needed to be calculated only once.
  e <- new.env();

  for (idx in new) {
    # Do actual reconstruction (depending on method, etc)
    .set.series(x,
                .elseries(x, idx, env = e), idx);
  }

  # Cleanup
  rm(list = ls(envir = e, all.names = TRUE),
     envir = e, inherits = FALSE);
}

cleanup <- function(x) {
  .remove(x, ls(.storage(x), pattern = "series:"));
}

reconstruct.ssa <- function(x, groups, ...,
                            drop.attributes = FALSE, cache = TRUE) {
  out <- list();

  if (missing(groups))
    groups <- as.list(1:min(nlambda(x), nu(x)));

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  # Grab indices of pre-cached values
  info <- .get.series.info(x);

  # Hack-hack-hack! Some routines will work much more efficiently if we'll
  # pass space to store some data which is needed to be calculated only once.
  e <- new.env();

  # Do actual reconstruction. Calculate the residuals on the way
  residuals <- .get(x, "F")
  for (i in seq_along(groups)) {
    group <- groups[[i]];
    new <- setdiff(group, info);
    cached <- intersect(group, info);

    if (length(new) == 0) {
      # Nothing to compute, just create zero output
      out[[i]] <- numeric(prod(x$length));
    } else {
      # Do actual reconstruction (depending on method, etc)
      out[[i]] <- .elseries(x, new, env = e);

      # Cache the reconstructed series, if this was requested
      if (cache && length(new) == 1)
        .set.series(x, out[[i]], new);
    }

    # Add pre-cached series
    out[[i]] <- out[[i]] + .get.series(x, cached);

    # Propagate attributes (e.g. dimension for 2d-SSA)
    attributes(out[[i]]) <- .get(x, "Fattr");
  }

  # Set names and drop the dimension, if necessary
  names(out) <- paste("F", 1:length(groups), sep="");

  # Calculate the residuals
  residuals <- .get(x, "F")
  rgroups <- unique(unlist(groups))
  info <- .get.series.info(x);
  rcached <- intersect(rgroups, info)
  rnew <- setdiff(rgroups, info)
  residuals <- residuals - .get.series(x, rcached)
  if (length(rnew))
    residuals <- residuals - .elseries(x, rnew, env = e)

  # Propagate attributes of residuals
  attributes(residuals) <- .get(x, "Fattr");
  F <- .get(x, "F")
  if (!drop.attributes)
    attributes(F) <- .get(x, "Fattr")

  # Cleanup
  rm(list = ls(envir = e, all.names = TRUE),
     envir = e, inherits = FALSE);

  attr(out, "residuals") <- residuals;
  attr(out, "series") <- F;

  # Reconstructed series can be pretty huge...
  class(out) <- paste(c(x$kind, "ssa"), "reconstruction", sep = ".")
  invisible(out);
}

residuals.ssa <- function(object, groups, ..., cache = TRUE) {
  groups <- list(if (missing(groups)) 1:min(nlambda(object), nu(object)) else unlist(groups))

  residuals(reconstruct(object, groups = groups, ..., cache = cache))
}

residuals.ssa.reconstruction <- function(object, ...) {
  attr(object, "residuals")
}

.elseries.default <- function(x, idx, ..., env = .GlobalEnv) {
  if (max(idx) > nlambda(x))
    stop("Too few eigentriples computed for this decomposition")

  lambda <- .get(x, "lambda");
  U <- .get(x, "U");

  res <- numeric(prod(x$length));

  for (i in idx) {
    if (nv(x) >= i) {
      # FIXME: Check, whether we have factor vectors for reconstruction
      # FIXME: Get rid of .get call
      V <- .get(x, "V")[, i];
    } else {
      # No factor vectors available. Calculate them on-fly.
      V <- calc.v(x, i, env = env);
    }

    res <- res + lambda[i] * .hankelize.one(x, U = U[, i], V = V);
  }

  res;
}

nu <- function(x) {
  ifelse(.exists(x, "U"), ncol(.get(x, "U")), 0);
}

nv <- function(x) {
  ifelse(.exists(x, "V"), ncol(.get(x, "V")), 0);
}

nlambda <- function(x) {
  ifelse(.exists(x, "lambda"), length(.get(x, "lambda")), 0);
}

clone.ssa <- function(x, copy.storage = TRUE, copy.cache = TRUE, ...) {
  obj <- .clone(x, copy.storage = copy.storage);
  if (copy.cache == FALSE)
    cleanup(obj);

  obj;
}

clusterify.ssa <- function(x, group, nclust = length(group) / 2,
                           ...,
                           type = c("wcor"), cache = TRUE) {
  type <- match.arg(type)

  if (missing(group))
    group <- as.list(1:nlambda(x));

  if (identical(type, "wcor")) {
    w <- wcor(x, groups = group, ..., cache = cache);
    g <- clusterify(w, nclust = nclust, ...);
    out <- lapply(g, function(idx) unlist(group[idx]));
  } else {
    stop("Unsupported clusterification method!");
  }
  out;
}

'$.ssa' <- function(x, name) {
  if (ind <- charmatch(name, names(x), nomatch = 0))
    return (x[[ind]]);

  if (.exists(x, name))
    return (.get(x, name));

  NULL;
}

.object.size <- function(x, pat = NULL) {
  env <- .storage(x);
  if (is.null(pat)) {
    members <- ls(envir = env, all.names = TRUE);
  } else {
    members <- ls(envir = env, pattern = pat);
  }

  l <- sapply(members, function(el) object.size(.get(x, el)))

  ifelse(length(l), sum(l), 0);
}

print.ssa <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n", sep="");
  cat("Series length:", paste(x$length, collapse = " x "));
  cat(",\tWindow length:", paste(x$window, collapse = " x "));
  cat(",\tSVD method:", x$svd.method);
  cat("\n\nComputed:\n");
  cat("Eigenvalues:", nlambda(x));
  cat(",\tEigenvectors:", nu(x));
  cat(",\tFactor vectors:", nv(x));
  cat("\n\nPrecached:",
      length(.get.series.info(x)),
      "elementary series (")
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

wnorm.ssa <- function(x, ...)
  stop("`wnorm' is not implemented for this kind of SSA")
