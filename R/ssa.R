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

determine.svd.method <- function(L, K, neig = NULL, ..., svd.method = "nutrlan") {
  truncated <- (identical(svd.method, "nutrlan") || identical(svd.method, "propack"))

  if (is.null(neig)) neig <- min(50, L, K)
  if (truncated) {
    # It's not wise to call truncated methods for small matrices at all
    if (L < 500) {
      truncated <- FALSE
      svd.method <- "eigen"
    } else if (neig > L /2) {
      # Check, whether desired eigentriples amount is too huge
      if (L < 500) {
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
                neig = NULL,
                mask = NULL,
                wmask = NULL,
                ...,
                kind = c("1d-ssa", "2d-ssa", "toeplitz-ssa", "mssa", "cssa"),
                svd.method = c("auto", "nutrlan", "propack", "svd", "eigen"),
                force.decompose = TRUE) {
  svd.method <- match.arg(svd.method)
  kind <- match.arg(kind)
  xattr <- attributes(x)
  iattr <- NULL
  # Grab class separately. This way we will capture the inherit class as well
  xclass <- class(x)

  # Do the fixups depending on the kind of SSA.
  if (identical(kind, "1d-ssa") || identical(kind, "toeplitz-ssa")) {
    # Coerce input to vector if necessary
    if (!is.vector(x))
      x <- as.vector(x)

    N <- length(x)

    if (is.null(neig))
      neig <- min(50, L, N - L + 1)

    # Fix svd method, if needed
    if (identical(svd.method, "auto"))
      svd.method <- determine.svd.method(L, N - L + 1, neig, ...)

    wmask <- fmask <- weights <- NULL
  } else if (identical(kind, "2d-ssa")) {
    # Coerce input to matrix if necessary
    if (!is.matrix(x))
      x <- as.matrix(x)

    N <- dim(x);

    if (is.null(mask)) {
      mask <- !is.na(x)
    } else {
      mask <- mask & !is.na(x)
    }

    wmask <- .fiface.eval(substitute(wmask),
                          envir = parent.frame(),
                          circle = circle.mask,
                          triangle = triangle.mask)
    if (is.null(wmask)) {
      wmask <- matrix(TRUE, L[1], L[2])
    } else {
      L <- dim(wmask)
    }

    if (is.null(neig))
      neig <- min(50, prod(L), prod(N - L + 1))

    # Fix SVD method.
    if (identical(svd.method, "auto"))
      svd.method <- determine.svd.method(prod(L), prod(N - L + 1), neig, ..., svd.method = "nutrlan")

    fmask <- factor.mask(mask, wmask)

    if (!all(wmask) || !all(fmask)) {
      weights <- field.weights(wmask, fmask)

      ommited <- sum(mask & (weights == 0))
      if (ommited > 0) {
        warning(sprintf("Some field elements were not covered by shaped window. %d elements will be ommited", ommited))
      }

      if (all(weights == 0)) {
        stop("Nothing to decompose: the given field shape is empty")
      }
    } else {
      weights <- NULL
    }

    if (all(wmask))
      wmask <- NULL
    if (all(fmask))
      fmask <- NULL
  } else if (identical(kind, "mssa")) {
    # We assume that we have mts-like object. With series in the columns.
    # Coerce input to series.list object
    # Note that this will correctly remove leading and trailing NA's
    x <- .to.series.list(x, na.rm = TRUE)
    # Grab the inner attributes, if any
    iattr <- lapply(x, attributes)

    N <- sapply(x, length)

    # If L is provided it should be length 1
    if (missing(L)) {
      L <- (min(N) + 1) %/% 2
    } else {
      if (length(L) > 1)
        warning("length of L is > 1, only the first element will be used")
      L <- L[1]
    }

    if (is.null(neig))
      neig <- min(50, L, sum(N - L + 1))

    # Fix SVD method.
    if (identical(svd.method, "auto"))
      svd.method <- determine.svd.method(L, sum(N - L + 1), neig, ...)

    wmask <- NULL
    if (!all(N == max(N))) {
      K <- N - L + 1

      weights <- matrix(0, max(N), length(N))
      fmask <- matrix(FALSE, max(K), length(N))
      for (idx in seq_along(N)) {
        weights[seq_len(N[idx]), idx] <- .hweights.default(N[idx], L)
        fmask[seq_len(K[idx]), idx] <- TRUE
      }
    } else {
      fmask <- weights <- NULL
    }
  } else if (identical(kind, "cssa")) {
    # Sanity check - the input series should be complex
    if (!is.complex(x))
      stop("complex SSA should be performed on complex time series")
    N <- length(x)

    if (is.null(neig))
      neig <- min(50, L, N - L + 1)

    # Fix SVD method.
    if (identical(svd.method, "auto"))
      svd.method <- determine.svd.method(L, N - L + 1, neig, ..., svd.method = "eigen")

    wmask <- fmask <- weights <- NULL
  }
  stopifnot(!is.null(neig))

  # Normalize the kind to be used
  kind <- sub("-", ".", kind, fixed = TRUE)

  # Create information body
  this <- list(length = N,
               window = L,
               call = match.call(),
               kind = kind,
               series = deparse(substitute(x)),
               svd.method = svd.method)

  # Create data storage
  this <- .create.storage(this);

  # Save series
  .set(this, "F", x);

  # Save masks and weights
  .set(this, "wmask", wmask)
  .set(this, "fmask", fmask)
  .set(this, "weights", weights)

  # Save attributes
  .set(this, "Fattr", xattr)
  .set(this, "Fclass", xclass)
  .set(this, "Iattr", iattr)

  # Make this S3 object
  class(this) <- c(paste(kind, svd.method, sep = "."), kind, "ssa");

  # Perform additional init steps, if necessary
  .init(this)

  # Decompose, if necessary
  if (force.decompose)
    this <- decompose(this, neig = neig, ...);

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
  if (desired > min(.traj.dim(x)))
    stop("Cannot decompose that much, desired elementary series index is too huge")

  # Continue decomposition, if necessary
  if (desired > min(nlambda(x), nu(x)))
    decompose(x, ..., neig = min(desired + 1, min(.traj.dim(x))))

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

  for (idx in new) {
    # Do actual reconstruction (depending on method, etc)
    .set.series(x,
                .elseries(x, idx), idx);
  }
}

cleanup <- function(x) {
  .remove(x, ls(.storage(x), pattern = "series:"));
}

.apply.attributes.default <- function(x, F,
                                      fixup = FALSE,
                                      only.new = TRUE, drop = FALSE) {
  a <- (if (drop) NULL else .get(x, "Fattr"))
  cls <- (if (drop) NULL else .get(x, "Fclass"))

  if (fixup) {
     # Try to guess the indices of known time series classes
    if ("ts" %in% cls) {
      tsp <- a$tsp
      return (ts(F,
                 start = if (only.new) tsp[2] + 1/tsp[3] else tsp[1],
                 frequency = tsp[3]))
    }
  } else {
    attributes(F) <- a
  }

  F
}

.group.names <- function(groups) {
  group.names <- names(groups)
  if (is.null(group.names)) group.names <- rep("", length(groups))

  ifelse(group.names != "", group.names, paste0("F", seq_along(groups)))
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

  # Do actual reconstruction. Calculate the residuals on the way
  residuals <- .get(x, "F")
  for (i in seq_along(groups)) {
    group <- groups[[i]];
    new <- setdiff(group, info);
    cached <- intersect(group, info);

    if (length(new)) {
      # Do actual reconstruction (depending on method, etc)
      out[[i]] <- .elseries(x, new);

      # Cache the reconstructed series, if this was requested
      if (cache && length(new) == 1)
        .set.series(x, out[[i]], new);

       # Add pre-cached series
      if (length(cached))
        out[[i]] <- out[[i]] + .get.series(x, cached);
    } else if (length(cached)) {
      out[[i]] <- .get.series(x, cached)
    } else {
      stop("group cannot be empty")
    }

    # Propagate attributes (e.g. dimension for 2d-SSA)
    out[[i]] <- .apply.attributes(x, out[[i]], fixup = FALSE, drop = drop.attributes)
  }

  # Set names
  names(out) <- .group.names(groups)

  # Calculate the residuals
  residuals <- .get(x, "F")
  rgroups <- unique(unlist(groups))
  info <- .get.series.info(x);
  rcached <- intersect(rgroups, info)
  rnew <- setdiff(rgroups, info)
  if (length(rcached))
    residuals <- residuals - .get.series(x, rcached)
  if (length(rnew))
    residuals <- residuals - .elseries(x, rnew)

  # Propagate attributes of residuals
  residuals <- .apply.attributes(x, residuals, fixup = FALSE, drop = drop.attributes)
  F <- .apply.attributes(x, .get(x, "F"), fixup = FALSE, drop = drop.attributes)

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

.elseries.default <- function(x, idx, ...) {
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
      V <- calc.v(x, i);
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
  clp <- (if (length(x$window) > 1) " x " else ", ")
  cat("\nCall:\n", deparse(x$call), "\n\n", sep="");
  cat("Series length:", paste(x$length, collapse = clp));
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
