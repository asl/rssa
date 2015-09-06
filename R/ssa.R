#   R package for Singular Spectrum Analysis
#   Copyright (c) 2008-2015 Anton Korobeynikov <asl@math.spbu.ru>
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

.default.neig <- function(x, ...)
  UseMethod(".default.neig")

.default.neig.ssa <- function(x, ...) {
  tjdim <- .traj.dim(x)

  min(50, tjdim)
}

.determine.svd.method <- function(x, kind, neig = NULL,
                                  ...,
                                  svd.method = (if (identical(kind, "cssa")) "eigen" else "nutrlan")) {
  tjdim <- .traj.dim(x)
  L <- tjdim[1]; K <- tjdim[2]

  truncated <- (identical(svd.method, "nutrlan") || identical(svd.method, "propack"))

  if (is.null(neig))
    neig <- .default.neig(x, ...)

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
                mask = NULL, wmask = NULL,
                column.projector = "none", row.projector = "none",
                column.oblique = "identical", row.oblique = "identical",
                ...,
                kind = c("1d-ssa", "2d-ssa", "nd-ssa", "toeplitz-ssa", "mssa", "cssa"),
                circular = FALSE,
                svd.method = c("auto", "nutrlan", "propack", "svd", "eigen"),
                force.decompose = TRUE) {
  svd.method <- match.arg(svd.method)

  # Squeeze the attributes
  xattr <- attributes(x)
  iattr <- NULL
  # Grab class separately. This way we will capture the inherit class as well
  xclass <- class(x)

  call <- match.call(); cargs <- as.list(call)[-1]
  ## wmask is special and will be treated separately later
  cargs$wmask <- NULL
  ecall <- do.call("call", c("ssa", lapply(cargs, eval, parent.frame())))

  ## Provide some sane defaults, e.g. complex inputs should default to cssa
  if (missing(kind)) {
    if (is.complex(x))
      kind <- "cssa"
    else if (inherits(x, "mts") || inherits(x, "data.frame") || inherits(x, "list") || inherits(x, "series.list"))
      kind <- "mssa"
    else if (is.matrix(x))
      kind <- "2d-ssa"
    else if (is.array(x))
      kind <- "nd-ssa"
    else
      kind <- "1d-ssa"
  }
  kind <- match.arg(kind)

  # Do the fixups depending on the kind of SSA.
  if (identical(kind, "1d-ssa") || identical(kind, "toeplitz-ssa")) {
    ## Nothing special here (yet!)
  } else if (identical(kind, "2d-ssa") || identical(kind, "nd-ssa")) {
    # 2d-SSA is just a special case of nd-ssa
    if (length(dim(x)) == 2)
      kind <- c("2d-ssa", "nd-ssa")
    else
      kind <- "nd-ssa"
  } else if (identical(kind, "mssa")) {
    ## Nothing special here (yet!)
  } else if (identical(kind, "cssa")) {
    ## Nothing special here (yet!)
  } else {
    N <- -1;
    fmask <- NULL
    stop("invalid SSA kind")
  }

  if (!identical(column.projector, "none") || !identical(row.projector, "none")) {
    # Add `pssa` class if appropriate implementation exists

    if (!any(match(kind, c("1d-ssa", "2d-ssa", "nd-ssa")))) {
      stop("SSA with projection is not implemented for such SSA kind yet")
    }

    kind <- c("pssa", paste("pssa", kind, sep = "-"), kind)
  }

  if (!identical(column.oblique, "identical") || !identical(row.oblique, "identical")) {
    # Add `wossa` class if appropriate implementation exists

    if (!any(match(kind, c("1d-ssa", "2d-ssa", "nd-ssa"))) || "pssa" %in% kind) {
      stop("SSA with weights is not implemented for such SSA kind yet")
    }

    kind <- c("wossa", paste("wossa", kind, sep = "-"), kind)
  }

  # Normalize the kind to be used
  kind <- gsub("-", ".", kind, fixed = TRUE)

  # Create information body
  this <- list(call = call, ecall = ecall,
               kind = kind,
               svd.method = svd.method)

  # Create data storage
  this <- .create.storage(this)

  # Save the names of the essential fields
  this$fields <- c("F",
                   "wmask", "fmask", "weights", "circular",
                   "Fattr", "Fclass", "Iattr",
                   "column.projector", "row.projector",
                   "column.oblique", "row.oblique")

  # Make this S3 object
  class(this) <- c(kind, "ssa")

  ## Perform additional init steps, if necessary. We cannot simply eval .init in
  ## the current environment because we're using S3 dispatch at the same
  ## time... UseMethod uses NSE.
  ## NOTE: This will modify the *current* environment (local vars of the function)
  parent.env <- parent.frame()
  eval(.init.fragment(this))

  # Save attributes
  .set(this, "Fattr", xattr)
  .set(this, "Fclass", xclass)
  .set(this, "Iattr", iattr)

  # Deprecated stuff
  .deprecate(this, "lambda", "sigma")

  ## Window and series length should be ready by this moment
  this$length <- N
  this$window <- L

  ## Save series
  .set(this, "F", x)

  ## Save masks, weights and topology
  .set(this, "wmask", wmask)
  .set(this, "fmask", fmask)
  .set(this, "weights", weights)
  .set(this, "circular", circular)

  ## Store projectors
  .set(this, "column.projector", column.projector)
  .set(this, "row.projector", row.projector)

  ## Store oblique matrices
  .set(this, "column.oblique", column.oblique)
  .set(this, "row.oblique", row.oblique)

  # Determine the desired number of eigentriples, if necessary
  if (is.null(neig))
    neig <- .default.neig(this, ...)

  # Fix SVD method
  if (identical(svd.method, "auto"))
    svd.method <- .determine.svd.method(this, kind = kind, neig = neig, ...)

  this$svd.method <- svd.method
  
  # Decompose, if necessary
  if (force.decompose) {
    if (!is.null(weights) && all(weights == 0))
      stop("Nothing to decompose: the given field shape is empty")

    this <- decompose(this, neig = neig, ...);
  }

  this;
}

.init.fragment.default <- function(x, ...) {
  # Do nothing
  x
}

.maybe.continue <- function(x, groups, ...) {
  L <- x$window
  K <- x$length - x$window + 1

  # Determine the upper bound of desired eigentriples
  desired <- max(unlist(groups), -Inf)

  # Sanity check
  if (desired > min(.traj.dim(x)))
    stop("Cannot decompose that much, desired elementary series index is too huge")

  # Continue decomposition, if necessary
  if (desired > min(nsigma(x), nu(x)))
    decompose(x, ..., neig = min(desired + 1 - nspecial(x), min(.traj.dim(x)))) #TODO: Fix it for PSSA

  desired
}

precache <- function(x, n, ...) {
  if (missing(n)) {
    warning("Amount of sub-series missed, precaching EVERYTHING",
            immediate. = TRUE);
    n <- nsigma(x)
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
                                      fixup = FALSE, only.new = TRUE,
                                      reverse = FALSE,
                                      drop = FALSE) {
  a <- (if (drop) NULL else .get(x, "Fattr"))
  cls <- (if (drop) NULL else .get(x, "Fclass"))

  if (fixup) {
     # Try to guess the indices of known time series classes
    if ("ts" %in% cls) {
      tsp <- a$tsp
      return (if (!reverse)
                ts(F,
                   start = if (only.new) tsp[2] + 1/tsp[3] else tsp[1],
                   frequency = tsp[3])
              else
                ts(F,
                   end = if (only.new) tsp[1] - 1/tsp[3] else tsp[2],
                   frequency = tsp[3])

              )
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
    groups <- as.list(1:min(nsigma(x), nu(x)));

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = groups, ...)

  # Grab indices of pre-cached values
  info <- .get.series.info(x);

  # Do actual reconstruction. Calculate the residuals on the way
  residuals <- .F(x)
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
      out[[i]] <- 0. * .F(x)
      if (!is.null(x$weights))
        out[[i]][x$weights == 0] <- NA
    }

    # Propagate attributes (e.g. dimension for 2d-SSA)
    out[[i]] <- .apply.attributes(x, out[[i]], fixup = FALSE, drop = drop.attributes)
  }

  # Set names
  names(out) <- .group.names(groups)

  # Calculate the residuals
  residuals <- .F(x)
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
  F <- .apply.attributes(x, .F(x), fixup = FALSE, drop = drop.attributes)

  attr(out, "residuals") <- residuals;
  attr(out, "series") <- F;

  # Reconstructed series can be pretty huge...
  class(out) <- paste(c(x$kind, "ssa"), "reconstruction", sep = ".")
  invisible(out);
}

residuals.ssa <- function(object, groups, ..., cache = TRUE) {
  groups <- list(if (missing(groups)) 1:min(nsigma(object), nu(object)) else unlist(groups))

  residuals(reconstruct(object, groups = groups, ..., cache = cache))
}

residuals.ssa.reconstruction <- function(object, ...) {
  attr(object, "residuals")
}

.elseries.default <- function(x, idx, ...) {
  if (max(idx) > nsigma(x))
    stop("Too few eigentriples computed for this decomposition")

  dec <- .decomposition(x)
  sigma <- .sigma(dec)
  U <- .U(dec)

  res <- numeric(prod(x$length));
  for (i in idx) {
    if (nv(x) >= i) {
      # FIXME: Check, whether we have factor vectors for reconstruction
      # FIXME: Get rid of .get call
      V <- .V(x)[, i];
    } else {
      # No factor vectors available. Calculate them on-fly.
      V <- calc.v(x, i);
    }

    res <- res + sigma[i] * .hankelize.one(x, U = U[, i], V = V);
  }

  res;
}

nu <- function(x) {
  res <- ncol(.U(x))
  ifelse(is.null(res), 0, res)
}

nv <- function(x) {
  res <- ncol(.V(x))
  ifelse(is.null(res), 0, res)
}

nlambda <- function(x) {
  warning("this function is deprecated, use `nsigma' instead")
  nsigma(x)
}

nsigma <- function(x) {
  length(.sigma(x))
}

is.shaped <- function(x) {
  ## Easy case: non-null masks in case of non-mssa
  if ((!is.null(x$wmask) || !is.null(x$fmask) || !is.null(x$weights)) && !inherits(x, "mssa"))
    return (TRUE)

  ## In case of mssa, check whether we have any zero meaningfull weights
  if (inherits(x, "mssa") && any(.hweights(x) == 0))
    return (TRUE)

  return (FALSE)
}

clone.ssa <- function(x, copy.storage = TRUE, copy.cache = TRUE, ...) {
  obj <- .clone(x, copy.storage = copy.storage)

  # We need to copy the "essential" fields
  if (copy.storage == FALSE)
    for (field in x$fields)
      .set(obj, field, .get(x, field))

  if (copy.cache == FALSE)
    cleanup(obj)

  obj;
}

'$.ssa' <- function(x, name) {
  # First, check the fields of the object itself
  if (ind <- charmatch(name, names(x), nomatch = 0))
    return (x[[ind]])

  # Now, check the fields of the storage
  res <- .get(x, name, allow.null = TRUE)
  if (!is.null(res)) {
     # Check for deprecation
    if (isTRUE(attr(res, "deprecated"))) {
      msg <- paste("the field `", name, "' is deprecated", sep = "")
      instead <- attr(res, "instead")

      # If no substitution is available, just stop here
      if (is.null(instead))
        stop(msg)

      # Otherwise, warn and fallback to new name
      warning(paste(msg, ". use `", instead, "' instead.", sep = ""))
      res <- Recall(x, instead)
    }
    return (res)
  }

  # Final special case: the fields of the decomposition
  .decomposition(x)[[name]]
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
  cat("\nSpecial triples: ", nspecial(x));
  cat("\n\nComputed:\n");
  cat("Eigenvalues:", nsigma(x));
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
