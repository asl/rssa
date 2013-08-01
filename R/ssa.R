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

determine.svd.method <- function(L, N, neig = NULL, ..., svd.method = "nutrlan") {
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
                neig = NULL,
                ...,
                kind = c("1d-ssa", "2d-ssa", "toeplitz-ssa", "pssa"),
                svd.method = c("auto", "nutrlan", "propack", "svd", "eigen"),
                left.projector = "constant",
                right.projector = "constant",
                force.decompose = TRUE) {
  svd.method <- match.arg(svd.method);
  kind <- match.arg(kind);
  xattr <- attributes(x);

  if (identical(kind, "1d-ssa") || identical(kind, "toeplitz-ssa") || identical(kind, "pssa")) {
    # Coerce input to vector if necessary
    if (!is.vector(x))
      x <- as.vector(x);

    N <- length(x);

    if (is.null(neig))
      neig <- min(50, L, N - L + 1)

    # Fix svd method, if needed
    if (identical(svd.method, "auto"))
      svd.method <- determine.svd.method(L, N, neig, ...)
  } else if (identical(kind, "2d-ssa")) {
    # Coerce input to matrix if necessary
    if (!is.matrix(x))
      x <- as.matrix(x);

    N <- dim(x);

    if (is.null(neig))
      neig <- min(50, prod(L), prod(N - L + 1))

    if (identical(svd.method, "auto"))
      svd.method <- "nutrlan"
  }
  stopifnot(!is.null(neig))

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

  # Compute and store projectors
  if (identical(kind, "pssa")) {
    K <- N - L + 1

    if (is.character(left.projector) || !is.matrix(left.projector))
      left.projector <- orthopoly(left.projector, L)
    if (is.character(right.projector) || !is.matrix(right.projector))
      right.projector <- orthopoly(right.projector, K)


    .set(this, "left.projector", left.projector)
    .set(this, "right.projector", right.projector)
  }

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

  for (idx in new) {
    # Do actual reconstruction (depending on method, etc)
    .set.series(x,
                .elseries(x, idx), idx);
  }
}

cleanup <- function(x) {
  .remove(x, ls(.storage(x), pattern = "series:"));
}

union.groups <- function(groups) {
  triple.kinds <- unique(do.call("c", lapply(groups, names)))
  res <- list()
  for (triple.kind in triple.kinds)
    res[[triple.kind]] <- unique(do.call("c", lapply(groups, function(group) group[[triple.kind]])))

  res
}

# Formula-like interface
.fiface.eval <- function(expr, envir = parent.frame(), ...) {
  env <- as.environment(list(...))
  parent.env(env) <- envir
  env$I <- function(expr) eval(substitute(expr), envir = envir)

  eval(expr, envir = env)
}

.prepare.groups <- function(x, groups.expression, envir = parent.frame(), ...) {
  .left <- function(idx = "all") {
    if (is.character(idx)) {
      proj.triples <- .get(x, "left.proj.triples", allow.null = TRUE)
      idx <- if (!is.null(proj.triples)) seq_len(length(proj.triples$lambda)) else c()
    }

    attr(idx, "triple.kind") <- "left"
    idx
  }

  .right <- function(idx = "all") {
    if (is.character(idx)) {
      proj.triples <- .get(x, "right.proj.triples", allow.null = TRUE)
      idx <- if (!is.null(proj.triples)) seq_len(length(proj.triples$lambda)) else c()
    }

    attr(idx, "triple.kind") <- "right"
    idx
  }

 .svd <- function(idx = "all") {
    if (is.character(idx)) {
      idx <- seq_len(min(nlambda(x), nu(x)))
    }

    attr(idx, "triple.kind") <- "svd"
    idx
  }

  .auto <- function() {
    stop("AutoSSA hasn't been implemented yet")
  }

  groups <- .fiface.eval(groups.expression, envir = envir,
                         left = .left,
                         right = .right,
                         svd = .svd,
                         auto = .auto)

  .get.triple.kind <- function(subgroup, default = "svd") {
    triple.kind.attr <- attr(subgroup, "triple.kind")
    if (!is.null(triple.kind.attr)) triple.kind.attr else default
  }

  groups <- lapply(groups, function(group) {
      if (is.numeric(group)) {
        res <- list(unique(group))
        names(res) <- .get.triple.kind(group)
        res
      } else if (is.list(group)) {
        names.group <- names(group)
        if (is.null(names.group)) names.group <- rep("", length(group))
        triple.kinds <- ifelse(names.group %in% c("left", "right", "svd"),
                               names.group,
                               sapply(group, .get.triple.kind))
        res <- list()
        for (name in unique(triple.kinds)) {
          res[[name]] <- unique(do.call("c", group[triple.kinds == name]))
        }
        res
      }
   })

  # Drop empty subgroups
  groups <- lapply(groups, function(group) group[sapply(group, length) > 0])

  # Assign names
  groups.names <- names(groups)
  if (is.null(groups.names)) groups.names <- rep("", length(groups))
  names(groups) <- ifelse(groups.names != "",
                          groups.names,
                          paste("F", seq_along(groups), sep = ""))

  groups
}

reconstruct.ssa <- function(x, groups, ...,
                            drop.attributes = FALSE, cache = TRUE) {
  out <- list();

  if (missing(groups))
    groups.expression <- expression(as.list(svd("all")))
  else
    groups.expression <- substitute(groups)

  groups <- .prepare.groups(x, groups.expression, envir = parent.frame())

  # Continue decomposition, if necessary
  svd.groups <- unique(do.call("c", lapply(groups, function(group) group$svd)))
  if (length(svd.groups) > 0)
    .maybe.continue(x, groups = svd.groups, ...)

  # Do actual reconstruction. Calculate the residuals on the way
  residuals <- .get(x, "F")
  for (i in seq_along(groups)) {
    group <- groups[[i]];
    out[[i]] <- .elseries(x, group, cache = cache);

    # Propagate attributes (e.g. dimension for 2d-SSA)
    attributes(out[[i]]) <- .get(x, "Fattr");
  }

  # Set names
  names(out) <- names(groups)

  # Calculate the residuals
  rgroups <- union.groups(groups)
  residuals <- .get(x, "F") - .elseries(x, rgroups, cache = cache)

  # Propagate attributes of residuals
  attributes(residuals) <- .get(x, "Fattr");
  F <- .get(x, "F")
  if (!drop.attributes)
    attributes(F) <- .get(x, "Fattr")

  attr(out, "residuals") <- residuals;
  attr(out, "series") <- F;

  # Reconstructed series can be pretty huge...
  class(out) <- paste(c(x$kind, "ssa"), "reconstruction", sep = ".")
  invisible(out);
}

residuals.ssa <- function(object, groups, ..., cache = TRUE) {
  if (missing(groups)) {
    groups.expression <- expression(list(left("all"), right("all"), svd("all")))
  } else {
    groups.expression <- substitute(groups)
  }
  groups <- .prepare.groups(object, groups.expression, envir = parent.frame())

  residuals(reconstruct(object, groups = groups, ..., cache = cache))
}

residuals.ssa.reconstruction <- function(object, ...) {
  attr(object, "residuals")
}

.elseries.default <- function(x, idx, ..., cache = TRUE) {
  if (is.numeric(idx))
    idx <- list(svd = idx)

  res <- numeric(prod(x$length))

  if ("svd" %in% names(idx)) {
    if (max(idx$svd) > nlambda(x))
      stop("Too few eigentriples computed for this decomposition")

    lambda <- .get(x, "lambda")
    U <- .get(x, "U")


    # Grab indices of pre-cached values
    info <- .get.series.info(x);

    for (i in idx$svd) {
      if (i %in% info) {
        elcomponent <- .get.series(x, i)
      } else {
        if (nv(x) >= i) {
          # FIXME: Check, whether we have factor vectors for reconstruction
          # FIXME: Get rid of .get call
          V <- .get(x, "V")[, i]
        } else {
          # No factor vectors available. Calculate them on-fly.
          V <- calc.v(x, i, env = env)
        }

        elcomponent <- lambda[i] * .hankelize.one(x, U = U[, i], V = V)
        if (cache) .set.series(x, elcomponent, i)
      }

      res <- res + elcomponent
    }
  }

  if ("left" %in% names(idx)) {
    triples <- .get(x, "left.proj.triples")
    lambda <- triples$lambda
    U <- triples$U
    V <- triples$V

    for (i in idx$left)
      res <- res + lambda[i] * .hankelize.one(x, U = U[, i], V = V[, i])
  }

  if ("right" %in% names(idx)) {
    triples <- .get(x, "right.proj.triples")
    lambda <- triples$lambda
    U <- triples$U
    V <- triples$V

    for (i in idx$right)
      res <- res + lambda[i] * .hankelize.one(x, U = U[, i], V = V[, i])
  }

  res
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
