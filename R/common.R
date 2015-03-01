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

# Convenient helpers to operate with data storage

.create.storage <- function(x) {
  attr(x, ".env") <- new.env();
  x;
}

.storage <- function(x)
  attr(x, ".env");

.get <- function(x, name, default,
                 allow.null = FALSE) {
  ret <- NULL

  # Make sure default is evaluated only when necessary
  if (.exists(x, name))
    ret <- get(name, envir = .storage(x))
  else if (!allow.null || !missing(default))
    ret <- default

  ret
}

.get.or.create <- function(x, name, default) {
  (if (.exists.non.null(x, name))
   get(name, envir = .storage(x))
  else
   assign(name, default, envir = .storage(x), inherits = FALSE))
}

.set <- function(x, name, value)
  assign(name, value, envir = .storage(x), inherits = FALSE);

.exists <- function(x, name)
  exists(name, envir = .storage(x), inherits = FALSE);

.exists.non.null <- function(x, name) {
  ret <- FALSE

  if (exists(name, envir = .storage(x), inherits = FALSE)) {
    val <- get(name, envir = .storage(x))
    ret <- !is.null(val) && (typeof(val) != "externalptr" || !.is.extptrnull(val))
  }
  ret
}

.remove <- function(x, name)
  rm(list = name, envir = .storage(x), inherits = FALSE);

.deprecate <- function(x, name, instead = NULL) {
  val <- NA
  attr(val, "deprecated") <- TRUE
  attr(val, "instead") <- instead
  .set(x, name, val)
}

.clone <- function(x, copy.storage = TRUE) {
  # Copy the information body
  obj <- x;

  # Make new storage
  obj <- .create.storage(obj);

  # Copy the contents of data storage
  if (copy.storage) {
    clone.env <- .storage(obj);
    x.env <- .storage(x);
    for (field in ls(envir = x.env, all.names = TRUE)) {
      value <- get(field, envir = x.env, inherits = FALSE);
      attr(value, "..cloned") <- NULL;
      assign(field, value, envir = clone.env, inherits = FALSE);
    }
  }

  obj
}

.get.series.info <- function(x) {
  if (.exists(x, "series:info"))
    return (.get(x, "series:info"));

  numeric(0);
}

.append.series.info <- function(x, index) {
  .set(x, "series:info",
           union(.get.series.info(x), index));
}

.set.series <- function(x, F, index) {
  name <- paste("series:", index, sep = "");
  .set(x, name, F);
  .append.series.info(x, index);
  index;
}

.get.series <- function(x, index) {
  name <- paste("series:", index[1], sep = "")
  F <- .get(x, name)
  for (i in index[-1]) {
    name <- paste("series:", i, sep = "")
    F <- F + .get(x, name)
  }
  F
}

.F <- function(x)
  .get(x, "F")

.decomposition <- function(x, field) {
  d <- .get(x, "decomposition", allow.null = TRUE)
  if (missing(field)) d else if (length(field) == 1) d[[field]] else d[field]
}

.set.decomposition <- function(x, ..., kind = "ssa.decomposition") {
  val <- list(...)
  class(val) <- kind

  d <- .get(x, "decomposition", allow.null = TRUE)
  if (is.null(d))
    d <- list()

  .set(x, "decomposition", modifyList(d, val))
}

.U.default <- function(x)
  x$U

.colspan.default <- function(x, idx)
  x$U[, idx, drop = FALSE]

.V.default <- function(x)
  x$V

.rowspan.default <- function(x, idx)
  x$V[, idx, drop = FALSE]

.sigma.default <- function(x)
  x$sigma

.U.ssa <- function(x)
  .decomposition(x, "U")

.colspan.ssa <- function(x, idx)
  .colspan(.decomposition(x), idx)

.V.ssa <- function(x)
  .decomposition(x, "V")

.rowspan.ssa <- function(x, idx)
  calc.v(x, idx)

.sigma.ssa <- function(x)
  .decomposition(x, "sigma")

nspecial.ssa <- function(x)
  return(0)

.is.extptrnull <- function(x)
  .Call("is_extptrnull", x)

.na.omit <- function(x, ...) {
  # Drop initial and final NAs
  good <- which(!is.na(x))
  if (!length(good))
    stop("all times contain an NA")

  omit <- integer()
  n <- length(x)
  st <- min(good)
  if (st > 1)
    omit <- c(omit, 1L:(st-1))
  en <- max(good)
  if (en < n)
    omit <- c(omit, (en+1):n)
  cls <- attr(x, "class")
  if (length(omit)) {
    if ("ts" %in% cls) {
      tm <- time(x)
      xfreq <- frequency(x)
    }
    x <- x[st:en]
    attr(omit, "class") <- "omit"
    attr(x, "na.action") <- omit
    if ("ts" %in% cls)
      tsp(x) <- c(tm[st], tm[en], xfreq)
    if (!is.null(cls)) class(x) <- cls
  }

  x
}

.to.series.list <- function(x, na.rm = TRUE, template = NULL) {
  if (inherits(x, "series.list"))
    return(x)

  # Coerce to list if necessary
  if (!is.list(x)) {
    if (!is.matrix(x))
      x <- as.matrix(x)

    x <- lapply(seq_len(ncol(x)), function(i) x[, i])
  }

  # Note that this will correctly remove leading and trailing NA, but will fail for internal NA's
  if (is.null(template)) {
    # If no template, remove leading and trailing NAs
    NA.fun <- (if (na.rm) .na.omit else identity)
    res <- lapply(x, NA.fun)
  } else {
    # Elsewise omit elements from the same places as in template series (it's necessary for correct conversion to inner format)
    res <- lapply(seq_along(x),
                  function(i) {
                    res <- x[[i]]
                    remove <- attr(res, "na.action") <- attr(template[[i]], "na.action")
                    if (!is.null(remove))
                      res <- res[-remove]

                    res
                  })
  }

  class(res) <- "series.list"

  res
}

.from.series.list <- function(x,
                              pad = c("none", "left", "right"),
                              simplify. = TRUE) {
  pad <- match.arg(pad)

  # First, get rid of "na.omit attribute"
  res <- sapply(x,
                function(x) {
                  removed <- attr(x, "na.action")
                  if (!is.null(removed) && length(removed) > 0) {
                    res <- numeric(length(x) + length(removed))
                    res[removed] <- NA
                    res[-removed] <- x
                    res
                  } else
                  x
                },
                simplify = (if (identical(pad, "none")) simplify. else FALSE))

  # If no padding or further simplification is required, return
  if (identical(pad, "none") || simplify. == FALSE)
    return(res)

  # Pad with NA's
  if (any(sapply(res, is.ts))) {
    do.call(ts.union, res)
  } else {
    ml <- max(sapply(res, length))
    sapply(res,
           function(x) {
             l <- length(x)
             if (identical(pad, "left"))
               c(rep.int(NA, ml - l), x)
             else
               c(x, rep.int(NA, ml - l))
           },
           simplify = simplify.)
  }
}

Ops.series.list <- function(e1, e2 = NULL) {
  unary <- nargs() == 1L
  lclass <- nzchar(.Method[1L])
  rclass <- !unary && (nzchar(.Method[2L]))

  FUN <- get(.Generic, envir = parent.frame(), mode = "function")

  if (lclass && rclass) {
    if (length(e1) != length(e2))
      stop("series list should have equal number of elements")

    res <- lapply(seq_len(length(e1)), function(i) FUN(e1[[i]], e2[[i]]))
  } else if (lclass) {
    res <- (if (unary) lapply(e1, FUN) else lapply(e1, function(x) FUN(x, e2)))
  } else {
    res <- (if (unary) lapply(e2, FUN) else lapply(e2, function(x) FUN(e1, x)))
  }

  class(res) <- "series.list"

  res
}

# Formula-like interface
.fiface.eval <- function(expr, envir = parent.frame(), ...) {
  env <- as.environment(list(...))
  parent.env(env) <- envir
  env$I <- function(expr) eval(substitute(expr), envir = envir)

  eval(expr, envir = env)
}

# Generics

# 'ssa' object
clone <- function(x, ...)
  UseMethod("clone");
reconstruct <- function(x, ...)
  UseMethod("reconstruct");
calc.v <- function(x, ...)
  UseMethod("calc.v");
wnorm <- function(x, ...)
  UseMethod("wnorm")

.hankelize.one <- function(x, ...)
  UseMethod(".hankelize.one")
.hankelize.multi <- function(x, ...)
  UseMethod(".hankelize.multi")
.U <- function(x, ...)
  UseMethod(".U")
.colspan <- function(x, ...)
  UseMethod(".colspan")
.V <- function(x, ...)
  UseMethod(".V")
.rowspan <- function(x, ...)
  UseMethod(".rowspan")
.sigma <- function(x, ...)
  UseMethod(".sigma")
nspecial <- function(x)
  UseMethod("nspecial")
.elseries <- function(x, ...)
  UseMethod(".elseries")
.init <- function(x, ...)
  UseMethod(".init")
.traj.dim <- function(x, ...)
  UseMethod(".traj.dim")
.apply.attributes <- function(x, ...)
  UseMethod(".apply.attributes")

# There is decompose() call in stats package, we need to take control over it
decompose <- function(x, ...) UseMethod("decompose");
decompose.default <- function(x, ...) stats::decompose(x, ...);
formals(decompose.default) <- c(formals(decompose.default), alist(... = ));
