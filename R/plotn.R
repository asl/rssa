#   R package for Singular Spectrum Analysis
#   Copyright (c) 2014 Alex Shlemov <shlemovalex@gmail.com>
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


.do.slice.array <- function(x, slice) {
  dim.names <- names(slice)
  if (is.null(dim.names) || any(dim.names == "")) {
    stop("All `slice' elements should be named")
    # TODO Maybe assign unnamed arguments by order like in function call?
  }

  parse.dim.index <- function(name) {
    if (grepl("^d?\\d+$", name)) {
      dim.index <- as.numeric(sub("d", "", name, fixed = TRUE))
    } else {
      dim.index <- switch(name,
                          x =, i = 1,
                          y =, j = 2,
                          z =, k = 3,
                          t = 4)
      if (is.null(dim.index)) {
        stop(sprintf("%s is not proper dimension name", name))
      }
    }

    dim.index
  }

  dim.indices <- sapply(dim.names, parse.dim.index)

  stopifnot(length(unique(dim.indices)) == length(dim.indices))

  # Coerce argument to array
  if (!is.array(x)) x <- as.array(x)
  rank <- length(dim(x))

  stopifnot(max(dim.indices) <= rank)

  args <- as.list(rep(TRUE, rank))
  for (i in seq_along(dim.indices)) {
    args[[dim.indices[i]]] <- slice[[i]]
  }

  res <- do.call("[", c(list(x, drop = TRUE), args))

  res
}

.do.slice.ssa.reconstruction <- function(x, slice) {
  for (i in seq_along(x)) {
    x[[i]] <- .do.slice.array(x[[i]], slice)
  }

  attr(x, "series") <- .do.slice.array(attr(x, "series"), slice)
  attr(x, "residuals") <- .do.slice.array(attr(x, "residuals"), slice)

  x
}

plot.nd.ssa.reconstruction <- function(x, slice, ...) {
  x <- .do.slice.ssa.reconstruction(x, slice)

  series <- attr(x, "series")
  rank <- length(dim(as.array(series)))

  if (rank == 1)
    plot.1d.ssa.reconstruction(x, ...)
  else if (rank == 2)
    plot.2d.ssa.reconstruction(x, ...)
  else
    stop("Cannot display array of rank higher than 2. Use `slice' argument")
}

.plot.ssa.vectors.nd.ssa <- function(x, slice, ...,
                                    what = c("eigen", "factor"),
                                    plot.contrib = FALSE,
                                    idx) {
  what <- match.arg(what)

  dots <- list(...)

  if (max(idx) > nsigma(x))
    stop("Too few eigentriples computed for this decomposition")

  N <- x$length
  L <- x$window
  K <- ifelse(x$circular, N, N - L + 1)

  if (identical(what, "eigen")) {
    mask <- x$wmask
    dimension <- L
    vmatrix <- .U(x)[, idx, drop = FALSE]
  } else if (identical(what, "factor")) {
    mask <- x$fmask
    dimension <- K
    vmatrix <- matrix(NA_real_, nrow = .traj.dim(x)[2], ncol = length(idx))

    vmatrix[, idx <= nv(x)] <- .V(x)[, idx[idx <= nv(x)]]
    if (any(idx > nv(x))) {
      # Some factor vectors are not available. Calculate them on-fly.
      vmatrix[, idx > nv(x)] <- calc.v(x, idx[idx > nv(x)])
    }
  }

  if (is.null(mask)) mask <- array(TRUE, dim = dimension)

  # We actually create a reconstruction object with all stuff there..
  res <- lapply(seq_len(ncol(vmatrix)),
                function(i) {
                  vec <- array(NA_real_, dim = dimension)
                  vec[mask] <- vmatrix[, i]

                  vec
                })

  # Make and set fake initial series and residuals
  fakeseries <- array(0, dim = dimension)
  fakeresiduals <- array(NA_real_, dim = dimension); fakeresiduals[mask] <- 0
  attr(res, "series") <- fakeseries
  attr(res, "residuals") <- fakeresiduals

  names(res) <- if (!plot.contrib) idx else paste(idx, " (", .contribution(x, idx, ...), "%)", sep = "")


  # Provide convenient defaults
  dots <- .defaults(dots,
                    xlab = "",
                    ylab =  "",
                    main = if (identical(what, "eigen")) "Eigenvectors" else "Factor vectors",
                    as.table = TRUE,
                    scales = list(draw = FALSE, relation = "free"),
                    aspect = 1,
                    plot.type = "l",
                    symmetric = TRUE,
                    ref = TRUE)

  do.call(plot.nd.ssa.reconstruction, c(list(res,
                                            slice = slice,
                                            add.original = FALSE,
                                            add.residuals = FALSE),
                                       dots))
}
