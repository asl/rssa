#   R package for Singular Spectrum Analysis
#   Copyright (c) 2013 Anton Korobeynikov <anton at korobeynikov dot info>
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


.traj.dim.mssa <- function(x) {
  c(x$window, sum(x$length - x$window + 1))
}

.hmat.striped <- function(x, fft.plan) {
  N <- x$length; L <- x$window

  F <- .F(x)
  field <- matrix(0., max(N), length(N))

  weights <- .get(x, "weights")
  if (!is.null(weights))
    mask <- weights > 0
  else
    mask <- matrix(TRUE, max(N), length(N))

  for (idx in seq_along(N)) {
    imask <- which(mask[seq_len(N[idx]), idx])
    field[imask, idx] <- F[[idx]][imask]
  }

  new.hbhmat(field, L = c(L, 1),
             wmask = NULL,
             fmask = .get(x, "fmask"),
             weights = weights)
}

.get.or.create.mhmat <- function(x) {
  .get.or.create(x, "hmat",
                 .hmat.striped(x))
}

decompose.mssa <- function(x,
                           neig = min(50, L, sum(K)),
                           ...,
                           force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1
  stop("Unsupported SVD method for MSSA!")
}

decompose.mssa.svd <- function(x,
                               neig = min(L, sum(K)),
                               ...,
                               force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Create circulant and convert it to ordinary matrix
  h <- as.matrix.hbhmat(.get.or.create.mhmat(x))

  # Do decomposition
  S <- svd(h, nu = neig, nv = neig)

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)

  x
}

decompose.mssa.eigen <- function(x, ...,
                                 neig = min(L, sum(K)),
                                 force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Create circulant and compute XX^T in form of ordinary matrix
  C <- tcrossprod.hbhmat(.get.or.create.mhmat(x))

  # Do decomposition
  S <- eigen(C, symmetric = TRUE)

  # Fix small negative values
  S$values[S$values < 0] <- 0

  # Save results
  .set.decomposition(x,
                     sigma = sqrt(S$values[1:neig]),
                     U = S$vectors[, 1:neig, drop = FALSE])

  x
}

decompose.mssa.propack <- function(x,
                                   neig = min(50, L, sum(K)),
                                   ...,
                                   force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.")

  h <- .get.or.create.mhmat(x)
  S <- propack.svd(h, neig = neig, ...)

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)

  x
}

decompose.mssa.nutrlan <- function(x,
                                   neig = min(50, L, sum(K)),
                                   ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  h <- .get.or.create.mhmat(x)

  sigma <- .sigma(x)
  U <- .U(x)

  S <- trlan.svd(h, neig = neig, ...,
                 lambda = sigma, U = U)

  # Save results
  .set.decomposition(x, sigma = S$d, U = S$u)

  x
}

calc.v.mssa<- function(x, idx, ...) {
  sigma <-.sigma(x)[idx]

  if (any(sigma <= .Machine$double.eps)) {
    sigma[sigma <= .Machine$double.eps] <- Inf
    warning("some sigmas are equal to zero. The corresponding vectors will be zeroed")
  }

  U <- .U(x)[, idx, drop = FALSE]
  h <- .get.or.create.mhmat(x)

  invisible(sapply(1:length(idx),
                   function(i) hbhmatmul(h, U[, i], transposed = TRUE) / sigma[i]))
}

.hankelize.one.mssa <- function(x, U, V) {
  h <- .get.or.create.hbhmat(x)
  storage.mode(U) <- storage.mode(V) <- "double"
  F <- .Call("hbhankelize_one_fft", U, V, h)

  ## FIXME: This is ugly
  N <- x$length; mN <- max(N)
  cidx <- unlist(lapply(seq_along(N), function(idx) seq_len(N[idx]) + mN * (idx - 1)))

  F[cidx]
}

.elseries.mssa <- function(x, idx, ...) {
  if (max(idx) > nsigma(x))
    stop("Too few eigentriples computed for this decomposition")

  N <- x$length
  sigma <- .sigma(x)
  U <- .U(x)
  F <- .F(x)

  res <- numeric(sum(N))
  for (i in idx) {
    if (nv(x) >= i) {
      # FIXME: Check, whether we have factor vectors for reconstruction
      # FIXME: Get rid of .get call
      V <- .V(x)[, i]
    } else {
      # No factor vectors available. Calculate them on-fly.
      V <- calc.v(x, i)
    }

    res <- res + sigma[i] * .hankelize.one(x, U = U[, i], V = V)
  }

  cN <- c(0, cumsum(N))
  sres <- list()
  for (i in seq_along(N)) {
    sres[[i]] <- res[(cN[i]+1):cN[i+1]]
    attr(sres[[i]], "na.action") <- attr(F[[i]], "na.action")
  }
  class(sres) <- "series.list"

  sres
}

.apply.attributes.mssa <- function(x, F,
                                   fixup = FALSE,
                                   only.new = TRUE, drop = FALSE) {
  ia <- (if (drop) NULL else .get(x, "Iattr"))

  # MSSA is a bit different from the default case. We need to convert (if
  # possible) to original object
  stopifnot(inherits(F, "series.list"))

  .restore.attr <- function(F, cls, attr, fixup = FALSE, only.new = TRUE) {
    if (fixup) {
      # Try to guess the indices of known time series classes
      if ("ts" %in% cls) {
        tsp <- attr$tsp
        F <- ts(F,
                start = if (only.new) tsp[2] + 1/tsp[3] else tsp[1],
                frequency = tsp[3])
        attr(F, "na.action") <- NULL
      }
      # It's safe to propagate dimnames in any case
      dimnames(F) <- attr$dimnames
    } else {
      ## Restore attributes

      ## HACK: We need to handle data frames as a special case, because we cannot
      ## convert to them just applying the attributes :(
      if ("data.frame" %in% cls)
        F <- as.data.frame(F)

      attributes(F) <- attr
    }

    F
  }

  # Restore the inner attributes, if any
  for (idx in seq_along(F)) {
    cia <- ia[[idx]]
    if (!is.null(cia))
      F[[idx]] <- .restore.attr(F[[idx]],
                                cia$class, cia,
                                fixup = fixup, only.new = only.new)
  }

  a <- (if (drop) NULL else .get(x, "Fattr"))
  cls <- (if (drop) NULL else .get(x, "Fclass"))

  # Pad with NA's if necessary and optionaly convert to matrix
  F <- .from.series.list(F, pad = "left",
                         simplify. = !("list" %in% cls))

  .restore.attr(F, cls, a,
                fixup = fixup, only.new = only.new)
}

plot.mssa.reconstruction <- function(x,
                                     slice = list(),
                                     ...,
                                     type = c("raw", "cumsum"),
                                     plot.method = c("native", "matplot", "xyplot"),
                                     na.pad = c("left", "right"),
                                     base.series = NULL,
                                     add.original = TRUE,
                                     add.residuals = TRUE) {
  type <- match.arg(type)
  plot.method <- match.arg(plot.method)
  na.pad <- match.arg(na.pad)
  original <- attr(x, "series")
  res <- attr(x, "residuals")

  # Handle base series, if any
  if (!is.null(base.series)) {
    stop("Base series are not supported yet")
    # stopifnot(inherits(base.series, "ssa.reconstruction"))
    # m0 <- matrix(unlist(base.series), ncol = length(base.series))
    # original <- attr(base.series, "series")
  }

  # Nifty defaults
  dots <- list(...)
  if (is.data.frame(original) && identical(plot.method, "native"))
    dots <- .defaults(dots,
                      main = "Reconstructed Series")
  else
    dots <- .defaults(dots,
                      main = "Reconstructed Series",
                      type = "l",
                      ylab = "")

  # Prepare the array with all the data
  m <- lapply(x, .to.series.list, na.rm = FALSE)
  # Fix the slices
  if (is.null(slice$series))
    slice$series <- seq_along(m[[1]])
  if (is.null(slice$component))
    slice$component <- seq_along(m)

  # Slice it
  m <- sapply(m,
              .from.series.list, pad = na.pad, simplify. = "array",
              simplify = "array")
  if (length(dim(m)) == 2)
    dim(m) <- c(NROW(m), 1, NCOL(m))

  m <- m[, slice$series, slice$component, drop = FALSE]

  # Transform the matrix, if necessary
  if (identical(type, "cumsum")) {
    m <- apply(m, c(1, 2), cumsum)
    if (length(dim(m)) == 2)
      dim(m) <- c(1, dim(m))

    m <- aperm(m, c(2, 3, 1))
  }

  m <- matrix(unlist(m), ncol = length(slice$series) * length(slice$component))

  # if (!is.null(base.series)) m <- cbind(m0, m)

  # Fill in the column names
  if (is.list(original)) {
    odimnames <- names(original)
    if (is.null(odimnames))
      odimnames <- paste0("F", seq_along(original))
  } else {
    odimnames <- colnames(original, do.NULL = FALSE, prefix = "F")
  }
  rdimnames <- odimnames[slice$series]

  # Fix the attributes
  a <- attributes(original)

  # Get rid of dim and dimnames attribute
  a$dim <- NULL
  a$dimnames <- NULL
  a$names <- NULL

  # Merge the attributes in
  # HACK: We need to handle data frames as a special case, because we cannot
  # convert to them just applying the attributes :(
  if (is.data.frame(original))
    m <- as.data.frame(m)
  attributes(m) <- append(attributes(m), a)

  cnames <- names(x)[slice$component]
  mnames <- sapply(seq_along(cnames), function(x) paste(rdimnames, if (cnames[x] != "") cnames[x] else x))
  if (add.original) {
    original <- .from.series.list(.to.series.list(original), pad = na.pad, simplify. = "array")
    if (is.ts(original))
      m <- ts.union(original[, slice$series], m)
    else
      m <- cbind(original[, slice$series], m)
    mnames <- c(paste("Original", rdimnames), mnames)
  }
  if (add.residuals) {
    res <- .from.series.list(.to.series.list(res), pad = na.pad, simplify. = "array")
    if (is.ts(original))
      m <- ts.union(m, res[, slice$series])
    else
      m <- cbind(m, res[, slice$series])
    mnames <- c(mnames, paste("Residuals", rdimnames))
  }
  colnames(m) <- mnames

  # Plot'em'all!
  if (identical(plot.method, "xyplot"))
    do.call(xyplot, c(list(m), dots))
  else if (identical(plot.method, "matplot") || !is.object(m))
    do.call(matplot, c(list(x = m), dots))
  else if (identical(plot.method, "native"))
    do.call(plot, c(list(m), dots))
  else
    stop("Unknown plot method")
}

xyplot.matrix <- function(x, ..., outer = TRUE) {
  dots <- list(...)
  x <- cbind(1:nrow(x), as.data.frame(x))
  stopifnot(nrow(x) > 1)
  nms <- sprintf("`%s`", colnames(x))
  form <- sprintf("%s ~ %s", paste(nms[-1], collapse = "+"), nms[1])
  # Provide convenient defaults
  dots <- .defaults(dots,
                    as.table = TRUE)
  do.call(xyplot, c(list(as.formula(form), data = x), dots, list(outer = outer)))
}

.plot.ssa.vectors.mssa <- function(x, ...,
                                   what = c("eigen", "factor"),
                                   plot.contrib, idx, plot.type = "l") {
  what <- match.arg(what)

  # Eigenvectors are same as for 1D-SSA case
  if (identical(what, "eigen")) {
    ## Do the call with the same set of arguments
    mplot <- match.call(expand.dots = TRUE)
    mplot[[1L]] <- as.name(".plot.ssa.vectors.1d.ssa")
    mplot[[2L]] <- x
    eval(mplot, parent.frame())
  } else {
    dots <- list(...)

    # Factor vectors are special...
    # We actually create a reconstruction object with all stuff there..
    if (max(idx) > nsigma(x))
      stop("Too few eigentriples computed for this decomposition")

    N <- x$length
    L <- x$window
    K <- N - L + 1
    F <- .F(x)

    cK <- c(0, cumsum(K))
    res <- list()
    for (i in idx) {
      if (nv(x) >= i) {
        # FIXME: Check, whether we have factor vectors for reconstruction
        # FIXME: Get rid of .get call
        V <- .V(x)[, i]
      } else {
        # No factor vectors available. Calculate them on-fly.
        V <- calc.v(x, i)
      }

      sres <- list()
      for (j in seq_along(K)) {
        sres[[j]] <- V[(cK[j]+1):cK[j+1]]
        attr(sres[[j]], "na.action") <- attr(F[[j]], "na.action")
      }
      class(sres) <- "series.list"
      res[[i]] <- .apply.attributes(x, sres,
                                    fixup = TRUE, only.new = FALSE, drop = FALSE)
    }

    # Build pseudo-original series
    oF <- lapply(seq_along(F),
                 function(idx) {
                   x <- F[[idx]]
                   removed <- attr(x, "na.action")
                   if (!is.null(removed) && length(removed) > 0) {
                     m <- max(setdiff(seq_along(length(x) + length(removed)), removed))
                     to.fix <- removed > m
                     removed[to.fix] <- removed[to.fix] - L + 1
                   }
                   res <- x[1:K[idx]]
                   attr(res, "na.action") <- removed
                   res
                 })
    class(oF) <- "series.list"
    attr(res, "series") <- .apply.attributes(x, oF,
                                             fixup = TRUE, only.new = FALSE, drop = FALSE)

    names(res) <- if (!plot.contrib) idx else paste(idx, " (", .contribution(x, idx, ...), "%)", sep = "")

    # Provide convenient defaults
    dots <- .defaults(dots,
                      xlab = "",
                      ylab =  "",
                      main = "Factor vectors",
                      as.table = TRUE,
                      scales = list(draw = FALSE, relation = "free"),
                      aspect = 1,
                      symmetric = TRUE,
                      ref = TRUE)

    do.call(plot.mssa.reconstruction, c(list(res,
                                             plot.method = "xyplot",
                                             add.original = FALSE,
                                             add.residuals = FALSE),
                                        dots))
  }
}
