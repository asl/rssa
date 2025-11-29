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

.traj.dim.cmssa <- function(x) {
  Ldim <- sum(x$wmask)
  if (Ldim == 0)
    Ldim <- x$window

  Kdim <- sum(x$fmask)
  if (Kdim == 0)
    Kdim <- sum(x$length - x$window + 1)

  c(Ldim, Kdim)
}

.chmat.striped <- function(x, fft.plan, re = TRUE) {
  N <- x$length; L <- x$window

  if (re)
    F <- lapply(.F(x), Re)
  else
    F <- lapply(.F(x), Im)

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

.get.or.create.cmhmat <- function(x, re = TRUE) {
  .get.or.create(x, paste0("hmat_", if (re) "r" else "i"),
                 .chmat.striped(x, fft.plan = .get.or.create.cfft.plan(x), re = re))
}

.get.or.create.trajmat.cmssa <- .get.or.create.cmhmat

# TODO: add diagonal averaging similarly to the hankel function.
# TODO: efficient implementation similarly to new.hbhmat.
mhankel <- function(X, L) {
  if (!inherits(X, "series.list"))
    X <- .to.series.list(X)

  Reduce(cbind, lapply(X, hankel, L = L))
}

decompose.cmssa <- function(x,
                            neig = NULL,
                            ...,
                            force.continue = FALSE) {
  # Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  if (is.null(neig))
    neig <- .default.neig(x, ...)

  if (identical(x$svd.method, "svd")) {
    S <- svd(mhankel(.F(x), L = x$window), nu = neig, nv = neig)
    .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else if (identical(x$svd.method, "eigen")) {
    h <- mhankel(.F(x), L = x$window)

    ## FIXME: Build the complex L-covariance matrix properly
    S <- eigen(tcrossprod(h, Conj(h)), symmetric = TRUE)

    ## Fix small negative values
    S$values[S$values < 0] <- 0

    ## Save results
    .set.decomposition(x,
                       sigma = sqrt(S$values[seq_len(neig)]),
                       U = S$vectors[, seq_len(neig), drop = FALSE])
  } else if (identical(x$svd.method, "primme")) {
    if (!requireNamespace("PRIMME", quietly = TRUE))
      stop("PRIMME package is required for SVD method `primme'")

    R <- .get.or.create.cmhmat(x, re = TRUE)
    I <- .get.or.create.cmhmat(x, re = FALSE)

    matmul <- function(x, y) {
      if (is.matrix(y))
        apply(y, 2, ematmul, emat = x, transposed = FALSE)
      else
        ematmul(x, y, transposed = FALSE)
    }

    tmatmul <- function(x, y) {
      if (is.matrix(y))
        apply(y, 2, ematmul, emat = x, transposed = TRUE)
      else
        ematmul(x, y, transposed = TRUE)
    }

    A <- function(x, trans) {
      rX <- Re(x)
      iX <- Im(x)
      if (identical(trans, "c")) {
        (tmatmul(R, rX) + tmatmul(I, iX)) +
          1i * (tmatmul(R, iX) - tmatmul(I, rX))
      } else {
        (matmul(R, rX) - matmul(I, iX)) +
          1i * (matmul(R, iX) + matmul(I, rX))
      }
    }

    S <- PRIMME::svds(A, NSvals = neig, m = nrow(R), n = ncol(R), isreal = FALSE, ...)
    .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else if (identical(x$svd.method, "irlba")) {
    if (!requireNamespace("irlba", quietly = TRUE))
      stop("irlba package is required for SVD method `irlba'")

    h <- mhankel(.F(x), L = x$window)
    S <- irlba::irlba(h, nv = neig, ...)
    .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else if (identical(x$svd.method, "rsvd")) {
    if (!requireNamespace("irlba", quietly = TRUE))
      stop("irlba package is required for SVD method `rsvd'")

    h <- mhankel(.F(x), L = x$window)
    S <- irlba::svdr(h, k = neig, ...)
    .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else
    stop("unsupported SVD method")

  x
}

calc.v.cmssa <- function(x, idx, ...) {
  nV <- nv(x)

  V <- matrix(NA_complex_, .traj.dim(x)[2], length(idx))
  idx.old <- idx[idx <= nV]
  idx.new <- idx[idx > nV]

  if (length(idx.old) > 0) {
    V[, idx <= nV] <- .V(x)[, idx.old]
  }

  if (length(idx.new) > 0) {
    sigma <- .sigma(x)[idx.new]

    if (any(sigma <= .Machine$double.eps)) {
      sigma[sigma <= .Machine$double.eps] <- Inf
      warning("some sigmas are equal to zero. The corresponding vectors will be zeroed")
    }

    U <- .U(x)[, idx.new, drop = FALSE]


    tmatmul <- function(x, y) {
      if (is.matrix(y))
        apply(y, 2, ematmul, emat = x, transposed = TRUE)
      else
        ematmul(x, y, transposed = TRUE)
    }

    Rh <- .get.or.create.cmhmat(x, re = TRUE)
    Ih <- .get.or.create.cmhmat(x, re = FALSE)

    hcrossprod <- function(x) {
      rx <- Re(x)
      ix <- Im(x)
      (tmatmul(Rh, rx) + tmatmul(Ih, ix)) +
        1i * (tmatmul(Rh, ix) - tmatmul(Ih, rx))
    }

    V[, idx > nV] <- apply(U, 2, hcrossprod) / sigma

  }

  invisible(V)
}

.hankelize.one.cmssa <- function(x, U, V, fft.plan = .get.or.create.cfft.plan(x)) {
  hankelize.part <- function(x, U, V) {
    h <- .get.or.create.cmhmat(x)
    storage.mode(U) <- storage.mode(V) <- "double"
    F <- .Call("hbhankelize_one_fft", U, V, h@.xData)

    ## FIXME: This is ugly
    N <- x$length
    mN <- max(N)
    cidx <- unlist(lapply(seq_along(N), function(idx)
      seq_len(N[idx]) + mN * (idx - 1)))

    F[cidx]
  }

  R1 <- hankelize.part(x, Re(U), Re(V))
  R2 <- hankelize.part(x, Im(U), Im(V))
  I1 <- hankelize.part(x, Re(U), Im(V))
  I2 <- hankelize.part(x, Im(U), Re(V))

  (R1 + R2) + 1i * (-I1 + I2)
}

.elseries.cmssa <- function(x, idx, ...) {
  if (max(idx) > nsigma(x))
    stop("Too few eigentriples computed for this decomposition")

  N <- x$length
  sigma <- .sigma(x)
  U <- .U(x)
  F <- .F(x)

  res <- complex(sum(N))
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

plot.cmssa.reconstruction <- function(x, ...) {
  rx = lapply(x, Re)
  ix = lapply(x, Im)

  a <- attributes(x)
  ar <- a
  ar$series <- Re(a$series)
  ar$residuals <- Re(a$residuals)
  attributes(rx) <- ar
  ai <- a
  ai$series <- Im(a$series)
  ai$residuals <- Im(a$residuals)
  attributes(ix) <- ai

  dots <- list(...)
  if (is.null(dots$main))
    title <- "Reconstructed Series"
  else
    title <- dots$main
  dots$main <- NULL

  do.call(plot.mssa.reconstruction, c(list(x = rx, main = paste(title, "(Real part)")), dots))
  do.call(plot.mssa.reconstruction, c(list(x = ix, main = paste(title, "(Imaginary part)")), dots))
}

.init.fragment.cmssa <- function(this)
  expression({
    if (any(circular))
      stop("Circular variant of complex multichannel SSA isn't implemented yet")

    # We assume that we have mts-like object. With series in the columns.
    # Coerce input to series.list object
    # Note that this will correctly remove leading and trailing NA's
    x <- .to.series.list(x, na.rm = TRUE)

    # Sanity check - the input series should be complex
    if (!all(sapply(x, is.complex)))
      stop("complex SSA should be performed on complex time series")

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

    wmask <- NULL
    if (!all(N == max(N)) || any(sapply(x, anyNA))) {
      K <- N - L + 1

      weights <- matrix(0, max(N), length(N))
      fmask <- matrix(FALSE, max(K), length(N))
      wmask <- rep(TRUE, L)
      for (idx in seq_along(N)) {
        mask <- !is.na(x[[idx]])
        fmask[seq_len(K[idx]), idx] <- .factor.mask.1d(mask, wmask)
        weights[seq_len(N[idx]), idx] <- .field.weights.1d(wmask, fmask[seq_len(K[idx]), idx])
      }
    } else {
      fmask <- weights <- NULL
    }

    column.projector <- row.projector <- NULL
  })
