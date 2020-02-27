#   R package for Singular Spectrum Analysis
#   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
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

#   Routines for normal hankel SSA

hankel <- function(X, L) {
  if (is.matrix(X) && nargs() == 1) {
     L <- nrow(X); K <- ncol(X); N <- K + L - 1
     left  <- c(1:L, L*(2:K))
     right <- c(1+L*(0:(K-1)), ((K-1)*L+2):(K*L))
     v <- sapply(1:N, function(i) mean(X[seq.int(left[i], right[i], by = L-1)]))
     return (v)
  }

  # Coerce output to vector, if necessary
  if (!is.vector(X))
    X <- as.vector(X)
  N <- length(X)
  if (missing(L))
    L <- (N + 1) %/% 2
  K <- N - L + 1
  outer(1:L, 1:K, function(x,y) X[x+y-1])
}

.convolve1 <- function(x, y, conj = TRUE, type = "circular") {
  if (length(type) > 1) {
    warning("Incorrect argument length: length(type) > 1, the first value will be used")
    type <- type[1]
  }

  type <- match.arg(type, choices = c("circular", "open", "filter"))

  convolution.size <- function(length.x, length.y, type) {
    switch(type,
           circular = list(input = length.x, output = length.x),
           open = list(input = length.x + length.y - 1, output = length.x + length.y - 1),
           filter = list(input = length.x, output = length.x - length.y + 1))
  }

  cs <- convolution.size(length(x), length(y), type)

  input.dim <- cs$input
  output.dim <- cs$output

  X <- Y <- rep(0, input.dim)

  X[seq_along(x)] <- x
  Y[seq_along(y)] <- y

  tmp <- Re(fft(fft(X) * if (conj) Conj(fft(Y)) else fft(Y), inverse = TRUE)) / input.dim
  tmp[seq_len(output.dim)]
}

.factor.mask.1d <- function(field.mask, window.mask, circular = FALSE) {
  field.mask[] <- as.numeric(field.mask)
  window.mask[] <- as.numeric(window.mask)
  tmp <- .convolve1(field.mask, window.mask, conj = TRUE,
                   type = ifelse(circular, "circular", "filter"))

  abs(tmp - sum(window.mask)) < 0.5 # ==0, but not exact in case of numeric error
}

.field.weights.1d <- function(window.mask, factor.mask, circular = FALSE) {
  window.mask[] <- as.numeric(window.mask)
  factor.mask[] <- as.numeric(factor.mask)
  res <- .convolve1(factor.mask, window.mask, conj = FALSE,
                   type = ifelse(circular, "circular", "open"))
  res[] <- as.integer(round(res))

  res
}

.get.or.create.fft.plan <- function(x) {
  .get.or.create(x, "fft.plan", fft.plan.1d(x$length, L = x$window, circular = x$circular,
                                            wmask = x$wmask, fmask = x$fmask, weights = x$weights))
}

.get.or.create.hmat <- function(x) {
  .get.or.create(x, "hmat",
                 new.hmat(.F(x), L = x$window, circular = x$circular,
                          wmask = x$wmask, fmask = x$fmask, weights = x$weights,
                          fft.plan = .get.or.create.fft.plan(x)))
}

.get.or.create.trajmat.1d.ssa <- .get.or.create.hmat
.get.or.create.trajmat.toeplitz.ssa <- .get.or.create.hmat

.hankelize.one.default <- function(U, V, fft.plan = NULL) {
  L <- length(U); K <- length(V); N = K + L - 1
  fft.plan <- (if (is.null(fft.plan)) fft.plan.1d(N) else fft.plan)
  storage.mode(U) <- storage.mode(V) <- "double"
  .Call("hankelize_one_fft", U, V, fft.plan)
}

.hankelize.one.1d.ssa <- function(x, U, V, fft.plan = NULL) {
  fft.plan <- (if (is.null(fft.plan)) .get.or.create.fft.plan(x) else fft.plan)
  storage.mode(U) <- storage.mode(V) <- "double"
  .Call("hankelize_one_fft", U, V, fft.plan)
}

.hankelize.multi.default <- function(U, V, fft.plan) {
  stopifnot(is.numeric(V))
  storage.mode(U) <- storage.mode(V) <- "double"
  .Call("hankelize_multi_fft", U, V, fft.plan)
}

fft.plan.1d <- function(N, L, circular = FALSE,
                        wmask = NULL, fmask = NULL, weights = NULL) {
  storage.mode(N) <- "integer"

  if (!is.null(wmask)) {
    storage.mode(wmask) <- "logical"
  }

  if (!is.null(fmask)) {
    storage.mode(fmask) <- "logical"
  }

  if (is.null(weights)) {
    if (!circular) {
      weights <- .hweights.default(N, L)
    } else {
      weights <- rep(L, N)
    }
  }
  storage.mode(weights) <- "integer"

  .Call("initialize_fft_plan", N, wmask, fmask, weights)
}

is.fft.plan <- function(fft.plan) {
  .Call("is_fft_plan", fft.plan)
}

new.hmat <- function(F, L = (N + 1) %/% 2, circular = FALSE,
                     wmask = NULL, fmask = NULL, weights = NULL,
                     fft.plan = NULL) {
  if (length(circular) > 1) {
    warning("Incorrect argument length: length(circular) > 1, the first value will be used")
    circular <- circular[1]
  }

  N <- length(F)

  if (!is.null(weights)) {
    mask <- weights > 0
    F[!mask] <- mean(F[mask]) # Improve FFT stability & remove NAs
  }

  storage.mode(F) <- "double"
  storage.mode(L) <- "integer"
  storage.mode(circular) <- "logical"


  h <- new("extmat",
           .Call("initialize_hmat", F, L, circular,
                 if (is.null(fft.plan)) fft.plan.1d(N, L, circular, wmask, fmask, weights) else fft.plan))
}

hcols <- function(h) {
  ncol(h)
}

hrows <- function(h) {
  nrow(h)
}

is.hmat <- function(h) {
  is.extmat(h) && .Call("is_hmat", h@.xData)
}

hmatmul <- function(hmat, v, transposed = FALSE) {
  ematmul(hmat, v, transposed = transposed)
}

.traj.dim.ssa <- function(x) {
  c(x$window, x$length - x$window + 1)
}

.traj.dim.1d.ssa <- function(x) {
  Ldim <- sum(x$wmask)
  if (Ldim == 0)
    Ldim <- x$window

  Kdim <- sum(x$fmask)
  if (Kdim == 0)
    Kdim <- x$length - ifelse(x$circular, 0, x$window - 1)

  c(Ldim, Kdim)
}

decompose.ssa <- function(x,
                          neig = NULL,
                          ...,
                          force.continue = FALSE) {
  ## Check, whether continuation of decomposition is requested
  if (!force.continue && nsigma(x) > 0 &&
      !capable(x, "decompose.continue"))
    stop("Continuation of decomposition is not yet implemented for this method.")

  if (is.null(neig))
    neig <- .default.neig(x, ...)

  if (identical(x$svd.method, "svd")) {
    S <- svd(as.matrix(.get.or.create.trajmat(x)), nu = neig, nv = neig)
    .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else if (identical(x$svd.method, "eigen")) {
    S <- eigen(tcrossprod(.get.or.create.trajmat(x)), symmetric = TRUE)

    ## Fix small negative values
    S$values[S$values < 0] <- 0

    .set.decomposition(x,
                       sigma = sqrt(S$values[1:neig]),
                       U = S$vectors[, 1:neig, drop = FALSE])
  } else if (identical(x$svd.method, "nutrlan")) {
    S <- trlan.svd( .get.or.create.trajmat(x), neig = neig, ...,
                   lambda = .sigma(x), U = .U(x))
    .set.decomposition(x, sigma = S$d, U = S$u)
  } else if (identical(x$svd.method, "propack")) {
    S <- propack.svd(.get.or.create.trajmat(x), neig = neig, ...)
    .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else if (identical(x$svd.method, "rspectra")) {
    if (!requireNamespace("RSpectra", quietly = TRUE))
        stop("RSpectra package is required for SVD method `rspectra'")
    h <- .get.or.create.trajmat(x)
    A <- function(x, args) ematmul(args, x)
    Atrans <- function(x, args) ematmul(args, x, transposed = TRUE)
    S <- RSpectra::svds(A, k = neig, Atrans = Atrans, dim = dim(h), args = h, ...)
    .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else if (identical(x$svd.method, "primme")) {
    if (!requireNamespace("PRIMME", quietly = TRUE))
        stop("PRIMME package is required for SVD method `primme'")
    h <- .get.or.create.trajmat(x)
    pA <-function(x, trans) if (identical(trans, "c")) crossprod(h, x) else h %*% x
    S <- PRIMME::svds(pA, NSvals = neig, m = nrow(h), n = ncol(h), isreal = TRUE, ...)
    .set.decomposition(x, sigma = S$d, U = S$u, V = S$v)
  } else
    stop("unsupported SVD method")

  x
}


calc.v.ssa <- function(x, idx, ...) {
  nV <- nv(x)

  V <- matrix(NA_real_, .traj.dim(x)[2], length(idx))
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

    h <- .get.or.create.trajmat(x)
    V[, idx > nV] <- crossprod(h, U) / rep(sigma, each = nrow(V))
  }

  invisible(V)
}

.init.fragment.1d.ssa <- function(this)
  expression({
  if (length(circular) > 1)
    warning("Incorrect argument length: length(circular) > 1, the first value will be used")
  if (length(circular) != 1)
    circular <- circular[1]

  ## Coerce input to vector (we have already saved attrs)
  x <- as.vector(x)
  N <- length(x)

  ## Calculate masks
  mask <- if (is.null(mask)) !is.na(x) else mask & !is.na(x)

  ecall$wmask <- wmask
  if (is.null(wmask)) {
    wmask <- rep(TRUE, L)
  } else {
    L <- length(wmask)
  }

  K <- if (circular) N else N - L + 1

  fmask <- .factor.mask.1d(mask, wmask, circular = circular)

  if (!all(wmask) || !all(fmask) || any(circular)) {
    weights <- .field.weights.1d(wmask, fmask, circular = circular)

    ommited <- sum(mask & (weights == 0))
    if (ommited > 0)
      warning(sprintf("Some field elements were not covered by shaped window. %d elements will be ommited", ommited))

    if (all(weights == 0))
      warning("Nothing to decompose: the given field shape is empty")
  } else {
    weights <- NULL
  }

  if (all(wmask))
    wmask <- NULL
  if (all(fmask))
    fmask <- NULL
})
