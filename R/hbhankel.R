#   R package for Singular Spectrum Analysis
#   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
#   Copyright (c) 2009 Konstantin Usevich <usevich.k.d@gmail.com>
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

#   Routines for hankel-block hankel (aka 2d) SSA


.convolution.dims <- function(x.dim, y.dim, type = "circular") {
  type <- sapply(type, match.arg, choices = c("circular", "open", "filter"))

  stopifnot(length(x.dim) == length(y.dim))
  rank <- length(x.dim)

  if (length(type) != rank) {
    # Use recycling

    if (rank %% length(type) != 0)
      warning("longer object length is not a multiple of shorter object length")

    type <- type[(seq_len(rank) - 1) %% length(type) + 1]
  }

  stopifnot(x.dim[type == "filter"] >= y.dim[type == "filter"])

  output.dim <- ifelse(type == "circular",
                       x.dim,
                       ifelse(type == "open",
                              x.dim + y.dim - 1,
                              x.dim - y.dim + 1))

  input.dim <- ifelse(type == "open",
                      x.dim + y.dim - 1,
                      x.dim)

  names(input.dim) <- names(output.dim) <- names(x.dim)

  list(input.dim = input.dim,
       output.dim = output.dim)
}

.convolven <- function(x, y, conj = TRUE, type = "circular") {
  if (is.null(dim(x))) dim(x) <- length(x)
  if (is.null(dim(y))) dim(y) <- length(y)

  io.dims <- .convolution.dims(dim(x), dim(y), type)
  input.dim <- io.dims$input.dim
  output.dim <- io.dims$output.dim

  storage.mode(x) <- storage.mode(y) <- "double"
  storage.mode(conj) <- "logical"
  storage.mode(input.dim) <- storage.mode(output.dim) <- "integer"

  res <- .Call("convolveN", x, y, input.dim, output.dim, conj)
  if (length(output.dim) == 1)
    dim(res) <- NULL
  else
    dim(res) <- output.dim

  res
}

.factor.mask.2d <- function(field.mask, window.mask, circular = FALSE) {
  field.mask[] <- as.numeric(field.mask)
  window.mask[] <- as.numeric(window.mask)
  tmp <- .convolven(field.mask, window.mask, conj = TRUE,
                   type = ifelse(circular, "circular", "filter"))

  abs(tmp - sum(window.mask)) < 0.5 # ==0, but not exact in case of numeric error
}

.field.weights.2d <- function(window.mask, factor.mask, circular = FALSE) {
  window.mask[] <- as.numeric(window.mask)
  factor.mask[] <- as.numeric(factor.mask)
  res <- .convolven(factor.mask, window.mask, conj = FALSE,
                   type = ifelse(circular, "circular", "open"))
  res[] <- as.integer(round(res))

  res
}

.ball.mask <- function(R, rank) {
  I <- array(seq_len(2*R - 1), dim = rep(2*R - 1, rank))

  Is <- lapply(seq_len(rank),
               function(i) {
                 perm <- seq_len(rank)
                 perm[1] <- i
                 perm[i] <- 1
                 aperm(I, perm)
               })

  dist2 <- array(0, dim = rep(2*R - 1, rank))
  for (I in Is) {
    dist2 <- dist2 + (I - R)^2
  }

  dist2 <= (R - 1)^2
}

.simplex.mask <- function(side, rank) {
  I <- array(seq_len(side), dim = rep(side, rank))

  Is <- lapply(seq_len(rank),
               function(i) {
                 perm <- seq_len(rank)
                 perm[1] <- i
                 perm[i] <- 1
                 aperm(I, perm)
               })

  dist <- array(0, dim = rep(side, rank))
  for (I in Is) {
    dist <- dist + I
  }

  dist <= side + rank - 1
}

new.hbhmat <- function(F, L = (N + 1) %/% 2,
                       wmask = NULL, fmask = NULL, weights = NULL,
                       circular = FALSE) {
  rank <- length(dim(F))
  if (length(circular) > rank)
    warning("Incorrect argument length: length(circular) > rank, two leading values will be used")
  if (length(circular) != rank)
    circular <- circular[(seq_len(rank) - 1) %% length(circular) + 1]

  N <- dim(F)

  storage.mode(F) <- "double"
  storage.mode(L) <- "integer"
  storage.mode(circular) <- "logical"

  if (!is.null(wmask)) {
    storage.mode(wmask) <- "logical"
  }

  if (!is.null(fmask)) {
    storage.mode(fmask) <- "logical"
  }

  if (!is.null(weights)) {
    mask <- weights > 0
    F[!mask] <- mean(F[mask]) # Improve FFT stability & remove NAs
  } else {
    weights <- .hweightsn(N, L)
  }
  storage.mode(weights) <- "integer"

  h <- new("extmat",
           .Call("initialize_hbhmat", F, L, wmask, fmask, weights, circular))
}

hbhcols <- function(h) {
  ncol(h)
}

hbhrows <- function(h) {
  nrow(h)
}

is.hbhmat <- function(h) {
  is.extmat(h) && .Call("is_hbhmat", h@.xData)
}

hbhmatmul <- function(hmat, v, transposed = FALSE) {
  ematmul(hmat, v, transposed = transposed)
}

.get.or.create.hbhmat <- function(x) {
  .get.or.create(x, "hmat",
                 new.hbhmat(.F(x), L = x$window, wmask = x$wmask, fmask = x$fmask,
                            weights = x$weights, circular = x$circular))
}

.get.or.create.trajmat.nd.ssa <- .get.or.create.hbhmat

.traj.dim.nd.ssa <- function(x) {
  Ldim <- sum(x$wmask)
  if (Ldim == 0)
    Ldim <- prod(x$window)

  Kdim <- sum(x$fmask)
  if (Kdim == 0)
    Kdim <- prod(x$length - ifelse(x$circular, 0, x$window - 1))

  c(Ldim, Kdim)
}

.hankelize.one.nd.ssa <- function(x, U, V) {
  h <- .get.or.create.hbhmat(x)
  storage.mode(U) <- storage.mode(V) <- "double"
  .Call("hbhankelize_one_fft", U, V, h@.xData)
}

.init.fragment.2d.ssa <- .init.fragment.nd.ssa <- function(this)
  expression({
    ## Coerce input to array if necessary
    if (!is.array(x))
      x <- as.array(x)
    N <- dim(x)

    rank <- length(dim(x))

    wmask <- .fiface.eval(substitute(wmask),
                          envir = parent.env,
                          circle = function(...) .ball.mask(..., rank = rank),
                          triangle = function(...) .simplex.mask(..., rank = rank))
    ecall$wmask <- wmask
    if (is.null(wmask)) {
      wmask <- array(TRUE, dim = L)
    } else {
      L <- dim(wmask)
    }

    # Fix rank (ndims) of x
    rank <- length(L)
    if (length(dim(x)) < rank)
      dim(x) <- c(dim(x), rep(1, rank - length(dim(x))))

    mask <- if (is.null(mask)) !is.na(x) else mask & !is.na(x)

    # Check `circular' argument
    if (length(circular) > rank)
      warning("Incorrect argument length: length(circular) > rank, the leading values will be used")
    if (length(circular) != rank)
      circular <- circular[(seq_len(rank) - 1) %% length(circular) + 1]

    K <- ifelse(circular, N, N - L + 1)

    fmask <- .factor.mask.2d(mask, wmask, circular = circular)

    if (!all(wmask) || !all(fmask) || any(circular)) {
      weights <- .field.weights.2d(wmask, fmask, circular = circular)

      ommited <- sum(mask & (weights == 0))
      if (ommited > 0) {
        warning(sprintf("Some field elements were not covered by shaped window. %d elements will be ommited", ommited))
      }

      if (all(weights == 0)) {
        warning("Nothing to decompose: the given field shape is empty")
      }
    } else {
      weights <- NULL
    }

    if (all(wmask))
      wmask <- NULL
    if (all(fmask))
      fmask <- NULL
  })
