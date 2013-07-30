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

.get.or.create.mfft.plan <- function(x) {
  .get.or.create(x, "fft.plan", lapply(x$length, fft.plan.1d))
}

.traj.dim.mssa <- function(x) {
  c(x$window, sum(x$length - x$window + 1))
}

.hmat.striped <- function(x, fft.plan) {
  # FIXME: think about NA's in the end
  N <- x$length; L <- x$window; K <- N - L + 1

  F <- .get(x, "F")
  h <- lapply(seq_along(N),
              function(idx) new.hmat(F[[idx]], L = L,
                                     fft.plan = fft.plan[[idx]]))
  b <- c(0, cumsum(K))
  matmul <- function(v) {
    res <- numeric(L)
    for (idx in seq_along(h)) {
      res <- res + hmatmul(h[[idx]], v[(b[idx]+1):b[idx+1]], transposed = FALSE)
    }
    res
  }
  tmatmul <- function(v) unlist(lapply(h, hmatmul, v = v, transposed = TRUE))

  extmat(matmul, tmatmul, nrow = L, ncol = sum(K))
}


.get.or.create.mhmat <- function(x) {
  fft.plan <- .get.or.create.mfft.plan(x)
  .get.or.create(x, "hmat",
                 .hmat.striped(x, fft.plan = fft.plan))
}

decompose.mssa.svd <- function(x,
                               neig = min(L, K),
                               ...,
                               force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  F <- .get(x, "F")
  h <- do.call(cbind, lapply(seq_along(N), function(idx) hankel(F[[idx]], L = L)))

  # Do decomposition
  S <- svd(h, nu = neig, nv = neig)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)
  if (!is.null(S$v))
    .set(x, "V", S$v)

  x
}

decompose.mssa.eigen <- function(x, ...,
                                 force.continue = FALSE) {
  N <- x$length; L <- x$window

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build L-cov matrix
  F <- .get(x, "F")
  fft.plan <- .get.or.create.mfft.plan(x)
  C <- matrix(0, L, L)
  for (idx in seq_along(N)) {
    C <- C + Lcov.matrix(F[[idx]], L = L, fft.plan = fft.plan[[idx]])
  }

  # Do decomposition
  S <- eigen(C, symmetric = TRUE)

  # Fix small negative values
  S$values[S$values < 0] <- 0

  # Save results
  .set(x, "lambda", sqrt(S$values[1:neig]))
  .set(x, "U", S$vectors[, 1:neig, drop = FALSE])

  x
}

decompose.mssa.propack <- function(x,
                                   neig = min(50, L, prod(K)),
                                   ...,
                                   force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.")

  h <- .get.or.create.mhmat(x)
  S <- propack.svd(h, neig = neig, ...)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)
  if (!is.null(S$v))
    .set(x, "V", S$v)

  x
}

decompose.mssa.nutrlan <- function(x,
                                   neig = min(50, L, prod(K)),
                                   ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  h <- .get.or.create.mhmat(x)

  lambda <- .get(x, "lambda", allow.null = TRUE)
  U <- .get(x, "U", allow.null = TRUE)

  S <- trlan.svd(h, neig = neig, ...,
                 lambda = lambda, U = U)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)

  x
}

.init.mssa <- function(x, ...) {
  # Initialize FFT plan
  .get.or.create.mfft.plan(x)
}

calc.v.mssa<- function(x, idx, env = .GlobalEnv, ...) {
  lambda <- .get(x, "lambda")[idx]
  U <- .get(x, "U")[, idx, drop = FALSE]
  h <- .get.or.create.mhmat(x)

  invisible(sapply(1:length(idx),
                   function(i) ematmul(h, U[, i], transposed = TRUE) / lambda[i]))
}

.hankelize.one.mssa <- function(x, U, V) {
  N <- x$length; L <- x$window; K <- N - L + 1

  b <- c(0, cumsum(K))
  fft.plan <- .get.or.create.mfft.plan(x)

  # FIXME: All these apply's are really ugly. Switch to C version...
  unlist(lapply(seq_along(K),
                function(idx) .hankelize.one.1d.ssa(x, U, V[(b[idx]+1):b[idx+1]], fft.plan[[idx]])))
}

.elseries.mssa <- function(x, idx, ..., env = .GlobalEnv) {
  if (max(idx) > nlambda(x))
    stop("Too few eigentriples computed for this decomposition")

  N <- x$length
  lambda <- .get(x, "lambda")
  U <- .get(x, "U")
  F <- .get(x, "F")

  res <- numeric(sum(N))
  for (i in idx) {
    if (nv(x) >= i) {
      # FIXME: Check, whether we have factor vectors for reconstruction
      # FIXME: Get rid of .get call
      V <- .get(x, "V")[, i]
    } else {
      # No factor vectors available. Calculate them on-fly.
      V <- calc.v(x, i, env = env)
    }

    res <- res + lambda[i] * .hankelize.one(x, U = U[, i], V = V)
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
  a <- (if (drop) NULL else .get(x, "Fattr"))
  cls <- (if (drop) NULL else .get(x, "Fclass"))

  # MSSA is a bit different from the default case. We need to convert (if
  # possible) to original object
  stopifnot(inherits(F, "series.list"))

  # First, pad with NA's if necessary
  F <- lapply(F,
              function(x) {
                removed <- attr(x, "na.action")
                if (!is.null(removed)) {
                  res <- numeric(length(x) + length(removed))
                  res[removed] <- NA
                  res[-removed] <- x
                  res
                } else
                x
              })
  # Optionaly convert to matrix
  if ("matrix" %in% cls)
    F <- simplify2array(F)

  if (fixup) {
     # Try to guess the indices of known time series classes
    if ("ts" %in% cls) {
      tsp <- a$tsp
      return (ts(F,
                 start = if (only.new) tsp[2] + 1/tsp[3] else tsp[1],
                 frequency = tsp[3]))
    } else if (!is.null(cls)) {
      warning("do not know how to fixup attributes for this input")
    }
  } else {
    # Restore attributes
    attributes(F) <- a
  }

  F
}


.get.or.create.cfft.plan <- function(x) {
  .get.or.create(x, "fft.plan", fft.plan.1d(x$length))
}

.traj.dim.cssa <- function(x) {
  c(x$window, sum(x$length - x$window + 1))
}

.chmat <- function(x, fft.plan) {
  N <- x$length; L <- x$window; K <- N - L + 1
  F <- .get(x, "F")

  R <- new.hmat(Re(F), L = L, fft.plan = fft.plan)
  I <- new.hmat(Im(F), L = L, fft.plan = fft.plan)

  matmul <- function(v)
    c(hmatmul(R, v[1:K], transposed = FALSE) - hmatmul(I, v[(K+1):(2*K)], transposed = FALSE),
      hmatmul(I, v[1:K], transposed = FALSE) + hmatmul(R, v[(K+1):(2*K)], transposed = FALSE))
  tmatmul <- function(v)
    c( hmatmul(R, v[1:L], transposed = TRUE) + hmatmul(I, v[(L+1):(2*L)], transposed = TRUE),
      -hmatmul(I, v[1:L], transposed = TRUE) + hmatmul(R, v[(L+1):(2*L)], transposed = TRUE))

  extmat(matmul, tmatmul, nrow = 2*L, ncol = 2*K)
}

.get.or.create.chmat <- function(x) {
  .get.or.create(x, "hmat",
                 .chmat(x, fft.plan = .get.or.create.cfft.plan(x)))
}

decompose.cssa.svd <- function(x,
                               neig = min(L, K),
                               ...,
                               force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  F <- .get(x, "F")
  h <- hankel(F, L = L)

  # Do decomposition
  S <- svd(h, nu = neig, nv = neig)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)
  if (!is.null(S$v))
    .set(x, "V", S$v)

  x
}

.traj.dim.cssa.svd <- function(x) {
  c(x$window, x$length - x$window + 1)
}

cssa.to.complex <- function(values, vectors) {
  # First, make sure values come into the pairs
  d1 <- values[c(TRUE, FALSE)]
  d2 <- values[c(FALSE, TRUE)]
  if (any((d1 - d2) / d2 > 1e-3 & d2 > 0))
    warning("Too big difference between consecutive eigenvalues. CSSA might not converge")

  # And vectors
  Y1 <- vectors[1:(nrow(vectors)/2),    c(TRUE,  FALSE)]
  Z1 <- vectors[-(1:(nrow(vectors)/2)), c(TRUE,  FALSE)]
  Y2 <- vectors[1:(nrow(vectors)/2),    c(FALSE, TRUE)]
  Z2 <- vectors[-(1:(nrow(vectors)/2)), c(FALSE, TRUE)]

  V1 <- Y1 + 1i*Z1
  V2 <- Y2 + 1i*Z2

  # Sanity check
  if (any((Mod(V1 - 1i*V2) > 1e-6) & (Mod(V1 + 1i*V2) > 1e-6)))
    warning("Too big difference between consecutive eigenvectors. CSSA might not converge")

  list(d = d2, u = V2
       # , vectors2 = V2,
       # dd = (d1 - d2) / d2 > 1e-3 & d2 > 0,
       # vv = (Mod(V1 - 1i*V2) > 1e-6) & (Mod(V1 + 1i*V2) > 1e-6)
       )
}

decompose.cssa.eigen <- function(x, ...,
                                 force.continue = FALSE) {
  N <- x$length; L <- x$window

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decomposition is not supported for this method.")

  # Build hankel matrix
  F <- .get(x, "F")

  R <- hankel(Re(F), L = L)
  I <- hankel(Im(F), L = L)
  h <- cbind(rbind(R, I), rbind(-I, R))

  # Do decomposition
  # FIXME: Build the L-covariance matrix properly
  S <- eigen(tcrossprod(h, h), symmetric = TRUE)

  # Fix small negative values
  S$values[S$values < 0] <- 0

  S <- cssa.to.complex(sqrt(S$values), S$vectors)

  # Save results
  .set(x, "lambda", sqrt(S$values[1:neig]))
  .set(x, "U", S$vectors[, 1:neig, drop = FALSE])

  x
}

decompose.cssa.propack <- function(x,
                                   neig = min(50, L, K),
                                   ...,
                                   force.continue = FALSE) {
  N <- x$length; L <- x$window; K <- N - L + 1

  # Check, whether continuation of decomposition is requested
  if (!force.continue && nlambda(x) > 0)
    stop("Continuation of decompostion is not yet implemented for this method.")

  h <- .get.or.create.chmat(x)
  S <- propack.svd(h, neig = 2*neig, ...)

  S <- cssa.to.complex(S$d, S$u)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)
  if (!is.null(S$v))
    .set(x, "V", S$v)

  x
}

decompose.cssa.nutrlan <- function(x,
                                   neig = min(50, L, K),
                                   ...) {
  N <- x$length; L <- x$window; K <- N - L + 1

  h <- .get.or.create.chmat(x)

  lambda <- .get(x, "lambda", allow.null = TRUE)
  U <- .get(x, "U", allow.null = TRUE)

  S <- trlan.svd(h, neig = 2*neig, ...,
                 lambda = lambda, U = U)

  S <- cssa.to.complex(S$d, S$u)

  # Save results
  .set(x, "lambda", S$d)
  if (!is.null(S$u))
    .set(x, "U", S$u)

  x
}

.traj.dim.cssa <- function(x) {
  c(2*x$window, 2*(x$length - x$window + 1))
}

.init.cssa <- function(x, ...) {
  # Initialize FFT plan
  .get.or.create.cfft.plan(x)
}

calc.v.cssa<- function(x, idx, env = .GlobalEnv, ...) {
  lambda <- .get(x, "lambda")[idx]
  U <- .get(x, "U")[, idx, drop = FALSE]
  h <- .get.or.create.chmat(x)

  invisible(sapply(1:length(idx),
                   function(i) {
                     v <- ematmul(h, c(Re(U[, i]), Im(U[, i])), transposed = TRUE) / lambda[i]
                     v[1:(length(v) / 2)] + 1i*v[-(1:(length(v) / 2))]
                   }))
}

.hankelize.one.cssa <- function(x, U, V, fft.plan = NULL) {
  fft.plan <- (if (is.null(fft.plan)) fft.plan.1d(x$length) else fft.plan)
  R1 <- .hankelize.one.default(Re(U), Re(V), fft.plan = fft.plan)
  R2 <- .hankelize.one.default(Im(U), Im(V), fft.plan = fft.plan)
  I1 <- .hankelize.one.default(Re(U), Im(V), fft.plan = fft.plan)
  I2 <- .hankelize.one.default(Im(U), Re(V), fft.plan = fft.plan)

  (R1 + R2) + 1i*(-I1 + I2)
}
