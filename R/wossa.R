#   R package for Singular Spectrum Analysis
#   Copyright (c) 2015 Alex Shlemov <shlemovalex@gmail.com>
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

#   Routines for weighted SSA

.identity.emat <- function(n) {
  extmat(identity, identity, n, n)
}

.mask.emat <- function(mask) {
  mask <- as.logical(as.vector(mask))

  matmul <- function(v) {
    v[mask]
  }

  tmatmul <- function(v) {
    res <- numeric(lenght(mask))
    res[mask] <- v
  }

  extmat(matmul, tmatmul, sum(mask), length(mask))
}

.prod.emat <- function(A, B) {
  stopifnot(ncol(A) == nrow(B))

  matmul <- function(v) {
    A %*% (B %*% v)
  }

  tmatmul <- function(v) {
    v %*% A %*% B
  }

  extmat(matmul, tmatmul, nrow(A), ncol(B))
}

.diag.emat <- function(d) {
  d <- as.vector(d)

  matmul <- function(v) {
    v * d
  }

  extmat(matmul, matmul, length(d), length(d))
}

.whmat <- function(x) {
  hmat <- .get.or.create.trajmat(x)
  column.oblique <- .get(x, "column.oblique")[[1]]
  row.oblique <- .get(x, "row.oblique")[[1]]

  matmul <- function(v) {
    v <- v %*% row.oblique
    v <- hmatmul(hmat, v, transposed = FALSE)
    v <- column.oblique %*% v

    v
  }

  tmatmul <- function(v) {
    v <- v %*% column.oblique
    v <- hmatmul(hmat, v, transposed = TRUE)
    v <- row.oblique %*% v

    v
  }

  if (length(.get(x, "column.oblique")) > 2 && length(.get(x, "row.oblique")) > 2) {
    column.oblique <- .get(x, "column.oblique")[[3]]
    row.oblique <- .get(x, "row.oblique")[[3]]

    matmul <- function(v) {
      v <- v * row.oblique
      v <- hmatmul(hmat, v, transposed = FALSE)
      v <- column.oblique * v

      v
    }

    tmatmul <- function(v) {
      v <- v * column.oblique
      v <- hmatmul(hmat, v, transposed = TRUE)
      v <- row.oblique * v

      v
    }
  }

  extmat(matmul, tmatmul, hrows(hmat), hcols(hmat))
}

.get.or.create.whmat <- function(x)
  .get.or.create(x, "whmat", .whmat(x))


decompose.wossa <- function(x,
                            neig = NULL,
                            ...,
                            force.continue = FALSE) {
  ## Check, whether continuation of decomposition is requested
  ## FIXME: Check the caps
  if (!force.continue && nsigma(x) > nspecial(x) &&
       !identical(x$svd.method, "nutrlan"))
    stop("Continuation of decomposition is not yet implemented for this method.")

  if (is.null(neig))
    neig <- .default.neig(x, ...)

  column.ioblique <- .get(x, "column.oblique")[[2]]
  row.ioblique <- .get(x, "row.oblique")[[2]]

  if (identical(x$svd.method, "svd")) {
    S <- svd(as.matrix(.get.or.create.whmat(x)), nu = neig, nv = neig)
    oU <- S$u; oV <- S$v; osigma <- S$d
  } else if (identical(x$svd.method, "eigen")) {
    S <- eigen(tcrossprod(.get.or.create.whmat(x)), symmetric = TRUE)

    ## Fix small negative values
    S$values[S$values < 0] <- 0
    oU <- S$vectors; oV <- NULL; osigma <- sqrt(S$values)
  } else if (identical(x$svd.method, "propack")) {
    S <- propack.svd(.get.or.create.whmat(x), neig = neig, ...)
    U <- S$u; V <- S$v; sigma <- S$d
  } else if (identical(x$svd.method, "nutrlan")) {
    S <- trlan.svd(.get.or.create.phmat(x), neig = neig, ...,
                   lambda = .decomposition(x)$osigma, U = .decomposition(x)$oU)
    oU <- S$u; oV <- NULL; osigma <- S$d
  } else
    stop("unsupported SVD method")

  if (is.null(oV)) {
    Z <- crossprod(.get.or.create.whmat(x), oU)
    oV <- sweep(Z, 2, osigma, FUN = "/")
  }

  U <- column.ioblique %*% oU
  V <- row.ioblique %*% oV

  sU <- apply(U, 2, function(x) sqrt(sum(x^2)))
  sV <- apply(V, 2, function(x) sqrt(sum(x^2)))
  U <- sweep(U, 2, sU, FUN = "/")
  V <- sweep(V, 2, sV, FUN = "/")
  sigma <- osigma * sU * sV

  # Precompute weights
  column.oblique <- .get(x, "column.oblique")[[3]]
  row.oblique <- .get(x, "row.oblique")[[3]]
  weights <- .hankelize.one(x, column.oblique^2, row.oblique^2)

  .set.decomposition(x,
                     sigma = sigma, U = U, V = V,
                     oU = oU, osigma = osigma,
                     weights = weights, # TODO Mb use it as genereal `weights`
                     kind = "weighted.oblique.decomposition")

  x
}

.init.fragment.wossa <- function(this)
  expression({
    ## First, initialize the main object
    ## We cannot use NextMethod here due to non-standard evaluation
    class.wo.wossa <- class(this)[!grepl("^wossa", class(this))]
    eval(getS3method(".init.fragment", class.wo.wossa)(this))

    .vector.pseudo.inverse <- function(v, eps = 1e-6) {
      iv <- 1 / v
      iv[abs(v) < eps] <- 0

      iv
    }

    .oblique.matrix.coerce <- function(A, n) {
      if (is.null(A) || identical(A, "identity")) {
        I <- .identity.emat(n)
        list(I, I, rep(1, n), rep(1, n))
      } else if (is.vector(A) && !is.list(A)) {
        A <- sqrt(A)
        stopifnot(length(A) == n)
        iA <- .vector.pseudo.inverse(A)
        list(.diag.emat(A), .diag.emat(iA), A, iA)
      } else if (is.list(A)) {
        stopifnot(ncol(A[[1]]) == n)
        stopifnot(nrow(A[[2]]) == n)
        stopifnot(ncol(A[[2]]) == nrow(A[[1]]))

        stop("Non-diagonal oblique SSA is not implemented yet")

        A
      } else if (is.matrix(A) || is.extmat(A)) {
        stopifnot(ncol(A) == n)
        iA = pseudo.inverse(as.matrix(A))

        stop("Non-diagonal oblique SSA is not implemented yet")

        list(A, iA)
      } else {
        stop("Unknown type of oblique input object")
      }
    }

    ## TODO Think about MSSA case
    column.oblique <- .oblique.matrix.coerce(column.oblique, prod(L))
    row.oblique <- .oblique.matrix.coerce(row.oblique, prod(K))

    ## Shape oblique matrices if needed
    if (!is.null(wmask)) {
      stop("Weighted-oblique decomposition is not implemented yet")
      # column.oblique[[1]] <- .prod.emat(column.oblique[[1]], gcc)
    }
    if (!is.null(fmask)) {
      stop("Weighted-oblique decomposition is not implemented yet")
      # column.weight <- column.weight[as.vector(fmask)]
    }
  })

.colspan.wossa <- function(x, idx) {
  qr.Q(qr(.U(x)[, idx, drop = FALSE]))
}

.rowspan.wossa <- function(x, idx) {
  qr.Q(qr(calc.v(x, idx)))
}

.elseries.wossa <- function(x, idx, ...) {
  if (max(idx) > nsigma(x))
    stop("Too few eigentriples computed for this decomposition")

  dec <- .decomposition(x)
  sigma <- .sigma(dec)
  U <- .U(dec)

  column.oblique <- .get(x, "column.oblique")[[3]]
  row.oblique <- .get(x, "row.oblique")[[3]]
  weights <- .decomposition(x, "weights")

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

    res <- res + sigma[i] * .hankelize.one(x, U = U[, i] * column.oblique^2, 
                                           V = V * row.oblique^2) / weights
  }

  res;
}
