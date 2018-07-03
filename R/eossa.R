#   R package for Singular Spectrum Analysis
#   Copyright (c) 2017-2018 Alex Shlemov <shlemovalex@gmail.com>
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

#   Routines for ESPRIT-based Oblique SSA


# TODO use QR instead of SVD and make basis real before clustering and joining
.clust.basis <- function(U, roots, k = 2, h = NULL, order = FALSE) {
  # Reorder roots by their freqs
  ord <- order(abs(Arg(roots)))
  roots <- roots[ord]
  U <- U[, ord, drop = FALSE]

  stopifnot(length(k) == 1)

  # Check for argument k, k <= the number of roots with nonegative imagine part
  maxk <- sum(Im(roots) >= -.Machine$double.eps) # Maybe just Im >= 0?
  if (k > maxk) {
    stop(sprintf("k exceeds the number of different ESPRIT roots with non-negative imaginary parts (%d%)", maxk))
  }

  d <- stats::dist(cbind(Re(roots), abs(Im(roots))), method = "euclidian") # TODO Use the proper distance from KDU

  hc <- hclust(d, method = "complete")

  idx <- cutree(hc, k = k, h = h)

  groups <- tapply(seq_along(idx), idx, identity)
  names(groups) <- paste("F", names(groups), sep = "")

  for (group in groups) {
    U[, group] <- svd(cbind(Re(U[, group, drop = FALSE]),
                                          Im(U[, group, drop = FALSE])),
                                    nu = length(group), nv = 0)$u
  }

  U <- Re(U)

  if (order) {
    U <- U[, unlist(groups), drop = FALSE]
    l <- sapply(groups, length)
    cl <- cumsum(c(0, l))
    groups <- lapply(seq_along(l), function(i) (cl[i] + 1) : cl[i+1])
  }

  list(basis = U, groups = groups)
}

eossa <- function(x, ...)
  UseMethod("eossa")

eossa.ssa <- function(x,
                      nested.groups, k = 2,
                      subspace = c("column", "row"),
                      dimensions = NULL,
                      solve.method = c("ls", "tls"),
                      beta = 8,
                      ...) {
  if (missing(nested.groups))
    nested.groups <- as.list(1:min(nsigma(x), nu(x)))

  subspace <- match.arg(subspace)
  solve.method <- match.arg(solve.method)

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = nested.groups, ...)

  idx <- sort(unique(unlist(nested.groups)))
  triples <- .get.orth.triples(x, idx, do.orthogonalize = FALSE)
  osigma <- triples$sigma; U <- triples$U; V <- triples$V

  if (identical(subspace, "column")) {
    vectors <- U
    V <- V * rep(osigma, each = nrow(V))
    wmask <- .wmask(x)
  } else if (identical(subspace, "row")) {
    vectors <- V
    U <- U * rep(sqrt(osigma), each = nrow(U))
    wmask <- .fmask(x)
  }
  sigma <- rep(1, length(idx))

  if (is.null(dimensions)) {
    dimensions <- seq_len(.dim(x))
  }

  d <- dim(wmask)

  if (max(dimensions) > length(d)) {
    stop(sprintf("some of input dimension indices exceed the actual number of object dimensions (%d)",
                 length(d)))
  }

  Zs <- lapply(dimensions,
               function(ndim) {
                 .shift.matrix(vectors,
                               wmask = wmask,
                               ndim = ndim,
                               circular = x$circular[ndim],
                               solve.method = solve.method)
               })

  Z <- .matrix.linear.combination(Zs, beta)
  Ze <- eigen(Z, symmetric = FALSE)

  # sm <- 0.5 * (Usm + t(Vsm)) # TODO implement two-sided ESPRIT????
  mb <- .clust.basis(Ze$vectors, Ze$values, k = k)
  C <- mb$basis
  nested.groups <- mb$groups

  U <- U %*% C
  # V <- V %*% solve(t(C))  # TODO Use qr.solve here
  V <- t(qr.solve(C, t(V)))

  x <- clone(x, copy.cache = FALSE) # TODO Maybe we should to preserve the relevant part of the cache?
  .save.oblique.decomposition(x, sigma, U, V, idx)

  # Return to real group numbers
  nested.groups <- lapply(nested.groups, function(group) idx[group])

  # Grab old iossa.groups.all value
  iossa.groups.all <- .get(x, "iossa.groups.all", allow.null = TRUE)
  if (is.null(iossa.groups.all)) {
    iossa.groups.all <- list()
  }

  valid.groups <- as.logical(sapply(iossa.groups.all,
                                    function(group) length(intersect(group, idx)) == 0))
  .set(x, "iossa.groups",  nested.groups)
  .set(x, "iossa.groups.all", c(nested.groups, iossa.groups.all[valid.groups]))

  # Save nested components
  .set(x, "ossa.set", idx)

  if (!is.null(.decomposition(x, "nPR"))) {
    if (any(idx <= .decomposition(x, "nPR"))) {
      .set.decomposition(x, nPR = 0, nPL = 0)
    } else if (any(idx <= sum(unlist(.decomposition(x, c("nPR", "nPL")))))){
      .set.decomposition(x, nPL = 0)
    }
  }

  if (!inherits(x, "ossa")) {
    class(x) <- c("ossa", class(x))
  }

  # Save call info
  x$call <- match.call()

  invisible(x)
}
