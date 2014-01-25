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


Cond <- function(A) {
  # Condition number for basis
  d <- svd(A)$d
  d[1] / d[length(d)]
}

high.rank.rate <- function(F, L, r, ...) {
  N <- length(F); K <- N - L + 1
  ss <- ssa(F, L, neig = min(r + 1, L, K), ...)

  1 - sum(.sigma(ss)[seq_len(r)]^2) / wnorm(F, L)^2
}

pseudo.inverse <- function(A) {
  # Moore-Penrose pseudo-inverse
  qrA <- qr(A)
  solve(qr.R(qrA), t(qr.Q(qrA)))
}

orthogonalize <- function(Y, Z, sigma, side = c("bi", "left", "right"), normalize = TRUE) {
  side <- match.arg(side)

  if (missing(sigma)) {
    sigma <- rep(1, ncol(Y))
  }

  rank <- length(sigma)
  # Check parameters' dims
  stopifnot(ncol(Y) == rank, ncol(Z) == rank)

  if (identical(side, "bi")) {
    # Low-rank SVD
    qrY <- qr(Y)
    qrZ <- qr(Z)
    dec <- svd(qr.R(qrY) %*% (sigma * t(qr.R(qrZ))))

    list(d = dec$d, u = qr.Q(qrY) %*% dec$u, v = qr.Q(qrZ) %*% dec$v)
  } else if (identical(side, "left")) {
    qrY <- qr(Y)
    u <- qr.Q(qrY)
    v <- Z %*% (sigma * t(qr.R(qrY)))
    if (normalize) {
      d <- sqrt(colSums(v^2))
      v <- v / matrix(d, nrow = nrow(v), ncol = rank, byrow = TRUE)
    } else {
      d <- rep(1, rank)
    }

    list(d = d, u = u, v = v)
  } else if (identical(side, "right")) {
    dec <- Recall(Y = Z, Z = Y, sigma = sigma, side = "left", normalize = normalize)

    list(d = dec$d, u = dec$v, v = dec$u)
  }
}

.gwcor <- function(X, L, LM, RM) {
  fft.plan <- fft.plan.1d(nrow(X))
  mx <- apply(X, 2,
              function(v) {
                h <- new.hmat(v, L = L, fft.plan = fft.plan)
                as.vector(LM %*% hmatmul.wrap(h, t(RM), transposed = FALSE))
              })

  # Compute covariations
  cov <- crossprod(mx)

  # Convert to correlations
  cor <- cov2cor(cov)

  # Fix possible numeric error
  cor[cor > 1] <- 1; cor[cor < -1] <- -1

  # Add class
  class(cor) <- "wcor.matrix"

  # Return
  cor
}

print.iossa.result <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  max.abs.nodiag <- function(mx) {
    diag(mx) <- 0
    mx[which.max(abs(mx))]
  }

  cat("\nI-OSSA result:\n")
  cat("\tConverged:             ", ifelse(x$converged, "yes", "no"), "\n", sep = "")
  cat("\tIterations:            ", x$iter, "\n", sep = "")
  cat("\tCondition numbers:     ", paste0(format(x$Cond, digits = digits), collapse = ", "), "\n", sep = "")
  cat("\tInitial mean(tau):     ", format(mean(x$initial.hrr), digits = digits), "\n", sep = "")
  cat("\tInitial tau:           ", paste0(format(x$initial.hrr, digits = digits), collapse = ", "), "\n", sep = "")
  cat("\tI-OSSA mean(tau):      ", format(mean(x$hrr), digits = digits), "\n", sep = "")
  cat("\tI-OSSA tau:            ", paste0(format(x$hrr, digits = digits), collapse = ", "), "\n", sep = "")
  cat("\tInitial max wcor:      ", format(max.abs.nodiag(x$initial.wcor), digits = digits), "\n", sep = "")
  cat("\tI-OSSA max wcor:       ", format(max.abs.nodiag(x$wcor), digits = digits), "\n", sep = "")
  cat("\n")
}

summary.iossa.result <- function(object, digits = max(3, getOption("digits") - 3), ...)
  print.iossa.result(x = object, digits = digits, ...)

print.ossa <- function(x, ...) {
  NextMethod()

  iossa.result <- .get(x, "iossa.result", allow.null = TRUE)
  if (!is.null(iossa.result))
    print(iossa.result)
}

.save.oblique.decomposition <- function(x, nosigma, Y, Z, idx) {
  sigma <- .sigma(x)
  U <- .U(x)
  V <- if (nv(x) < max(idx)) calc.v.ossa(x, seq_len(max(idx))) else .V(x)

  ynorms <- sqrt(colSums(Y^2))
  znorms <- sqrt(colSums(Z^2))

  nosigma <- nosigma * ynorms * znorms
  Y <- Y / rep(ynorms, each = nrow(Y))
  Z <- Z / rep(znorms, each = nrow(Z))

  for (i in seq_along(idx)) {
    sigma[idx[i]] <- nosigma[i]
    U[, idx[i]] <- Y[, i]
    V[, idx[i]] <- Z[, i]
  }

  .set.decomposition(x, sigma = sigma, U = U, V = V)
}

.make.ossa.result <- function(x, Fs, ranks, IBL, IBR, initial.Fs, svd.method = "auto") {
  L <- x$window

  hrr <- sapply(seq_along(Fs),
                function(i) {
                  high.rank.rate(Fs[[i]], L, ranks[i], svd.method = svd.method)
                })

  initial.hrr <- sapply(seq_along(initial.Fs),
                        function(i) {
                          high.rank.rate(initial.Fs[[i]], L, ranks[i], svd.method = svd.method)
                        })

  names(Fs) <- paste("F", seq_along(Fs), sep = "")
  out <- list(F = .F(x),
              Fs = Fs,
              Cond = c(Cond(IBL), Cond(IBR)),
              gwcor = .gwcor(do.call(cbind, Fs), L, IBL, IBR),
              wcor = wcor(do.call(cbind, Fs), L),
              initial.wcor = wcor(do.call(cbind, initial.Fs), L),
              hrr = hrr,
              initial.hrr = initial.hrr,
              initial.rec = initial.Fs)
  class(out) <- "iossa.result"

  invisible(out)
}

.component.grouping <- function(elem.components, groups) {
  lapply(groups,
         function(group) {
           rowSums(elem.components[, group, drop = FALSE])
         })
}

svd2LRsvd <- function(d, u, v, basis.L, basis.R, need.project = TRUE, fast = TRUE) {
  rank <- length(d)
  # Check parameters' dims
  stopifnot(ncol(u) == rank, ncol(v) == rank)

  ub.L <- if (is.null(basis.L)) {
        basis.L <- u
        diag(nrow = rank, ncol = rank)
      } else {
        crossprod(u, basis.L)
      }

  vb.R <- if (is.null(basis.R)) {
        basis.R <- v
        diag(nrow = rank, ncol = rank)
      } else {
        crossprod(v, basis.R)
      }

  if (fast) {
    dec <- svd(t(ub.L) %*% (vb.R / d))
    dec$d <- 1 / dec$d
    # Reverse order
    dec$d <- rev(dec$d)
    dec$u <- dec$u[, rank:1, drop = FALSE]
    dec$v <- dec$v[, rank:1, drop = FALSE]
  } else {
    dec <- svd(solve(ub.L) %*% (d * t(solve(vb.R))))
  }

  if (need.project) {
    basis.L <- u %*% ub.L
    basis.R <- v %*% vb.R
  }

  list(sigma = dec$d, Y = basis.L %*% dec$u, Z = basis.R %*% dec$v)
}

.get.orth.triples <- function(x, idx, eps = sqrt(.Machine$double.eps)) {
  # Determine the upper bound of desired eigentriples
  desired <- max(idx)

  # Continue decomposition, if necessary
  if (desired > min(nsigma(x), nu(x)))
    decompose(x, neig = desired)

  sigma <- .sigma(x)[idx]
  U <- .U(x)[, idx, drop = FALSE]
  V <- if (nv(x) < desired) calc.v.ossa(x, idx) else .V(x)[, idx, drop = FALSE]

  # TODO Perform orthogonolize if it's only needed
  dec <- orthogonalize(U, V, sigma, side = "bi")
  sigma <- dec$d; U <- dec$u; V <- dec$v

  if (min(sigma) < eps)
    warning("Decomposition isn't minimal. Some singular values equal to zero")

  list(sigma = sigma, U = U, V = V)
}

iossa <- function(x, nested.groups, ..., tol = 1e-5, kappa = 1.2,
                  maxiter = 1000,
                  initial.approx = reconstruct(x, nested.groups),
                  norm = function(x) sqrt(mean(x^2)),
                  trace = FALSE,
                  kappa.balance = 0.5,
                  svd.method = "auto") {
  N <- x$length; L <- x$window; K <- N - L + 1

  if (missing(nested.groups))
    nested.groups <- as.list(1:min(nsigma(x), nu(x)))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = nested.groups, ...)

  rec <- reconstruct(x, groups = nested.groups)

  Fs <- initial.approx

  idx <- sort(unique(unlist(nested.groups)))
  triples <- .get.orth.triples(x, idx)
  lambda <- triples$sigma; U <- triples$U; V <- triples$V

  # Replace (in nested.groups) ET's numbers for their ranks (order numbers) in set of signal ETs
  nested.groups <- lapply(nested.groups, function(group) match(group, idx))
  ranks <- sapply(nested.groups, length)

  # If we use reordering, reorder components
  if (!is.null(kappa)) {
    cumranks <- cumsum(ranks)
    for (i in seq_along(nested.groups)) {
      nested.groups[[i]] <- seq(to = cumranks[i], length.out = ranks[i])
    }
  }

  # Grab the FFT plan
  fft.plan <- .get.or.create.fft.plan(x)

  converged <- FALSE
  for (iter in seq_len(maxiter)) {
    lrss <- numeric(length(nested.groups))
    triples.list <- list()
    for (i in seq_along(nested.groups)) {
      cur.ssa <- ssa(Fs[[i]], L, svd.method = svd.method, neig = min(ranks[i] + 1, L, K))
      lrss[i] <- sum(cur.ssa$sigma[-seq_len(ranks[i])]^2)
      triples.list[[i]] <- .get.orth.triples(cur.ssa, seq_len(ranks[i]))
    }

    if (trace) cat(sprintf("LRSS(%d): %s\n", iter, paste0(sqrt(lrss), collapse = " ")))

    lambdas <- lapply(triples.list, function(el) el$sigma)
    Us <- lapply(triples.list, function(el) el$U)
    Vs <- lapply(triples.list, function(el) el$V)

    if (!is.null(kappa)) {
      # If we use reordering, force separability
      for (i in seq_along(nested.groups[-length(nested.groups)])) {
        div <- lambdas[[i]][ranks[i]] / lambdas[[i + 1]][1]
        if (div < kappa) {
          mul <- kappa / div
          lambdas[[i + 1]] <- lambdas[[i + 1]] / mul
          Us[[i + 1]] <- Us[[i + 1]] * mul^kappa.balance
          Vs[[i + 1]] <- Vs[[i + 1]] * mul^(1 - kappa.balance)
        }
      }
    }

    basU <- do.call(cbind, Us)
    basV <- do.call(cbind, Vs)

    # Oblique SVD
    dec <- svd2LRsvd(lambda, U, V, basU, basV, need.project = TRUE)
    sigma <- dec$sigma; Y <- dec$Y; Z <- dec$Z

    elem.components <- .hankelize.multi.default(t(sigma * t(Y)), Z, fft.plan)
    Fs.new <- .component.grouping(elem.components, nested.groups)

    if (trace) {
      Xs <- lapply(nested.groups,
                   function(group) Y[, group, drop = FALSE] %*% (sigma[group] * t(Z[, group, drop = FALSE])))
      svddist <- sapply(seq_along(Xs),
                        function(i) sum((Xs[[i]] - hankel(Fs[[i]], L)) ^ 2))
      cat(sprintf("SVDD(%d): %s\n", iter, paste0(sqrt(svddist), collapse = " ")))

      nonhankeleness <- sapply(seq_along(Xs),
                               function(i) sum((Xs[[i]] - hankel(Fs.new[[i]], L)) ^ 2))
      cat(sprintf("NHNK(%d): %s\n", iter, paste0(sqrt(nonhankeleness), collapse = " ")))
    }

    deltas <- sapply(seq_along(nested.groups), function(i) norm(Fs.new[[i]] - Fs[[i]]))
    if (max(deltas) < tol) {
      converged <- TRUE
      Fs <- Fs.new
      break
    }

    Fs <- Fs.new
  }

  IBL <- pseudo.inverse(Y)
  IBR <- pseudo.inverse(Z)

  x <- clone(x, copy.cache = FALSE) # TODO Maybe preserve relevant part of cache?
  .save.oblique.decomposition(x, sigma, Y, Z, idx)

  # Update class for x
  if (!inherits(x, "ossa")) {
    class(x) <- c("ossa", class(x))
  }

  # Save call info
  x$call <- match.call()

  out <- .make.ossa.result(x, Fs, ranks, IBL, IBR, rec, svd.method = svd.method)
  out[c("iter", "converged", "kappa", "maxiter", "tol", "call")] <-
      list(iter, converged, kappa, maxiter, tol, match.call())

  .set(x, "iossa.result", out)

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

  invisible(x)
}

fssa <- function(x, nested.groups, FILTER = diff, gamma = 1, ...) {
  if (!is.function(FILTER)) {
    FILTER.coeffs <- FILTER
    FILTER <- function(x) {
      out <- filter(x, FILTER.coeffs)
      out[!is.na(out)]
    }
  }

  N <- x$length; L <- x$window; K <- N - L + 1

  if (missing(nested.groups))
    nested.groups <- as.list(1:min(nsigma(x), nu(x)))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = nested.groups, ...)

  idx <- sort(unique(unlist(nested.groups)))
  triples <- .get.orth.triples(x, idx)
  lambda <- triples$sigma; U <- triples$U; V <- triples$V

  Z <- V * rep(lambda, each = nrow(V))

  fZ <- apply(Z, 2, FILTER)
  dec <- eigen(crossprod(rbind(Z, gamma * fZ)), symmetric = TRUE)
  U <- U %*% dec$vectors
  Z <- Z %*% dec$vectors
  sigma <- rep(1, ncol(U))

  x <- clone(x, copy.cache = FALSE) # TODO Maybe preserve relevant part of cache?
  .save.oblique.decomposition(x, sigma, U, Z, idx)

  if (!is.null(.get(x, "iossa.groups", allow.null = TRUE))) {
    .set(x, "iossa.groups", .fix.iossa.groups(.get(x, "iossa.groups", allow.null = TRUE), idx))
    .set(x, "iossa.groups.all", .fix.iossa.groups(.get(x, "iossa.groups.all", allow.null = TRUE), idx))
  }

  if (!inherits(x, "ossa")) {
    class(x) <- c("ossa", class(x))
  }

  # Save call info
  x$call <- match.call()

  invisible(x)
}

.fix.iossa.groups <- function(iossa.groups, group) {
  if (length(iossa.groups) == 0) {
    return(list())
  }

  touched <- sapply(iossa.groups, function(g) length(intersect(g, group)) > 0)

  if (any(touched)) {
    iossa.groups <- c(iossa.groups[!touched], union(group, unlist(iossa.groups[touched])))
  }

  iossa.groups
}

decompose.ossa <- function(x, ...) {
  # We can simply implement continuation if we store reference to initial decomposition
  stop("Continuation of decomposition is impossible for ObliqueSSA")
}

.colspan.ossa <- function(x, idx) {
  qr.Q(qr(.U(x)[, idx, drop = FALSE]))
}

.rowspan.ossa <- function(x, idx) {
  qr.Q(qr(calc.v(x, idx)))
}

calc.v.ossa <- function(x, idx, ...) {
  N <- x$length; L <- x$window; K <- N - L + 1
  nV <- nv(x)

  V <- matrix(NA_real_, K, length(idx))
  idx.old <- idx[idx <= nV]
  idx.new <- idx[idx > nV]

  if (length(idx.old) > 0) {
    V[, idx <= nV] <- .V(x)[, idx.old]
  }

  if (length(idx.new) > 0) {
    sigma <- .sigma(x)[idx.new]
    U <- .U(x)[, idx.new, drop = FALSE]
    h <- .get.or.create.hmat(x)
    V[, idx > nV] <- sapply(seq_along(idx.new),
                            function(i) hmatmul(h, U[, i], transposed = TRUE) / sigma[i])
  }

  invisible(V)
}

genwcor <- function(x, groups, basis, ..., cache = TRUE) {
  if (missing(groups)) {
    groups <- x$iossa.groups.all
    if (is.null(groups)) {
      groups <- as.list(seq_len(min(nu(x), nsigma(x))))
    }
  }

  if (missing(basis)) {
    basis <- unlist(.get(x, "iossa.groups.all", allow.null = TRUE))
    if (is.null(basis)) {
      basis <- seq_len(min(nu(x), nsigma(x)))
    }
  }
  basis <- unique(basis)

  # Compute reconstruction.
  F <- reconstruct(x, groups, ..., cache = cache)

  LM <- pseudo.inverse(.U(x)[, basis, drop = FALSE]) # No colspan here!!!!
  RM <- pseudo.inverse(calc.v(x, basis)) # No rowspan here!!!!

  # Compute generalized w-correlations and return
  res <- .gwcor(do.call(cbind, F), L = x$window, LM, RM)
  colnames(res) <- rownames(res) <- names(F)

  res
}
