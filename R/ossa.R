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

# TODO move this functionality to `clone.ssa()`
.remove.ossa.pssa <- function(class) {
  class <- class[(class != "ossa") & !grepl("^pssa", class)]
}

.clone.with.new.series <- function(x, F, ...) {
  # Get conversion
  conversion <- .inner.fmt.conversion(x)

  s <- clone(x, copy.cache = FALSE, copy.storage = FALSE)
  .set(s, "F", conversion(F))
  class(s) <- .remove.ossa.pssa(class(s))

  s
}

high.rank.rate <- function(F, rank, ssaobj = ssa(F), ...) {
  s <- .clone.with.new.series(ssaobj, F)

  # Decompose for computation leading singular values
  .maybe.continue(s, groups = seq_len(rank), ...)

  1 - sum(.sigma(s)[seq_len(rank)]^2) / wnorm(s)^2
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

# TODO Make emat really polymorphic class. Such wrappings are very ugly
.traj.matrix <- function(x, ...)
  UseMethod(".traj.matrix")

.traj.matrix.1d.ssa <- .traj.matrix.toeplitz.ssa <- function(x, ...)
  .get.or.create.hmat(x, ...)

.traj.matrix.nd.ssa <- function(x, ...)
  .get.or.create.hbhmat(x, ...)

.traj.matrix.mssa <- function(x, ...)
  .get.or.create.mhmat(x, ...)

.traj.matrix.cssa <- function(x, ...)
  .get.or.create.chmat(x, ...)

.traj.matrix.ssa <- function(x, ...)
  stop("`.traj.matrix' is not implemented for this kind of SSA")

.frob.dist.series.lowrank <- function(F, sigma, Y, Z, ssaobj) {
  # This is efficient implementation of |T(F) - Y diag(sigma) Z^T|^2_F
  s <- .clone.with.new.series(ssaobj, F)
  TF <- .traj.matrix(s)

  sZ <- Z * rep(sigma, each = nrow(Z))

  res <- sum((sZ %*% crossprod(Y)) * sZ) + wnorm(s)^2 - 2*sum(crossprod(TF, Y) * sZ)

  # Fix possible numeric error, force dist be nonnegative
  if (res < 0)
    res <- 0

  res
}

.owcor <- function(Fs, LM, RM, ssaobj) {
  mx <- lapply(Fs,
               function(F) {
                 s <- .clone.with.new.series(ssaobj, F)
                 TF <- .traj.matrix(s)
                 as.vector(LM %*% tcrossprod(TF, RM))
              })
  mx <- do.call(cbind, mx)

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

  cat("\nIterative O-SSA result:\n")
  cat("\tConverged:             ", ifelse(x$converged, "yes", "no"), "\n", sep = "")
  cat("\tIterations:            ", x$iter, "\n", sep = "")
  cat("\tInitial mean(tau):     ", format(mean(x$initial.tau), digits = digits), "\n", sep = "")
  cat("\tInitial tau:           ", paste0(format(x$initial.tau, digits = digits), collapse = ", "), "\n", sep = "")
  cat("\tI. O-SSA mean(tau):    ", format(mean(x$tau), digits = digits), "\n", sep = "")
  cat("\tI. O-SSA tau:          ", paste0(format(x$tau, digits = digits), collapse = ", "), "\n", sep = "")
  cat("\tInitial max wcor:      ", format(max.abs.nodiag(x$initial.wcor), digits = digits), "\n", sep = "")
  cat("\tI. O-SSA max wcor:     ", format(max.abs.nodiag(x$wcor), digits = digits), "\n", sep = "")
  cat("\tI. O-SSA max owcor:    ", format(max.abs.nodiag(x$owcor), digits = digits), "\n", sep = "")
  cat("\n")
}

summary.iossa.result <- function(object, digits = max(3, getOption("digits") - 3), ...)
  print.iossa.result(x = object, digits = digits, ...)

print.ossa <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  NextMethod()

  iossa.result <- .get(x, "iossa.result", allow.null = TRUE)
  if (!is.null(iossa.result))
    print(iossa.result, digits = digits)
}

summary.ossa <- function(object, digits = max(3, getOption("digits") - 3), ...)
  print(x = object, digits = digits, ...)

.save.oblique.decomposition <- function(x, nosigma, Y, Z, idx) {
  sigma <- .sigma(x)
  U <- .U(x)
  V <- if (nv(x) < max(idx)) calc.v(x, seq_len(max(idx))) else .V(x)

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

.make.ossa.result <- function(x, Fs, ranks, idx, initial.Fs) {
  IBL <- pseudo.inverse(.U(x)[, idx, drop = FALSE])
  IBR <- pseudo.inverse(calc.v(x, idx))

  tau <- sapply(seq_along(Fs),
                function(i) {
                  high.rank.rate(Fs[[i]], ranks[i], ssaobj = x)
                })

  initial.tau <- sapply(seq_along(initial.Fs),
                        function(i) {
                          high.rank.rate(initial.Fs[[i]], ranks[i], ssaobj = x)
                        })

  names(Fs) <- paste("F", seq_along(Fs), sep = "")
  out <- list(cond = c(Cond(IBL), Cond(IBR)),
              owcor = .owcor(Fs, IBL, IBR, ssaobj = x), #MB just call owcor(x) here?
              wcor = wcor(x, Fs = Fs),
              initial.wcor = wcor(x, Fs = initial.Fs),
              tau = tau,
              initial.tau = initial.tau,
              initial.rec = initial.Fs)
  class(out) <- "iossa.result"

  invisible(out)
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
  V <- if (nv(x) < desired) calc.v(x, idx) else .V(x)[, idx, drop = FALSE]

  # TODO Perform orthogonolize if it's only needed
  dec <- orthogonalize(U, V, sigma, side = "bi")
  sigma <- dec$d; U <- dec$u; V <- dec$v

  if (min(sigma) < eps)
    warning("Decomposition isn't minimal. Some singular values equal to zero")

  list(sigma = sigma, U = U, V = V)
}

iossa <- function(x, nested.groups, ..., tol = 1e-5, kappa = 2,
                  maxiter = 100,
                  norm = function(x) sqrt(mean(x^2)),
                  trace = FALSE,
                  kappa.balance = 0.5) {
  if (inherits(x, "cssa"))
    stop("I-OSSA is not implemented for Complex SSA yet")

  # Get mask
  mask <- .hweights(x) > 0

  if (missing(nested.groups))
    nested.groups <- as.list(1:min(nsigma(x), nu(x)))

  # Continue decomposition, if necessary
  .maybe.continue(x, groups = nested.groups, ...)

  nested.groups <- lapply(nested.groups, unique)
  if (length(unique(unlist(nested.groups))) != sum(sapply(nested.groups, length)))
    stop("Intersected nested groups")

  Fs <- rec <- reconstruct(x, groups = nested.groups) # Get initial approx

  idx <- sort(unique(unlist(nested.groups)))
  triples <- .get.orth.triples(x, idx)
  osigma <- triples$sigma; U <- triples$U; V <- triples$V

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

  converged <- FALSE
  for (iter in seq_len(maxiter)) {
    lrss <- numeric(length(nested.groups))
    triples.list <- list()
    for (i in seq_along(nested.groups)) {
      cur.ssa <- .clone.with.new.series(x, Fs[[i]])

      .maybe.continue(cur.ssa, seq_len(ranks[i]))
      lrss[i] <- sum(cur.ssa$sigma[-seq_len(ranks[i])]^2)
      triples.list[[i]] <- .get.orth.triples(cur.ssa, seq_len(ranks[i]))
    }

    if (trace) cat(sprintf("LRSS(%d): %s\n", iter, paste0(sqrt(lrss), collapse = " ")))

    osigmas <- lapply(triples.list, function(el) el$sigma)
    Us <- lapply(triples.list, function(el) el$U)
    Vs <- lapply(triples.list, function(el) el$V)

    if (!is.null(kappa)) {
      # If we use reordering, force separability
      for (i in seq_along(nested.groups[-length(nested.groups)])) {
        div <- osigmas[[i]][ranks[i]] / osigmas[[i + 1]][1]
        if (div < kappa) {
          mul <- kappa / div
          osigmas[[i + 1]] <- osigmas[[i + 1]] / mul
          Us[[i + 1]] <- Us[[i + 1]] * mul^kappa.balance
          Vs[[i + 1]] <- Vs[[i + 1]] * mul^(1 - kappa.balance)
        }
      }
    }

    basU <- do.call(cbind, Us)
    basV <- do.call(cbind, Vs)

    # Oblique SVD
    dec <- svd2LRsvd(osigma, U, V, basU, basV, need.project = TRUE)
    sigma <- dec$sigma; Y <- dec$Y; Z <- dec$Z

    # Clone is used for hankelization
    # We copy storage for saving quasi-hankel circulant, it is used for stored fft_plan
    # TODO Separate circulants and plans just like in 1d-SSA
    s <- clone(x, copy.cache = FALSE, copy.storage = TRUE)

    # Add `ossa` class for prevent redecomposition
    class(s) <- c("ossa", class(s)[class(s) != "ossa"])

    .set.decomposition(s, sigma = sigma, U = Y, V = Z)
    Fs.new <- reconstruct(s, groups = nested.groups)

    if (trace) {
      svddist <- sapply(seq_along(nested.groups),
                        function(i) .frob.dist.series.lowrank(Fs[[i]],
                                                              sigma[nested.groups[[i]]],
                                                              Y[, nested.groups[[i]], drop = FALSE],
                                                              Z[, nested.groups[[i]], drop = FALSE],
                                                              x))

      nonhank <- sapply(seq_along(nested.groups),
                        function(i) .frob.dist.series.lowrank(Fs.new[[i]],
                                                              sigma[nested.groups[[i]]],
                                                              Y[, nested.groups[[i]], drop = FALSE],
                                                              Z[, nested.groups[[i]], drop = FALSE],
                                                              x))
      cat(sprintf("SVDD(%d): %s\n", iter, paste0(sqrt(svddist), collapse = " ")))
      cat(sprintf("NHNK(%d): %s\n", iter, paste0(sqrt(nonhank), collapse = " ")))
    }

    deltas <- sapply(seq_along(nested.groups), function(i) .series.dist(Fs.new[[i]], Fs[[i]], norm, mask))
    if (max(deltas) < tol) {
      converged <- TRUE
      Fs <- Fs.new
      break
    }

    Fs <- Fs.new
  }

  x <- clone(x, copy.cache = FALSE) # TODO Maybe preserve relevant part of cache?
  .save.oblique.decomposition(x, sigma, Y, Z, idx)

  # Update class for x
  if (!inherits(x, "ossa")) {
    class(x) <- c("ossa", class(x))
  }

  # Save call info
  x$call <- match.call()

  out <- .make.ossa.result(x, Fs, ranks, idx, rec)
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

  # Save nested components
  .set(x, "ossa.set", idx)

  if (!is.null(.decomposition(x, "nPR"))) {
    if (any(idx <= .decomposition(x, "nPR"))) {
      .set.decomposition(x, nPR = 0, nPL = 0)
    } else if (any(idx <= sum(unlist(.decomposition(x, c("nPR", "nPL")))))){
      .set.decomposition(x, nPL = 0)
    }
  }

  invisible(x)
}

fossa <- function(x, nested.groups, FILTER = diff, gamma = 1, normalize = FALSE, ...) {
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
  osigma <- triples$sigma; U <- triples$U; V <- triples$V

  Z <- V * rep(osigma, each = nrow(V))

  Y <- if (normalize) V else Z

  fY <- apply(Y, 2, FILTER)
  dec <- eigen(crossprod(rbind(Y, gamma * fY)), symmetric = TRUE)
  U <- U %*% dec$vectors
  Z <- Z %*% dec$vectors
  sigma <- rep(1, ncol(U))

  x <- clone(x, copy.cache = FALSE) # TODO Maybe preserve relevant part of cache?
  .save.oblique.decomposition(x, sigma, U, Z, idx)

  if (!is.null(.get(x, "iossa.groups", allow.null = TRUE))) {
    .set(x, "iossa.groups", .fix.iossa.groups(.get(x, "iossa.groups", allow.null = TRUE), idx))
    .set(x, "iossa.groups.all", .fix.iossa.groups(.get(x, "iossa.groups.all", allow.null = TRUE), idx))
  }

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

owcor <- function(x, groups, ..., cache = TRUE) {
  # Check class
  stopifnot(inherits(x, "ossa"))

  # Get basis
  basis <- .get(x, "ossa.set")

  if (missing(groups)) {
    groups <- x$iossa.groups
    if (is.null(groups)) {
      groups <- as.list(basis)
    }
  }

  if (!all(unlist(groups) %in% basis))
    stop(sprintf("groups in `groups' must consist of components from current OSSA set (%s)",
                 paste0(basis, collapse = ", ")))

  # Compute reconstruction.
  Fs <- reconstruct(x, groups, ..., cache = cache)

  LM <- pseudo.inverse(.U(x)[, basis, drop = FALSE]) # No colspan here!!!!
  RM <- pseudo.inverse(calc.v(x, basis)) # No rowspan here!!!!

  # Compute oblique w-correlations and return
  res <- .owcor(Fs, LM, RM, ssaobj = x)
  colnames(res) <- rownames(res) <- names(Fs)

  res
}
