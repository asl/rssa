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

high.rank.rate <- function(F, L, r) {
  ss <- ssa(F, L, svd.method = "propack", neig = r + 1)

  1 - sum(ss$sigma[(1:r)]^2) / wnorm(F, L)^2
}

pseudo.inverse <- function(A) {
  # Moore-Penrose pseudo-inverse MB Use MASS::ginv()? But ginv uses SVD instead of QR.
  qrA <- qr(A)
  solve(qr.R(qrA), t(qr.Q(qrA)))
}

orthogonalize <- function(Y, Z, sigma, side = c("bi", "left", "right"), normalize = TRUE) {
  # Orthtogonalize. From some minimal decomposition to SVD or one side othodecomposition
  if (missing(Z)) {
    # Just QR-orthogonalization
    return(qr.Q(qr(Y)))
  }

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

gwcor <- function(X, L, LM, RM) {
  # TODO Optimize it
  mx <- apply(X, 2, function(x) as.vector(LM %*% hankel(x, L) %*% t(RM)))

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

plot.ossa.result <- function(x, ...) {
  if (is.null(x$converged))
    x$converged <- TRUE

  if (is.null(x$iter))
    x$iter <- NA

  main <- sprintf("M = %3.2f iter = %d hrr: %4e -> %4e converged = %s", x$Cond, x$iter,
      x$mean.initial.hrr, x$mean.hrr,
      ifelse(x$converged, "T", "F"))

  col <- c("grey", rainbow(length(x$Fs)))
  matplot(do.call(cbind, c(list(x$F), x$Fs)), type = "l", main = main, lty = "solid",
      col = col, ylab = "F", xlab = "x", ...)
  legend("topright", legend = c("Initial", names(x$Fs)),
      col = col, lty = "solid")
}

summary.ossa.result <- function(x, ...) {
  cat("\nOSSA result:\n")
  cat("\tConverged:             ", ifelse(x$converged, "yes", "no"), "\n", sep = "")
  cat("\tCondition number:      ", x$Cond, "\n", sep = "")
  cat("\tIterations:            ", x$iter, "\n", sep = "")
  cat("\tInitial mean(tau):     ", x$mean.initial.hrr, "\n", sep = "")
  cat("\tInitial tau:           ", paste0(x$initial.hrr, collapse = ", "), "\n", sep = "")
  cat("\tOSSA mean(tau):        ", x$mean.hrr, "\n", sep = "")
  cat("\tOSSA tau:              ", paste0(x$hrr, collapse = ", "), "\n", sep = "")

  wcor <- x$wcor
  diag(wcor) <- 0
  max.wcor <- wcor[which.max(abs(wcor))]

  initial.wcor <- x$initial.wcor
  diag(initial.wcor) <- 0
  max.initial.wcor <- initial.wcor[which.max(abs(initial.wcor))]

  cat("\tInitial max wcor:      ", max.initial.wcor, "\n", sep = "")
  cat("\tOSSA max wcor:         ", max.wcor, "\n", sep = "")
  cat("\n")
}

.save.oblique.decomposition <- function(x, nosigma, Y, Z, idx, IBL, IBR) {
  U <- .get(x, "U")
  sigma <- .get(x, "sigma")
  V <- if (nv(x) < max(idx)) calc.v.ossa(x, seq_len(max(idx))) else .get(x, "V")

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

  .set(x, "sigma", sigma)
  .set(x, "U", U)
  .set(x, "V", V)

  .set(x, "L", IBL)
  .set(x, "R", IBR)

  # Save indexes of oblique triples (we need orthogonalization for them, in forecast, in calc.v)
  .set(x, "oblique.triples.idx", union(.get(x, "oblique.triples.idx", allow.null = TRUE), idx))
}

.make.ossa.result <- function(x, Fs, ranks, IBL, IBR, initial.Fs) {
  L <- x$window

  hrr <- sapply(seq_along(Fs),
                function(i) {
                  high.rank.rate(Fs[[i]], L, ranks[i])
                })

  initial.hrr <- sapply(seq_along(initial.Fs),
                        function(i) {
                          high.rank.rate(initial.Fs[[i]], L, ranks[i])
                        })

  names(Fs) <- paste("F", seq_along(Fs), sep = "")
  out <- list(
      F = .get(x, "F"),
      Fs = Fs,
      Cond = Cond(IBL) * Cond(IBR),
      gwcor = gwcor(do.call(cbind, Fs), L, IBL, IBR),
      wcor = wcor(do.call(cbind, Fs), L),
      initial.wcor = wcor(do.call(cbind, initial.Fs), L),
      hrr = hrr, mean.hrr = mean(hrr),
      initial.hrr = initial.hrr, mean.initial.hrr = mean(initial.hrr),
      initial.rec = initial.Fs)
  class(out) <- "ossa.result"

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

.get.orth.triples <- function(x, idx) {
  # Determine the upper bound of desired eigentriples
  desired <- max(idx)

  # Continue decomposition, if necessary
  if (desired > min(nsigma(x), nu(x)))
    decompose(x, neig = desired)

  lambda <- .get(x, "sigma")[idx]
  U <- .get(x, "U")[, idx, drop = FALSE]
  V <- if (nv(x) < desired) calc.v.ossa(x, idx) else .get(x, "V")[, idx, drop = FALSE]

  # # Orthogonolize decomposition if needed
  # oblique.triples.idx <- .get(x, "oblique.triples.idx", allow.null = TRUE)
  # if (any(idx %in% oblique.triples.idx) || inherits(x, "toeplitz.ssa")) {
  #   dec <- orthogonalize(U, V, lambda, side = "bi")
  #   lambda <- dec$d; U <- dec$u; V <- dec$v
  # }

  dec <- orthogonalize(U, V, lambda, side = "bi")
  lambda <- dec$d; U <- dec$u; V <- dec$v
  list(lambda = lambda, U = U, V = V)
}

assa <- function(x, groups, ..., tol = 1e-5, kappa = 1.2,
                 maxIter = 1000,
                 initial.approx,
                 norm = function(x) {
                   sqrt(mean(x^2))
                 }, trace = FALSE) {
  x <- clone(x, copy.cache = FALSE)
  N <- x$length; L <- x$window; K <- N - L + 1

  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))

  rec <- reconstruct(x, groups = groups)

  if (missing(initial.approx))
    initial.approx <- rec

  Fs <- initial.approx

  ind <- sort(unique(unlist(groups)))
  triples <- .get.orth.triples(x, ind)
  lambda <- triples$lambda; U <- triples$U; V <- triples$V

  # Replace (in groups) ET's numbers for their ranks (order numbers) in set of signal ETs
  groups <- lapply(groups, function(group) match(group, ind))
  ranks <- sapply(groups, length)

  # If we use reordering, reorder components
  if (!is.null(kappa)) {
    cumranks <- cumsum(ranks)
    for (i in seq_along(groups)) {
      groups[[i]] <- seq(to = cumranks[i], length.out = ranks[i])
    }
  }

  # Grab the FFT plan
  fft.plan <- .get.or.create.fft.plan(x)

  converged <- FALSE
  for (iter in seq_len(maxIter)) {
    lrss <- numeric(length(groups))
    triples.list <- list()
    for (i in seq_along(groups)) {
      cur.ssa <- ssa(Fs[[i]], L, svd.method = "propack", neig = 2 * ranks[i])
      lrss[i] <- sum(cur.ssa$sigma[-seq_len(ranks[i])]^2)
      triples.list[[i]] <- .get.orth.triples(cur.ssa, 1:ranks[i])
    }

    if (trace) cat(sprintf("LRSS(%d): %s\n", iter, paste0(sqrt(lrss), collapse = " ")))

    lambdas <- lapply(triples.list, function(el) el$lambda)
    Us <- lapply(triples.list, function(el) el$U)
    Vs <- lapply(triples.list, function(el) el$V)

    if (!is.null(kappa)) {
      # If we use reordering, force separability
      for (i in seq_along(groups[-length(groups)])) {
        div <- lambdas[[i]][ranks[i]] / lambdas[[i + 1]][1]
        if (div < kappa) {
          mul <- kappa / div
          lambdas[[i + 1]] <- lambdas[[i + 1]] / mul
          Us[[i + 1]] <- Us[[i + 1]] * sqrt(mul)
          Vs[[i + 1]] <- Vs[[i + 1]] * sqrt(mul)
        }
      }
    }

    basU <- do.call(cbind, Us)
    basV <- do.call(cbind, Vs)

    # Oblique SVD
    dec <- svd2LRsvd(lambda, U, V, basU, basV, need.project = TRUE)
    sigma <- dec$sigma; Y <- dec$Y; Z <- dec$Z

    elem.components <- .hankelize.multi.default(t(sigma * t(Y)), Z, fft.plan)
    Fs.new <- .component.grouping(elem.components, groups)

    if (trace) {
      Xs <- lapply(groups,
                   function(group) Y[, group, drop = FALSE] %*% (sigma[group] * t(Z[, group, drop = FALSE])))
      svddist <- sapply(seq_along(Xs),
                        function(i) sum((Xs[[i]] - hankel(Fs[[i]], L)) ^ 2))
      cat(sprintf("SVDD(%d): %s\n", iter, paste0(sqrt(svddist), collapse = " ")))

      nonhankeleness <- sapply(seq_along(Xs),
                               function(i) sum((Xs[[i]] - hankel(Fs.new[[i]], L)) ^ 2))
      cat(sprintf("NHNK(%d): %s\n", iter, paste0(sqrt(nonhankeleness), collapse = " ")))
    }

    deltas <- sapply(seq_along(groups), function(i) norm(Fs.new[[i]] - Fs[[i]]))
    if (max(deltas) < tol) {
      converged <- TRUE
      Fs <- Fs.new
      break
    }

    Fs <- Fs.new
  }

  IBL <- pseudo.inverse(Y)
  IBR <- pseudo.inverse(Z)

  .save.oblique.decomposition(x, sigma, Y, Z, ind, IBL, IBR)

  # Update class for x
  if (!inherits(x, "ossa")) {
    class(x) <- c("ossa", class(x))
  }

  out <- .make.ossa.result(x, Fs, ranks, IBL, IBR, rec)
  out[c("iter", "converged", "kappa", "maxIter", "tol", "call")] <-
      list(iter, converged, kappa, maxIter, tol, match.call())

  class(out) <- c(class(out), "assa.result")

  .set(x, "ossa.result", out)
  invisible(x)
}

fssa <- function(x, groups, FILTER = diff, gamma = 0, ...) {
  x <- clone(x, copy.cache = FALSE)

  if (!is.function(FILTER)) {
    FILTER.coeffs <- FILTER
    FILTER <- function(x) {
      out <- filter(x, FILTER.coeffs)
      out[!is.na(out)]
    }
  }

  N <- x$length; L <- x$window; K <- N - L + 1

  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))

  rec <- reconstruct(x, groups = groups)

  ind <- sort(unique(unlist(groups)))
  triples <- .get.orth.triples(x, ind)
  lambda <- triples$lambda; U <- triples$U; V <- triples$V

  # Replace (in groups) ET's numbers for their ranks (order numbers) in set of signal ETs
  groups <- lapply(groups, function(group) match(group, ind))
  ranks <- sapply(groups, length)

  # Grab the FFT plan
  fft.plan <- .get.or.create.fft.plan(x)

  fV <- apply(V, 2, FILTER)
  # We transform right basis as that it became orthogonal in special inner product
  qrR <- qr.R(qr(rbind(V, gamma * fV)))

  general.way <- TRUE
  if (general.way) {
    # GENERALLY
    pbasV <- V %*% solve(qrR)

    # We can do it more optimal here. But this is more general
    dec <- svd2LRsvd(lambda, U, V, NULL, pbasV, need.project = FALSE)
    Y <- dec$Y; Z <- dec$Z; sigma <- dec$sigma
  } else {
    # FASTER
    dec <- svd(lambda * t(qrR))
    Y <- U %*% dec$u
    Z <- V %*% solve(qrR, dec$v)
    sigma <- dec$d
  }

  elem.components <- .hankelize.multi.default(t(sigma * t(Y)), Z, fft.plan)
  Fs <- .component.grouping(elem.components, groups)

  IBL <- t(Y)
  IBR <- pseudo.inverse(Z)

  .save.oblique.decomposition(x, sigma, Y, Z, ind, IBL, IBR)

  # Update class for x
  if (!inherits(x, "ossa")) {
    class(x) <- c("ossa", class(x))
  }

  out <- .make.ossa.result(x, Fs, ranks, IBL, IBR, rec)
  out[c("filter", "gamma", "call")] <-
      list(FILTER, gamma, match.call())

  class(out) <- c(class(out), "fssa.result")

  .set(x, "ossa.result", out)
  invisible(x)
}

decompose.ossa <- function(x, ...) {
  stop("Continuation of decomposition is impossible for ObliqueSSA. Yet.")
}

calc.v.ossa <- function(x, idx, ...) {
  if (length(idx) == 0)
    return(matrix(0, x$length - x$window + 1, 0))

  oblique.triples.idx <- .get(x, "oblique.triples.idx", allow.null = TRUE)
  oblique.idx <- intersect(idx, oblique.triples.idx)
  orth.idx <- setdiff(idx, oblique.triples.idx)

  sigma <- .get(x, "sigma")[idx]
  U <- .get(x, "U")[ , idx, drop = FALSE]

  h <- .get.or.create.hmat(x)

  if (length(oblique.idx) > 0) {
    iU <- t(pseudo.inverse(.get(x, "U")[, oblique.triples.idx, drop = FALSE]))
    U[, match(oblique.idx, idx)] <- iU[, match(oblique.idx, oblique.triples.idx)]
  }

  invisible(sapply(seq_along(idx),
          function(i) hmatmul(h, U[, i], transposed = TRUE) / sigma[i]))
}

# calc.v.ossa.prototype <- function(x, idx) {
#   cached.idx <- idx[idx <= nv(x)]
#   new.idx <- idx[idx > nv (x)]
# 
#   V <- matrix(NA, x$length - x$window + 1, length(idx))
#   V[, idx <= nv(x)] <- x$V[, cached.idx]
# 
#   if (length(new.idx) > 0) {
#     sigma <- .get(x, "sigma")[new.idx]
#     U <- .get(x, "U")[ , new.idx, drop = FALSE]
# 
#     h <- .get.or.create.hmat(x)
# 
#     V[, idx > nv(x)] <- sapply(seq_along(idx),
#                                function(i) hmatmul(h, U[, i], transposed = TRUE) / sigma[i])
#   }
# 
#   invisible(V)
# }

mcor <- function(x, groups) {
  L <- x$window; N <- x$length; K <- N - L + 1
  if (missing(groups))
    groups <- as.list(1:nlambda(x))

  rank <- max(unlist(groups))
  U <- x$U[, seq_len(rank), drop = FALSE]
  V <- calc.v(x, seq_len(rank))

  cov.el <- crossprod(U) * crossprod(V)

  cov <- outer(seq_along(groups), seq_along(groups),
               Vectorize(function(i, j) sum(cov.el[groups[[i]], groups[[j]]])))

  # Convert to correlations
  cor <- cov2cor(cov)

  # Fix possible numeric error
  cor[cor > 1] <- 1; cor[cor < -1] <- -1

  # Add class
  class(cor) <- "wcor.matrix"

  # Set names
  colnames(cor) <- rownames(cor) <- paste("F", seq_along(groups), sep = "")

  # Return
  cor
}

wcor.ossa <- function(x, groups, type = c("vanilla", "generalized", "matrix"), ..., cache = TRUE) {
  L <- x$window; N <- x$length; K <- N - L + 1
  if (missing(groups))
    groups <- as.list(1:nlambda(x))

  # Compute reconstruction.
  F <- reconstruct(x, groups, ..., cache = cache)
  mx <- matrix(unlist(F), nrow = N, ncol = length(groups))
  colnames(mx) <- names(F)

  type <- match.arg(type)

  if (identical(type, "generalized")) {
    # Get inner product matrices
    L.mx <- .get(x, "L", allow.null = TRUE)
    if (is.null(L.mx)) {
      # Use vanilla inner product
      # This is not optimal, big matrix
      L.mx <- diag(nrow = L, ncol = L)
    }

    R.mx <- .get(x, "R", allow.null = TRUE)
    if (is.null(R.mx)) {
      # Use vanilla inner product
      # This is not optimal, big matrix
      R.mx <- diag(nrow = K, ncol = K)
    }

    # Finally, compute generalized w-correlations and return
    res <- gwcor(mx, L = L, L.mx, R.mx)
    colnames(res) <- rownames(res) <- paste("F", seq_along(groups), sep = "")

    res
  } else if (identical(type, "vanilla")) {
    # Compute standard w-correlations
    res <- wcor(mx, L = L)
    colnames(res) <- rownames(res) <- paste("F", seq_along(groups), sep = "")

    res
  } else if (identical(type, "matrix")) {
    mcor(x, groups)
  } else {
    stop("Unsupported type of wcor!")
  }
}
