#   R package for Singular Spectrum Analysis
#   Copyright (c) 2012 Alexander Shlemov <shlemovalex@gmail.com>
#   Copyright (c) 2012 Anton Korobeynikov <asl@math.spbu.ru>
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

basis2lrf <- function(U) {
  N <- nrow(U);
  lpf <- U %*% t(U[N, , drop = FALSE]);

  (lpf[-N]) / (1 - lpf[N]);
}

"lrf.1d-ssa" <- function(this, group, ...) {
  # Determine the upper bound of desired eigentriples
  desired <- max(group);

  # Continue decomposition, if necessary
  if (desired > nu(this))
    decompose(this, ..., neig = desired);

  U <- .get(this, "U")[, group, drop = FALSE];

  res <- basis2lrf(U);
  class(res) <- "lrf";

  res;
}

companion.matrix.lrf <- function(x) {
  n <- length(x) + 1;
  res <- matrix(0, n, n);
  res[-n, n] <- x;
  res[seq(from = 2, by = n + 1, length.out = n - 1)] <- 1;

  res;
}

roots.lrf <- function(x) {
  # polyroot(c(-x, 1));
  # Much more complicated but much more stable
  eigen(companion.matrix.lrf(x), only.values = TRUE)$values;
}

plot.lrf <- function(x, ..., raw = FALSE) {
  r <- roots(x);
  if (raw) {
    plot(r, ...);
  } else {
    xlim <- range(c(Re(r), +1, -1));
    ylim <- range(c(Im(r), +1, -1));

    plot(r, ...,
         xlim = xlim, ylim = ylim,
         main = "Roots of Linear Recurrence Formula",
         xlab = "Real Part",
         ylab = "Imaginary Part",
         asp = 1);
    symbols(0, 0, circles = 1, add = TRUE, inches = FALSE);
  }
}

apply.lrf <- function(F, lrf, len = 1) {
  # TODO Rewrite this function on C
  N <- length(F);
  r <- length(lrf);

  # Sanity check of inputs
  if (r > N)
    stop("Wrong length of LRF");

  # Run the actual LRF
  F <- c(F, rep(NA, len));

  for (i in 1:len)
    F[N+i] <- sum(F[(N+i-r) : (N+i-1)]*lrf);

  F;
}

"rforecast.1d-ssa" <- function(this, groups, len = 1,
                               base = c("reconstructed", "original"),
                               ..., cache = TRUE) {
  base <- match.arg(base);
  if (missing(groups))
    groups <- as.list(1:min(nlambda(this), nu(this)));

  # Determine the upper bound of desired eigentriples
  desired <- max(unlist(groups));

  # Continue decomposition, if necessary
  if (desired > min(nlambda(this), nu(this)))
    decompose(this, ..., neig = desired);

  # Grab the reconstructed series if we're basing on them
  if (identical(base, "reconstructed"))
    r <- reconstruct(this, groups = groups, ..., cache = cache);

  out <- list();
  for (i in seq_along(groups)) {
    group <- groups[[i]];

    # Calculate the LRF corresponding to group
    lrf <- lrf(this, group)

    # Calculate the forecasted values
    out[[i]] <- apply.lrf(if (identical(base, "reconstructed")) r[[i]] else .get(this, "F"),
                          lrf, len)
    # FIXME: try to fixup the attributes
  }

  names(out) <- paste(sep = "", "F", 1:length(groups));

  out
}


"vforecast.1d-ssa" <- function(this, groups, len = 1,
                               ...) {
  L <- this$window
  K <- this$length - L + 1;
  L.s <- min(L, K);
  N <- K + L - 1 + len + L - 1;
  N.res <- K + L - 1 + len;

  dv <- c(1:(L.s-1), rep(L.s, N-2*L.s+2), (L.s-1):1);

  convolve.open <- function(F, G) {
    NF <- length(F)
    NG <- length(G)

    NN <- nextn(NF+NG-1)
    ZFZ <- c(rep(0, NG-1), F, rep(0, NN-(NF+NG-1)));
    GZ <- c(G, rep(0, NN-NG));

    res <- fft(fft(ZFZ)*Conj(fft(GZ)), inverse = TRUE)/NN;

    Re(res)[1:(NF+NG-1)];
  }

  if (missing(groups))
    groups <- as.list(1:min(nlambda(this), nu(this)))

  # Determine the upper bound of desired eigentriples
  desired <- max(unlist(groups));

  # Continue decomposition, if necessary
  if (desired > min(nlambda(this), nu(this)))
    decompose(this, ..., neig = desired);

  lambda <- .get(this, "lambda");
  U <- .get(this, "U");

  V <- if (nv(this) >= desired) .get(this, "V") else NULL

  out <- list()
  for (i in seq_along(groups)) {
    ET <- unique(groups[[i]])

    Vl <- matrix(NA, N.res, length(ET));
    Uet <- U[ , ET, drop=FALSE];

    Vl[1:K, ] <- (if (is.null(V))
                    calc.v(this, idx = ET)
                  else V[ , ET]) %*% diag(lambda[ET], nrow = length(ET))

    U.head <- Uet[-L, , drop=FALSE];
    U.tail <- Uet[-1, , drop=FALSE];
    P <- solve(t(U.head) %*% U.head, t(U.head) %*% U.tail);

    for (j in (K+1):(K+len+L-1)) {
      Vl[j, ] <- P %*% Vl[j-1, ];
    }

    res <- rep(0, N);
    for (j in 1:length(ET)) {
      res <- res + convolve.open(Vl[ , j], rev(Uet[ , j]));
    }

    out[[i]] <- (res/dv)[1:N.res];
  };

  names(out) <- paste(sep = "", "F", 1:length(groups));

  out
}

"lrf.toeplitz-ssa" <- `lrf.1d-ssa`;
"vforecast.toeplitz-ssa" <- `vforecast.1d-ssa`;
"rforecast.toeplitz-ssa" <- `rforecast.1d-ssa`;

rforecast.ssa <- function(x, groups, len = 1,
                          base = c("reconstructed", "original"),
                          ..., cache = TRUE) {
  stop("generic recurrent forecast not implemented yet!")
}

lrf.ssa <- function(x, group) {
  stop("generic LRF calculation not implemented yet!")
}

vforecast.ssa <- function(x, groups, len = 1,
                          ...) {
  stop("generic vector forecast not implemented yet!")
}

lrf <- function(this, ...)
  UseMethod("lrf")
roots <- function(x)
  UseMethod("roots")
rforecast <- function(this, ...)
  UseMethod("rforecast")
vforecast <- function(this, ...)
  UseMethod("vforecast")
