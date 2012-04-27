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
  N <- length(F);
  r <- length(lrf);

  # Sanity check of inputs
  if (r > N)
    stop("Wrong length of LRF");

  # Run the actual LRF
  F <- c(F, rep(NA, len));
  
  # TODO Rewrite this function on C
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
    lf <- lrf(this, group);

    # Calculate the forecasted values
    out[[i]] <- apply.lrf(if (identical(base, "reconstructed")) r[[i]] else .get(this, "F"),
                          lf, len);
    # FIXME: try to fixup the attributes
    # FIXME: ASH It's impossible --- different lengths =(
  }

  names(out) <- paste(sep = "", "F", 1:length(groups));

  # Forecasted series can be pretty huge...
  invisible(out);
}



"vforecast.1d-ssa" <- function(this, groups, len = 1,
                               ...,
                               cache = TRUE) {
  L <- this$window;
  K <- this$length - L + 1;
  N <- K + L - 1 + len + L - 1;
  N.res <- K + L - 1 + len;

  if (missing(groups))
    groups <- as.list(1:min(nlambda(this), nu(this)));

  # Determine the upper bound of desired eigentriples
  desired <- max(unlist(groups));

  # Continue decomposition, if necessary
  if (desired > min(nlambda(this), nu(this)))
    decompose(this, ..., neig = desired);

  lambda <- .get(this, "lambda");
  U <- .get(this, "U");

  if (nv(this) >= desired) {
    V <- .get(this, "V");
  } else {
    if (cache){
      V <- cbind(.get(this, "V", allow.null = TRUE),
                 calc.v(this, idx = (nv(this) + 1):desired));
      .set(this, "V", V);
    } else {
      V <- NULL;
    }
  }
  
  # Make hankel matrix for fast hankelization (we use it for plan) 
  h <- new.hmat(double(N), L);
  
  out <- list();
  for (i in seq_along(groups)) {
    group <- unique(groups[[i]]);
    
    Uet <- U[, group, drop = FALSE];
    Vet <- if (is.null(V)) calc.v(this, idx = group) else V[, group, drop = FALSE];
    
    Z <- rbind(t(lambda[group] * t(Vet)), matrix(NA, len + L - 1, length(group)));
    
    U.head <- Uet[-L, , drop = FALSE];
    U.tail <- Uet[-1, , drop = FALSE];
    Pi <- Uet[L, ];
    tUhUt <- t(U.head) %*% U.tail;
    P <- tUhUt + 1 / (1 - sum(Pi^2)) * Pi %*% (t(Pi) %*% tUhUt);
    
    for(j in (K + 1):(K + len + L - 1)) {
      Z[j, ] <- P %*% Z[j - 1, ];
    }
     
    res <- double(N);
    for(j in seq_along(group)) {
      res <- res + .hankelize.one.hankel(Uet[ , j], Z[ , j], h);
    }
    
    out[[i]] <- res[1:N.res];
  }
  
  names(out) <- paste(sep = "", "F", 1:length(groups));
  
  # Forecasted series can be pretty huge...
  invisible(out);
}

"lrf.toeplitz-ssa" <- get("lrf.1d-ssa");
"vforecast.toeplitz-ssa" <- get("vforecast.1d-ssa");
"rforecast.toeplitz-ssa" <- get("rforecast.1d-ssa");

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
