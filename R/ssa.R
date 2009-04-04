#   R package for Singular Spectrum Analysis
#   Copyright (c) 2008 Anton Korobeynikov <asl@math.spbu.ru>
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

new.ssa <- function(x,
                    L = (N - 1) %/% 2,
                    ..., 
                    method = c("hankel", "toeplitz"),
                    force.decompose = TRUE) {
  method <- match.arg(method);
  N <- length(x);

  # Create information body
  this <- list(length = N,
               window = L,
               call = match.call(),
               series = deparse(substitute(x)),
               method = method);

  # Create data storage
  this <- .create.storage(this);

  # Save series
  .set(this, "F", x);
  
  # Make this S3 object
  class(this) <- "ssa";

  # Decompose, if necessary
  if (force.decompose)
    decompose(this, ...);

  this;
}

.decompose.ssa.toeplitz <- function(this, ...) {
  stop("Unimplemented!")
}

decompose.ssa <- function(this, ...) {
  method <- this$method;

  if (identical(method, "hankel")) {
    .decompose.ssa.hankel(this, ...);
  } else if (identical(method, "toeplitz")) {
    .decompose.ssa.toeplitz(this, ...);
  } else {
    stop("Unknown method in SSA")
  }
}

precache.ssa <- function(this, n, ...) {
  if (missing(n)) {
    warning("Amount of sub-series missed, precaching EVERYTHING",
            immediate. = TRUE);
    n <- nlambda(this);
  }

  # Calculate numbers of sub-series to be calculated
  info <- .get.series.info(this);
  new <- setdiff(1:n, info);

  # We're supporting only 'full' data for now
  U <- .get(this, "U");
  V <- .get(this, "V");
  lambda <- .get(this, "lambda");

  F <- .hankelize.multi(U[,new], V[,new]);

  # Return numbers of sub-series cached
  invisible(sapply(new,
                   function(i) {
                     .set.series(this,
                                 lambda[i] * F[,i],
                                 i)}));
}

.get.series.info <- function(this) {
  if (.exists(this, "series:info"))
    return (.get(this, "series:info"));

  numeric(0);
}

.append.series.info <- function(this, index) {
  .set(this, "series:info",
           union(.get.series.info(this), index));
}

.set.series <- function(this, F, index) {
  name <- paste("series:", index, sep = "");
  .set(this, name, F);
  .append.series.info(this, index);
  index;
}

.get.series <- function(this, index) {
  F <- numeric(this$length);
  for (i in index) {
    name <- paste("series:", i, sep = "");
    F <- F + .get(this, name);
  }
  F;
}

cleanup.ssa <- function(this, ...) {
  .remove(this, ls(.storage(this), pattern = "series:"));
  gc();
}

reconstruct.ssa <- function(this, groups, ..., cache = TRUE) {
  out <- list();
  nu <- nu(this); nv <- nv(this);

  if (missing(groups))
    groups <- as.list(1:nlambda(this));

  # We're supporting only 'full' data for now
  U <- .get(this, "U");
  V <- .get(this, "V");
  lambda <- .get(this, "lambda");

  # Grab indices of pre-cached values
  info <- .get.series.info(this);

  # Do actual reconstruction
  for (i in seq_along(groups)) {
    group <- groups[[i]];
    new <- setdiff(group, info);
    cached <- intersect(group, info);

    if (length(new) == 0) {
      # Nothing to compute, just create zero output
      out[[i]] <- numeric(this$length);
    } else if (length(new) == 1) {
      # Special case for rank one reconstruction
      out[[i]] <- lambda[new] * .hankelize.one(U[, new], V[, new]);

      # Cache the reconstructed series, if this was requested
      if (cache)
        .set.series(this, out[[i]], new);
    } else {
      out[[i]] <- hankel(U[, new] %*%
                         diag(lambda[new], nrow = length(group)) %*%
                         t(V[, new]));
    }
    # Add pre-cached series
    out[[i]] <- out[[i]] + .get.series(this, cached);
  }

  names(out) <- paste("F", 1:length(groups), sep="");

  # Reconstructed series can be pretty huge...
  invisible(out);
}

nu.ssa <- function(this, ...) {
  ifelse(.exists(this, "U"), ncol(.get(this, "U")), 0);
}

nv.ssa <- function(this, ...) {
  ifelse(.exists(this, "V"), ncol(.get(this, "V")), 0);
}

nlambda.ssa <- function(this, ...) {
  ifelse(.exists(this, "lambda"), length(.get(this, "lambda")), 0);
}

clone.ssa <- function(this, copy.cache = TRUE, ...) {
  obj <- .clone(this);
  if (copy.cache == FALSE)
    cleanup(obj);

  obj;
}

clusterify.ssa <- function(this, groups, nclust = length(groups) / 2,
                           ...,
                           type = c("wcor"), cache = TRUE) {
  type <- match.arg(type)

  if (missing(groups))
    groups <- as.list(1:nlambda(this));

  if (identical(type, "wcor")) {
    w <- wcor(this, groups = groups, ..., cache = cache);
    g <- clusterify(w, nclust = nclust, ...);
    out <- lapply(g, function(idx) unlist(groups[idx]));
  } else {
    stop("Unsupported clusterification method!");
  }
  out;
}

'$.ssa' <- function(this, name) {
  if (ind <- charmatch(name, names(this), nomatch = 0))
    return (this[[ind]]);

  if (.exists(this, name))
    return (.get(this, name));

  NULL;
}

#.F <- function(x) exp(-.01 * x)*cos(x/100) + 0.05*rnorm(length(x));
# F <- .F(1:5000);
