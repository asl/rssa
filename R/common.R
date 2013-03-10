#   R package for Singular Spectrum Analysis
#   Copyright (c) 2008, 2009 Anton Korobeynikov <asl@math.spbu.ru>
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

# Convenient helpers to operate with data storage

.create.storage <- function(x) {
  attr(x, ".env") <- new.env();
  x;
}

.storage <- function(x)
  attr(x, ".env");

.get <- function(x, name, allow.null = FALSE) {
  ret <- NULL;
  if (!allow.null || .exists(x, name))
    ret <- get(name, envir = .storage(x));

  ret;
}

.set <- function(x, name, value)
  assign(name, value, envir = .storage(x), inherits = FALSE);

.exists <- function(x, name)
  exists(name, envir = .storage(x), inherits = FALSE);

.remove <- function(x, name)
  rm(list = name, envir = .storage(x), inherits = FALSE);

.clone <- function(x, copy.storage = TRUE) {
  # Copy the information body
  obj <- x;

  # Make new storage
  obj <- .create.storage(obj);

  # Copy the contents of data storage
  if (copy.storage) {
    clone.env <- .storage(obj);
    x.env <- .storage(x);
    for (field in ls(envir = x.env, all.names = TRUE)) {
      value <- get(field, envir = x.env, inherits = FALSE);
      attr(value, "..cloned") <- NULL;
      assign(field, value, envir = clone.env, inherits = FALSE);
    }
  }

  .init(obj)
}

.get.series.info <- function(x) {
  if (.exists(x, "series:info"))
    return (.get(x, "series:info"));

  numeric(0);
}

.append.series.info <- function(x, index) {
  .set(x, "series:info",
           union(.get.series.info(x), index));
}

.set.series <- function(x, F, index) {
  name <- paste("series:", index, sep = "");
  .set(x, name, F);
  .append.series.info(x, index);
  index;
}

.get.series <- function(x, index) {
  F <- numeric(prod(x$length));
  for (i in index) {
    name <- paste("series:", i, sep = "");
    F <- F + .get(x, name);
  }
  F;
}

# Generics

# 'ssa' object
clone <- function(x, ...)
  UseMethod("clone");
reconstruct <- function(x, ...)
  UseMethod("reconstruct");
clusterify <- function(x, ...)
  UseMethod("clusterify");
calc.v <- function(x, ...)
  UseMethod("calc.v");

.hankelize.one <- function(x, ...)
  UseMethod(".hankelize.one")
.elseries <- function(x, ...)
  UseMethod(".elseries")
.init <- function(x, ...)
  UseMethod(".init")

# There is decompose() call in stats package, we need to take control over it
decompose <- function(x, ...) UseMethod("decompose");
decompose.default <- function(x, ...) stats::decompose(x, ...);
formals(decompose.default) <- c(formals(decompose.default), alist(... = ));
