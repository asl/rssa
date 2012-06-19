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

.create.storage <- function(this) {
  attr(this, ".env") <- new.env();
  this;
}

.storage <- function(this)
  attr(this, ".env");

.get <- function(this, name, allow.null = FALSE) {
  ret <- NULL;
  if (!allow.null || .exists(this, name))
    ret <- get(name, envir = .storage(this));

  ret;
}

.set <- function(this, name, value)
  assign(name, value, envir = .storage(this), inherits = FALSE);

.exists <- function(this, name)
  exists(name, envir = .storage(this), inherits = FALSE);

.remove <- function(this, name)
  rm(list = name, envir = .storage(this), inherits = FALSE);

.clone <- function(this, copy.storage = TRUE) {
  # Copy the information body
  obj <- this;

  # Make new storage
  obj <- .create.storage(obj);

  # Copy the contents of data storage
  if (copy.storage) {
    clone.env <- .storage(obj);
    this.env <- .storage(this);
    for (field in ls(envir = this.env, all.names = TRUE)) {
      value <- get(field, envir = this.env, inherits = FALSE);
      attr(value, "..cloned") <- NULL;
      assign(field, value, envir = clone.env, inherits = FALSE);
    }
  }

  obj;
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
  F <- numeric(prod(this$length));
  for (i in index) {
    name <- paste("series:", i, sep = "");
    F <- F + .get(this, name);
  }
  F;
}

# Generics

# 'ssa' object
clone <- function(this, ...)
  UseMethod("clone");
reconstruct <- function(this, ...)
  UseMethod("reconstruct");
nu <- function(this, ...)
  UseMethod("nu");
nv <- function(this, ...)
  UseMethod("nv");
nlambda <- function(this, ...)
  UseMethod("nlambda");
precache <- function(this, ...)
  UseMethod("precache");
cleanup <- function(this, ...)
  UseMethod("cleanup");
clusterify <- function(this, ...)
  UseMethod("clusterify");
calc.v <- function(this, ...)
  UseMethod("calc.v");

.object.size <- function(this, ...)
  UseMethod(".object.size")

.hankelize.one <- function(this, ...)
  UseMethod(".hankelize.one")

# There is decompose() call in stats package, we need to take control over it
decompose <- function(x, ...) UseMethod("decompose");
decompose.default <- stats::decompose;
formals(decompose.default) <- c(formals(decompose.default), alist(... = ));
