#   R package for Singular Spectrum Analysis
#   Copyright (c) 2015-2016 Anton Korobeynikov <anton at korobeynikov dot info>
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

.capabilities <- list()

.register.capability <- function(name, fun,
                                 pred = function(...) TRUE,
                                 alias) {
  if (missing(alias))
    alias <- fun
  .capabilities[[alias]] <<- list(fun = fun, pred = pred, name = name)
}

capable <- function(x, capname) {
  cap.entry <- .capabilities[[capname]]
  stopifnot(!is.null(cap.entry))
  any(sapply(class(x), .check.caps, cap = cap.entry, x = x))
}

.check.caps <- function(cap, klass, x)
  !is.null(getS3method(f = cap$fun, class = klass, optional = TRUE)) && cap$pred(x)

ssa.capabilities <- function(x) {
  res <- list()
  for (klass in class(x))
    res[[klass]] <- sapply(.capabilities, .check.caps, klass = klass, x = x)
  res <- apply(simplify2array(res), 1, any)
  names(res) <- sapply(.capabilities[names(res)], function(x) x$name)
  res
}

.register.capability("Decomposition", "decompose")
.caps.continue <- function(x) {
  ## Only nu-trlan SVD method is capable of continuation
  identical(x$svd.method, "nutrlan")
}
.register.capability("Continuation of a decomposition", "decompose", .caps.continue, "decompose.continue")

.register.capability("Reconstruction", "reconstruct")
.register.capability("Plotting", "plot")

.caps.vforecast <- function(x) {
  ## No forecast in shaped and circular case
  !is.shaped(x) && !x$circular
}

.caps.rforecast <- function(x) {
  ## No forecast in shaped case
  !is.shaped(x) || x$circular
}
.register.capability("Recurrent forecast", "rforecast", .caps.rforecast)
.register.capability("Vector forecast", "vforecast", .caps.vforecast)

.caps.gapfill <- function(x) {
  ## Gapfilling should always start from shaped object
  is.shaped(x)
}
.register.capability("Gapfilling via forecast", "gapfill", .caps.gapfill)
.register.capability("Iterative gapfilling", "igapfill", .caps.gapfill)

.register.capability("Cadzow iterations", "cadzow")
.register.capability("W-correlations", "wnorm")

.caps.lrr <- function(x) {
  ## We don't support LRR in shaped case
  !is.shaped(x) || x$circular
}
.register.capability("LRR", "lrr", .caps.lrr)

.caps.iossa <- function(x) {
  ## No complex case
  !inherits(x, "cssa")
}
.register.capability("Iterative O-SSA nested decomposition", "iossa", .caps.iossa)
.register.capability("Filter-adjusted O-SSA nested decomposition", "fossa")

.register.capability("Automatic grouping via w-correlations", "grouping.auto.wcor")

.caps.autossa <- function(x) {
  ## No periodogram in shaped case
  !is.shaped(x)
}
.register.capability("Automatic grouping via periodogram", "grouping.auto.pgram", .caps.autossa)

