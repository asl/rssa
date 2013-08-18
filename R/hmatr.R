#   R package for Singular Spectrum Analysis
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

hmatr <- function(F, ...,
                  B = N %/% 4, T = N %/% 4, L = B %/% 2,
                  neig = 10) {
  N <- length(F)

  # Pre-calculate embedding vectors and their squared norms
  th <- t(hankel(F, L = L))
  cth2 <- c(0, cumsum(rowSums(th^2)))
  cth2 <- (cth2[1:(N-T)+(T-L+1)] - cth2[1:(N-T)])

  hc <- function(idx) {
    Fb <- F[idx:(idx+B)]   # Form a basis subspace
    s <- ssa(Fb, L = L, ..., neig = min(2*neig, 50))

    # Calculate the distance
    U <- s$U[, 1:neig, drop = FALSE]
    # FIXME: Can we use FFT stuff here somehow?
    cXU2 <- c(0, cumsum(rowSums((th %*% U)^2)))

    1 - (cXU2[1:(N-T)+(T-L+1)] - cXU2[1:(N-T)]) / cth2
  }

  h <- sapply(1:(N-B), hc)
  class(h) <- "hmatr"

  invisible(h)
}

plot.hmatr <- function(x,
                       col = rev(heat.colors(256)),
                       main = "Heterogeneity Matrix",
                       xlab = "", ylab = "",
                       ...) {
  image(1:nrow(x), 1:ncol(x), x,
        col = col, xlab = xlab, ylab = ylab, main = main, ...)
}
