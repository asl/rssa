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

#   Routines for operations with hankel matrices.

tcirc.old <- function(F, L = (N - 1) %/% 2) {
  N <- length(F); K = N - L + 1;
  .res <- list();
  .res$C <- as.vector(fft(c(F[K:N], F[1:(K-1)])));
  .res$L <- L;
  return (.res);
}

hmatmul.old <- function(C, v) {
  v <- as.vector(fft(C$C * fft(c(rev(v), rep(0, C$L-1))), inverse = TRUE));
  Re((v/length(C$C))[1:C$L]);
}

new.hmat <- function(F,
                     L = (N - 1) %/% 2) {
  N <- length(F);
  storage.mode(F) <- "double";
  storage.mode(L) <- "integer";
  .Call("initialize_hmat", F, L);
}

hmatmul <- function(hmat, v, transposed = FALSE) {
  storage.mode(v) <- "double";
  storage.mode(transposed) <- "logical";
  .Call("hmatmul", hmat, v, transposed);
}

mes <- function(N = 1000, L = (N %/% 2), n = 50) {
  F <- rnorm(N);
  v <- rnorm(N - L + 1);
  C <- tcirc.old(F, L = L);
  X <- hankel(F, L = L);
  h <- new.hmat(F, L = L);
  st1 <- system.time(for (i in 1:n) X %*% v);
  st2 <- system.time(for (i in 1:n) hmatmul.old(C, v));
  st3 <- system.time(for (i in 1:n) hmatmul(h, v));
  c(st1[["user.self"]], st2[["user.self"]], st3[["user.self"]]);
}


#Rprof();
#for (i in 1:250) {
#  r1 <- X %*% v;
#  r2 <- hmul(C, v);
#}
#Rprof(NULL);
#print(max(abs(r1-r2)));
#summaryRprof();
