#   R package for Singular Spectrum Analysis
#   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
#   Copyright (c) 2009 Konstantin Usevich <usevich.k.d@gmail.com>
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

#   Routines for hankel-block hankel (aka 2d) SSA

tcircpile <- function(F, Lx = (Nx - 1) %/% 2, Ly = (Ny - 1) %/% 2) {
  Nx <- nrow(F); Ny <- ncol(F);
  Kx <- Nx - Lx + 1; Ky <- Ny - Ly + 1;

  .res <- list(Lx = Lx, Ly = Ly,
               Kx = Kx, Ky = Ky);

  TF <- cbind(F[,Ky:Ny],F[,1:(Ky-1)]);
  .res$Cblock <- fft(rbind(TF[Kx:Nx,],TF[1:(Kx-1),]));

  .res;
}

hbhmatmul.old <- function(C, v) {
  revv <- matrix(c(rev(v), rep(0, C$Kx*(C$Ly-1))), C$Kx, ncol(C$Cblock));
  revv <- rbind(revv, matrix(0, (C$Lx-1), ncol(revv)));

  mult <- fft(C$Cblock * fft(revv), inverse = TRUE);

  Re((mult/(prod(dim(C$Cblock))))[1:C$Lx,1:C$Ly]);
}

thbhmatmul.old <- function(C, v) {
  revv <- matrix(c(rep(0, C$Lx*(C$Ky-1)), rev(v)), C$Lx, ncol(C$Cblock));
  revv <- rbind(matrix(0, (C$Kx-1), ncol(revv)), revv);

  mult <- fft(C$Cblock * fft(revv), inverse = TRUE);

  Re((mult/(prod(dim(C$Cblock))))[C$Lx:(C$Lx+C$Kx-1),C$Ly:(C$Ly+C$Ky-1)]);
}

new.hbhmat <- function(F,
                       Lx = (Nx - 1) %/% 2,
                       Ly = (Ny - 1) %/% 2) {
  Nx <- nrow(F); Ny <- ncol(F);
  storage.mode(F) <- "double";
  storage.mode(Lx) <- "integer";
  storage.mode(Ly) <- "integer";
  h <- .Call("initialize_hbhmat", F, Lx, Ly);
}

hbhcols <- function(h) {
  .Call("hbhankel_cols", h)
}

hbhrows <- function(h) {
  .Call("hbhankel_rows", h)
}

is.hbhmat <- function(h) {
  .Call("is_hbhmat", h)
}

hbhmatmul <- function(hmat, v, transposed = FALSE) {
  storage.mode(v) <- "double";
  storage.mode(transposed) <- "logical";
  .Call("hbhmatmul", hmat, v, transposed);
}

mes <- function(Nx = 200, Ny = 90, Lx = 100, Ly = 50, n = 50) {
  Kx <- Nx - Lx +1;
  Ky <- Ny - Ly +1;
  F <- matrix(rnorm(Nx*Ny),Nx,Ny);
  #F <- matrix(1:(Nx*Ny),Nx,Ny);
  C <- tcircpile(F, Lx,Ly);
  X <- outer(1:(Lx*Ly), 1:(Kx*Ky), function(x,z) { F[(((x-1)%%Lx)+((z-1)%%Kx))+(((x-1)%/%Lx)+((z-1)%/%Kx))*Nx+1] });

  #X <- matrix(X, Lx*Ly, Kx*Ky);

  v <- rnorm(Kx*Ky);
  st1 <- system.time(for (i in 1:n) X %*% v);
  st2 <- system.time(for (i in 1:n) hbhmatmul(C, matrix(v,Kx,Ky)));

  print(c(st1[["user.self"]],st2[["user.self"]]),5);

  #v1 <- X %*% v;
  #v2 <- hbhmatmul(C, matrix(v,Kx,Ky));
  #print(max(abs(v1-v2)));
}
