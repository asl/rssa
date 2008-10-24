#   R package for Static Singular Analysis
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

# FIXME: more checks
.wcor <- function(X, L) {
  if (missing(L))
    L <- dim(X)[2];
  N <- dim(X)[1];
  K <- N - L + 1;
  Ls <- min(L, K); Ks <- max(L, K);

  # Compute weights
  w <- c(1:Ls, rep(Ls, Ks-Ls), seq(from = N-Ks, to = 1, by = -1));

  # Compute w-covariation
  cov <- crossprod(sqrt(w) * X);

  # Convert to correlations
  Is <- 1/sqrt(diag(cov));
  Is * cov * rep(Is, each = nrow(cov));
}

N = 399;
a = 1.005;
T = 200;
F1 <- (1/a)^(1:N);
F2 <- (a^(1:N))*cos(2*pi*(1:N)/T);
.wcor(cbind(F1, F2));
