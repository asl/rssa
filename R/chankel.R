# R package for Singular Spectrum Analysis
# Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
# Copyright (c) 2013 Alexander Shlemov <shlemovalex@gmail.com>
#
# This program is free software; you can redistribute it
# and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation;
# either version 2 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the
# Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
# MA 02139, USA.

# Routines for SSA with centering

chankel <- function(X, L, lcv, rcv) {
  N <- length(X);
  K <- N - L + 1;
  h <- outer(1:L, 1:K, function(x, y) X[x + y - 1]);

  h - outer(lcv, rep(1, K)) - outer(rep(1, L), rcv);
}

new.chmat <- function(F, L, lcv, rcv) {
  storage.mode(F) <- "double";
  storage.mode(L) <- "integer";
  storage.mode(lcv) <- "double";
  storage.mode(rcv) <- "double";

  .Call("initialize_hmat", F, L, lcv, rcv);
}
