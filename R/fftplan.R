#   R package for Singular Spectrum Analysis
#   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
#   Copyright (c) 2013 Alexander Shlemov <shlemovalex@gmail.com>
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

#   Routines for FFT plans

.get.or.create.fft.plan <- function(x) {
  .get.or.create(x, "fft.plan", fft.plan.1d(x$length))
}

fft.plan.1d <- function(N) {
  storage.mode(N) <- "integer"
  .Call("initialize_fft_plan", N)
}

is.fft.plan <- function(fft.plan) {
  .Call("is_fft_plan", fft.plan)
}

dim.fft.plan <- function(fft.plan) {
  .Call("dim_fft_plan", fft.plan)
}
