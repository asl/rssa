/*
 *   R package for Singular Spectrum Analysis
 *   Copyright (c) 2009-2010 Anton Korobeynikov <asl@math.spbu.ru>
 *   Copyright (c) 2013 Konstantin Usevich <konstantin.usevich@statmod.ru>
 *   Copyright (c) 2014 Alex Shlemov <shlemovalex@gmail.com>
 *
 *   This program is free software; you can redistribute it
 *   and/or modify it under the terms of the GNU General Public
 *   License as published by the Free Software Foundation;
 *   either version 2 of the License, or (at your option)
 *   any later version.
 *
 *   This program is distributed in the hope that it will be
 *   useful, but WITHOUT ANY WARRANTY; without even the implied
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *   PURPOSE.  See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public
 *   License along with this program; if not, write to the
 *   Free Software Foundation, Inc., 675 Mass Ave, Cambridge,
 *   MA 02139, USA.
 */

#ifndef __MASKS_H__
#define __MASKS_H__

#include <R.h>
#include <Rinternals.h>

typedef struct {
  R_len_t num;
  R_len_t *ind; /* Indices in an N array or an Nx x Ny array */
} area_indices;

void free_area(area_indices *area);
unsigned *alloc_weights(SEXP weights);
#endif /* __MASKS_H__ */
