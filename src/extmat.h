/*
 *   R package for Singular Spectrum Analysis
 *   Copyright (c) 2009 Anton Korobeynikov <asl@math.spbu.ru>
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

#ifndef __EXTMAT_H__
#define __EXTMAT_H__

/* External matrix structure */
typedef void (*mulfn) (double* out, const double* v, const void* matrix);
typedef unsigned (*infofn) (const void* matrix);

typedef struct {
  const char* type;
  void* matrix;
  mulfn mulfn;
  mulfn tmulfn;
  infofn ncol;
  infofn nrow;
} ext_matrix;

#endif /* __EXTMAT_H__ */
