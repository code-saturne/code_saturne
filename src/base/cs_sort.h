#ifndef __CS_SORT_H__
#define __CS_SORT_H__

/*============================================================================
 * Functions related to in-place sorting of arrays.
 *===========================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *===========================================================================*/

/*============================================================================
 * Type definitions
 *===========================================================================*/

/*=============================================================================
 * Static global variables
 *===========================================================================*/

/*=============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Sort an array "a" between its left bound "l" and its right bound "r"
 * thanks to a shell sort (Knuth algorithm).
 * Index location of the sorted array are stored in loc. a is unchanged.
 *
 * parameters:
 *   l    <-- left bound
 *   r    <-- right bound
 *   a    <-> array to sort (not modified)
 *   loc  <-> position by increasing order (size = r-l)
 *---------------------------------------------------------------------------*/

void
cs_sort_shell_inplace(cs_lnum_t        l,
                      cs_lnum_t        r,
                      const cs_lnum_t  a[],
                      cs_lnum_t        loc[]);

/*----------------------------------------------------------------------------
 * Sort an array "a" between its left bound "l" and its right bound "r"
 * thanks to a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l <-- left bound
 *   r <-- right bound
 *   a <-> array to sort
 *---------------------------------------------------------------------------*/

void
cs_sort_shell(cs_lnum_t  l,
              cs_lnum_t  r,
              cs_lnum_t  a[]);

/*----------------------------------------------------------------------------
 * Sort an array of integers "a" between its left bound "l" and its
 * right bound "r" using a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l <-- left bound
 *   r <-- right bound
 *   a <-> array to sort
 *---------------------------------------------------------------------------*/

void
cs_sort_int_shell(cs_lnum_t  l,
                  cs_lnum_t  r,
                  int        a[]);

/*----------------------------------------------------------------------------
 * Sort a global array "a" between its left bound "l" and its right bound "r"
 * thanks to a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l <-- left bound
 *   r <-- right bound
 *   a <-> array to sort
 *---------------------------------------------------------------------------*/

void
cs_sort_gnum_shell(cs_lnum_t  l,
                   cs_lnum_t  r,
                   cs_gnum_t  a[]);

/*----------------------------------------------------------------------------
 * Sort an array "a" and apply the sort to its associated array "b" (local
 * numbering)
 * Sort is realized thanks to a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l     -->   left bound
 *   r     -->   right bound
 *   a     <->   array to sort
 *   b     <->   associated array
 *---------------------------------------------------------------------------*/

void
cs_sort_coupled_shell(cs_lnum_t  l,
                      cs_lnum_t  r,
                      cs_lnum_t  a[],
                      cs_lnum_t  b[]);

/*----------------------------------------------------------------------------
 * Sort an array "a" and apply the sort to its associated array "b" (local
 * numbering)
 * Sort is realized thanks to a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l     -->   left bound
 *   r     -->   right bound
 *   a     <->   array to sort
 *   b     <->   associated array
 *---------------------------------------------------------------------------*/

void
cs_sort_dcoupled_shell(int     l,
                       int     r,
                       int     a[],
                       double  b[]);

/*----------------------------------------------------------------------------
 * Sort an array "a" and apply the sort to its associated array "b" (local
 * numbering)
 * Sort is realized thanks to a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l     -->   left bound
 *   r     -->   right bound
 *   a     <->   array to sort
 *   b     <->   associated array
 *---------------------------------------------------------------------------*/

void
cs_sort_sicoupled_shell(int        l,
                        int        r,
                        cs_lnum_t  a[],
                        short int  b[]);

/*----------------------------------------------------------------------------
 * Sort an array "a" and apply the sort to its associated array "b" (local
 * numbering)
 * Sort is realized thanks to a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l     -->   left bound
 *   r     -->   right bound
 *   a     <->   array to sort
 *   b     <->   associated array
 *---------------------------------------------------------------------------*/

void
cs_sort_coupled_gnum_shell(cs_lnum_t  l,
                           cs_lnum_t  r,
                           cs_gnum_t  a[],
                           cs_gnum_t  b[]);

/*----------------------------------------------------------------------------
 * Order an array of local numbers.
 *
 * parameters:
 *   number   <-> array of numbers to sort
 *   n_elts   <-- number of elements considered
 *----------------------------------------------------------------------------*/

void
cs_sort_lnum(cs_lnum_t  number[],
             size_t     n_elts);

/*----------------------------------------------------------------------------
 * Sort rows of an indexed structure.
 *
 * parameters:
 *   n_elts  <-- number of indexed elements
 *   elt_idx <-- element index (size: n_elts+1)
 *   elts    <-> indexed values
 *
 * returns:
 *   true if no values were encountered multiple times in a given row
 *---------------------------------------------------------------------------*/

bool
cs_sort_indexed(cs_lnum_t        n_elts,
                const cs_lnum_t  elt_idx[],
                cs_lnum_t        elts[]);

/*----------------------------------------------------------------------------
 * Sort rows of an indexed structure of global ids
 *
 * parameters:
 *   n_elts  <-- number of indexed elements
 *   elt_idx <-- element index (size: n_elts+1)
 *   elts    <-> indexed values
 *
 * returns:
 *   true if no values were encountered multiple times in a given row
 *---------------------------------------------------------------------------*/

bool
cs_sort_indexed_gnum(cs_lnum_t        n_elts,
                     const cs_lnum_t  elt_idx[],
                     cs_gnum_t        elts[]);

/*----------------------------------------------------------------------------
 * Sort an array of global number and remove duplicate entries.
 *
 * The calling code may choose to reallocate the array to the new, compacted
 * size; this is not done automatically, to avoid the overhead of memory
 * management in cases where it is not useful (i.e. when the array is
 * discarded soon after use).
 *
 * parameters:
 *   n_elts <-- initial number of elements
 *   elts   <-> elements to sort
 *
 * returns:
 *   number of compacted elements
 *---------------------------------------------------------------------------*/

cs_lnum_t
cs_sort_and_compact_gnum(cs_lnum_t  n_elts,
                         cs_gnum_t  elts[]);

/*----------------------------------------------------------------------------
 * Sort an array of global number couples and remove duplicate entries.
 *
 * Lexicographical ordering is considered.
 *
 * The calling code may choose to reallocate the array to the new, compacted
 * size; this is not done automatically, to avoid the overhead of memory
 * management in cases where it is not useful (i.e. when the array is
 * discarded soon after use).
 *
 * parameters:
 *   n_elts <-- initial number of elements
 *   elts   <-> elements to sort (size: n_elts*2, interleaved)
 *
 * returns:
 *   number of compacted elements
 *---------------------------------------------------------------------------*/

cs_lnum_t
cs_sort_and_compact_gnum_2(cs_lnum_t  n_elts,
                           cs_gnum_t  elts[]);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SORT_H__ */
