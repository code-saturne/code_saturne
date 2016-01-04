/*============================================================================
 * Functions related to in-place sorting of arrays.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sort.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

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
                      cs_lnum_t        loc[])
{
  int i, j, h;
  cs_lnum_t  v;

  const cs_lnum_t  range = r-l;

  /* Compute stride */
  for (h = 1; h <= range/9; h = 3*h+1) ;

  /* Initialize loc */
  for (i = 0; i < range; i++)
    loc[i] = l+i;

  /* Sort array */
  for (; h > 0; h /= 3) {
    for (i = h; i < range; i++) {

      v = a[loc[i]], j = i;
      while ((j >= h) && (v < a[loc[j-h]]))
        loc[j] = loc[j-h], j -= h;
      loc[j] = loc[i];

    } /* Loop on array elements */
  } /* End of loop on stride */

}

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
              cs_lnum_t  a[])
{
  int i, j, h;
  cs_lnum_t  v;

  /* Compute stride */
  for (h = 1; h <= (r-l)/9; h = 3*h+1) ;

  /* Sort array */
  for (; h > 0; h /= 3) {
    for (i = l+h; i < r; i++) {

      cs_lnum_t  v = a[i];

      j = i;
      while ((j >= l+h) && (v < a[j-h]))
        a[j] = a[j-h], j -= h;
      a[j] = v;

    } /* Loop on array elements */
  } /* End of loop on stride */

}

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
                   cs_gnum_t  a[])
{
  int i, j, h;

  /* Compute stride */
  for (h = 1; h <= (r-l)/9; h = 3*h+1) ;

  /* Sort array */
  for (; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      cs_gnum_t  v = a[i];

      j = i;
      while ((j >= l+h) && (v < a[j-h])) {
        a[j] = a[j-h];
        j -= h;
      }
      a[j] = v;

    } /* Loop on array elements */

  } /* End of loop on stride */

}

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
cs_sort_coupled_shell(cs_lnum_t   l,
                      cs_lnum_t   r,
                      cs_lnum_t   a[],
                      cs_lnum_t   b[])
{
  int  i, j, h;
  cs_lnum_t  size = r - l;

  if (size == 0)
    return;

  /* Compute stride */
  for (h = 1; h <= size/9; h = 3*h+1) ;

  /* Sort array */
  for ( ; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      cs_lnum_t  va = a[i];
      cs_lnum_t  vb = b[i];

      j = i;
      while ( (j >= l+h) && (va < a[j-h]) ) {
        a[j] = a[j-h];
        b[j] = b[j-h];
        j -= h;
      }
      a[j] = va;
      b[j] = vb;

    } /* Loop on array elements */

  } /* End of loop on stride */

}

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
                       double  b[])
{
  int  i, j, h;

  int  size = r - l;

  if (size == 0)
    return;

  /* Compute stride */
  for (h = 1; h <= size/9; h = 3*h+1) ;

  /* Sort array */
  for ( ; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      int  va = a[i];
      double  vb = b[i];

      j = i;
      while ( (j >= l+h) && (va < a[j-h]) ) {
        a[j] = a[j-h];
        b[j] = b[j-h];
        j -= h;
      }
      a[j] = va;
      b[j] = vb;

    } /* Loop on array elements */

  } /* End of loop on stride */

}

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
                        int        a[],
                        short int  b[])
{
  int  i, j, h;

  int  size = r - l;

  if (size == 0)
    return;

  /* Compute stride */
  for (h = 1; h <= size/9; h = 3*h+1) ;

  /* Sort array */
  for ( ; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      int  va = a[i];
      short int  vb = b[i];

      j = i;
      while ( (j >= l+h) && (va < a[j-h]) ) {
        a[j] = a[j-h];
        b[j] = b[j-h];
        j -= h;
      }
      a[j] = va;
      b[j] = vb;

    } /* Loop on array elements */

  } /* End of loop on stride */

}

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
cs_sort_coupled_gnum_shell(cs_lnum_t   l,
                           cs_lnum_t   r,
                           cs_gnum_t   a[],
                           cs_gnum_t   b[])
{
  int  i, j, h;
  cs_lnum_t  size = r - l;

  if (size == 0)
    return;

  /* Compute stride */
  for (h = 1; h <= size/9; h = 3*h+1) ;

  /* Sort array */
  for ( ; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      cs_gnum_t  va = a[i];
      cs_gnum_t  vb = b[i];

      j = i;
      while ( (j >= l+h) && (va < a[j-h]) ) {
        a[j] = a[j-h];
        b[j] = b[j-h];
        j -= h;
      }
      a[j] = va;
      b[j] = vb;

    } /* Loop on array elements */

  } /* End of loop on stride */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS

