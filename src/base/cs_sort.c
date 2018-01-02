/*============================================================================
 * Functions related to in-place sorting of arrays.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of a cs_lnum_t (integer) array.
 *
 * parameters:
 *   number    <-> pointer to elements that should be ordered
 *   level     <-- level of the binary tree to descend
 *   n_elts    <-- number of elements in the binary tree to descend
 *----------------------------------------------------------------------------*/

static inline void
_sort_descend_tree(cs_lnum_t  number[],
                   size_t     level,
                   size_t     n_elts)
{
  size_t lv_cur;
  cs_lnum_t num_save;

  num_save = number[level];

  while (level <= (n_elts/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n_elts - 1)
      if (number[lv_cur+1] > number[lv_cur]) lv_cur++;

    if (lv_cur >= n_elts) break;

    if (num_save >= number[lv_cur]) break;

    number[level] = number[lv_cur];
    level = lv_cur;

  }

  number[level] = num_save;
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of a cs_gnum_t (integer) array.
 *
 * parameters:
 *   number    <-> pointer to elements that should be ordered
 *   level     <-- level of the binary tree to descend
 *   n_elts    <-- number of elements in the binary tree to descend
 *----------------------------------------------------------------------------*/

static inline void
_sort_descend_tree_gnum(cs_gnum_t  number[],
                        size_t     level,
                        size_t     n_elts)
{
  size_t lv_cur;
  cs_gnum_t num_save;

  num_save = number[level];

  while (level <= (n_elts/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n_elts - 1)
      if (number[lv_cur+1] > number[lv_cur]) lv_cur++;

    if (lv_cur >= n_elts) break;

    if (num_save >= number[lv_cur]) break;

    number[level] = number[lv_cur];
    level = lv_cur;

  }

  number[level] = num_save;
}

/*----------------------------------------------------------------------------
 * Order an array of global numbers.
 *
 * parameters:
 *   number   <-> array of numbers to sort
 *   n_elts   <-- number of elements considered
 *----------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER)
#pragma optimization_level 0 /* Bug above this with icc 15.0.1 20141023
                                for similar sort_lnum ? */
#endif

static void
_sort_gnum(cs_gnum_t  number[],
           size_t     n_elts)
{
  if (n_elts < 2)
    return;

  /* Use shell sort for short arrays */

  if (n_elts < 50) {

    size_t inc;

    /* Compute increment */
    for (inc = 1; inc <= n_elts/9; inc = 3*inc+1);

    /* Sort array */
    while (inc > 0) {
      for (size_t i = inc; i < n_elts; i++) {
        cs_gnum_t num_save = number[i];
        size_t j = i;
        while (j >= inc && number[j-inc] > num_save) {
          number[j] = number[j-inc];
          j -= inc;
        }
        number[j] = num_save;
      }
      inc = inc / 3;
    }

  }

  else {

    /* Create binary tree */

    size_t i = (n_elts / 2);
    do {
      i--;
      _sort_descend_tree_gnum(number, i, n_elts);
    } while (i > 0);

    /* Sort binary tree */

    for (size_t j = n_elts - 1; j > 0; j--) {
      cs_gnum_t num_save   = number[0];
      number[0] = number[j];
      number[j] = num_save;
      _sort_descend_tree_gnum(number, 0, j);
    }
  }
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the lexicographical ordering of a cs_gnum_t
 * couples array.
 *
 * parameters:
 *   number    <-> pointer to elements that should be ordered
 *   level     <-- level of the binary tree to descend
 *   n_elts    <-- number of elements in the binary tree to descend
 *----------------------------------------------------------------------------*/

static inline void
_sort_descend_tree_gnum_2(cs_gnum_t  number[],
                          size_t     level,
                          size_t     n_elts)
{
  size_t lv_cur;

  cs_gnum_t num_save[2] = {number[level*2], number[level*2+1]};

  while (level <= (n_elts/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n_elts - 1)
      if (   number[(lv_cur+1)*2] > number[lv_cur*2]
          || (   number[(lv_cur+1)*2] == number[lv_cur*2]
              && number[(lv_cur+1)*2+1] > number[lv_cur*2+1])) lv_cur++;

    if (lv_cur >= n_elts) break;

    if (! (  num_save[0] < number[lv_cur*2]
           || (   num_save[0] == number[lv_cur*2]
               && num_save[1] < number[lv_cur*2+1]))) break;

    number[level*2] = number[lv_cur*2];
    number[level*2+1] = number[lv_cur*2+1];
    level = lv_cur;

  }

  number[level*2] = num_save[0];
  number[level*2+1] = num_save[1];
}

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
  cs_lnum_t i, j, h;

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

/*----------------------------------------------------------------------------
 * Order an array of local numbers.
 *
 * parameters:
 *   number   <-> array of numbers to sort
 *   n_elts   <-- number of elements considered
 *----------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER)
#pragma optimization_level 0 /* Bug above this with icc 15.0.1 20141023 ? */
#endif

void
cs_sort_lnum(cs_lnum_t  number[],
             size_t     n_elts)
{
  if (n_elts < 2)
    return;

  /* Use shell sort for short arrays */

  if (n_elts < 50) {

    size_t inc;

    /* Compute increment */
    for (inc = 1; inc <= n_elts/9; inc = 3*inc+1);

    /* Sort array */
    while (inc > 0) {
      for (size_t i = inc; i < n_elts; i++) {
        cs_lnum_t num_save = number[i];
        size_t j = i;
        while (j >= inc && number[j-inc] > num_save) {
          number[j] = number[j-inc];
          j -= inc;
        }
        number[j] = num_save;
      }
      inc = inc / 3;
    }

  }

  else {

    /* Create binary tree */

    size_t i = (n_elts / 2);
    do {
      i--;
      _sort_descend_tree(number, i, n_elts);
    } while (i > 0);

    /* Sort binary tree */

    for (size_t j = n_elts - 1; j > 0; j--) {
      cs_lnum_t num_save   = number[0];
      number[0] = number[j];
      number[j] = num_save;
      _sort_descend_tree(number, 0, j);
    }
  }
}

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
                cs_lnum_t        elts[])
{
  bool retval = true;

  /* Sort line elements by column id (for better access patterns) */

# pragma omp parallel for  if(n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++) {
    cs_lnum_t *restrict _elts = elts + elt_idx[i];
    cs_lnum_t _n_elts = elt_idx[i+1] - elt_idx[i];
    cs_lnum_t id_prev = -1;
    cs_sort_lnum(_elts, _n_elts);
    for (cs_lnum_t j = 0; j < _n_elts; j++) {
      if (_elts[j] == id_prev)
        retval = false;
      id_prev = _elts[j];
    }
  }

  return retval;
}

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
                     cs_gnum_t        elts[])
{
  bool retval = true;

  /* Sort line elements by column id (for better access patterns) */

# pragma omp parallel for  if(n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++) {
    cs_gnum_t *restrict _elts = elts + elt_idx[i];
    cs_lnum_t _n_elts = elt_idx[i+1] - elt_idx[i];
    _sort_gnum(_elts, _n_elts);
    for (cs_lnum_t j = 1; j < _n_elts; j++) {
      if (_elts[j] == _elts[j-1])
        retval = false;
    }
  }

  return retval;
}

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
                         cs_gnum_t  elts[])
{
  if (n_elts < 2)
    return n_elts;

  /* Check if operations are needed, as overhead of checking should
     be significantly less than overhead of uneeded operations */

  bool no_need = true;
  for (cs_lnum_t i = 1; i < n_elts && no_need; i++) {
    if (elts[i] <= elts[i-1])
      no_need = false;
  }
  if (no_need)
    return n_elts;

  /* Use shell sort for short arrays */

  if (n_elts < 50) {

    cs_lnum_t inc;

    /* Compute increment */
    for (inc = 1; inc <= n_elts/9; inc = 3*inc+1);

    /* Sort array */
    while (inc > 0) {
      for (cs_lnum_t i = inc; i < n_elts; i++) {
        cs_gnum_t num_save = elts[i];
        cs_lnum_t j = i;
        while (j >= inc && elts[j-inc] > num_save) {
          elts[j] = elts[j-inc];
          j -= inc;
        }
        elts[j] = num_save;
      }
      inc = inc / 3;
    }

  }

  else {

    /* Create binary tree */

    cs_lnum_t i = (n_elts / 2);
    do {
      i--;
      _sort_descend_tree_gnum(elts, i, n_elts);
    } while (i > 0);

    /* Sort binary tree */

    for (cs_lnum_t j = n_elts - 1; j > 0; j--) {
      cs_gnum_t num_save   = elts[0];
      elts[0] = elts[j];
      elts[j] = num_save;
      _sort_descend_tree_gnum(elts, 0, j);
    }
  }

  /* Now compact array */

  cs_lnum_t j = 1;
  cs_gnum_t e_prev = elts[0];

  for (cs_lnum_t i = 1; i < n_elts; i++) {
    if (elts[i] != e_prev) {
      elts[j++] = elts[i];
      e_prev = elts[i];
    }
  }

  return j;
}

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
                           cs_gnum_t  elts[])
{
  if (n_elts < 2)
    return n_elts;

  /* Check if operations are needed, as overhead of checking should
     be significantly less than overhead of uneeded operations */

  bool no_need = true;
  for (cs_lnum_t i = 1; i < n_elts && no_need; i++) {
    if (   elts[i*2] <= elts[(i-1)*2]
        || (   elts[i*2] == elts[(i-1)*2]
            && elts[i*2+1] <= elts[(i-1)*2+1]))
      no_need = false;
  }
  if (no_need)
    return n_elts;

  /* Use shell sort for short arrays */

  if (n_elts < 50) {

    cs_lnum_t inc;

    /* Compute increment */
    for (inc = 1; inc <= n_elts/9; inc = 3*inc+1);

    /* Sort array */
    while (inc > 0) {
      for (cs_lnum_t i = inc; i < n_elts; i++) {
        cs_gnum_t num_save[2] = {elts[i*2], elts[i*2+1]};
        cs_lnum_t j = i;
        while (j >= inc
               && (   elts[(j-inc)*2] > num_save[0]
                   || (   elts[(j-inc)*2] == num_save[0]
                       && elts[(j-inc)*2+1] > num_save[1]))) {
          elts[j*2] = elts[(j-inc)*2];
          elts[j*2+1] = elts[(j-inc)*2+1];
          j -= inc;
        }
        elts[j*2] = num_save[0];
        elts[j*2+1] = num_save[1];
      }
      inc = inc / 3;
    }

  }

  else {

    /* Create binary tree */

    cs_lnum_t i = (n_elts / 2);
    do {
      i--;
      _sort_descend_tree_gnum_2(elts, i, n_elts);
    } while (i > 0);

    /* Sort binary tree */

    for (cs_lnum_t j = n_elts - 1; j > 0; j--) {
      cs_gnum_t num_save[2] = {elts[0], elts[1]};
      elts[0] = elts[j*2];
      elts[1] = elts[j*2+1];
      elts[j*2]   = num_save[0];
      elts[j*2+1] = num_save[1];
      _sort_descend_tree_gnum_2(elts, 0, j);
    }
  }

  /* Now compact array */

  cs_lnum_t j = 1;
  cs_gnum_t e_prev[2] = {elts[0], elts[1]};

  for (cs_lnum_t i = 1; i < n_elts; i++) {
    if (elts[i*2+1] != e_prev[1] || elts[i*2] != e_prev[0]) {
      elts[j*2]   = elts[i*2];
      elts[j*2+1] = elts[i*2+1];
      j++;
      e_prev[0] = elts[i*2];
      e_prev[1] = elts[i*2+1];
    }
  }

  return j;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

