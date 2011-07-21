/*============================================================================
 * Functions related to the ordering of local arrays of global numbers and
 * calculation of global ranks in parallel mode.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_order.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of a fvm_gnum (integer) array.
 *
 * parameters:
 *   number    <-- pointer to numbers of entities that should be ordered.
 *                 (if NULL, a default 1 to n numbering is considered)
 *   level     <-- level of the binary tree to descend
 *   nb_ent    <-- number of entities in the binary tree to descend
 *   order     <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_gnum_descend_tree(const fvm_gnum_t  number[],
                         size_t            level,
                         const size_t      nb_ent,
                         fvm_lnum_t        order[])
{
  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (nb_ent/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < nb_ent - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      if (number[i1] > number[i2]) lv_cur++;
    }

    if (lv_cur >= nb_ent) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (number[i1] >= number[i2]) break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order an array of global numbers.
 *
 * parameters:
 *   number   <-- array of entity numbers (if NULL, a default 1 to n
 *                numbering is considered)
 *   order    <-- pre-allocated ordering table
 *   nb_ent   <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_local(const fvm_gnum_t  number[],
             fvm_lnum_t        order[],
             const size_t      nb_ent)
{
  size_t i;
  fvm_lnum_t o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < nb_ent ; i++)
    order[i] = i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2) ;
  do {
    i--;
    _order_gnum_descend_tree(number, i, nb_ent, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = nb_ent - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_gnum_descend_tree(number, 0, i, order);
  }
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the lexicographical ordering of a strided
 * fvm_gnum array.
 *
 * parameters:
 *   number    <-- pointer to numbers of entities that should be ordered.
 *                 (if NULL, a default 1 to n numbering is considered)
 *   stride    <-- stride of array (number of values to compare)
 *   level     <-- level of the binary tree to descend
 *   nb_ent    <-- number of entities in the binary tree to descend
 *   order     <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_gnum_descend_tree_s(const fvm_gnum_t  number[],
                           size_t            stride,
                           size_t            level,
                           const size_t      nb_ent,
                           fvm_lnum_t        order[])
{
  size_t i_save, i1, i2, j, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (nb_ent/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < nb_ent - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      for (j = 0; j < stride; j++) {
        if (number[i1*stride + j] != number[i2*stride + j])
          break;
      }

      if (j < stride) {
        if (number[i1*stride + j] > number[i2*stride + j])
          lv_cur++;
      }

    }

    if (lv_cur >= nb_ent) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    for (j = 0; j < stride; j++) {
      if (number[i1*stride + j] != number[i2*stride + j])
        break;
    }

    if (j == stride) break;
    if (number[i1*stride + j] >= number[i2*stride + j]) break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order a strided array of global numbers lexicographically.
 *
 * parameters:
 *   number   <-- array of entity numbers (if NULL, a default 1 to n
 *                numbering is considered)
 *   stride   <-- stride of array (number of values to compare)
 *   order    --> pre-allocated ordering table
 *   nb_ent   <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_local_s(const fvm_gnum_t  number[],
               size_t            stride,
               fvm_lnum_t        order[],
               const size_t      nb_ent)
{
  size_t i;
  fvm_lnum_t o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < nb_ent ; i++)
    order[i] = i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2) ;
  do {
    i--;
    _order_gnum_descend_tree_s(number, stride, i, nb_ent, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = nb_ent - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_gnum_descend_tree_s(number, stride, 0, i, order);
  }
}

/*----------------------------------------------------------------------------
 * Indicate if element i1 from an indexed list is lexicographically
 * greater than or equal to element i2.
 *
 * parameters:
 *   i1        <-- position in index for the first element
 *   i2        <-- position in index for the second element
 *   index     <-- number of values to compare for each entity
 *   number    <-- pointer to numbers of entities that should be ordered.
 *                 (if NULL, a default 1 to n numbering is considered)
 *
 * returns:
 *   true if element i1 is greater or equal than element i2, false otherwise
 *----------------------------------------------------------------------------*/

inline static _Bool
_indexed_is_greater_or_equal(size_t            i1,
                             size_t            i2,
                             const fvm_lnum_t  index[],
                             const fvm_gnum_t  number[])
{
  fvm_lnum_t  i;

  fvm_lnum_t  i1_s = index[i1], i1_e = index[i1+1], s1 = i1_e - i1_s;
  fvm_lnum_t  i2_s = index[i2], i2_e = index[i2+1], s2 = i2_e - i2_s;

  if (s1 >= s2) {

    for (i = 0; i < s2; i++) {
      if (number[i1_s + i] > number[i2_s + i])
        return true;
      else if (number[i1_s + i] < number[i2_s + i])
        return false;
    }

    return true;
  }
  else { /* s1 < s2 */

    for (i = 0; i < s1; i++) {
      if (number[i1_s + i] > number[i2_s + i])
        return true;
      else if (number[i1_s + i] < number[i2_s + i])
        return false;
    }

    return false;
  }

}

/*----------------------------------------------------------------------------
 * Indicate if element i1 from an indexed list is lexicographically
 * strictly greater than element i2.
 *
 * parameters:
 *   i1        <-- position in index for the first element
 *   i2        <-- position in index for the second element
 *   index     <-- number of values to compare for each entity
 *   number    <-- pointer to numbers of entities that should be ordered.
 *                 (if NULL, a default 1 to n numbering is considered)
 *
 * returns:
 *   true if element i1 is strictly greater than element i2, false otherwise
 *----------------------------------------------------------------------------*/

inline static _Bool
_indexed_is_greater(size_t            i1,
                    size_t            i2,
                    const fvm_lnum_t  index[],
                    const fvm_gnum_t  number[])
{
  fvm_lnum_t  i;

  fvm_lnum_t  i1_s = index[i1], i1_e = index[i1+1], s1 = i1_e - i1_s;
  fvm_lnum_t  i2_s = index[i2], i2_e = index[i2+1], s2 = i2_e - i2_s;

  if (s1 > s2) {

    for (i = 0; i < s2; i++) {
      if (number[i1_s + i] > number[i2_s + i])
        return true;
      else if (number[i1_s + i] < number[i2_s + i])
        return false;
    }

    return true;
  }
  else { /* s1 <= s2 */

    for (i = 0; i < s1; i++) {
      if (number[i1_s + i] > number[i2_s + i])
        return true;
      else if (number[i1_s + i] < number[i2_s + i])
        return false;
    }

    return false;
  }

}

/*----------------------------------------------------------------------------
 * Descend binary tree for the lexicographical ordering of an indexed
 * array of global numbers.
 *
 * parameters:
 *   number    <-- pointer to numbers of entities that should be ordered.
 *                 (if NULL, a default 1 to n numbering is considered)
 *   index     <-- number of values to compare for each entity
 *   level     <-- level of the binary tree to descend
 *   nb_ent    <-- number of entities in the binary tree to descend
 *   order     <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_gnum_descend_tree_i(const fvm_gnum_t  number[],
                           const fvm_lnum_t  index[],
                           size_t            level,
                           const size_t      nb_ent,
                           fvm_lnum_t        order[])
{
  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (nb_ent/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < nb_ent - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      /* Test if element in position i1 is greater than element in
         position i2 */

      if (_indexed_is_greater(i1, i2, index, number))
        lv_cur++;

    }

    if (lv_cur >= nb_ent) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (_indexed_is_greater_or_equal(i1, i2, index, number))
      break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order an indexed array of global numbers lexicographically.
 *
 * parameters:
 *   number   <-- array of entity numbers (if NULL, a default 1 to n
 *                numbering is considered)
 *   index    <-- number of values to compare for each entity
 *   order    --> pre-allocated ordering table
 *   nb_ent   <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_local_i(const fvm_gnum_t  number[],
               const fvm_lnum_t  index[],
               fvm_lnum_t        order[],
               const size_t      nb_ent)
{
  size_t i;
  fvm_lnum_t o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < nb_ent ; i++)
    order[i] = i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2);
  do {
    i--;
    _order_gnum_descend_tree_i(number, index, i, nb_ent, order);
  } while (i > 0);


  /* Sort binary tree */

  for (i = nb_ent - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_gnum_descend_tree_i(number, index, 0, i, order);
  }

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Test if an array of global numbers is ordered.
 *
 * parameters:
 *   list     <-- optional list (1 to n numbering) of selected entities
 *                (or NULL if all nb_ent are selected). This list may
 *                contain element numbers in any order
 *   number   <-- array of all entity numbers (number of entity i
 *                given by number[i] or number[list[i] - 1]) if list exists
 *                (if NULL, a default 1 to n numbering is considered)
 *   nb_ent   <-- number of entities considered
 *
 * returns:
 *   1 if ordered, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
fvm_order_local_test(const fvm_lnum_t  list[],
                     const fvm_gnum_t  number[],
                     const size_t      nb_ent)
{
  size_t i = 0;

  /* If numbering is explicit */

  if (number != NULL) {

    if (list != NULL) {
      for (i = 1 ; i < nb_ent ; i++) {
        if (number[list[i] - 1] < number[list[i-1] - 1])
          break;
      }
    }
    else {
      for (i = 1 ; i < nb_ent ; i++) {
        if (number[i] < number[i-1])
          break;
      }
    }

  /* If numbering is implicit */

  }
  else {

    if (list != NULL) {
      for (i = 1 ; i < nb_ent ; i++) {
        if (list[i] < list[i-1])
          break;
      }
    }
    else
      i = nb_ent;
  }

  if (i == nb_ent || nb_ent == 0)
    return 1;
  else
    return 0;
}

/*----------------------------------------------------------------------------
 * Test if an array of global numbers is lexicographically ordered.
 *
 * parameters:
 *   list     <-- optional list (1 to n numbering) of selected entities
 *                (or NULL if all nb_ent are selected). This list may
 *                contain element numbers in any order
 *   number   <-- array of all entity numbers (number of entity i
 *                given by number[i] or number[list[i] - 1]) if list exists
 *                (if NULL, a default 1 to n numbering is considered)
 *   stride   <-- stride of number array (number of values to compare)
 *   nb_ent   <-- number of entities considered
 *
 * returns:
 *   1 if ordered, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
fvm_order_local_test_s(const fvm_lnum_t  list[],
                       const fvm_gnum_t  number[],
                       size_t            stride,
                       size_t            nb_ent)
{
  size_t j;
  size_t i = 0;

  /* If numbering is explicit */

  if (number != NULL) {

    if (list != NULL) {
      for (i = 1 ; i < nb_ent ; i++) {
        size_t j_prev, k;
        bool unordered = false;
        j_prev = list[i-1] - 1;
        j = list[i] - 1;
        for (k = 0; k < stride; k++) {
          if (number[j_prev*stride + k] < number[j*stride + k])
            break;
          else if (number[j_prev*stride + k] > number[j*stride + k])
            unordered = true;
        }
        if (unordered == true)
          break;
      }
    }
    else {
      for (i = 1 ; i < nb_ent ; i++) {
        size_t i_prev, k;
        bool unordered = false;
        i_prev = i-1;
        for (k = 0; k < stride; k++) {
          if (number[i_prev*stride + k] < number[i*stride + k])
            break;
          else if (number[i_prev*stride + k] > number[i*stride + k])
            unordered = true;
        }
        if (unordered == true)
          break;
      }
    }

  /* If numbering is implicit */

  }
  else {

    if (list != NULL) {
      for (i = 1 ; i < nb_ent ; i++) {
        if (list[i] < list[i-1])
          break;
      }
    }
    else
      i = nb_ent;
  }

  if (i == nb_ent || nb_ent == 0)
    return 1;
  else
    return 0;
}

/*----------------------------------------------------------------------------
 * Return an ordering table associated with an array of global numbers.
 *
 * parameters:
 *   list     <-- optional list (1 to n numbering) of selected entities
 *                (or NULL if all nb_ent are selected). This list may
 *                contain element numbers in any order
 *   number   <-- array of all entity numbers (number of entity i
 *                given by number[i] or number[list[i] - 1]) if list exists
 *                (if NULL, a default 1 to n numbering is considered)
 *   nb_ent   <-- number of entities considered
 *
 * returns:
 *   pointer to list of nb_ent entities (0 to n-1 numbering) ordered by
 *   increasing associated number. The calling code is responsible for
 *   freeing this array when it is not needed anymore.
 *----------------------------------------------------------------------------*/

fvm_lnum_t *
fvm_order_local(const fvm_lnum_t  list[],
                const fvm_gnum_t  number[],
                size_t            nb_ent)
{
  fvm_lnum_t *order;

  BFT_MALLOC(order, nb_ent, fvm_lnum_t);

  fvm_order_local_allocated(list,
                            number,
                            order,
                            nb_ent);

  return order;
}

/*----------------------------------------------------------------------------
 * Return a lexicographical ordering table associated with a strided array
 * of global numbers.
 *
 * parameters:
 *   list     <-- optional list (1 to n numbering) of selected entities
 *                (or NULL if all nb_ent are selected). This list may
 *                contain element numbers in any order
 *   number   <-- array of all entity numbers (number of entity i
 *                given by number[i] or number[list[i] - 1]) if list exists
 *                (if NULL, a default 1 to n numbering is considered)
 *   stride   <-- stride of number array (number of values to compare)
 *   nb_ent   <-- number of entities considered
 *
 * returns:
 *   pointer to list of nb_ent entities (0 to n-1 numbering) ordered by
 *   increasing associated number. The calling code is responsible for
 *   freeing this array when it is not needed anymore.
 *----------------------------------------------------------------------------*/

fvm_lnum_t *
fvm_order_local_s(const fvm_lnum_t  list[],
                  const fvm_gnum_t  number[],
                  size_t            stride,
                  size_t            nb_ent)
{
  fvm_lnum_t *order;

  BFT_MALLOC(order, nb_ent, fvm_lnum_t);

  fvm_order_local_allocated_s(list,
                              number,
                              stride,
                              order,
                              nb_ent);

  return order;
}

/*----------------------------------------------------------------------------
 * Return a lexicographical ordering table associated with an indexed array
 * of global numbers.
 *
 * parameters:
 *   list     <-- optional list (1 to n numbering) of selected entities
 *                (or NULL if all nb_ent are selected). This list may
 *                contain element numbers in any order
 *   number   <-- array of all entity numbers (numbers of entity i start
 *                at index[i] or _index[i] (reduced index) if list exists).
 *                If list = NULL, a default 1 to n numbering is considered)
 *   index    <-- number of values to compare for each entity
 *   nb_ent   <-- number of entities considered
 *
 * returns:
 *   pointer to list of nb_ent entities (0 to n-1 numbering) ordered by
 *   increasing associated number. The calling code is responsible for
 *   freeing this array when it is not needed anymore.
 *----------------------------------------------------------------------------*/

fvm_lnum_t *
fvm_order_local_i(const fvm_lnum_t  list[],
                  const fvm_gnum_t  number[],
                  const fvm_lnum_t  index[],
                  size_t            nb_ent)
{
  fvm_lnum_t *order = NULL;

  BFT_MALLOC(order, nb_ent, fvm_lnum_t);

  fvm_order_local_allocated_i(list, number, index, order, nb_ent);

  return order;
}

/*----------------------------------------------------------------------------
 * Compute an ordering table associated with an array of global numbers.
 *
 * parameters:
 *   list     <-- optional list (1 to n numbering) of selected entities
 *                (or NULL if all nb_ent are selected). This list may
 *                contain element numbers in any order
 *   number   <-- array of all entity numbers (number of entity i
 *                given by number[i] or number[list[i] - 1]) if list exists
 *                (if NULL, a default 1 to n numbering is considered)
 *   order    --> pointer to pre-allocated ordering table
 *   nb_ent   <-- number of entities considered
 *----------------------------------------------------------------------------*/

void
fvm_order_local_allocated(const fvm_lnum_t  list[],
                          const fvm_gnum_t  number[],
                          fvm_lnum_t        order[],
                          const size_t      nb_ent)
{
  size_t i;
  fvm_gnum_t *number_list;

  /* Explicit numbering */

  if (number != NULL) {

    if (list != NULL) {
      BFT_MALLOC(number_list, nb_ent, fvm_gnum_t);
      for (i = 0 ; i < nb_ent ; i++)
        number_list[i] = number[list[i] - 1];
      _order_local(number_list,
                   order,
                   nb_ent);
      BFT_FREE(number_list);
    }
    else
      _order_local(number,
                   order,
                   nb_ent);

  }

  /* Implicit numbering */

  else {

    if (list != NULL) {
      BFT_MALLOC(number_list, nb_ent, fvm_gnum_t);
      for (i = 0 ; i < nb_ent ; i++)
        number_list[i] = (fvm_gnum_t)(list[i]);
      _order_local(number_list,
                   order,
                   nb_ent);
      BFT_FREE(number_list);
    }
    else {
      for (i = 0 ; i < nb_ent ; i++)
        order[i] = i;
    }

  }

}

/*----------------------------------------------------------------------------
 * Compute a lexicographical ordering table associated with an array of
 * strided global numbers.
 *
 * parameters:
 *   list     <-- optional list (1 to n numbering) of selected entities
 *                (or NULL if all nb_ent are selected). This list may
 *                contain element numbers in any order
 *   number   <-- array of all entity numbers (numbers of entity i start
 *                at number[i*stride] or number[(list[i] - 1)*stride]) if
 *                list exists (if NULL, a default 1 to n numbering is
 *                considered)
 *   stride   <-- stride of number array (number of values to compare)
 *   order    --> pointer to pre-allocated ordering table
 *   nb_ent   <-- number of entities considered
 *----------------------------------------------------------------------------*/

void
fvm_order_local_allocated_s(const fvm_lnum_t  list[],
                            const fvm_gnum_t  number[],
                            size_t            stride,
                            fvm_lnum_t        order[],
                            const size_t      nb_ent)
{
  size_t i, j;
  fvm_gnum_t *number_list;

  /* Explicit numbering */

  if (number != NULL) {

    if (list != NULL) {
      BFT_MALLOC(number_list, nb_ent*stride, fvm_gnum_t);
      for (i = 0 ; i < nb_ent ; i++) {
        for (j = 0; j < stride; j++)
          number_list[i*stride + j] = number[(list[i] - 1)*stride + j];
      }
      _order_local_s(number_list,
                     stride,
                     order,
                     nb_ent);
      BFT_FREE(number_list);
    }
    else
      _order_local_s(number,
                     stride,
                     order,
                     nb_ent);

  }

  /* Implicit numbering */

  else

    fvm_order_local_allocated(list,
                              NULL,
                              order,
                              nb_ent);

}

/*----------------------------------------------------------------------------
 * Compute a lexicographical ordering table associated with an indexed array
 * of global numbers.
 *
 * parameters:
 *   list     <-- optional list (1 to n numbering) of selected entities
 *                (or NULL if all nb_ent are selected). This list may
 *                contain element numbers in any order
 *   number   <-- array of all entity numbers (numbers of entity i start
 *                at index[i] or _index[i] (reduced index) if list exists).
 *                If list = NULL, a default 1 to n numbering is considered)
 *   index    <-- number of values to compare for each entity (from 0)
 *   order    --> pointer to pre-allocated ordering table
 *   nb_ent   <-- number of entities considered
 *----------------------------------------------------------------------------*/

void
fvm_order_local_allocated_i(const fvm_lnum_t  list[],
                            const fvm_gnum_t  number[],
                            const fvm_lnum_t  index[],
                            fvm_lnum_t        order[],
                            const size_t      nb_ent)
{
  /* Explicit numbering */

  if (number != NULL) {

    if (list != NULL) {

      size_t  i, j, k, ent_id, _shift;

      fvm_lnum_t  *_index = NULL;
      fvm_gnum_t  *number_list = NULL;

      BFT_MALLOC(_index, nb_ent + 1, fvm_lnum_t);

      /* Count reduced size */

      for (i = 0; i < nb_ent; i++) {
        ent_id = list[i]-1;
        _index[i+1] = index[ent_id+1] - index[ent_id];
      }

      _index[0] = 0;
      for (i = 0; i < nb_ent; i++)
        _index[i+1] += _index[i];

      BFT_MALLOC(number_list, _index[nb_ent], fvm_gnum_t);

      /* Define reduced index and adjacency */

      for (i = 0; i < nb_ent; i++) {

        ent_id = list[i]-1;
        _shift = _index[i];

        for (j = index[ent_id], k = 0;
             (fvm_lnum_t)j < index[ent_id+1]; j++, k++)
          number_list[_shift + k] = number[j];

      }

      /* Local ordering */

      _order_local_i(number_list, _index, order, nb_ent);

      BFT_FREE(_index);
      BFT_FREE(number_list);

    }
    else { /* Local ordering */

      _order_local_i(number, index, order, nb_ent);

    }

  } /* number != NULL */

  /* Implicit numbering */

  else

    fvm_order_local_allocated(list, NULL, order, nb_ent);

}

/*----------------------------------------------------------------------------
 * Build local renumbering array based on ordering of entities.
 *
 * parameters:
 *   order    <-- 0 to n-1 ordering of entities by increasing attribute
 *   nb_ent   <-- number of entities considered
 *   base     <-- renumbering starting number (typically 0 or 1)
 *
 * returns:
 *   pointer to renumbering array (0 to n-1 numbering) indicating the new
 *   index of renumbered entities; The calling code is responsible for
 *   freeing this array when it is not needed anymore
 *----------------------------------------------------------------------------*/

fvm_lnum_t *
fvm_order_local_renumbering(const fvm_lnum_t  order[],
                            const size_t      nb_ent)
{
  size_t i;
  fvm_lnum_t *number;

  if (nb_ent < 1)
    return NULL;

  assert(order != NULL);

  BFT_MALLOC(number, nb_ent, fvm_lnum_t);

#if defined(DEBUG) && !defined(NDEBUG)
  /* Initialize with "impossible" number (so as to have a reproducible
     and detectable error in the renumbering array in case of
     incorrect order[] values) */
  for (i = 0 ; i < nb_ent ; i++)
    number[i] = - 1;
#endif

  /* Populate renumbering array */

  for (i = 0 ; i < nb_ent ; i++)
    number[order[i]] = i;

#if defined(DEBUG) && !defined(NDEBUG)
  /* Check renumbering array */
  for (i = 0 ; i < nb_ent ; i++)
    assert(number[i] >= 0);
#endif

  return number;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
