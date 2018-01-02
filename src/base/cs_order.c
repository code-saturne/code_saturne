/*============================================================================
 * \file Functions related to the ordering of local arrays.
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
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_order.h"

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
 * Descend binary tree for the ordering of a cs_gnum_t (integer) array.
 *
 * parameters:
 *   number    <-- pointer to numbers of entities that should be ordered.
 *                 (if NULL, a default 1 to n numbering is considered)
 *   level     <-- level of the binary tree to descend
 *   nb_ent    <-- number of entities in the binary tree to descend
 *   order     <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_gnum_descend_tree(const cs_gnum_t   number[],
                         size_t            level,
                         const size_t      nb_ent,
                         cs_lnum_t         order[])
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
_order_gnum(const cs_gnum_t   number[],
            cs_lnum_t         order[],
            const size_t      nb_ent)
{
  size_t i;
  cs_lnum_t o_save;

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
 * cs_gnum_t array.
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
_order_gnum_descend_tree_s(const cs_gnum_t   number[],
                           size_t            stride,
                           size_t            level,
                           const size_t      nb_ent,
                           cs_lnum_t         order[])
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
_order_gnum_s(const cs_gnum_t   number[],
              size_t            stride,
              cs_lnum_t         order[],
              const size_t      nb_ent)
{
  size_t i;
  cs_lnum_t o_save;

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
                             const cs_lnum_t   index[],
                             const cs_gnum_t   number[])
{
  cs_lnum_t   i;

  cs_lnum_t   i1_s = index[i1], i1_e = index[i1+1], s1 = i1_e - i1_s;
  cs_lnum_t   i2_s = index[i2], i2_e = index[i2+1], s2 = i2_e - i2_s;

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
                    const cs_lnum_t   index[],
                    const cs_gnum_t   number[])
{
  cs_lnum_t   i;

  cs_lnum_t   i1_s = index[i1], i1_e = index[i1+1], s1 = i1_e - i1_s;
  cs_lnum_t   i2_s = index[i2], i2_e = index[i2+1], s2 = i2_e - i2_s;

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
_order_gnum_descend_tree_i(const cs_gnum_t   number[],
                           const cs_lnum_t   index[],
                           size_t            level,
                           const size_t      nb_ent,
                           cs_lnum_t         order[])
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
_order_gnum_i(const cs_gnum_t   number[],
              const cs_lnum_t   index[],
              cs_lnum_t         order[],
              const size_t      nb_ent)
{
  size_t i;
  cs_lnum_t o_save;

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

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of a cs_lnum_t (integer) array.
 *
 * parameters:
 *   number    <-- pointer to numbers of entities that should be ordered.
 *                 (if NULL, a default 1 to n numbering is considered)
 *   level     <-- level of the binary tree to descend
 *   nb_ent    <-- number of entities in the binary tree to descend
 *   order     <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_lnum_descend_tree(const cs_lnum_t   number[],
                         size_t            level,
                         const size_t      nb_ent,
                         cs_lnum_t         order[])
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
 * Order an array of local numbers.
 *
 * parameters:
 *   number   <-- array of entity numbers (if NULL, a default 1 to n
 *                numbering is considered)
 *   order    <-- pre-allocated ordering table
 *   nb_ent   <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_lnum(const cs_lnum_t   number[],
            cs_lnum_t         order[],
            const size_t      nb_ent)
{
  size_t i;
  cs_lnum_t o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < nb_ent ; i++)
    order[i] = i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2) ;
  do {
    i--;
    _order_lnum_descend_tree(number, i, nb_ent, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = nb_ent - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_lnum_descend_tree(number, 0, i, order);
  }
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the lexicographical ordering of a strided
 * cs_lnum_t array.
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
_order_lnum_descend_tree_s(const cs_lnum_t   number[],
                           size_t            stride,
                           size_t            level,
                           const size_t      nb_ent,
                           cs_lnum_t         order[])
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
_order_lnum_s(const cs_lnum_t   number[],
              size_t            stride,
              cs_lnum_t         order[],
              const size_t      nb_ent)
{
  size_t i;
  cs_lnum_t o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < nb_ent ; i++)
    order[i] = i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2) ;
  do {
    i--;
    _order_lnum_descend_tree_s(number, stride, i, nb_ent, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = nb_ent - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_lnum_descend_tree_s(number, stride, 0, i, order);
  }
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of a cs_real_t array.
 *
 * parameters:
 *   value    <-- pointer to values of entities that should be ordered.
 *   level    <-- level of the binary tree to descend
 *   nb_ent   <-- number of entities in the binary tree to descend
 *   order    <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_real_descend_tree(const cs_real_t   value[],
                         size_t            level,
                         const size_t      nb_ent,
                         cs_lnum_t         order[])
{
  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (nb_ent/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < nb_ent - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      if (value[i1] > value[i2]) lv_cur++;
    }

    if (lv_cur >= nb_ent) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (value[i1] >= value[i2]) break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order an array of local values.
 *
 * parameters:
 *   value   <-- array of entity values
 *   order   <-- pre-allocated ordering table
 *   nb_ent  <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_real(const cs_real_t   value[],
            cs_lnum_t         order[],
            const size_t      nb_ent)
{
  size_t i;
  cs_lnum_t o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < nb_ent ; i++)
    order[i] = i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2) ;
  do {
    i--;
    _order_real_descend_tree(value, i, nb_ent, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = nb_ent - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_real_descend_tree(value, 0, i, order);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if an array of global numbers is ordered.
 *
 * \param[in]  list    optional list (1 to n numbering) of selected entities
 *                     (or NULL if all nb_ent are selected). This list may
 *                     contain element numbers in any order
 * \param[in]  number  array of all entity numbers (number of entity i given
 *                     by number[i] or number[list[i] - 1]) if list exists
 *                     (if NULL, a default 1 to n numbering is considered)
 * \param[in]  nb_ent  number of entities considered
 *
 * \return  1 if ordered, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_order_gnum_test(const cs_lnum_t  list[],
                   const cs_gnum_t  number[],
                   size_t           nb_ent)
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return an ordering table associated with an array of global numbers.
 *
 * \param[in]  list    optional list (1 to n numbering) of selected entities
 *                     (or NULL if all nb_ent are selected). This list may
 *                     contain element numbers in any order
 * \param[in]  number  array of all entity numbers (number of entity i given
 *                     by number[i] or number[list[i] - 1]) if list exists
 *                     (if NULL, a default 1 to n numbering is considered)
 * \param[in]  nb_ent  number of entities considered
 *
 * \return  pointer to list of nb_ent entities (0 to n-1 numbering) ordered by
 *          increasing associated number. The calling code is responsible for
 *          freeing this array when it is not needed anymore.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t *
cs_order_gnum(const cs_lnum_t  list[],
              const cs_gnum_t  number[],
              size_t           nb_ent)
{
  cs_lnum_t *order;

  BFT_MALLOC(order, nb_ent, cs_lnum_t);

  cs_order_gnum_allocated(list,
                            number,
                            order,
                            nb_ent);

  return order;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a lexicographical ordering table associated with a strided
 * array of global numbers.
 *
 * \param[in]  list    optional list (1 to n numbering) of selected entities
 *                     (or NULL if all nb_ent are selected). This list may
 *                     contain element numbers in any order
 * \param[in]  number  array of all entity numbers (number of entity i
 *                     given by number[i] or number[list[i] - 1]) if list
 *                     exists (if NULL, a default 1 to n numbering is
 *                     considered)
 * \param[in]  stride  stride of number array (number of values to compare)
 * \param[in]  nb_ent  number of entities considered
 *
 * \return  pointer to list of nb_ent entities (0 to n-1 numbering) ordered by
 *          increasing associated number. The calling code is responsible for
 *          freeing this array when it is not needed anymore.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t *
cs_order_gnum_s(const cs_lnum_t  list[],
                const cs_gnum_t  number[],
                size_t           stride,
                size_t           nb_ent)
{
  cs_lnum_t *order;

  BFT_MALLOC(order, nb_ent, cs_lnum_t);

  cs_order_gnum_allocated_s(list,
                              number,
                              stride,
                              order,
                              nb_ent);

  return order;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a lexicographical ordering table associated with an indexed
 * array of global numbers.
 *
 * \param[in]  list    optional list (1 to n numbering) of selected entities
 *                     (or NULL if all nb_ent are selected). This list may
 *                     contain element numbers in any order
 * \param[in]  number  array of all entity numbers (numbers of entity i start
 *                     at index[i] or _index[i] (reduced index) if list exists).
 *                     If list = NULL, a default 1 to n numbering is considered)
 * \param[in]  index   number of values to compare for each entity
 * \param[in]  nb_ent  number of entities considered
 *
 * \return  pointer to list of nb_ent entities (0 to n-1 numbering) ordered by
 *          increasing associated number. The calling code is responsible for
 *          freeing this array when it is not needed anymore.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t *
cs_order_gnum_i(const cs_lnum_t  list[],
                const cs_gnum_t  number[],
                const cs_lnum_t  index[],
                size_t           nb_ent)
{
  cs_lnum_t *order = NULL;

  BFT_MALLOC(order, nb_ent, cs_lnum_t);

  cs_order_gnum_allocated_i(list, number, index, order, nb_ent);

  return order;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute an ordering table associated with an array of global numbers.
 *
 * \param[in]   list    optional list (1 to n numbering) of selected entities
 *                      (or NULL if all nb_ent are selected). This list may
 *                      contain element numbers in any order
 * \param[in]   number  array of all entity numbers (number of entity i given
 *                      by number[i] or number[list[i] - 1]) if list exists
 *                      (if NULL, a default 1 to n numbering is considered)
 * \param[out]  order   pointer to pre-allocated ordering table
 * \param[in]   nb_ent  number of entities considered
 */
/*----------------------------------------------------------------------------*/

void
cs_order_gnum_allocated(const cs_lnum_t  list[],
                        const cs_gnum_t  number[],
                        cs_lnum_t        order[],
                        size_t           nb_ent)
{
  size_t i;
  cs_gnum_t *number_list;

  /* Explicit numbering */

  if (number != NULL) {

    if (list != NULL) {
      BFT_MALLOC(number_list, nb_ent, cs_gnum_t);
      for (i = 0 ; i < nb_ent ; i++)
        number_list[i] = number[list[i] - 1];
      _order_gnum(number_list,
                   order,
                   nb_ent);
      BFT_FREE(number_list);
    }
    else
      _order_gnum(number,
                   order,
                   nb_ent);

  }

  /* Implicit numbering */

  else {

    if (list != NULL) {
      BFT_MALLOC(number_list, nb_ent, cs_gnum_t);
      for (i = 0 ; i < nb_ent ; i++)
        number_list[i] = (cs_gnum_t)(list[i]);
      _order_gnum(number_list,
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a lexicographical ordering table associated with an array of
 * strided global numbers.
 *
 * \param[in]   list    optional list (1 to n numbering) of selected entities
 *                      (or NULL if all nb_ent are selected). This list may
 *                      contain element numbers in any order
 * \param[in]   number  array of all entity numbers (numbers of entity i start
 *                      at number[i*stride] or number[(list[i] - 1)*stride])
 *                      if list exists (if NULL, a default 1 to n numbering is
 *                      considered)
 * \param[in]   stride  stride of number array (number of values to compare)
 * \param[out]  order   pointer to pre-allocated ordering table
 * \param[in]   nb_ent  number of entities considered
 */
/*----------------------------------------------------------------------------*/

void
cs_order_gnum_allocated_s(const cs_lnum_t  list[],
                          const cs_gnum_t  number[],
                          size_t           stride,
                          cs_lnum_t        order[],
                          const size_t     nb_ent)
{
  size_t i, j;
  cs_gnum_t *number_list;

  /* Explicit numbering */

  if (number != NULL) {

    if (list != NULL) {
      BFT_MALLOC(number_list, nb_ent*stride, cs_gnum_t);
      for (i = 0 ; i < nb_ent ; i++) {
        for (j = 0; j < stride; j++)
          number_list[i*stride + j] = number[(list[i] - 1)*stride + j];
      }
      _order_gnum_s(number_list,
                     stride,
                     order,
                     nb_ent);
      BFT_FREE(number_list);
    }
    else
      _order_gnum_s(number,
                     stride,
                     order,
                     nb_ent);

  }

  /* Implicit numbering */

  else

    cs_order_gnum_allocated(list,
                              NULL,
                              order,
                              nb_ent);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a lexicographical ordering table associated with an indexed
 *        array of global numbers.
 *
 * \param[in]   list    optional list (1 to n numbering) of selected entities
 *                      (or NULL if all nb_ent are selected). This list may
 *                      contain element numbers in any order
 * \param[in]   number  array of all entity numbers (numbers of entity i start
 *                      at index[i] or _index[i] (reduced index) if list
 *                      exists). If list = NULL, a default 1 to n numbering
 *                      is considered)
 * \param[in]   index   number of values to compare for each entity (from 0)
 * \param[out]  order   pointer to pre-allocated ordering table
 * \param[in]   nb_ent  number of entities considered
 */
/*----------------------------------------------------------------------------*/

void
cs_order_gnum_allocated_i(const cs_lnum_t  list[],
                          const cs_gnum_t  number[],
                          const cs_lnum_t  index[],
                          cs_lnum_t        order[],
                          size_t           nb_ent)
{
  /* Explicit numbering */

  if (number != NULL) {

    if (list != NULL) {

      size_t  i, j, k, ent_id, _shift;

      cs_lnum_t   *_index = NULL;
      cs_gnum_t   *number_list = NULL;

      BFT_MALLOC(_index, nb_ent + 1, cs_lnum_t);

      /* Count reduced size */

      for (i = 0; i < nb_ent; i++) {
        ent_id = list[i]-1;
        _index[i+1] = index[ent_id+1] - index[ent_id];
      }

      _index[0] = 0;
      for (i = 0; i < nb_ent; i++)
        _index[i+1] += _index[i];

      BFT_MALLOC(number_list, _index[nb_ent], cs_gnum_t);

      /* Define reduced index and adjacency */

      for (i = 0; i < nb_ent; i++) {

        ent_id = list[i]-1;
        _shift = _index[i];

        for (j = index[ent_id], k = 0;
             (cs_lnum_t)j < index[ent_id+1]; j++, k++)
          number_list[_shift + k] = number[j];

      }

      /* Local ordering */

      _order_gnum_i(number_list, _index, order, nb_ent);

      BFT_FREE(_index);
      BFT_FREE(number_list);

    }
    else { /* Local ordering */

      _order_gnum_i(number, index, order, nb_ent);

    }

  } /* number != NULL */

  /* Implicit numbering */

  else

    cs_order_gnum_allocated(list, NULL, order, nb_ent);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute an ordering table associated with an array of local numbers.
 *
 * \param[in]   list    optional list (1 to n numbering) of selected entities
 *                      (or NULL if all nb_ent are selected). This list may
 *                      contain element numbers in any order
 * \param[in]   number  array of all entity numbers (number of entity i given
 *                      by number[i] or number[list[i] - 1]) if list exists
 *                      (if NULL, a default 1 to n numbering is considered)
 * \param[out]  order   pointer to pre-allocated ordering table
 * \param[in]   nb_ent  number of entities considered
 */
/*----------------------------------------------------------------------------*/

void
cs_order_lnum_allocated(const cs_lnum_t  list[],
                        const cs_lnum_t  number[],
                        cs_lnum_t        order[],
                        size_t           nb_ent)
{
  size_t i;
  cs_lnum_t *number_list;

  /* Explicit numbering */

  if (number != NULL) {

    if (list != NULL) {
      BFT_MALLOC(number_list, nb_ent, cs_lnum_t);
      for (i = 0 ; i < nb_ent ; i++)
        number_list[i] = number[list[i] - 1];
      _order_lnum(number_list,
                  order,
                  nb_ent);
      BFT_FREE(number_list);
    }
    else
      _order_lnum(number,
                  order,
                  nb_ent);

  }

  /* Implicit numbering */

  else {

    if (list != NULL) {
      BFT_MALLOC(number_list, nb_ent, cs_lnum_t);
      for (i = 0 ; i < nb_ent ; i++)
        number_list[i] = (cs_gnum_t)(list[i]);
      _order_lnum(number_list,
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a lexicographical ordering table associated with an array of
 * strided local numbers.
 *
 * \param[in]   list    optional list (1 to n numbering) of selected entities
 *                      (or NULL if all nb_ent are selected). This list may
 *                      contain element numbers in any order
 * \param[in]   number  array of all entity numbers (numbers of entity i start
 *                      at number[i*stride] or number[(list[i] - 1)*stride])
 *                      if list exists (if NULL, a default 1 to n numbering is
 *                      considered)
 * \param[in]   stride  stride of number array (number of values to compare)
 * \param[out]  order   pointer to pre-allocated ordering table
 * \param[in]   nb_ent  number of entities considered
 */
/*----------------------------------------------------------------------------*/

void
cs_order_lnum_allocated_s(const cs_lnum_t  list[],
                          const cs_lnum_t  number[],
                          size_t           stride,
                          cs_lnum_t        order[],
                          const size_t     nb_ent)
{
  size_t i, j;
  cs_lnum_t *number_list;

  /* Explicit numbering */

  if (number != NULL) {

    if (list != NULL) {
      BFT_MALLOC(number_list, nb_ent*stride, cs_lnum_t);
      for (i = 0 ; i < nb_ent ; i++) {
        for (j = 0; j < stride; j++)
          number_list[i*stride + j] = number[(list[i] - 1)*stride + j];
      }
      _order_lnum_s(number_list,
                     stride,
                     order,
                     nb_ent);
      BFT_FREE(number_list);
    }
    else
      _order_lnum_s(number,
                     stride,
                     order,
                     nb_ent);

  }

  /* Implicit numbering */

  else

    cs_order_lnum_allocated(list,
                              NULL,
                              order,
                              nb_ent);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute an ordering table associated with an array of local values.
 *
 * \param[in]   list    optional list (1 to n numbering) of selected entities
 *                      (or NULL if all nb_ent are selected). This list may
 *                      contain element numbers in any order
 * \param[in]   val     array of all entity values (value of entity i given
 *                      by value[i] or value[list[i] - 1]) if list exists
 * \param[out]  order   pointer to pre-allocated ordering table
 * \param[in]   nb_ent  number of entities considered
 */
/*----------------------------------------------------------------------------*/

void
cs_order_real_allocated(const cs_lnum_t  list[],
                        const cs_real_t  val[],
                        cs_lnum_t        order[],
                        size_t           nb_ent)
{
  size_t i;
  cs_real_t *val_list;

  /* Explicit numbering */

  if (list != NULL) {
    BFT_MALLOC(val_list, nb_ent, cs_real_t);
    for (i = 0 ; i < nb_ent ; i++)
      val_list[i] = val[list[i] - 1];
    _order_real(val_list,
                order,
                nb_ent);
    BFT_FREE(val_list);
  }
  else
    _order_real(val,
                order,
                nb_ent);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build local renumbering array based on ordering of entities.
 *
 * \param[in]  order   0 to n-1 ordering of entities by increasing attribute
 * \param[in]  nb_ent  number of entities considered
 *
 * \return  pointer to renumbering array (0 to n-1 numbering) indicating the
 *          new index of renumbered entities; The calling code is responsible
 *          for freeing this array when it is not needed anymore.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t *
cs_order_renumbering(const cs_lnum_t  order[],
                     size_t           nb_ent)
{
  size_t i;
  cs_lnum_t *number;

  if (nb_ent < 1)
    return NULL;

  assert(order != NULL);

  BFT_MALLOC(number, nb_ent, cs_lnum_t);

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
/*!
 * \brief Reorder data based on ordering array.
 *
 * \param[in]      n_elts      number of elements
 * \param[in]      elt_size    element size
 * \param[in]      order       reordering array
 * \param[in,out]  data        data
 *
 * \return  new size of data
 */
/*----------------------------------------------------------------------------*/

void
cs_order_reorder_data(cs_lnum_t         n_elts,
                      size_t            elt_size,
                      const cs_lnum_t   order[],
                      void             *data)
{
  unsigned char *tmp;
  unsigned char *_data = data;

  BFT_MALLOC(tmp, n_elts*elt_size, unsigned char);

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    cs_lnum_t j = order[i];
    const unsigned char *src = _data + j*elt_size;
    unsigned char *dest = tmp + i*elt_size;
    for (size_t k = 0; k < elt_size; k++)
      dest[k] = src[k];
  }
  memcpy(data, tmp, n_elts*elt_size);

  BFT_FREE(tmp);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
