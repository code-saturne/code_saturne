#ifndef __FVM_ORDER_H__
#define __FVM_ORDER_H__

/*============================================================================
 * Functions related to the ordering of local arrays of global numbers and
 * calculation of global ranks in parallel mode.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2008  EDF

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
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
                     size_t            nb_ent);

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
                       size_t            nb_ent);

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
 *   freeing this array when it is not needed anymore
 *----------------------------------------------------------------------------*/

fvm_lnum_t *
fvm_order_local(const fvm_lnum_t  list[],
                const fvm_gnum_t  number[],
                size_t            nb_ent);

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
                  size_t            nb_ent);

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
                  size_t            nb_ent);

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
                          const size_t      nb_ent);

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
                            const size_t      nb_ent);

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
                            const size_t      nb_ent);

/*----------------------------------------------------------------------------
 * Build local renumbering array based on ordering of entities.
 *
 * parameters:
 *   order    <-- 0 to n-1 ordering of entities by increasing attribute
 *   nb_ent   <-- number of entities considered
 *
 * returns:
 *   pointer to renumbering array (0 to n-1 numbering) indicating the new
 *   index of renumbered entities; The calling code is responsible for
 *   freeing this array when it is not needed anymore
 *----------------------------------------------------------------------------*/

fvm_lnum_t *
fvm_order_local_renumbering(const fvm_lnum_t  order[],
                            const size_t      nb_ent);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_ORDER_H__ */
