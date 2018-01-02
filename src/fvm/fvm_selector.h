#ifndef __FVM_SELECTOR_H__
#define __FVM_SELECTOR_H__

/*============================================================================
 * Mechanism for entity selection based on groups or attributes
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_group.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _fvm_selector_t  fvm_selector_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a selector object.
 *
 * parameters:
 *   dim                  <-- spatial dimension (coordinates and normals)
 *   n_elements           <-- number of selectable elements
 *   group_class_set      <-- pointer to group class set definition
 *   group_class_id       <-- group class id associated with each element
 *                            (size: n_elements)
 *   group_class_id_base; <-- Starting group class id base (usually 0 or 1)
 *   coords               <-- coordinates (interlaced) associated with each
 *                            element, whether vertex, face or cell center, ...
 *                            (size: n_elements * dim)
 *   normals              <-- normals (interlaced) associated with each element
 *                            if applicable (such as for face normals), or NULL
 *
 * returns:
 *   pointer to new selector
 *----------------------------------------------------------------------------*/

fvm_selector_t *
fvm_selector_create(int                           dim,
                    cs_lnum_t                     n_elements,
                    const fvm_group_class_set_t  *group_class_set,
                    const int                     group_class_id[],
                    int                           group_class_id_base,
                    const double                  coords[],
                    const double                  normals[]);

/*----------------------------------------------------------------------------
 * Destruction of a selector structure.
 *
 * parameters:
 *   this_selector <-> selector to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_selector_t *
fvm_selector_destroy(fvm_selector_t  *this_selector);

/*----------------------------------------------------------------------------
 * Define the list of the elements verifying the criteria described
 * by a character string
 *
 * The selected_element[] array must be pre-allocated, and be of sufficient
 * size to contain all elements associated with the selector.
 *
 * parameters:
 *   this_selector       <-> pointer to selector
 *   str                 <-- string defining selection criteria
 *   elt_id_base         <-- element id base (usually 0 or 1)
 *   n_selected_elements <-- number of elements selected
 *   selected_elements   <-> selected elements list (1 to n numbering)
 *
 * returns:
 *   criteria id associated by selector with str
 *----------------------------------------------------------------------------*/

int
fvm_selector_get_list(fvm_selector_t  *this_selector,
                      const char      *str,
                      cs_lnum_t        elt_id_base,
                      cs_lnum_t       *n_selected_elements,
                      cs_lnum_t       *selected_elements);

/*----------------------------------------------------------------------------
 * Define the list of group classes verifying the criteria described
 * by a character string.
 *
 * The selected_gc[] array must be pre-allocated, and be of sufficient
 * size to contain all elements associated with the selector.
 *
 * parameters:
 *   this_selector  <-> pointer to selector
 *   str            <-- string defining selection criteria
 *   n_selected_gcs <-- number of group classes selected
 *   selected_gcs   <-> selected group class list (0 to n numbering,
 *                      as group class "zero" may exist)
 *
 * returns:
 *   criteria id associated by selector with str
 *----------------------------------------------------------------------------*/

int
fvm_selector_get_gc_list(fvm_selector_t  *this_selector,
                         const char      *str,
                         int             *n_selected_gcs,
                         int              selected_gcs[]);

/*----------------------------------------------------------------------------
 * Return the number of operands associated with a selection criteria
 * which are missing in the selector's associated group class set.
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *   criteria_id   <-- id of criteria returned by fvm_selector_get_list()
 *
 * returns:
 *   number of missing operands
 *----------------------------------------------------------------------------*/

int
fvm_selector_n_missing(const fvm_selector_t  *this_selector,
                       int                    criteria_id);

/*----------------------------------------------------------------------------
 * Return a pointer to the name of an of operand associated with a selection
 * criteria which is missing in the selector's associated group class set.
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *   criteria_id   <-- id of criteria returned by fvm_selector_get_list()
 *   missing_id    <-- id of missing operand for this criteria
 *
 * returns:
 *   pointer to name of missing operand
 *----------------------------------------------------------------------------*/

const char *
fvm_selector_get_missing(const fvm_selector_t  *this_selector,
                         int                    criteria_id,
                         int                    missing_id);

/*----------------------------------------------------------------------------
 * Get statistics on selector usage
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *   n_evals       <-> number of evaluations, or NULL
 *   eval_wtime    <-> evaluation wall-clock time, or NULL
 *----------------------------------------------------------------------------*/

void
fvm_selector_get_stats(const fvm_selector_t  *this_selector,
                       int                   *n_evals,
                       double                *eval_wtime);

/*----------------------------------------------------------------------------
 * Dump the contents of a selector structure in human readable form
 *
 * parameters:
 *   this_selector <-- pointer to selector
 *----------------------------------------------------------------------------*/

void
fvm_selector_dump(const fvm_selector_t  *this_selector);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_SELECTOR_H__ */
