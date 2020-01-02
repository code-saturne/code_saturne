#ifndef __FVM_SELECTOR_POSTFIX_H__
#define __FVM_SELECTOR_POSTFIX_H__

/*============================================================================
 * Expression handling for entity selection based on groups or attributes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

typedef struct _fvm_selector_postfix_t  fvm_selector_postfix_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a postfix expression from an infix expression
 *
 * parameters:
 *   infix        <-- infix expression
 *   n_groups     <-- number of groups
 *   n_attributes <-- number of attributes
 *   group_name   <-- array group names (sorted)
 *   attribute    <-- array of attribute numbers (sorted)
 *
 * returns:
 *   pointer to created postfix structure
 *----------------------------------------------------------------------------*/

fvm_selector_postfix_t *
fvm_selector_postfix_create(const char  *infix,
                            int          n_groups,
                            int          n_attributes,
                            const char  *group_name[],
                            const int    attribute[]);

/*----------------------------------------------------------------------------
 * Destroy a postfix expression
 *
 * parameters:
 *   pf <-> pointer to postfix structure pointer
 *----------------------------------------------------------------------------*/

void
fvm_selector_postfix_destroy(fvm_selector_postfix_t  **postfix);

/*----------------------------------------------------------------------------
 * Return a pointer to the infix string associated with a postfix expression
 *
 * parameters:
 *   pf <-- pointer to postfix structure
 *
 * returns:
 *   pointer to original infix string
 *----------------------------------------------------------------------------*/

const char *
fvm_selector_postfix_get_infix(const fvm_selector_postfix_t  *pf);

/*----------------------------------------------------------------------------
 * Indicate if a postfix expression depends on coordinates
 *
 * parameters:
 *   pf <-- pointer to postfix structure
 *
 * returns:
 *   true if expression depends on coordinates, false otherwise
 *----------------------------------------------------------------------------*/

bool
fvm_selector_postfix_coords_dep(const fvm_selector_postfix_t  *pf);

/*----------------------------------------------------------------------------
 * Indicate if a postfix expression depends on normals
 *
 * parameters:
 *   pf <-- pointer to postfix structure
 *
 * returns:
 *   true if expression depends on normals, false otherwise
 *----------------------------------------------------------------------------*/

bool
fvm_selector_postfix_normals_dep(const fvm_selector_postfix_t  *pf);

/*----------------------------------------------------------------------------
 * Return the number of operands associated with a postfix expression
 * missing in the associated group class set
 *
 * parameters:
 *   pf <-- pointer to postfix structure
 *
 * returns:
 *   number of missing operands
 *----------------------------------------------------------------------------*/

int
fvm_selector_postfix_n_missing(const fvm_selector_postfix_t  *pf);

/*----------------------------------------------------------------------------
 * Return a pointer to the name of an of operand associated with a postfix
 * expression but missing in the associated group class set
 *
 * parameters:
 *   pf <-- pointer to postfix structure
 *   id <-- id of missing operand (0 to fvm_selector_postfix_n_missing())
 *
 * returns:
 *   pointer to name of missing operand
 *----------------------------------------------------------------------------*/

const char *
fvm_selector_postfix_get_missing(const fvm_selector_postfix_t  *pf,
                                 int                            id);

/*----------------------------------------------------------------------------
 * Evaluate a postfix expression
 *
 * parameters:
 *   pf           <-- pointer to postfix structure
 *   n_groups     <-- number of groups associated with group class
 *   n_attributes <-- number of attributes associated with group class
 *   group_id     <-- array group ids associated with group class
 *   attribute_id <-- array of attribute ids associated with group class
 *   coords       <-- coordinates associated with evaluation, or NULL
 *   normal       <-- normal associated with evaluation, or NULL
 *
 * returns:
 *   true or false base on expression evaluation
 *----------------------------------------------------------------------------*/

bool
fvm_selector_postfix_eval(const fvm_selector_postfix_t  *pf,
                          int                            n_groups,
                          int                            n_attributes,
                          const int                      group_id[],
                          const int                      attribute_id[],
                          const double                   coords[],
                          const double                   normal[]);

/*----------------------------------------------------------------------------
 * Dump the contents of a postfix structure in human readable form
 *
 * parameters:
 *   pf           <-> pointer to postfix structure
 *   n_groups     <-- number of groups
 *   n_attributes <-- number of attributes
 *   group_name   <-- array group names (sorted)
 *   attribute    <-- array of attribute numbers (sorted)
 *----------------------------------------------------------------------------*/

void
fvm_selector_postfix_dump(const fvm_selector_postfix_t  *pf,
                          int                            n_groups,
                          int                            n_attributes,
                          const char                    *group_name[],
                          const int                      attribute[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_SELECTOR_POSTFIX_H__ */
