#ifndef __FVM_GROUP_H__
#define __FVM_GROUP_H__

/*============================================================================
 * Definition of entity groups
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2007  EDF

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
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining (equivalence) classes of groups
 *----------------------------------------------------------------------------*/

/*
  Pointer to opaque group class and group class set objects.
  The structure contents are private, and are defined in fvm_group.c
*/

typedef struct _fvm_group_class_t     fvm_group_class_t;
typedef struct _fvm_group_class_set_t fvm_group_class_set_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the number of groups of a group class.
 *
 * parameters:
 *   this_group_class <-- pointer to group class structure
 *
 * returns:
 *   number of groups in group class
 *----------------------------------------------------------------------------*/

int
fvm_group_class_get_n_groups(const fvm_group_class_t  *this_group_class);

/*----------------------------------------------------------------------------
 * Return the array of group names of a group class.
 *
 * parameters:
 *   this_group_class <-- pointer to group class structure
 *
 * returns:
 *   pointer to array of group names in group class
 *----------------------------------------------------------------------------*/

const char **
fvm_group_class_get_group_names(const fvm_group_class_t  *this_group_class);

/*----------------------------------------------------------------------------
 * Return the number of attributes of a group class.
 *
 * parameters:
 *   this_group_class <-- pointer to group class structure
 *
 * returns:
 *   number of attributes in group class
 *----------------------------------------------------------------------------*/

int
fvm_group_class_get_n_attributes(const fvm_group_class_t  *this_group_class);

/*----------------------------------------------------------------------------
 * Return the array of attributes of a group class.
 *
 * parameters:
 *   this_group_class <-- pointer to group class structure
 *
 * returns:
 *   pointer to array of attributes in group class
 *----------------------------------------------------------------------------*/

const int *
fvm_group_class_get_attributes(const fvm_group_class_t  *this_group_class);

/*----------------------------------------------------------------------------
 * Creation of a group class set structure.
 *
 * returns:
 *   pointer to created group class set structure
 *----------------------------------------------------------------------------*/

fvm_group_class_set_t *
fvm_group_class_set_create(void);

/*----------------------------------------------------------------------------
 * Add a group class to a set.
 *
 * parameters:
 *   this_group_class_set <-> pointer to group class set structure
 *   n_groups             <-- number of groups in class
 *   group_names          <-- array of group names
 *   n_attributes         <-- number of attributes in class
 *   attributes           <-- list of attributes
 *----------------------------------------------------------------------------*/

void
fvm_group_class_set_add(fvm_group_class_set_t   *this_group_class_set,
                        int                      n_groups,
                        int                      n_attributes,
                        const char             **group_names,
                        const int               *attributes);

/*----------------------------------------------------------------------------
 * Destruction of a group class set structure.
 *
 * parameters:
 *   this_class_set <-- pointer to structure which should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_group_class_set_t *
fvm_group_class_set_destroy(fvm_group_class_set_t  *this_class_set);

/*----------------------------------------------------------------------------
 * Return number of classes in a group class set.
 *
 * parameters:
 *   this_group_class_set <-> pointer to group class set structure
 *
 * returns:
 *   number of classes in a group class set
 *----------------------------------------------------------------------------*/

int
fvm_group_class_set_size(const fvm_group_class_set_t  *this_group_class_set);

/*----------------------------------------------------------------------------
 * Return pointer to a given group class in a group class set.
 *
 * parameters:
 *   this_group_class_set <-- pointer to group class set structure
 *   group_class_id       <-- index of group class in set (0 to n-1)
 *
 * returns:
 *   pointer to group class structure
 *----------------------------------------------------------------------------*/

const fvm_group_class_t *
fvm_group_class_set_get(const fvm_group_class_set_t  *this_group_class_set,
                        int                           group_class_id);

/*----------------------------------------------------------------------------
 * Dump printout of a group class set
 *
 * parameters:
 *   this_class_set      <-- pointer to group class set to be dumped
 *----------------------------------------------------------------------------*/

void
fvm_group_class_set_dump(const fvm_group_class_set_t  *this_group_class_set);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_GROUP_H__ */
