#ifndef __FVM_GROUP_H__
#define __FVM_GROUP_H__

/*============================================================================
 * Definition of entity groups
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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 *----------------------------------------------------------------------------*/

void
fvm_group_class_set_add(fvm_group_class_set_t   *this_group_class_set,
                        int                      n_groups,
                        const char             **group_names);

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
 * Copy a group class set, optionally with a renumbering array.
 *
 * parameters:
 *   this_class_set <-- pointer to group class set to be copied
 *   n_renums       <-- number of group classes
 *   renum          <-- group class renumbering (1 to n), or NULL
 *
 * returns:
 *   pointer to new copy of group class set
 *----------------------------------------------------------------------------*/

fvm_group_class_set_t  *
fvm_group_class_set_copy(const fvm_group_class_set_t  *this_group_class_set,
                         int                           n_renums,
                         int                           renum[]);

/*----------------------------------------------------------------------------
 * Dump printout of a group class set
 *
 * parameters:
 *   this_class_set      <-- pointer to group class set to be dumped
 *----------------------------------------------------------------------------*/

void
fvm_group_class_set_dump(const fvm_group_class_set_t  *this_group_class_set);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_GROUP_H__ */
