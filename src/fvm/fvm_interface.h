#ifndef __FVM_INTERFACE_H__
#define __FVM_INTERFACE_H__

/*============================================================================
 * Main structure for handling of interfaces associating mesh entities
 * (such as inter-processor or periodic connectivity between cells, faces,
 * or vertices);
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2006-2007  EDF

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
#include "fvm_periodicity.h"

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

/*----------------------------------------------------------------------------
 * Structure defining an I/O numbering scheme
 *----------------------------------------------------------------------------*/

/*
  Pointer to structures representing an interface and a list of interfaces.
  The structures themselves are private, and is defined in fvm_interface.c
*/

typedef struct _fvm_interface_t     fvm_interface_t;
typedef struct _fvm_interface_set_t fvm_interface_set_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return process rank associated with an interface's distant entities.
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   process rank associated with the interface's distant entities
 *----------------------------------------------------------------------------*/

int
fvm_interface_rank(const fvm_interface_t  *this_interface);

/*----------------------------------------------------------------------------
 * Return number of local and distant entities defining an interface.
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   number of local and distant entities defining the interface
 *----------------------------------------------------------------------------*/

fvm_lnum_t
fvm_interface_size(const fvm_interface_t  *this_interface);

/*----------------------------------------------------------------------------
 * Return pointer to array of local entity numbers defining an interface.
 *
 * The size of the array may be obtained by fvm_interface_size().
 * The array is owned by the interface structure, and is not copied
 * (hence the constant qualifier for the return value).
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   pointer to array of local entity numbers defining the interface
 *----------------------------------------------------------------------------*/

const fvm_lnum_t *
fvm_interface_get_local_num(const fvm_interface_t  *this_interface);

/*----------------------------------------------------------------------------
 * Return pointer to array of distant entity numbers defining an interface.
 *
 * The size of the array may be obtained by fvm_interface_size().
 * The array is owned by the interface structure, and is not copied
 * (hence the constant qualifier for the return value).
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   pointer to array of local entity numbers defining the interface
 *----------------------------------------------------------------------------*/

const fvm_lnum_t *
fvm_interface_get_distant_num(const fvm_interface_t  *this_interface);

/*----------------------------------------------------------------------------
 * Return size of index of sub-sections for different transformations.
 *
 * The index is applicable to both local_num and distant_num arrays,
 * with purely parallel equivalences appearing at position 0, and
 * equivalences through periodic transform i at position i+1;
 * Its size should thus be equal to 1 + number of periodic transforms + 1,
 * In absence of periodicity, it may be 0, as the index is not needed.
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   transform index size for the interface
 *----------------------------------------------------------------------------*/

fvm_lnum_t
fvm_interface_get_tr_index_size(const fvm_interface_t  *this_interface);

/*----------------------------------------------------------------------------
 * Return pointer to index of sub-sections for different transformations.
 *
 * The index is applicable to both local_num and distant_num arrays,
 * with purely parallel equivalences appearing at position 0, and
 * equivalences through periodic transform i at position i+1;
 * In absence of periodicity, it may be NULL, as it is not needed.
 *
 * parameters:
 *   this_interface <-- pointer to interface structure
 *
 * returns:
 *   pointer to transform index for the interface
 *----------------------------------------------------------------------------*/

const fvm_lnum_t *
fvm_interface_get_tr_index(const fvm_interface_t  *this_interface);

/*----------------------------------------------------------------------------
 * Creation of a list of interfaces between entities of a same type.
 *
 * These interfaces may be used to identify equivalent vertices or faces using
 * domain splitting, as well as periodic entities (on the same or on
 * distant ranks).
 *
 * Note that periodicity information will be completed and made consistent
 * based on the input, so that if a periodic couple is defined on a given rank,
 * the reverse couple wil be defined, whether it is also defined on the same
 * or a different rank.
 *
 * In addition, multiple periodicity interfaces will be built automatically
 * if the periodicity structure provides for composed periodicities, so they
 * need not be defined prior to this function.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   n_ent                <-- number of local entities considered
 *                            (size of parent_entity_number[]
 *   parent_entity_number <-- pointer to list of selected entities local
 *                            numbers (1 to n), or NULL if all first n_ent
 *                            entities are used
 *   global_number        <-- pointer to list of global (i.e. domain splitting
 *                            independent) entity numbers
 *   periodicity          <-- periodicity information (NULL if none)
 *   n_periodic_lists     <-- number of periodic lists (may be local)
 *   periodicity_num      <-- periodicity number (1 to n) associated with
 *                            each periodic list (primary periodicities only)
 *   n_periodic_couples   <-- number of periodic couples associated with
 *                            each periodic list
 *   periodic_couples     <-- array indicating periodic couples (using
 *                            global numberings) for each list
 *
 * returns:
 *  pointer to list of interfaces (possibly NULL in serial mode)
 *----------------------------------------------------------------------------*/

fvm_interface_set_t *
fvm_interface_set_create(fvm_lnum_t                n_ent,
                         const fvm_lnum_t          parent_entity_number[],
                         const fvm_gnum_t          global_number[],
                         const fvm_periodicity_t  *periodicity,
                         int                       n_periodic_lists,
                         const int                 periodicity_num[],
                         const fvm_lnum_t          n_periodic_couples[],
                         const fvm_gnum_t   *const periodic_couples[]);

/*----------------------------------------------------------------------------
 * Destruction of an interface list.
 *
 * parameters:
 *   this_interface_set <-- pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_interface_set_t *
fvm_interface_set_destroy(fvm_interface_set_t  *this_interface_set);

/*----------------------------------------------------------------------------
 * Return number of interfaces associated with an interface set.
 *
 * parameters:
 *   this_interface_set <-- pointer to interface set structure
 *
 * returns:
 *   number of interfaces in set
 *----------------------------------------------------------------------------*/

int
fvm_interface_set_size(const fvm_interface_set_t  *this_interface_set);

/*----------------------------------------------------------------------------
 * Return pointer to a given interface in an interface set.
 *
 * parameters:
 *   this_interface_set <-- pointer to interface set structure
 *   interface_id       <-- index of interface in set (0 to n-1)
 *
 * returns:
 *   pointer to interface structure
 *----------------------------------------------------------------------------*/

const fvm_interface_t *
fvm_interface_set_get(const fvm_interface_set_t  *this_interface_set,
                      int                         interface_id);

/*----------------------------------------------------------------------------
 * Return pointer to the periocicity structure associated of an interface set.
 *
 * parameters:
 *   this_interface_set <-- pointer to interface set structure
 *
 * returns:
 *   pointer to periodicity structure, or NULL
 *----------------------------------------------------------------------------*/

const fvm_periodicity_t *
fvm_interface_set_periodicity(const fvm_interface_set_t  *this_interface_set);

/*----------------------------------------------------------------------------
 * Dump printout of an interface list.
 *
 * parameters:
 *   this_interface_set <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvm_interface_set_dump(const fvm_interface_set_t  *this_interface_set);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_INTERFACE_H__ */
