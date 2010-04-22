#ifndef __FVM_IO_NUM_H__
#define __FVM_IO_NUM_H__

/*============================================================================
 * Main structure for an I/O numbering scheme associated with mesh entities
 * (such as cells, faces, and vertices);
 *
 * In parallel mode, such a scheme is important so as to redistribute
 * locally numbered entities on n processes to files written by p
 * processes, with p <= n.
 *
 * Only the case where p = 1 is presently implemented, so the numbering
 * scheme is simply based on entity's global labels.
 *
 * For p > 1, it would probably be necessary to extend the numbering
 * schemes so as to account for the fact that a given entity may have
 * a main index on its main associated domain, but may be present
 * as a ghost entity with another index on neighboring domains.
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

/*----------------------------------------------------------------------------
 * Structure defining an I/O numbering scheme
 *----------------------------------------------------------------------------*/

/*
  Pointer to an I/O numbering scheme structure. The structure
  itself is private, and is defined in fvm_io_num.c
*/

typedef struct _fvm_io_num_t fvm_io_num_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   parent_entity_number <-- pointer to list of selected entitie's parent's
 *                            numbers, or NULL if all first nb_ent entities
 *                            are used
 *   parent_global_number <-- pointer to list of global (i.e. domain splitting
 *                            independent) parent entity numbers
 *   n_entities           <-- number of entities considered
 *   share_parent_global  <-- if non zero, try to share parent_global_number
 *                            instead of using a local copy
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create(const fvm_lnum_t  parent_entity_number[],
                  const fvm_gnum_t  parent_global_number[],
                  const size_t      n_entities,
                  const int         share_parent_global);

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure,
 * sharing a given global numbering array.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   global_number <-- pointer to list of global (i.e. domain splitting
 *                     independent) entity numbers
 *   global_count  <-- global number of entities
 *   n_entities    <-- number of local entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_shared(const fvm_gnum_t  global_number[],
                         fvm_gnum_t        global_count,
                         size_t            n_entities);

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on an an initial
 * I/O numbering and a number of new entities per base entity.
 *
 * This is useful for example to create an I/O numbering for
 * triangles based on split polygons, whose I/O numbering is defined.
 *
 * parameters:
 *   base_io_num    <-- pointer to base I/O numbering structure
 *   n_sub_entities <-- number of new entities per base entity
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_sub(const fvm_io_num_t  *base_io_num,
                           const fvm_lnum_t     n_sub_entities[]);

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on a strided adjacency.
 *
 * The corresponding entities must be locally ordered.
 *
 * parameters:
 *   parent_entity_number <-- pointer to list of selected entitie's parent's
 *                            numbers, or NULL if all first n_ent entities
 *                            are used
 *   adjacency            <-- entity adjacency (1 to n global numbering)
 *   n_entities           <-- number of entities considered
 *   stride               <-- values per entity
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_adj_s(const fvm_lnum_t  parent_entity_number[],
                             const fvm_gnum_t  adjacency[],
                             size_t            n_entities,
                             size_t            stride);

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on an indexed adjacency.
 *
 * The corresponding entities do not need to be locally ordered.
 *
 * parameters:
 *  parent_entity_number <-- pointer to list of selected entitie's parent's
 *                           numbers, or NULL if all first n_ent entities
 *                           are used
 *  index                <-- index on entities for adjacency
 *  adjacency            <-- entity adjacency (1 to n global numbering)
 *  n_entities           <-- number of entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_adj_i(const fvm_lnum_t  parent_entity_number[],
                             const fvm_lnum_t  index[],
                             const fvm_gnum_t  adjacency[],
                             fvm_lnum_t        n_entities);

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on coordinates.
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be determined by lexicographical ordering of coordinates.
 *
 * parameters:
 *   coords     <-- pointer to entity coordinates (interlaced)
 *   dim        <-- spatial dimension
 *   n_entities <-- number of entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_coords(const fvm_coord_t  coords[],
                              int                dim,
                              size_t             n_entities);

/*----------------------------------------------------------------------------
 * Creation of an I/O numbering structure based on a simple accumulation
 * (i.e. scan) of counts on successive ranks.
 *
 * parameters:
 *   n_entities <-- number of entities considered
 *
 * returns:
 *  pointer to I/O numbering structure
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_create_from_scan(size_t  n_entities);

/*----------------------------------------------------------------------------
 * Destruction of a I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_io_num_t *
fvm_io_num_destroy(fvm_io_num_t  * this_io_num);

/*----------------------------------------------------------------------------
 * Return local number of entities associated with an I/O numbering
 * structure.
 *
 * parameters:
 *   this_io_num <-- pointer to I/O/ numbering structure
 *
 * returns:
 *  local number of associated entities
 *----------------------------------------------------------------------------*/

fvm_lnum_t
fvm_io_num_get_local_count(const fvm_io_num_t  *const this_io_num);

/*----------------------------------------------------------------------------
 * Return global number of entities associated with an I/O numbering
 * structure.
 *
 * parameters:
 *   this_io_num <-- pointer to I/O/ numbering structure
 *
 * returns:
 *  global number of associated entities
 *----------------------------------------------------------------------------*/

fvm_gnum_t
fvm_io_num_get_global_count(const fvm_io_num_t  *const this_io_num);

/*----------------------------------------------------------------------------
 * Return global numbering associated with an I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to I/O/ numbering structure
 *
 * returns:
 *  pointer to array of global numbers associated with local entities
 *  (1 to n numbering)
 *----------------------------------------------------------------------------*/

const fvm_gnum_t *
fvm_io_num_get_global_num(const fvm_io_num_t  *const this_io_num);

/*----------------------------------------------------------------------------
 * Return the global number of sub-entities associated with an initial
 * entity whose global numbering is known, given the number of
 * sub-entities per initial entity.
 *
 * parameters:
 *   this_io_num    <-- pointer to base io numbering
 *   n_sub_entities <-- number of sub-entities per initial entity
 *   comm           <-- associated MPI communicator
 *
 * returns:
 *   global number of sub-entities
 *----------------------------------------------------------------------------*/

fvm_gnum_t
fvm_io_num_global_sub_size(const fvm_io_num_t  *this_io_num,
                           const fvm_lnum_t     n_sub_entities[]);

/*----------------------------------------------------------------------------
 * Dump printout of a I/O numbering structure.
 *
 * parameters:
 *   this_io_num <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvm_io_num_dump(const fvm_io_num_t  *const this_io_num);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_IO_NUM_H__ */
