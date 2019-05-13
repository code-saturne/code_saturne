#ifndef __FVM_NODAL_H__
#define __FVM_NODAL_H__

/*============================================================================
 * Main structure for a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_mesh.h"
#include "fvm_defs.h"
#include "fvm_group.h"
#include "fvm_io_num.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining a mesh in nodal definition
 *----------------------------------------------------------------------------*/

typedef struct _fvm_nodal_t fvm_nodal_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Number of vertices associated with each "nodal" element type */

extern const int  fvm_nodal_n_vertices_element[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a nodal mesh representation structure.
 *
 * parameters:
 *   name <-- name that should be assigned to the nodal mesh
 *   dim  <-- spatial dimension
 *
 * returns:
 *  pointer to created nodal mesh representation structure
 *----------------------------------------------------------------------------*/

fvm_nodal_t *
fvm_nodal_create(const char  *name,
                 int          dim);

/*----------------------------------------------------------------------------
 * Destruction of a nodal mesh representation structure.
 *
 * parameters:
 *   this_nodal  <-> pointer to structure that should be destroyed
 *
 * returns:
 *  NULL pointer
 *----------------------------------------------------------------------------*/

fvm_nodal_t *
fvm_nodal_destroy(fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Copy a nodal mesh representation structure, sharing arrays with the
 * original structure.
 *
 * parameters:
 *   this_nodal  <-> pointer to structure that should be copied
 *
 * returns:
 *   pointer to created nodal mesh representation structure
 *----------------------------------------------------------------------------*/

fvm_nodal_t *
fvm_nodal_copy(const fvm_nodal_t *this_nodal);

/*----------------------------------------------------------------------------
 * Reduction of a nodal mesh representation structure: only the associations
 * (numberings) necessary to redistribution of fields for output are
 * conserved, the full connectivity being in many cases no longer useful
 * once it has been output. If the del_vertex_num value is set
 * to true, vertex-based values may not be output in parallel mode
 * after this function is called.
 *
 * parameters:
 *   this_nodal        <-> pointer to structure that should be reduced
 *   del_vertex_num    <-- indicates if vertex parent indirection and
 *                         I/O numbering are destroyed (1) or not (0)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_reduce(fvm_nodal_t  *this_nodal,
                 int           del_vertex_num);

/*----------------------------------------------------------------------------
 * Change entity parent numbering; this is useful when entities of the
 * parent mesh have been renumbered after a nodal mesh representation
 * structure's creation.
 *
 * parameters:
 *   this_nodal          <-- nodal mesh structure
 *   new_parent_num      <-- pointer to local parent renumbering array
 *                           ({1, ..., n} <-- {1, ..., n})
 *   entity_dim          <-- 3 for cells, 2 for faces, 1 for edges,
 *                           and 0 for vertices
 *----------------------------------------------------------------------------*/

void
fvm_nodal_change_parent_num(fvm_nodal_t        *this_nodal,
                            const cs_lnum_t     new_parent_num[],
                            int                 entity_dim);

/*----------------------------------------------------------------------------
 * Remove entity parent numbering; this is useful for example when we
 * want to assign coordinates or fields to an extracted mesh using
 * arrays relative to the mesh, and not to its parent.
 *
 * This is equivalent to calling fvm_nodal_change_parent_num(), with
 * 'trivial' (1 o n) new_parent_num[] values.
 *
 * parameters:
 *   this_nodal          <-- nodal mesh structure
 *   entity_dim          <-- 3 for cells, 2 for faces, 1 for edges,
 *                           and 0 for vertices
 *----------------------------------------------------------------------------*/

void
fvm_nodal_remove_parent_num(fvm_nodal_t  *this_nodal,
                            int           entity_dim);

/*----------------------------------------------------------------------------
 * Build external numbering for entities based on global numbers.
 *
 * parameters:
 *   this_nodal           <-- nodal mesh structure
 *   parent_global_number <-- pointer to list of global (i.e. domain splitting
 *                            independent) parent entity numbers
 *   entity_dim           <-- 3 for cells, 2 for faces, 1 for edges,
 *                            and 0 for vertices
 *----------------------------------------------------------------------------*/

void
fvm_nodal_init_io_num(fvm_nodal_t        *this_nodal,
                      const cs_gnum_t     parent_global_numbers[],
                      int                 entity_dim);

/*----------------------------------------------------------------------------
 * Transfer an existing global numbering structure to a nodal mesh.
 *
 * This assumes the local number of vertices is identical to that of
 * the numberign structure.
 *
 * The argument pointer to the global numbering structure is set to NULL
 * so the caller does not spuriously modify the structure after this call.
 *
 * A call to this function may be used instead of fvm_nodal_init_io_num
 * (for vertices only).
 *
 * parameters:
 *   this_nodal    <-- nodal mesh structure
 *   io_num        <-> pointer to external numbering structure
 *                     whose property is transferred to the mesh
 *                     independent) parent entity numbers
 *----------------------------------------------------------------------------*/

void
fvm_nodal_transfer_vertex_io_num(fvm_nodal_t    *this_nodal,
                                 fvm_io_num_t  **io_num);

/*----------------------------------------------------------------------------
 * Set entity tags (for non-vertex entities).
 *
 * The number of entities of the given dimension may be obtained
 * through fvm_nodal_get_n_entities(), the tag[] array is populated
 * in local section order, section by section).
 *
 * parameters:
 *   this_nodal <-- nodal mesh structure
 *   tag        <-- tag values to assign
 *   entity_dim <-- 3 for cells, 2 for faces, 1 for edges
 *----------------------------------------------------------------------------*/

void
fvm_nodal_set_tag(fvm_nodal_t  *this_nodal,
                  const int     tag[],
                  int           entity_dim);

/*----------------------------------------------------------------------------
 * Remove entity tags.
 *
 * parameters:
 *   this_nodal <-- nodal mesh structure
 *   entity_dim <-- 3 for cells, 2 for faces, 1 for edges
 *----------------------------------------------------------------------------*/

void
fvm_nodal_remove_tag(fvm_nodal_t  *this_nodal,
                     int           entity_dim);

/*----------------------------------------------------------------------------
 * Preset number and list of vertices to assign to a nodal mesh.
 *
 * If the parent_vertex_num argument is NULL, the list is assumed to
 * be {1, 2, ..., n}. If parent_vertex_num is given, it specifies a
 * list of n vertices from a larger set (1 to n numbering).
 *
 * Ownership of the given parent vertex numbering array is
 * transferred to the nodal mesh representation structure.
 *
 * This function should be called before fvm_nodal_set_shared_vertices()
 * or fvm_nodal_transfer_vertices() if we want to force certain
 * vertices to appear in the mesh (especially if we want to define
 * a mesh containing only vertices).
 *
 * parameters:
 *   this_nodal        <-> nodal mesh structure
 *   n_vertices        <-- number of vertices to assign
 *   parent_vertex_num <-- parent numbers of vertices to assign
 *----------------------------------------------------------------------------*/

void
fvm_nodal_define_vertex_list(fvm_nodal_t  *this_nodal,
                             cs_lnum_t     n_vertices,
                             cs_lnum_t     parent_vertex_num[]);

/*----------------------------------------------------------------------------
 * Assign shared vertex coordinates to an extracted nodal mesh,
 * renumbering vertex numbers based on those really referenced,
 * and updating connectivity arrays in accordance.
 *
 * This function should be called once all element sections have
 * been added to a nodal mesh representation.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   vertex_coords   <-- coordinates of parent vertices (interlaced)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_set_shared_vertices(fvm_nodal_t       *this_nodal,
                              const cs_coord_t   vertex_coords[]);

/*----------------------------------------------------------------------------
 * Assign private vertex coordinates to a nodal mesh,
 * renumbering vertex numbers based on those really referenced,
 * and updating connectivity arrays in accordance.
 *
 * Ownership of the given coordinates array is transferred to
 * the nodal mesh representation structure.
 *
 * This function should only be called once all element sections
 * have been added to a nodal mesh representation.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   vertex_coords   <-- coordinates of parent vertices (interlaced)
 *
 * returns:
 *   updated pointer to vertex_coords (may be different from initial
 *   argument if vertices were renumbered).
 *----------------------------------------------------------------------------*/

cs_coord_t *
fvm_nodal_transfer_vertices(fvm_nodal_t  *this_nodal,
                            cs_coord_t    vertex_coords[]);

/*----------------------------------------------------------------------------
 * Make vertex coordinates of a nodal mesh private.
 *
 * If vertex coordinates were previously shared, those coordinates that
 * are actually references are copied, and the relation to parent vertices
 * is discarded.
 *
 * If vertices were already private, the mesh is not modified.
 *
 * parameters:
 *   this_nodal <-> nodal mesh structure
 *----------------------------------------------------------------------------*/

void
fvm_nodal_make_vertices_private(fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Assign group class set descriptions to a nodal mesh.
 *
 * The structure builds its own copy of the group class sets,
 * renumbering them so as to discard those not referenced.
 * Empty group classes are also renumbered to zero.
 *
 * This function should only be called once all element sections
 * have been added to a nodal mesh representation.
 *
 * parameters:
 *   this_nodal <-> nodal mesh structure
 *   gc_set     <-- group class set descriptions
 *----------------------------------------------------------------------------*/

void
fvm_nodal_set_group_class_set(fvm_nodal_t                  *this_nodal,
                              const fvm_group_class_set_t  *gc_set);

/*----------------------------------------------------------------------------
 * Assign global vertex labels to a nodal mesh.
 *
 * As these are expected to be used only for small sets (i.e. probes)
 * where the point set is built from a global definition and data movement
 * would add complexity and overhead, the labels refer to a global view
 * on rank 0.
 *
 * The size of the labels pointers array should be the same as that
 * returned by fvm_nodal_n_g_vertices();
 *
 * This function should only be called once the nodal mesh representation
 * has been completed, as most functions modifying its vertex definitions
 * will remove these labels.
 *
 * parameters:
 *   this_nodal <-> nodal mesh structure
 *   g_labels   <-- global vertex labels, or NULL
 *----------------------------------------------------------------------------*/

void
fvm_nodal_transfer_global_vertex_labels(fvm_nodal_t  *this_nodal,
                                        char         *g_labels[]);

/*----------------------------------------------------------------------------
 * Obtain the name of a nodal mesh.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *
 * returns:
 *   pointer to constant string containing the mesh name
 *----------------------------------------------------------------------------*/

const char *
fvm_nodal_get_name(const fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Return spatial dimension of the nodal mesh.
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *
 * returns:
 *  spatial dimension
 *----------------------------------------------------------------------------*/

int
fvm_nodal_get_dim(const fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Return maximum dimension of entities in a nodal mesh.
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *
 * returns:
 *  maximum dimension of entities in mesh (0 to 3)
 *----------------------------------------------------------------------------*/

int
fvm_nodal_get_max_entity_dim(const fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Return number of entities of a given dimension in a nodal mesh.
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *   entity_dim <-- dimension of entities we want to count (0 to 3)
 *
 * returns:
 *  number of entities of given dimension in mesh
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_nodal_get_n_entities(const fvm_nodal_t  *this_nodal,
                         int                 entity_dim);

/*----------------------------------------------------------------------------
 * Return global number of vertices associated with nodal mesh.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *
 * returns:
 *   global number of vertices associated with nodal mesh
 *----------------------------------------------------------------------------*/

cs_gnum_t
fvm_nodal_get_n_g_vertices(const fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Return global number of elements of a given type associated with nodal mesh.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   element_type         <-- type of elements for query
 *
 * returns:
 *   global number of elements of the given type associated with nodal mesh
 *----------------------------------------------------------------------------*/

cs_gnum_t
fvm_nodal_get_n_g_elements(const fvm_nodal_t  *this_nodal,
                           fvm_element_t       element_type);

/*----------------------------------------------------------------------------
 * Return local number of elements of a given type associated with nodal mesh.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *   element_type         <-- type of elements for query
 *
 * returns:
 *   local number of elements of the given type associated with nodal mesh
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_nodal_get_n_elements(const fvm_nodal_t  *this_nodal,
                         fvm_element_t       element_type);

/*----------------------------------------------------------------------------
 * Return local parent numbering array for all entities of a given
 * dimension in a nodal mesh.
 *
 * The number of entities of the given dimension may be obtained
 * through fvm_nodal_get_n_entities(), the parent_num[] array is populated
 * with the parent entity numbers of those entities, in order (i.e. in
 * local section order, section by section).
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *   entity_dim <-- dimension of entities we are interested in (0 to 3)
 *   parent_num --> entity parent numbering (array must be pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_parent_num(const fvm_nodal_t  *this_nodal,
                         int                 entity_dim,
                         cs_lnum_t           parent_num[]);

/*----------------------------------------------------------------------------
 * Return pointer to global vertex labels of a nodal mesh.
 *
 * As these are expected to be used only for small sets (i.e. probes)
 * where the point set is built from a global definition and data movement
 * would adds complexity and overhead, the labels refer to a global view
 * on rank 0; for the same reason, only shared labels are needed.
 *
 * The size of the labels pointers array returned should be the same
 * as that returned by fvm_nodal_n_g_vertices();
 *
 * parameters:
 *   this_nodal <-> nodal mesh structure
 *
 * returns:
 *   pointer to global vertex labels, or NULL
 *----------------------------------------------------------------------------*/

const char **
fvm_nodal_get_global_vertex_labels(const fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Return a const pointer to a nodal mesh representation's parent mesh.
 *
 * parameters:
 *   this_nodal <-> pointer to structure that should be reduced
 *
 * return:
 *   const pointer to parent mesh
 *----------------------------------------------------------------------------*/

const cs_mesh_t  *
fvm_nodal_get_parent(const fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Associate a parent mesh to a nodal mesh representation structure.
 *
 * parameters:
 *   this_nodal <-> pointer to structure that should be reduced
 *   parent     <-- const pointer to parent mesh
 *----------------------------------------------------------------------------*/

void
fvm_nodal_set_parent(fvm_nodal_t      *this_nodal,
                     const cs_mesh_t  *parent);

/*----------------------------------------------------------------------------
 * Compute tesselation a a nodal mesh's sections of a given type, and add the
 * corresponding structure to the mesh representation.
 *
 * If global element numbers are used (i.e. in parallel mode), this function
 * should be only be used after calling fvm_nodal_init_io_num().
 *
 * If some mesh sections have already been tesselated, their tesselation
 * is unchanged.
 *
 * parameters:
 *   this_nodal  <-> pointer to nodal mesh structure
 *   type        <-> element type that should be tesselated
 *   error_count --> number of elements with a tesselation error
 *                   counter (optional)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_tesselate(fvm_nodal_t    *this_nodal,
                    fvm_element_t   type,
                    cs_lnum_t      *error_count);

/*----------------------------------------------------------------------------
 * Build a nodal representation structure based on extraction of a
 * mesh's edges.
 *
 * parameters:
 *   name        <-- name to assign to extracted mesh
 *   this_nodal  <-> pointer to nodal mesh structure
 *----------------------------------------------------------------------------*/

fvm_nodal_t *
fvm_nodal_copy_edges(const char         *name,
                     const fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Dump printout of a nodal representation structure.
 *
 * parameters:
 *   this_nodal <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvm_nodal_dump(const fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_NODAL_H__ */
