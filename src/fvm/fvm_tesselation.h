#ifndef __FVM_TESSELATION_H__
#define __FVM_TESSELATION_H__

/*============================================================================
 * Structure describing a mesh section's subdivision into simpler elements
 *
 * This is mostly useful to replace polygons or polyhedra by simpler
 * elements such as triangles, tetrahedra, and prisms upon data export.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include "fvm_io_num.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define FVM_TESSELATION_N_SUB_TYPES_MAX 2

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining a tesselation of a mesh section.
 *----------------------------------------------------------------------------*/

/*
  Pointer to tesselation structure. The structure
  itself is private, and is defined in fvm_tesselation.c
*/

typedef struct _fvm_tesselation_t fvm_tesselation_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a mesh section's subdivision into simpler elements.
 *
 * The structure contains pointers to the mesh section's connectivity,
 * (passed as arguments), which is not copied. This structure should thus
 * always be destroyed before the mesh section to which it relates.
 *
 * Unused connectivity array arguments should be set to NULL (such as
 * face_index[] and face_num[] for 2D or regular (strided) elements,
 * and vertex_index[] for strided elements.
 *
 * At this stage, the structure does not yet contain tesselation information.
 *
 * parameters:
 *   element_type       <-- type of elements considered
 *   n_elements         <-- number of elements
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          dimension [n_elements + 1]
 *   face_num           <-- element -> face numbers (1 to n, signed,
 *                           > 0 for outwards pointing face normal
 *                           < 0 for inwards pointing face normal);
 *                          dimension: [face_index[n_elements]], or NULL
 *   vertex_index       <-- element face -> vertices index (O to n-1);
 *                          dimension: [n_cell_faces + 1], [n_elements + 1],
 *                          or NULL depending on face_index and vertex_index
 *   vertex_num         <-- element -> vertex connectivity (1 to n)
 *   global_element_num <-- global element numbers (NULL in serial mode)
 *
 * returns:
 *  pointer to created mesh section tesselation structure
 *----------------------------------------------------------------------------*/

fvm_tesselation_t *
fvm_tesselation_create(fvm_element_t        element_type,
                       cs_lnum_t            n_elements,
                       const cs_lnum_t      face_index[],
                       const cs_lnum_t      face_num[],
                       const cs_lnum_t      vertex_index[],
                       const cs_lnum_t      vertex_num[],
                       const fvm_io_num_t  *global_element_num);

/*----------------------------------------------------------------------------
 * Destruction of a mesh section tesselation structure.
 *
 * parameters:
 *   this_tesselation <-> pointer to structure that should be destroyed
 *
 * returns:
 *  NULL pointer
 *----------------------------------------------------------------------------*/

fvm_tesselation_t *
fvm_tesselation_destroy(fvm_tesselation_t  * this_tesselation);

/*----------------------------------------------------------------------------
 * Tesselate a mesh section referred to by an fvm_tesselation_t structure.
 *
 * parameters:
 *   this_tesselation   <-> partially initialized tesselation structure
 *   dim                <-- spatial dimension
 *   vertex_coords      <-- associated vertex coordinates array
 *   parent_vertex_num  <-- optional indirection to vertex coordinates
 *   error_count        --> number of elements with a tesselation error
 *                          counter (optional)
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_init(fvm_tesselation_t  *this_tesselation,
                     int                 dim,
                     const cs_coord_t    vertex_coords[],
                     const cs_lnum_t     parent_vertex_num[],
                     cs_lnum_t          *error_count);

/*----------------------------------------------------------------------------
 * Reduction of a nodal mesh polygon splitting representation structure;
 * only the associations (numberings) necessary to redistribution of fields
 * for output are conserved, the full connectivity being no longer useful
 * once it has been output.
 *
 * parameters:
 *   this_tesselation <-> pointer to structure that should be reduced
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_reduce(fvm_tesselation_t  * this_tesselation);

/*----------------------------------------------------------------------------
 * Return number of parent elements of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   number of parent elements
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_tesselation_n_elements(const fvm_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------
 * Return global number of added vertices associated with a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   global number of added vertices associated with the tesselation
 *----------------------------------------------------------------------------*/

cs_gnum_t
fvm_tesselation_n_g_vertices_add(const fvm_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------
 * Return (local) number of added vertices associated with a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   global number of added vertices associated with the tesselation
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_tesselation_n_vertices_add(const fvm_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------
 * Return number of resulting sub-types of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   number of resulting sub-types of the tesselation
 *----------------------------------------------------------------------------*/

int
fvm_tesselation_n_sub_types(const fvm_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------
 * Return given sub-types of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   sub_type_id      <-- index of sub-type in tesselation (0 to n-1)
 *
 * returns:
 *   sub-types of the tesselation with the given index
 *----------------------------------------------------------------------------*/

fvm_element_t
fvm_tesselation_sub_type(const fvm_tesselation_t  *this_tesselation,
                         int                       sub_type_id);

/*----------------------------------------------------------------------------
 * Return number of elements of a given sub-type of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   sub_type_id      <-- index of sub-type in tesselation (0 to n-1)
 *
 * returns:
 *   sub-types of the tesselation with the given index
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_tesselation_n_sub_elements(const fvm_tesselation_t  *this_tesselation,
                               fvm_element_t             sub_type);

/*----------------------------------------------------------------------------
 * Obtain the global and maximum number of elements of a given sub-type
 * of a tesselation.
 *
 * parameters:
 *   this_tesselation    <-- tesselation structure
 *   sub_type_id         <-- index of sub-type in tesselation (0 to n-1)
 *   n_sub_elements_glob --> global number of sub-elements of the given type
 *   n_sub_elements_max  --> maximum number of sub-elements per element
 *                           of the given type (for all ranks)
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_get_global_size(const fvm_tesselation_t  *this_tesselation,
                                fvm_element_t             sub_type,
                                cs_gnum_t                *n_sub_elements_glob,
                                cs_lnum_t                *n_sub_elements_max);

/*----------------------------------------------------------------------------
 * Return global numbering of added vertices associated with a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   pointer to global numbering of added vertices for this tesselation,
 *   or NULL if no added vertices are present.
 *----------------------------------------------------------------------------*/

const fvm_io_num_t *
fvm_tesselation_global_vertex_num(const fvm_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------
 * Compute coordinates of added vertices for a tesselation of polyhedra.
 *
 * One additional vertex is added near the center of each polyhedra.
 * For element types other than polyhedra, there is no need for added
 * vertices, so this function returns immediately.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   vertex_coords      --> coordinates of added vertices
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_vertex_coords(const fvm_tesselation_t  *this_tesselation,
                              cs_coord_t                vertex_coords[]);

/*----------------------------------------------------------------------------
 * Return index of sub-elements associated with each element of a given
 * sub-type of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   sub_type_id      <-- index of sub-type in tesselation (0 to n-1)
 *
 * returns:
 *   index of sub-elements associated with each element (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

const cs_lnum_t *
fvm_tesselation_sub_elt_index(const fvm_tesselation_t  *this_tesselation,
                              fvm_element_t             sub_type);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Decode tesselation to a connectivity buffer.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   extra_vertex_base  <-- starting number for added vertices
 *   global_vertex_num  <-- global vertex numbering
 *   extra_vertex_base  <-- starting number for added vertices
 *   vertex_num         --> sub-element (global) vertex connectivity
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_decode_g(const fvm_tesselation_t  *this_tesselation,
                         fvm_element_t             connect_type,
                         const fvm_io_num_t       *global_vertex_num,
                         cs_gnum_t                 extra_vertex_base,
                         cs_gnum_t                 vertex_num[]);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Decode tesselation to a connectivity buffer.
 *
 * To avoid requiring huge buffers and computing unneeded element
 * connectivities, this function may decode a partial connectivity range,
 * starting at polygon index start_id and ending either when the indicated
 * buffer size or the last polygon is attained.
 * It returns the effective polygon index end.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   start_id           <-- start index of polygons subset in parent section
 *   buffer_limit       <-- maximum number of sub-elements of destination
 *                          element type allowable for vertex_num[] buffer
 *   extra_vertex_base  <-- starting number for added vertices
 *   vertex_num         --> sub-element (global) vertex connectivity
 *
 * returns:
 *   polygon index corresponding to end of decoded range
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_tesselation_decode(const fvm_tesselation_t  *this_tesselation,
                       fvm_element_t             connect_type,
                       cs_lnum_t                 start_id,
                       cs_lnum_t                 buffer_limit,
                       cs_lnum_t                 extra_vertex_base,
                       cs_lnum_t                 vertex_num[]);

/*----------------------------------------------------------------------------
 * Distribute "per element" data from the base elements to their tesselation.
 *
 * The same data array is used for input and output, so as to avoid requiring
 * excess allocation in typical use cases (extracting data from a parent mesh
 * to a buffer and distributing it as per its tesselation).
 * The data array should be at least of size:
 * [sub_elt_index[end_id] - sub_elt_index[start_id] * size
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   start_id           <-- start index of elements subset in parent section
 *   end_id             <-- end index of elements subset in parent section
 *   size               <-- data size for each element (sizeof(type)*stride)
 *   data               <-> undistributed data in, distributed data out
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_distribute(const fvm_tesselation_t  *this_tesselation,
                           fvm_element_t             connect_type,
                           cs_lnum_t                 start_id,
                           cs_lnum_t                 end_id,
                           size_t                    size,
                           void                     *data);

/*----------------------------------------------------------------------------
 * Compute field values at added vertices for a tesselation of polyhedra.
 *
 * One additional vertex is added near the center of each polyhedra.
 * For element types other than polyhedra, there is no need for added
 * vertices, so this function returns immediately.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   vertex_coords    <-- coordinates of added vertices
 *   src_dim          <-- dimension of source data
 *   src_dim_shift    <-- source data dimension shift (start index)
 *   dest_dim         <-- destination data dimension (1 if non interlaced)
 *   start_id         <-- added vertices start index
 *   end_id           <-- added vertices past the end index
 *   src_interlace    <-- indicates if source data is interlaced
 *   src_datatype     <-- source data type (float, double, or int)
 *   dest_datatype    <-- destination data type (float, double, or int)
 *   n_parent_lists   <-- number of parent lists (if parent_num != NULL)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   parent_num       <-- if n_parent_lists > 0, parent entity numbers
 *   src_data         <-- array of source arrays (at least one, with one per
 *                        source dimension if non interlaced, times one per
 *                        parent list if multiple parent lists, with
 *                        x_parent_1, y_parent_1, ..., x_parent_2, ...) order
 *   dest_data        --> destination buffer
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_vertex_values(const fvm_tesselation_t  *this_tesselation,
                              int                       src_dim,
                              int                       src_dim_shift,
                              int                       dest_dim,
                              cs_lnum_t                 start_id,
                              cs_lnum_t                 end_id,
                              cs_interlace_t            src_interlace,
                              cs_datatype_t             src_datatype,
                              cs_datatype_t             dest_datatype,
                              int                       n_parent_lists,
                              const cs_lnum_t           parent_num_shift[],
                              const cs_lnum_t           parent_num[],
                              const void         *const src_data[],
                              void               *const dest_data);

/*----------------------------------------------------------------------------
 * Dump printout of a mesh section tesselation structure.
 *
 * parameters:
 *   this_tesselation  <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_dump(const fvm_tesselation_t  *this_tesselation);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_TESSELATION_H__ */
