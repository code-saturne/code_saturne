#ifndef __FVM_WRITER_HELPER_H__
#define __FVM_WRITER_HELPER_H__

/*============================================================================
 * Helper types and functions for mesh and field writers
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

#include "cs_block_dist.h"
#include "cs_part_to_block.h"
#include "fvm_defs.h"
#include "fvm_nodal.h"
#include "fvm_writer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * FVM nodal to writer section translation list
 *----------------------------------------------------------------------------*/

typedef struct _fvm_writer_section_t {

  struct _fvm_writer_section_t  *next;  /* Pointer to next element
                                           in list (NULL at end) */

  const fvm_nodal_section_t  *section;  /* Corresponding section in mesh */

  cs_gnum_t   extra_vertex_base;        /* Start global number of added
                                           vertices (for tesselation) */

  cs_lnum_t   num_shift;                /* Element number shift when no
                                           parent lists are used */
  fvm_element_t  type;                  /* Corresponding element type (may
                                           differ from  section->type when
                                           using tesselations) */

  bool    continues_previous;           /* Indicates if the corresponding FVM
                                           nodal section should be appended
                                           to the previous one on output */

} fvm_writer_section_t;

/*----------------------------------------------------------------------------
 * FVM nodal to writer field output helper
 *----------------------------------------------------------------------------*/

/*
  Pointer to structure keeping track of the status of a writer's field
  output state. The structure itself is private, and is defined in fvm_writer.c
*/

typedef struct _fvm_writer_field_helper_t fvm_writer_field_helper_t;

/*----------------------------------------------------------------------------
 * Function pointer for output of field values by a writer helper
 *
 * parameters:
 *   context      <-> pointer to writer and field context
 *   datatype     <-- output datatype
 *   dimension    <-- output field dimension
 *   component_id <-- output component id (if non-interleaved)
 *   block_start  <-- start global number of element for current block
 *   block_end    <-- past-the-end global number of element for current block
 *   buffer       <-> associated output buffer
 *----------------------------------------------------------------------------*/

typedef void
(fvm_writer_field_output_t) (void           *context,
                             cs_datatype_t   datatype,
                             int             dimension,
                             int             component_id,
                             cs_gnum_t       block_start,
                             cs_gnum_t       block_end,
                             void           *buffer);

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Build list of sections to output
 *
 * Depending on whether multiple sections of a given element type
 * are allowed or not, sections may be ordered in different ways.
 * Discarded element types are not added to the list.
 *
 * parameters:
 *   mesh                 <-- pointer to nodal mesh structure
 *   group_by_type        <-- if true, group sections of same type
 *   group_all            <-- if true, all sections continue previous ones
 *   min_export_dim       <-- minimum dimension of sections to export
 *   discard_polygons     <-- ignore polygonal sections
 *   discard_polyhedra    <-- ignore polyhedral sections
 *   divide_polygons      <-- tesselate polygonal sections
 *   divide_polyhedra     <-- tesselate polyhedral sections
 *
 * returns:
 *   array of section translations (must be freed by caller),
 *   or NULL if section list is completely empty
 *----------------------------------------------------------------------------*/

fvm_writer_section_t *
fvm_writer_export_list(const fvm_nodal_t          *mesh,
                       int                         min_export_dim,
                       bool                        group_by_type,
                       bool                        group_all,
                       bool                        discard_polygons,
                       bool                        discard_polyhedra,
                       bool                        divide_polygons,
                       bool                        divide_polyhedra);

/*----------------------------------------------------------------------------
 * Count number of extra vertices when tesselations are present
 *
 * parameters:
 *   mesh               <-- pointer to nodal mesh structure
 *   divide_polyhedra   <-- true if polyhedra are tesselated
 *   n_extra_vertices_g --> global number of extra vertices (optional)
 *   n_extra_vertices   --> local number of extra vertices (optional)
 *----------------------------------------------------------------------------*/

void
fvm_writer_count_extra_vertices(const fvm_nodal_t  *mesh,
                                bool                divide_polyhedra,
                                cs_gnum_t          *n_extra_vertices_g,
                                cs_lnum_t          *n_extra_vertices);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Build block info and part to block distribution helper for vertices.
 *
 * parameters:
 *   min_rank_step    <-- minimum step between output ranks
 *   min_block_size   <-- minimum block buffer size
 *   n_g_add_vertices <-- global number of vertices due to tesselated polyhedra
 *   n_add_vertices   <-- local number of vertices due to tesselated polyhedra
 *   mesh             <-- pointer to nodal mesh structure
 *   bi               --> block information structure
 *   d                --> part to block distributor
 *   comm             <-- associated communicator
 *----------------------------------------------------------------------------*/

void
fvm_writer_vertex_part_to_block_create(int                     min_rank_step,
                                       cs_lnum_t               min_block_size,
                                       cs_gnum_t               n_g_add_vertices,
                                       cs_lnum_t               n_add_vertices,
                                       const fvm_nodal_t      *mesh,
                                       cs_block_dist_info_t   *bi,
                                       cs_part_to_block_t    **d,
                                       MPI_Comm                comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Return extra vertex coordinates when tesselations are present
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   n_extra_vertices <-- number of extra vertices
 *
 * returns:
 *   array containing all extra vertex coordinates
 *----------------------------------------------------------------------------*/

cs_coord_t *
fvm_writer_extra_vertex_coords(const fvm_nodal_t  *mesh,
                               cs_lnum_t           n_extra_vertices);

/*----------------------------------------------------------------------------
 * Create field writer helper structure.
 *
 * Local values are initialized, ang lobal values are set to zero
 * (they may be initialized by calling fvm_writer_field_helper_init_g()).
 *
 * parameters:
 *   mesh                 <-- pointer to nodal mesh structure
 *   section_list         <-- point to export section list helper structure
 *   field_dim            <-- indicates output field dimension
 *   interlace            <-- indicates if output is interlaced
 *   location             <-- indicates if field is cell or node based
 *   datatype             <-- datatype of destination buffers
 *
 * returns:
 *   pointer to allocated and initialized field writer helper
 *----------------------------------------------------------------------------*/

fvm_writer_field_helper_t *
fvm_writer_field_helper_create(const fvm_nodal_t          *mesh,
                               const fvm_writer_section_t *section_list,
                               int                         field_dim,
                               cs_interlace_t              interlace,
                               cs_datatype_t               datatype,
                               fvm_writer_var_loc_t        location);

/*----------------------------------------------------------------------------
 * Destroy FVM field writer helper.
 *
 * parameters: *   helper <-> pointer to pointer to structure that should be destroyed
 *----------------------------------------------------------------------------*/

void
fvm_writer_field_helper_destroy(fvm_writer_field_helper_t **helper);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Set MPI info for an fvm_writer_field_helper structure.
 *
 * parameters:
 *   helper         <-> pointer to structure that should be initialized
 *   min_rank_step  <-- minimum step between output ranks
 *   min_block_size <-- minimum block buffer size
 *   comm           <-- associated MPI communicator
 *
 * returns:
 *   pointer to allocated and initialized field writer helper
 *----------------------------------------------------------------------------*/

void
fvm_writer_field_helper_init_g(fvm_writer_field_helper_t   *helper,
                               int                          min_rank_step,
                               int                          min_block_size,
                               MPI_Comm                     comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Return sizes associated with a writer field helper.
 *
 * parameters:
 *   helper                   <-- pointer to helper structure
 *   input_size               --> Total field locations in input (or NULL)
 *   output_size              --> Total field locations in output (or NULL)
 *   min_output_buffer_size   --> Minimum required buffer size (or NULL)
 *----------------------------------------------------------------------------*/

void
fvm_writer_field_helper_get_size(const fvm_writer_field_helper_t  *helper,
                                 size_t  *input_size,
                                 size_t  *output_size,
                                 size_t  *min_output_buffer_size);

/*----------------------------------------------------------------------------
 * Return the output dimension associated with a writer field helper.
 *
 * parameters:
 *   helper                   <-- pointer to helper structure
 *
 * returns:
 *   field dimension associated with helper
 *----------------------------------------------------------------------------*/

int
fvm_writer_field_helper_field_dim(const fvm_writer_field_helper_t  *helper);

/*----------------------------------------------------------------------------
 * Return the output datatype associated with a writer field helper.
 *
 * parameters:
 *   helper                   <-- pointer to helper structure
 *
 * returns:
 *   output datatype associated with helper
 *----------------------------------------------------------------------------*/

cs_datatype_t
fvm_writer_field_helper_datatype(const fvm_writer_field_helper_t  *helper);

/*----------------------------------------------------------------------------
 * Partially distribute field values to a local output buffer.
 *
 * parameters:
 *   helper             <-> pointer to helper structure
 *   export_section     <-- pointer to section helper structure
 *   src_dim            <-- dimension of source data
 *   src_dim_shift      <-- source data dimension shift (start index)
 *   src_interlace      <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- indicates the data type of (source) field values
 *   field_values       <-- array of associated field value arrays
 *   datatype           <-- input data type
 *   field_values       <-- pointer to input array
 *   output_buffer      <-- pointer to output buffer
 *                          (working array only for ranks > 0)
 *   output_buffer_size <-- size of output buffer (in datatype units)
 *   output_size        --> size of output upon return (in datatype units)
 *
 * returns:
 *   0 if values were distributed to the output buffer, 1 if the end of the
 *   section has already been reached and no values were left to distribute.
 *----------------------------------------------------------------------------*/

int
fvm_writer_field_helper_step_el(fvm_writer_field_helper_t   *helper,
                                const fvm_writer_section_t  *export_section,
                                int                          src_dim,
                                int                          src_dim_shift,
                                cs_interlace_t               src_interlace,
                                int                          n_parent_lists,
                                const cs_lnum_t              parent_num_shift[],
                                cs_datatype_t                datatype,
                                const void            *const field_values[],
                                void                        *output_buffer,
                                size_t                       output_buffer_size,
                                size_t                      *output_size);

/*----------------------------------------------------------------------------
 * Partially distribute per node field values to a local output buffer.
 *
 * parameters:
 *   helper             <-> pointer to helper structure
 *   mesh               <-- pointer to associated mesh
 *   src_dim            <-- dimension of source data
 *   src_dim_shift      <-- source data dimension shift (start index)
 *   src_interlace      <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- indicates the data type of (source) field values
 *   field_values       <-- array of associated field value arrays
 *   datatype           <-- input data type
 *   field_values       <-- pointer to input array
 *   output_buffer      <-- pointer to output buffer
 *                          (working array only for ranks > 0)
 *   output_buffer_size <-- size of output buffer (in datatype units)
 *   output_size        --> size of output upon return (in datatype units)
 *
 * returns:
 *   0 if values were distributed to the output buffer, 1 if the end of the
 *   section has already been reached and no values were left to distribute.
 *----------------------------------------------------------------------------*/

int
fvm_writer_field_helper_step_nl(fvm_writer_field_helper_t   *helper,
                                const fvm_nodal_t           *mesh,
                                int                          src_dim,
                                int                          src_dim_shift,
                                cs_interlace_t               src_interlace,
                                int                          n_parent_lists,
                                const cs_lnum_t              parent_num_shift[],
                                cs_datatype_t                datatype,
                                const void            *const field_values[],
                                void                        *output_buffer,
                                size_t                       output_buffer_size,
                                size_t                      *output_size);

/*----------------------------------------------------------------------------
 * Output per element field values, using writer-specific function
 * and context structure pointers.
 *
 * Note that if the output data is not interleaved, for multidimensional data,
 * the output function is called once per component, using the same buffer.
 * This is a good fit for most options, but if a format requires writing
 * additional buffering may be required in the context.
 *
 * parameters:
 *   helper           <-> pointer to helper structure
 *   context          <-> pointer to writer context
 *   export_section   <-- pointer to section helper structure
 *   src_dim          <-- dimension of source data
 *   src_interlace    <-- indicates if field in memory is interlaced
 *   comp_order       <-- field component reordering array, or NULL
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   field_values     <-- array of associated field value arrays
 *   output_func      <-- pointer to output function
 *
 * returns:
 *   pointer to next section helper structure in list
 *----------------------------------------------------------------------------*/

const fvm_writer_section_t *
fvm_writer_field_helper_output_e(fvm_writer_field_helper_t   *helper,
                                 void                        *context,
                                 const fvm_writer_section_t  *export_section,
                                 int                          src_dim,
                                 cs_interlace_t               src_interlace,
                                 const int                   *comp_order,
                                 int                          n_parent_lists,
                                 const cs_lnum_t              parent_num_shift[],
                                 cs_datatype_t                datatype,
                                 const void            *const field_values[],
                                 fvm_writer_field_output_t   *output_func);

/*----------------------------------------------------------------------------
 * Output per node field values, using writer-specific function
 * and context structure pointers.
 *
 * Note that if the output data is not interleaved, for multidimensional data,
 * the output function is called once per component, using the same buffer.
 * This is a good fit for most options, but if a format requires writing
 * additional buffering may be required in the context.
 *
 * parameters:
 *   helper             <-> pointer to helper structure
 *   context            <-> pointer to writer context
 *   mesh               <-- pointer to nodal mesh
 *   divide_polyhedra   <-- tesselate polyhedral sections
 *   src_dim            <-- dimension of source data
 *   src_dim_shift      <-- source data dimension shift (start index)
 *   src_interlace      <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- indicates the data type of (source) field values
 *   field_values       <-- array of associated field value arrays
 *   datatype           <-- input data type
 *   field_values       <-- pointer to input array
 *   output_func        <-- pointer to output function
 *----------------------------------------------------------------------------*/

void
fvm_writer_field_helper_output_n(fvm_writer_field_helper_t  *helper,
                                 void                       *context,
                                 const fvm_nodal_t          *mesh,
                                 int                         src_dim,
                                 cs_interlace_t              src_interlace,
                                 const int                  *comp_order,
                                 int                         n_parent_lists,
                                 const cs_lnum_t             parent_num_shift[],
                                 cs_datatype_t               datatype,
                                 const void           *const field_values[],
                                 fvm_writer_field_output_t  *output_func);

/*----------------------------------------------------------------------------
 * Set string representing a field component's name based on its id.
 *
 * parameters:
 *   s                  --> destination string
 *   s_size             <-- maximum string size
 *   lowercase          <-- true if lowercase is required
 *   dimension          <-- field dimension
 *   component_id       <-- field component id
 *----------------------------------------------------------------------------*/

void
fvm_writer_field_component_name(char    *s,
                                size_t   s_size,
                                bool     lowercase,
                                int      dimension,
                                int      component_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_WRITER_HELPER_H__ */
