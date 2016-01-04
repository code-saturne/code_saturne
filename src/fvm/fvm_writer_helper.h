#ifndef __FVM_WRITER_HELPER_H__
#define __FVM_WRITER_HELPER_H__

/*============================================================================
 * Helper types and functions for mesh and field writers
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include "fvm_gather.h"
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

  _Bool   continues_previous;           /* Indicates if the corresponding FVM
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
 *   group_same_type      <-- group sections of the same type
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
fvm_writer_export_list(const fvm_nodal_t  *mesh,
                       int                 min_export_dim,
                       _Bool               group_same_type,
                       _Bool               discard_polygons,
                       _Bool               discard_polyhedra,
                       _Bool               divide_polygons,
                       _Bool               divide_polyhedra);

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
 * parameters:
 *   helper <-> pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_writer_field_helper_t *
fvm_writer_field_helper_destroy(fvm_writer_field_helper_t *helper);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize global values for an fvm_writer_field_helper structure.
 *
 * Internal buffers for gathering of data to rank 0 of the given
 * communicator are also allocated.
 *
 * parameters:
 *   helper        <-> pointer to structure that should be initialized
 *   section_list  <-- point to export section list helper structure
 *   mesh          <-- pointer to nodal mesh structure
 *   comm          <-- associated MPI communicator
 *
 * returns:
 *   pointer to allocated and initialized field writer helper
 *----------------------------------------------------------------------------*/

void
fvm_writer_field_helper_init_g(fvm_writer_field_helper_t   *helper,
                               const fvm_writer_section_t  *section_list,
                               const fvm_nodal_t           *mesh,
                               MPI_Comm                     comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Return sizes associated with a writer field helper.
 *
 * parameters:
 *   helper                   <-- pointer to helper structure
 *   input_size               --> Total field locations in input (or NULL)
 *   output_size              --> Total field locations in output (or NULL)
 *   max_grouped_elements_out --> Max. field locations in a single group
 *                                (elements of a given type if sections are
 *                                grouped, elements of a given section
 *                                otherwise; NULL if unused)
 *   min_output_buffer_size   --> Minimum required buffer size (or NULL)
 *----------------------------------------------------------------------------*/

void
fvm_writer_field_helper_get_size(const fvm_writer_field_helper_t  *helper,
                                 size_t  *input_size,
                                 size_t  *output_size,
                                 size_t  *max_grouped_elements_out,
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
 * Partially distribute field values to an output buffer.
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
fvm_writer_field_helper_step_e(fvm_writer_field_helper_t   *helper,
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
 * Partially distribute per node field values to an output buffer.
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
fvm_writer_field_helper_step_n(fvm_writer_field_helper_t   *helper,
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

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_WRITER_HELPER_H__ */
