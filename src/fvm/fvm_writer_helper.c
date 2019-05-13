/*============================================================================
 * Helper types and functions for mesh and field writers
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_convert_array.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_writer.h"
#include "fvm_writer_helper.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define FVM_WRITER_MIN_ELEMENTS     32
#define FVM_WRITER_MIN_SUB_ELEMENTS 32

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * FVM nodal to writer field output helper
 *----------------------------------------------------------------------------*/

struct _fvm_writer_field_helper_t {

  int                         field_dim;        /* Field dimension */
  cs_interlace_t              interlace;        /* Field interlaced or not */
  cs_datatype_t               datatype;         /* Output datatype */
  fvm_writer_var_loc_t        location;         /* Variable location */

  cs_gnum_t                   input_size;       /* Total input location size
                                                   (field values / dimension) */
  cs_gnum_t                   output_size;      /* Total output location size
                                                   (field values / dimension) */

  /* Global dimensions */

  cs_gnum_t   n_g_vertices_add;            /* Global number of added vertices */

  /* Local dimensions */

  cs_lnum_t   n_vertices_add;              /* Local number of added vertices */

  cs_lnum_t   n_sub_elements_max;          /* Maximum number of sub-elements
                                              per element (1 if no tesselation
                                              is used, global maximum in
                                              parallel mode) */

  /* State related fields */

  cs_lnum_t                   start_id;      /* Local section start */
  const fvm_writer_section_t *last_section;  /* Current section pointer */

  int         n_ranks;                       /* Number of ranks
                                                in communicator */

#if defined(HAVE_MPI)

  /* Additionnal parallel state */

  MPI_Comm     comm;                       /* Associated MPI communicator */
  int          rank;                       /* Local Rank in communicator */
  int          min_rank_step;              /* Minimum rank step for output */
  int          min_block_size;             /* Minimum block size for output */

#endif
};

/*============================================================================
 * Static and constant variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Reorder data components for output if needed.
 *
 * parameters:
 *   n_ent      <-- number of entities
 *   n_comp     <-- number of components
 *   datatype   <-- associated datatype
 *   comp_order <-- component reordering array
 *   val        <-> values buffer
 *----------------------------------------------------------------------------*/

static void
_reorder_components(size_t          n_ent,
                    size_t          n_comp,
                    cs_datatype_t   datatype,
                    const int       comp_order[],
                    unsigned char   val[])
{
  if (comp_order != NULL && n_comp > 1) {
    size_t datatype_size = cs_datatype_size[datatype];
    assert(datatype_size <= 8);
    assert(n_comp <= 9);
    size_t comp_size = datatype_size*n_comp;
    for (size_t j = 0; j < n_ent; j++) {
      unsigned char swap_buf[72];
      for (size_t k = 0; k < comp_size; k++)
        swap_buf[k] = val[j*comp_size + k];
      for (size_t k = 0; k < n_comp; k++) {
        const size_t k0 = k*datatype_size;
        const size_t k1 = comp_order[k]*datatype_size;
        for (size_t l = 0; l < datatype_size; l++)
          val[j*comp_size + k0 + l] = swap_buf[k1 + l];
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Zero array values for output if needed.
 *
 * parameters:
 *   n_ent      <-- number of entities
 *   val        <-> values buffer
 *----------------------------------------------------------------------------*/

static void
_zero_values(size_t          n_ent,
             cs_datatype_t   datatype,
             void           *val)
{
  switch(datatype) {
  case CS_FLOAT:
    {
      float *val_p = val;
      for (size_t i = 0; i < n_ent; i++)
        val_p[i] = 0.0;
    }
    break;
  case CS_DOUBLE:
    {
      double *val_p = val;
      for (size_t i = 0; i < n_ent; i++)
        val_p[i] = 0.0;
    }
    break;
  case CS_INT32:
    {
      int32_t *val_p = val;
      for (size_t i = 0; i < n_ent; i++)
        val_p[i] = 0.0;
    }
    break;
  case CS_INT64:
    {
      int64_t *val_p = val;
      for (size_t i = 0; i < n_ent; i++)
        val_p[i] = 0.0;
    }
    break;
  case CS_UINT32:
    {
      uint32_t *val_p = val;
      for (size_t i = 0; i < n_ent; i++)
        val_p[i] = 0.0;
    }
    break;
  case CS_UINT64:
    {
      uint64_t *val_p = val;
      for (size_t i = 0; i < n_ent; i++)
        val_p[i] = 0.0;
    }
    break;
  default:
    assert(0);
  }
}

/*----------------------------------------------------------------------------
 * Compare sections by associated type
 *
 * parameters:
 *   s1 <-- pointer to first section
 *   s2 <-- pointer to second section
 *
 * returns:
 *   difference of section types
 *----------------------------------------------------------------------------*/

static int
_compare_sections(const void *s1_p,
                  const void *s2_p)
{
  const fvm_writer_section_t *s1 = s1_p;
  const fvm_writer_section_t *s2 = s2_p;

  return ((int)(s1->type) - (int)(s2->type));
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get extra vertex global numbers when tesselations are present
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   n_extra_vertices <-- number of extra vertices
 *   vtx_gnum         --> extra vertex global numbers (size: n_extra_vertices)
 *----------------------------------------------------------------------------*/

static void
_extra_vertex_get_gnum(const fvm_nodal_t  *mesh,
                       cs_lnum_t           n_extra_vertices,
                       cs_gnum_t           vtx_gnum[])
{
  int  i;
  cs_lnum_t   j = 0;
  cs_lnum_t   start_id = 0;
  cs_gnum_t   gnum_shift
    = fvm_io_num_get_global_count(mesh->global_vertex_num);

  if (n_extra_vertices > 0) { /* Implies divide_polyhedra */

    for (i = 0; i < mesh->n_sections; i++) {

      const fvm_nodal_section_t  *const  section = mesh->sections[i];

      if (   section->type == FVM_CELL_POLY
          && section->tesselation != NULL) {

        cs_lnum_t n_extra_vertices_section
          = fvm_tesselation_n_vertices_add(section->tesselation);

        if (n_extra_vertices_section > 0) {

          const fvm_io_num_t *extra_vertex_num
            = fvm_tesselation_global_vertex_num(section->tesselation);
          const cs_gnum_t *extra_gnum
            = fvm_io_num_get_global_num(extra_vertex_num);

          for (j = 0; j < n_extra_vertices_section; j++)
            vtx_gnum[start_id + j] = extra_gnum[j] + gnum_shift;

          start_id += n_extra_vertices_section;

        }

        gnum_shift
          = fvm_tesselation_n_g_vertices_add(section->tesselation);

      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Output per-element field values in parallel mode.
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

#if defined(__INTEL_COMPILER) && defined(__KNC__)
#pragma optimization_level 1 /* Crash with O2 on KNC with icc 14.0.0 20130728 */
#endif

static const fvm_writer_section_t *
_field_helper_output_eg(fvm_writer_field_helper_t          *helper,
                        void                               *context,
                        const fvm_writer_section_t         *export_section,
                        int                                 src_dim,
                        cs_interlace_t                      src_interlace,
                        const int                          *comp_order,
                        int                                 n_parent_lists,
                        const cs_lnum_t                     parent_num_shift[],
                        cs_datatype_t                       datatype,
                        const void                   *const field_values[],
                        fvm_writer_field_output_t          *output_func)
{
  fvm_writer_field_helper_t *h = helper;

  cs_block_dist_info_t  bi;
  cs_part_to_block_t  *d = NULL;

  int         n_sections = 0;
  bool        have_tesselation = false;
  cs_lnum_t   part_size = 0, block_size = 0;
  cs_gnum_t   block_sub_size = 0, block_start = 0, block_end = 0;
  cs_gnum_t   n_g_elements = 0;

  int  *part_n_sub = NULL, *block_n_sub = NULL;
  unsigned char  *part_values = NULL;
  unsigned char  *block_values = NULL, *_block_values = NULL;

  cs_gnum_t         *_g_elt_num = NULL;
  const cs_gnum_t   *g_elt_num
    = fvm_io_num_get_global_num(export_section->section->global_element_num);

  const fvm_writer_section_t  *current_section = NULL;

  const size_t stride = (h->interlace == CS_INTERLACE) ? h->field_dim : 1;
  const size_t elt_size = cs_datatype_size[h->datatype];
  const size_t min_block_size =   h->min_block_size
                                / (elt_size*stride);

  /* Loop on sections to count output size */

  current_section = export_section;
  do {

    const fvm_nodal_section_t  *section = current_section->section;

    n_sections += 1;
    n_g_elements += fvm_io_num_get_global_count(section->global_element_num);
    part_size += fvm_io_num_get_local_count(section->global_element_num);
    if (current_section->type != section->type)
      have_tesselation = true;

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Build global numbering if necessary */

  if (n_sections > 1) {

    cs_lnum_t start_id = 0;
    cs_gnum_t gnum_shift = 0;

    BFT_MALLOC(_g_elt_num, part_size, cs_gnum_t);
    g_elt_num = _g_elt_num;

    /* loop on sections which should be appended */

    current_section = export_section;
    do {

      const fvm_nodal_section_t  *section = current_section->section;
      const cs_lnum_t section_size
        = fvm_io_num_get_local_count(section->global_element_num);

      const cs_gnum_t * s_gnum
        = fvm_io_num_get_global_num(section->global_element_num);

      for (cs_lnum_t j = 0, k = start_id; j < section_size; j++, k++)
        _g_elt_num[k] = s_gnum[j] + gnum_shift;

      start_id += section_size;
      gnum_shift += fvm_io_num_get_global_count(section->global_element_num);

      current_section = current_section->next;


    } while (   current_section != NULL
             && current_section->continues_previous == true);
  }

  /* Build sub-element count if necessary */

  if (have_tesselation) {

    cs_lnum_t start_id = 0;

    BFT_MALLOC(part_n_sub, part_size, int);

    current_section = export_section;
    do {

      const fvm_nodal_section_t  *section = current_section->section;
      const cs_lnum_t section_size
        = fvm_io_num_get_local_count(section->global_element_num);

      if (current_section->type != section->type) {
        const cs_lnum_t   *sub_element_idx
          = fvm_tesselation_sub_elt_index(section->tesselation,
                                          current_section->type);
        for (cs_lnum_t j = 0; j < section_size; j++)
          part_n_sub[start_id + j] = sub_element_idx[j+1] - sub_element_idx[j];
      }
      else {
        for (cs_lnum_t j = 0; j < section_size; j++)
          part_n_sub[start_id + j] = 1;
      }
      start_id += section_size;

      current_section = current_section->next;

    } while (   current_section != NULL
             && current_section->continues_previous == true);
  }

  /* Build distribution structures */

  bi = cs_block_dist_compute_sizes(h->rank,
                                   h->n_ranks,
                                   h->min_rank_step,
                                   min_block_size,
                                   n_g_elements);

  block_size = bi.gnum_range[1] - bi.gnum_range[0];

  d = cs_part_to_block_create_by_gnum(h->comm, bi, part_size, g_elt_num);

  if (_g_elt_num != NULL)
    cs_part_to_block_transfer_gnum(d, _g_elt_num);

  g_elt_num = NULL;
  _g_elt_num = NULL;

  /* Distribute sub-element info in case of tesselation */

  if (have_tesselation) {

    BFT_MALLOC(block_n_sub, block_size, int);

    cs_part_to_block_copy_array(d,
                                CS_INT_TYPE,
                                1,
                                part_n_sub,
                                block_n_sub);
    BFT_FREE(part_n_sub);

    for (cs_lnum_t j = 0; j < block_size; j++)
      block_sub_size += block_n_sub[j];

  }
  else
    block_sub_size = block_size;

  /* Number of loops on dimension and conversion output dimension */

  const int n_dim_loops = (h->interlace == CS_INTERLACE) ? 1 : h->field_dim;
  const int convert_dim = (h->interlace == CS_INTERLACE) ? h->field_dim : 1;

  /* To save space, in case of tesselation, part_values and _block_values
     point to the same memory space, as they are not needed simultaneously.
     Without tesselation, _block_values simply points to block_values */

  BFT_MALLOC(block_values, block_size*elt_size*stride, unsigned char);

  if (have_tesselation) {
    BFT_MALLOC(part_values,
               (  CS_MAX(part_size, (cs_lnum_t)block_sub_size)
                * elt_size*convert_dim),
               unsigned char);
    MPI_Scan(&block_sub_size, &block_end, 1, CS_MPI_GNUM, MPI_SUM, h->comm);
    block_end += 1;
    block_start = block_end - block_sub_size;
    _block_values = part_values;
  }
  else {
    BFT_MALLOC(part_values, part_size*elt_size*convert_dim, unsigned char);
    block_start = bi.gnum_range[0];
    block_end = bi.gnum_range[1];
    _block_values = block_values;
  }

  /* Loop on dimension (for non-interlaced output) */

  for (int comp_id = 0; comp_id < n_dim_loops; comp_id++) {

    /* Distribute partition to block values */

    if (comp_id < src_dim) {

      cs_lnum_t start_id = 0;
      cs_lnum_t src_shift = 0;

      const int comp_id_in = comp_order != NULL ? comp_order[comp_id] : comp_id;

      /* loop on sections which should be appended */

      current_section = export_section;
      do {

        unsigned char *_part_values =   part_values
                                      + (size_t)start_id*elt_size*convert_dim;

        const fvm_nodal_section_t  *section = current_section->section;

        if (n_parent_lists == 0)
          src_shift = current_section->num_shift;

        fvm_convert_array(src_dim,
                          comp_id_in,
                          convert_dim,
                          src_shift,
                          section->n_elements + src_shift,
                          src_interlace,
                          datatype,
                          h->datatype,
                          n_parent_lists,
                          parent_num_shift,
                          section->parent_element_num,
                          field_values,
                          _part_values);

        start_id += fvm_io_num_get_local_count(section->global_element_num);

        current_section = current_section->next;

      } while (   current_section != NULL
               && current_section->continues_previous == true);

      /* Reorder components if required
         (for interlaced output; done though dim_loops if non-interlaced) */

      if (comp_order != NULL && convert_dim > 1)
        _reorder_components(part_size,
                            convert_dim,
                            h->datatype,
                            comp_order,
                            part_values);

      /* Distribute part values */

      cs_part_to_block_copy_array(d,
                                  h->datatype,
                                  convert_dim,
                                  part_values,
                                  block_values);

      /* Scatter values to sub-elements in case of tesselation */

      if (have_tesselation) {
        size_t i = 0;
        size_t comp_size = elt_size * convert_dim;
        for (size_t j = 0; j < (size_t)block_size; j++) {
          for (int k = 0; k < block_n_sub[j]; k++) {
            for (size_t l = 0; l < comp_size; l++)
              _block_values[i++] = block_values[j*comp_size + l];
          }
        }
      }

    }

    /* Zero extra dimensions
       (for non-interlaced output; done in fvm_convert_array if interlaced) */

    else
      _zero_values(block_sub_size, h->datatype, _block_values);

    /* Write block values */

    output_func(context,
                h->datatype,
                h->field_dim,
                comp_id,
                block_start,
                block_end,
                _block_values);

  } /* end of loop on spatial dimension */

  BFT_FREE(block_values);
  BFT_FREE(part_values);

  cs_part_to_block_destroy(&d);

  if (block_n_sub != NULL)
    BFT_FREE(block_n_sub);

  /* Return pointer to next section */

  return current_section;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Output per-element field values in serial mode.
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

static const fvm_writer_section_t *
_field_helper_output_el(fvm_writer_field_helper_t          *helper,
                        void                               *context,
                        const fvm_writer_section_t         *export_section,
                        int                                 src_dim,
                        cs_interlace_t                      src_interlace,
                        const int                          *comp_order,
                        int                                 n_parent_lists,
                        const cs_lnum_t                     parent_num_shift[],
                        cs_datatype_t                       datatype,
                        const void                   *const field_values[],
                        fvm_writer_field_output_t          *output_func)
{
  fvm_writer_field_helper_t *h = helper;

  int         n_sections = 0;
  cs_lnum_t   sub_size = 0;
  cs_lnum_t   n_elements = 0;

  unsigned char  *values = NULL;

  const fvm_writer_section_t  *current_section = NULL;

  const size_t elt_size = cs_datatype_size[h->datatype];

  /* Loop on sections to count output size */

  current_section = export_section;
  do {

    const fvm_nodal_section_t  *section = current_section->section;

    n_sections += 1;
    n_elements += section->n_elements;
    if (current_section->type != section->type)
      sub_size += fvm_tesselation_n_sub_elements(section->tesselation,
                                                 current_section->type);
    else
      sub_size += section->n_elements;

    current_section = current_section->next;

  } while (   current_section != NULL
           && current_section->continues_previous == true);

  /* Number of loops on dimension and conversion output dimension */

  const int n_dim_loops = (h->interlace == CS_INTERLACE) ? 1 : h->field_dim;
  const int convert_dim = (h->interlace == CS_INTERLACE) ? h->field_dim : 1;

  BFT_MALLOC(values,
             (  CS_MAX(n_elements, (cs_lnum_t)sub_size)
              * elt_size*convert_dim),
             unsigned char);

  /* Loop on dimension (for non-interlaced output) */

  for (int comp_id = 0; comp_id < n_dim_loops; comp_id++) {

    /* Extract values */

    if (comp_id < src_dim) {

      cs_lnum_t start_id = 0;
      cs_lnum_t src_shift = 0;

      const int comp_id_in = comp_order != NULL ? comp_order[comp_id] : comp_id;

      /* loop on sections which should be appended */

      current_section = export_section;
      do {

        unsigned char *_values = values + (size_t)start_id*elt_size*convert_dim;

        const fvm_nodal_section_t  *section = current_section->section;

        if (n_parent_lists == 0)
          src_shift = current_section->num_shift;

        fvm_convert_array(src_dim,
                          comp_id_in,
                          convert_dim,
                          src_shift,
                          section->n_elements + src_shift,
                          src_interlace,
                          datatype,
                          h->datatype,
                          n_parent_lists,
                          parent_num_shift,
                          section->parent_element_num,
                          field_values,
                          _values);

        /* Scatter values to sub-elements in case of tesselation */

        if (current_section->type != section->type) {
          fvm_tesselation_distribute(section->tesselation,
                                     export_section->type,
                                     0,
                                     section->n_elements,
                                     elt_size*convert_dim,
                                     _values);
          start_id += fvm_tesselation_n_sub_elements(section->tesselation,
                                                     current_section->type);
        }
        else
          start_id += section->n_elements;

        current_section = current_section->next;

      } while (   current_section != NULL
               && current_section->continues_previous == true);

      /* Reorder components if required
         (for interlaced output; done though dim_loops if non-interlaced) */

      if (comp_order != NULL && convert_dim > 1)
        _reorder_components(sub_size,
                            convert_dim,
                            h->datatype,
                            comp_order,
                            values);

    }

    /* Zero extra dimensions
       (for non-interlaced output; done in fvm_convert_array if interlaced) */

    else
      _zero_values(sub_size, h->datatype, values);

    /* Write block values */

    {
      int eo = 0;
      for (cs_lnum_t ii = 0; ii < convert_dim*sub_size; ii++) {
        if (values[ii] != 0)
          eo++;
      }
    }
    output_func(context,
                h->datatype,
                h->field_dim,
                comp_id,
                1,             /* block_start, */
                sub_size + 1,  /* block_end */
                values);

  } /* end of loop on spatial dimension */

  BFT_FREE(values);

  /* Return pointer to next section */

  return current_section;
}

#if defined(HAVE_MPI )

/*----------------------------------------------------------------------------
 * Output per-node field values in parallel mode.
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

static void
_field_helper_output_ng(fvm_writer_field_helper_t        *helper,
                        void                             *context,
                        const fvm_nodal_t                *mesh,
                        int                               src_dim,
                        cs_interlace_t                    src_interlace,
                        const int                        *comp_order,
                        int                               n_parent_lists,
                        const cs_lnum_t                   parent_num_shift[],
                        cs_datatype_t                     datatype,
                        const void                 *const field_values[],
                        fvm_writer_field_output_t        *output_func)
{
  fvm_writer_field_helper_t *h = helper;

  cs_block_dist_info_t  bi;

  cs_lnum_t       part_size = 0, block_size = 0;
  unsigned char  *part_values = NULL, *block_values = NULL;
  cs_part_to_block_t  *d = NULL;

  const size_t stride = (h->interlace == CS_INTERLACE) ? h->field_dim : 1;
  const size_t elt_size = cs_datatype_size[h->datatype];
  const size_t min_block_size =   h->min_block_size
                                / (elt_size*stride);

  /* Initialize distribution info */

  fvm_writer_vertex_part_to_block_create(h->min_rank_step,
                                         min_block_size,
                                         helper->n_g_vertices_add,
                                         helper->n_vertices_add,
                                         mesh,
                                         &bi,
                                         &d,
                                         h->comm);

  part_size = cs_part_to_block_get_n_part_ents(d);
  block_size = bi.gnum_range[1] - bi.gnum_range[0];

  /* Number of loops on dimension and conversion output dimension */

  const int n_dim_loops = (h->interlace == CS_INTERLACE) ? 1 : h->field_dim;
  const int convert_dim = (h->interlace == CS_INTERLACE) ? h->field_dim : 1;

  BFT_MALLOC(part_values, part_size*elt_size*convert_dim, unsigned char);
  BFT_MALLOC(block_values, block_size*elt_size*convert_dim, unsigned char);

  /* Loop on dimension (for non-interlaced output) */

  for (int comp_id = 0; comp_id < n_dim_loops; comp_id++) {

    /* Distribute partition to block values */

    if (comp_id < src_dim) {

      cs_lnum_t start_id = 0;
      cs_lnum_t end_id = mesh->n_vertices;

      const int comp_id_in = comp_order != NULL ? comp_order[comp_id] : comp_id;

      /* Distribute partition to block values */

      /* Main vertices */

      fvm_convert_array(src_dim,
                        comp_id_in,
                        convert_dim,
                        start_id,
                        end_id,
                        src_interlace,
                        datatype,
                        h->datatype,
                        n_parent_lists,
                        parent_num_shift,
                        mesh->parent_vertex_num,
                        field_values,
                        part_values);

      /* Additional vertices in case of tesselation
         (end_id == part_size with no tesselation or if all tesselated
         sections have been accounted for).*/

      if (helper->n_vertices_add > 0) {

        for (cs_lnum_t j = 0; end_id < part_size && j < mesh->n_sections; j++) {

          const fvm_nodal_section_t  *section = mesh->sections[j];

          if (section->type == FVM_CELL_POLY && section->tesselation != NULL) {

            cs_lnum_t   n_extra_vertices
              = fvm_tesselation_n_vertices_add(section->tesselation);

            start_id = end_id;
            end_id = start_id + n_extra_vertices;

            size_t part_values_shift =   (size_t)start_id
                                       * (size_t)convert_dim
                                       * elt_size;

            fvm_tesselation_vertex_values(section->tesselation,
                                          src_dim,
                                          comp_id_in,
                                          convert_dim,
                                          0,
                                          n_extra_vertices,
                                          h->interlace,
                                          datatype,
                                          h->datatype,
                                          n_parent_lists,
                                          parent_num_shift,
                                          mesh->parent_vertex_num,
                                          field_values,
                                          part_values + part_values_shift);

          }

        } /* End of loops on tesselated sections */

        assert(end_id == part_size);

      }

      /* Reorder components if required
         (for interlaced output; done though dim_loops if non-interlaced) */

      if (comp_order != NULL && convert_dim > 1)
        _reorder_components(part_size,
                            convert_dim,
                            h->datatype,
                            comp_order,
                            part_values);

      /* Distribute part values */

      cs_part_to_block_copy_array(d,
                                  h->datatype,
                                  convert_dim,
                                  part_values,
                                  block_values);

    }

    /* Zero extra dimensions
       (for non-interlaced output; done in fvm_convert_array if interlaced) */

    else
      _zero_values(part_size, h->datatype, block_values);

    output_func(context,
                h->datatype,
                h->field_dim,
                comp_id,
                bi.gnum_range[0],
                bi.gnum_range[1],
                block_values);

  }

  BFT_FREE(block_values);
  BFT_FREE(part_values);

  cs_part_to_block_destroy(&d);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Output per-node field values in serial mode.
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

static void
_field_helper_output_nl(fvm_writer_field_helper_t        *helper,
                        void                             *context,
                        const fvm_nodal_t                *mesh,
                        int                               src_dim,
                        cs_interlace_t                    src_interlace,
                        const int                        *comp_order,
                        int                               n_parent_lists,
                        const cs_lnum_t                   parent_num_shift[],
                        cs_datatype_t                     datatype,
                        const void                 *const field_values[],
                        fvm_writer_field_output_t        *output_func)
{
  fvm_writer_field_helper_t *h = helper;

  cs_lnum_t       n_vertices = mesh->n_vertices + helper->n_vertices_add;
  unsigned char  *values = NULL;

  const size_t elt_size = cs_datatype_size[h->datatype];

  /* Number of loops on dimension and conversion output dimension */

  const int n_dim_loops = (h->interlace == CS_INTERLACE) ? 1 : h->field_dim;
  const int convert_dim = (h->interlace == CS_INTERLACE) ? h->field_dim : 1;

  /* Allocate buffer */

  BFT_MALLOC(values, n_vertices*elt_size*convert_dim, unsigned char);

  /* Loop on dimension (for non-interlaced output) */

  for (int comp_id = 0; comp_id < n_dim_loops; comp_id++) {

    /* Distribute partition to block values */

    if (comp_id < src_dim) {

      cs_lnum_t start_id = 0;
      cs_lnum_t end_id = mesh->n_vertices;

      const int comp_id_in = comp_order != NULL ? comp_order[comp_id] : comp_id;

      /* Distribute partition to block values */

      /* Main vertices */

      fvm_convert_array(src_dim,
                        comp_id_in,
                        convert_dim,
                        start_id,
                        end_id,
                        src_interlace,
                        datatype,
                        h->datatype,
                        n_parent_lists,
                        parent_num_shift,
                        mesh->parent_vertex_num,
                        field_values,
                        values);

      /* Additional vertices in case of tesselation
         (end_id == part_size with no tesselation or if all tesselated
         sections have been accounted for).*/

      if (helper->n_vertices_add > 0) {

        for (int i = 0; i < mesh->n_sections; i++) {

          const fvm_nodal_section_t  *section = mesh->sections[i];

          if (section->type == FVM_CELL_POLY && section->tesselation != NULL) {

            cs_lnum_t n_extra_vertices
              = fvm_tesselation_n_vertices_add(section->tesselation);

            start_id = end_id;
            end_id = start_id + n_extra_vertices;

            size_t values_shift =   (size_t)start_id
                                  * (size_t)convert_dim
                                  * elt_size;

            fvm_tesselation_vertex_values(section->tesselation,
                                          src_dim,
                                          comp_id_in,
                                          convert_dim,
                                          0,
                                          n_extra_vertices,
                                          h->interlace,
                                          datatype,
                                          h->datatype,
                                          n_parent_lists,
                                          parent_num_shift,
                                          mesh->parent_vertex_num,
                                          field_values,
                                          values + values_shift);

          }

        } /* End of loops on tesselated sections */

        assert(end_id == n_vertices);

      }

      /* Reorder components if required
         (for interlaced output; done though dim_loops if non-interlaced) */

      if (comp_order != NULL && convert_dim > 1)
        _reorder_components(n_vertices,
                            convert_dim,
                            h->datatype,
                            comp_order,
                            values);

    }

    /* Zero extra dimensions
       (for non-interlaced output; done in fvm_convert_array if interlaced) */

    else
      _zero_values(n_vertices, h->datatype, values);

    output_func(context,
                h->datatype,
                h->field_dim,
                comp_id,
                1,
                n_vertices + 1,
                values);

  }

  BFT_FREE(values);
}

/*============================================================================
 * Semi-private function definitions (prototypes in fvm_writer_helper.h)
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
                       bool                        divide_polyhedra)
{
  int  i, j;
  int  n_sections = 0;
  int  n_sub_types = 0;
  fvm_element_t sub_type[FVM_TESSELATION_N_SUB_TYPES_MAX];
  fvm_writer_section_t *export_list = NULL;

  cs_lnum_t   num_shift = 0;
  cs_gnum_t   extra_vertex_base = fvm_nodal_n_g_vertices(mesh) + 1;

  /* Initial count and allocation */

  n_sections = 0;

  for (i = 0 ; i < mesh->n_sections ; i++) {

    const fvm_nodal_section_t  *const  section = mesh->sections[i];

    /* Output if entity dimension equal to highest in mesh
       (i.e. no output of faces if cells present, or edges
       if cells or faces) */

    if (section->entity_dim >= min_export_dim) {

      /* Optionally discard polygons or polyhedra */
      if (   (section->type == FVM_FACE_POLY && discard_polygons == true)
          || (section->type == FVM_CELL_POLY && discard_polyhedra == true))
        continue;

      /* Possibly add multiple sections when dividing polygons or polyhedra
         (ignore those sections if no tesselation present) */
      else if (   (section->type == FVM_FACE_POLY && divide_polygons == true)
               || (section->type == FVM_CELL_POLY && divide_polyhedra == true)) {
        if (section->tesselation != NULL)
          n_sections += fvm_tesselation_n_sub_types(section->tesselation);
      }

      else
        n_sections += 1;

    }

  }

  /* If no sections are present no list is returned */

  if (n_sections == 0)
    return NULL;

  BFT_MALLOC(export_list, n_sections, fvm_writer_section_t);

  for (i = 0 ; i < n_sections - 1 ; i++)
    (export_list[i]).next = export_list + i + 1;
  (export_list[n_sections - 1]).next = NULL;

  /* Build unsorted list */

  n_sections = 0;

  for (i = 0 ; i < mesh->n_sections ; i++) {

    const fvm_tesselation_t    *tesselation = NULL;
    const fvm_nodal_section_t  *const  section = mesh->sections[i];

    /* Ignore sections with entity dimension other than the highest */
    if (section->entity_dim < min_export_dim)
      continue;

    /* Ignore polygonal or polyhedra sections if they should be discarded */
    if (   (section->type == FVM_FACE_POLY && discard_polygons == true)
        || (section->type == FVM_CELL_POLY && discard_polyhedra == true))
      continue;

    /* Depending on section type and tesselation,
       we may have multiple sub-types */

    n_sub_types = 1;
    sub_type[0] = section->type;

    if (   (section->type == FVM_FACE_POLY && divide_polygons == true)
        || (section->type == FVM_CELL_POLY && divide_polyhedra == true)) {

      if (section->tesselation != NULL) {
        tesselation = section->tesselation;
        n_sub_types = fvm_tesselation_n_sub_types(section->tesselation);
        for (j = 0; j < n_sub_types; j++)
          sub_type[j] = fvm_tesselation_sub_type(section->tesselation, j);
      }
      else
        continue; /* ignore section if no tesselation present */

    }

    for (j = 0; j < n_sub_types; j++) {

      (export_list[n_sections]).section = section;
      (export_list[n_sections]).type = sub_type[j];
      (export_list[n_sections]).continues_previous = false;
      if (tesselation == NULL)
        (export_list[n_sections]).extra_vertex_base = 0;
      else
        (export_list[n_sections]).extra_vertex_base = extra_vertex_base;
      (export_list[n_sections]).num_shift = num_shift;

      n_sections++;

    }

    if (tesselation != NULL)
      extra_vertex_base += fvm_tesselation_n_g_vertices_add(tesselation);

    num_shift += section->n_elements;

  }

  /* Now order list so as to merge sections of similar type */

  if (group_by_type == true && n_sections > 1) {

    qsort(export_list, n_sections, sizeof(fvm_writer_section_t),
          _compare_sections);

    for (i = 1; i < n_sections; i++) {
      if ((export_list[i-1]).type == (export_list[i]).type)
        (export_list[i]).continues_previous = true;
    }

  }

  /* Also group all sections if required */

  if (group_all) {
    for (i = 1; i < n_sections; i++)
      (export_list[i]).continues_previous = true;
  }

  for (i = 0; i < n_sections - 1; i++)
    (export_list[i]).next = &(export_list[i+1]);
  export_list[n_sections - 1].next = NULL;

  return export_list;
}

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
                                cs_lnum_t          *n_extra_vertices)
{
  int  i;

  const int  export_dim = fvm_nodal_get_max_entity_dim(mesh);

  /* Initial count and allocation */

  if (n_extra_vertices_g != NULL)
    *n_extra_vertices_g = 0;
  if (n_extra_vertices != NULL)
    *n_extra_vertices   = 0;

  if (divide_polyhedra) {

    for (i = 0; i < mesh->n_sections; i++) {

      const fvm_nodal_section_t  *section = mesh->sections[i];

      /* Output if entity dimension equal to highest in mesh
         (i.e. no output of faces if cells present, or edges
         if cells or faces) */

      if (   section->entity_dim == export_dim
          && section->type == FVM_CELL_POLY
          && section->tesselation != NULL) {

        if (n_extra_vertices_g != NULL)
          *n_extra_vertices_g
            += fvm_tesselation_n_g_vertices_add(section->tesselation);

        if (n_extra_vertices != NULL)
          *n_extra_vertices
            += fvm_tesselation_n_vertices_add(section->tesselation);

      }
    }
  }
}

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
                                       MPI_Comm                comm)
{
  cs_gnum_t    n_g_vertices_tot = 0;
  cs_lnum_t    n_vertices_tot = 0;

  cs_gnum_t   *_g_num = NULL;

  int                    rank, n_ranks;
  cs_block_dist_info_t  _bi;
  cs_part_to_block_t   *_d;

  const cs_lnum_t   n_vertices
    = fvm_io_num_get_local_count(mesh->global_vertex_num);
  cs_gnum_t   n_g_vertices
    = fvm_io_num_get_global_count(mesh->global_vertex_num);
  const cs_gnum_t   *g_num
    = fvm_io_num_get_global_num(mesh->global_vertex_num);

  /* Communicator info */

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &n_ranks);

  /* Vertex counts */

  n_vertices_tot = n_vertices + n_add_vertices;
  n_g_vertices_tot = n_g_vertices + n_g_add_vertices;

  _bi = cs_block_dist_compute_sizes(rank,
                                    n_ranks,
                                    min_rank_step,
                                    min_block_size,
                                    n_g_vertices_tot);

  /* Global vertex numbers */

  if (n_add_vertices > 0) {

    BFT_MALLOC(_g_num, n_vertices_tot, cs_gnum_t);

    memcpy(_g_num, g_num, n_vertices*sizeof(cs_gnum_t));
    _extra_vertex_get_gnum(mesh, n_add_vertices, _g_num + n_vertices);

    g_num = _g_num;

  }

  /* Build distribution structures */

  _d = cs_part_to_block_create_by_gnum(comm, _bi, n_vertices_tot, g_num);

  if (n_add_vertices > 0)
    cs_part_to_block_transfer_gnum(_d, _g_num);

  /* Return initialized structures */

  if (bi != NULL)
    *bi = _bi;

  if (d != NULL)
    *d = _d;
}

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
                               cs_lnum_t           n_extra_vertices)
{
  int  i;
  cs_lnum_t   n_extra_vertices_section;

  size_t  coord_shift = 0;
  cs_coord_t  *coords = NULL;

  if (n_extra_vertices > 0) { /* This implies divide_polyhedra = true */

    BFT_MALLOC(coords, n_extra_vertices * 3, cs_coord_t);

    for (i = 0; i < mesh->n_sections; i++) {

      const fvm_nodal_section_t  *const  section = mesh->sections[i];

      if (   section->type == FVM_CELL_POLY
          && section->tesselation != NULL) {

        n_extra_vertices_section
          = fvm_tesselation_n_vertices_add(section->tesselation);

        if (n_extra_vertices_section > 0) {

          fvm_tesselation_vertex_coords(section->tesselation,
                                        coords + coord_shift);

          coord_shift += n_extra_vertices_section * 3;

        }

      }
    }
  }

  return coords;
}

/*----------------------------------------------------------------------------
 * Create field writer helper structure.
 *
 * Local values are initialized, ang lobal values are set to zero
 * (they may be initialized by calling fvm_writer_field_helper_init_g()).
 *
 * The mesh argument is not used when location is FVM_WRITER_PER_ELEMENT,
 * so NULL may be given instead in this case. The section_list
 * argument is used even when location is FVM_WRITER_PER_NODE to determine
 * if polyhedra are are divided (using extra vertices). Thus, NULL
 * may be given instead only if we are sure the corresponding exported mesh
 * does not contain divided polyhedra.
 *
 * parameters:
 *   mesh            <-- pointer to nodal mesh structure
 *   section_list    <-- point to export section list helper structure
 *   field_dim       <-- indicates output field dimension
 *   interlace       <-- indicates if output is interlaced
 *   location        <-- indicates if field is cell or node based
 *   datatype        <-- datatype of destination buffers
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
                               fvm_writer_var_loc_t        location)
{
  fvm_writer_field_helper_t *h = NULL;

  /* Initialize structure */

  BFT_MALLOC(h, 1, fvm_writer_field_helper_t);

  /* Initialize structure */

  h->field_dim = field_dim;
  h->interlace = interlace;
  h->datatype = datatype;
  h->location = location;

  /* Zero fields that will be computed later */

  h->input_size = 0;
  h->output_size = 0;

  h->n_g_vertices_add = 0;
  h->n_vertices_add = 0;

  h->n_sub_elements_max = 1;

  /* State values */

  h->start_id = 0;
  h->last_section = NULL;

  h->n_ranks = 1;

#if defined(HAVE_MPI)

  h->comm = MPI_COMM_NULL;
  h->rank = -1;
  h->min_rank_step = 1;
  h->min_block_size = 0;

#endif /* defined(HAVE_MPI) */

  /* Compute local dimensions */

  if (location == FVM_WRITER_PER_ELEMENT) {

    const fvm_writer_section_t *export_section = section_list;

    while (export_section != NULL) {

      const fvm_nodal_section_t *section = export_section->section;

      cs_lnum_t n_sub_elements;
      cs_lnum_t n_sub_elements_max = 1;

      if (export_section->type == section->type) { /* Regular section */
        n_sub_elements = section->n_elements;
      }
      else { /* Tesselated section */
        fvm_tesselation_get_global_size(section->tesselation,
                                        export_section->type,
                                        NULL,
                                        &n_sub_elements_max);
        n_sub_elements = fvm_tesselation_n_sub_elements(section->tesselation,
                                                        export_section->type);
      }

      /* Update global helper values */

      h->input_size  += section->n_elements;
      h->output_size += n_sub_elements;

      h->n_sub_elements_max = CS_MAX(h->n_sub_elements_max,
                                     n_sub_elements_max);

      /* continue with next section */

      export_section = export_section->next;

    } /* End of loop on sections */

  }

  else if (location == FVM_WRITER_PER_NODE) {

    int n_added_vertex_sections = 0;

    const fvm_writer_section_t *export_section;

    h->input_size  = mesh->n_vertices;
    h->output_size = mesh->n_vertices;

    /* Determine if polyhedra are tesselated */

    for (export_section = section_list;
         export_section != NULL;
         export_section = export_section->next) {

      const fvm_nodal_section_t *section = export_section->section;

      if (   export_section->type != section->type
          && section->type == FVM_CELL_POLY)
        n_added_vertex_sections += 1;

    }

    /* In case of polyhedra tesselation, count for added vertices;
       at this stage, we know that polyhedra are tesselated */

    if (n_added_vertex_sections > 0) {

      int ii, jj;

      for (ii = 0, jj = 0; ii < mesh->n_sections; ii++) {

        const fvm_nodal_section_t *section = mesh->sections[ii];

        if (section->type == FVM_CELL_POLY) {

          cs_lnum_t n_vertices_add
            = fvm_tesselation_n_vertices_add(section->tesselation);

          h->output_size += n_vertices_add;
          h->n_g_vertices_add
            = fvm_tesselation_n_g_vertices_add(section->tesselation);
          h->n_vertices_add += n_vertices_add;

          jj++;

        }

      }
    }

  }

  return h;
}

/*----------------------------------------------------------------------------
 * Destroy FVM field writer helper.
 *
 * parameters:
 *   helper <-> pointer to pointer to structure that should be destroyed
 *----------------------------------------------------------------------------*/

void
fvm_writer_field_helper_destroy(fvm_writer_field_helper_t **helper)
{
  if (helper != NULL)
    BFT_FREE(*helper);
}

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
                               MPI_Comm                     comm)
{
  fvm_writer_field_helper_t *h = helper;

  h->input_size = 0; /* Previously computed, so reset */
  h->output_size = 0;

  /* Get info on the current MPI communicator */

  if (comm != MPI_COMM_NULL) {
    h->min_rank_step = min_rank_step;
    h->min_block_size = min_block_size;
    h->comm = comm;
    MPI_Comm_rank(comm, &(h->rank));
    MPI_Comm_size(comm, &(h->n_ranks));
  }

  if (h->n_ranks < 2)
    h->rank = -1;
}

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
                                 size_t  *min_output_buffer_size)
{
  const fvm_writer_field_helper_t *h = helper;

  assert(h != NULL);

  if (input_size != NULL)
    *input_size = h->input_size;
  if (output_size != NULL)
    *output_size = h->output_size;

  if (min_output_buffer_size != NULL) {

    size_t min_size = 0;

    if (h->n_sub_elements_max > 1) {
      min_size = h->n_sub_elements_max * FVM_WRITER_MIN_SUB_ELEMENTS;
      min_size = CS_MIN(min_size, h->output_size);
    }

    if (h->output_size > 0)
      min_size = CS_MAX(FVM_WRITER_MIN_ELEMENTS, min_size);

    if (min_size > h->output_size)
      min_size = h->output_size;

    /* With interlaced multidimensional arrays,
       multiply buffer sizes by the stride */

    if (h->field_dim > 1 && h->interlace == CS_INTERLACE)
      min_size *= h->field_dim;

    *min_output_buffer_size = min_size;

  }

}

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
fvm_writer_field_helper_field_dim(const fvm_writer_field_helper_t  *helper)
{
  return helper->field_dim;
}

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
fvm_writer_field_helper_datatype(const fvm_writer_field_helper_t  *helper)
{
  return helper->datatype;
}

/*----------------------------------------------------------------------------
 * Partially distribute per element field values to an output buffer.
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
                                size_t                      *output_size)
{
  fvm_writer_field_helper_t *h = helper;

  int  retval = 0;

  cs_gnum_t   slice_output_size = 0;

  int  stride = 1;
  cs_lnum_t   end_id = 0;
  cs_lnum_t   num_shift = 0;

  size_t output_buffer_base_size = output_buffer_size;

  const fvm_nodal_section_t  *section = export_section->section;
  const cs_lnum_t   *parent_entity_num = section->parent_element_num;

  /* If output data is interlaced, set stride */

  if (h->field_dim > 1 && h->interlace == CS_INTERLACE) {
    stride = h->field_dim;
    output_buffer_base_size /= h->field_dim;
  }

  if (n_parent_lists == 0)
    num_shift = export_section->num_shift;

  /* If slice end has not been attained yet, extract values,
     update, and exit */

  if (h->start_id < section->n_elements) {

    /* Standard section */
    /*------------------*/

    if (export_section->type == section->type) {

      end_id = CS_MIN(h->start_id + (cs_lnum_t)output_buffer_base_size,
                      section->n_elements);

      fvm_convert_array(src_dim,
                        src_dim_shift,
                        stride,
                        h->start_id + num_shift,
                        end_id + num_shift,
                        src_interlace,
                        datatype,
                        h->datatype,
                        n_parent_lists,
                        parent_num_shift,
                        parent_entity_num,
                        field_values,
                        output_buffer);

      slice_output_size = end_id - h->start_id;

    }

    /* Tesselated section */
    /*--------------------*/

    else if (export_section->type != section->type) {

      cs_lnum_t   n_sub_elements_max = 0;

      const fvm_tesselation_t  *tesselation = section->tesselation;
      const cs_lnum_t *sub_element_idx
        = fvm_tesselation_sub_elt_index(section->tesselation,
                                        export_section->type);

      cs_lnum_t   output_buffer_size_min
        = fvm_tesselation_n_sub_elements(section->tesselation,
                                         export_section->type);

      fvm_tesselation_get_global_size(section->tesselation,
                                      export_section->type,
                                      NULL,
                                      &n_sub_elements_max);

      output_buffer_size_min = CS_MIN(output_buffer_size_min,
                                      (  n_sub_elements_max
                                       * FVM_WRITER_MIN_SUB_ELEMENTS));

      /* The buffer passed to this function should not be too small,
         as its size should have been computed in a consistent manner,
         but we check here in case it was not defined correctly by the
         calling writer. */

      if ((size_t)output_buffer_size_min > output_buffer_base_size) {
        bft_error(__FILE__, __LINE__, 0,
                  _("Output buffer too small:\n"
                    "Current size = %lu, minimum size required = %lu."),
                  (unsigned long)output_buffer_size,
                  (unsigned long)(output_buffer_size_min * stride));
      }

      /* compute index end to fill output_buffer after distribution */

      for (end_id = h->start_id;
           (   end_id < section->n_elements
            &&  (sub_element_idx[end_id] <   (cs_lnum_t)output_buffer_base_size
                                           + sub_element_idx[h->start_id]));
           end_id++);
      if (  sub_element_idx[end_id] - sub_element_idx[h->start_id]
          > (cs_lnum_t)output_buffer_base_size)
        end_id--;

      /* Extract data from parent; note that the following call to
         fvm_tesselation_distribute() will modify the var_buffer
         array (which avoids requiring an additional buffer) */

      fvm_convert_array(src_dim,
                        src_dim_shift,
                        stride,
                        h->start_id + num_shift,
                        end_id + num_shift,
                        src_interlace,
                        datatype,
                        h->datatype,
                        n_parent_lists,
                        parent_num_shift,
                        parent_entity_num,
                        field_values,
                        output_buffer);

      /* distribute data to tesselation and write values */

      fvm_tesselation_distribute(tesselation,
                                 export_section->type,
                                 h->start_id,
                                 end_id,
                                 cs_datatype_size[h->datatype] * stride,
                                 output_buffer);

      slice_output_size = (  sub_element_idx[end_id]
                           - sub_element_idx[h->start_id]);

    } /* End for tesselated section */

    h->start_id = end_id;

  } /* End for this slice */

  /* If we have reached the end of the section, reset indexes */

  else {

    h->last_section = export_section;
    h->start_id     = 0;

    retval = 1;

  }

  /* Set return values */

  *output_size = slice_output_size * stride;

  return retval;
}

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
                                size_t                      *output_size)
{
  fvm_writer_field_helper_t *h = helper;

  int  retval = 0;

  cs_gnum_t   slice_output_size = 0;
  cs_lnum_t   end_id = 0;

  int  stride = 1;

  const cs_lnum_t   *parent_entity_num = mesh->parent_vertex_num;

  /* If output data is interlaced, set stride */

  if (h->field_dim > 1 && h->interlace == CS_INTERLACE)
    stride = h->field_dim;

  /* Main vertices */
  /*---------------*/

  if (h->start_id < mesh->n_vertices) {

    cs_lnum_t n_entities = mesh->n_vertices;

    end_id = h->start_id + (output_buffer_size / stride);
    end_id = CS_MIN(end_id, n_entities);

    fvm_convert_array(src_dim,
                      src_dim_shift,
                      stride,
                      h->start_id,
                      end_id,
                      src_interlace,
                      datatype,
                      h->datatype,
                      n_parent_lists,
                      parent_num_shift,
                      parent_entity_num,
                      field_values,
                      output_buffer);

    slice_output_size = end_id - h->start_id;

    h->start_id = end_id;

  }

  /* Additional vertices in case of tesselation */
  /*--------------------------------------------*/

  else if (h->start_id < mesh->n_vertices + h->n_vertices_add) {

    cs_lnum_t start_id_prev = h->start_id;

    for (int i = 0; i < mesh->n_sections; i++) {

      const fvm_nodal_section_t  *const  section = mesh->sections[i];

      if (   section->type == FVM_CELL_POLY
          && section->tesselation != NULL) {

        cs_lnum_t n_vertices_section
          = fvm_tesselation_n_vertices_add(section->tesselation);

        /* Find current section */

        if (   n_vertices_section > 0
            && h->start_id < start_id_prev + n_vertices_section) {

          /* Temporarily shift h->start_id to section value */

          h->start_id -= start_id_prev;

          end_id = h->start_id + (output_buffer_size / stride);
          end_id = CS_MIN(end_id, h->start_id + n_vertices_section);

          slice_output_size = end_id - h->start_id;

          if (   (h->datatype == CS_DOUBLE || h->datatype == CS_FLOAT)
              && (datatype == CS_DOUBLE || datatype == CS_FLOAT)) {

            fvm_tesselation_vertex_values(section->tesselation,
                                          src_dim,
                                          src_dim_shift,
                                          stride,
                                          h->start_id,
                                          end_id,
                                          src_interlace,
                                          datatype,
                                          h->datatype,
                                          n_parent_lists,
                                          parent_num_shift,
                                          mesh->parent_vertex_num,
                                          field_values,
                                          output_buffer);

          }
          else
            _zero_values(slice_output_size*stride, datatype, output_buffer);

          /* Advance to next section if necessary */

          h->start_id = start_id_prev + end_id;
          start_id_prev = h->start_id;

        }

      }

    }

  }
  else {

    /* If we have reached the end of the section, reset start_id */

    h->start_id = 0;

  }

  /* Set return values */

  *output_size = slice_output_size * stride;

  if (slice_output_size == 0)
    retval = 1;

  return retval;
}

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
                                 fvm_writer_field_output_t   *output_func)
{
  assert(helper != NULL);

#if defined(HAVE_MPI)

  if (helper->n_ranks > 1)
    return _field_helper_output_eg(helper,
                                   context,
                                   export_section,
                                   src_dim,
                                   src_interlace,
                                   comp_order,
                                   n_parent_lists,
                                   parent_num_shift,
                                   datatype,
                                   field_values,
                                   output_func);

#endif /* defined(HAVE_MPI) */

  return _field_helper_output_el(helper,
                                 context,
                                 export_section,
                                 src_dim,
                                 src_interlace,
                                 comp_order,
                                 n_parent_lists,
                                 parent_num_shift,
                                 datatype,
                                 field_values,
                                 output_func);
}

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
                                 fvm_writer_field_output_t  *output_func)
{
  assert(helper != NULL);

#if defined(HAVE_MPI)

  if (helper->n_ranks > 1)
    _field_helper_output_ng(helper,
                            context,
                            mesh,
                            src_dim,
                            src_interlace,
                            comp_order,
                            n_parent_lists,
                            parent_num_shift,
                            datatype,
                            field_values,
                            output_func);

#endif /* defined(HAVE_MPI) */

  if (helper->n_ranks  == 1)
    _field_helper_output_nl(helper,
                            context,
                            mesh,
                            src_dim,
                            src_interlace,
                            comp_order,
                            n_parent_lists,
                            parent_num_shift,
                            datatype,
                            field_values,
                            output_func);
}

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
                                int      component_id)
{
  static const char *comp3[] = {"X", "Y", "Z"};
  static const char *comp6[] = {"XX", "YY", "ZZ", "XY", "XZ", "YZ"};
  static const char *comp9[] = {"XX", "XY", "XZ",
                                "YX", "YY", "YZ",
                                "ZX", "ZY", "ZZ"};

  s[0] = '\0';

  if (   dimension > 1 && s_size > 1
      && component_id > -1
      && component_id < dimension) {

    if (dimension == 3)
      strcpy(s, comp3[component_id]);

    else if (s_size > 2) {
      if (dimension == 6)
        strcpy(s, comp6[component_id]);
      else if (dimension == 9)
        strcpy(s, comp9[component_id]);
    }

    /* Fallback */

    if (s[0] == '\0') {
      snprintf(s, s_size, "%d", component_id);
      s[s_size -1] = '\0';
    }

    if (lowercase) {
      size_t l = strlen(s);
      for (size_t i = 0; i < l; i++)
        s[i] = tolower(s[i]);
    }
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

END_C_DECLS
