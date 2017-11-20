/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to CoProcessing stats daemon output
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Statistics library header
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

#undef BUILD_WITH_MPI
#define BUILD_WITH_MPI 1

#include "melissa_api.h"

#else

#include "melissa_api_no_mpi.h"

#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_defs.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_writer_helper.h"
#include "fvm_writer_priv.h"

#include "cs_block_dist.h"
#include "cs_file.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_melissa.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Melissa output writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char        *name;               /* Writer name */

  int          rank;               /* Rank of current process in communicator */
  int          n_ranks;            /* Number of processes in communicator */
  int          parameters[2];      /* array of parameters */

  int          nt_cur;             /* current time step */
  int          time_stamp;         /* current time step */

#if defined(HAVE_MPI)
  int          min_rank_step;      /* Minimum rank step */
  int          min_block_size;     /* Minimum block buffer size */
  MPI_Comm     block_comm;         /* Associated MPI block communicator */
  MPI_Comm     comm;               /* Associated MPI communicator */
#endif

} fvm_to_melissa_writer_t;

/*----------------------------------------------------------------------------
 * Context structure for fvm_writer_field_helper_output_* functions.
 *----------------------------------------------------------------------------*/

typedef struct {

  fvm_to_melissa_writer_t  *writer;      /* Pointer to writer structure */
  const char               *name;        /* current field name */

} _melissa_context_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Write block of a vector of doubles to Melissa output
 *
 * parameters:
 *   n_values <-- number of values to write
 *   num_end  <-- global number of past the last element for this block
 *   values   <-- pointer to values block array
 *   f        <-- file to write to
 *----------------------------------------------------------------------------*/

static void
_write_block_doubles_l(size_t             n_values,
                      const double        values[],
                      _melissa_context_t *c)
{
  /* TODO send(values, n_values); */

  static int is_init = 0;
  int n = (int)n_values;
  int coupling = 1;

  if (is_init == 0) {
#if defined(HAVE_MPI)
    melissa_init(&n,
                 &c->writer->n_ranks,
                 &c->writer->rank,
                 &c->writer->parameters[0],
                 &c->writer->parameters[1],
                 &c->writer->comm,
                 &coupling);
#else
    melissa_init_no_mpi (&n);
#endif
    is_init = 1;
  }

  if (c->writer->time_stamp > 0) {
    melissa_send(&(c->writer->time_stamp),
                 c->name,
                 values,
                 &c->writer->rank,
                 &c->writer->parameters[0],
                 &c->writer->parameters[1]);
  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Output function for field values.
 *
 * This function is passed to fvm_writer_field_helper_output_* functions.
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

static void
_field_output_g(void           *context,
                cs_datatype_t   datatype,
                int             dimension,
                int             component_id,
                cs_gnum_t       block_start,
                cs_gnum_t       block_end,
                void           *buffer)
{
  CS_UNUSED(dimension);
  CS_UNUSED(component_id);

  _melissa_context_t *c = context;

  int local_vect_size = block_end - block_start;

  assert(datatype == CS_DOUBLE);

  static int is_init = 0;
  int coupling = 1;

  if (is_init == 0) {
#if defined(HAVE_MPI)
    melissa_init(&local_vect_size,
                 &c->writer->n_ranks,
                 &c->writer->rank,
                 &c->writer->parameters[0],
                 &c->writer->parameters[1],
                 &c->writer->comm,
                 &coupling);
#else
    melissa_init_no_mpi(&local_vect_size);
#endif
    is_init = 1;
  }

  if (c->writer->time_stamp > 0) {
    melissa_send(&(c->writer->time_stamp),
                 c->name,
                 buffer,
                 &c->writer->rank,
                 &c->writer->parameters[0],
                 &c->writer->parameters[1]);
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write field values associated with nodal values of a nodal mesh to
 * an Melissa output in serial mode.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   n_entities         <-- number of entities
 *   input_dim          <-- input field dimension
 *   output_dim         <-- output field dimension
 *   interlace          <-- indicates if field in memory is interlaced
 *   n_parent_lists     <-- indicates if field values are to be obtained
 *                          directly through the local entity index (when 0) or
 *                          through the parent entity numbers (when 1 or more)
 *   parent_num_shift   <-- parent list to common number index shifts;
 *                          size: n_parent_lists
 *   datatype           <-- input data type (output is real)
 *   field_values       <-- array of associated field value arrays
 *   f                  <-- associated file handle
 *----------------------------------------------------------------------------*/

static void
_export_field_values_nl(const fvm_nodal_t           *mesh,
                        fvm_writer_field_helper_t   *helper,
                        int                          input_dim,
                        cs_interlace_t               interlace,
                        int                          n_parent_lists,
                        const cs_lnum_t              parent_num_shift[],
                        cs_datatype_t                datatype,
                        const void            *const field_values[],
                        void                        *context)
{
  int  i;
  size_t  output_size;
  double  *output_buffer;

  int output_dim = fvm_writer_field_helper_field_dim(helper);

  const size_t  output_buffer_size
    = mesh->n_vertices >  16 ? (mesh->n_vertices / 4) : mesh->n_vertices;

  BFT_MALLOC(output_buffer, output_buffer_size, double);

  for (i = 0; i < output_dim; i++) {

    const int i_in = i;

    while (fvm_writer_field_helper_step_nl(helper,
                                           mesh,
                                           input_dim,
                                           i_in,
                                           interlace,
                                           n_parent_lists,
                                           parent_num_shift,
                                           datatype,
                                           field_values,
                                           output_buffer,
                                           output_buffer_size,
                                           &output_size) == 0)
      _write_block_doubles_l(output_buffer_size, output_buffer, context);
  }

  BFT_FREE(output_buffer);
}

/*----------------------------------------------------------------------------
 * Write field values associated with element values of a nodal mesh to
 * an Melissa output.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   export_section   <-- pointer to EnSight section helper structure
 *   helper           <-- pointer to general writer helper structure
 *   input_dim        <-- input field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   field_values     <-- array of associated field value arrays
 *   f                <-- associated file handle
 *----------------------------------------------------------------------------*/

static const fvm_writer_section_t *
_export_field_values_el(const fvm_writer_section_t      *export_section,
                        fvm_writer_field_helper_t       *helper,
                        int                              input_dim,
                        cs_interlace_t                   interlace,
                        int                              n_parent_lists,
                        const cs_lnum_t                  parent_num_shift[],
                        cs_datatype_t                    datatype,
                        const void                *const field_values[],
                        void                            *context)
{
  int  i;
  size_t  input_size = 0, output_size = 0;
  size_t  min_output_buffer_size = 0, output_buffer_size = 0;
  double  *output_buffer = NULL;

  const fvm_writer_section_t  *current_section = NULL;

  int output_dim = fvm_writer_field_helper_field_dim(helper);

  /* Blocking for arbitrary buffer size, but should be small enough
     to add little additional memory requirement (in proportion), large
     enough to limit number of write calls. */

  fvm_writer_field_helper_get_size(helper,
                                   &input_size,
                                   &output_size,
                                   &min_output_buffer_size);

  output_buffer_size = input_size;
  output_buffer_size = CS_MAX(output_buffer_size, min_output_buffer_size);
  output_buffer_size = CS_MAX(output_buffer_size, 128);
  output_buffer_size = CS_MIN(output_buffer_size, output_size);

  BFT_MALLOC(output_buffer, output_buffer_size, double);

  /* Loop on dimension (de-interlace vectors, always 3D for EnSight) */

  for (i = 0; i < output_dim; i++) {

    bool loop_on_sections = true;

    const int i_in = i;

    current_section = export_section;

    while (loop_on_sections == true) {

      while (fvm_writer_field_helper_step_el(helper,
                                             current_section,
                                             input_dim,
                                             i_in,
                                             interlace,
                                             n_parent_lists,
                                             parent_num_shift,
                                             datatype,
                                             field_values,
                                             output_buffer,
                                             output_buffer_size,
                                             &output_size) == 0)
        _write_block_doubles_l(output_buffer_size, output_buffer, context);

      current_section = current_section->next;

      if (   current_section == NULL
          || current_section->continues_previous == false)
        loop_on_sections = false;

    } /* while (loop on sections) */

  } /* end of loop on spatial dimension */

  BFT_FREE(output_buffer);

  return current_section;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to Melissa output writer.
 *
 * No options are currently handled
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque Melissa output writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_melissa_init_writer(const char             *name,
                           const char             *path,
                           const char             *options,
                           fvm_writer_time_dep_t   time_dependency,
                           MPI_Comm                comm)
#else
void *
fvm_to_melissa_init_writer(const char             *name,
                           const char             *path,
                           const char             *options,
                           fvm_writer_time_dep_t   time_dependency)
#endif
{
  CS_UNUSED(path);
  CS_UNUSED(time_dependency);

  fvm_to_melissa_writer_t  *w = NULL;

  /* Initialize writer */

  BFT_MALLOC(w, 1, fvm_to_melissa_writer_t);

  BFT_MALLOC(w->name, strlen(name) + 1, char);
  strcpy(w->name, name);

  w->rank = 0;
  w->n_ranks = 1;

#if defined(HAVE_MPI)
  {
    int mpi_flag, rank, n_ranks, min_rank_step, min_block_size;
    MPI_Comm w_block_comm, w_comm;
    w->min_rank_step = 1;
    w->min_block_size = 1024*1024*8;
    w->block_comm = MPI_COMM_NULL;
    w->comm = MPI_COMM_NULL;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag && comm != MPI_COMM_NULL) {
      w->comm = comm;
      MPI_Comm_rank(w->comm, &rank);
      MPI_Comm_size(w->comm, &n_ranks);
      w->rank = rank;
      w->n_ranks = n_ranks;
      cs_file_get_default_comm(&min_rank_step, &min_block_size,
                               &w_block_comm, &w_comm);
      if (comm == w_comm) {
        w->min_rank_step = min_rank_step;
        w->min_block_size = min_block_size;
        w->block_comm = w_block_comm;
      }
      w->comm = comm;
    }
  }
#endif /* defined(HAVE_MPI) */

  /* Parse options */

  w->parameters[0] = 0;
  w->parameters[1] = 0;

  if (options != NULL) {

    size_t l = strlen(options);

    if (l > 0) {

      /* Add metadata */

      char       *options_c = NULL;
      char       *temp_char = NULL;
      const char  s[2] = " ";

      BFT_MALLOC(options_c, l+1, char);
      strncpy(options_c, options, l);
      options_c[l] = '\0';

      /* get the first token */
      temp_char = strtok(options_c, s);

      /* walk through other tokens */
      for (int i = 0; i < 2 && temp_char != NULL; i++) {
        w->parameters[i] = atoi(temp_char);
        temp_char = strtok(NULL, s);
      }

      BFT_FREE(options_c);

    }
  }

  w->nt_cur = -1;
  w->time_stamp = -1;

  /* Return writer */
  return w;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to Melissa output writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque Ensight Gold writer structure.
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_melissa_finalize_writer(void  *this_writer_p)
{
  fvm_to_melissa_writer_t *w = (fvm_to_melissa_writer_t *)this_writer_p;
  static int disconnected = 0;

  if (disconnected == 0) {
    BFT_FREE(w->name);
    melissa_finalize();
    BFT_FREE(w);
    disconnected = 1;
  }
  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with an EnSight geometry.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_melissa_set_mesh_time(void          *this_writer_p,
                             const int      time_step,
                             const double   time_value)
{
  fvm_to_melissa_writer_t  *w = (fvm_to_melissa_writer_t *)this_writer_p;
  if (time_step > w->nt_cur) {
    w->nt_cur = time_step;
    w->time_stamp += 1;
  }
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a an Melissa output
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *----------------------------------------------------------------------------*/

void
fvm_to_melissa_export_nodal(void               *this_writer_p,
                            const fvm_nodal_t  *mesh)
{

}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a Melissa output.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *   mesh             <-- pointer to associated nodal mesh structure
 *   name             <-- variable name
 *   location         <-- variable definition location (nodes or elements)
 *   dimension        <-- variable dimension (0: constant, 1: scalar,
 *                        3: vector, 6: sym. tensor, 9: asym. tensor)
 *   interlace        <-- indicates if variable in memory is interlaced
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   time_step        <-- number of the current time step
 *   time_value       <-- associated time value
 *   field_values     <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
fvm_to_melissa_export_field(void                  *this_writer_p,
                            const fvm_nodal_t     *mesh,
                            const char            *name,
                            fvm_writer_var_loc_t   location,
                            int                    dimension,
                            cs_interlace_t         interlace,
                            int                    n_parent_lists,
                            const cs_lnum_t        parent_num_shift[],
                            cs_datatype_t          datatype,
                            int                    time_step,
                            double                 time_value,
                            const void      *const field_values[])
{
  int   output_dim;

  const fvm_writer_section_t  *export_section = NULL;
  fvm_writer_field_helper_t   *helper = NULL;
  fvm_writer_section_t        *export_list = NULL;
  fvm_to_melissa_writer_t     *w = (fvm_to_melissa_writer_t *)this_writer_p;

  const int  rank = w->rank;
  const int  n_ranks = w->n_ranks;

  if (time_step > w->nt_cur) {
    w->nt_cur = time_step;
    w->time_stamp += 1;
  }
  else if (time_step < w->nt_cur)
    return;

  /* Initialization */
  /*----------------*/

  output_dim = dimension;
  if (dimension != 1)
    return; /* currently handled only for scalars */

  /* Initialize writer helper */
  /*--------------------------*/

  /* Build list of sections that are used here, in order of output */
  export_list = fvm_writer_export_list(mesh,
                                       fvm_nodal_get_max_entity_dim(mesh),
                                       true,
                                       false,
                                       true,
                                       true,
                                       false,
                                       false);

  helper = fvm_writer_field_helper_create(mesh,
                                          export_list,
                                          output_dim,
                                          CS_NO_INTERLACE,
                                          CS_DOUBLE,
                                          location);

#if defined(HAVE_MPI)

  if (n_ranks > 1)
    fvm_writer_field_helper_init_g(helper,
                                   w->min_rank_step,
                                   w->min_block_size,
                                   w->comm);

#endif

  _melissa_context_t c;
  c.writer = w;
  c.writer->time_stamp = time_step;
  c.name = name;

  /* Per node variable */
  /*-------------------*/

  if (location == FVM_WRITER_PER_NODE) {

#if defined(HAVE_MPI)

    if (n_ranks > 1)
      fvm_writer_field_helper_output_n(helper,
                                       &c,
                                       mesh,
                                       dimension,
                                       interlace,
                                       NULL,
                                       n_parent_lists,
                                       parent_num_shift,
                                       datatype,
                                       field_values,
                                       _field_output_g);

#endif /* defined(HAVE_MPI) */

    if (n_ranks == 1)
      _export_field_values_nl(mesh,
                              helper,
                              dimension,
                              interlace,
                              n_parent_lists,
                              parent_num_shift,
                              datatype,
                              field_values,
                              &c);
  }

  /* Per element variable */
  /*----------------------*/

  else if (location == FVM_WRITER_PER_ELEMENT) {
    export_section = export_list;

    while (export_section != NULL) {

      /* Output per grouped sections */

#if defined(HAVE_MPI)

      if (n_ranks > 1)
        export_section = fvm_writer_field_helper_output_e(helper,
                                                          &c,
                                                          export_section,
                                                          dimension,
                                                          interlace,
                                                          NULL,
                                                          n_parent_lists,
                                                          parent_num_shift,
                                                          datatype,
                                                          field_values,
                                                          _field_output_g);

#endif /* defined(HAVE_MPI) */

      if (n_ranks == 1)
        export_section = _export_field_values_el(export_section,
                                                 helper,
                                                 dimension,
                                                 interlace,
                                                 n_parent_lists,
                                                 parent_num_shift,
                                                 datatype,
                                                 field_values,
                                                 &c);

    } /* End of loop on sections */

  } /* End for per element variable */

  /* Free helper structures */
  /*------------------------*/
  fvm_writer_field_helper_destroy(&helper);
  BFT_FREE(export_list);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
