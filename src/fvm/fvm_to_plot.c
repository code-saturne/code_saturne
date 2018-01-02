/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to plot files
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_defs.h"
#include "fvm_convert_array.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "bft_printf.h"
#include "fvm_writer_helper.h"
#include "fvm_writer_priv.h"

#include "cs_block_dist.h"
#include "cs_file.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_plot.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Type of 1D plot file format */

typedef enum {
  CS_PLOT_DAT,  /* .dat file (usable by Qtplot or Grace) */
  CS_PLOT_CSV   /* .csv file (readable by ParaView or spreadsheat) */
} cs_plot_format_t;

/*----------------------------------------------------------------------------
 * Plot writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char        *name;               /* Writer name */
  char        *path;               /* Path prefix */

  int          rank;               /* Rank of current process in communicator */
  int          n_ranks;            /* Number of processes in communicator */

  cs_plot_format_t  format;        /* Plot format */

  int               nt;            /* Time step */
  double            t;             /* Time value */

  int               n_cols;        /* Number of columns */
  int               n_cols_max;    /* Max. number of columns */
  int               n_rows;        /* Number of rows */

  cs_real_t        *buffer;        /* Values buffer */

  char             *file_name;     /* File name */
  FILE             *f;             /* Associated plots */

#if defined(HAVE_MPI)
  MPI_Comm     comm;               /* Associated MPI communicator */
#endif

} fvm_to_plot_writer_t;

/*----------------------------------------------------------------------------
 * Context structure for fvm_writer_field_helper_output_* functions.
 *----------------------------------------------------------------------------*/

typedef struct {

  fvm_to_plot_writer_t  *writer;    /* Pointer to writer structure */

  const char            *name;      /* Field name */

} _plot_context_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
_field_output(void           *context,
              cs_datatype_t   datatype,
              int             dimension,
              int             component_id,
              cs_gnum_t       block_start,
              cs_gnum_t       block_end,
              void           *buffer)
{
  CS_UNUSED(datatype);

  _plot_context_t *c = context;

  fvm_to_plot_writer_t  *w = c->writer;

  if (w->rank > 0)
    return;

  const int n_rows = (block_end > block_start) ? block_end - block_start : 0;

  if (w->n_cols == 0)
    w->n_rows = n_rows;
  else if (w->n_rows != n_rows) {
    const char e[] = "";
    const char *name = (c->name != NULL) ? c->name : e;
    bft_printf(_("Warning: inconsistent data size for plot \"%s\" between\n"
                 "field \"%s\" and previous outputs; values dropped.\n"),
               w->name, name);
    return;
  }

  /* Create file and write headers on first output */

  if (w->f == NULL) {

    char t_stamp[32];
    if (w->nt >= 0)
      sprintf(t_stamp, "_%.4i", w->nt);
    else
      t_stamp[0] = '\0';
    size_t l =   strlen(w->path) + strlen(w->name)
               + strlen(t_stamp) + 4 + 1;
    BFT_REALLOC(w->file_name, l, char);

    if (w->format == CS_PLOT_DAT)
      sprintf(w->file_name, "%s%s%s.dat", w->path, w->name, t_stamp);
    else
      sprintf(w->file_name, "%s%s%s.csv", w->path, w->name, t_stamp);

    w->f = fopen(w->file_name, "w");

    if (w->f ==  NULL) {
      bft_error(__FILE__, __LINE__, errno,
                _("Error opening file: \"%s\""), w->file_name);
      return;
    }

    if (w->format == CS_PLOT_DAT) {
      fprintf(w->f, "# Code_Saturne plot output\n#\n");
      if (w->nt < 0)
        fprintf(w->f, "# time independent\n");
      else {
        fprintf(w->f, "# time step id: %i\n", w->nt);
        fprintf(w->f, "# time:         %12.5e\n#\n", w->t);
      }
      fprintf(w->f, "#COLUMN_TITLES: ");
    }

  }

  assert(component_id == 0);

  if (w->n_cols + dimension > w->n_cols_max) {
    while (w->n_cols + dimension > w->n_cols_max) {
      if (w->n_cols_max == 0)
        w->n_cols_max = 4;
      else
        w->n_cols_max *= 2;
    }
    BFT_REALLOC(w->buffer, w->n_rows*w->n_cols_max, cs_real_t);
  }

  for (int i = 0; i < dimension; i++) {

    /* Add header info */

    char name_buf[64];

    if (c->name != NULL)
      strncpy(name_buf, c->name, 63);
    else
      name_buf[0] = '\0';
    name_buf[63] = '\0';

    if (dimension > 1) {
      size_t l = strlen(name_buf);
      if (l > 59)
        l = 59;
      if (l > 0)
        name_buf[l++] = '_';
      fvm_writer_field_component_name(name_buf+l, 3, true, dimension, i);
    }

    if (w->format == CS_PLOT_DAT) {
      if (w->n_cols > 0)
        fprintf(w->f, " | %s", name_buf);
      else
        fprintf(w->f, " %s", name_buf);
    }
    else if (w->format == CS_PLOT_CSV) {
      if (w->n_cols > 0)
        fprintf(w->f, ", %s", name_buf);
      else
        fprintf(w->f, "%s", name_buf);
    }

    /* Update buffer */

    const cs_real_t *src = buffer;
    cs_real_t *dest = w->buffer + w->n_rows*w->n_cols;

    for (cs_lnum_t j = 0; j < w->n_rows; j++)
      dest[j] = src[j*dimension + i];
    w->n_cols++;

  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to plot file writer.
 *
 * Options are:
 *   csv                 output CSV (comma-separated-values) files
 *   dat                 output dat (space-separated) files
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque plot writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_plot_init_writer(const char             *name,
                        const char             *path,
                        const char             *options,
                        fvm_writer_time_dep_t   time_dependency,
                        MPI_Comm                comm)
#else
void *
fvm_to_plot_init_writer(const char             *name,
                        const char             *path,
                        const char             *options,
                        fvm_writer_time_dep_t   time_dependency)
#endif
{
  CS_UNUSED(time_dependency);

  fvm_to_plot_writer_t  *w = NULL;

  /* Initialize writer */

  BFT_MALLOC(w, 1, fvm_to_plot_writer_t);

  BFT_MALLOC(w->name, strlen(name) + 1, char);
  strcpy(w->name, name);

  BFT_MALLOC(w->path, strlen(path) + 1, char);
  strcpy(w->path, path);

  w->rank = 0;
  w->n_ranks = 1;

#if defined(HAVE_MPI)
  {
    int mpi_flag, rank, n_ranks;
    w->comm = MPI_COMM_NULL;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag && comm != MPI_COMM_NULL) {
      w->comm = comm;
      MPI_Comm_rank(w->comm, &rank);
      MPI_Comm_size(w->comm, &n_ranks);
      w->rank = rank;
      w->n_ranks = n_ranks;
    }
  }
#endif /* defined(HAVE_MPI) */

  /* Defaults */

  w->format = CS_PLOT_CSV;

  w->nt = -1;
  w->t = -1;

  w->n_cols = 0;
  w->n_cols_max = 0;
  w->n_rows = 0;

  w->buffer = NULL;

  w->file_name = NULL;
  w->f = NULL;

  /* Parse options */

  if (options != NULL) {

    int i1, i2, l_opt;
    int l_tot = strlen(options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {

      for (i2 = i1; i2 < l_tot && options[i2] != ' '; i2++);
      l_opt = i2 - i1;

      if ((l_opt == 3) && (strncmp(options + i1, "csv", l_opt) == 0))
        w->format = CS_PLOT_CSV;
      else if ((l_opt == 3) && (strncmp(options + i1, "dat", l_opt) == 0))
        w->format = CS_PLOT_DAT;

      for (i1 = i2 + 1 ; i1 < l_tot && options[i1] == ' ' ; i1++);

    }

  }

  /* Return writer */

  return w;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to plot file writer.
 *
 * parameters:
 *   writer <-- pointer to opaque plot Gold writer structure.
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_plot_finalize_writer(void  *writer)
{
  fvm_to_plot_writer_t  *w
    = (fvm_to_plot_writer_t *)writer;

  BFT_FREE(w->name);
  BFT_FREE(w->path);

  fvm_to_plot_flush(writer);

  BFT_FREE(w->file_name);

  BFT_FREE(w);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with an plot geometry.
 *
 * parameters:
 *   writer     <-- pointer to associated writer
 *   time_step  <-- time step number
 *   time_value <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_plot_set_mesh_time(void    *writer,
                          int      time_step,
                          double   time_value)
{
  fvm_to_plot_writer_t  *w = (fvm_to_plot_writer_t *)writer;

  w->nt = time_step;
  w->t = time_value;

  if (w->n_cols > 0)
    fvm_to_plot_flush(writer);
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a an plot file
 *
 * parameters:
 *   writer <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *----------------------------------------------------------------------------*/

void
fvm_to_plot_export_nodal(void               *writer,
                         const fvm_nodal_t  *mesh)
{
#if 0
  if (fvm_nodal_get_max_entity_dim(mesh) == 0) {

    fvm_to_plot_writer_t  *w
      = (fvm_to_plot_writer_t *)writer;

    fvm_writer_field_helper_t  *helper
      = fvm_writer_field_helper_create(mesh,
                                       NULL, /* section list */
                                       mesh->dim,
                                       CS_INTERLACE,
                                       CS_REAL_TYPE,
                                       FVM_WRITER_PER_NODE);

#if defined(HAVE_MPI)
    if (w->n_ranks > 1)
      fvm_writer_field_helper_init_g(helper,
                                     w->n_ranks,
                                     0,
                                     w->comm);
#endif

    _plot_context_t c = {.writer = w};

    int n_parent_lists = (mesh->parent_vertex_num != NULL) ? 1 : 0;
    cs_lnum_t parent_num_shift[1] = {0};
    const cs_real_t  *coo_ptr[1] = {mesh->vertex_coords};

    fvm_writer_field_helper_output_n(helper,
                                     &c,
                                     mesh,
                                     mesh->dim,
                                     CS_INTERLACE,
                                     NULL,
                                     n_parent_lists,
                                     parent_num_shift,
                                     CS_COORD_TYPE,
                                     (const void **)coo_ptr,
                                     _field_output);

  }

  /* Free helper structures */

  fvm_writer_field_helper_destroy(&helper);
#endif
}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to an plot file.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   writer    <-- pointer to associated writer
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
fvm_to_plot_export_field(void                  *writer,
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
  fvm_to_plot_writer_t  *w = (fvm_to_plot_writer_t *)writer;

  /* If time step changes, update it */

  if (time_step != w->nt)
    fvm_to_plot_set_mesh_time(writer,
                              time_step,
                              time_value);

  /* Initialize writer helper */

  cs_datatype_t  dest_datatype = CS_REAL_TYPE;

  if (datatype >= CS_INT32 && datatype <= CS_UINT64)
    dest_datatype = CS_INT64;

  fvm_writer_field_helper_t  *helper
    = fvm_writer_field_helper_create(mesh,
                                     NULL, /* section list */
                                     dimension,
                                     CS_NO_INTERLACE,
                                     dest_datatype,
                                     location);

#if defined(HAVE_MPI)

  if (w->n_ranks > 1)
    fvm_writer_field_helper_init_g(helper,
                                   w->n_ranks,
                                   0,
                                   w->comm);

#endif

  /* Per node variable only */

  if (location == FVM_WRITER_PER_NODE) {

    _plot_context_t c = {.writer = w, .name = name};

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
                                     _field_output);

  }

  /* Free helper structures */

  fvm_writer_field_helper_destroy(&helper);
}

/*----------------------------------------------------------------------------
 * Flush files associated with a given writer.
 *
 * In this case, the effective writing to file is done.
 *
 * parameters:
 *   writer <-- pointer to associated writer
 *----------------------------------------------------------------------------*/

void
fvm_to_plot_flush(void  *writer)
{
  fvm_to_plot_writer_t  *w = (fvm_to_plot_writer_t *)writer;

  if (w->f != NULL && w->buffer != NULL) {

    /* Transpose output on write */

    int n_cols = w->n_cols;
    int n_rows = w->n_rows;

    if (w->format == CS_PLOT_DAT) {

      fprintf(w->f, "\n");

      for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols - 1; j++)
          fprintf(w->f, "%12.5e ", w->buffer[w->n_rows*j + i]);
        if (n_cols > 0) {
          int j = n_cols -1;
          fprintf(w->f, "%12.5e\n", w->buffer[w->n_rows*j + i]);
        }
      }

    }
    else if (w->format == CS_PLOT_CSV) {

      fprintf(w->f, "\n");

      for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols-1; j++)
          fprintf(w->f, "%12.5e, ", w->buffer[w->n_rows*j + i]);
        if (n_cols > 0) {
          int j = n_cols -1;
          fprintf(w->f, "%12.5e\n", w->buffer[w->n_rows*j + i]);
        }
      }

    }

    w->n_rows = 0;
    w->n_cols = 0;
    w->n_cols_max = 0;

    if (fclose(w->f) != 0)
      bft_error(__FILE__, __LINE__, errno,
                _("Error closing file: \"%s\""), w->file_name);

    w->f = NULL;

  }

  BFT_FREE(w->buffer);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
