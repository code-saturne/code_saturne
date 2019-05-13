/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to time plot files
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
#include "fvm_writer_helper.h"
#include "fvm_writer_priv.h"

#include "cs_block_dist.h"
#include "cs_file.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"
#include "cs_time_plot.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_time_plot.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * time plot writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char        *name;               /* Writer name */
  char        *prefix;             /* Plot file prefix */

  int          rank;               /* Rank of current process in communicator */
  int          n_ranks;            /* Number of processes in communicator */

  cs_time_plot_format_t  format;   /* Time plot format */

  float             flush_wtime;   /* Elapsed time interval between file
                                      flushes (if < 0, no forced flush) */
  int               n_buf_steps;   /* Number of time steps in output
                                      buffer if file is not to be kept open */

  bool              use_iteration; /* Should we use the iteration number
                                      instead of the physical time ? */

  int               nt;            /* Time step */
  double            t;             /* Time value */

  int               n_plots;       /* Number of associated plots */

  cs_map_name_to_id_t   *f_map;    /* field names to plots mapping */
  cs_time_plot_t  **tp;            /* Associated plots */

#if defined(HAVE_MPI)
  MPI_Comm     comm;               /* Associated MPI communicator */
#endif

} fvm_to_time_plot_writer_t;

/*----------------------------------------------------------------------------
 * Context structure for fvm_writer_field_helper_output_* functions.
 *----------------------------------------------------------------------------*/

typedef struct {

  fvm_to_time_plot_writer_t  *writer;    /* Pointer to writer structure */

  const fvm_nodal_t          *mesh;      /* Pointer to mesh */
  const char                 *name;      /* Field name */

} _time_plot_context_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Output function for coordinates files
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
_coords_output(void           *context,
               cs_datatype_t   datatype,
               int             dimension,
               int             component_id,
               cs_gnum_t       block_start,
               cs_gnum_t       block_end,
               void           *buffer)
{
  if (dimension > 3 || buffer == NULL)
    return;

  _time_plot_context_t *c = context;

  fvm_to_time_plot_writer_t  *w = c->writer;

  assert(datatype == CS_REAL_TYPE);
  assert(component_id == 0);

  char *file_name;
  FILE *_f;

  const cs_real_t *coords = buffer;
  const cs_lnum_t n_coords = block_end - block_start;

  char t_stamp[32];
  if (w->nt >= 0)
    sprintf(t_stamp, "_%.4i", w->nt);
  else
    t_stamp[0] = '\0';

  size_t l = strlen(w->prefix) + strlen("coords") + strlen(t_stamp) + 4 + 1;

  BFT_MALLOC(file_name, l, char);

  if (w->format == CS_TIME_PLOT_DAT)
    sprintf(file_name, "%scoords%s.dat", w->prefix, t_stamp);
  else if (w->format == CS_TIME_PLOT_CSV)
    sprintf(file_name, "%scoords%s.csv", w->prefix, t_stamp);

  _f = fopen(file_name, "w");
  if (_f == NULL) {
    bft_error(__FILE__, __LINE__, errno,
              _("Error opening file: \"%s\""), file_name);
    return;
  }

  /* whitespace-separated format */

  if (w->format == CS_TIME_PLOT_DAT) {

    const char **probe_names = fvm_nodal_get_global_vertex_labels(c->mesh);

    if (probe_names != NULL) {
      fprintf(_f, "# Monitoring point names:\n");
      for (cs_lnum_t i = 0; i < n_coords; i++)
        fprintf(_f, "#   %6i %16s\n",
                i+1, probe_names[i]);
      fprintf(_f, "#\n");
    }

    fprintf(_f, _("# Monitoring point coordinates:\n"));

    switch(dimension) {
    case 3:
      for (cs_lnum_t i = 0; i < n_coords; i++)
        fprintf(_f, "# %6i %14.7e %14.7e %14.7e\n",
                i + 1,
                coords[i*3],
                coords[i*3 + 1],
                coords[i*3 + 2]);
      break;
    case 2:
      for (cs_lnum_t i = 0; i < n_coords; i++)
        fprintf(_f, "# %6i %14.7e %14.7e\n",
                i + 1,
                coords[i*2],
                coords[i*2 + 1]);
      break;
    case 1:
      for (cs_lnum_t i = 0; i < n_coords; i++)
        fprintf(_f, "# %6i %14.7e\n",
                i + 1,
                coords[i]);
      break;
    }
    fprintf(_f, "#\n");

  }

  /* CSV format */

  else if (w->format == CS_TIME_PLOT_CSV) {

    switch(dimension) {
    case 3:
      fprintf(_f, "x, y, z\n");
      for (cs_lnum_t i = 0; i < n_coords; i++) {
        fprintf(_f, "%14.7e, %14.7e, %14.7e\n",
                coords[i*3], coords[i*3 + 1], coords[i*3 + 2]);
      }
      break;
    case 2:
      fprintf(_f, "x, y\n");
      for (cs_lnum_t i = 0; i < n_coords; i++) {
        fprintf(_f, "%14.7e, %14.7e\n",
                coords[i*2], coords[i*2 + 1]);
      }
      break;
    case 1:
      fprintf(_f, "x\n");
      for (cs_lnum_t i = 0; i < n_coords; i++) {
        fprintf(_f, "%14.7e\n", coords[i]);
      }
      break;
    }

  }

  /* Close file */

  if (fclose(_f) != 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Error closing file: \"%s\""), file_name);

  BFT_FREE(file_name);
}

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

  cs_time_plot_t *p;

  _time_plot_context_t       *c = context;
  fvm_to_time_plot_writer_t  *w = c->writer;

  if (buffer == NULL) return;

  assert(component_id == 0);

  int n_vals = block_end - block_start;

  cs_real_t *_vals = NULL;

  if (dimension > 1)
    BFT_MALLOC(_vals, n_vals, cs_real_t);

  for (int _component_id = 0; _component_id < dimension; _component_id++) {

    /* Build plot name */

    char tmpn[128], tmpe[6];

    char *plot_name = tmpn;

    fvm_writer_field_component_name(tmpe, 6, false, dimension, _component_id);

    size_t lce = strlen(tmpe);
    size_t l =  strlen(c->name) + 1;

    if (lce > 0)
      l += 2 + lce;

    if (l > 128)
      BFT_MALLOC(plot_name, l, char);

    if (lce > 0)
      sprintf(plot_name, "%s[%s]", c->name, tmpe);
    else
      strcpy(plot_name, c->name);

    int p_id = cs_map_name_to_id(w->f_map, plot_name);

    if (p_id >= w->n_plots) {

      w->n_plots += 1;
      BFT_REALLOC(w->tp, w->n_plots, cs_time_plot_t *);

      int n_probes = (block_end > block_start) ? block_end - block_start : 0;

      const char **probe_names = fvm_nodal_get_global_vertex_labels(c->mesh);

      w->tp[p_id] = cs_time_plot_init_probe(plot_name,
                                            w->prefix,
                                            w->format,
                                            w->use_iteration,
                                            w->flush_wtime,
                                            w->n_buf_steps,
                                            n_probes,
                                            NULL,
                                            NULL, /* probe_coords */
                                            probe_names);

    }

    if (plot_name != tmpn)
      BFT_FREE(plot_name);

    p = w->tp[p_id];

    if (p != NULL) {
      const cs_real_t *vals = (const cs_real_t *)buffer;
      if (dimension > 1) {
        for (int i = 0; i < n_vals; i++)
          _vals[i] = vals[i*dimension + _component_id];
        vals = _vals;
      }
      cs_time_plot_vals_write(p,
                              w->nt,
                              w->t,
                              n_vals,
                              vals);
    }
  }

  BFT_FREE(_vals);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to time plot file writer.
 *
 * Options are:
 *   csv                 output CSV (comma-separated-values) files
 *   dat                 output dat (space-separated) files
 *   use_iteration       use time step id instead of time value for
 *                       first column
 *   flush_wtime=<wt>    flush output file every 'wt' seconds
 *   n_buf_steps=<n>     write output to file every 'n' output steps
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque time plot writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_time_plot_init_writer(const char             *name,
                             const char             *path,
                             const char             *options,
                             fvm_writer_time_dep_t   time_dependency,
                             MPI_Comm                comm)
#else
void *
fvm_to_time_plot_init_writer(const char             *name,
                             const char             *path,
                             const char             *options,
                             fvm_writer_time_dep_t   time_dependency)
#endif
{
  CS_NO_WARN_IF_UNUSED(time_dependency);

  fvm_to_time_plot_writer_t  *w = NULL;

  /* Initialize writer */

  BFT_MALLOC(w, 1, fvm_to_time_plot_writer_t);

  BFT_MALLOC(w->name, strlen(name) + 1, char);
  strcpy(w->name, name);

  if (strlen(name) > 0) {
    BFT_MALLOC(w->prefix, strlen(path) + 1 + strlen(name) + 1, char);
    sprintf(w->prefix, "%s%s_", path, name);
  }
  else {
    BFT_MALLOC(w->prefix, strlen(path) + 1, char);
    strcpy(w->prefix, path);
  }

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

  w->format = CS_TIME_PLOT_CSV;

  cs_time_plot_get_flush_default(&(w->flush_wtime),
                                 &(w->n_buf_steps));

  w->use_iteration = false;

  w->nt = -1;
  w->t = -1;

  w->n_plots = 0;
  w->f_map = (w->rank > 0) ? NULL : cs_map_name_to_id_create();
  w->tp = NULL;

  /* Parse options */

  if (options != NULL) {

    int i1, i2, l_opt;
    int l_tot = strlen(options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {

      for (i2 = i1; i2 < l_tot && options[i2] != ' '; i2++);
      l_opt = i2 - i1;

      if ((l_opt == 3) && (strncmp(options + i1, "csv", l_opt) == 0))
        w->format = CS_TIME_PLOT_CSV;
      else if ((l_opt == 3) && (strncmp(options + i1, "dat", l_opt) == 0))
        w->format = CS_TIME_PLOT_DAT;
      else if ((l_opt == 13) && (strcmp(options + i1, "use_iteration") == 0))
        w->use_iteration = true;
      else if (strncmp(options + i1, "n_buf_steps=", 12) == 0) {
        const char *s = options + i1 + 12;
        int nb, nr;
        nr = sscanf(s, "%d", &nb);
        if (nr == 1)
          w->n_buf_steps = nb;
      }
      else if (strncmp(options + i1, "flush_wtime=", 12) == 0) {
        const char *s = options + i1 + 12;
        int nr;
        float wt;
        nr = sscanf(s, "%g", &wt);
        if (nr == 1)
          w->flush_wtime = wt;
      }

      for (i1 = i2 + 1 ; i1 < l_tot && options[i1] == ' ' ; i1++);

    }

  }

  /* Return writer */

  return w;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to time plot file writer.
 *
 * parameters:
 *   writer <-- pointer to opaque time plot writer structure.
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_time_plot_finalize_writer(void  *writer)
{
  fvm_to_time_plot_writer_t  *w  = (fvm_to_time_plot_writer_t *)writer;

  BFT_FREE(w->name);
  BFT_FREE(w->prefix);

  if (w->rank <= 0) {
    for (int i = 0; i < w->n_plots; i++)
      cs_time_plot_finalize(&(w->tp[i]));
    BFT_FREE(w->tp);
    cs_map_name_to_id_destroy(&(w->f_map));
  }

  BFT_FREE(w);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with a time plot geometry.
 *
 * parameters:
 *   writer     <-- pointer to associated writer
 *   time_step  <-- time step number
 *   time_value <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_time_plot_set_mesh_time(void          *writer,
                               const int      time_step,
                               const double   time_value)
{
  fvm_to_time_plot_writer_t  *w = (fvm_to_time_plot_writer_t *)writer;

  w->nt = time_step;
  w->t = time_value;
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a an time plot file
 *
 * parameters:
 *   writer <-- pointer to associated writer
 *   mesh   <-- pointer to nodal mesh structure that should be written
 *----------------------------------------------------------------------------*/

void
fvm_to_time_plot_export_nodal(void               *writer,
                              const fvm_nodal_t  *mesh)
{
  if (fvm_nodal_get_max_entity_dim(mesh) == 0) {

    fvm_to_time_plot_writer_t  *w
      = (fvm_to_time_plot_writer_t *)writer;

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

    _time_plot_context_t c = {.writer = w,
                              .mesh = mesh,
                              .name = NULL};

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
                                     _coords_output);

    fvm_writer_field_helper_destroy(&helper);
  }
}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to an time plot file.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   writer           <-- pointer to associated writer
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
fvm_to_time_plot_export_field(void                  *writer,
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
  fvm_to_time_plot_writer_t  *w = (fvm_to_time_plot_writer_t *)writer;

  /* If time step changes, update it */

  if (time_step != w->nt)
    fvm_to_time_plot_set_mesh_time(writer,
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
                                     CS_INTERLACE,
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

    _time_plot_context_t c = {.writer = w,
                              .mesh = mesh,
                              .name = name};

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

/*----------------------------------------------------------------------------*/

END_C_DECLS
