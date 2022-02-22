/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to CoProcessing stats daemon output
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

#if defined(HAVE_MELISSA_MPI_05)
  #include "melissa_api.h"
  #define HAVE_MELISSA_MPI 1
#elif defined(HAVE_MELISSA_MPI)
  #include "melissa/api.h"
#endif

#if defined(HAVE_MELISSA_NO_MPI)
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

#include "cs_map.h"
#include "cs_parall.h"

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

  bool         dry_run;            /* If true, do not connect to Melissa */
  FILE        *tracefile;          /* optional file for tracing */

  int          rank;               /* Rank of current process in communicator */
  int          n_ranks;            /* Number of processes in communicator */

  size_t       buffer_size;        /* buffer size required */

  cs_map_name_to_id_t  *f_map;     /* field names mapping */
  int         *f_ts;               /* last field output time step */

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

  fvm_to_melissa_writer_t  *writer;      /* pointer to writer structure */

  const char               *name;        /* current field name */
  int                       time_step;   /* current_time_step */
  bool                      call_init;   /* call melissa_init ? */

} _melissa_context_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Since the Melissa API specifies the field name (and communicator)
   only, only one Melissa writer can be active a a given time. */

static fvm_to_melissa_writer_t  *_writer = NULL;

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
_field_c_output(void           *context,
                cs_datatype_t   datatype,
                int             dimension,
                int             component_id,
                cs_gnum_t       block_start,
                cs_gnum_t       block_end,
                void           *buffer)
{
  _melissa_context_t *c = context;
  fvm_to_melissa_writer_t  *w = c->writer;

  cs_gnum_t n_values = block_end - block_start;

  assert(datatype == CS_DOUBLE);
  double *values = (double*)buffer;

  /* Build component name */

  char tmpn[128], tmpe[6];
  char *c_name = tmpn;

  fvm_writer_field_component_name(tmpe, 6, false, dimension, component_id);

  size_t lce = strlen(tmpe);
  size_t l =  strlen(c->name) + 1;

  if (lce > 0)
    l += 2 + lce;

  if (l > 128)
    BFT_MALLOC(c_name, l, char);

  if (lce > 0)
    sprintf(c_name, "%s[%s]", c->name, tmpe);
  else
    strcpy(c_name, c->name);

  /* Check which API we use */

#if defined(HAVE_MELISSA_NO_MPI)
  bool use_melissa_mpi = false;
#else
  bool use_melissa_mpi = true;
#endif

#if defined(HAVE_MELISSA_MPI)
  if (w->block_comm != MPI_COMM_NULL)
    use_melissa_mpi = true;
#endif

  /* Case of first call for this variable */

  if (c->call_init) {

    if (w->tracefile != NULL)
      fprintf(w->tracefile, "[melissa_init] name: %s (time_step: %d)\n",
              c_name, c->time_step);

    if (w->dry_run == false) {
#if defined(HAVE_MELISSA_MPI)
      if (use_melissa_mpi == true)
        melissa_init(c_name, n_values, w->block_comm);
#endif
#if defined(HAVE_MELISSA_NO_MPI)
      if (use_melissa_mpi == false)
        melissa_init_no_mpi(c_name, n_values);
#endif
    }

  }

  /* Output values */

  if (w->tracefile != NULL)
    fprintf(w->tracefile, "[melissa_send] name: %s (time_step: %d)\n",
            c_name, c->time_step);

  if (w->dry_run == false) {
#if defined(HAVE_MELISSA_MPI)
    if (use_melissa_mpi == true)
      melissa_send(c_name, values);
#endif
#if defined(HAVE_MELISSA_NO_MPI)
    if (use_melissa_mpi == false)
      melissa_send_no_mpi(c_name, values);
#endif
  }

  /* Local cleanup */
  if (c_name != tmpn)
    BFT_FREE(c_name);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to Melissa output writer.
 *
 * Options are:
 *   dry_run               trace output to <name>.log file, but do not
 *                         actually communicate with Melissa server
 *   rank_step=<integer>   MPI rank step
 *   trace                 trace output to <name>.log file
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
  CS_UNUSED(time_dependency);

  fvm_to_melissa_writer_t  *w = NULL;

  /* Parse options */

  int rank_step = 1;
  bool trace = false;
  bool dry_run = false;

  if (options != NULL) {

    int i1 = 0, i2 = 0;
    int l_tot = strlen(options);

    const char rs[] = "rank_step=";
    const int l_rs = strlen(rs);

    while (i1 < l_tot) {

      for (i2 = i1; i2 < l_tot && options[i2] != ' '; i2++);
      int l_opt = i2 - i1;

      if ((l_opt == 5) && (strncmp(options + i1, "trace", l_opt) == 0))
        trace = true;

      else if ((l_opt == 7) && (strncmp(options + i1, "dry_run", l_opt) == 0)) {
        dry_run = true;
        trace = true;
      }

      else if ((strncmp(options + i1, rs, l_rs) == 0)) {
        if (l_opt < l_rs+32) { /* 32 integers more than enough
                                  for maximum integer string */
          char options_c[32];
          strncpy(options_c, options+i1+l_rs, l_opt-l_rs);
          options_c[l_opt-l_rs] = '\0';
          rank_step = atoi(options_c);

        }
      }

      for (i1 = i2 + 1; i1 < l_tot && options[i1] == ' '; i1++);

    }
  }

  if (dry_run == false && _writer != NULL) {
    bft_error(__FILE__, __LINE__, errno,
              _("Error creating Melissa writer: \"%s\":\n"
                "only one Melissa server may be used, and is already used by\n"
                "writer: \"%s\"."), name, _writer->name);
    return NULL;
  }

  /* Initialize writer */

  BFT_MALLOC(w, 1, fvm_to_melissa_writer_t);

  BFT_MALLOC(w->name, strlen(name) + 1, char);
  strcpy(w->name, name);

  w->dry_run = dry_run;
  w->tracefile = NULL;

  w->rank = 0;
  w->n_ranks = 1;

  w->buffer_size = 0;

  w->f_map = cs_map_name_to_id_create();
  w->f_ts = NULL;

#if defined(HAVE_MPI)
  {
    int mpi_flag, rank, n_ranks;
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
      if (rank_step < 1)
        rank_step = 1;
      else if (rank_step > n_ranks)
        rank_step = n_ranks;

#if defined(HAVE_MELISSA_MPI)
      w->min_rank_step = rank_step;
      if (rank_step > 1) {
        w->block_comm = cs_base_get_rank_step_comm(rank_step);
      }
      else
        w->block_comm = comm;
#else
      w->min_rank_step = n_ranks;
      w->block_comm = MPI_COMM_NULL;
#endif

      w->comm = comm;
    }
  }
#endif /* defined(HAVE_MPI) */

  if (trace && w->rank < 1) {

    int  path_len = 0;
    if (path != NULL)
      path_len = strlen(path) + 1;

    int tracefile_path_l = path_len + strlen(name) + strlen(".log") + 1;
    char *tracefile_path;
    BFT_MALLOC(tracefile_path, tracefile_path_l, char);
    if (path != NULL)
      strcpy(tracefile_path, path);
    else
      tracefile_path[0] = '\0';
    strcat(tracefile_path, name);
    strcat(tracefile_path, ".log");

    w->tracefile = fopen(tracefile_path, "w");
    if (w->tracefile ==  NULL)
      bft_error(__FILE__, __LINE__, errno,
                _("Error opening file: \"%s\""), tracefile_path);

    BFT_FREE(tracefile_path);

  }

  if (w->dry_run == false)
    _writer = w;

  /* Return writer */
  return w;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to Melissa output writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque Melissa writer structure.
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_melissa_finalize_writer(void  *this_writer_p)
{
  fvm_to_melissa_writer_t *w = (fvm_to_melissa_writer_t *)this_writer_p;

  if (_writer == w)
    _writer = NULL;

  if (w->tracefile != NULL) {
    if (fclose(w->tracefile) != 0)
      bft_error(__FILE__, __LINE__, errno,
                _("Error closing file: \"%s.log\""), w->name);
  }
  BFT_FREE(w->name);

  cs_map_name_to_id_destroy(&(w->f_map));
  BFT_FREE(w->f_ts);

  if (w->dry_run == false)
    melissa_finalize();

  BFT_FREE(w);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with a Melissa geometry.
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
  CS_UNUSED(this_writer_p);
  CS_UNUSED(time_value);
  CS_UNUSED(time_step);
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a Melissa output
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *----------------------------------------------------------------------------*/

void
fvm_to_melissa_export_nodal(void               *this_writer_p,
                            const fvm_nodal_t  *mesh)
{
  CS_UNUSED(this_writer_p);
  CS_UNUSED(mesh);
}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a Melissa output.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field.
 *
 * Fields already exported by another Melissa writer will be ignored,
 * as will fields already exported for the same or greater time step,
 * or exported with a different time dependency.
 * Use a different field name if you really need to output fields in such
 * a configuration.
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
 *   time_value       <-- associated time value (ignored)
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
  CS_UNUSED(time_value);

  fvm_to_melissa_writer_t  *w = (fvm_to_melissa_writer_t *)this_writer_p;

  _melissa_context_t  c;
  c.writer = w;
  c.name = name;
  c.time_step = time_step;
  c.call_init = false;

  fvm_writer_field_helper_t   *helper = NULL;
  fvm_writer_section_t        *export_list = NULL;

  const int  n_ranks = w->n_ranks;

  int f_id = cs_map_name_to_id_try(w->f_map, name);
  if (f_id < 0) {
    f_id = cs_map_name_to_id(w->f_map, name);
    int n_fields = cs_map_name_to_id_size(w->f_map);
    BFT_REALLOC(w->f_ts, n_fields, int);
    c.call_init = true;
    if (w->tracefile != NULL) {
      fprintf
        (w->tracefile,
         "[init] name: %s (time step: %d; location: %d; dimension: %d)\n",
         name, time_step, location, dimension);
    }
  }
  else {
    int prev_ts = w->f_ts[f_id];
    if (prev_ts < 0 || prev_ts >= time_step) {
      if (w->tracefile != NULL) {
        fprintf
          (w->tracefile,
           "[ignore] name: %s (time step: %d; location: %d; dimension: %d)\n",
           name, time_step, location, dimension);
      }
      return;
    }
  }

  w->f_ts[f_id] = time_step;

  /* Initialize writer helper */
  /*--------------------------*/

  int export_dim = fvm_nodal_get_max_entity_dim(mesh);

  /* Build list of sections that are used here, in order of output */
  export_list = fvm_writer_export_list(mesh,
                                       export_dim,
                                       export_dim,
                                       -1,
                                       false, /* group by type */
                                       true,  /* group all */
                                       false,
                                       false,
                                       false,
                                       false);

  helper = fvm_writer_field_helper_create(mesh,
                                          export_list,
                                          dimension,
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

  /* Per node variable */
  /*-------------------*/

  if (location == FVM_WRITER_PER_NODE) {

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
                                     _field_c_output);

  }

  /* Per element variable */
  /*----------------------*/

  else if (location == FVM_WRITER_PER_ELEMENT) {

    /* Output per grouped sections */

    fvm_writer_field_helper_output_e(helper,
                                     &c,
                                     export_list,
                                     dimension,
                                     interlace,
                                     NULL,
                                     n_parent_lists,
                                     parent_num_shift,
                                     datatype,
                                     field_values,
                                     _field_c_output);

  } /* End for per element variable */

  /* Free helper structures */

  fvm_writer_field_helper_destroy(&helper);
  BFT_FREE(export_list);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
