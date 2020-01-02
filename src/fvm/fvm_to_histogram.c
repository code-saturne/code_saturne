/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to histogram files
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include <errno.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "fvm_defs.h"

#include "fvm_convert_array.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_to_vtk_histogram.h"
#include "fvm_writer_helper.h"
#include "fvm_writer_priv.h"

#include "cs_block_dist.h"
#include "cs_file.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_histogram.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Context structure for fvm_writer_field_helper_output_* functions.
 *----------------------------------------------------------------------------*/

typedef struct {

  fvm_to_histogram_writer_t  *writer;    /* Pointer to writer structure */

  const char                 *name;      /* Field name */

} _histogram_context_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static void * _catalyst_plugin = NULL;

static fvm_to_histogram_display_t  *_fvm_to_vtk_display_histogram_png = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_PLUGIN_CATALYST)

/*----------------------------------------------------------------------------
 * Load Catalyst plugin.
 *----------------------------------------------------------------------------*/

static void
_load_plugin_catalyst(void)
{
  /* Open from shared library */

  _catalyst_plugin = cs_base_dlopen_plugin("fvm_catalyst");

  /* Load symbols from shared library */

  /* Function pointers need to be double-casted so as to first convert
     a (void *) type to a memory address and then convert it back to the
     original type. Otherwise, the compiler may issue a warning.
     This is a valid ISO C construction. */

  _fvm_to_vtk_display_histogram_png = (fvm_to_histogram_display_t *) (intptr_t)
    cs_base_get_dl_function_pointer(_catalyst_plugin,
                                    "fvm_to_vtk_display_histogram_png",
                                    true);
}

#endif /* defined(HAVE_PLUGIN_CATALYST)*/

/*----------------------------------------------------------------------------*/
/*!
 * Display histograms in txt format.
 *
 * parameters:
 *  var_min  <--  minimum variable value
 *  var_max  <--  maximum variable value
 *  count    <--  count for each histogram slice
 *  w        <--> histogram writer
 *  var_name <--  name of the variable
 */
/*----------------------------------------------------------------------------*/

static void
_display_histogram_txt(cs_real_t                   var_min,
                       cs_real_t                   var_max,
                       cs_gnum_t                   count[],
                       fvm_to_histogram_writer_t  *w,
                       char                       *var_name)
{
  double var_step;
  int i, j;

  /* Open the txt file */
  w->f = fopen(w->file_name, "w");

  if (w->f ==  NULL) {
    bft_error(__FILE__, __LINE__, errno,
              _("Error opening file: \"%s\""), w->file_name);
    return;
  }

  /* Print header */
  fprintf(w->f, "# Code_Saturne histograms\n#\n");
  if (w->nt < 0)
    fprintf(w->f, "# time independent\n");
  else {
    fprintf(w->f, "# time step id: %i\n", w->nt);
    fprintf(w->f, "# time:         %12.5e\n#\n", w->t);
  }

  /* Print base min, max, and increment */
  fprintf(w->f, "# Variable : %s\n\n",var_name);

  fprintf(w->f, _("    minimum value =         %10.5e\n"), (double)var_min);
  fprintf(w->f, _("    maximum value =         %10.5e\n\n"), (double)var_max);

  var_step = CS_ABS(var_max - var_min) / w->n_sub;

  if (CS_ABS(var_max - var_min) > 0.) {

    /* Number of elements in each subdivision */

    for (i = 0, j = 1; i < w->n_sub - 1; i++, j++)
      fprintf(w->f, "    %3d : [ %10.5e ; %10.5e [ = %10llu\n",
                 i+1, var_min + i*var_step, var_min + j*var_step,
                 (unsigned long long)(count[i]));

    fprintf(w->f, "    %3d : [ %10.5e ; %10.5e ] = %10llu\n",
               w->n_sub,
               var_min + (w->n_sub - 1)*var_step,
               var_max,
               (unsigned long long)(count[w->n_sub - 1]));

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * Display histograms in tex format.
 *
 * parameters:
 *  var_min  <--  minimum variable value
 *  var_max  <--  maximum variable value
 *  count    <--  count for each histogram slice
 *  w        <--> histogram writer
 *  var_name <--  name of the variable
 */
/*----------------------------------------------------------------------------*/

static void
_display_histogram_tex(cs_real_t                   var_min,
                       cs_real_t                   var_max,
                       cs_gnum_t                   count[],
                       fvm_to_histogram_writer_t  *w,
                       char                       *var_name)
{
  double var_step = CS_ABS(var_max - var_min) / w->n_sub;
  int i;

  if (var_step > 0) {
    /* Open the tex file if non-zero histogram */
    w->f = fopen(w->file_name, "w");

    if (w->f ==  NULL) {
      bft_error(__FILE__, __LINE__, errno,
                _("Error opening file: \"%s\""), w->file_name);
      return;
    }

    fprintf(w->f, "\\begin{center}\n");
    fprintf(w->f, "\\begin{tikzpicture}[font=\\footnotesize]\n");

    fprintf(w->f, "  \\begin{axis}[\n");
    fprintf(w->f, "    ybar,\n");
    fprintf(w->f, "    bar width=18pt,\n");
    fprintf(w->f, "    xlabel={%s},\n",var_name);
    fprintf(w->f, "    ylabel={Number of elements},\n");
    fprintf(w->f, "    ymin=0,\n");
    fprintf(w->f, "    ytick=\\empty,\n");
    fprintf(w->f, "    xtick=data,\n");
    fprintf(w->f, "    axis x line=bottom,\n");
    fprintf(w->f, "    axis y line=left,\n");
    fprintf(w->f, "    enlarge x limits=0.06,\n");

    fprintf(w->f, "    symbolic x coords={");
    for (i = 0 ; i < w->n_sub-1 ; i++)
      fprintf(w->f, "%.3e,",var_min + (i+0.5)*var_step);
    fprintf(w->f, "%.3e},\n",var_min + (w->n_sub-0.5)*var_step);

    fprintf(w->f, "    xticklabel style={rotate=45,font=\\scriptsize},\n");
    fprintf(w->f,
            "    nodes near coords={\\pgfmathprintnumber\\pgfplotspointmeta}\n");
    fprintf(w->f, "  ]\n");
    fprintf(w->f, "    \\addplot[fill=blue] coordinates {\n");
    for (i = 0 ; i < w->n_sub ; i++)
      fprintf(w->f, "        (%.3e,%llu)\n",
              var_min + (i+0.5)*var_step, (unsigned long long)count[i]);
    fprintf(w->f, "    };\n");
    fprintf(w->f, "  \\end{axis}\n");

    fprintf(w->f, "\\end{tikzpicture}\n");
    fprintf(w->f, "\\end{center}\n");
  }
}

/*----------------------------------------------------------------------------
 * Compute the minimum and the maximum of a vector (locally).
 *
 * parameters:
 *   n_vals    <-- local number of elements
 *   var       <-- pointer to vector
 *   min       --> minimum
 *   max       --> maximum
 *----------------------------------------------------------------------------*/

static void
_compute_local_minmax(cs_lnum_t        n_vals,
                      const cs_real_t  var[],
                      cs_real_t       *min,
                      cs_real_t       *max)
{
  cs_lnum_t  i;
  cs_real_t  _min = DBL_MAX, _max = -DBL_MAX;

  for (i = 0; i < n_vals; i++) {
    _min = CS_MIN(_min, var[i]);
    _max = CS_MAX(_max, var[i]);
  }

  *min = _min;
  *max = _max;
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a real vector.
 *
 * parameters:
 *   var_min      <--  minimum variable value
 *   var_max      <--  maximum variable value
 *   count        <->  count for each histogram slice (size: w->n_sub)
 *                     local values in, global values out
 *   display_func <--  function pointer to display the histogram
 *   w            <--> histogram writer
 *   var_name     <--  name of the variable
 *----------------------------------------------------------------------------*/

static void
_display_histograms(cs_real_t                   var_min,
                    cs_real_t                   var_max,
                    cs_gnum_t                   count[],
                    fvm_to_histogram_display_t  *display_func,
                    fvm_to_histogram_writer_t   *w,
                    char                        *var_name)
{
  int  i;

#if defined(HAVE_MPI)

  if (w->n_ranks > 1) {

    cs_gnum_t *g_count = NULL;

    BFT_MALLOC(g_count, w->n_sub, cs_gnum_t);

    MPI_Allreduce(count, g_count, w->n_sub, CS_MPI_GNUM, MPI_SUM,
                  w->comm);

    for (i = 0; i < w->n_sub; i++)
      count[i] = g_count[i];

    BFT_FREE(g_count);
  }

#endif

  if (w->rank == 0)
    display_func(var_min, var_max, count, w, var_name);
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a real vector on cells or
 * boundary faces.
 *
 * parameters:
 *   n_vals       <--  number of values
 *   var          <--  variable values
 *   display_func <--  function pointer to display the histogram
 *   w            <--> histogram writer
 *   var_name     <--  name of the variable
 *----------------------------------------------------------------------------*/

static void
_histogram(cs_lnum_t                    n_vals,
           const cs_real_t              var[],
           fvm_to_histogram_display_t  *display_func,
           fvm_to_histogram_writer_t   *w,
           char                        *var_name)
{
  cs_lnum_t  i;
  int        j, k;

  cs_real_t  step;
  cs_real_t  max, min, _max, _min;

  cs_gnum_t *count = NULL;

  BFT_MALLOC(count, w->n_sub, cs_gnum_t);

  assert (sizeof(double) == sizeof(cs_real_t));

  /* Compute global min and max */

  _compute_local_minmax(n_vals, var, &_min, &_max);

  /* Default initialization */

  min = _min;
  max = _max;

#if defined(HAVE_MPI)

  if (w->n_ranks > 1) {
    MPI_Allreduce(&_min, &min, 1, CS_MPI_REAL, MPI_MIN,
                  w->comm);

    MPI_Allreduce(&_max, &max, 1, CS_MPI_REAL, MPI_MAX,
                  w->comm);
  }

#endif

  /* Define axis subdivisions */

  for (j = 0; j < w->n_sub; j++)
    count[j] = 0;

  if (CS_ABS(max - min) > 0.) {

    step = CS_ABS(max - min) / w->n_sub;

    /* Loop on values */

    for (i = 0; i < n_vals; i++) {

      /* Associated subdivision */

      for (j = 0, k = 1; k < w->n_sub; j++, k++) {
        if (var[i] < min + k*step)
          break;
      }
      count[j] += 1;

    }

  }

  _display_histograms(min, max, count, display_func, w, var_name);

  BFT_FREE(count);
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

  _histogram_context_t *c = (_histogram_context_t *)context;

  fvm_to_histogram_writer_t  *w = c->writer;

  const int n_vals = (block_end > block_start) ? block_end - block_start : 0;

  char tmpn[128], tmpe[6];

  char *var_name = tmpn;

  fvm_writer_field_component_name(tmpe, 6, false, dimension, component_id);

  size_t lce = strlen(tmpe);
  size_t lv =  strlen(c->name) + 1;

  if (lce > 0)
    lv += 2 + lce;

  if (lv > 128)
    BFT_MALLOC(var_name, lv, char);

  if (lce > 0)
    sprintf(var_name, "%s[%s]", c->name, tmpe);
  else
    strcpy(var_name, c->name);

  /* File name */
  char t_stamp[32];
  sprintf(t_stamp, "_%s_%.4i", var_name, w->nt);
  size_t l =   strlen(w->path) + strlen(w->name)
             + strlen(t_stamp) + 4 + 1;

  BFT_REALLOC(w->file_name, l, char);

  const cs_real_t *vals = (const cs_real_t *)buffer;

  if (w->format == CS_HISTOGRAM_TXT) {

    sprintf(w->file_name, "%s%s%s.txt", w->path, w->name, t_stamp);
    _histogram(n_vals, vals, _display_histogram_txt, w, var_name);

  } else if (w->format == CS_HISTOGRAM_TEX) {

    sprintf(w->file_name, "%s%s%s.tex", w->path, w->name, t_stamp);
    _histogram(n_vals, vals, _display_histogram_tex, w, var_name);

#if defined(HAVE_CATALYST)
  } else if (w->format == CS_HISTOGRAM_PNG) {

    sprintf(w->file_name, "%s%s%s.png", w->path, w->name, t_stamp);
    _histogram(n_vals, vals, _fvm_to_vtk_display_histogram_png, w, var_name);

#endif
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to histogram file writer.
 *
 * Options are:
 *   txt                 output txt (space-separated) files
 *   tex                 output TeX (TixZ) files
 *   png                 output PNG files
 *   [n_sub]             number of subdivisions
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque histogram writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_histogram_init_writer(const char             *name,
                             const char             *path,
                             const char             *options,
                             fvm_writer_time_dep_t   time_dependency,
                             MPI_Comm                comm)
#else
void *
fvm_to_histogram_init_writer(const char             *name,
                             const char             *path,
                             const char             *options,
                             fvm_writer_time_dep_t   time_dependency)
#endif
{
  CS_UNUSED(time_dependency);

  fvm_to_histogram_writer_t  *w = NULL;

  /* Initialize writer */

  BFT_MALLOC(w, 1, fvm_to_histogram_writer_t);

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

  w->format = CS_HISTOGRAM_TXT;

  w->nt = -1;
  w->t = -1;
  w->buffer = NULL;

  w->file_name = NULL;
  w->f = NULL;

  w->n_sub = 5; /* default */

  /* Parse options */

  if (options != NULL) {

    int i1, i2, l_opt;
    int l_tot = strlen(options);
    bool n_sub_read = false;

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {

      for (i2 = i1; i2 < l_tot && options[i2] != ' '; i2++);
      l_opt = i2 - i1;

      if ((l_opt == 3) && (strncmp(options + i1, "txt", l_opt) == 0))
        w->format = CS_HISTOGRAM_TXT;
      else if ((l_opt == 3) && (strncmp(options + i1, "tex", l_opt) == 0)) {
        w->format = CS_HISTOGRAM_TEX;
        if (!n_sub_read)
          w->n_sub = 10;
      }
#if defined(HAVE_CATALYST)
      else if ((l_opt == 3) && (strncmp(options + i1, "png", l_opt) == 0)) {
        w->format = CS_HISTOGRAM_PNG;
        if (!n_sub_read)
          w->n_sub = 10;
#if defined(HAVE_PLUGIN_CATALYST)
        _load_plugin_catalyst();
#else
        _fvm_to_vtk_display_histogram_png = fvm_to_vtk_display_histogram_png;
#endif
      }
#endif
      else {
        const char *n_sub_opt = options+i1;
        while (*n_sub_opt != '\0' && !isdigit(*n_sub_opt))
          n_sub_opt++;
        if (atoi(n_sub_opt) > 0) {
          w->n_sub = atoi(n_sub_opt);
          n_sub_read = true;
        }
      }

      for (i1 = i2 + 1 ; i1 < l_tot && options[i1] == ' ' ; i1++);

    }

  }

  /* Return writer */

  return w;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to histogram file writer.
 *
 * parameters:
 *   writer <-- pointer to opaque histogram Gold writer structure.
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_histogram_finalize_writer(void  *writer)
{
  fvm_to_histogram_writer_t  *w
    = (fvm_to_histogram_writer_t *)writer;

  BFT_FREE(w->name);
  BFT_FREE(w->path);

  fvm_to_histogram_flush(writer);

  BFT_FREE(w->file_name);

#if defined(HAVE_PLUGIN_CATALYST)
  if (w->format == CS_HISTOGRAM_PNG)
    cs_base_dlclose("fvm_catalyst",
                    _catalyst_plugin); /* decrease reference count or free */
#endif

  BFT_FREE(w);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with a histogram geometry.
 *
 * parameters:
 *   writer     <-- pointer to associated writer
 *   time_step  <-- time step number
 *   time_value <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_histogram_set_mesh_time(void    *writer,
                               int      time_step,
                               double   time_value)
{
  fvm_to_histogram_writer_t  *w = (fvm_to_histogram_writer_t *)writer;

  w->nt = time_step;
  w->t = time_value;

  fvm_to_histogram_flush(writer);
}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a histogram file.
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
fvm_to_histogram_export_field(void                  *writer,
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
  fvm_to_histogram_writer_t  *w = (fvm_to_histogram_writer_t *)writer;

  /* If time step changes, update it */

  if (time_step != w->nt)
    fvm_to_histogram_set_mesh_time(writer,
                                   time_step,
                                   time_value);

  /* Initialize writer helper */

  cs_datatype_t  dest_datatype = CS_REAL_TYPE;

  if (datatype >= CS_INT32 && datatype <= CS_UINT64)
    dest_datatype = CS_INT64;

  fvm_writer_section_t  *export_list
    = fvm_writer_export_list(mesh,
                             fvm_nodal_get_max_entity_dim(mesh),
                             true,
                             true,
                             false,
                             false,
                             false,
                             true);

  fvm_writer_field_helper_t  *helper
    = fvm_writer_field_helper_create(mesh,
                                     export_list,
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

  _histogram_context_t c = {.writer = w, .name = name};

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
                                   _field_output);

  BFT_FREE(export_list);

  /* Free helper structures */

  fvm_writer_field_helper_destroy(&helper);
}

/*----------------------------------------------------------------------------
 * Flush files associated with a given writer.
 *
 * parameters:
 *   writer <-- pointer to associated writer
 *----------------------------------------------------------------------------*/

void
fvm_to_histogram_flush(void  *writer)
{
  fvm_to_histogram_writer_t  *w = (fvm_to_histogram_writer_t *)writer;

  if (w->f != NULL && w->buffer != NULL) {

    if (fclose(w->f) != 0)
      bft_error(__FILE__, __LINE__, errno,
                _("Error closing file: \"%s\""), w->file_name);

    w->f = NULL;

  }

  BFT_FREE(w->buffer);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
