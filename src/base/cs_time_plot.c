/*============================================================================
 * Writer helper for time-varying 1d data
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_base.h"
#include "cs_timer.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_time_plot.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Basic file output structure */
/*-----------------------------*/

/* This structure should be non-null only for ranks
   on which it is active */

struct _cs_time_plot_t {

  char       *plot_name;        /* Associated plot name */
  char       *file_name;        /* Associated file name */

  FILE       *f;                /* Associated file */

  cs_time_plot_format_t   format;  /* Associated format */

  bool        use_iteration;    /* Use iteration number instead of
                                 * physical time ? */
  bool        write_time_value; /* Output time value ? */

  double      flush_times[2];   /* 0: File flush time interval
                                 *    (if < 0, no forced flush)
                                 * 1: last write time
                                 *    (-2 before first write) */
  double      buffer_steps[2];  /* 0: Maximum number of time steps in
                                 *    output buffer (if > 0, file is
                                 *    not kept open)
                                 * 1: current number of time steps in
                                 *    output buffer */

  size_t      buffer_size;      /* Buffer size if required */
  size_t      buffer_end;       /* Current buffer end */
  char       *buffer;           /* Associated buffer if required */

  struct _cs_time_plot_t  *prev;  /* Previous in flush list */
  struct _cs_time_plot_t  *next;  /* Next in flush list */

};

/*============================================================================
 *  Global variables
 *============================================================================*/

static cs_time_plot_t   *_plots_head = NULL, *_plots_tail = NULL;

static float             _flush_wtime_default = -1;
static int               _n_buffer_steps_default = -1;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Ensure buffer is large enough for upcoming data.
 *
 * The buffer is reallocated if necessary.
 *
 * parameters:
 *   p        <-> time plot values file handler
 *   min_size <-- minimum buffer size
 *----------------------------------------------------------------------------*/

static void
_ensure_buffer_size(cs_time_plot_t  *p,
                    size_t           min_size)
{
  /* Write data to line buffer */

  if (p->buffer_size < min_size) {
    p->buffer_size = CS_MAX(1, p->buffer_size);
    while (p->buffer_size < min_size)
      p->buffer_size *= 2;
    BFT_REALLOC(p->buffer, p->buffer_size, char);
  }
}

/*----------------------------------------------------------------------------
 * Write file header for xmgrace/qsplotlib readable .dat files
 *
 * parameters:
 *   p                <-> time plot values file handler
 *   n_probes         <-- number of probes associated with this variable ?
 *   probe_list       <-- numbers (1 to n) of probes if filtered, or NULL
 *   probe_coords     <-- probe coordinates
 *   probe_names      <-- probe names, or NULL
 *----------------------------------------------------------------------------*/

static void
_write_probe_header_dat(cs_time_plot_t    *p,
                        int                n_probes,
                        const int         *probe_list,
                        const cs_real_t    probe_coords[],
                        const char        *probe_names[])
{
  int i, probe_id;
  int col_id = 0;
  FILE *_f = p->f;

  if (_f != NULL) {
    fclose(_f);
    p->f = NULL;
  }

  _f = fopen(p->file_name, "w");
  if (_f == NULL) {
    bft_error(__FILE__, __LINE__, errno,
              _("Error opening file: \"%s\""), p->file_name);
    return;
  }

  fprintf(_f, _("# Time varying values for: %s\n"
                "#\n"), p->plot_name);

  if (probe_coords != NULL) {
    fprintf(_f, _("# Monitoring point coordinates:\n"));
    for (i = 0; i < n_probes; i++) {
      probe_id = i;
      if (probe_list != NULL)
        probe_id = probe_list[i] - 1;
      if (probe_names != NULL)
        fprintf(_f, "# %16s [%14.7e, %14.7e, %14.7e]\n",
                probe_names[i],
                probe_coords[probe_id*3],
                probe_coords[probe_id*3 + 1],
                probe_coords[probe_id*3 + 2]);
      else
        fprintf(_f, "#   %6i [%14.7e, %14.7e, %14.7e]\n",
                probe_id + 1,
                probe_coords[probe_id*3],
                probe_coords[probe_id*3 + 1],
                probe_coords[probe_id*3 + 2]);
    }
    fprintf(_f, "#\n");
  }
  else if (probe_names != NULL) {
    fprintf(_f, _("# Monitoring points:\n"));
    for (i = 0; i < n_probes; i++)
      fprintf(_f, "# %s\n", probe_names[i]);
    fprintf(_f, "#\n");
  }

  fprintf(_f, _("# Columns:\n"));
  if (p->use_iteration)
    fprintf(_f, _("#   %d:     Time step number\n"), col_id++);
  else
    fprintf(_f, _("#   %d:     Physical time\n"), col_id++);
  fprintf(_f, _("#   %d - :  Values at monitoring points\n"), col_id);

  fprintf(_f,
          "#\n"
          "#TITLE: %s\n"
          "#COLUMN_TITLES: ", p->plot_name);
  if (p->use_iteration)
    fprintf(_f, " nt");
  else
    fprintf(_f, " t");
  for (i = 0; i < n_probes; i++) {
    probe_id = i;
    if (probe_list != NULL)
      probe_id = probe_list[i] - 1;
    if (probe_names != NULL) {
      if (probe_coords != NULL)
        fprintf(_f, " | %s [%9.5e, %9.5e, %9.5e]",
                probe_names[i],
                probe_coords[probe_id*3],
                probe_coords[probe_id*3 + 1],
                probe_coords[probe_id*3 + 2]);
      else
        fprintf(_f, " | %s", probe_names[i]);
    }
    else {
      if (probe_coords != NULL)
        fprintf(_f, " | %d [%9.5e, %9.5e, %9.5e]",
                probe_id + 1,
                probe_coords[probe_id*3],
                probe_coords[probe_id*3 + 1],
                probe_coords[probe_id*3 + 2]);
      else
        fprintf(_f, " | %d", probe_id + 1);
    }
  }
  fprintf(_f, "\n");

  fprintf(_f, "#COLUMN_UNITS: ");
  if (p->use_iteration)
    fprintf(_f, " iter");
  else
    fprintf(_f, " s");
  for (probe_id = 0; probe_id < n_probes; probe_id++)
    fprintf(_f, " -");
  fprintf(_f, "\n#\n");

  /* Close file or assign it to handler depending on options */

  if (p->buffer_steps[0] > 0) {
    if (fclose(_f) != 0)
      bft_error(__FILE__, __LINE__, errno,
                _("Error closing file: \"%s\""), p->file_name);
  }
  else
    p->f = _f;
}

/*----------------------------------------------------------------------------
 * Write probe coordinates for CSV files
 *
 * parameters:
 *   file_prefix      <-- file name prefix
 *   plot_name        <-- plot name
 *   n_probes         <-- number of probes associated with this plot
 *   probe_list       <-- numbers (1 to n) of probes if filtered, or NULL
 *   probe_coords     <-- probe coordinates
 *----------------------------------------------------------------------------*/

static void
_write_probe_coords_csv(const char        *file_prefix,
                        const char        *plot_name,
                        int                n_probes,
                        const int         *probe_list,
                        const cs_real_t    probe_coords[])
{
  int i, probe_id;
  char *file_name;
  FILE *_f;

  BFT_MALLOC(file_name,
             strlen(file_prefix) + strlen(plot_name) + strlen("_coords") + 4 + 1,
             char);

  if (probe_coords != NULL) {
    sprintf(file_name, "%s%s%s.csv", file_prefix, plot_name, "_coords");

    _f = fopen(file_name, "w");
    if (_f == NULL) {
      bft_error(__FILE__, __LINE__, errno,
          _("Error opening file: \"%s\""), file_name);
      return;
    }

    fprintf(_f, "x, y, z\n");
    for (i = 0; i < n_probes; i++) {
      probe_id = i;
      if (probe_list != NULL)
        probe_id = probe_list[i] - 1;
      fprintf(_f, "%14.7e, %14.7e, %14.7e\n",
          probe_coords[probe_id*3],
          probe_coords[probe_id*3 + 1],
          probe_coords[probe_id*3 + 2]);
    }

    /* Close file or assign it ot handler depending on options */

    if (fclose(_f) != 0)
      bft_error(__FILE__, __LINE__, errno,
          _("Error closing file: \"%s\""), file_name);

  }

  BFT_FREE(file_name);
}

/*----------------------------------------------------------------------------
 * Write file header for CSV files
 *
 * parameters:
 *   p                <-> time plot values file handler
 *   n_probes         <-- number of probes associated with this variable ?
 *   probe_list       <-- numbers (1 to n) of probes if filtered, or NULL
 *   probe_coords     <-- probe coordinates
 *   probe_names      <-- probe names, or NULL
 *----------------------------------------------------------------------------*/

static void
_write_probe_header_csv(cs_time_plot_t    *p,
                        int                n_probes,
                        const int         *probe_list,
                        const cs_real_t    probe_coords[],
                        const char        *probe_names[])
{
  int i, probe_id;
  FILE *_f = p->f;

  if (_f != NULL) {
    fclose(_f);
    p->f = NULL;
  }

  _f = fopen(p->file_name, "w");
  if (_f == NULL) {
    bft_error(__FILE__, __LINE__, errno,
              _("Error opening file: \"%s\""), p->file_name);
    return;
  }

  if (p->use_iteration)
    fprintf(_f, " iteration");
  else
    fprintf(_f, "t");
  for (i = 0; i < n_probes; i++) {
    probe_id = i;
    if (probe_list != NULL)
      probe_id = probe_list[i] - 1;
    if (probe_coords != NULL) {
      if (probe_names != NULL)
        fprintf(_f, ", %s [%9.5e; %9.5e; %9.5e]",
                probe_names[i],
                probe_coords[probe_id*3],
                probe_coords[probe_id*3 + 1],
                probe_coords[probe_id*3 + 2]);
      else
        fprintf(_f, ", %d [%9.5e; %9.5e; %9.5e]",
                probe_id + 1,
                probe_coords[probe_id*3],
                probe_coords[probe_id*3 + 1],
                probe_coords[probe_id*3 + 2]);
    }
    else if (probe_names != NULL)
      fprintf(_f, ", %s", probe_names[i]);
    else
      fprintf(_f, ", %d", probe_id + 1);
  }
  fprintf(_f, "\n");

  /* Close file or assign it to handler depending on options */

  if (p->buffer_steps[0] > 0) {
    if (fclose(_f) != 0)
      bft_error(__FILE__, __LINE__, errno,
                _("Error closing file: \"%s\""), p->file_name);
  }
  else
    p->f = _f;
}

/*----------------------------------------------------------------------------
 * Write file header for xmgrace/qsplotlib readable .dat files
 *
 * parameters:
 *   p                  <-> time plot values file handler
 *   n_structures       <-- number of structures associated with this plot
 *   mass_matrixes      <-- mass matrix coefficients (3x3 blocks)
 *   damping_matrixes   <-- damping matrix coefficients (3x3 blocks)
 *   stiffness_matrixes <-- stiffness matrix coefficients (3x3 blocks)
 *----------------------------------------------------------------------------*/

static void
_write_struct_header_dat(cs_time_plot_t    *p,
                         int                n_structures,
                         const cs_real_t    mass_matrixes[],
                         const cs_real_t    damping_matrixes[],
                         const cs_real_t    stiffness_matrixes[])
{
  int i, struct_id;
  int col_id = 0;
  FILE *_f = p->f;

  const int perm_id[9] = {0, 3, 6, 1, 4, 7, 2, 5, 8};

  if (_f != NULL) {
    fclose(_f);
    p->f = NULL;
  }

  _f = fopen(p->file_name, "w");
  if (_f == NULL) {
    bft_error(__FILE__, __LINE__, errno,
              _("Error opening file: \"%s\""), p->file_name);
    return;
  }

  fprintf(_f, _("# Time varying values for: %s\n"
                "#\n"), p->plot_name);

  fprintf(_f, _("# Number of structures: %d\n"
                "#\n"), n_structures);

  for (struct_id = 0; struct_id < n_structures; struct_id++) {
    double m_tmp[9], d_tmp[9], s_tmp[9];
    for (i = 0; i < 9; i++) {
      m_tmp[i] = mass_matrixes[perm_id[i] + struct_id*9];
      d_tmp[i] = damping_matrixes[perm_id[i] + struct_id*9];
      s_tmp[i] = stiffness_matrixes[perm_id[i] + struct_id*9];
    }
    fprintf(_f, _("# Structure: %i\n"
                  "#\n"), struct_id + 1);
    fprintf(_f, _("# Mass:       [%14.7e, %14.7e, %14.7e]\n"
                  "#             [%14.7e, %14.7e, %14.7e]\n"
                  "#             [%14.7e, %14.7e, %14.7e]\n\n"),
            m_tmp[0], m_tmp[1], m_tmp[2],
            m_tmp[3], m_tmp[4], m_tmp[5],
            m_tmp[6], m_tmp[7], m_tmp[8]);
    fprintf(_f, _("# Damping:    [%14.7e, %14.7e, %14.7e]\n"
                  "#             [%14.7e, %14.7e, %14.7e]\n"
                  "#             [%14.7e, %14.7e, %14.7e]\n\n"),
            d_tmp[0], d_tmp[1], d_tmp[2],
            d_tmp[3], d_tmp[4], d_tmp[5],
            d_tmp[6], d_tmp[7], d_tmp[8]);
    fprintf(_f, _("# Stiffness:  [%14.7e, %14.7e, %14.7e]\n"
                  "#             [%14.7e, %14.7e, %14.7e]\n"
                  "#             [%14.7e, %14.7e, %14.7e]\n\n"),
            s_tmp[0], s_tmp[1], s_tmp[2],
            s_tmp[3], s_tmp[4], s_tmp[5],
            s_tmp[6], s_tmp[7], s_tmp[8]);
   }
  fprintf(_f,
          _("# (when structure characteristics are variable, the values\n"
            "# above are those at the computation initialization.\n\n"));

  fprintf(_f, _("# Columns:\n"));
  if (p->use_iteration)
    fprintf(_f, _("#   %d:     Time step number\n"), col_id++);
  else
    fprintf(_f, _("#   %d:     Physical time\n"), col_id++);
  fprintf(_f, _("#   %d - :  Values for each structure\n"), col_id);


  fprintf(_f,
          "#\n"
          "#TITLE: %s\n"
          "#COLUMN_TITLES: ", p->plot_name);
  if (p->use_iteration)
    fprintf(_f, " nt");
  else
    fprintf(_f, " t");
  for (i = 0; i < n_structures; i++)
    fprintf(_f, " | %d", i + 1);
  fprintf(_f, "\n");

  fprintf(_f, "#COLUMN_UNITS: ");
  if (p->use_iteration)
    fprintf(_f, " iter");
  else
    fprintf(_f, " s");
  for (i = 0; i < n_structures; i++)
    fprintf(_f, " -");
  fprintf(_f, "\n#\n");

  /* Close file or assign it ot handler depending on options */

  if (p->buffer_steps[0] > 0) {
    if (fclose(_f) != 0)
      bft_error(__FILE__, __LINE__, errno,
                _("Error closing file: \"%s\""), p->file_name);
  }
  else
    p->f = _f;
}

/*----------------------------------------------------------------------------
 * Write file header for CSV files
 *
 * parameters:
 *   p                  <-> time plot values file handler
 *   n_structures       <-- number of structures associated with this plot
 *----------------------------------------------------------------------------*/

static void
_write_struct_header_csv(cs_time_plot_t  *p,
                         int              n_structures)
{
  int i;
  FILE *_f = p->f;

  if (_f != NULL) {
    fclose(_f);
    p->f = NULL;
  }

  _f = fopen(p->file_name, "w");
  if (_f == NULL) {
    bft_error(__FILE__, __LINE__, errno,
              _("Error opening file: \"%s\""), p->file_name);
    return;
  }

  if (p->use_iteration)
    fprintf(_f, " iteration");
  else
    fprintf(_f, "t");
  for (i = 0; i < n_structures; i++) {
    fprintf(_f, ",%d", i + 1);
  }
  fprintf(_f, "\n");

  /* Close file or assign it ot handler depending on options */

  if (p->buffer_steps[0] > 0) {
    if (fclose(_f) != 0)
      bft_error(__FILE__, __LINE__, errno,
                _("Error closing file: \"%s\""), p->file_name);
  }
  else
    p->f = _f;
}

/*----------------------------------------------------------------------------
 * Add a time plot to the global time plots array.
 *----------------------------------------------------------------------------*/

static void
_time_plot_register(cs_time_plot_t  *p)
{
  p->prev = _plots_tail;
  p->next = NULL;

  if (_plots_head == NULL)
    _plots_head = p;
  else if (_plots_head->next == NULL)
    _plots_head->next = p;

  if (_plots_tail != NULL)
    _plots_tail->next = p;

  _plots_tail = p;
}

/*----------------------------------------------------------------------------
 * Remove a time plot from the global time plots array.
 *----------------------------------------------------------------------------*/

static void
_time_plot_unregister(cs_time_plot_t  *p)
{
  if (_plots_head == p)
    _plots_head = p->next;
  if (_plots_tail == p)
    _plots_tail = p->prev;

  if (p->prev != NULL)
    p->prev->next = p->next;
  if (p->next != NULL)
    p->next->prev = p->prev;
}

/*----------------------------------------------------------------------------
 * Create time plot writer for a given variable
 *
 * parameters:
 *   plot_name        <-- plot (variable) name
 *   file_prefix      <-- file name prefix
 *   format           <-- associated file format
 *   use_iteration    <-- should we use the iteration number instead of the
 *                        physical time ?
 *   flush_wtime      <-- elapsed time interval between file flushes
 *                        (if < 0, no forced flush)
 *   n_buffer_steps   <-- number of time steps in output buffer if
 *                        file is not to be kept open
 *
 * returns
 *   pointer to corresponding probe writer
 *----------------------------------------------------------------------------*/

static cs_time_plot_t *
_plot_file_create(const char             *plot_name,
                  const char             *file_prefix,
                  cs_time_plot_format_t   format,
                  bool                    use_iteration,
                  double                  flush_wtime,
                  int                     n_buffer_steps)
{
  size_t i;

  cs_time_plot_t *p = NULL;

  BFT_MALLOC(p, 1, cs_time_plot_t);
  BFT_MALLOC(p->plot_name, strlen(plot_name) + 1, char);
  BFT_MALLOC(p->file_name,
             strlen(file_prefix) + strlen(plot_name) + 4 + 1,
             char);

  strcpy(p->plot_name, plot_name);
  switch (format) {
  case CS_TIME_PLOT_DAT:
    sprintf(p->file_name, "%s%s.dat", file_prefix, plot_name);
    break;
  case CS_TIME_PLOT_CSV:
    sprintf(p->file_name, "%s%s.csv", file_prefix, plot_name);
    break;
  default:
    break;
  }

  for (i = strlen(file_prefix); p->file_name[i] != '\0'; i++) {
    if (isspace(p->file_name[i]))
      p->file_name[i] = '_';
  }

  p->f = NULL;

  p->format = format;
  p->use_iteration = use_iteration;

  p->flush_times[0] = flush_wtime;
  p->flush_times[1] = -2;

  p->buffer_steps[0] = n_buffer_steps;
  p->buffer_steps[1] = 0;

  p->buffer_size = 256;
  p->buffer_end = 0;

  BFT_MALLOC(p->buffer, p->buffer_size, char);

  _time_plot_register(p);

  return p;
}

/*----------------------------------------------------------------------------
 * Write buffered values to file if applicable
 *
 * parameters:
 *   p <-> time plot values file handler
 *----------------------------------------------------------------------------*/

static void
_plot_file_check_or_write(cs_time_plot_t  *p)
{
  size_t n_written;

  /* Return immediately if we are buffering and not writing now */

  if (   p->buffer_steps[0] > 0
      && p->buffer_steps[1] < p->buffer_steps[0]) {
    p->buffer_steps[1] += 1;
      return;
  }

  /* Ensure file is open */

  if (p->f == NULL) {
    p->f = fopen(p->file_name, "a");
    if (p->f == NULL) {
      bft_error(__FILE__, __LINE__, errno,
                _("Error re-opening file: \"%s\""), p->file_name);
      p->buffer_end = 0;
      return;
    }
  }

  /* Write buffer contents */
  n_written = fwrite(p->buffer, 1, p->buffer_end, p->f);

  if (n_written < p->buffer_end)
    bft_error(__FILE__, __LINE__, ferror(p->f),
              _("Error writing file: \"%s\""), p->file_name);

  p->buffer_end = 0;

  /* Close or flush file depending on options */

  if (p->buffer_steps[0] > 0) {

    assert(p->buffer_steps[1] >= p->buffer_steps[0]);
    if (fclose(p->f) != 0)
      bft_error(__FILE__, __LINE__, errno,
                  _("Error closing file: \"%s\""), p->file_name);
    p->f = NULL;
    p->buffer_steps[1] = 0;

  }
  else {
    double cur_time = cs_timer_wtime();
    if (   p->flush_times[0] > 0
        && (cur_time - p->flush_times[1]) > p->flush_times[0]) {
      p->flush_times[1] = cur_time;
      fflush(p->f);
    }
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a plot file writer for probe-type plots
 *
 * This function should only be called by one rank for a given data series.
 *
 * parameters:
 *   plot_name        <-- plot (variable) name
 *   file_prefix      <-- file name prefix
 *   format           <-- associated file format
 *   use_iteration    <-- should we use the iteration number instead of the
 *                        physical time ?
 *   flush_wtime      <-- elapsed time interval between file flushes
 *                        (if < 0, no forced flush)
 *   n_buffer_steps   <-- number of time steps in output buffer if
 *                        file is not to be kept open
 *   n_probes         <-- number of probes associated with this plot
 *   probe_list       <-- numbers (1 to n) of probes if filtered, or NULL
 *   probe_coords     <-- probe coordinates, or NULL
 *   probe_names      <-- probe names, or NULL
 *
 * returns:
 *   pointer to new time plot writer
 *----------------------------------------------------------------------------*/

cs_time_plot_t *
cs_time_plot_init_probe(const char             *plot_name,
                        const char             *file_prefix,
                        cs_time_plot_format_t   format,
                        bool                    use_iteration,
                        double                  flush_wtime,
                        int                     n_buffer_steps,
                        int                     n_probes,
                        const int              *probe_list,
                        const cs_real_t         probe_coords[],
                        const char             *probe_names[])
{
  cs_time_plot_t  *p = _plot_file_create(plot_name,
                                         file_prefix,
                                         format,
                                         use_iteration,
                                         flush_wtime,
                                         n_buffer_steps);

  switch (format) {
  case CS_TIME_PLOT_DAT:
    _write_probe_header_dat(p, n_probes, probe_list, probe_coords, probe_names);
    break;
  case CS_TIME_PLOT_CSV:
    _write_probe_coords_csv(file_prefix,
                            plot_name,
                            n_probes,
                            probe_list,
                            probe_coords);
    _write_probe_header_csv(p, n_probes, probe_list, probe_coords, probe_names);
    break;
  default:
    break;
  }

  cs_time_plot_flush(p); /* Flush after creation so plot tools can start
                            analyzing file immediately, while computation
                            is running */

  return p;
}

/*----------------------------------------------------------------------------
 * Initialize a plot file writer for structure-type plots
 *
 * This function should only be called by one rank for a given data series.
 *
 * parameters:
 *   plot_name          <-- plot (variable) name
 *   file_prefix        <-- file name prefix
 *   format             <-- associated file format
 *   use_iteration      <-- should we use the iteration number instead of the
 *                          physical time ?
 *   flush_wtime        <-- elapsed time interval between file flushes
 *                          (if < 0, no forced flush)
 *   n_buffer_steps     <-- number of time steps in output buffer if
 *                          file is not to be kept open
 *   n_structures       <-- number of structures associated with this plot
 *   mass_matrixes      <-- mass matrix coefficients (3x3 blocks)
 *   damping_matrixes   <-- damping matrix coefficients (3x3 blocks)
 *   stiffness_matrixes <-- stiffness matrix coefficients (3x3 blocks)
 *
 * returns:
 *   pointer to new time plot writer
 *----------------------------------------------------------------------------*/

cs_time_plot_t *
cs_time_plot_init_struct(const char             *plot_name,
                         const char             *file_prefix,
                         cs_time_plot_format_t   format,
                         bool                    use_iteration,
                         double                  flush_wtime,
                         int                     n_buffer_steps,
                         int                     n_structures,
                         const cs_real_t         mass_matrixes[],
                         const cs_real_t         damping_matrixes[],
                         const cs_real_t         stiffness_matrixes[])
{
  cs_time_plot_t  *p = _plot_file_create(plot_name,
                                         file_prefix,
                                         format,
                                         use_iteration,
                                         flush_wtime,
                                         n_buffer_steps);

  switch (format) {
  case CS_TIME_PLOT_DAT:
    _write_struct_header_dat(p, n_structures,
                             mass_matrixes, damping_matrixes, stiffness_matrixes);
    break;
  case CS_TIME_PLOT_CSV:
    _write_struct_header_csv(p, n_structures);
  break;
  default:
    break;
  }

  return p;
}

/*----------------------------------------------------------------------------
 * Finalize time plot writer for a given variable
 *
 * This function should only be called by one rank for a given data series.
 *
 * parameters:
 *   p <-> time plot values file handler
 *----------------------------------------------------------------------------*/

void
cs_time_plot_finalize(cs_time_plot_t  **p)
{
  if (p != NULL) {

    cs_time_plot_t *_p = *p;

    _time_plot_unregister(_p);

    if (_p->buffer_steps[0] > 0)
      _p->buffer_steps[1] = _p->buffer_steps[0] + 1;

    _plot_file_check_or_write(_p);

    if (_p->f != NULL) {
      if (fclose(_p->f) != 0)
        bft_error(__FILE__, __LINE__, errno,
                  _("Error closing file: \"%s\""), _p->file_name);
    }

    BFT_FREE(_p->buffer);
    BFT_FREE(_p->file_name);
    BFT_FREE(_p->plot_name);

    BFT_FREE(*p);
  }
}

/*----------------------------------------------------------------------------
 * Write time plot values
 *
 * This function should only be called by one rank for a given data series.
 *
 * parameters:
 *   p      <-- pointer to associated plot structure
 *   tn     <-- associated time step number
 *   t      <-- associated time value
 *   n_vals <-- number of associated time values
 *   vals   <-- associated time values
 *----------------------------------------------------------------------------*/

void
cs_time_plot_vals_write(cs_time_plot_t  *p,
                        int              tn,
                        double           t,
                        int              n_vals,
                        const cs_real_t  vals[])
{
  int i;

  if (p == NULL)
    return;

  /* Write data to line buffer */

  _ensure_buffer_size(p, p->buffer_end + 64);

  switch (p->format) {

  case CS_TIME_PLOT_DAT:

    if (p->use_iteration)
      p->buffer_end += sprintf(p->buffer + p->buffer_end, " %8d", tn);
    else
      p->buffer_end += sprintf(p->buffer + p->buffer_end, " %14.7e", t);

    for (i = 0; i < n_vals; i++) {
      _ensure_buffer_size(p, p->buffer_end + 64);
      p->buffer_end += sprintf(p->buffer + p->buffer_end, " %14.7e",
                               (double)(vals[i]));
    }

    p->buffer_end += sprintf(p->buffer + p->buffer_end, "\n");

    break;

  case CS_TIME_PLOT_CSV:

    if (p->use_iteration)
      p->buffer_end += sprintf(p->buffer + p->buffer_end, "%8d", tn);
    else
      p->buffer_end += sprintf(p->buffer + p->buffer_end, "%14.7e", t);

    for (i = 0; i < n_vals; i++) {
      _ensure_buffer_size(p, p->buffer_end + 64);
      p->buffer_end += sprintf(p->buffer + p->buffer_end, ", %14.7e",
                               (double)(vals[i]));
    }

    p->buffer_end += sprintf(p->buffer + p->buffer_end, "\n");

    break;

  default:
    break;
  }

  /* Output buffer if necessary */

  _plot_file_check_or_write(p);
}

/*----------------------------------------------------------------------------
 * Flush buffered values to file if applicable
 *
 * parameters:
 *   p <-> time plot values file handler
 *----------------------------------------------------------------------------*/

void
cs_time_plot_flush(cs_time_plot_t  *p)
{
  /* Force buffered variant output */

  if (p->buffer_end > 0) {
    if (p->buffer_steps[0] > 0)
      p->buffer_steps[1] = p->buffer_steps[0];
    _plot_file_check_or_write(p);
  }

  if (p->f != NULL) {
    if (p->flush_times[0] > 0)
      p->flush_times[1] = cs_timer_wtime();
    fflush(p->f);
  }
}

/*----------------------------------------------------------------------------
 * flush all time plots
 *----------------------------------------------------------------------------*/

void
cs_time_plot_flush_all(void)
{
  for (cs_time_plot_t *p = _plots_head; p != NULL; p = p->next)
    cs_time_plot_flush(p);
}

/*----------------------------------------------------------------------------
 * Set time plot file writer flush behavior defaults.
 *
 * parameters:
 *   flush_wtime     <-- elapsed time interval between file flushes;
 *                       if < 0, no forced flush
 *   n_buffer_steps  <-- number of time steps in output buffer if
 *                       file is not to be kept open
 *----------------------------------------------------------------------------*/

void
cs_time_plot_set_flush_default(float  flush_wtime,
                               int    n_buffer_steps)
{
  _flush_wtime_default    = flush_wtime;
  _n_buffer_steps_default = n_buffer_steps;
}

/*----------------------------------------------------------------------------
 * Return time plot file writer flush behavior defaults.
 *
 * parameters:
 *   flush_wtime     --> elapsed time interval between file flushes;
 *                       if < 0, no forced flush (NULL if not queried)
 *   n_buffer_steps  <-- number of time steps in output buffer if
 *                       file is not to be kept open (NULL if not queried)
 *----------------------------------------------------------------------------*/

void
cs_time_plot_get_flush_default(float  *flush_wtime,
                               int    *n_buffer_steps)
{
  if (flush_wtime != NULL)
    *flush_wtime = _flush_wtime_default;

  if (n_buffer_steps != NULL)
    *n_buffer_steps = _n_buffer_steps_default;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
