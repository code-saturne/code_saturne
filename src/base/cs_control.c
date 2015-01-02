/*============================================================================
 * Interactive control management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_UNISTD_H) && defined(HAVE_ACCESS)
#include <unistd.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_file.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_restart.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_control.h"

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_control.c
 *
 *  \brief Handle control file usable for interactive change of stop,
 *         post-processing or checkpoint behavior.
 */

/*============================================================================
 * Static global variables
 *============================================================================*/

static double  _control_file_wt_interval = 0.;
static double  _control_file_wt_last = -1.;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read next value, expecting integer
 *
 * parameters:
 *   cur_line <-- pointer to line buffer
 *   s        <-> current string position
 *   val      --> integer read
 *
 * returns:
 *   number of integers read (1 for success, 0 otherwise)
 *----------------------------------------------------------------------------*/

static int
_read_next_int(const char   *cur_line,
               const char  **s,
               int          *val)
{
  int n_val = 0;

  const char *p = *s;
  while (*p != '\0' && *p != ' ' && *p != '\t')
    p++;
  while (*p != '\0' && (*p == ' ' || *p == '\t'))
    p++;
  *s = p;

  n_val = sscanf(*s, "%i", val);

  if (n_val == 0)
    bft_printf(_("   ignored: \"%s\"\n"), cur_line);

  return n_val;
}

/*----------------------------------------------------------------------------
 * Read next optional value, expecting integer
 *
 * parameters:
 *   s        <-> current string position
 *   val     --> integer read
 *
 * returns:
 *   number of integers read (1 for success, 0 otherwise)
 *----------------------------------------------------------------------------*/

static int
_read_next_opt_int(const char  **s,
                   int          *val)
{
  int n_val = 0;

  const char *p = *s;
  while (*p != '\0' && *p != ' ' && *p != '\t')
    p++;
  while (*p != '\0' && (*p == ' ' || *p == '\t'))
    p++;
  *s = p;

  n_val = sscanf(*s, "%i", val);

  return n_val;
}

/*----------------------------------------------------------------------------
 * Read next value, expecting double precision floating point value
 *
 * parameters:
 *   cur_line <-- pointer to line buffer
 *   s        <-> current string position
 *   val      --> value read
 *
 * returns:
 *   number of values read (1 for success, 0 otherwise)
 *----------------------------------------------------------------------------*/

static int
_read_next_double(const char   *cur_line,
                  const char  **s,
                  double       *val)
{
  int n_val = 0;

  const char *p = *s;
  while (*p != '\0' && *p != ' ' && *p != '\t')
    p++;
  while (*p != '\0' && (*p == ' ' || *p == '\t'))
    p++;
  *s = p;

  n_val = sscanf(*s, "%lg", val);

  if (n_val == 0)
    bft_printf(_("   ignored: \"%s\"\n"), cur_line);

  return n_val;
}

/*----------------------------------------------------------------------------
 * Handle command file line relative to checkpointing
 *
 * parameters:
 *   cur_line <-- pointer to line buffer
 *   s        <-> pointer to current position in line
 *----------------------------------------------------------------------------*/

static void
_control_checkpoint(const char   *cur_line,
                    const char  **s)
{
  *s += 11; /* shift in string by lenght of "checkpoint_" part */

  if (strncmp(*s, "time_step ", 10) == 0) {
    int nt;
    if (_read_next_int(cur_line, s, &nt) > 0) {
      cs_restart_checkpoint_set_next_ts(nt);
      bft_printf("  %-32s %12d\n",
                 "checkpoint_time_step", nt);
    }
  }
  else if (strncmp(*s, "time_value ", 11) == 0) {
    double t;
    if (_read_next_double(cur_line, s, &t) > 0) {
      cs_restart_checkpoint_set_next_tv(t);
      bft_printf("  %-32s %12.5g\n",
                 "checkpoint_time_value", t);
    }
  }
  else if (strncmp(*s, "wall_time ", 10) == 0) {
    double wt;
    if (_read_next_double(cur_line, s, &wt) > 0) {
      cs_restart_checkpoint_set_next_wt(wt);
      bft_printf("  %-32s %12.5g\n",
                 "checkpoint_wall_time", wt);
    }
  }
  else if (strncmp(*s, "time_step_interval ", 19) == 0) {
    int nt;
    if (_read_next_int(cur_line, s, &nt) > 0) {
      cs_restart_checkpoint_set_defaults(nt, -1., -1.);
      bft_printf("  %-32s %12d\n",
                 "checkpoint_time_step_interval", nt);
    }
  }
  else if (strncmp(*s, "time_value_interval ", 20) == 0) {
    double t;
    if (_read_next_double(cur_line, s, &t) > 0) {
      if (t > 0) {
        cs_restart_checkpoint_set_defaults(-1, t, -1.);
        bft_printf("  %-32s %12.5g\n",
                   "checkpoint_time_value_interval", t);
      }
      else
        bft_printf("  %-32s %12.5g %s\n",
                   "checkpoint_time_value_interval", t, _("ignored"));
    }
  }
  else if (strncmp(*s, "wall_time_interval ", 19) == 0) {
    double wt;
    if (_read_next_double(cur_line, s, &wt) > 0) {
      if (wt > 0) {
        cs_restart_checkpoint_set_defaults(-1, -1., wt);
        bft_printf("  %-32s %12.5g\n",
                   "checkpoint_wall_time_interval", wt);
      }
      else
        bft_printf("  %-32s %12.5g %s\n",
                   "checkpoint_wall_time_interval", wt, _("ignored"));
    }
  }
  else
    bft_printf(_("   ignored: \"%s\"\n"), cur_line);

}

/*----------------------------------------------------------------------------
 * Handle command file line relative to postprocessing
 *
 * parameters:
 *   ts       <-- pointer to time step status
 *   s        <-> pointer to current position in line
 *----------------------------------------------------------------------------*/

static void
_control_postprocess(const cs_time_step_t   *ts,
                     char                   *cur_line,
                     const char            **s)
{
  *s += 12; /* shift in string by lenght of "postprocess_" part */

  if (strncmp(*s, "time_step ", 10) == 0) {
    int nt, writer_id;
    if (_read_next_int(cur_line, s, &nt) > 0) {
      if (_read_next_opt_int(s, &writer_id) == 0)
        writer_id = 0;
      if (nt >= 0)
        nt = CS_MAX(nt, ts->nt_cur);
      else
        nt = CS_MAX(nt, -ts->nt_cur);
      cs_post_add_writer_t_step(writer_id, nt);
      bft_printf("  %-32s %12d %12d\n",
                 "postprocess_time_step", nt, writer_id);
    }
  }
  else if (strncmp(*s, "time_value ", 11) == 0) {
    int writer_id;
    double t;
    if (_read_next_double(cur_line, s, &t) > 0) {
      if (_read_next_opt_int(s, &writer_id) == 0)
        writer_id = 0;
      if (t >= 0)
        t = CS_MAX(t, ts->t_cur);
      else
        t = CS_MAX(t, -ts->t_cur);
      bft_printf("  %-32s %12.5g %12d\n",
                 "postprocess_time_value", t, writer_id);
      cs_post_add_writer_t_value(writer_id, t);
    }
  }
  else
    bft_printf(_("   ignored: \"%s\"\n"), cur_line);

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the presence of a control file and deal with the interactive
 *        control.
 */
/*----------------------------------------------------------------------------*/

void
cs_control_check_file(void)
{
  int nt_max;

  long f_size = 0;
  char *s;
  char *buffer = NULL, *cur_line = NULL;
  const cs_time_step_t  *ts = cs_glob_time_step;

  const char path[] = "control_file";

  /* Test existence and size of file */

  if (cs_glob_rank_id <= 0) {

    if (   _control_file_wt_interval <= 0.
        ||(    cs_timer_wtime() - _control_file_wt_last
           >= _control_file_wt_interval)) {

#if defined(HAVE_UNISTD_H) && defined(HAVE_ACCESS)

      /* Test existence of file using access() before stat(),
         as this may be less costly on some filesytems
         (such as on LUSTRE, due to metadata handling aspects). */

      if (access(path, F_OK) == 0)
        f_size = cs_file_size(path);

#else

      f_size = cs_file_size(path);

#endif

    }

  }

#if defined(HAVE_MPI)
  if (cs_glob_rank_id >= 0)
    MPI_Bcast(&f_size,  1, MPI_LONG,  0, cs_glob_mpi_comm);
#endif

  /* If no control file is present, we are done */

  if (f_size == 0)
    return;

  /* If file exists, handle it */

  BFT_MALLOC(buffer, f_size + 1, char);

  if (cs_glob_rank_id <= 0) {

    size_t r_size = 0;
    FILE *control_file = fopen("control_file", "r");

    if (control_file != NULL) {
      r_size = fread(buffer, 1, f_size, control_file);
      buffer[r_size] = '\0'; /* precaution for partial read */
      fclose(control_file);
      remove("control_file");
    }
    else
      bft_printf
        (_("\n"
           " Warning: error opening %s (ignored):\n"
           " --------\n"
           "   \"%s\"\n\n"), path, strerror(errno));


    _control_file_wt_last = cs_timer_wtime();

  }

#if defined(HAVE_MPI)
  if (cs_glob_rank_id >= 0)
    MPI_Bcast(buffer, f_size + 1, MPI_CHAR, 0, cs_glob_mpi_comm);
#endif

  /* Now all ranks have required buffer */

  cur_line = buffer;

  bft_printf
    (_("\n"
       " Options set or changed by \"control_file\":\n"
       " -----------------------------------------\n\n"));

  /* Loop on buffer's lines */

  /* Note that when parsing lines, we do not use strtok() type functions, to
     avoid modifying a line (so that log/error/warning messages may simply use
     that line); hence also tests using strncp on strings appended with a
     whitespace (always present in case of arguments) rather than strcmp. */

  while (cur_line != NULL) {

    /* Prepare current and next line for parsing */

    char *next_line = cur_line;

    while (   *next_line != '\0'
           && *next_line != '\n' && *next_line != '\r' && *next_line != '\f')
      next_line++;

    *next_line = '\0'; /* marks end of cur_line */
    next_line += 1;

    if (next_line >= (buffer + f_size))
      next_line = NULL;
    else {
      while (    *next_line != '\0'
             && (*next_line == '\n' || *next_line == '\r' || *next_line == '\f'))
        next_line++;
    }

    /* Check for keywords given in control_file and store the related values */

    size_t l = strlen(cur_line);

    for (size_t i = 0; i < l; i++) {
      if (cur_line[i] == '#') {
        cur_line[i] = '\0';
        break;
      }
      else
        cur_line[i] = tolower(cur_line[i]);
    }

    for (s = cur_line; *s == ' ' || *s == '\t'; s++)
      *s = ' ';

    if (*s == '\0') {
      cur_line = next_line;
      continue;
    }

    /* Calculation end options
       default with no keyword; max_time_step */

    nt_max = -1;
    if (sscanf(s, "%i", &nt_max) > 0)
      nt_max = CS_MAX(nt_max, 0);
    else if (strncmp(s, "max_time_step ", 14) == 0) {
      if (_read_next_int(cur_line, (const char **)&s, &nt_max) > 0)
        nt_max = CS_MAX(nt_max, 0);
    }

    if (nt_max > -1) {
      nt_max = CS_MAX(nt_max, ts->nt_cur);
      cs_time_step_define_nt_max(nt_max);
      bft_printf("  %-32s %12d (%s %d)\n",
                 "max_time_step", ts->nt_max, _("current:"), ts->nt_cur);
    }
    else if (strncmp(s, "max_time_value ", 15) == 0) {
      double t_max;
      if (_read_next_double(cur_line, (const char **)&s, &t_max) > 0)
        t_max = CS_MAX(t_max, ts->t_cur);
      cs_time_step_define_t_max(t_max);
      bft_printf("  %-32s %12.5g (%s %12.5g)\n",
                 "max_time_value", ts->t_max, _("current:"), ts->t_cur);
    }

    /* Control file check interval */

    else if (strncmp(s, "control_file_wtime_interval ", 28) == 0) {
      double wt;
      if (_read_next_double(cur_line, (const char **)&s, &wt) > 0)
        _control_file_wt_interval = wt;
    }

    /* Checkpointing options */

    else if (strncmp(s, "checkpoint_", 11) == 0)
      _control_checkpoint(cur_line, (const char **)&s);

    /* Postprocessing options */

    else if (strncmp(s, "postprocess_", 12) == 0)
      _control_postprocess(ts, cur_line, (const char **)&s);

    /* Unhandled lines */

    else
      bft_printf(_("   ignored: \"%s\"\n"), cur_line);

    /* Prepare for next line */

    cur_line = next_line;

  } /* End of loop on lines */

  bft_printf
    (_("\n"
       " Finished reading \"control_file\".\n\n"));

  BFT_FREE(buffer);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
