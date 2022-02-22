/*============================================================================
 * Tracing utility functions for profiling and debugging
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

#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_mem_usage.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_trace.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

#if defined(HAVE_MPI)

typedef struct
{
  double val;
  int    rank;
} _fvm_trace_mpi_double_int_t;

#endif

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Print memory usage status.
 *
 * If no associated description is given, a call number will be issued.
 *
 * parameters:
 *   descr <-- optional description, or NULL
 *----------------------------------------------------------------------------*/

void
fvm_trace_mem_status(const char  *descr)
{
  int    i, itot;
  double valreal[4];

#if defined(HAVE_MPI)
  MPI_Comm comm = cs_glob_mpi_comm;
  int rank_id = cs_glob_rank_id;
  int n_ranks = cs_glob_n_ranks;
  int  imax = 0, imin = 0;
  int  id_min[4];
  _fvm_trace_mpi_double_int_t  val_in[4], val_min[4], val_max[4];
#else
  int n_ranks = 1;
#endif

  int   val_flag[4] = {1, 1, 1, 1};
  char  unit[]    = {'k', 'm', 'g', 't', 'p'};

  static int call_id = 0;

  const char  *type_str[] = {"max. measured       ",
                             "max. instrumented   ",
                             "current measured    ",
                             "current instrumented"};

  /* Memory summary */

  if (descr != NULL)
    bft_printf(_("\nMemory use summary: %s\n\n"), descr);
  else
    bft_printf(_("\nMemory use summary (call %d):\n\n"), call_id);

  valreal[0] = (double)bft_mem_usage_max_pr_size();
  valreal[1] = (double)bft_mem_size_max();
  valreal[2] = (double)bft_mem_usage_pr_size();
  valreal[3] = (double)bft_mem_size_current();

  /* Ignore inconsistent measurements */

  for (i = 0; i < 4; i++) {
    if (valreal[i] < 1.0)
      val_flag[i] = 0;
  }

#if defined(HAVE_MPI)

  if (n_ranks > 1) {

    MPI_Reduce(val_flag, id_min, 4, MPI_INT, MPI_MIN, 0, comm);
    for (i = 0; i < 4; i++) {
      val_in[i].val = valreal[i];
      val_in[i].rank = rank_id;
    }
    MPI_Reduce(&val_in, &val_min, 4, MPI_DOUBLE_INT, MPI_MINLOC, 0, comm);
    MPI_Reduce(&val_in, &val_max, 4, MPI_DOUBLE_INT, MPI_MAXLOC, 0, comm);
    if (rank_id == 0) {
      for (i = 0; i < 4; i++) {
        val_flag[i]  = id_min[i];
        valreal[i] = val_max[i].val;
      }
    }
  }
#endif

  /* Similar handling of several instrumentation methods */

  for (i = 0; i < 4; i++) {

    /* If an instrumentation method returns an apparently consistent
       result, print it. */

    if (val_flag[i] == 1) {

      for (itot = 0;
           valreal[i] > 1024. && unit[itot] != 'p';
           itot++)
        valreal[i] /= 1024.;
#if defined(HAVE_MPI)
      if (n_ranks > 1 && rank_id == 0) {
        for (imin = 0;
             val_min[i].val > 1024. && unit[imin] != 'p';
             imin++)
          val_min[i].val /= 1024.;
        for (imax = 0;
             val_max[i].val > 1024. && unit[imax] != 'p';
             imax++)
          val_max[i].val /= 1024.;
      }
#endif

      /* Print to log file */

#if defined(HAVE_MPI)
      if (n_ranks > 1 && rank_id == 0) {
        bft_printf(_("  %s : %10.3f %cb min (rank %d), "
                     " %10.3f %cb max (rank %d)\n"),
                   type_str[i], val_min[i].val, unit[imin], val_min[i].rank,
                   val_max[i].val, unit[imax], val_max[i].rank);
      }
#endif
      if (n_ranks == 1)
        bft_printf(_("  %s : %12.3f %cb\n"),
                   type_str[i], valreal[i], unit[itot]);
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
