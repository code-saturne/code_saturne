/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *============================================================================*/

/* VERS */

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
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function)
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  /* Variables declaration */

  /* Handle to the moy.dat file that will be filled with Tav;
   declared as static so as to keep its value between calls. */
  static FILE *file = NULL;

  /* Get pointers to the mesh and mesh quantities structures */
  const cs_mesh_t *m= cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  /* Number of cells */
  const int n_cells = m->n_cells;

  /* Cell volumes */
  const cs_real_t *cell_vol = fvq->cell_vol;

  /* Get the temperature field */
  const cs_field_t *temp = cs_field_by_name_try("temperature");

  /* Total domain volume */
  cs_real_t voltot = fvq->tot_vol;

  /* Temp * volumes sum and Tavg */
  cs_real_t temptot =0., Tavg =0.;

  /* Compute the sum T*vol */
  for (int ii = 0 ; ii < n_cells ; ii++)
    temptot += temp->val[ii]*cell_vol[ii];

  /* Parallel sums */
  cs_parall_sum(1, CS_REAL_TYPE, &temptot);

  /* Compute Tavg */
  Tavg = temptot / voltot;

  /* Open the file moy.dat at the first iteration
     and write the first comment line only on the
     first processor (0 in parallel, -1 in serial) */
  if (cs_glob_time_step->nt_cur == 1 && cs_glob_rank_id <= 0) {
    file = fopen("moy.dat", "a");
    fprintf(file, "#Time (s)   Average Temperature (C) \n");
  }

/* Print the average temperature at the current time step
   on the first processor only */
  if (cs_glob_rank_id <= 0)
    fprintf(file, "%.6f  %.6f\n",cs_glob_time_step->t_cur, Tavg);

/* Close the file moy.dat at the last iteration
   on the first processor only */
  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_max
   && cs_glob_rank_id <= 0)
    fclose(file);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
