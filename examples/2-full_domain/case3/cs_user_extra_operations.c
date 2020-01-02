/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *============================================================================*/

/* VERS */

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
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

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
 */
/*----------------------------------------------------------------------------*/

/* Global declaration of the moy.dat file that will be filled with Tav */
static FILE *file = NULL;

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  /* Variables declaration */
  /* Get pointers to the mesh and mesh quantities structures */
  const cs_mesh_t *m= cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  /* Number of cells */
  const int n_cells = m->n_cells;

  /* Cell volumes */
  const cs_real_t *cell_vol = fvq->cell_vol;

  /* Get the temperature field */
  const cs_field_t *temp = cs_field_by_name_try("temperature");

  /* Cell volumes sum, Temp * volumes sum and Tavg */
  cs_real_t voltot = 0., temptot =0., Tavg =0.;

  /* Compute the sum of the cell volumes */
  for (int ii = 0 ; ii < n_cells ; ii++)
    voltot += cell_vol[ii];

  /* Compute the sum T*vol */
  for (int ii = 0 ; ii < n_cells ; ii++)
    temptot += temp->val[ii]*cell_vol[ii];

  /* Parallel sums */
  cs_parall_sum(1, CS_DOUBLE, &voltot);
  cs_parall_sum(1, CS_DOUBLE, &temptot);

  /* Compute Tavg */
  Tavg = temptot / voltot ;

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
