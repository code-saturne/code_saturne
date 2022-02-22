/*============================================================================
 * User initialization prior to solving time steps.
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
 * \file cs_user_initialization.c
 *
 * \brief Initialization prior to solving time steps.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization.c
 *
 * \brief Initialize variables.
 *
 * This function is called at beginning of the computation
 * (restart or not) before the time step loop.
 *
 * This is intended to initialize or modify (when restarted)
 * variable and time step values.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t     *domain)
{
  CS_UNUSED(domain);

  /*! [remapper_init] */
  /* Apply only on computation start, not on restarts */
  if (domain->time_step->nt_prev > 0)
    return;

#if defined(HAVE_MEDCOUPLING_LOADER)

  double t1 = cs_timer_wtime();

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  bft_printf("In cs_user_initialization HAVE_MEDCOUPLING_LOADER\n");

  /* Number of fields to interpolate from the MED file */
  const int  nremapper_fields = 2;

  /* Names of the nremapper_fields fields to read */
  const char  *field_names[] = {"p", "U"};

  /* We request a remapper with a given name. If it does not exist,
   * the function returns a NULL pointer. */
  cs_medcoupling_remapper_t *r
    = cs_medcoupling_remapper_by_name_try("init");

  if (r == NULL) {

    /* Space dimension of the elements (2 for faces, 3 for cells) */
    int elts_dim = 3;

    /* Indexes needed to read the time step from the
     * file (0, 0 if only one exists) */
    int it0 = 0, it1 = 0;

    /* Path to file */
    const char file_name[] = "../../../MESH/channel395_OF_4.1_5000.med";

    /* The remapper is created. We retrieve its id from the function.
       * The function inputs are:
       * 1) Name of the remapper
       * 2) dimension of the mesh elements
       * 3) selection criteria
       * 4) path to the med file
       * 5) number of fields to interpolate
       * 6) names of the fields to interpolate
       * 7 + 8) time iteration index and order */
    int r_id = cs_medcoupling_remapper_initialize("init",
                                                  elts_dim,
                                                  "all[]",
                                                  file_name,
                                                  nremapper_fields,
                                                  field_names,
                                                  it0,
                                                  it1);

    /* Retrieve the pointer */
    r = cs_medcoupling_remapper_by_id(r_id);

    /* Set options */
    cs_medcoupling_remapper_set_options(r,
                                        "IntersectionType",
                                        "PointLocator");

    /* We create the interpolation matrix => Here it is only called once
     * since the mesh is not moving */
    cs_medcoupling_remapper_setup(r);

  }

  /* We retrieve arrays containing the interpolated values.
   * Inputs are:
   * 1) remapper
   * 2) id of the field to interpolate
   * 3) a default value (if no intersection is obtained) */

  /* We loop in the fields matching field_names defined above */

  for (int f_id = 0; f_id < nremapper_fields; f_id++) {

    cs_real_t *m_vals = cs_medcoupling_remapper_copy_values(r, f_id, 1.);

    /* For pressure */

    if (f_id == 0) { /* or (strcmp(field_names[f_id], "p") == 0) */

      cs_real_t *cvar_p = CS_F_(p)->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cvar_p[c_id] = m_vals[c_id];
      }

    }

    /* For velocity */

    else if (f_id == 1) { /* or (strcmp(field_names[f_id], "U") == 0) */

      cs_real_3_t *cvar_vel = (cs_real_3_t *)CS_F_(vel)->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cvar_vel[c_id][0] = m_vals[c_id*3];
        cvar_vel[c_id][1] = m_vals[c_id*3 + 1];
        cvar_vel[c_id][2] = m_vals[c_id*3 + 2];
      }
    }

    BFT_FREE(m_vals);

  } /* End of loop on fields */

  double t2 = cs_timer_wtime();

  bft_printf("Time in MEDCoupling: %f seconds \n", t2 - t1);
  bft_printf_flush();

#endif

  /*! [remapper_init] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
