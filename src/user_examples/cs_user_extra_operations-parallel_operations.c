/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations-parallel_operations.c
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

  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
  const cs_lnum_t n_i_faces = domain->mesh->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)domain->mesh->i_face_cells;

  const cs_real_t *cell_vol = domain->mesh_quantities->cell_vol;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)domain->mesh_quantities->cell_cen;

  const cs_real_3_t *cvar_vel = (const cs_real_3_t *)CS_F_(vel)->val;

  /* Example use of parallel utility functions for several operations
   * ================================================================ */

  /* Sum of an integer counter 'ii', here the number of cells */

  /* ![example_1] */
  cs_gnum_t g_ii = n_cells;

  cs_parall_sum(1, CS_GNUM_TYPE, &g_ii);

  bft_printf("%s: total number of cells = %ld\n", __func__, (long)g_ii);
  /* ![example_1] */

  /* Maximum of an integer counter 'ii', here the number of cells */

  /* ![example_2] */
  cs_lnum_t ii = n_cells;

  cs_parall_max(1, CS_LNUM_TYPE, &ii);

  bft_printf("%s: max. number of cells per rank = %d\n", __func__, ii);
  /* ![example_2] */

  /* Sum of a real 'rrr', here the volume; */

  /* ![example_3] */
  cs_real_t rrr = cs_sum(n_cells, cell_vol);

  cs_parall_sum(1, CS_REAL_TYPE, &rrr);

  bft_printf("%s: total domain volume = %14.5e\n", __func__, rrr);
  /* ![example_3] */

  /* Maximum of a real 'rrr', here the volume */

  /* ![example_4] */
  rrr = 0;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (rrr < cell_vol[c_id])
      rrr = cell_vol[c_id];
  }

  cs_parall_max(1, CS_REAL_TYPE, &rrr);

  bft_printf("%s: max cell volume = %14.5e\n", __func__, rrr);
  /* ![example_4] */

  /* Minimum of a real 'rrr', here the volume */

  /* ![example_5] */
  rrr = cs_math_big_r;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (rrr > cell_vol[c_id])
      rrr = cell_vol[c_id];
  }

  cs_parall_min(1, CS_REAL_TYPE, &rrr);

  bft_printf("%s: min cell volume = %14.5e\n", __func__, rrr);
  /* ![example_5] */

 /* Maximum of a real and associated real values;
  * here the volume and its location (3 coordinates) */

  /* ![example_6] */
  rrr = -1;
  cs_real_t xyz[3] = {0, 0, 0};

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (rrr < cell_vol[c_id]) {
      rrr = cell_vol[c_id];
      for (int i = 0; i < 3; i++)
        xyz[i] = cell_cen[c_id][i];
    }
  }

  cs_parall_max_loc_vals(3, &rrr, xyz);

  bft_printf("%s: Max. volume = %14.5e.\n", __func__, rrr);
  bft_printf("Location(x,y,z) = %14.5e, %14.5e, %14.5e\n",
             xyz[0], xyz[1], xyz[2]);
  /* ![example_6] */

  /* Minimum of a real and associated real values;
   * here the volume and its location (3 coordinates) */

  /* ![example_7] */
  rrr = 1e30;
  xyz[0] = 0;
  xyz[1] = 0;
  xyz[2] = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (rrr > cell_vol[c_id]) {
      rrr = cell_vol[c_id];
      for (int i = 0; i < 3; i++)
        xyz[i] = cell_cen[c_id][i];
    }
  }

  cs_parall_min_loc_vals(3, &rrr, xyz);

  bft_printf("%s: Min. volume = %14.5e.\n ", __func__, rrr);
  bft_printf("           Location (x,y,z) = %14.5e, %14.5e, %14.5e\n",
             xyz[0], xyz[1], xyz[2]);
  /* ![example_7] */

  /* Sum of an array of integers;
   * here, the number of cells, faces, and boundary faces

   * local values; note that to avoid counting interior faces on
   * parallel boundaries twice, we check if "i_face_cells[f_id][0] <= n_cells",
   * as on a parallel boundary, this is always true for one domain
   * and false for the other. */

  /* ![example_8] */
  cs_gnum_t g_itab[3] = {n_cells, 0, n_b_faces};

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    if (i_face_cells[f_id][0] <= n_cells)
      g_itab[1]++;

  cs_parall_sum(3, CS_GNUM_TYPE, g_itab);

  bft_printf("%s:\n"
             "Number of cells = %ld\n"
             "Number of interior faces = %ld\n"
             "Number of boundary faces = %ld\n\n",
             __func__, (long)g_itab[0], (long)g_itab[1], (long)g_itab[2]);
  /* ![example_8] */

  /* Maxima from an array of integers;
   * here, the number of cells, faces, and boundary faces */

  /* ![example_9] */
  cs_lnum_t itab[3];
  itab[0] = n_cells;
  itab[1] = n_i_faces;
  itab[2] = n_b_faces;

  /* global maxima */
  cs_parall_max(3, CS_LNUM_TYPE, itab);

  bft_printf("%s:\n"
             " Max. number of cells per rank = %d\n"
             " Max. number of interior faces per rank = %d\n"
             " Max. number of boundary faces per rank = %d\n\n",
             __func__, (int)itab[0], (int)itab[1], (int)itab[2]);

  /* ![example_9] */

  /* Minima from an array of integers;
   * here, the number of cells, faces, and boundary faces */

  /* ![example_10] */
  itab[0] = n_cells;
  itab[1] = n_i_faces;
  itab[2] = n_b_faces;

  /* global minima */
  cs_parall_min(3, CS_LNUM_TYPE, itab);

  bft_printf("%s:\n"
             " Min. number of cells per rank = %d\n"
             " Min. number of interior faces per rank = %d\n"
             " Min. number of boundary faces per rank = %d\n\n",
             __func__, (int)itab[0], (int)itab[1], (int)itab[2]);
  /* ![example_10] */

  /* Sum of an array of reals;
   * here, the 3 velocity components (so as to compute a mean for example);
   * note that a simple loop is a bad way of computing a large sum;
   * hierarchical or other precise BLAS-type operations should be used for
   * better numerical precision. */

  /* ![example_11] */
  xyz[0] = 0;
  xyz[1] = 0;
  xyz[2] = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (cs_lnum_t i = 0; i < 3; i++)
      xyz[i] += cvar_vel[c_id][i];
  }

  /* global sum */
  cs_parall_sum(3, CS_REAL_TYPE, xyz);

  bft_printf("%s:\n"
             " Sum of U on the domain = %14.5e\n"
             " Sum of V on the domain = %14.5e\n"
             " Sum of V on the domain = %14.5e\n\n",
             __func__, xyz[0], xyz[1], xyz[2]);
  /* ![example_11] */

  /* Maximum of an array of reals;
   * here, the 3 velocity components */

  /* ![example_12] */
  xyz[0] = cvar_vel[0][0];
  xyz[1] = cvar_vel[0][1];
  xyz[2] = cvar_vel[0][2];

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (cs_lnum_t i = 0; i < 3; i++)
      xyz[i] = cs_math_fmax(cvar_vel[c_id][i], xyz[i]);
  }

  /* global maximum */
  cs_parall_max(3, CS_REAL_TYPE, xyz);

  bft_printf("%s:\n"
             " Max. of U on the domain = %14.5e\n"
             " Max. of V on the domain = %14.5e\n"
             " Max. of V on the domain = %14.5e\n\n",
             __func__, xyz[0], xyz[1], xyz[2]);
  /* ![example_12] */

  /* Minimum of an array of reals;
   * here, the 3 velocity components */

  /* ![example_13] */
  xyz[0] = cvar_vel[0][0];
  xyz[1] = cvar_vel[0][1];
  xyz[2] = cvar_vel[0][2];

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (int i = 0; i < 3; i++)
      xyz[i] = cs_math_fmin(cvar_vel[c_id][i], xyz[i]);
  }

  /* global maximum */
  cs_parall_min(3, CS_REAL_TYPE, xyz);

  bft_printf("%s:\n"
             " Min. of U on the domain = %14.5e\n"
             " Min. of V on the domain = %14.5e\n"
             " Min. of V on the domain = %14.5e\n\n",
             __func__, xyz[0], xyz[1], xyz[2]);
  /* ![example_13] */

  /* Broadcast an array of local integers to other ranks;
   * in this example, we use the number of cells, interior faces, and boundary
   * faces from process rank 0 (root_rank). */

  /* ![example_14] */
  int root_rank = 0;
  itab[0] = n_cells;
  itab[1] = n_i_faces;
  itab[2] = n_b_faces;

  /* broadcast from root_rank to all others */
  cs_parall_bcast(root_rank, 3, CS_LNUM_TYPE, itab);

  bft_printf("%s: On rank %d\n"
             "              Number of cells          = %d\n"
             "              Number of interior faces = %d\n"
             "              Number of boundary faces = %d\n\n",
             __func__, root_rank, (int)itab[0], (int)itab[1], (int)itab[2]);
  /* ![example_14] */

  /* Broadcast an array of local reals to other ranks;
   * in this example, we use 3 velocity values from process root_rank = 0. */

  /* ![example_15] */
  xyz[0] = cvar_vel[0][0];
  xyz[1] = cvar_vel[0][1];
  xyz[2] = cvar_vel[0][2];

  /* broadcast from rank irangv to all others */
  cs_parall_bcast(root_rank, 3, CS_REAL_TYPE, xyz);

  bft_printf("%s: On rank %d\n"
             "          Velocity U in first cell = %14.5e\n"
             "          Velocity V in first cell = %14.5e\n"
             "          Velocity W in first cell = %14.5e\n\n",
             __func__, root_rank, xyz[0], xyz[1], xyz[2]);
  /* ![example_15] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
