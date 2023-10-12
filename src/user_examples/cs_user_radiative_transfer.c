/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <string.h>
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

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_user_radiative_transfer.c
 *
 * \brief User function for input of radiative transfer parameters:
 *        absorption coefficient and net radiation flux.
 *
 *  See \ref cs_user_radiative_transfer for examples.
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Absorption coefficient for radiative module
 *
 * It is necessary to define the value of the fluid's absorption coefficient Ck.
 *
 * This value is defined automatically for specific physical models, such
 * as gas and coal combustion, so this function is not used by these models.
 *
 * For a transparent medium, the coefficient should be set to 0.
 *
 * In the case of the P-1 model, we check that the optical length is at
 * least of the order of 1.
 *
 * \param[in]   bc_type  boundary face types
 * \param[out]  ck       medium's absorption coefficient (zero if transparent)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_rad_transfer_absorption(const int  bc_type[],
                                cs_real_t  ck[])
{
  /* Note: ck argument could be removed, access using field ids */
  {
    cs_real_t *cpro_cak0 = CS_FI_(rad_cak, 0)->val;
    assert(cpro_cak0 == ck);
  }

  /*< [abso_coeff] */
  /*
   * Absorption coefficient of the medium (m-1).
   */

  for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++)
    ck[cell_id] = 0.;

  /*< [abso_coeff]*/

  /*< [abso_coeff_meteo] */

  {
    /* Read file the first time step */

    if (cs_glob_time_step->nt_cur <= cs_glob_time_step->nt_prev + 1) {

      /* Read from file */

      int n_v_segs = 77;
      cs_real_t zray2[78], k_up[78], k_down[78];

      FILE *f;
      float t;

      f = fopen("zray.dat", "r");
      for (int i = 0; i < n_v_segs+1; i++) {
        fscanf(f, "%g", &t);
        zray2[i] = t;
      }
      fclose(f);
      f = fopen("kup.dat", "r");
      for (int i = 0; i < n_v_segs+1; i++) {
        fscanf(f, "%g", &t);
        k_up[i] = t;
      }
      fclose(f);
      f = fopen("kdown.dat", "r");
      for (int i = 0; i < n_v_segs+1; i++) {
        fscanf(f, "%g", &t);
        k_down[i] = t;
      }
      fclose(f);

      cs_field_t *f_ck_u = cs_field_by_name("rad_absorption_coeff_up");
      cs_field_t *f_ck_d = cs_field_by_name("rad_absorption_coeff_down");

      for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++) {

        cs_real_t z = cs_glob_mesh_quantities->cell_cen[cell_id*3 + 2];

        /* Search for 1D cell */
        int iz = 0;

        cs_real_t tol = 1.0;
        /* Search loop */
        for (; iz < n_v_segs && z > (zray2[iz+1]-tol); iz++);

        f_ck_u->val[cell_id] = fabs(k_up[iz+1]);
        f_ck_d->val[cell_id] = fabs(k_down[iz+1]);

      }
    }

  }
  /*< [abso_coeff_meteo] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the net radiation flux.
 *
 * The density of net radiation flux must be calculated
 * consistently with the boundary conditions of the intensity.
 * The density of net flux is the balance between the radiative
 * emiting part of a boundary face (and not the reflecting one)
 * and the radiative absorbing part.
 *
 * \param[in]   bc_type   boundary face types
 * \param[in]   twall     inside current wall temperature (K)
 * \param[in]   qincid    radiative incident flux  (W/m2)
 * \param[in]   xlam      conductivity (W/m/K)
 * \param[in]   epa       thickness (m)
 * \param[in]   eps       emissivity (>0)
 * \param[in]   ck        absorption coefficient
 * \param[out]  net_flux  net flux (W/m2)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_rad_transfer_net_flux(const int        bc_type[],
                              const cs_real_t  twall[],
                              const cs_real_t  qincid[],
                              const cs_real_t  xlam[],
                              const cs_real_t  epa[],
                              const cs_real_t  eps[],
                              const cs_real_t  ck[],
                              cs_real_t        net_flux[])
{
  /*< [loc_var_dec_2]*/
  const cs_real_t stephn = cs_physical_constants_stephan;
  /*< [loc_var_dec_2]*/

  /*< [init]*/
  /* Initializations */
  /* get BC coeffs for radiance (explicit part of the first band) */
  cs_real_t *coefap = CS_FI_(radiance, 0)->bc_coeffs->a;

  /*< [init]*/

  /*< [net_flux]*/
  /* Net flux dendity for the boundary faces
   * The provided examples are sufficient in most of cases.*/

  /* If the boundary conditions given above have been modified
   *   it is necessary to change the way in which density is calculated from
   *   the net radiative flux consistently.*/

  /* The rule is:
   *   the density of net flux is a balance between the emitting energy from a
   *   boundary face (and not the reflecting energy) and the absorbing radiative
   *   energy. Therefore if a wall heats the fluid by radiative transfer, the
   *   net flux is negative */

  for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_b_faces; ifac++) {

    /* Wall faces */
    if (   bc_type[ifac] == CS_SMOOTHWALL
        || bc_type[ifac] == CS_ROUGHWALL)
      net_flux[ifac] = eps[ifac] * (qincid[ifac] - stephn * pow(twall[ifac], 4));

    /* Symmetry   */
    else if (bc_type[ifac] == CS_SYMMETRY)
      net_flux[ifac] = 0.0;

    /* Inlet/Outlet    */
    else if (   bc_type[ifac] == CS_INLET
             || bc_type[ifac] == CS_CONVECTIVE_INLET
             || bc_type[ifac] == CS_OUTLET
             || bc_type[ifac] == CS_FREE_INLET) {
      if (cs_glob_rad_transfer_params->type == CS_RAD_TRANSFER_DOM)
        net_flux[ifac] = qincid[ifac] - cs_math_pi * coefap[ifac];
      else if (cs_glob_rad_transfer_params->type == CS_RAD_TRANSFER_P1)
        net_flux[ifac] = 0.0;
    }
    /* Stop if there are forgotten faces   */
    else
      bft_error
        (__FILE__, __LINE__, 0,
         "In %s:\n"
         "  non-handled boundary faces for net flux calculation\n\n"
         "  Last face: %10ld; zone = %d; nature = %d\n",
         __func__,
         (long)ifac,
         cs_boundary_zone_face_zone_id()[ifac],
         bc_type[ifac]);

  }
  /*< [net_flux]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
