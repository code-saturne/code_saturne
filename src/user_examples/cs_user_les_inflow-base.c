
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
#include <string.h>

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
 * \file cs_user_les_inflow-base.c
 *
 * \brief Generation of synthetic turbulence at LES inlets initialization.
 *
 * See \ref les_inflow for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define parameters of synthetic turbulence at LES inflow.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_les_inflow_define(void)
{
  /*! [set_restart] */
  cs_les_inflow_set_restart(false,   /* allow_read */
                            false);  /* allow_write */
  /*! [set_restart] */

  /* First synthetic turbulence inlet: the Batten Method is used
     for boundary faces of zone "INLET_1" */

  /*! [init_1] */
  {
    /* Number of turbulence structures */
    int n_entities = 50;

    /* Velocity, turbulent kinetic energy and dissipation rate */
    cs_real_t vel_r[3] = {18.0, 0, 0};
    cs_real_t k_r = 4.0;
    cs_real_t eps_r = 4.0;

    cs_les_inflow_add_inlet(CS_INFLOW_BATTEN,
                            false,
                            cs_boundary_zone_by_name("INLET_1"),
                            n_entities,
                            0,  /* verbosity */
                            vel_r,
                            k_r,
                            eps_r);
  }
  /*! [init_1] */

  /* Second synthetic turbulence inlet: the Synthetic Eddy Method
     Batten Method is used for boundary faces of zone "INLET_2" */

  /*! [init_2] */
  {
    /* Number of turbulence structures */
    int n_entities = 50;

    /* Velocity, turbulent kinetic energy and dissipation rate */
    cs_real_t vel_r[3] = {12.0, 0, 0};
    cs_real_t k_r = 3.0;
    cs_real_t eps_r = 3.0;

    cs_les_inflow_add_inlet(CS_INFLOW_SEM,
                            false,  /* volume mode */
                            cs_boundary_zone_by_name("INLET_2"),
                            n_entities,
                            0,  /* verbosity */
                            vel_r,
                            k_r,
                            eps_r);
  }
  /*! [init_2] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update of the characteristics of a given synthetic turbulence inlet.
 *
 * \param[in]   zone       pointer to associated boundary zone
 * \param[out]  vel_r      reference mean velocity
 * \param[out]  k_r        reference turbulent kinetic energy
 * \param[out]  eps_r      reference turbulent dissipation
 */
/*----------------------------------------------------------------------------*/

void
cs_user_les_inflow_update(const cs_zone_t  *zone,
                          cs_real_t         vel_r[3],
                          cs_real_t        *k_r,
                          cs_real_t        *eps_r)
{
  /*! [update_1] */
  if (strcmp(zone->name, "INLET_1") == 0) {
    /* Velocity, turbulent kinetic energy and dissipation rate */
    vel_r[0] = 19.0;
    vel_r[1] = 0.0;
    vel_r[2] = 0.0;
    *k_r = 5.0;
    *eps_r = 5.0;
  }
  /*! [update_1] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Definition of mean velocity, Reynolds stresses and dissipation rate
 *        for each boundary face of the given synthetic turbulence inlet.
 *
 * Rij components are ordered as usual: XX, YY, ZZ, XY, YZ, XZ
 *
 * Arrays are pre-initialized before this function is called
 * (see \ref cs_user_les_inflow_define).
 *
 * \param[in]       zone    pointer to associated boundary zone
 * \param[in, out]  vel_l   velocity a zone faces
 * \param[in, out]  rij_l   Reynolds stresses at zone faces
 * \param[in, out]  eps_l   reference turbulent dissipation
 */
/*----------------------------------------------------------------------------*/

void
cs_user_les_inflow_advanced(const cs_zone_t  *zone,
                            cs_real_3_t       vel_l[],
                            cs_real_6_t       rij_l[],
                            cs_real_t         eps_l[])
{
  /*  Example 1:
   *   - mean velocity, Reynolds stresses an dissipation are deduced
   *     from a wall law for the first synthetic turbulence inlet,
   *   - no modification of the statistics of the flow is provided for the
   *     second synthetic turbulence inlet */

  /*! [example_1] */
  if (strcmp(zone->name, "INLET_1") == 0) {

    const cs_lnum_t n_faces = zone->n_elts;
    const cs_lnum_t *face_ids = zone->elt_ids;

    const cs_real_t d2_s3 = 2.0/3.0;

    const cs_mesh_t   *m = cs_glob_mesh;
    const cs_real_3_t *b_face_cog
      = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;

    const cs_real_t *cpro_mu = CS_F_(mu)->val;

    const cs_real_t uref = cs_glob_turb_ref_values->uref;

    /* Approximation of the friction velocity */
    const cs_real_t utau = uref/20.0;

    /* Reference length scale */
    const cs_real_t h_ref = 1.0;

    const cs_real_t kappa = cs_turb_xkappa;

    for (cs_lnum_t i = 0; i < n_faces; i++) {

      const cs_lnum_t f_id = face_ids[i];
      const cs_lnum_t c_id = m->b_face_cells[f_id];

      const cs_real_t x_mu = cpro_mu[c_id];
      const cs_real_t r_fro = utau * h_ref / x_mu;

      /* Dimensionless wall distance */
      const cs_real_t yy = 1.0 - b_face_cog[f_id][1];
      const cs_real_t yplus = yy / h_ref * r_fro;

      /* Reichart laws (dimensionless) */
      const cs_real_t uplus = log(1. + 0.4*yplus)/kappa
                             + 7.8*(1. - exp(-yplus/11.0)
                                    - yplus/11.0*exp(-0.33*yplus));
      const cs_real_t kplus =   0.07*yplus*yplus*exp(-yplus/8.0)
                              + (1.0 - exp(-yplus/20.0)) *4.5
                              / (1. + 4.0*yplus/r_fro);
      const cs_real_t eplus =  (1.0/kappa)
                              / pow((  cs_math_pow4(yplus)
                                     + cs_math_pow4(15.)), 0.25);

      /* Arrays are filled with dimensionnal stats */
      vel_l[i][0] = uplus * utau;
      vel_l[i][1] = 0.0;
      vel_l[i][2] = 0.0;

      rij_l[i][0] = d2_s3 * kplus * cs_math_pow2(utau);
      rij_l[i][1] = d2_s3 * kplus * cs_math_pow2(utau);
      rij_l[i][2] = d2_s3 * kplus * cs_math_pow2(utau);
      rij_l[i][3] = 0;
      rij_l[i][4] = 0;
      rij_l[i][5] = 0;

      eps_l[i] = eplus * cs_math_pow4(utau) / x_mu;
    }
  }
  /*! [example_1] */

  /*  Example 2:
   *  - Reynolds stresses and dissipation at the inlet are computed
   *    using the turbulence intensity and standard laws for
   *    a circular pipe for the first synthetic turbulence inlet,
   *   - no modification of the statistics of the flow is provided for the
   *     second synthetic turbulence inlet */


  /*! [example_2] */
  if (strcmp(zone->name, "INLET_1") == 0) {

    const cs_lnum_t n_faces = zone->n_elts;

    const cs_real_t d2_s3 = 2.0/3.0;

    for (cs_lnum_t i = 0; i < n_faces; i++) {

      vel_l[i][0] = 1.1;
      vel_l[i][1] = 1.1;
      vel_l[i][2] = 1.1;

      cs_real_t uref2 = cs_math_3_square_norm(vel_l[i]);
      uref2 = cs_math_fmax(uref2, 1e-12);

      /* Hydraulic diameter */
      const cs_real_t x_hd = 0.075;

      /* Turbulence intensity */
      const cs_real_t x_ti = 0.02;

      cs_real_t x_k = 1e-12;
      cs_real_t x_eps = 1e-12;

      cs_turbulence_bc_ke_turb_intensity(uref2,
                                         x_ti,
                                         x_hd,
                                         &x_k,
                                         &x_eps);

      rij_l[i][0] = d2_s3 * x_k;
      rij_l[i][1] = d2_s3 * x_k;
      rij_l[i][2] = d2_s3 * x_k;
      rij_l[i][3] = 0;
      rij_l[i][4] = 0;
      rij_l[i][5] = 0;

      eps_l[i] = x_eps;
    }
  }
  /*! [example_2] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
