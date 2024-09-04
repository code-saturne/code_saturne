/*============================================================================
 * Compute the flux at the boundary.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"

#include "cs_cf_boundary_conditions.h"
#include "cs_cf_thermo.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cf_boundary_flux.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cf_boundary_flux.c

  Compute the flux at the boundary.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes the analytical flux at the boundary for Euler and Energy
 *
 * The Euler equations used to compute the flux are:
 * \f{eqnarray*}{
 *    \der{\rho}{t} + \divs \left(\rho\vect{u}\right) &=&0
 * \\ \der{\rho \vect{u}}{t} + \divv \left(\vect{u}\otimes\rho\vect{u}\right)&=&0
 * \\ \der{E}{t} + \divs \left(\rho\vect{u} E\right) &=&0
 * \f}
 *
 * \param[in]  f_id        face id
 * \param[in]  val_ext_en  Dirichlet value for the total energy
 * \param[in]  val_ext_p   Dirichlet value for the pressure
 * \param[in]  val_ext_v   Dirichlet value for the velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_boundary_analytical_flux(const cs_lnum_t    f_id,
                               const cs_real_t   *val_ext_en,
                               const cs_real_t   *val_ext_p,
                               const cs_real_3_t *val_ext_v)
{
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t *b_f_face_surf = fvq->b_f_face_surf;
  const cs_real_3_t *b_face_u_normal
    = (const cs_real_3_t *)fvq->b_face_u_normal;

  cs_real_3_t *cofacv = (cs_real_3_t *)CS_F_(vel)->bc_coeffs->ac;
  cs_real_t *coface = CS_F_(e_tot)->bc_coeffs->ac;
  const cs_real_t *rho_b = CS_F_(rho_b)->val;

  /* Compute values needed for analytical flux
     ----------------------------------------- */

  const cs_real_t und = cs_math_3_dot_product(val_ext_v[f_id],
                                              b_face_u_normal[f_id]);

  const cs_real_t rund = rho_b[f_id] * und;

  /* Convective analytical flux
     -------------------------- */

  /* Tag the faces where an analytical flux is computed
     The tag will be used in bilsc2 to retrieve the faces
     where an analytical flux has to be imposed */

  int *icvfli = cs_cf_boundary_conditions_get_icvfli();
  icvfli[f_id] = 1;

  /* Momentum flux (the centered pressure contribution is
                    directly taken into account
                    in the pressure BC) */

  cofacv[f_id][0] = b_f_face_surf[f_id] * rund * val_ext_v[f_id][0];
  cofacv[f_id][1] = b_f_face_surf[f_id] * rund * val_ext_v[f_id][1];
  cofacv[f_id][2] = b_f_face_surf[f_id] * rund * val_ext_v[f_id][2];

  /* Total energy flux */
  coface[f_id] =   b_f_face_surf[f_id]
                 * (rund * val_ext_en[f_id] + und * val_ext_p[f_id]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes the Rusanov flux at the boundary for Euler and Energy
 *
 * d rho   /dt + div rho u             = 0
 * d rho u /dt + div rho u u + grad  P = 0
 * d E     /dt + div rho u E + div u P = 0
 *
 * \param[in]       f_id        face id
 * \param[in]       val_ext_en  Dirichlet value for the total energy
 * \param[in, out]  val_ext_p   Dirichlet value for the pressure
 * \param[in, out]  val_ext_v   Dirichlet value for the velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_boundary_rusanov(const cs_lnum_t  f_id,
                       const cs_real_t *val_ext_en,
                       cs_real_t       *val_ext_p,
                       cs_real_3_t     *val_ext_v)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_real_t *restrict b_f_face_surf = fvq->b_f_face_surf;
  const cs_real_3_t *b_face_u_normal = (const cs_real_3_t *)fvq->b_face_u_normal;

  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  const int icp = fluid_props->icp;
  const int icv = fluid_props->icv;

  /* Initialization
     -------------- */

  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_e_tot = CS_F_(e_tot);
  cs_real_3_t *vel = (cs_real_3_t *)f_vel->val;
  cs_real_t *pres = CS_F_(p)->val;
  cs_real_t *energy = f_e_tot->val;

  cs_real_3_t *cofacv = (cs_real_3_t *)f_vel->bc_coeffs->ac;
  cs_real_t *coface = f_e_tot->bc_coeffs->ac;

  cs_real_t *rho = CS_F_(rho)->val;
  cs_real_t *rho_b = CS_F_(rho_b)->val;

  const cs_lnum_t c_id = b_face_cells[f_id];

  /* Initialize local specific heat values and handle uniform cases */

  cs_real_t c_cp = (icp >= 0) ? CS_F_(cp)->val[c_id] : 0.;
  cs_real_t c_cv = (icv >= 0) ? CS_F_(cv)->val[c_id] : 0.;

  /* Compute values needed for Rusanov scheme
     ---------------------------------------- */

  const cs_real_t *n = b_face_u_normal[f_id];

  const cs_real_t b_vel_n = cs_math_3_dot_product(val_ext_v[f_id], n);
  const cs_real_t c_vel_n = cs_math_3_dot_product(vel[c_id], n);
  const cs_real_t b_mass_flux = rho_b[f_id] * b_vel_n;
  const cs_real_t c_mass_flux = rho[c_id] * c_vel_n;

  cs_real_t b_c2;
  cs_cf_thermo_c_square(&c_cp,
                        &c_cv,
                        &val_ext_p[f_id],
                        &rho_b[f_id],
                        NULL,
                        NULL,
                        NULL,
                        &b_c2,
                        1);

  cs_real_t c_c2;
  cs_cf_thermo_c_square(&c_cp,
                        &c_cv,
                        &pres[c_id],
                        &rho[c_id],
                        NULL,
                        NULL,
                        NULL,
                        &c_c2,
                        1);

  const cs_real_t b_c = sqrt(b_c2);
  const cs_real_t c_c = sqrt(c_c2);
  const cs_real_t rrus = cs_math_fmax(cs_math_fabs(b_vel_n) + b_c,
                                      cs_math_fabs(c_vel_n) + c_c);

  /* boundary mass flux computed with Rusanov scheme */

  const cs_real_t r_b_masfl
    = 0.5*(b_mass_flux+c_mass_flux) - 0.5*rrus*(rho_b[f_id]-rho[c_id]);

  /* Rusanov normal velocity (computed using boundary density) */

  const cs_real_t r_b_vel_n = r_b_masfl / rho_b[f_id];

  /* Update velocity boundary condition */

  for (cs_lnum_t i = 0; i < 3; i++)
    val_ext_v[f_id][i] += (r_b_vel_n - b_vel_n) * n[i];

  /* Convective Rusanov flux
     ----------------------- */

  /* Tag the faces where a Rusanov flux is computed
     The tag will be used in bilsc2 to retrieve the faces
     where a Rusanov flux has to be imposed */

   int *icvfli = cs_cf_boundary_conditions_get_icvfli();
   icvfli[f_id] = 1;

   /* Momentum flux (the centered pressure contribution is directly
      taken into account in the pressure BC) */

   for (cs_lnum_t i = 0; i < 3; i++) {
     const cs_real_t sum_flux
       = b_mass_flux * val_ext_v[f_id][i] + c_mass_flux * vel[c_id][i];
     const cs_real_t diff_mom
       = rho_b[f_id] * val_ext_v[f_id][i] - rho[c_id] * vel[c_id][i];

     cofacv[f_id][i] = 0.5 * b_f_face_surf[f_id] *(sum_flux - rrus * diff_mom);
   }

   /* BC for the pressure gradient in the momentum balance */
   val_ext_p[f_id] = 0.5 * (val_ext_p[f_id] + pres[c_id]);

   /* Total energy flux */

   const cs_real_t sum_e
     = b_mass_flux * val_ext_en[f_id] + c_mass_flux * energy[c_id];

   const cs_real_t sum_p
     = b_vel_n * val_ext_p[f_id] + c_vel_n * pres[c_id];

   const cs_real_t diff_e
       = rho_b[f_id] * val_ext_en[f_id] - rho[c_id] * energy[c_id];

   coface[f_id] = 0.5 * b_f_face_surf[f_id] * (sum_e + sum_p - rrus * diff_e);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
