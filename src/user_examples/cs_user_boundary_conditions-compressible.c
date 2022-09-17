/*============================================================================
 * User definition of boundary conditions.
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
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_boundary_conditions.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*=============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  bc_type  boundary face types
 *
 * The icodcl and rcodcl arrays are pre-initialized based on default
 * and GUI-defined definitions, and may be modified here.
 *
 * For a given variable field f, and a given face "face_id", these arrays
 * may be used as follows:
 *
 * - Boundary condition type code given at:
 *   f->bc_coeffs->icodcl[face_id]
 *
 * - Dirichlet value defined at:
 *   f->bc_coeffs->rcodcl1[face_id]
 *
 * - Interior exchange coefficient (infinite if no exchange) at:
 *   f->bc_coeffs->rcodcl2[face_id]
 *
 * - Flux density defined at:
 *   f->bc_coeffs->rcodcl3[face_id]
 *
 * For vector or tensor fields, these arrays are not interleaved,
 * so for a given face "face_id" and field component "comp_id", acess
 * is as follows (where n_b_faces is domain->mesh->n_b_faces):
 *
 *   f->bc_coeffs->rcodcl1[n_b_faces*comp_id + face_id]
 *   f->bc_coeffs->rcodcl2[n_b_faces*comp_id + face_id]
 *   f->bc_coeffs->rcodcl3[n_b_faces*comp_id + face_id]
 *
 * Only the icodcl code values from the first component are used in the case
 * of vector or tensor fields, so the icodcl values can be defined as for
 * a scalar.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(cs_domain_t  *domain,
                            int           bc_type[])
{
  /*![loc_var_dec]*/
  /* Initialization */

  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)domain->mesh_quantities->b_face_normal;

  const cs_zone_t  *zn = NULL;
  const int n_fields = cs_field_n_fields();
  const cs_real_t *cpro_rho = CS_F_(rho)->val;
  const int keysca = cs_field_key_id("scalar_id");

  const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  /*==========================================================================
   * Assign boundary conditions to boundary faces here

   * For each subset:
   * - use selection criteria to filter boundary faces of a given subset
   * - loop on faces from a subset
   *   - set the boundary condition for each face
   *==========================================================================*/

  /* Example of inlet/outlet for which everything is known

   * Without assuming the subsonic or supersonic nature of the inlet
   * the user wishes to impose all the characteristics of the flow,
   * a supersonic inlet is a particular case.

   * The turbulence and the user scalars take a zero flux if the
   * velocity is outward. ---*/

  /*! [example_1] */

  zn = cs_boundary_zone_by_name("1 and X <= 1.0");
  for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {

    const cs_lnum_t face_id = zn->elt_ids[ilelt];
    const cs_lnum_t c_id = b_face_cells[face_id];

    bc_type[face_id] = CS_ESICF;

    /* Velocity */
    CS_F_(vel)->bc_coeffs->rcodcl1[face_id] = 5.0;
    CS_F_(vel)->bc_coeffs->rcodcl1[n_b_faces*1 + face_id] = 0;
    CS_F_(vel)->bc_coeffs->rcodcl1[n_b_faces*2 + face_id] = 0;

    /* Pressure, Density, Temperature, Total Specific Energy

     *     Only 2 of the 4 variables are independant,
     *     hence one can impose values for any couple of variables
     *     (except Temperature-Energy) and the two other variables
     *     will be computed automatically */

    /* Choose a couple of variables which values are to be imposed
     * and delete the others (that will be computed with the help of
     * the thermodynamic laws in cs_cf_thermo.c). */

    /* Pressure (in Pa) */
    CS_F_(p)->bc_coeffs->rcodcl1[face_id] = 5.5;

    /* Temperature (in K) */
    CS_F_(t_kelvin)->bc_coeffs->rcodcl1[face_id] = 300.0;

    /* Specific Total Energy (in J/kg) */
    CS_F_(e_tot)->bc_coeffs->rcodcl1[face_id] = 355.e3;

    /* Turbulence */

    cs_real_t uref2 =  0;
    for(int ii = 0; ii < CS_F_(vel)->dim; ii++)
      uref2 += cs_math_pow2(CS_F_(vel)->bc_coeffs->rcodcl1[n_b_faces*ii +face_id]);
    cs_parall_sum(1, CS_REAL_TYPE, &uref2);
    uref2 = cs_math_fmax(uref2, 1.e-12);

    /* Turbulence example computed using equations valid for a pipe.

     * We will be careful to specify a hydraulic diameter adapted
     * to the current inlet.

     * We will also be careful if necessary to use a more precise
     * formula for the dynamic viscosity use in the calculation of
     * the Reynolds number (especially if it is variable, it may be
     * useful to take the law from 'cs_user_physical_properties'.
     * Here, we use by default the 'viscl0" value.
     * Regarding the density, we have access to its value at boundary
     * faces (romb) so this value is the one used here (specifically,
     * it is consistent with the processing in 'cs_user_physical_properties',
     * in case of variable density) */

    /* Hydraulic diameter */
    const cs_real_t dhy = 0.075;

    const cs_real_t rhomoy = cpro_rho[c_id];

    cs_turbulence_bc_inlet_hyd_diam(face_id,
                                    uref2,
                                    dhy,
                                    rhomoy,
                                    viscl0);

    /*  Handle scalars (avoid modifying rho and energy) */
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *fld  = cs_field_by_id(f_id);
      if ((fld == CS_F_(rho)) || (fld == CS_F_(e_tot)))
        continue;
      int sc_id = -1;
      if (fld->type & CS_FIELD_VARIABLE)
        sc_id = cs_field_get_key_int(fld, keysca) - 1;
      if (sc_id < 0)
        continue;
      fld->bc_coeffs->rcodcl1[face_id] = 1.0;
    }
  }
  /*! [example_1] */

  /* Supersonic outlet example

   * All the characteristics are outward,
   * nothing needs to be imposed (only internal values are used
   * to compute the boundary flux).

   * for the turbulence and the scalar, if values of rcodcl are
   * provided here, we impose them as Dirichlet if the mass flux is
   * inward ; otherwise a zero flux is imposed (outward mass flux or
   * rcodcl values given here).
   * Note that for turbulence, RCODCL has to be filled in for all
   * turbulent variable (otherwise a zero flux is imposed). */

  /*! [example_2] */
   zn = cs_boundary_zone_by_name("2");

   for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {
     const cs_lnum_t face_id = zn->elt_ids[ilelt];
     bc_type[face_id] = CS_SSPCF;
   }
   /*! [example_2] */

   /* Example of subsonic inlet (total pressure, total enthalpy) */

   /*! [example_3] */
   zn = cs_boundary_zone_by_name("3");

   for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {

     const cs_lnum_t face_id = zn->elt_ids[ilelt];

     bc_type[face_id] = CS_EPHCF;

     /* Total pressure (in Pa) */
     CS_F_(p)->bc_coeffs->rcodcl1[face_id] = 1e5;

     /* Total energy */
     CS_F_(e_tot)->bc_coeffs->rcodcl1[face_id] = 294465.0;

     /* Direction of the velocity: normal to inlet faces */
     for(int ii = 0; ii < CS_F_(vel)->dim; ii++)
       CS_F_(vel)->bc_coeffs->rcodcl1[n_b_faces*ii +face_id]
         = -b_face_normal[face_id][ii];

     /* Turbulence (no turbulence)*/

     /*  Handle scalars
         (do not loop on nscal to avoid modifying rho and energy) */
     for (int f_id = 0; f_id < n_fields; f_id++) {
       cs_field_t *fld  = cs_field_by_id(f_id);
       if ((fld == CS_F_(rho)) || (fld == CS_F_(e_tot)))
         continue;
       int sc_id = -1;
       if (fld->type & CS_FIELD_VARIABLE)
         sc_id = cs_field_get_key_int(fld, keysca) - 1;
       if (sc_id < 0)
         continue;
       fld->bc_coeffs->rcodcl1[face_id] = 1.0;
     }

   }
   /*! [example_3] */

   /* Subsonic outlet example */

   /*! [example_4] */
   zn = cs_boundary_zone_by_name("4");

   for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {
     const cs_lnum_t face_id = zn->elt_ids[ilelt];
     bc_type[face_id] = CS_SOPCF;

     /* Pressure (in Pa) */
     CS_F_(p)->bc_coeffs->rcodcl1[face_id] = 5.e5;
   }
   /*! [example_4] */

   /* Wall example */

   /*! [example_5] */
   zn = cs_boundary_zone_by_name("5");

   for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {
     const cs_lnum_t face_id = zn->elt_ids[ilelt];
     bc_type[face_id] = CS_SMOOTHWALL;

     /* Sliding wall */

     /* By default, the wall does not slide.
      * If the wall slides, define the nonzero components of its velocity.
      * The velocity will be projected in a plane tangent to the wall.
      * In the following example, we prescribe Ux = 1. */

     CS_F_(vel)->bc_coeffs->rcodcl1[face_id] = 1.0;

     /* Prescribed temperature
      *
      * By default, the wall is adiabatic.
      * If the wall has a prescribed temperature, indicate it by setting
      * icodcl = 5 and define a value in Kelvin in rcodcl1
      * In the following example, we prescribe T = 293.15 K
      * (example activated if needed=1) */

     CS_F_(t_kelvin)->bc_coeffs->icodcl[face_id] = 5;
     CS_F_(t_kelvin)->bc_coeffs->rcodcl1[face_id] = 20. + 273.15;

     /* Prescribed flux

      * By default, the wall is adiabatic.
      * If the wall has a prescribed flux, indicate it by setting
      * icodcl = 3 and define the value in Watt/m2 in rcodcl3
      * In the following example, we prescribe a flux of 1000 W/m2
      * - a midday in the summer - (example is activated if needed=1) */

     CS_F_(t_kelvin)->bc_coeffs->icodcl[face_id] = 3;
     CS_F_(t_kelvin)->bc_coeffs->rcodcl3[face_id] = 1000.0;
   }
   /*! [example_5] */

   /*![loc_var_dec]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
