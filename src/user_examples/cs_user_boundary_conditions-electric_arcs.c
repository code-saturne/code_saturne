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
  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;

  const cs_real_3_t *cdgfbo
    = (const cs_real_3_t *)domain->mesh_quantities->b_face_cog;

  /*  Initialization */

  const int ielcor = cs_glob_elec_option->ielcor;

  cs_field_t *potr = cs_field_by_name("elec_pot_r");;
  cs_field_t *poti = NULL;
  cs_field_t *potva = NULL;
  if (cs_glob_physical_model_flag[CS_JOULE_EFFECT] > 1) {
    poti = cs_field_by_name("elec_pot_i");
  }
  if (cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > 1) {
    potva = cs_field_by_name("vec_potential");
  }

  /* Assign boundary conditions to boundary faces here

   * For each subset:
   * - use selection criteria to filter boundary faces of a given subset
   * - loop on faces from a subset
   *   - set the boundary condition for each face

   * For boundary faces of zone 1 assign an inlet
   *   and assign a cathode for "electric" variables. */

  /*! [example_1] */
  {
    const cs_zone_t *zone = cs_boundary_zone_by_name("1");
    cs_field_t *th_f = cs_thermal_model_field();
    const cs_real_t *bfpro_rom =  CS_F_(rho_b)->val;

    cs_real_t *rcodcl1_vel = CS_F_(vel)->bc_coeffs->rcodcl1;

    for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

      const cs_lnum_t face_id = zone->elt_ids[ilelt];

      bc_type[face_id] = CS_INLET;

      /* Velocity */
      rcodcl1_vel[n_b_faces*0 + face_id] = 0;
      rcodcl1_vel[n_b_faces*1 + face_id] = 0;
      rcodcl1_vel[n_b_faces*2 + face_id] = 0;

      /* Turbulence */

      if (   (cs_glob_turb_model->itytur == 2)
          || (cs_glob_turb_model->itytur == 3)
          || (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI)
          || (cs_glob_turb_model->iturb == CS_TURB_K_OMEGA)
          || (cs_glob_turb_model->iturb == CS_TURB_SPALART_ALLMARAS)) {

        cs_real_t uref2 =  0;
        for(int ii = 0; ii < 3; ii++)
          uref2 += cs_math_pow2(rcodcl1_vel[n_b_faces*ii +face_id]);
        uref2 = cs_math_fmax(uref2, 1.e-12);

        /* Turbulence example computed using equations valid for a pipe.

         * We will be careful to specify a hydraulic diameter adapted
         * to the current inlet.

         * We will also be careful if necessary to use a more precise
         * formula for the dynamic viscosity use in the calculation of
         * the Reynolds number (especially if it is variable, it may be
         * useful to take the law from 'cs_user_physical_properties'.
         * Here, we use by default the 'viscl0" value.

         * Regarding the density, we have acess to its value at boundary
         * faces (romb) so this value is the one used here (specifically,
         * it is consistent with the processing in 'cs_user_physical_properties',
         * in case of variable density)

         * Hydraulic diameter */
        const cs_real_t dhy = 0.075;

        const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

        const cs_real_t rhomoy = bfpro_rom[face_id];

        cs_turbulence_bc_inlet_hyd_diam(face_id,
                                        uref2,
                                        dhy,
                                        rhomoy,
                                        viscl0);

        /* Handle Scalars

         * Enthalpy in J/kg (iscalt)
         * On this example we impose the value of the enthalpy
         * the arbitrary value of 1.d6 corresponds to a temperature
         * of 2200 Kelvin for argon at atmospheric pressure (see dp_ELE) */

        th_f->bc_coeffs->icodcl[face_id] = 1;
        th_f->bc_coeffs->rcodcl1[face_id] = 1e6;

        /* Specific model for Electric arc:
         * Vector Potential: Zero flux by default beacuse we don't a lot about
         * vector potential (what we know, is that A is equal to zero at
         * and infinite distance)

         * All the boundary conditions for A are zero flux, except on some
         * selected faces where we need to impose a value in order to have
         * a stable calculation (well defined problem).
         * These faces are selected where we are sure that the electrical current
         * density remains very low generally far from the center of the electric
         * arc and from the electrodes (see above) */

        if (potva != NULL) {
          potva->bc_coeffs->icodcl[face_id] = 3;
          potva->bc_coeffs->rcodcl3[face_id] = 0;
          potva->bc_coeffs->icodcl[n_b_faces + face_id] = 3;
          potva->bc_coeffs->rcodcl3[n_b_faces + face_id] = 0;
          potva->bc_coeffs->icodcl[n_b_faces*2 + face_id] = 3;
          potva->bc_coeffs->rcodcl3[n_b_faces*2 + face_id] = 0;
        }
      }
    }
    /*! [example_1] */
  }

  /* For boundary faces of zone 2 assign a free outlet
   * =================================================
   * and example of electrode for Joule Effect by direct conduction.
   * ============================================================== */

  /*! [example_2] */
  {
    const cs_zone_t *zone = cs_boundary_zone_by_name("2");
    for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

      const cs_lnum_t face_id = zone->elt_ids[ilelt];

      bc_type[face_id] = CS_OUTLET;

      /* Handle Scalars

       * Enthalpy in J/kg  (By default zero flux with CS_OUTLET)
       *    Nothing to do

       * Mass fraction of the (n-1) gas mixture components
       * (Zero flux by defaut with ISOLIB)
       *    Nothing to do

       * Specific model for Joule Effect by direct conduction:

       * If you want to make a simulation with an imposed Power PUISIM
       * (you want to get PUISIM imposed in cs_user_parameters.c and
       * PUISIM = Amp x Volt) you need to impose IELCOR=1 in cs_user_parameters.c
       * The boundary conditions will be scaled by COEJOU coefficient
       * for example the electrical potential will be multiplied bu COEJOU
       * (both real and imaginary part of the electrical potential if needed)

       * COEJOU is automatically defined in order that the calculated dissipated
       * power by Joule effect (both real and imaginary part if needed) is equal
       * to PUISIM

       * At the beginning of the calculation, COEJOU ie equal to 1;
       * COEJOU is writing and reading in the result files.

       * If you don't want to calculate with by scaling,
       * you can impose directly the value. */

      if (cs_glob_physical_model_flag[CS_JOULE_EFFECT] > -1) {
        potr->bc_coeffs->icodcl[face_id] = 1;
        if (ielcor == 1) {
          const cs_real_t coejou = cs_glob_elec_option->coejou;
          potr->bc_coeffs->rcodcl1[face_id] = 500*coejou;
        }
        else {
          potr->bc_coeffs->rcodcl1[face_id] = 500.;
        }
      }

      if ( cs_glob_physical_model_flag[CS_JOULE_EFFECT] > 1) {
        poti->bc_coeffs->icodcl[face_id] = 1;
        if (ielcor == 1) {
          const cs_real_t coejou = cs_glob_elec_option->coejou;
          poti->bc_coeffs->rcodcl1[face_id] = sqrt(3)*500*coejou;
        }
        else {
          poti->bc_coeffs->rcodcl1[face_id] = sqrt(3)*500.;
        }
      }

    }
    /*! [example_2] */
  }

  /* For boundary faces of zone 2 assign a free outlet
   * =================================================================
   * and example of anode for electric arc.
   * ===================================== */

  /*! [example_3] */
  {
    const cs_zone_t *zone = cs_boundary_zone_by_name("2");
    for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

      const cs_lnum_t face_id = zone->elt_ids[ilelt];

      bc_type[face_id] = CS_OUTLET;

      /* Handle scalars

       * Enthalpy in J/kg  (Zero flux by default with ISOLIB)
       *    Nothing to do

       * Real component of the electrical potential */

      /*  For electric arc model,
       * ======================
       * We generally calculate the "electric variables" assuming that the total
       * intensity of the electrical current is imposed (COUIMP is the value of
       * the imposed total current).

       * In that case, you need to impose IELCOR=1 in cs_user_parameters.c
       * The "electrical variables" will be scaled by COEPOT coefficient :
       * for example the electrical potential will be multiplied by COEPOT,
       * Joule effect will be multipied by COEPOT * COEPOT and so on (see
       * cs_user_electric_scaling.c).

       * COEJOU is defined in cs_user_electric_scaling.c: different possibilities
       * are described in cs_user_electric_scaling.c, depending on the different
       * physics you want to simulate (scaling from current, from power, special
       * model for restriking ...).

       * The variable DPOT is defined: it corresponds to the electrical potential
       * difference between the electrodes (Anode potential - cathode Potential).
       * DPOT is calculated in cs_user_electric_scaling.c. DPOT is saved at each
       * time step, and for a following calculation

       * DPOT is the value of the boundary condition on anode assuming that
       * the cathode potential is equel to zero.

       * It is also possible to fix the value of the potential on the anode.
       * (for example, 1000 Volts).*/

      potr->bc_coeffs->icodcl[face_id] = 1;

      if ((cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > -1) && (ielcor == 1))
        potr->bc_coeffs->rcodcl1[face_id] = cs_glob_elec_option->pot_diff;
      else
        potr->bc_coeffs->rcodcl1[face_id] = 1000.;

      /* Mass fraction of the (n-1) gas mixture components
       *  zero flux by default with CS_OUTLET
       *  nothing to do

       * vector Potential
       *  zero flux by default with CS_OUTLET
       *  nothing to do */
    }
    /*! [example_3] */
  }

  /* For boundary faces of zone 2 assign a wall
   * ==========================================
   * and example of potential vector Dirichlet condition */

  /*! [example_4] */
  {
    const cs_real_3_t *cvara_potva
      = (const cs_real_3_t *)cs_field_by_name("vec_potential")->val_pre;

    const cs_zone_t *zone = cs_boundary_zone_by_name("2");
    for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

      const cs_lnum_t face_id = zone->elt_ids[ilelt];

      bc_type[face_id] = CS_SMOOTHWALL;

      /* Wall: zero flow (zero flux for pressure)
       *       friction for velocities (+ turbulent variables)
       *       zero flux for scalars

       * Handle scalars
       * Enthalpy in J/kg  (Zero flux by default)
       *   Nothing to do

       * Real component of the electrical potential
       *   Zero flux by default
       *   Nothing to do

       * Specific model for Electric arc:
       * ===============================

       * Vector potential  A (Ax, Ay, Az)

       * Zero flux by default because we don't a lot about vector potential
       * (what we know, is that A is equal to zero at the infinite)

       * All the boundary conditions for A are zero flux, except on some chosen
       * faces where we need to impose a value in order to have a stable
       * calculation. These faces are chosen where we are sure that the
       * electrical current density remains very low generally far from the
       * center of the electric arc and from the electrodes:

       * On the following example, we choose to impose a "dirichlet" value for
       * the 3 components of A on a small zone of the boundary located near the
       * vertical free outlet of the computation domain.

       * In this example, the electric arc is at the center of the computational
       * domain, located on z axis (near x = 0 and y = 0).
       * The x (1st) and y (the 3rd) coordinates are contained between
       * -2.5 cm nd 2.5 cm:

       *    Ax(t, x,y,z) = Ax(t-dt, x=2.5cm, y=2.5cm, z)
       *    Ay(t, x,y,z) = Ay(t-dt, x=2.5cm, y=2.5cm, z)
       *    Az(t, x,y,z) = Az(t-dt, x=2.5cm, y=2.5cm, z) */

      if (cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > 1) {
        if (   (cdgfbo[face_id][0] <= 249e-2)
            || (cdgfbo[face_id][0] >= 249e-2)
            || (cdgfbo[face_id][2] <= 249e-2)
            || (cdgfbo[face_id][2] >= 249e-2)) {

          const cs_lnum_t c_id = b_face_cells[face_id];
          potva->bc_coeffs->icodcl[face_id] = 1;
          potva->bc_coeffs->rcodcl1[face_id] = cvara_potva[c_id][0];

          potva->bc_coeffs->icodcl[n_b_faces + face_id] = 1;
          potva->bc_coeffs->rcodcl1[n_b_faces + face_id] = cvara_potva[c_id][1];

          potva->bc_coeffs->icodcl[n_b_faces*2 + face_id] = 1;
          potva->bc_coeffs->rcodcl1[n_b_faces*2 + face_id]
            = cvara_potva[c_id][2];

        }
      }
    }
    /*![example_4]*/
  }

  /* For boundary faces of zone 51 assign a wall
   * ============================================
   * and restriking model for electric arc (anode boundaray condition)
   * ================================================================= */

  /*![example_5]*/
  {
    cs_field_t *th_f = cs_thermal_model_field();

    const cs_zone_t *zone = cs_boundary_zone_by_name("51");
    for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

      const cs_lnum_t face_id = zone->elt_ids[ilelt];

      bc_type[face_id] = CS_SMOOTHWALL;

      /* Enthalpy (J/kg) :
       *   imposed heat transfer coefficient */

      th_f->bc_coeffs->icodcl[face_id] = 1;
      th_f->bc_coeffs->rcodcl1[face_id] = 2.e4;
      th_f->bc_coeffs->rcodcl2[face_id] = 1.e5;

      /* Real electrical potential: anode boundary condition;
       * pot_diff calculated in cs_user_electric_scaling.c */

      potr->bc_coeffs->icodcl[face_id] = 1;

      if ((cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > 1) && (ielcor == 1)) {
        if (   (cs_glob_elec_option->irestrike == 1)
            && (cs_glob_time_step->nt_cur <= cs_glob_elec_option->ntdcla+30)) {
          cs_real_t z1 = cs_glob_elec_option->restrike_point[2] - 2.e-4;
          if (z1 <= 0.0)
            z1 = 0.0;
          cs_real_t z2 = cs_glob_elec_option->restrike_point[2] + 2.e-4;
          if (z2 >= 2e-2)
            z2 = 2e-2;

          if ((cdgfbo[face_id][2] >= z1) && (cdgfbo[face_id][2] <= z2)) {
            th_f->bc_coeffs->icodcl[face_id] = 1;
            th_f->bc_coeffs->rcodcl1[face_id] = cs_glob_elec_option->pot_diff;
          }
          else {
            th_f->bc_coeffs->icodcl[face_id] = 3;
            th_f->bc_coeffs->rcodcl3[face_id] = 0;
          }
        }
      }

      /* Vector potential : Zero flux */

      if (potva != NULL) {
        potva->bc_coeffs->icodcl[n_b_faces*0 + face_id] = 3;
        potva->bc_coeffs->icodcl[n_b_faces*1 + face_id] = 3;
        potva->bc_coeffs->icodcl[n_b_faces*2 + face_id] = 3;
        potva->bc_coeffs->rcodcl3[n_b_faces*0 + face_id] = 0;
        potva->bc_coeffs->rcodcl3[n_b_faces*1 + face_id] = 0;
        potva->bc_coeffs->rcodcl3[n_b_faces*2 + face_id] = 0;
      }
    }
    /*![example_5]*/
  }

  /* For boundary faces of zone 4 assign a symetry
   * ============================================= */

  /*![example_6]*/
  {
    const cs_zone_t *zone = cs_boundary_zone_by_name("4");
    for (cs_lnum_t ilelt = 0; ilelt < zone->n_elts; ilelt++) {

      const cs_lnum_t face_id = zone->elt_ids[ilelt];

      bc_type[face_id] = CS_SYMMETRY;
      /* For all scalars, by default a zero flux condition is assumed
       * (for potentials also)

       * In Joule effect direct conduction,
       * we can use an anti-symetry condition for the imaginary component of the
       * electrical potential depending on the electrode configuration: */

      if (poti != NULL) {
        poti->bc_coeffs->icodcl[face_id] = 1;
        poti->bc_coeffs->rcodcl1[face_id] = 0.;
      }
    }
    /*![example_6]*/
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
