/*============================================================================
 * Compressible flow boundary conditions.
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

/*----------------------------------------------------------------------------
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_array.h"
#include "bft_mem.h"
#include "cs_cf_model.h"
#include "cs_cf_boundary_flux.h"
#include "cs_cf_thermo.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cf_boundary_conditions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Prototypes for Fortran functions and variables.
 *============================================================================*/

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cf_boundary_conditions.c
        Compressible flow boundary conditions.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Imposed thermal flux indicator at the boundary
  (some boundary contributions of the total energy equation
  have to be cancelled) */
static int *_ifbet = nullptr;

/*! boundary convection flux indicator of a Rusanov or an analytical flux
  (some boundary contributions of the momentum eq. have to be cancelled) */
static int *_icvfli = nullptr;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize coal model.
 *
 * \pram[in, out]  cm  pointer to coal model pointer to destroy.
 */
/*----------------------------------------------------------------------------*/

static void
_cf_boundary_conditions_finalize(void)
{
  BFT_FREE(_ifbet);
  BFT_FREE(_icvfli);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic boundary condition for compressible flows
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_boundary_conditions(int  bc_type[])
{
  /* Initializations
     --------------- */

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_3_t *b_face_cog
    = reinterpret_cast<const cs_real_3_t *>(fvq->b_face_cog);

  const cs_real_3_t *cell_cen
    = reinterpret_cast<const cs_real_3_t *>(fvq->cell_cen);

  const cs_real_t *b_dist = fvq->b_dist;

  const cs_fluid_properties_t  *fluid_props = cs_glob_fluid_properties;
  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  const cs_real_t cv0 = fluid_props->cv0;

  const int icfgrp = cs_glob_cf_model->icfgrp;

  // Threshold to detect whether some BC values were set or are at default.
  const cs_real_t r_inf_05 = cs_math_infinite_r * 0.5;

  cs_real_t *w5, *w7, *wbfa, *wbfb;
  cs_real_3_t *bc_vel;
  cs_real_t *bc_en, *bc_pr, *bc_tk;
  BFT_MALLOC(w5, n_cells_ext, cs_real_t);
  BFT_MALLOC(w7, n_b_faces, cs_real_t);
  BFT_MALLOC(wbfa, n_b_faces, cs_real_t);
  BFT_MALLOC(wbfb, n_b_faces, cs_real_t);
  BFT_MALLOC(bc_vel, n_b_faces, cs_real_3_t);
  BFT_MALLOC(bc_en, n_b_faces, cs_real_t);
  BFT_MALLOC(bc_pr, n_b_faces, cs_real_t);
  BFT_MALLOC(bc_tk, n_b_faces, cs_real_t);

  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_pr = CS_F_(p);
  cs_field_t *f_en = CS_F_(e_tot);
  cs_field_t *f_tk = CS_F_(t);

  auto *vel_icodcl = f_vel->bc_coeffs->icodcl;
  auto *vel_rcodcl1 = f_vel->bc_coeffs->rcodcl1;
  auto *vel_rcodcl3 = f_vel->bc_coeffs->rcodcl3;

  auto *pr_icodcl  = f_pr->bc_coeffs->icodcl;
  auto *pr_rcodcl1 = f_pr->bc_coeffs->rcodcl1;
  auto *pr_rcodcl2 = f_pr->bc_coeffs->rcodcl2;
  auto *pr_rcodcl3 = f_pr->bc_coeffs->rcodcl3;

  auto *tk_icodcl  = f_tk->bc_coeffs->icodcl;
  auto *tk_rcodcl1 = f_tk->bc_coeffs->rcodcl1;
  auto *tk_rcodcl3 = f_tk->bc_coeffs->rcodcl3;

  auto *en_icodcl  = f_en->bc_coeffs->icodcl;
  auto *en_rcodcl1 = f_en->bc_coeffs->rcodcl1;
  auto *en_rcodcl3 = f_en->bc_coeffs->rcodcl3;

  const cs_real_3_t *vel = (const cs_real_3_t *)(CS_F_(vel)->val);
  const cs_real_t *crom = (const cs_real_t *)(CS_F_(rho)->val);
  const cs_real_t *brom = (const cs_real_t *)(CS_F_(rho_b)->val);
  const cs_real_t *cpro_cv = nullptr;
  const cs_real_t *dt = CS_F_(dt)->val;

  if (fluid_props->icv >= 0)
    cpro_cv = static_cast<const cs_real_t *>
                (cs_field_by_id(fluid_props->icv)->val);

  cs_real_t *bc_fracv = nullptr, *bc_fracm = nullptr, *bc_frace = nullptr;
  int *fracv_icodcl = nullptr;
  int *fracm_icodcl = nullptr;
  int *frace_icodcl = nullptr;
  cs_real_t *fracv_rcodcl1 = nullptr;
  cs_real_t *fracm_rcodcl1 = nullptr;
  cs_real_t *frace_rcodcl1 = nullptr;

  // Mixture fractions for the homogeneous two-phase flows
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 2) {
    BFT_MALLOC(bc_fracv, n_b_faces, cs_real_t);
    BFT_MALLOC(bc_fracm, n_b_faces, cs_real_t);
    BFT_MALLOC(bc_frace, n_b_faces, cs_real_t);
    fracv_icodcl = CS_F_(volume_f)->bc_coeffs->icodcl;
    fracm_icodcl = CS_F_(mass_f)->bc_coeffs->icodcl;
    frace_icodcl = CS_F_(energy_f)->bc_coeffs->icodcl;
    fracv_rcodcl1 = CS_F_(volume_f)->bc_coeffs->rcodcl1;
    fracm_rcodcl1 = CS_F_(mass_f)->bc_coeffs->rcodcl1;
    frace_rcodcl1 = CS_F_(energy_f)->bc_coeffs->rcodcl1;
  }

  cs_array_copy(n_b_faces, f_pr->bc_coeffs->b, wbfb);

  /* Computation of epsilon_sup = e - CvT
     Needed if walls with imposed temperature are set. */

  int icalep = 0;
  for (auto face_id = 0; face_id < n_b_faces; face_id++) {
    if (tk_icodcl[face_id] == 5)
      icalep = 1;
  }

  if (icalep > 0) { // Local only, no need for parallel sync of icalep.
    cs_cf_thermo_eps_sup(crom, w5, m->n_cells);
    cs_cf_thermo_eps_sup(brom, w7, n_b_faces);
  }

  int *ifbet = cs_cf_boundary_conditions_get_ifbet();

  /* Main loop on boundary faces for series BC computation
     ----------------------------------------------------- */

  for (auto face_id = 0; face_id < n_b_faces; face_id++) {

    auto cell_id = b_face_cells[face_id];

    switch(bc_type[face_id]) {

    /* Wall faces
       --------- */

    case CS_SMOOTHWALL:
      [[fallthrough]];
    case CS_ROUGHWALL:
      {
        /* Pressure:
           if the gravity is prevailing: hydrostatic pressure
           (warning: the density is here explicit and the term is
           an approximation). */

        if (icfgrp == 1) {
          pr_icodcl[face_id] = 15;
          cs_real_t hint = dt[cell_id] / b_dist[face_id];
          pr_rcodcl3[face_id]
            = -hint * cs_math_3_distance_dot_product(b_face_cog[face_id],
                                                     cell_cen[cell_id],
                                                     grav)
              * crom[cell_id];
        }
        else {
          /* generally proportional to the bulk value
             (Pboundary = COEFB*Pi)
             The part deriving from pinf in stiffened gas is explicit for now */

          cs_cf_thermo_wall_bc(wbfa, wbfb, face_id);

          if (wbfb[face_id] < r_inf_05 && wbfb[face_id] > 0.) {
            pr_icodcl[face_id] = 12;
            pr_rcodcl1[face_id] = wbfa[face_id];
            pr_rcodcl2[face_id] = wbfb[face_id];
          }
          else {
            // If rarefaction is too strong: homogeneous Dirichlet
            pr_icodcl[face_id] = 13;
            pr_rcodcl1[face_id] = 0.;
          }
        }

        /* Velocity and turbulence are treated in a standard manner
           in cs_boundary_condition_set_coeffs.

           For thermal B.C., a pre-processing has be done here since the solved
           variable is the total energy
           (internal energy + epsilon_sup + kinetic energy).
           Especially, when a temperature is imposed on a wall, clptur treatment
           has to be prepared. Except for the solved energy all the variables rho
           and s will take arbitrarily a zero flux B.C. (their B.C.s are only
           used for gradient reconstructions and imposing something other than
           zero flux could bring out spurious values near the boundary layer). */

        // Adiabatic by default
        if (tk_icodcl[face_id] == 0 && en_icodcl[face_id] == 0) {
          tk_icodcl[face_id] = 3;
          tk_rcodcl3[face_id] = 0.;
        }

        /* Imposed temperature */
        if (tk_icodcl[face_id] == 5) {

          /* The value of the energy that leads to the right flux is imposed.
             However it should be noted that it is the B.C. for the diffusion
             flux. For the gradient reconstruction, something else will be
             needed. For example, a zero flux or an other B.C. respecting a
             profile: it may be possible to treat the total energy as the
             temperature, keeping in mind that the total energy contains
             the cinetic energy, which could make the choice of the profile more
             difficult. */

          en_icodcl[face_id] = 5;
          if (cpro_cv == nullptr)
            en_rcodcl1[face_id] = cv0 * tk_rcodcl1[face_id];
          else
            en_rcodcl1[face_id] = cpro_cv[cell_id] * tk_rcodcl1[face_id];
          en_rcodcl1[face_id] +=   0.5*cs_math_3_square_norm(vel[cell_id])
                                 + w5[cell_id];
          // w5 contains epsilon_sup

          /* fluxes in grad(epsilon_sup and kinetic energy) have to be zero
             since they are already accounted for in the energy diffusion term
             ifbet[face_id] = 1;

             Dirichlet condition on the temperature for gradient reconstruction
             used only in post-processing (typically Nusselt computation). */

          tk_icodcl[face_id] = 1;

        }

        /* Imposed flux */
        else if (tk_icodcl[face_id] == 3) {

          // zero flux on energy
          en_icodcl[face_id] = 3;
          en_rcodcl3[face_id] = tk_rcodcl3[face_id];

          // Fluxes in grad(epsilon_sup and cinetic energy) have to be zero
          // since they are already accounted for in the energy diffusion term
          ifbet[face_id] = 1;

          // zero flux for the possible temperature reconstruction
          tk_icodcl[face_id] = 3;
          tk_rcodcl3[face_id] = 0.;

        }
      }
      break;

    /* Imposed Inlet/outlet (for example: supersonic inlet)
       ---------------------------------------------------- */

    case CS_ESICF:
      {
        // We have
        //   - velocity,
        //   - 2 variables among P, rho, T, E (but not the couple (T,E)),
        //   - turbulence variables
        //   - scalars

        // We look for the variable to be initialized
        // (if a zero value has been given, it is not adapted, so it will
        // be considered as not initialized and the computation will stop
        // displaying an error message. The boundary density may
        // be pre-initialized to the cell density also, so is tested last.

        cs_real_t drom = abs(crom[cell_id] - brom[face_id]);

        int level = 0;
        int iccfth = 10000;
        if (pr_rcodcl1[face_id] < r_inf_05) {
          iccfth = 2*iccfth;
          level += 1;
        }
        if (tk_rcodcl1[face_id] < r_inf_05) {
          iccfth = 5*iccfth;
          level += 1;
        }
        if (en_rcodcl1[face_id] < r_inf_05) {
          iccfth = 7*iccfth;
          level += 1;
        }

        if (brom[face_id] > 0.  &&  (level < 2  ||  drom > cs_math_epzero)) {
          iccfth = 3*iccfth;
          level += 1;
        }

        if (level != 2) {
          bft_error
            (__FILE__, __LINE__, 0,
             _("Prescribed inlet conditions for compressible flow (CS_ESICF)\n"
               "   Two and only two independant variables among\n"
               "   P, rho, T and E must be given. Here %d variables have been\n"
               "   given; check the user-defined BC definitions."), level);
        }

        iccfth += 900;

        // Compute missing thermo variables among P, rho, T, E.
        bc_en[face_id] = en_rcodcl1[face_id];
        bc_pr[face_id] = pr_rcodcl1[face_id];
        bc_tk[face_id] = tk_rcodcl1[face_id];
        for (cs_lnum_t i = 0; i < 3; i++)
          bc_vel[face_id][i] = vel_rcodcl1[n_b_faces*i + face_id];

        cs_cf_thermo(iccfth, face_id, bc_en, bc_pr, bc_tk, bc_vel);
      }
      break;

    /* Outlet with imposed pressure
       ---------------------------- */

    case CS_SOPCF:
      {
        // If no value was given for P or if its value is negative,
        // abort the computation (a negative value could be possible,
        // but in most cases it would be an error).
        if (pr_rcodcl1[face_id] < -r_inf_05) {
          bft_error
            (__FILE__, __LINE__, 0,
             _("The pressure was not set at an outlet with imposed pressure."
               " (CS_SOPCF);\n"
               " check the user-defined BC definitions."));
        }

        bc_en[face_id] = en_rcodcl1[face_id];
        bc_pr[face_id] = pr_rcodcl1[face_id];
        bc_tk[face_id] = tk_rcodcl1[face_id];
        for (cs_lnum_t i = 0; i < 3; i++)
          bc_vel[face_id][i] = vel_rcodcl1[n_b_faces*i + face_id];

        cs_cf_thermo_subsonic_outlet_bc(bc_en, bc_pr, bc_vel, face_id);
      }
      break;

    /* Inlet with Ptot, Htot imposed (reservoir boundary conditions)
       ------------------------------------------------------------- */

    case CS_EPHCF:
      {
        // If values for Ptot and Htot were not given, the computation stops.

        // en_rcodcl1 contains the boundary total enthalpy values
        // prescribed by the user.

        if (pr_rcodcl1[face_id] < -r_inf_05 || en_rcodcl1[face_id] < -r_inf_05) {
          bft_error
            (__FILE__, __LINE__, 0,
             _("The total pressure or total enthalpy were not provided\n"
               "at inlet with total pressure and total enthalpy imposed"
               " (CS_EPHCF);\n"
               " check the user-defined BC definitions."));
        }

        bc_en[face_id] = en_rcodcl1[face_id];
        bc_pr[face_id] = pr_rcodcl1[face_id];
        bc_tk[face_id] = tk_rcodcl1[face_id];
        for (cs_lnum_t i = 0; i < 3; i++)
          bc_vel[face_id][i] = vel_rcodcl1[n_b_faces*i + face_id];

        cs_cf_thermo_ph_inlet_bc(bc_en, bc_pr, bc_vel, face_id);
      }
      break;

    /* Inlet  with imposed rho*U and rho*U*H
       ------------------------------------- */

    case CS_EQHCF:
      {
        //! TODO to be implemented
        bft_error
          (__FILE__, __LINE__, 0,
           _("Inlet with mass and enthalpy flow rate (CS_EQHCF)\n"
             "prescribed, but feature not iplemented/available yet."));

        // Use a scenario in which we have a 2-contact and a 3-relaxation
        // entering the domain. We determine the conditions on the interface
        // based on thermodynalics and use Rusanov for smoothing.

        // Both rho and vel must be given.
        if (en_rcodcl1[face_id] < -r_inf_05) {
          bft_error
            (__FILE__, __LINE__, 0,
             _("The mass or enthalpy flow rate were not specified at an inlet\n"
               "with imposed mass and enthalpy flow rate (CS_IEQHCF)\n"
               " check the user-defined BC definitions."));
        }
      }
      break;

    default:
      break;

    } // End for switch/case on bc type

    /* Complete the treatment for inlets and outlets:
     *  - boundary convective fluxes computation
     *    (analytical or Rusanov) if needed
     *  - B.C. code (Dirichlet or Neumann)
     *  - Dirichlet values
     *----------------------------------------------*/

    if (   bc_type[face_id] == CS_ESICF || bc_type[face_id] == CS_SSPCF
        || bc_type[face_id] == CS_EPHCF || bc_type[face_id] == CS_SOPCF
        || bc_type[face_id] == CS_EQHCF) {

      /* Boundary convective fluxes computation (analytical or Rusanov) if needed
         (gamma has already have been computed if Rusanov fluxes are computed) */

      // Compute Rusanov fluxes only for the imposed inlet for stability reasons.

      if (bc_type[face_id] == CS_ESICF) {
        // Dirichlet for velocity and pressure are computed in order to
        // impose the Rusanov fluxes in mass, momentum and energy balance.
        cs_cf_boundary_rusanov(face_id, bc_en, bc_pr, bc_vel);
      }

      // For the other types of inlets/outlets (subsonic outlet, QH inlet,
      // PH inlet), analytical fluxes are computed

      else if (bc_type[face_id] != CS_SSPCF) {
        // the pressure part of the boundary analytical flux is not added here,
        // but set through the pressure gradient boundary conditions (Dirichlet).
        cs_cf_boundary_analytical_flux(face_id, bc_en, bc_pr, bc_vel);
      }

      /* Copy of boundary values into the Dirichlet values array */

      if (bc_type[face_id] != CS_SSPCF) {
        en_rcodcl1[face_id] = bc_en[face_id];
        pr_rcodcl1[face_id] = bc_pr[face_id];
        tk_rcodcl1[face_id] = bc_tk[face_id];
        for (cs_lnum_t i = 0; i < 3; i++)
          vel_rcodcl1[n_b_faces*i + face_id]  = bc_vel[face_id][i];
        if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 2) {
          // FIXME fill bc_frac...
          assert(0);
          fracv_rcodcl1[face_id] = bc_fracv[face_id];
          fracm_rcodcl1[face_id] = bc_fracm[face_id];
          frace_rcodcl1[face_id] = bc_frace[face_id];
        }
      }
      else { // supersonic outlet
        en_rcodcl3[face_id] = 0.;
        pr_rcodcl3[face_id] = 0.;
        tk_rcodcl3[face_id] = 0.;
        for (cs_lnum_t i = 0; i < 3; i++)
          vel_rcodcl3[n_b_faces*i + face_id] = 0.;
        if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 2) {
          fracv_rcodcl1[face_id] = 0.;
          fracm_rcodcl1[face_id] = 0.;
          frace_rcodcl1[face_id] = 0.;
        }
      }

      /* Boundary conditions codes (Dirichlet or Neumann) */

      // P               : Neumann but pressure part of momentum flux is imposed
      //                   as a Dirichlet BC for the pressure gradient (code 13)
      // rho, U, E, T    : Dirichlet
      // k, R, eps, scal : Dirichlet/Neumann depending on the flux mass value

      if (bc_type[face_id] != CS_SSPCF) {
        // Pressure : - Dirichlet for the gradient computation, allowing to have
        //              the pressure part of the convective flux at the boundary
        //            - Homogeneous Neumann for the diffusion
        pr_icodcl[face_id]   = 13;
        // velocity
        vel_icodcl[face_id]  = 1;
        // total energy
        en_icodcl[face_id]   = 1;
        // temperature
        tk_icodcl[face_id]   = 1;
        // mixture fractions
        if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 2) {
          fracv_icodcl[face_id] = 1;
          fracm_icodcl[face_id] = 1;
          frace_icodcl[face_id] = 1;
        }
      }
      else { // supersonic outlet
        pr_icodcl[face_id]   = 3;
        vel_icodcl[face_id]  = 3;
        en_icodcl[face_id]   = 3;
        tk_icodcl[face_id]   = 3;
        // mixture fractions
        if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 2) {
          fracv_icodcl[face_id] = 3;
          fracm_icodcl[face_id] = 3;
          frace_icodcl[face_id] = 3;
        }
      }

      /* Turbulence and passive scalars:
         Dirichlet / Neumann depending on the mass flux,
         handled later in common part of the code
         (not specific to compressible flow). */

    } /* End case for open BC's */

  } // Loop on boundary faces.

  /* Free work arrays */

  BFT_FREE(w5);
  BFT_FREE(w7);
  BFT_FREE(wbfa);
  BFT_FREE(wbfb);
  BFT_FREE(bc_vel);
  BFT_FREE(bc_en);
  BFT_FREE(bc_pr);
  BFT_FREE(bc_tk);

  BFT_FREE(bc_fracv);
  BFT_FREE(bc_fracm);
  BFT_FREE(bc_frace);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate boundary flux indicator arrays.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_boundary_conditions_init(void)
{
  cs_base_at_finalize(_cf_boundary_conditions_finalize);

  BFT_REALLOC(_ifbet, cs_glob_mesh->n_b_faces, int);
  BFT_REALLOC(_icvfli, cs_glob_mesh->n_b_faces, int);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to boundary face indicator array of convection flux
 *        - 0 upwind scheme
 *        - 1 imposed flux
 */
/*----------------------------------------------------------------------------*/

int *
cs_cf_boundary_conditions_get_icvfli(void)
{
  return _icvfli;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to imposed thermal flux indicator at the boundary
 *        (some boundary contributions of the total energy eq. have to be
 *         cancelled)
 */
/*----------------------------------------------------------------------------*/

int *
cs_cf_boundary_conditions_get_ifbet(void)
{
  return _ifbet;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare (reset) condition coefficients specific to compressible flows.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_boundary_conditions_reset(void)
{
  const cs_lnum_t n_b_faces
    = cs_mesh_location_get_n_elts(CS_MESH_LOCATION_BOUNDARY_FACES)[0];

  /* Reset markers for:
     - use of Rusanov at boundary
     - forced convective flux at boundary
  */

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    _icvfli[i] = 0;
    _ifbet[i] = 0;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
