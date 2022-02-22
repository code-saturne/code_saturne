/*============================================================================
 * This function is called each time step to define physical properties
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
#include <string.h>

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
 * \file cs_user_physical_properties-scalar-drift.c
 *
 * \brief User definition of physical properties for scalars with a drift.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties(cs_domain_t *domain)
{
  /* Common variables and initializations
   * ------------------------------------ */

  /*! [init] */
  const cs_lnum_t n_cells = domain->mesh->n_cells;

  cs_real_t *cpro_viscl = CS_F_(mu)->val;

  /* Key id for drift scalar */
  const int keydri = cs_field_key_id("drift_scalar_model");

  /* Key id for diffusivity id */
  const int kivisl = cs_field_key_id("diffusivity_id");

  /* Number of fields */
  const int nfld = cs_field_n_fields();
  /*! [init] */

  /* The following examples should be adapted by the user
   * ==================================================== */

  /*
   *  Example: If thermophoresis is required, one MUST set the diffusivity
   *  -------
   *  (Brownian motion) to be variable in space and set the proper relation
   *  between the molecular diffusivity and T:
   *  ex: Kb x T x cuning /(3*pi*diamp(iscal)*cpro_viscl(iel))
   *
   *  Excluding:
   *  - temperature, enthalpy (handled above)
   *  - fluctuation variances (property equal to that of the associated scalar)
   *
   *  Below, we define the same diffusivity law for all scalars (except the
   *  ones excluded above).
   *  Values of this property must be defined at cell centers
   */

  /*! [example_1] */

  /* Loop over fields which are scalar with a drift */

  for (int iflid = 0; iflid < nfld; iflid++) {

    cs_field_t  *f = cs_field_by_id(iflid);

    /*! We only handle here scalars with a drift */

    if (! (f->type & CS_FIELD_VARIABLE))
      continue;

    int drift_flag = cs_field_get_key_int(f, keydri);

    if (drift_flag & CS_DRIFT_SCALAR_ADD_DRIFT_FLUX) {

      /* Position of variables, coefficients
       * ----------------------------------- */

      cs_real_t *cpro_taup = NULL, *cpro_taufpt = NULL, *cpro_viscls = NULL;

      /* Scalar's diffusivity (Brownian motion) */

      int ifcvsl = cs_field_get_key_int(f, kivisl);
      if (ifcvsl > -1)
        cpro_viscls = cs_field_by_id(ifcvsl)->val;

      /* Coefficients of drift scalar CHOSEN BY THE USER
         Values given here are fictitious */

      const cs_real_t diamp = 1.e-4;   /* particle diameter */
      const cs_real_t cuning = 1.;     /* Cuningham correction factor */
      const cs_real_t rhop = 1.e4;     /* particle density */

      /* Get corresponding relaxation time (cpro_taup) */

      char df_name[128]; df_name[127] = '\0';
      snprintf(df_name, 127, "drift_tau_%s", f->name);

      cpro_taup = cs_field_by_name(df_name)->val;

      /* Corresponding interaction time particle--eddies */

      if (drift_flag & CS_DRIFT_SCALAR_TURBOPHORESIS) {
        snprintf(df_name, 127, "drift_turb_tau_%s", f->name);
        cpro_taufpt = cs_field_by_name(df_name)->val;
      }

      /* Computation of the relaxation time of the particles
       * --------------------------------------------------- */

      const cs_real_t diamp2 = cs_math_pow2(diamp);

      if (diamp <= 1.e-6) {
        /* Cuningham's correction for submicronic particules */
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cpro_taup[c_id] =   cuning*diamp2*rhop / (18.*cpro_viscl[c_id]);
        }
      }
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cpro_taup[c_id] = diamp2*rhop / (18.*cpro_viscl[c_id]);
        }
      }

      /* Compute the interaction time particle--eddies (tau_fpt)
       * ------------------------------------------------------- */

      if (drift_flag & CS_DRIFT_SCALAR_TURBOPHORESIS) {

        const cs_turb_model_t *t_mdl = cs_glob_turb_model;
        const int ksigmas = cs_field_key_id("turbulent_schmidt");

        /* k-epsilon or v2-f models */
        if (t_mdl->itytur == 2 || t_mdl->itytur == 5) {
          const cs_real_t *cvar_k = CS_F_(k)->val;
          const cs_real_t *cvar_eps = CS_F_(eps)->val;
          const cs_real_t turb_schmidt = cs_field_get_key_double(f, ksigmas);
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
            cs_real_t xk = cvar_k[c_id];
            cs_real_t xeps = cvar_eps[c_id];
            cpro_taufpt[c_id] = (3./2.)*(cs_turb_cmu/turb_schmidt)*xk/xeps;
          }
        }

        /* Rij-epsilon models */
        else if (t_mdl->itytur == 3) {
          const cs_real_6_t *cvar_rij = (const cs_real_6_t *)(CS_F_(rij)->val);
          const cs_real_t *cvar_eps = CS_F_(eps)->val;
          cs_real_t beta1 = 0.5 + 3.0/(4.*cs_turb_xkappa);
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
            cs_real_t xk = 0.5 * (  cvar_rij[c_id][0]
                                  + cvar_rij[c_id][1]
                                  + cvar_rij[c_id][2]);
            cs_real_t xeps = cvar_eps[c_id];
            cpro_taufpt[c_id] = xk/xeps/beta1;
          }
        }

        /* k-omega models */
        if (t_mdl->itytur == 6) {
          const cs_real_t *cvar_omg = CS_F_(omg)->val;
          const cs_real_t turb_schmidt = cs_field_get_key_double(f, ksigmas);
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
            cs_real_t xomg = cvar_omg[c_id];
            cpro_taufpt[c_id] = (3./2.)*(1./turb_schmidt)/xomg;
          }
        }

      }

      /* Brownian diffusion at cell centers
       * ---------------------------------- */

      /* Stop if the diffusivity is not variable */
      if (ifcvsl < 0)
        bft_error(__FILE__, __LINE__, 0,
                  "The diffusivity is uniform while a variable diffusivity\n"
                  "is computed.");

      /* Temperature and density */
      if (CS_F_(t) == NULL)
        bft_error(__FILE__, __LINE__, 0,
                  "The temperature field on which physical properties depend\n"
                  "does not seem to be present.");

      const cs_real_t *cvar_t = CS_F_(t)->val;
      const cs_real_t *cpro_rom = CS_F_(rho)->val;

      /* Homogeneous to a dynamic viscosity */
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t xvart = cvar_t[c_id];
        cs_real_t rho = cpro_rom[c_id];
        cs_real_t viscl = cpro_viscl[c_id];
        cpro_viscls[c_id] =   rho*cs_physical_constants_kb*xvart
                            * cuning/(3.*cs_math_pi*diamp*viscl);
      }

    } /* End of test on drift */

  } /* End of loop on fields */

  /*! [example_1] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
