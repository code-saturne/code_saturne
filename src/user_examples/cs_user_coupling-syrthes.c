/*============================================================================
 * Code couplings definition with SYRTHES and code_saturne.
 *
 * 1) Define conjuguate heat transfer couplings with the SYRTHES code
 * 2) Define couplings with other instances of code_saturne
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
 * \file cs_user_coupling-syrthes.c
 *
 * \brief Code couplings definition with SYRTHES and code_saturne.
 *
 * See \ref user_coupling for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with SYRTHES code.
 *
 * This is done by calling the \ref cs_syr_coupling_define function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_syrthes_coupling(void)
{
  /*! [coupling_syrthes_1] */
  {
    int  verbosity = 1, plot = 1;
    float tolerance = 0.1;
    bool allow_nonmatching = false;

    /*-------------------------------------------------------------------------
     * Example 1:
     *
     * Boundary faces of group '3' coupled with instance named 'SYRTHES_01'.
     *-------------------------------------------------------------------------*/

    cs_syr_coupling_define("SYRTHES_01",
                           "3",               /* boundary criteria */
                           NULL,              /* volume_criteria */
                           ' ',               /* projection_axis */
                           allow_nonmatching,
                           tolerance,
                           verbosity,
                           plot);

  }
  /*! [coupling_syrthes_1] */

  /*! [coupling_syrthes_2] */
  {
    int  verbosity = 1, plot = 1;
    float tolerance = 0.1;
    bool allow_nonmatching = false;

    /*-------------------------------------------------------------------------
     * Example 2:
     *
     * Boundary faces of group 'Wall' coupled with 2D SYRTHES instance
     * named 'SYRTHES_02'.
     *-------------------------------------------------------------------------*/

    cs_syr_coupling_define("SYRTHES_02",
                           "Wall",            /* boundary criteria */
                           NULL,              /* volume_criteria */
                           'z',               /* projection_axis */
                           allow_nonmatching,
                           tolerance,
                           verbosity,
                           plot);

  }
  /*! [coupling_syrthes_2] */

  /*! [coupling_syrthes_3] */
  {
    int  verbosity = 1, plot = 1;
    float tolerance = 0.1;
    bool allow_nonmatching = false;

    /*-------------------------------------------------------------------------
     * Example 3:
     *
     * Cells in box with corners (0, 0, 0) and (1, 1, 1) coupled with
     * SYRTHES instance named 'Solid' (volume coupling).
     *-------------------------------------------------------------------------*/

    cs_syr_coupling_define("Solid",
                           NULL,                          /* boundary */
                           "box[0., 0., 0., 1., 1., 1.]", /* volume */
                           ' ',                           /* projection */
                           allow_nonmatching,
                           tolerance,
                           verbosity,
                           plot);
  }
  /*! [coupling_syrthes_3] */

  /* By default, conservativity forcing flag is switched off (value 0)
     If one wants to switch on the conservativity forcing flag:

     cs_syr_coupling_set_conservativity(1);
  */

  /* Only for a volume coupling:
      By default, implicit treatment is done. You can switch to
      an explicit treatment by using the following function:

     cs_syr_coupling_set_explicit_treatment();
  */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a volume exchange coefficient for SYRTHES couplings.
 *
 * \param[in]   coupling_id   Syrthes coupling id
 * \param[in]   syrthes_name  name of associated Syrthes instance
 * \param[in]   n_elts        number of associated cells
 * \param[in]   elt_ids       associated cell ids
 * \param[out]  h_vol         associated exchange coefficient (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_syrthes_coupling_volume_h(int               coupling_id,
                                  const char       *syrthes_name,
                                  cs_lnum_t         n_elts,
                                  const cs_lnum_t   elt_ids[],
                                  cs_real_t         h_vol[])
{
  CS_NO_WARN_IF_UNUSED(coupling_id);
  CS_NO_WARN_IF_UNUSED(syrthes_name);

  /* Example 1 of the computation of a constant volumic exchange coefficient
     ----------------------------------------------------------------------- */

  /*! [example_1] */
  cs_real_t hvol_cst = 1.0e6;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    h_vol[i] = hvol_cst;
  }
  /*! [example_1] */

  /* Example 2 of the computation of a volumic exchange coefficient
   * --------------------------------------------------------------
   *
   *  hvol(iel) =  hsurf(iel) * exchange_surface_by_unit_vol
   *
   *  with: hsurf = Nusselt * lambda / L
   *
   *  lambda is the thermal conductivity coefficient
   *   L is a characteristic length
   *
   *  Nusselt is computed by means of the Colburn correlation
   *
   *  Nu = 0.023 * Re^(0.8) * Pr^(1/3)
   *
   *  Re is the Reynolds number and Pr is the Prandtl number
   */

  /*! [example_2_init] */
  const cs_real_3_t  *cvar_vel = (const cs_real_3_t *)CS_F_(vel)->val;
  const cs_real_t  *cpro_rom = (const cs_real_t *)CS_F_(rho)->val;
  const cs_real_t  *cpro_mu = (const cs_real_t *)CS_F_(mu)->val;

  const cs_real_t  cp0 = cs_glob_fluid_properties->cp0;
  cs_real_t  visls_0 = -1;

  const cs_real_t  *cpro_cp = NULL, *cpro_viscls = NULL;
  cs_lnum_t  cp_step = 0, viscls_step = 0;

  if (CS_F_(cp) != NULL) {
    cpro_cp = (const cs_real_t *)CS_F_(cp)->val;
    cp_step = 1;
  }
  else {
    cpro_cp = &cp0;
  }

  /* Get thermal field and associated diffusivity
     (temperature only handled here) */

  const cs_field_t *fth = cs_thermal_model_field();

  const int viscl_id = cs_field_get_key_int(fth,
                                            cs_field_key_id("diffusivity_id"));

  if (viscl_id > -1) {
    cpro_viscls = (const cs_real_t *)cs_field_by_id(viscl_id)->val;
    viscls_step = 1;
  }
  else {
    visls_0 = cs_field_get_key_double(fth, cs_field_key_id("diffusivity_ref"));
    cpro_viscls = &visls_0;
  }

  const int is_temperature
    = cs_field_get_key_int(fth,
                           cs_field_key_id("is_temperature"));
  /*! [example_2_init] */

  /*! [example_2] */
  cs_real_t sexcvo = 36.18;  /* Surface area where exchanges take
                                place by unit of volume */
  cs_real_t l0 = 0.03;       /* Characteristic length */

  for (cs_lnum_t i = 0; i < n_elts; i++) {

    cs_lnum_t c_id = elt_ids[i];

    cs_real_t rho = cpro_rom[c_id];
    cs_real_t mu = cpro_mu[c_id];
    cs_real_t cp = cpro_cp[c_id*cp_step];

    cs_real_t lambda, lambda_over_cp;

    if (is_temperature) {
      lambda = cpro_viscls[c_id*viscls_step];
      lambda_over_cp = lambda / cp;
    }
    else {
      lambda_over_cp = cpro_viscls[c_id*viscls_step];
      lambda =  lambda_over_cp * cp;
    }

    /* Compute a local molecular Prandtl **(1/3) */

    cs_real_t  pr = mu / lambda_over_cp;

    /* Compute a local Reynolds number */

    cs_real_t uloc = cs_math_3_norm(cvar_vel[c_id]);
    cs_real_t re = fmax(uloc*rho*l0/mu, 1.);  /* To avoid division by zero */

    /* Compute Nusselt number using Colburn correlation */

    cs_real_t nu = 0.023 * pow(re, 0.8) * pow(pr, 1./3.);
    cs_real_t h_corr = nu * lambda / l0;

    /* Compute hvol */
    h_vol[i] = h_corr * sexcvo;
  }
  /*! [example_2] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
