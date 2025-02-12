/*============================================================================
 * Libby-Williams gas combustion model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_array_reduce.h"
#include "base/cs_base.h"
#include "base/cs_dispatch.h"
#include "cdo/cs_equation_param.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_physical_constants.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_parall.h"
#include "base/cs_physical_properties.h"
#include "base/cs_restart.h"
#include "base/cs_restart_default.h"
#include "mesh/cs_mesh.h"
#include "rayt/cs_rad_transfer.h"
#include "turb/cs_turbulence_model.h"

#include "pprt/cs_combustion_model.h"
#include "pprt/cs_physical_model.h"

#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_combustion_boundary_conditions.h"
#include "cogz/cs_combustion_ht_convert.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_combustion_ebu.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_lw.cpp
        Libby-Williams gas combustion model.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize specific fields for Libby-Williams gas combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_lw_fields_init(void)
{
  // Only when not a restart
  if (cs_restart_present())
    return;

  // Local variables

  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  cs_real_t *cvar_ygfm = cm->ygfm->val;

  cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  for (int igg = 0; igg < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; igg++)
    coefg[igg] = 0;

  // Initializations with air at tinitk
  // ----------------------------------

  // Mass fraction of fresh gas
  cs_array_real_set_scalar(n_cells_ext, 1.0, cvar_ygfm);

  // Mixture enthalpy
  if (cm->type % 2 == 1) {
    // mixture temperature: air at tinitk
    cs_real_t tinitk = cs_glob_fluid_properties->t0;

    // Air enthalpy at tinitk
    coefg[0] = 0.;
    coefg[1] = 1.;
    coefg[2] = 0.;
    cs_real_t hair = cs_gas_combustion_t_to_h(coefg, tinitk);

    // Mixture enthalpy
    cs_real_t *cvar_scalt = CS_F_(h)->val;
    cs_array_real_set_scalar(n_cells_ext, hair, cvar_scalt);
  }

  // Mass fraction
  cs_array_real_set_scalar(n_cells_ext, cm->lw.fmax, cm->yfm->val);

  // Mixture fraction
  cs_array_real_set_scalar(n_cells_ext, cm->lw.fmax, cm->fm->val);

  // No need to set yfp2m, fp2m, or coyfp to 0,
  // as this is the default for all fields.
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute physical properties for premixed flame Libby-Williams
 *        combustion model.
 *
 * Define the source terms for a given scalar over one time step.
 *
 * The equations read: \f$ rovsdt \delta a = smbrs \f$
 *
 * \f$ rovsdt \f$ et \f$ smbrs \f$ may already contain source term
 * so must not be overwritten, but incremented.
 *
 * For stability, only positive terms should be add in \f$ rovsdt \f$.
 * There is no constraint for \f$ smbrs \f$.
 * For a source term written \f$ S_{exp} + S_{imp} a \f$, source terms are:
 *           \f$ smbrs  = smbrs  + S_{exp} + S_{imp} a \f$
 *           \f$ rovsdt = rovsdt + \max(-S_{imp},0) \f$
 *
 * Here we set \f$ rovsdt \f$ and \f$ smbrs \f$ containing \f$ \rho \Omega \f$
 *   - \f$ smbrs \f$ in \f$ kg_a.s^{-1} \f$ (ex: for velocity:
 *     \f$ kg.m.s^{-2} \f$, for temperature: \f$ kg.C.s^{-1} \f$,
 *     for enthalpy: \f$ J.s^{-1} \f$)
 *   - \f$ rovsdt \f$ in \f$ kg.s^{-1} \f$
 *
 * \param[in]      f_sc          pointer to scalar field
 * \param[in,out]  smbrs         explicit right hand side
 * \param[in,out]  rovsdt        implicit terms
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_lw_source_terms(cs_field_t  *f_sc,
                              cs_real_t    smbrs[],
                              cs_real_t    rovsdt[])
{
  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  if (f_sc != cm->ygfm)
    return;

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f_sc);
  if (eqp->verbosity >= 1)
    bft_printf(_("Source terms for variable %s\n\n"), f_sc->name);

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  // Get variables and coefficients

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t *volume = cs_glob_mesh_quantities->cell_vol;

  const cs_real_t *cvara_k = nullptr;
  const cs_real_t *cvara_ep = nullptr;
  const cs_real_t *cvara_omg = nullptr;
  const cs_real_6_t *cvara_rij = nullptr;

  const cs_real_t *cvara_scal = f_sc->val_pre;
  const cs_real_t *cvara_yfm = cm->yfm->val_pre;
  const cs_real_t *cvara_fm = cm->fm->val_pre;

  cs_real_t *cpro_fmel[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_fmal[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_tscl[CS_COMBUSTION_GAS_MAX_DIRAC];
  cs_real_t *cpro_rhol[CS_COMBUSTION_GAS_MAX_DIRAC];

  const int n_dirac = cm->lw.n_dirac;
  for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
    cpro_fmel[dirac_id] = cm->lw.fmel[dirac_id]->val;
    cpro_fmal[dirac_id] = cm->lw.fmal[dirac_id]->val;
    cpro_rhol[dirac_id] = cm->lw.rhol[dirac_id]->val;
    cpro_tscl[dirac_id] = cm->lw.tscl[dirac_id]->val;
  }

  const cs_real_t epsi = 1.0e-10;

  cs_host_context ctx;

  /* Source term for mean fuel mass fraction
     --------------------------------------- */

  if (f_sc == cm->yfm) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      cs_real_t sum = 0.;
      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
        sum +=   cpro_rhol[dirac_id][c_id]
               * cpro_tscl[dirac_id][c_id] * volume[c_id];
      }

      // implicit term
      if (cvara_scal[c_id] > epsi)
        rovsdt[c_id] += fmax(-sum/cvara_scal[c_id], 0.0);

      // explicit term
      smbrs[c_id] += sum;
    });

  }

  /* Source term for variance of mean fuel mass fraction
     --------------------------------------------------- */

  if (f_sc == cm->yfp2m) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      cs_real_t sum = 0.;
      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
        sum  +=  (cpro_tscl[dirac_id][c_id]*volume[c_id]
                *(cpro_fmal[dirac_id][c_id] - cvara_yfm[c_id])
                *cpro_rhol[dirac_id][c_id]);
      }
      smbrs[c_id] += sum;
    });

  }

  /* Covariance source term
     ---------------------- */

  if (f_sc == cm->coyfp) {

    // Allocate work arrays

    cs_real_3_t *gradf, *grady;
    CS_MALLOC(gradf, n_cells_ext, cs_real_3_t);
    CS_MALLOC(grady, n_cells_ext, cs_real_3_t);

    // Gradient of F
    cs_field_gradient_scalar(cm->fm, true, 1, gradf);

    // Gradient of Yfuel
    cs_field_gradient_scalar(cm->yfm, true, 1, grady);

    // Access or reconstruct k and Epsilon based on turbulence model.

    cs_real_t  *w1 = nullptr;

    if (   cs_glob_turb_model->itytur == 2
        || cs_glob_turb_model->itytur == 5) {
      cvara_k = CS_F_(k)->val_pre;
      cvara_ep = CS_F_(eps)->val_pre;
    }

    else if (cs_glob_turb_model->itytur == 3) {
      cvara_rij = (const cs_real_6_t *)(CS_F_(k)->val_pre);
      cvara_ep = CS_F_(eps)->val_pre;

      CS_MALLOC(w1, n_cells_ext, cs_real_t);
      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
        w1[c_id] = 0.5 * cs_math_6_trace(cvara_rij[c_id]);
      });
      cvara_k = w1;
    }

    else if (cs_glob_turb_model->model == CS_TURB_K_OMEGA) {
      cvara_k = CS_F_(k)->val_pre;
      cvara_omg = CS_F_(omg)->val_pre;

      const cs_real_t cmu = cs_turb_cmu;
      CS_MALLOC(w1, n_cells_ext, cs_real_t);
      ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
        w1[c_id] = cmu * cvara_k[c_id]* cvara_omg[c_id];
      });
      cvara_ep = w1;
    }

    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    const int krvarfl = cs_field_key_id("variance_dissipation");
    cs_real_t turb_schmidt = cs_field_get_key_double(f_sc, ksigmas);
    cs_real_t rvarfl = cs_field_get_key_double(f_sc, krvarfl);

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {

      // To confirm:
      // The dissipation term should be implicit.
      // In the dissipation term, a constant Cf is missing.
      //   Can it be conidered == 1 ?
      // Check the sign of the production term.

      // Implicit term

      cs_real_t w11 =   cvara_ep[c_id] / (cvara_k[c_id] * rvarfl)
                      * volume[c_id] * crom[c_id];
      rovsdt[c_id] = rovsdt[c_id] + fmax(w11, 0.);

      // Gradient term

      cs_real_t tsgrad =  (  2.0
                           * visct[c_id]/(turb_schmidt)
                           * cs_math_3_dot_product(gradf[c_id], grady[c_id]))
                         * volume[c_id];

      // Dissipation term

      cs_real_t tsdiss = -w11 * cvara_scal[c_id];

      // Chemical term

      cs_real_t tschim = 0.;
      for (int dirac_id = 0; dirac_id < n_dirac; dirac_id++) {
        tschim += ( cpro_tscl[dirac_id][c_id]
                   *(cpro_fmel[dirac_id][c_id]-cvara_fm[c_id])
                   *volume[c_id]) * cpro_rhol[dirac_id][c_id];
      }

      // Sum of sources

      smbrs[c_id] += tschim + tsgrad + tsdiss;
    });

    CS_FREE(w1);
    CS_FREE(grady);
    CS_FREE(gradf);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
