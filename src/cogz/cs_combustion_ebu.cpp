/*============================================================================
 * EBU (Eddy Break-Up) gas combustion model.
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
 * Standard C library headers
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
  \file cs_combustion_ebu.cpp
        Eddy-Break-Up gas combustion model.
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize specific fields for EBU gas combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_ebu_fields_init(void)
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

  // No need to set fm to 0, as this is the default for all fields.
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute physical properties for EBU combustion model.
 *
 * \param[in, out]   mbrom    filling indicator of romb
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_ebu_physical_prop(int  *mbrom)
{
  static int n_calls = 0;  // Count calls
  n_calls++;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_rad_transfer_model_t rt_model = cs_glob_rad_transfer_params->type;
  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  const int sub_type = cm->type % 100;

  /* Initializations
     --------------- */

  // Allocate temporary arrays

  cs_real_t  *yfuegf, *yoxygf, *yprogf, *yfuegb, *yoxygb, *yprogb;
  cs_real_t  *temp, *masmel;
  CS_MALLOC(yfuegf, n_cells_ext, cs_real_t);
  CS_MALLOC(yoxygf, n_cells_ext, cs_real_t);
  CS_MALLOC(yprogf, n_cells_ext, cs_real_t);
  CS_MALLOC(yfuegb, n_cells_ext, cs_real_t);
  CS_MALLOC(yoxygb, n_cells_ext, cs_real_t);
  CS_MALLOC(yprogb, n_cells_ext, cs_real_t);
  CS_MALLOC(temp, n_cells_ext, cs_real_t);
  CS_MALLOC(masmel, n_cells_ext, cs_real_t);

  // Get variables and coefficients

  cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *cpro_temp = CS_F_(t)->val;
  cs_real_t *cpro_ym1 = cm->ym[0]->val;
  cs_real_t *cpro_ym2 = cm->ym[1]->val;
  cs_real_t *cpro_ym3 = cm->ym[2]->val;

  cs_real_t *cpro_ckabs = nullptr, *cpro_t4m = nullptr, *cpro_t3m = nullptr;
  if (rt_model != CS_RAD_TRANSFER_NONE) {
    cpro_ckabs = cm->ckabs->val;
    cpro_t4m = cm->t4m->val;
    cpro_t3m = cm->t3m->val;
  }

  cs_real_t *cvar_ygfm = cm->ygfm->val;
  cs_real_t *cvar_fm = nullptr, *cvar_scalt = nullptr;
  if (sub_type == 2 || sub_type == 3)
    cvar_fm = cm->fm->val;
  if (cm->type %2 == 1)
    cvar_scalt = CS_F_(h)->val;

  bool update_rho = false;
  if (n_calls > 1 || cs_restart_get_field_read_status(CS_F_(rho)->id) == 1)
    update_rho = true;

  /*
   * Determine thermochemical quantities
   * ----------------------------------- */

  // Fresh gas

  // User-defined
  //   fmel        --> Mixture fraction
  //                   constant for options 0 et 1, variable otherwise
  //   tgf         --> Fresh gas temperature in identical K
  //                   for fresh premix and dilution

  // Deducted
  //   yfuegf      --> Fresh gas fuel mass fraction
  //   yoxygf      --> Fresh gas oxidant mass fraction
  //   yprogf      --> Fresh gas products mass fraction
  //   hgf         --> Identical fresh gas mass enthaly
  //                   for fresh premix and dilution
  //   masmgf      --> Fresh gas molar mass
  //   ckabgf      --> Absorption coefficient

  // Burned gas

  // Deducted
  //   tgb         --> Burned gases temperature in K
  //   yfuegb      --> Burned gases fuel mass fraction
  //   yoxygb      --> Burned gases oxidant mass fraction
  //   yprogb      --> Burned gases products mass fraction
  //   masmgb      --> Burned gases molar mass
  //   ckabgb      --> Absorption coefficient

  // Mixture quantities

  //   masmel      --> Mixture Molar mass
  //   cpro_temp   --> Mixture temperature
  //   crom        --> Mixture density
  //   cpro_.f,o,p --> Mass fractions in F, O, P
  //   cpro_ckabs  --> Absorption coefficient
  //   cpro_t4m    --> T^4 term
  //   cpro_t3m    --> T^3 term

  // Mass fractions of fresh and burned gases in F, O, P

  cs_host_context ctx;

  cs_real_t frmel = cm->frmel;
  const cs_real_t fs[1] = {cm->fs[0]};

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    cs_real_t fmel = (cvar_fm != nullptr) ? cvar_fm[c_id] : frmel;

    yfuegf[c_id] = fmel;
    yoxygf[c_id] = 1.0 - fmel;
    yprogf[c_id] = 0.0;

    yfuegb[c_id] = fmax(0., (fmel-fs[0]) / (1.0-fs[0]));
    yprogb[c_id] = (fmel-yfuegb[c_id]) / fs[0];
    yoxygb[c_id] = 1.0 - yfuegb[c_id] - yprogb[c_id];
  });

  const cs_real_t epsi = 1e-6;
  const cs_real_t tgf = cm->tgf;
  const cs_real_t srrom = cm->srrom;
  const cs_real_t c_r = cs_physical_constants_r;
  const cs_real_t p0 = cs_glob_fluid_properties->p0;

  const cs_real_t ckabsg[3] = {cm->ckabsg[0], cm->ckabsg[1], cm->ckabsg[2]};
  const cs_real_t wmolg[3] = {cm->wmolg[0], cm->wmolg[1], cm->wmolg[2]};

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    cs_real_t ckabgf, ckabgb;
    if (cpro_ckabs != nullptr) {
      ckabgf =    yfuegf[c_id]*ckabsg[0] + yoxygf[c_id]*ckabsg[1]
                + yprogf[c_id]*ckabsg[2];
      ckabgb  =   yfuegb[c_id]*ckabsg[0] + yoxygb[c_id]*ckabsg[1]
                + yprogb[c_id]*ckabsg[2];
    }
    else {
      ckabgf = 0.;
      ckabgb = 0.;
    }

    // Molar mass of fresh gases
    cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
    coefg[0] = yfuegf[c_id];
    coefg[1] = yoxygf[c_id];
    coefg[2] = yprogf[c_id];
    cs_real_t nbmol = 0.0;
    for (int igg = 0; igg < 3; igg++)
      nbmol += coefg[igg] / wmolg[igg];
    cs_real_t masmgf = 1. / nbmol;

    // Enthalpy of fresh gases
    cs_real_t hgf = cs_gas_combustion_t_to_h(coefg, tgf);

    // Molar masse of burned gasses

    coefg[0] = yfuegb[c_id];
    coefg[1] = yoxygb[c_id];
    coefg[2] = yprogb[c_id];
    nbmol = 0.0;
    for (int igg = 0; igg < 3; igg++)
      nbmol += coefg[igg] / wmolg[igg];
    cs_real_t masmgb = 1.0 / nbmol;

    cs_real_t ygfm = cvar_ygfm[c_id];
    cs_real_t ygbm = 1.0 - ygfm;

    // Molar mass of mixture
    masmel[c_id] = 1.0 / (ygfm/masmgf + ygbm/masmgb);

    // Temperature of burned gases
    cs_real_t hgb = hgf;
    if (sub_type%2 == 1 && ygbm > epsi)
      hgb = (cvar_scalt[c_id] - hgf*ygfm) / ygbm;
    cs_real_t tgb = cs_gas_combustion_h_to_t(coefg, hgb);

    // Mixture temperature
    //   Rq PPl: it would be better to weight by the CP (GF et GB)
    cpro_temp[c_id] = ygfm*tgf + ygbm*tgb;

    // Temperature / molar mass
    cs_real_t temsmm = ygfm*tgf/masmgf + ygbm*tgb/masmgb;

    // Mixture density

    if (update_rho) {
      crom[c_id] =   srrom*crom[c_id]
                   + (1.0-srrom) * (p0/(c_r*temsmm));
    }

    // Mass fractions of global species
    cpro_ym1[c_id] = yfuegf[c_id]*ygfm + yfuegb[c_id]*ygbm;
    cpro_ym2[c_id] = yoxygf[c_id]*ygfm + yoxygb[c_id]*ygbm;
    cpro_ym3[c_id] = yprogf[c_id]*ygfm + yprogb[c_id]*ygbm;

    // Radiation quantities
    if (cpro_ckabs != nullptr) {
      cpro_ckabs[c_id] = ygfm*ckabgf + ygbm*ckabgb;
      cpro_t4m[c_id]  = ygfm*cs_math_pow4(tgf) + ygbm*cs_math_pow4(tgb);
      cpro_t3m[c_id]  = ygfm*cs_math_pow3(tgf) + ygbm*cs_math_pow3(tgb);
    }
  });

  /*
   * Compute rho and mass fractions of global species at boundaries
   * -------------------------------------------------------------- */

  cs_combustion_boundary_conditions_density_ebu_lw();

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

  for (int igg = 0; igg < cm->n_gas_species; igg++) {
    cs_real_t *bsval = cm->bym[igg]->val;
    cs_real_t *cpro_ymgg = cm->ym[igg]->val;

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST (cs_lnum_t f_id) {
      cs_lnum_t c_id = b_face_cells[f_id];
      bsval[f_id] = cpro_ymgg[c_id];
    });
  }

  CS_FREE(yfuegf);
  CS_FREE(yoxygf);
  CS_FREE(yprogf);
  CS_FREE(yfuegb);
  CS_FREE(yoxygb);
  CS_FREE(yprogb);
  CS_FREE(temp);
  CS_FREE(masmel);

  *mbrom = 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute physical properties for premixed flame EBU combustion model.
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
cs_combustion_ebu_source_terms(cs_field_t  *f_sc,
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
  const cs_real_t *cvara_scal = f_sc->val_pre;
  const cs_real_t *volume = cs_glob_mesh_quantities->cell_vol;

  const cs_real_t *cvara_k = nullptr;
  const cs_real_t *cvara_ep = nullptr;
  const cs_real_t *cvara_omg = nullptr;
  const cs_real_6_t *cvara_rij = nullptr;

  cs_real_t  *w1 = nullptr;

  cs_host_context ctx;

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

  const cs_real_t epzero = cs_math_epzero;
  const cs_real_t cebu = cm->cebu;

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    if (cvara_k[c_id] > epzero && cvara_ep[c_id] > epzero) {
      cs_real_t w3 =   cebu*cvara_ep[c_id] / cvara_k[c_id]
                     * crom[c_id] * volume[c_id] * (1.0 - cvara_scal[c_id]);
      smbrs[c_id] = smbrs[c_id] - cvara_scal[c_id]*w3;
      rovsdt[c_id] = rovsdt[c_id] + fmax(w3, 0.0);
    }
  });

  CS_FREE(w1);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
