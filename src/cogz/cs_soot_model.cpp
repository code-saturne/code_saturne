/*============================================================================
 * Soot production models.
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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "base/cs_array.h"
#include "pprt/cs_combustion_model.h"
#include "cogz/cs_cogz_headers.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cogz/cs_soot_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_soot_model.cpp
        Soot related models.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

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

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Specific physical models: soot production model.
 *
 * Define the source terms for the soot mass fraction
 * and the precursor number for soot model.
 * The equations read: \f$ rovsdt \delta a = smbrs \f$
 *
 * \f$ rovsdt \f$ et \f$ smbrs \f$ may already contain source term
 * so must not be overwritten, but incremented.
 *
 * For stability only positive terms should be add in \f$ rovsdt \f$.
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
 * \param[in]     f_id           scalar field id
 * \param[out]    smbrs          explicit right hand side
 * \param[out]    rovsdt         implicit terms
 */
/*----------------------------------------------------------------------------*/

void
cs_soot_production(int        f_id,
                   cs_real_t  smbrs[],
                   cs_real_t  rovsdt[])
{
  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  const int isoot = cm->isoot;
  const int n_gas_e = cm->n_gas_el_comp;
  const int n_gas_g = cm->n_gas_species;
  const cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t epsi = 1.e-6;

  /* Get fields */
  cs_field_t *f = cs_field_by_id(f_id);

  cs_real_t *cvar_temp = cs_field_by_name_try("temperature")->val;
  cs_real_t *cvara_scal = f->val_pre;
  cs_real_t *cvar_ym1 = cs_field_by_name("ym_fuel")->val;
  cs_real_t *cvar_ym2 = cs_field_by_name("ym_oxyd")->val;
  cs_real_t *cvar_ym3 = cs_field_by_name("ym_prod")->val;
  cs_real_t *cvara_ys = CS_F_(fsm)->val_pre;
  cs_real_t *cvara_yp = CS_F_(npm)->val_pre;

  /* Equation parameters */
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  /* Mesh quantities */
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const int n_cells = mesh->n_cells;
  const cs_real_t *volume = mq->cell_vol;

  /* Logging  */
  if (eqp->verbosity >= 1)
    bft_printf("Soot production calculation.\n");

  /*--- Moss et al.:
    zeta_s (ifsm) soot mass fraction zeta_s = (rho_s/rho).f_v
    zeta_n (inpm) precursor density  zeta_n = n / (rho.No)
    ----------------------------------------------------------------*/

  if (isoot == 1) {

    /* To be changed for other combustible !FIXME
       Methane CH4 (Syed, Stewart and Moss Symposium 1990) */
    cs_real_t caa = 6.54e4; // m^3/kg^2.K^0.5.s
    cs_real_t cbb = 1.3e7;  // m^3.K^-1/2.s^-1
    cs_real_t ccc = 0.1;    // m^3.kg^-2/3.K^-1/2.s^-1
    cs_real_t taa = 46.1e3; // K
    cs_real_t tcc = 12.6e3; // K
    cs_real_t d1s3 = 1./3.;

    # pragma omp parallel for  if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t po2, ka, kb, kt, kz, zetas, zetan;
      cs_real_t chi, wox, aa, bb, cc, dd;

      cs_real_t rho = crom[c_id];
      cs_real_t temp = cvar_temp[c_id];

      cs_real_t cexp = 0.;
      cs_real_t cimp = 0.;
      cs_real_t nn0 = 6.0223e23;

      cs_real_t xm = 1./ (  cvar_ym1[c_id] / cm->wmolg[0]
                          + cvar_ym2[c_id] / cm->wmolg[1]
                          + cvar_ym3[c_id] / cm->wmolg[2]);

      /* Fuel molar fraction */
      cs_real_t xfu = cvar_ym1[c_id] * xm / cm->wmolg[0];

      /* Rate of particule nucleation */
      aa = caa * cs_math_pow2(rho) * sqrt(temp) * xfu * exp(-taa/temp);

      /* Coagulation */
      bb = cbb * sqrt(temp);

      /* Surface growth of soot */
      cc = ccc * rho * sqrt(temp) * xfu * exp(-tcc/temp);
      po2 = cvar_ym2[c_id] * xm / cm->wmolg[1] * 1./4.76;

      /* Oxidation */
      ka = 20. * exp(-15098./temp);
      kb = 4.46e-3 * exp(-7650./temp);
      kt = 1.51e5 * exp(-48817./temp);
      kz = 21.3 * exp(2063./temp);

      chi = kb * po2 / (kb * po2 + kt);

      wox =   1.2e2 * (  (ka * po2 * chi) / (1. + kz * po2)
                       + kb * po2 * (1. - chi));

      dd = pow(36. * acos(-1.) / cs_math_pow2(cm->rosoot), d1s3);

      zetas = cvara_ys[c_id];
      zetan = cvara_yp[c_id];

      if (f_id == CS_F_(fsm)->id) {
        /* Surface growth : quadratic */
        if (zetas >= epsi) {
          cimp =   volume[c_id] *(pow(nn0, d1s3) * rho * cc * pow(zetas, -d1s3)
                 * pow(zetan, d1s3) - rho * dd * pow(nn0, d1s3)
                 * pow(zetan, d1s3) * pow(zetas, -d1s3) * wox);
        }
        cexp = volume[c_id] * (144. * aa);
      }
      if (f_id == CS_F_(npm)->id) {
        cimp = volume[c_id] * (-cs_math_pow2(rho) * bb * zetan);
        cexp = volume[c_id] * aa;
      }
      smbrs[c_id] += cexp + cimp * cvara_scal[c_id];
      rovsdt[c_id] += cs::max(-cimp, 0.);
    }

  }
  else if (isoot == 2) {

    const cs_real_t lsp_c2h4 = 0.106, af_c2h4 = 4.e-5;
    cs_real_t fst = cm->fs[0];
    cs_real_t a0 = 120., as = 160.e3;
    cs_real_t gama = 2.25, tact = 2000., yft = 1.;

    /* lsp_fuel set by the user */
    const cs_real_t lsp_fuel = cm->lsp_fuel;
    cs_real_t *cvar_ys  = CS_F_(fsm)->val;
    const cs_real_t zso = 1.25 * fst;
    const cs_real_t zsf   = 2.50 * fst;
    const cs_real_t afuel = af_c2h4 * lsp_c2h4 / lsp_fuel;

    if (f_id == CS_F_(fsm)->id) {

      cs_real_t *wsf, *wso;
      CS_MALLOC(wsf, n_cells, cs_real_t);
      CS_MALLOC(wso, n_cells, cs_real_t);

      cs_array_real_fill_zero(n_cells, wso);
      cs_array_real_fill_zero(n_cells, wsf);

      const cs_real_t *cvar_fm = CS_F_(fm)->val;

      #pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        const cs_real_t fm    = cvar_fm[c_id];
        const cs_real_t rho   = crom[c_id];      // Mixture density (kg/m3)
        const cs_real_t temp  = cvar_temp[c_id]; // Temperature
        const cs_real_t zetas = cvar_ys[c_id];   // Soot Mass Fraction

        cs_real_t yg[n_gas_g], ye[n_gas_e], xe[n_gas_e];

        yg[0] = cvar_ym1[c_id];
        yg[1] = cvar_ym2[c_id];
        yg[2] = cvar_ym3[c_id];

        cs_combustion_gas_yg2xye(yg, ye, xe);

        const cs_real_t cxo2 = rho * ye[1] / cm->wmole[1];

        if ((fm >= zso) && (fm <= zsf)) {
          wsf[c_id] =   afuel * pow(rho, 2) * (yft *(fm - fst)
                      / (1. - fst)) * pow(temp, gama) * exp(-tact / temp);
          wso[c_id] = -rho * a0 * cxo2 * sqrt(temp) * exp(-19670./temp) * as;
        }
        if (fm < zso) {
          wso[c_id] = -rho * a0 * cxo2 * sqrt(temp) * exp(-19670./temp) * as;
        }
        smbrs[c_id]  += (wsf[c_id] + wso[c_id] *zetas) * volume[c_id];
        rovsdt[c_id] += cs::max(-wso[c_id], 0.) * volume[c_id];
      }

      /* Free memory */
      CS_FREE(wsf);
      CS_FREE(wso);

    }

    else if (f_id == CS_F_(npm)->id) {
      cs_array_real_fill_zero(n_cells, smbrs);
      cs_array_real_fill_zero(n_cells, rovsdt);
    }

  } /* isoot = 2 */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
