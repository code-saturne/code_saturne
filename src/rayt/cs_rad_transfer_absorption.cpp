/*============================================================================
 * Radiative transfer module, absorption coefficient.
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
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "base/cs_assert.h"
#include "base/cs_base.h"
#include "comb/cs_coal.h"
#include "cogz/cs_combustion_gas.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "pprt/cs_combustion_model.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_time_step.h"

#include "rayt/cs_rad_transfer.h"
#include "rayt/cs_rad_transfer_modak.h"
#include "rayt/cs_rad_transfer_adf_models.h"
#include "rayt/cs_rad_transfer_fsck.h"
#include "rayt/cs_rad_transfer_rcfsk.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "rayt/cs_rad_transfer_absorption.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_rad_transfer_absorption.cpp
        Absorption coefficient computation for radiative transfer
        with specific physical models.

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
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
 * \brief Compute absorption coefficient for gas mix and particles for
 *        pulverized coal or other specific physical models.
 *
 * For the P-1 model, this function also checks whether the medium's optical
 * length is at least of the order of unity.
 *
 * \param[in]   tempk      gas phase temperature at cells (in Kelvin)
 * \param[out]  cpro_cak0  medium (gas) Absorption coefficient
 * \param[out]  kgas       radiation coefficients of the gray gases at cells
 *                         (per gas)
 * \param[out]  agas       weights of the gray gases at cells (per gas)
 * \param[out]  agasb      weights of the gray gases at boundary faces (per gas)
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_absorption(const cs_real_t  tempk[],
                           cs_real_t        cpro_cak0[],
                           cs_real_t        kgas[],
                           cs_real_t        agas[],
                           cs_real_t        agasb[])
{
  cs_real_t *w1 = nullptr, *w2 = nullptr, *w3 = nullptr;

  const cs_mesh_t *m = cs_glob_mesh;
  const int n_cells = m->n_cells;
  const int n_cells_ext = m->n_cells_with_ghosts;
  const int *pm_flag = cs_glob_physical_model_flag;

  cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;

  /* Initialization */

  if (   rt_params->imgrey >= 1
      || rt_params->imoadf >= 1
      || rt_params->imfsck >= 1) {
    CS_MALLOC(w1, n_cells_ext, cs_real_t);
    CS_MALLOC(w2, n_cells_ext, cs_real_t);
    CS_MALLOC(w3, n_cells_ext, cs_real_t);
  }

  cs_real_t *crom = CS_F_(rho)->val;

  /* Absorption coefficient of gas mix (m-1)
     --------------------------------------- */

  /* Gas combustion:
     - diffusion flame
     - premixed flame (EBU model) */

  if (   pm_flag[CS_COMBUSTION_3PT] >= 0
      || pm_flag[CS_COMBUSTION_SLFM] >= 0
      || pm_flag[CS_COMBUSTION_EBU] >= 0) {

    cs_combustion_gas_model_t  *cm = cs_glob_combustion_gas_model;

    if (rt_params->imgrey > 0) {

      const int n_gas_e = cm->n_gas_el_comp;
      const int n_gas_g = cm->n_gas_species;
      cs_real_t xpro;

      cs_real_t *_w;
      CS_MALLOC(_w, n_gas_e + n_gas_g + n_gas_e, cs_real_t);

      cs_real_t *xk = _w;
      cs_real_t *yi = _w + n_gas_e;
      cs_real_t *yk = _w + (n_gas_e + n_gas_g);

      const cs_real_t xsoot = cm->xsoot;
      const cs_real_t rosoot = cm->rosoot;

      const cs_real_t *cvar_fsm = nullptr;
      if (cm->isoot >= 1)
        cvar_fsm = CS_F_(fsm)->val;

      // If we are in multiphase, we get the first temperature field
      const cs_real_t *cpro_temp;
      if (CS_F_(t) != nullptr)
        cpro_temp = CS_F_(t)->val;
      else
        cpro_temp = CS_FI_(t, 0)->val;

      const cs_real_t *cpro_ym1 = cs_field_by_name("ym_fuel")->val;
      const cs_real_t *cpro_ym2 = cs_field_by_name("ym_oxyd")->val;
      const cs_real_t *cpro_ym3 = cs_field_by_name("ym_prod")->val;

      const double *restrict wmolg = cm->wmolg;

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

        cs_real_t xm = 1./ (  cpro_ym1[cell_id]/wmolg[0]
                            + cpro_ym2[cell_id]/wmolg[1]
                            + cpro_ym3[cell_id]/wmolg[2]);
        w1[cell_id] = cpro_ym3[cell_id]*xm/wmolg[2]*cm->xco2;
        w2[cell_id] = cpro_ym3[cell_id]*xm/wmolg[2]*cm->xh2o;

        if (cm->isoot >= 0) {   /* Soot model */
          cs_real_t ys;
          if (cm->isoot == 0 && cm->iic > 0)
            ys = cpro_ym3[cell_id]*cm->coefeg[2][cm->iic-1];
          else if (cm->isoot == 0)
            ys = xsoot * cpro_ym3[cell_id];
          else // (cm->isoot >= 1)
            ys = cvar_fsm[cell_id];

          w3[cell_id] = ys * crom[cell_id] / rosoot;
        }
        else
          w3[cell_id] = 0;

        if (rt_params->imgrey == 2) {
          yi[0] = cpro_ym1[cell_id];
          yi[1] = cpro_ym2[cell_id];
          yi[2] = cpro_ym3[cell_id];
          cs_assert(n_gas_g <= 3);  /* Otherwise fill values */

          cs_combustion_gas_yg2xye(yi, yk, xk);
          xpro = (xk[2] + xk[3]);
          cpro_cak0[cell_id]
            = 1225. * w3[cell_id] * cpro_temp[cell_id] + 0.1 * xpro;
        }
      }

      if (rt_params->imgrey == 1) {
        cs_rad_transfer_modak(cpro_cak0, w1, w2, w3, cpro_temp);
      }

      CS_FREE(_w);
    }
    else if (rt_params->imfsck == 2) {

      for (int gg_id = 0; gg_id < rt_params->nwsgg; gg_id++) {
        char f_name[64];
        snprintf(f_name, 63, "spectral_absorption_coeff_%2d", gg_id + 1);
        cs_field_t *f_kgabs = cs_field_by_name_try(f_name);

        if (f_kgabs != nullptr)
          for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
            kgas[n_cells*gg_id + cell_id] = f_kgabs->val[cell_id];
      }
    }
    else { /* if (rt_params->imgrey != 1) */
      const cs_real_t *cpro_ckabs = cs_field_by_name("kabs")->val;
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        cpro_cak0[cell_id] = cpro_ckabs[cell_id];
    }
  }

  else if (cs_glob_coal_model != nullptr) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    cs_real_t *cpro_temp = cs_field_by_name("temperature")->val;
    cs_real_t *cpro_yco2 = cs_field_by_name("ym_co2")->val;
    cs_real_t *cpro_yh2o = cs_field_by_name("ym_h2o")->val;
    cs_real_t *cpro_mmel = cs_field_by_name("xm")->val;

    const double *restrict wmole = cm->wmole;
    const int ico2 = cm->ico2 - 1;
    const int ih2o = cm->ih2o - 1;

    if (   rt_params->imgrey == 1
        || rt_params->imoadf >= 1
        || rt_params->imfsck >= 1) {

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        /* CO2 volume concentration */
        w1[cell_id] = cpro_mmel[cell_id]/wmole[ico2] * cpro_yco2[cell_id];
        /* H2O volume concentration */
        w2[cell_id] = cpro_mmel[cell_id]/wmole[ih2o] * cpro_yh2o[cell_id];
        /* soot volume fraction */
        w3[cell_id] = 0;
      }

      if (rt_params->imgrey == 1)
        cs_rad_transfer_modak(cpro_cak0, w1, w2, w3, cpro_temp);

      else if (rt_params->imoadf == 1)
        cs_rad_transfer_adf08(w1, w2, tempk, kgas, agas, agasb);

      else if (rt_params->imoadf == 2)
        cs_rad_transfer_adf50(w1, w2, tempk, kgas, agas, agasb);

      else if (rt_params->imfsck >= 1)
        cs_rad_transfer_fsck(w1, w2, tempk, kgas, agas, agasb);
    }

    else {
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        cpro_cak0[cell_id] = cm->ckabs0;
    }

  }

  /* Free temporary memory */

  CS_FREE(w1);
  CS_FREE(w2);
  CS_FREE(w3);

  /* Absorption coefficient of particles per class k2/x2 (m-1)
     --------------------------------------------------------- */

  /* Coal combustion */

  if (cs_glob_coal_model != nullptr) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    for (int icla = 0; icla < cm->nclacp; icla++) {

      char s[64];

      int icha = cm->ichcor[icla] - 1;

      snprintf(s, 63, "diam_p_%02d", icla+1); s[63] = '\0';
      cs_real_t *cpro_diam2 = cs_field_by_name(s)->val;

      snprintf(s, 63, "rho_p_%02d", icla+1); s[63] = '\0';
      cs_real_t *cpro_rom2 = cs_field_by_name(s)->val;

      cs_real_t *cpro_cak = CS_FI_(rad_cak, 1 + icla)->val;

      double  c0 = cm->xashch[icha] * cs_math_sq(cm->diam20[icla]);
      double  c1 = (1 - cm->xashch[icha]);

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

        /* particles diameter */
        double dd2 = sqrt(c0 + c1 * cs_math_sq(cpro_diam2[cell_id]));

        /* absorption coefficient */
        cpro_cak[cell_id] = 1.5 *crom[cell_id] / (cpro_rom2[cell_id]*dd2);

      }

    }

  }

  /* Electric arcs*/

  if (pm_flag[CS_ELECTRIC_ARCS] >= 0) {

    assert(rt_params->nrphas == 1);

    /* If the "radiation_source" field is not present,
       then assume we use directly the standard radiative absorption
       coefficient field, whose values are read by the electric model
       directly when reading "dp_elec" (or read no such values, in
       which case the default initialization to zero stands). */

    if (CS_F_(radsc) != nullptr) {

      cs_real_t *cpro_radsc = CS_F_(radsc)->val;

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        cpro_cak0[cell_id] = cpro_radsc[cell_id];

    }

  }

  /* Clip absorption coefficient with P-1 approximation
     -------------------------------------------------- */

  /* Standard check of absorption coefficient values. This coefficient should
     ensure an optical length of the order of unity.

     Do not use with ADF model !!! */

  if (rt_params->type == CS_RAD_TRANSFER_P1 && rt_params->imoadf == 0) {

    CS_MALLOC(w3, n_cells_ext, cs_real_t);

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      w3[cell_id] = cpro_cak0[cell_id];

    /* Absorption coefficient for a gas/coal particles mix */

    if (cs_glob_coal_model != nullptr) {

      cs_coal_model_t  *cm = cs_glob_coal_model;

      for (int icla = 0; icla < cm->nclacp; icla++) {

        char s[64];

        snprintf(s, 63, "x_p_%02d", icla+1); s[63] = '\0';
        cs_real_t *cpro_x_p = cs_field_by_name(s)->val;  /* property here */
        cs_real_t *cpro_cak = CS_FI_(rad_cak, 1 + icla)->val;

        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
          w3[cell_id] += (cpro_x_p[cell_id] * cpro_cak[cell_id]);

      }

    }

    cs_rad_transfer_absorption_check_p1(w3);

    CS_FREE(w3);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute absorption coefficient for the case of the RCFSK model.
 *
 * \param[in]   tempk      gas phase temperature at cells (in Kelvin)
 * \param[out]  kgas       radiation coefficients of the gray gases at cells
 *                         (per gas)
 * \param[out]  agas       weights of the gray gases at cells (per gas)
 * \param[out]  agasb      weights of the gray gases at boundary faces (per gas)
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_rcfsk_absorption(const cs_real_t  tempk[],
                                 cs_real_t        kgas[],
                                 cs_real_t        agas[],
                                 cs_real_t        agasb[])
{
  cs_real_t *w1 = nullptr, *w2 = nullptr, *w3 = nullptr;

  const cs_mesh_t *m           = cs_glob_mesh;
  const int        n_cells     = m->n_cells;
  const int        n_cells_ext = m->n_cells_with_ghosts;
  const int       *pm_flag     = cs_glob_physical_model_flag;

  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  /* Initialization */

  CS_MALLOC(w1, n_cells_ext, cs_real_t);
  CS_MALLOC(w2, n_cells_ext, cs_real_t);
  CS_MALLOC(w3, n_cells_ext, cs_real_t);

  cs_real_t *crom = CS_F_(rho)->val;

  const cs_real_t xsoot  = cm->xsoot;
  const cs_real_t rosoot = cm->rosoot;

  cs_real_t *cvar_fsm = nullptr;
  if (cm->isoot >= 1)
    cvar_fsm = CS_F_(fsm)->val;

  /* Absorption coefficient of gas mix (m-1)
     --------------------------------------- */

  /* Gas combustion:
     - diffusion flame */

  if (pm_flag[CS_COMBUSTION_3PT] >= 0) {

    const cs_real_t *cpro_ym1 = cs_field_by_name("ym_fuel")->val;
    const cs_real_t *cpro_ym2 = cs_field_by_name("ym_oxyd")->val;
    const cs_real_t *cpro_ym3 = cs_field_by_name("ym_prod")->val;

    const double *restrict wmolg = cm->wmolg;

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      cs_real_t xm
        = 1. / (  cpro_ym1[cell_id] / wmolg[0] + cpro_ym2[cell_id] / wmolg[1]
                + cpro_ym3[cell_id] / wmolg[2]);
      w1[cell_id] = cpro_ym3[cell_id] * xm / wmolg[2] * cm->xco2;
      w2[cell_id] = cpro_ym3[cell_id] * xm / wmolg[2] * cm->xh2o;

      if (cm->isoot >= 0) { /* Soot model */
        cs_real_t ys;
        if (cm->isoot == 0 && cm->iic > 0)
          ys = cpro_ym3[cell_id] * cm->coefeg[2][cm->iic - 1];
        else if (cm->isoot == 0)
          ys = xsoot * cpro_ym3[cell_id];
        else // if (cm->isoot >= 1)
          ys = cvar_fsm[cell_id];

        /* Calculation of soot volume fraction */
        w3[cell_id] = ys * crom[cell_id] / rosoot;
      }
      else
        w3[cell_id] = 0.;
    }

    /* Calculation of the absorption coefficient using the RCFSK scheme */
    cs_rad_transfer_rcfsk(w1, w2, w3, tempk, kgas, agas, agasb);

  }

  /* Free temporary memory */

  CS_FREE(w1);
  CS_FREE(w2);
  CS_FREE(w3);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the absorption validity fo the P-1 approximation.
 *
 * For the P-1 model, the medium's optical length should be at least of
 * the order of unity.
 *
 * \param[in]  cpro_cak  absorption coefficient values (at cells)
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_absorption_check_p1(const cs_real_t  cpro_cak[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_real_t *b_face_surf = cs_glob_mesh_quantities->b_face_surf;

  cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;

  cs_real_t s[2] = {0, 0};

  /* Compute the characteristic length of the computational domain */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
    s[0] += b_face_surf[face_id];

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    s[1] += cell_vol[cell_id];

  cs_parall_sum(2, CS_REAL_TYPE, s);

  cs_real_t xlc = 3.6 * s[1] / s[0];

  /* Clipping on ck  */

  cs_real_t xkmin = 1.0 / xlc;
  cs_gnum_t iok = 0;
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    if (cpro_cak[cell_id] < xkmin)
      iok++;
  }

  cs_parall_sum(1, CS_GNUM_TYPE, &iok);

  /* Warning if the optical thickness is too big   */

  cs_real_t pp = rt_params->xnp1mx / 100.0;
  if (iok > pp * cs_glob_mesh->n_g_cells) {

    if (   rt_params->iwrp1t < 2
        || cs_glob_time_step->nt_max < cs_glob_time_step->nt_cur - 2)
      bft_printf
        (_(" Warning: P-1 radiative model (in %s)\n"
           " --------\n"
           "   The optical length of the semi-transparent medium must be\n"
           "   at least of the order of unity to be in the application\n"
           "   domain of the P-1 approximation.\n"
           "   This does not seem to be the cas here.\n\n"
           "   The minimum absorption coefficient required to ensure\n"
           "   this optical length is xkmin = %11.4e.\n"
           "   This value is not reached for %11.4e%% of mesh cells.\n\n"
           "   The percentage of cells for which we allow this condition\n"
           "   not to be reached is currently set to:\n"
           "   \"cs_glob_rad_transfer_params->xnp1mx\" = %11.4e.\n\n"),
         __func__, xkmin, iok/cs_glob_mesh->n_g_cells*100.,
         rt_params->xnp1mx);
    rt_params->iwrp1t += 1;

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
