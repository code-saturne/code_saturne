/*============================================================================
 * Main functions dedicated to soil management in groundwater flows
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_array.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_param_types.h"
#include "cs_post.h"
#include "cs_reco.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_soil.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_gwf_soil.c

  \brief Main functions dedicated to soil management in groundwater flows when
         using CDO schemes

*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_GWF_SOIL_DBG  0

/*============================================================================
 * Structure definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_soil[] =
  " Stop execution. The structure related to a soil is empty.\n"
  " Please check your settings.\n";

static int  _n_soils = 0;
static cs_gwf_soil_t  **_soils = NULL;

/* The following array enables to get the soil id related to each cell.
   The array size is equal to n_cells */

static short int *_cell2soil_ids = NULL;

/* Array storing the dual volume of each vertex weighted by the soil
   porosity */

static double  *_dual_porous_volume = NULL;

/* Array storing the soil state */

static int  *_soil_state_array = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the (effective) liquid saturation and the liquid capacity
 *        (dSl/dPc) for a soil associated to a Van Genuchten model and a
 *        two-phase flow model.
 *
 * \param[in]  sp       set of modelling parameters
 * \param[in]  pc       capillarity pressure
 * \param[out] sl_e     effective liquid saturation
 * \param[out] sl       liquid saturation
 * \param[out] lcap     liquid capacity
 */
/*----------------------------------------------------------------------------*/

static inline void
_set_sle_lcap_vg(const cs_gwf_soil_vgm_tpf_param_t    *sp,
                 const double                          pc,
                 double                               *sl_e,
                 double                               *sl,
                 double                               *lcap)
{
  const double  pr_e = pc*sp->inv_pr_r;
  const double  pr_en = pow(pr_e, sp->n);
  const double  lcap_coef = sp->sl_range * sp->inv_pr_r * (1 - sp->n);
  const double  _sle = pow(pr_en + 1, -sp->m);

  *sl_e = _sle;
  *sl = sp->sl_r + sp->sl_range * _sle;
  *lcap = lcap_coef * pr_en/pr_e * _sle/(1 + pr_en);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the relative permeabilities for a soil associated to a Van
 *        Genuchten-Mualem model and a two-phase flow model.
 *
 * \param[in]  sp        set of modelling parameters
 * \param[in]  sl_e      effective liquid saturation
 * \param[out] krl       relative permeability for the liquid phase
 * \param[out] krg       relative permeability for the gas phase
 */
/*----------------------------------------------------------------------------*/

static inline void
_set_kr_vgm(const cs_gwf_soil_vgm_tpf_param_t    *sp,
            const double                          sl_e,
            double                               *krl,
            double                               *krg)
{
  const double  sl_e_coef = 1 - pow(sl_e, sp->inv_m);
  const double  krl_coef = 1 - pow(sl_e_coef, sp->m);

  *krl = sqrt(sl_e) * krl_coef * krl_coef;
  *krg = sqrt(1 - sl_e) * pow(sl_e_coef, 2*sp->m);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the relative permeabilities for a soil associated to a Van
 *        Genuchten-Mualem model and a two-phase flow model.
 *        Case of a joining function relying on a 2nd order polynomial for the
 *        relative permeability in the gas phase
 *
 * \param[in]  sp        set of modelling parameters
 * \param[in]  sl_e      effective liquid saturation
 * \param[out] krl       relative permeability for the liquid phase
 * \param[out] krg       relative permeability for the gas phase
 */
/*----------------------------------------------------------------------------*/

static inline void
_set_kr_vgm_poly2(const cs_gwf_soil_vgm_tpf_param_t    *sp,
                  const double                          sl_e,
                  double                               *krl,
                  double                               *krg)
{
  const double  sl_e_coef = 1 - pow(sl_e, sp->inv_m);
  const double  krl_coef = 1 - pow(sl_e_coef, sp->m);
  const double  delta_s = sl_e - sp->sle_thres;

  *krl = sqrt(sl_e) * krl_coef * krl_coef;
  *krg = sp->krg_alpha*delta_s*delta_s + sp->krg_beta*delta_s + sp->krg_star;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the joining parameters for a soil associated to a Van
 *        Genuchten-Mualem model.
 *
 * \param[in, out] sp         set of modelling parameters
 */
/*----------------------------------------------------------------------------*/

static void
_joining_param_vgm(cs_gwf_soil_vgm_tpf_param_t    *sp)
{
  const double  sle = sp->sle_thres, sle_conj = 1./(1 - sle);

  sp->pc_star = sp->pr_r * pow( pow(sle,-sp->inv_m)-1, 1./sp->n );

  const double  ratio = sp->pc_star*sp->inv_pr_r;
  const double  ratio_n = pow(ratio, sp->n);
  const double  ratio_coef = pow(ratio_n + 1, -sp->m-1);

  sp->dsldpc_star = sp->inv_pr_r * (1-sp->n) * ratio_coef * ratio_n/ratio;

  /* Joining coefficients */

  if (sp->sle_jtype == CS_GWF_SOIL_JOIN_C1_HYPERBOLIC) {
    sp->sle_alpha = sle_conj*sle_conj * sp->dsldpc_star;
    sp->sle_beta = sle_conj;
  }
  else if (sp->sle_jtype == CS_GWF_SOIL_JOIN_C1_EXPONENTIAL) {
    sp->sle_alpha = (1 - sle);
    sp->sle_beta = -sp->dsldpc_star*sle_conj;
  }

  if (sp->krg_jtype == CS_GWF_SOIL_JOIN_C1_POLY_ORDER2) {

    const double  sle_coef = 1 - pow(sle, sp->inv_m);

    sp->krg_star = sqrt(1 - sle) * pow(sle_coef, 2*sp->m);
    sp->dkrgdsl_star = -sqrt(1 - sle) * pow(sle_coef, 2*sp->m-1) *
      ( 0.5*sle_coef*sle_conj + 2*pow(sle, (1-sp->m)*sp->inv_m) );

    sp->krg_beta = 1/sle*sp->dkrgdsl_star;
    sp->krg_alpha = sle_conj * ( sp->krg_star*sle_conj - sp->krg_beta );

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function that compute the new values of the properties related to
 *        a soil with a Van Genuchten-Mualem.
 *        Case of an isotropic permeability and an unsteady Richards eq.
 *
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone          pointer to a cs_zone_t
 * \param[in, out] soil          pointer to a soil structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_soil_vgm_spf(const cs_real_t              t_eval,
                         const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *cdoq,
                         const cs_zone_t             *zone,
                         cs_gwf_soil_t               *soil)
{
  CS_NO_WARN_IF_UNUSED(t_eval);
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(cdoq);

  if (soil == NULL)
    return;

  assert(soil->hydraulic_model ==  CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE);

  /* Retrieve the soil parameters */

  cs_gwf_soil_vgm_spf_param_t  *sp = soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_uspf_t  *hc = soil->hydraulic_context;

  /* Only isotropic values are considered in this case */

  const double  iso_satval = soil->abs_permeability[0][0];
  const double  delta_m = soil->porosity - sp->residual_moisture;
  const cs_real_t  *head = hc->head_in_law;

  /* Retrieve field values associated to properties to update */

  cs_real_t  *permeability = hc->permeability_field->val;
  cs_real_t  *moisture = hc->moisture_field->val;
  cs_real_t  *capacity = hc->capacity_field->val;

  assert(capacity != NULL && permeability != NULL && moisture != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)                \
  shared(head, zone, sp, permeability, moisture, capacity)              \
  firstprivate(iso_satval, delta_m)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];
    const cs_real_t  h = head[c_id];

    if (h < 0) { /* S_e(h) = [1 + |alpha*h|^n]^(-m) */

      const double  coef = pow(fabs(sp->scale * h), sp->n);
      const double  se = pow(1 + coef, -sp->m);
      const double  se_pow_overm = pow(se, 1/sp->m);
      const double  coef_base = 1 - pow(1 - se_pow_overm, sp->m);

      /* Set the permeability value : abs_perm * rel_perm */

      permeability[c_id] =
        iso_satval * pow(se, sp->tortuosity) * coef_base*coef_base;

      /* Set the moisture content (or liquid saturation) */

      moisture[c_id] = se*delta_m + sp->residual_moisture;

      /* Set the soil capacity = \frac{\partial S_l}{partial h} */

      const double  ccoef = -sp->n * sp->m * delta_m;
      const double  se_m1 = se/(1. + coef);

      capacity[c_id] = ccoef * coef/h * se_m1;

    }
    else {

      /* Set the permeability value to the saturated values */

      permeability[c_id] = iso_satval;

      /* Set the moisture content (Sle = 1 in this case)*/

      moisture[c_id] = delta_m + sp->residual_moisture;

      /* Set the soil capacity */

      capacity[c_id] = 0.;

    }

  } /* Loop on selected cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new values of the properties related to a soil with a Van
 *        Genuchten-Mualem.
 *
 *        Case relying on the interpolation of the capillarity pressure at cell
 *        centers.
 *        Case of an isotropic permeability and a two-phase flow in porous
 *        media.
 *
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone          pointer to a cs_zone_t
 * \param[in, out] soil          pointer to a soil structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_soil_vgm_tpf_pc(const cs_real_t              t_eval,
                            const cs_mesh_t             *mesh,
                            const cs_cdo_connect_t      *connect,
                            const cs_cdo_quantities_t   *cdoq,
                            const cs_zone_t             *zone,
                            cs_gwf_soil_t               *soil)
{
  CS_NO_WARN_IF_UNUSED(t_eval);
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(cdoq);

  if (soil == NULL)
    return;

  /* Retrieve the soil parameters */

  cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_tpf_t  *hc = soil->hydraulic_context;

  /* Retrieve arrays to update */

  cs_real_t  *sliq = hc->l_saturation->val;
  cs_real_t  *lcap = hc->l_capacity;
  cs_real_t  *krl = hc->l_rel_permeability;
  cs_real_t  *krg = hc->g_rel_permeability;

  assert(sliq != NULL && lcap != NULL && krl != NULL && krg != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];
    const cs_real_t  pc = hc->c_pressure_cells[c_id];

    if (pc > 0) {

      double  sl_e;

      _set_sle_lcap_vg(sp, pc, &sl_e, &(sliq[c_id]), &(lcap[c_id]));
      _set_kr_vgm(sp, sl_e, &(krl[c_id]), &(krg[c_id]));

    }
    else { /* Saturated case */

      sliq[c_id] = sp->sl_s;
      lcap[c_id] = 0.;
      krl[c_id] = 1;
      krg[c_id] = 0.;

    }

  } /* Loop on selected cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new values of the properties related to a soil with a Van
 *        Genuchten-Mualem.
 *
 *        Case relying on the interpolation of the liquid saturation at cell
 *        centers.
 *        Case of an isotropic permeability and a two-phase flow in porous
 *        media.
 *
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone          pointer to a cs_zone_t
 * \param[in, out] soil          pointer to a soil structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_soil_vgm_tpf_sl(const cs_real_t              t_eval,
                            const cs_mesh_t             *mesh,
                            const cs_cdo_connect_t      *connect,
                            const cs_cdo_quantities_t   *cdoq,
                            const cs_zone_t             *zone,
                            cs_gwf_soil_t               *soil)
{
  CS_NO_WARN_IF_UNUSED(t_eval);
  CS_NO_WARN_IF_UNUSED(mesh);

  if (soil == NULL)
    return;

  const cs_adjacency_t  *c2v = connect->c2v;

  /* Retrieve the soil parameters */

  cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_tpf_t  *hc = soil->hydraulic_context;

  const cs_real_t  *pc_val = hc->c_pressure->val;

  /* Retrieve arrays to update */

  cs_real_t  *sliq = hc->l_saturation->val;
  cs_real_t  *lcap = hc->l_capacity;
  cs_real_t  *krl = hc->l_rel_permeability;
  cs_real_t  *krg = hc->g_rel_permeability;

  assert(sliq != NULL && lcap != NULL && krl != NULL && krg != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];

    double  sliq_csum = 0.;
    double  lcap_csum = 0.;
    double  krl_csum = 0.;
    double  krg_csum = 0.;

    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

      /* capillarity pressure at vertices */

      const cs_real_t  pc = pc_val[c2v->ids[j]];
      const cs_real_t  pvc = cdoq->pvol_vc[j];

      if (pc > 0) {

        double  sle_v, sliq_v, lcap_v, krl_v, krg_v;

        _set_sle_lcap_vg(sp, pc, &sle_v, &sliq_v, &lcap_v);
        _set_kr_vgm(sp, sle_v, &krl_v, &krg_v);

        sliq_csum += pvc*sliq_v;
        lcap_csum += pvc*lcap_v;
        krl_csum += pvc*krl_v;
        krg_csum += pvc*krg_v;

      }
      else { /* Saturated case */

        sliq_csum += pvc*sp->sl_s;
        lcap_csum += 0.;
        krl_csum += pvc;
        krg_csum += 0;

      }

    } /* Loop on cell vertices */

    /* Compute the cell values */

    const double  inv_volc = 1./cdoq->cell_vol[c_id];

    sliq[c_id] = inv_volc * sliq_csum;
    lcap[c_id] = inv_volc * lcap_csum;
    krl[c_id] = inv_volc * krl_csum;
    krg[c_id] = inv_volc * krg_csum;

  } /* Loop on selected cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new property values related to a soil with a Van
 *        Genuchten-Mualem.
 *
 *        C1 join relying on a hyperbolic function
 *        Case relying on the interpolation of the capillarity pressure at cell
 *        centers.
 *        Case of an isotropic permeability and a two-phase flow in porous
 *        media.
 *
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone          pointer to a cs_zone_t
 * \param[in, out] soil          pointer to a soil structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_soil_vgm_tpf_c1_hyperbolic_pc(const cs_real_t              t_eval,
                                          const cs_mesh_t             *mesh,
                                          const cs_cdo_connect_t      *connect,
                                          const cs_cdo_quantities_t   *cdoq,
                                          const cs_zone_t             *zone,
                                          cs_gwf_soil_t               *soil)
{
  CS_NO_WARN_IF_UNUSED(t_eval);
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(cdoq);

  if (soil == NULL)
    return;

  /* Retrieve the soil parameters */

  cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_tpf_t  *hc = soil->hydraulic_context;

  /* Retrieve arrays to update */

  cs_real_t  *sliq = hc->l_saturation->val;
  cs_real_t  *lcap = hc->l_capacity;
  cs_real_t  *krl = hc->l_rel_permeability;
  cs_real_t  *krg = hc->g_rel_permeability;

  assert(sliq != NULL && lcap != NULL && krl != NULL && krg != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    double  sl_e;

    const cs_lnum_t  c_id = zone->elt_ids[i];
    const cs_real_t  pc = hc->c_pressure_cells[c_id];

    if (pc < sp->pc_star) { /* Joining function */

      const double  denum = sp->sle_alpha*(pc - sp->pc_star) + sp->sle_beta;

      sl_e = 1 - 1/denum;
      sliq[c_id] = sp->sl_r + sp->sl_range*sl_e;
      lcap[c_id] = sp->sle_alpha / (denum*denum);

    }
    else { /* pc >= pc_star */

      _set_sle_lcap_vg(sp, pc, &sl_e, &(sliq[c_id]), &(lcap[c_id]));

    }

    _set_kr_vgm(sp, sl_e, &(krl[c_id]), &(krg[c_id]));

  } /* Loop on selected cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new property values related to a soil with a Van
 *        Genuchten-Mualem.
 *
 *        C1 join relying on a hyperbolic function
 *        Case relying on the interpolation of the liquid saturation at cell
 *        centers.
 *        Case of an isotropic permeability and a two-phase flow in porous
 *        media.
 *
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone          pointer to a cs_zone_t
 * \param[in, out] soil          pointer to a soil structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_soil_vgm_tpf_c1_hyperbolic_sl(const cs_real_t              t_eval,
                                          const cs_mesh_t             *mesh,
                                          const cs_cdo_connect_t      *connect,
                                          const cs_cdo_quantities_t   *cdoq,
                                          const cs_zone_t             *zone,
                                          cs_gwf_soil_t               *soil)
{
  CS_NO_WARN_IF_UNUSED(t_eval);
  CS_NO_WARN_IF_UNUSED(mesh);

  if (soil == NULL)
    return;

  const cs_adjacency_t  *c2v = connect->c2v;

  /* Retrieve the soil parameters */

  cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_tpf_t  *hc = soil->hydraulic_context;

  const cs_real_t  *pc_val = hc->c_pressure->val;

  /* Retrieve arrays to update */

  cs_real_t  *sliq = hc->l_saturation->val;
  cs_real_t  *lcap = hc->l_capacity;
  cs_real_t  *krl = hc->l_rel_permeability;
  cs_real_t  *krg = hc->g_rel_permeability;

  assert(sliq != NULL && lcap != NULL && krl != NULL && krg != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];

    double  sliq_csum = 0.;
    double  lcap_csum = 0.;
    double  krl_csum = 0.;
    double  krg_csum = 0.;

    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

      double  sle_v, sliq_v, lcap_v, krl_v, krg_v;

      /* capillarity pressure at vertices */

      const cs_real_t  pc = pc_val[c2v->ids[j]];
      const cs_real_t  pvc = cdoq->pvol_vc[j];

      if (pc < sp->pc_star) { /* Joining function */

        const double  denum = sp->sle_alpha*(pc - sp->pc_star) + sp->sle_beta;

        sle_v = 1 - 1/denum;
        sliq_v = sp->sl_r + sp->sl_range*sle_v;
        lcap_v = sp->sle_alpha / (denum*denum);

      }
      else { /* pc >= pc_star */

        _set_sle_lcap_vg(sp, pc, &sle_v, &sliq_v, &lcap_v);

      }

      _set_kr_vgm(sp, sle_v, &krl_v, &krg_v);

      sliq_csum += pvc * sliq_v;
      lcap_csum += pvc * lcap_v;
      krl_csum += pvc * krl_v;
      krg_csum += pvc * krg_v;

    } /* Loop on cell vertices */

    /* Compute the cell values */

    const double  inv_volc = 1./cdoq->cell_vol[c_id];

    sliq[c_id] = inv_volc * sliq_csum;
    lcap[c_id] = inv_volc * lcap_csum;
    krl[c_id] = inv_volc * krl_csum;
    krg[c_id] = inv_volc * krg_csum;

  } /* Loop on selected cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new property values related to a soil with a Van
 *        Genuchten-Mualem.
 *
 *        C1 join relying on a hyperbolic function for Sl
 *        C1 join relying on a 2nd order polynomial function for krg
 *        Case relying on the interpolation of the liquid saturation at cell
 *        centers.
 *        Case of an isotropic permeability and a two-phase flow in porous
 *        media.
 *
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone          pointer to a cs_zone_t
 * \param[in, out] soil          pointer to a soil structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_soil_vgm_tpf_c1_hyperbolic_sl_krg(const cs_real_t           t_eval,
                                              const cs_mesh_t          *mesh,
                                              const cs_cdo_connect_t   *connect,
                                              const cs_cdo_quantities_t   *cdoq,
                                              const cs_zone_t          *zone,
                                              cs_gwf_soil_t            *soil)
{
  CS_NO_WARN_IF_UNUSED(t_eval);
  CS_NO_WARN_IF_UNUSED(mesh);

  if (soil == NULL)
    return;

  const cs_adjacency_t  *c2v = connect->c2v;

  /* Retrieve the soil parameters */

  cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_tpf_t  *hc = soil->hydraulic_context;

  const cs_real_t  *pc_val = hc->c_pressure->val;

  /* Retrieve arrays to update */

  cs_real_t  *sliq = hc->l_saturation->val;
  cs_real_t  *lcap = hc->l_capacity;
  cs_real_t  *krl = hc->l_rel_permeability;
  cs_real_t  *krg = hc->g_rel_permeability;

  assert(sliq != NULL && lcap != NULL && krl != NULL && krg != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];

    double  sliq_csum = 0.;
    double  lcap_csum = 0.;
    double  krl_csum = 0.;
    double  krg_csum = 0.;

    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

      double  sle_v, sliq_v, lcap_v, krl_v, krg_v;

      /* capillarity pressure at vertices */

      const cs_real_t  pc = pc_val[c2v->ids[j]];
      const cs_real_t  pvc = cdoq->pvol_vc[j];

      if (pc < sp->pc_star) { /* Joining function */

        const double  denum = sp->sle_alpha*(pc - sp->pc_star) + sp->sle_beta;

        sle_v = 1 - 1/denum;
        sliq_v = sp->sl_r + sp->sl_range*sle_v;
        lcap_v = sp->sle_alpha / (denum*denum);

        _set_kr_vgm_poly2(sp, sle_v, &krl_v, &krg_v);

      }
      else { /* pc >= pc_star */

        _set_sle_lcap_vg(sp, pc, &sle_v, &sliq_v, &lcap_v);
        _set_kr_vgm(sp, sle_v, &krl_v, &krg_v);

      }

      sliq_csum += pvc * sliq_v;
      lcap_csum += pvc * lcap_v;
      krl_csum += pvc * krl_v;
      krg_csum += pvc * krg_v;

    } /* Loop on cell vertices */

    /* Compute the cell values */

    const double  inv_volc = 1./cdoq->cell_vol[c_id];

    sliq[c_id] = inv_volc * sliq_csum;
    lcap[c_id] = inv_volc * lcap_csum;
    krl[c_id] = inv_volc * krl_csum;
    krg[c_id] = inv_volc * krg_csum;

  } /* Loop on selected cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new property values related to a soil with a Van
 *        Genuchten-Mualem.
 *
 *        C1 join relying on an exponential function.
 *        Case relying on the interpolation of the capillarity pressure at cell
 *        centers.
 *        Case of an isotropic permeability and a two-phase flow in porous
 *        media.
 *
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone          pointer to a cs_zone_t
 * \param[in, out] soil          pointer to a soil structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_soil_vgm_tpf_c1_exponential_pc(const cs_real_t              t_eval,
                                           const cs_mesh_t             *mesh,
                                           const cs_cdo_connect_t      *connect,
                                           const cs_cdo_quantities_t   *cdoq,
                                           const cs_zone_t             *zone,
                                           cs_gwf_soil_t               *soil)
{
  CS_NO_WARN_IF_UNUSED(t_eval);
  CS_NO_WARN_IF_UNUSED(mesh);
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(cdoq);

  if (soil == NULL)
    return;

  /* Retrieve the soil parameters */

  cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_tpf_t  *hc = soil->hydraulic_context;

  /* Retrieve arrays to update */

  cs_real_t  *sliq = hc->l_saturation->val;
  cs_real_t  *lcap = hc->l_capacity;
  cs_real_t  *krl = hc->l_rel_permeability;
  cs_real_t  *krg = hc->g_rel_permeability;

  assert(sliq != NULL && lcap != NULL && krl != NULL && krg != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    double  sl_e;

    const cs_lnum_t  c_id = zone->elt_ids[i];
    const cs_real_t  pc = hc->c_pressure_cells[c_id];

    if (pc < sp->pc_star) { /* Joining function */

      const double  dpc = pc - sp->pc_star;

      sl_e = 1 - sp->sle_alpha*exp(sp->sle_beta*dpc);
      sliq[c_id] = sp->sl_r + sp->sl_range*sl_e;
      lcap[c_id] = sp->dsldpc_star * exp(sp->sle_beta*dpc);

    }
    else { /* pc >= pc_star */

      _set_sle_lcap_vg(sp, pc, &sl_e, &(sliq[c_id]), &(lcap[c_id]));

    }

    _set_kr_vgm(sp, sl_e, &(krl[c_id]), &(krg[c_id]));

  } /* Loop on selected cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new property values related to a soil with a Van
 *        Genuchten-Mualem.
 *
 *        C1 join relying on an exponential function for Sl
 *        Case relying on the interpolation of the liquid saturation at cell
 *        centers.
 *        Case of an isotropic permeability and a two-phase flow in porous
 *        media.
 *
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone          pointer to a cs_zone_t
 * \param[in, out] soil          pointer to a soil structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_soil_vgm_tpf_c1_exponential_sl(const cs_real_t             t_eval,
                                           const cs_mesh_t            *mesh,
                                           const cs_cdo_connect_t     *connect,
                                           const cs_cdo_quantities_t  *cdoq,
                                           const cs_zone_t            *zone,
                                           cs_gwf_soil_t              *soil)
{
  CS_NO_WARN_IF_UNUSED(t_eval);
  CS_NO_WARN_IF_UNUSED(mesh);

  if (soil == NULL)
    return;

  const cs_adjacency_t  *c2v = connect->c2v;

  /* Retrieve the soil parameters */

  cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_tpf_t  *hc = soil->hydraulic_context;

  const cs_real_t  *pc_val = hc->c_pressure->val;

  /* Retrieve arrays to update */

  cs_real_t  *sliq = hc->l_saturation->val;
  cs_real_t  *lcap = hc->l_capacity;
  cs_real_t  *krl = hc->l_rel_permeability;
  cs_real_t  *krg = hc->g_rel_permeability;

  assert(sliq != NULL && lcap != NULL && krl != NULL && krg != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];

    double  sliq_csum = 0.;
    double  lcap_csum = 0.;
    double  krl_csum = 0.;
    double  krg_csum = 0.;

    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

      double  sle_v, sliq_v, lcap_v, krl_v, krg_v;

      /* capillarity pressure at vertices */

      const cs_real_t  pc = pc_val[c2v->ids[j]];
      const cs_real_t  pvc = cdoq->pvol_vc[j];

      if (pc < sp->pc_star) { /* Joining function */

        const double  dpc = pc - sp->pc_star;

        sle_v = 1 - sp->sle_alpha*exp(sp->sle_beta*dpc);
        sliq_v = sp->sl_r + sp->sl_range*sle_v;
        lcap_v = sp->dsldpc_star * exp(sp->sle_beta*dpc);

      }
      else { /* pc >= pc_star */

        _set_sle_lcap_vg(sp, pc, &sle_v, &sliq_v, &lcap_v);

      }

      _set_kr_vgm(sp, sle_v, &krl_v, &krg_v);

      sliq_csum += pvc * sliq_v;
      lcap_csum += pvc * lcap_v;
      krl_csum += pvc * krl_v;
      krg_csum += pvc * krg_v;

    } /* Loop on cell vertices */

    /* Compute the cell values */

    const double  inv_volc = 1./cdoq->cell_vol[c_id];

    sliq[c_id] = inv_volc * sliq_csum;
    lcap[c_id] = inv_volc * lcap_csum;
    krl[c_id] = inv_volc * krl_csum;
    krg[c_id] = inv_volc * krg_csum;

  } /* Loop on selected cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new property values related to a soil with a Van
 *        Genuchten-Mualem.
 *
 *        C1 join relying on an exponential function for Sl
 *        C1 join relying on a 2nd order polynomial function for krg
 *        Case relying on the interpolation of the liquid saturation at cell
 *        centers.
 *        Case of an isotropic permeability and a two-phase flow in porous
 *        media.
 *
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone          pointer to a cs_zone_t
 * \param[in, out] soil          pointer to a soil structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_soil_vgm_tpf_c1_exponential_sl_krg(const cs_real_t          t_eval,
                                               const cs_mesh_t         *mesh,
                                               const cs_cdo_connect_t  *connect,
                                               const cs_cdo_quantities_t  *cdoq,
                                               const cs_zone_t         *zone,
                                               cs_gwf_soil_t           *soil)
{
  CS_NO_WARN_IF_UNUSED(t_eval);
  CS_NO_WARN_IF_UNUSED(mesh);

  if (soil == NULL)
    return;

  const cs_adjacency_t  *c2v = connect->c2v;

  /* Retrieve the soil parameters */

  cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_tpf_t  *hc = soil->hydraulic_context;

  const cs_real_t  *pc_val = hc->c_pressure->val;

  /* Retrieve arrays to update */

  cs_real_t  *sliq = hc->l_saturation->val;
  cs_real_t  *lcap = hc->l_capacity;
  cs_real_t  *krl = hc->l_rel_permeability;
  cs_real_t  *krg = hc->g_rel_permeability;

  assert(sliq != NULL && lcap != NULL && krl != NULL && krg != NULL);

  /* Main loop on cells belonging to this soil */

# pragma omp parallel for if (zone->n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];

    double  sliq_csum = 0.;
    double  lcap_csum = 0.;
    double  krl_csum = 0.;
    double  krg_csum = 0.;

    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

      double  sle_v, sliq_v, lcap_v, krl_v, krg_v;

      /* capillarity pressure at vertices */

      const cs_real_t  pc = pc_val[c2v->ids[j]];
      const cs_real_t  pvc = cdoq->pvol_vc[j];

      if (pc < sp->pc_star) { /* Joining function */

        const double  dpc = pc - sp->pc_star;

        sle_v = 1 - sp->sle_alpha*exp(sp->sle_beta*dpc);
        sliq_v = sp->sl_r + sp->sl_range*sle_v;
        lcap_v = sp->dsldpc_star * exp(sp->sle_beta*dpc);

        _set_kr_vgm_poly2(sp, sle_v, &krl_v, &krg_v);

      }
      else { /* pc >= pc_star */

        _set_sle_lcap_vg(sp, pc, &sle_v, &sliq_v, &lcap_v);
        _set_kr_vgm(sp, sle_v, &krl_v, &krg_v);

      }

      sliq_csum += pvc * sliq_v;
      lcap_csum += pvc * lcap_v;
      krl_csum += pvc * krl_v;
      krg_csum += pvc * krg_v;

    } /* Loop on cell vertices */

    /* Compute the cell values */

    const double  inv_volc = 1./cdoq->cell_vol[c_id];

    sliq[c_id] = inv_volc * sliq_csum;
    lcap[c_id] = inv_volc * lcap_csum;
    krl[c_id] = inv_volc * krl_csum;
    krg[c_id] = inv_volc * krg_csum;

  } /* Loop on selected cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build an array storing the associated soil for each cell
 *         The lifecycle of this array is managed by the code.
 *
 * \param[in] n_cells      number of cells
 */
/*----------------------------------------------------------------------------*/

static void
_build_cell2soil(cs_lnum_t    n_cells)
{
  BFT_MALLOC(_cell2soil_ids, n_cells, short int);

  if (_n_soils == 1)
    memset(_cell2soil_ids, 0, sizeof(short int)*n_cells);

  else {

    assert(_n_soils > 1);
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t j = 0; j < n_cells; j++)
      _cell2soil_ids[j] = -1; /* unset by default */

    for (int soil_id = 0; soil_id < _n_soils; soil_id++) {

      const cs_gwf_soil_t  *soil = _soils[soil_id];
      const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

      assert(z != NULL);

#     pragma omp parallel for if (z->n_elts > CS_THR_MIN)
      for (cs_lnum_t j = 0; j < z->n_elts; j++)
        _cell2soil_ids[z->elt_ids[j]] = soil_id;

    } /* Loop on soils */

    /* Check if every cells is associated to a soil */

    for (cs_lnum_t j = 0; j < n_cells; j++)
      if (_cell2soil_ids[j] == -1)
        bft_error(__FILE__, __LINE__, 0,
                  " %s: At least cell %ld has no related soil.\n",
                  __func__, (long)j);

  } /* n_soils > 1 */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that at least one soil has been defined and the model of soil
 *        exists.
 *        Raise an error if a problem is encoutered.
 */
/*----------------------------------------------------------------------------*/

static void
_check_soil_settings(void)
{
  if (_n_soils < 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Groundwater module is activated but no soil is defined.",
              __func__);
  if (_soils == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The soil structure is not allocated whereas %d soils"
              " have been added.\n", __func__, _n_soils);

  for (int i = 0; i < _n_soils; i++) {

    const cs_zone_t  *z = cs_volume_zone_by_id(_soils[i]->zone_id);
    assert(z != NULL);

    if (_soils[i]->model == CS_GWF_SOIL_N_HYDRAULIC_MODELS)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid model of soil attached to zone %s\n",
                __func__, z->name);

    if (z->n_g_elts < 1) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(" %s: The soil \"%s\" is defined but associated to no cell.\n"
                 " Please check your settings.\n",
                 __func__, z->name);
    }

    if (z->n_elts > 0)
      if (z->elt_ids == NULL)
        bft_error(__FILE__, __LINE__, 0,
                  " %s: One assumes that z->elt_ids != NULL.\n"
                  " This is not the case for the soil \"%s\"\n",
                  __func__, z->name);

  } /* Loop on soils */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the number of allocated soils
 *
 * \return the number of allocated soils
 */
/*----------------------------------------------------------------------------*/

int
cs_gwf_get_n_soils(void)
{
  return _n_soils;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new cs_gwf_soil_t structure and add it to the array of
 *        soils. An initialization by default of all members is performed.
 *
 * \param[in] zone                pointer to a volume zone structure
 * \param[in] hydraulic_model     main hydraulic model for the module
 * \param[in] model               type of model for the soil behavior
 * \param[in] perm_type           type of permeability (iso/anisotropic)
 * \param[in] k_abs               absolute (intrisic) permeability
 * \param[in] porosity            porosity or max. moisture content
 * \param[in] bulk_density        value of the mass density
 * \param[in] hydraulic_context   pointer to the context structure
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_create(const cs_zone_t                 *zone,
                   cs_gwf_model_type_t              hydraulic_model,
                   cs_gwf_soil_model_t              model,
                   cs_property_type_t               perm_type,
                   double                           k_abs[3][3],
                   double                           porosity,
                   double                           bulk_density,
                   void                            *hydraulic_context)
{
  cs_gwf_soil_t  *soil = NULL;

  BFT_MALLOC(soil, 1, cs_gwf_soil_t);

  soil->id = _n_soils;

  /* Attached a volume zone to the current soil */

  assert(zone != NULL);
  soil->zone_id = zone->id;

  /* Members related to the hydraulic model */

  soil->hydraulic_model = hydraulic_model;
  soil->hydraulic_context = hydraulic_context;

  /* Members related to the soil parameters/model */

  soil->model = model;
  soil->model_param = NULL;

  soil->bulk_density = bulk_density;
  soil->porosity = porosity;

  for (int ki = 0; ki < 3; ki++)
    for (int kj = 0; kj < 3; kj++)
      soil->abs_permeability[ki][kj] = k_abs[ki][kj];

  switch (perm_type) {

  case CS_PROPERTY_ISO:
    soil->abs_permeability_dim = 1;
    break;

  case CS_PROPERTY_ANISO:
    soil->abs_permeability_dim = 9;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of absolute permeability.\n", __func__);

  } /* Switch on the type of absolute permeability */

  /* Initialize function pointers */

  soil->update_properties = NULL;
  soil->free_model_param = NULL;

  /* Initialization which are specific to a soil model */

  switch (model) {

  case CS_GWF_SOIL_SATURATED:
    if (hydraulic_model != CS_GWF_MODEL_SATURATED_SINGLE_PHASE)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid type of soil with the general hydraulic model.\n"
                " In a saturated single-phase model, all soils have to be"
                " of type CS_GWF_SOIL_SATURATED.\n", __func__);
    break;

  case CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_SINGLE_PHASE:
    {
      cs_gwf_soil_vgm_spf_param_t  *sp = NULL;

      BFT_MALLOC(sp, 1, cs_gwf_soil_vgm_spf_param_t);

      sp->residual_moisture = 0.;

      sp->n = 1.25;
      sp->m = 1 - 1./sp->n;
      sp->scale = 1.;
      sp->tortuosity = 1.;

      soil->model_param = sp;

      if (perm_type & CS_PROPERTY_ISO)
        if (hydraulic_model == CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE)
          soil->update_properties = _update_iso_soil_vgm_spf;
        else
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Invalid type of hydraulic model.\n"
                    " Please check your settings.", __func__);
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of property for the permeability.\n"
                  " Please check your settings.", __func__);
    }
    break;

  case CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_TWO_PHASE:
    {
      cs_gwf_soil_vgm_tpf_param_t  *sp = NULL;

      BFT_MALLOC(sp, 1, cs_gwf_soil_vgm_tpf_param_t);
      soil->model_param = sp;

      /* Default values */

      cs_gwf_soil_set_vgm_tpf_param(soil,
                                    1.7,  /* n */
                                    1e6,  /* pr_r */
                                    0,    /* sl_r */
                                    1);   /* sl_s */
    }
    break;

  case CS_GWF_SOIL_USER:
    break; /* All has to be done by the user */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of soil model\n", __func__);
    break; /* Nothing to do */

  } /* Switch on the soil model */

  /* Store the new soils in the soil array */

  _n_soils++;
  BFT_REALLOC(_soils, _n_soils, cs_gwf_soil_t *);
  _soils[soil->id] = soil;

  return soil;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a soil structure from its id
 *
 * \param[in]  id      id to look for
 *
 * \return a pointer to a cs_gwf_soil_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_by_id(int   id)
{
  if (id > -1 && id < _n_soils)
    return _soils[id];
  else
    return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a soil structure from its name
 *
 * \param[in]  name      name to look for
 *
 * \return a pointer to a cs_gwf_soil_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_by_name(const char    *name)
{
  if (name == NULL)
    return NULL;

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *s = _soils[i];
    const cs_zone_t  *zone = cs_volume_zone_by_id(s->zone_id);

    if (strcmp(zone->name, name) == 0)
      return s;
  }

  /* Not found among the list */

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_free_all(void)
{
  if (_n_soils < 1)
    return;
  assert(_soils != NULL);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    if (soil->free_model_param != NULL)
      soil->free_model_param(&(soil->model_param));

    if (soil->model_param != NULL) {

      switch (soil->model) {

      case CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_SINGLE_PHASE:
        {
          cs_gwf_soil_vgm_spf_param_t  *sp = soil->model_param;

          BFT_FREE(sp);
          sp = NULL;
        }
        break;

      case CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_TWO_PHASE:
        {
          cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

          BFT_FREE(sp);
          sp = NULL;
        }
        break;

      default:
        cs_base_warn(__FILE__, __LINE__);
        bft_printf("%s: The context structure of a soil may not be freed.\n",
                   __func__);
        break;

      } /* Switch on predefined soil context */

    }

    /* The hydraulic context is shared and thus is freed during the free of the
       cs_gwf_t structure */

    BFT_FREE(soil);

  } /* Loop on soils */

  BFT_FREE(_soils);
  BFT_FREE(_cell2soil_ids);
  BFT_FREE(_dual_porous_volume);
  BFT_FREE(_soil_state_array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last initialization step for the soil structures/parameters
 *
 * \param[in] gwf_model      modelling used for the GWF module
 * \param[in] post_flag      which post-processing to do
 * \param[in] n_cells        number of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_finalize_setup(cs_gwf_model_type_t    gwf_model,
                           cs_flag_t              post_flag,
                           cs_lnum_t              n_cells)
{
  /* Check the settings */

  _check_soil_settings();

  /* Store the soil id for each cell */

  _build_cell2soil(n_cells);

  /* Allocate if needed the soil state */

  if ((post_flag & CS_GWF_POST_SOIL_STATE) ||
      (gwf_model != CS_GWF_MODEL_SATURATED_SINGLE_PHASE)) {

    BFT_MALLOC(_soil_state_array, n_cells, int);

    /* Default initialization */

    for (cs_lnum_t i = 0; i < n_cells; i++)
      _soil_state_array[i] = CS_GWF_SOIL_STATE_SATURATED;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the settings related to all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_log_setup(void)
{
  cs_log_printf(CS_LOG_SETUP, "  * GWF | Number of soils: %d\n", _n_soils);

  char  id[64];
  for (int i = 0; i < _n_soils; i++) {

    const cs_gwf_soil_t  *soil = _soils[i];
    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    sprintf(id, "        Soil.%d |", soil->id);

    cs_log_printf(CS_LOG_SETUP, "\n%s Zone: %s\n", id, z->name);
    cs_log_printf(CS_LOG_SETUP, "%s Bulk.density: %.1e\n",
                  id, soil->bulk_density);
    cs_log_printf(CS_LOG_SETUP, "%s Max.Porosity: %.3e (=saturated_moisture)\n",
                  id, soil->porosity);
    cs_log_printf(CS_LOG_SETUP, "%s Absolute permeability\n", id);
    cs_log_printf(CS_LOG_SETUP, "%s [%-4.2e %4.2e %4.2e;\n", id,
                  soil->abs_permeability[0][0],
                  soil->abs_permeability[0][1],
                  soil->abs_permeability[0][2]);
    cs_log_printf(CS_LOG_SETUP, "%s  %-4.2e %4.2e %4.2e;\n", id,
                  soil->abs_permeability[1][0],
                  soil->abs_permeability[1][1],
                  soil->abs_permeability[1][2]);
    cs_log_printf(CS_LOG_SETUP, "%s  %-4.2e %4.2e %4.2e]\n", id,
                  soil->abs_permeability[2][0],
                  soil->abs_permeability[2][1],
                  soil->abs_permeability[2][2]);

    /* Display the model parameters */

    switch (soil->model) {

    case CS_GWF_SOIL_SATURATED:
        cs_log_printf(CS_LOG_SETUP, "%s Model: *Saturated*\n", id);
      break;

    case CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_SINGLE_PHASE:
      {
        const cs_gwf_soil_vgm_spf_param_t  *sp = soil->model_param;

        cs_log_printf(CS_LOG_SETUP, "%s Model: "
                      "*Single_phase_Van_Genuchten_Mualem*\n", id);
        cs_log_printf(CS_LOG_SETUP, "%s Parameters:", id);
        cs_log_printf(CS_LOG_SETUP,
                      " residual_moisture %5.3e\n", sp->residual_moisture);
        cs_log_printf(CS_LOG_SETUP, "%s Parameters:", id);
        cs_log_printf(CS_LOG_SETUP, " n= %f, scale= %f, tortuosity= %f\n",
                      sp->n, sp->scale, sp->tortuosity);
      }
      break;

    case CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_TWO_PHASE:
      {
        const cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

        cs_log_printf(CS_LOG_SETUP, "%s Model: "
                      "*Two_phase_Van_Genuchten_Mualem*\n", id);
        cs_log_printf(CS_LOG_SETUP, "%s Parameters:", id);
        cs_log_printf(CS_LOG_SETUP, " residual_saturation  %5.3e\n", sp->sl_r);
        cs_log_printf(CS_LOG_SETUP, "%s Parameters:", id);
        cs_log_printf(CS_LOG_SETUP, " saturated_saturation %5.3e\n", sp->sl_s);
        cs_log_printf(CS_LOG_SETUP, "%s Parameters:", id);
        cs_log_printf(CS_LOG_SETUP, " n %f; m= %f; pr_r= %f\n",
                      sp->n, sp->m, sp->pr_r);

        switch(sp->sle_jtype) {
        case CS_GWF_SOIL_JOIN_NOTHING:
          cs_log_printf(CS_LOG_SETUP, "%s No joining function for Sl\n", id);
          break;

        case CS_GWF_SOIL_JOIN_C1_HYPERBOLIC:
          cs_log_printf(CS_LOG_SETUP,
                        "%s C1 hyperbolic joining function for Sl\n", id);
          cs_log_printf(CS_LOG_SETUP, "%s Joining parameters:", id);
          cs_log_printf(CS_LOG_SETUP, " sle %8.6e pc %5.3e\n",
                        sp->sle_thres, sp->pc_star);
          break;

        case CS_GWF_SOIL_JOIN_C1_EXPONENTIAL:
          cs_log_printf(CS_LOG_SETUP,
                        "%s C1 exponential joining function for Sl\n", id);
          cs_log_printf(CS_LOG_SETUP, "%s Joining parameters:", id);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Invalid joining function.", __func__);
        }

        switch(sp->krg_jtype) {
        case CS_GWF_SOIL_JOIN_NOTHING:
        case CS_GWF_SOIL_JOIN_C1_EXPONENTIAL:
        case CS_GWF_SOIL_JOIN_C1_HYPERBOLIC:
          cs_log_printf(CS_LOG_SETUP, "%s No joining function for krg\n", id);
          break;

        case CS_GWF_SOIL_JOIN_C1_POLY_ORDER2:
          cs_log_printf(CS_LOG_SETUP,
                        "%s C1 2nd order poly. joining function for krg\n", id);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Invalid joining function.", __func__);
        }
      }
      break;

    case CS_GWF_SOIL_USER:
      cs_log_printf(CS_LOG_SETUP, "%s Model: *User-defined*\n", id);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid model for groundwater module.\n"
                " Please check your settings.");

    } /* Switch model */

  } /* Loop on soils */

  cs_log_printf(CS_LOG_SETUP, "\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the soil properties
 *
 * \param[in] time_eval      time at which one evaluates properties
 * \param[in] mesh           pointer to the mesh structure
 * \param[in] connect        pointer to the cdo connectivity
 * \param[in] cdoq           pointer to the cdo quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_update(cs_real_t                     time_eval,
                   const cs_mesh_t              *mesh,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *cdoq)
{
  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];
    if (soil == NULL)
      continue;

    switch (soil->model) {

    case CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_SINGLE_PHASE:
    case CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_TWO_PHASE:
    case CS_GWF_SOIL_USER:
      {
        assert(soil->update_properties != NULL);

        const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);

        soil->update_properties(time_eval,
                                mesh, connect, cdoq,
                                zone,
                                soil);
      }
      break;

    default:
      break; /* Do nothing (for instance in the case of a saturated soil which
                is constant (steady and uniform) */

    } /* Switch on the soil model */

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the soil state associated to each cell w.r.t. the given
 *        liquid saturation
 *
 * \param[in] n_cells      number of mesh cells
 * \param[in] sliq         values of the liquid saturation in each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_update_soil_state(cs_lnum_t            n_cells,
                              const cs_real_t     *sliq)
{
  assert(_soil_state_array != NULL);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    if (sliq[c_id] < FLT_MIN)
      _soil_state_array[c_id] = CS_GWF_SOIL_STATE_DRY;
    else if (sliq[c_id] > 1 - FLT_MIN)
      _soil_state_array[c_id] = CS_GWF_SOIL_STATE_SATURATED;
    else
      _soil_state_array[c_id] = CS_GWF_SOIL_STATE_UNSATURATED;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the definition of the soil porosity and absolute permeability
 *        (which are properties always defined in the GWF module). One relies
 *        on the definition of these properties in each soil.
 *
 * \param[in, out] abs_permeability    pointer to a cs_property_t structure
 * \param[in, out] soil_porosity       pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_define_shared_properties(cs_property_t   *abs_permeability,
                                     cs_property_t   *soil_porosity)
{
  assert(abs_permeability != NULL && soil_porosity != NULL);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    /* Define the absolute permeability */

    if (abs_permeability->type & CS_PROPERTY_ISO) {

      assert(soil->abs_permeability_dim == 1);
      cs_property_def_iso_by_value(abs_permeability,
                                   z->name,
                                   soil->abs_permeability[0][0]);

    }
    else if (abs_permeability->type & CS_PROPERTY_ANISO) {

      cs_property_def_aniso_by_value(abs_permeability,
                                     z->name,
                                     soil->abs_permeability);

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of property.\n", __func__);

    /* Set the soil porosity */

    cs_property_def_iso_by_value(soil_porosity,
                                 z->name,
                                 soil->porosity);

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the definition of the soil porosity and absolute porosity (which
 *        are properties always defined). This relies on the definition of
 *        each soil.
 *
 * \param[in, out] moisture_content   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_define_sspf_property(cs_property_t   *moisture_content)
{
  assert(moisture_content != NULL);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    if (soil->model != CS_GWF_SOIL_SATURATED)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid way of setting soil parameter.\n"
                " All soils are not considered as saturated.", __func__);

    /* Set the moisture content. In this case, one set the moisture content to
       the soil porosity since one considers that the soil is fully
       saturated */

    const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

    cs_property_def_iso_by_value(moisture_content,
                                 z->name,
                                 soil->porosity);

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if one needs to compute the capillarity at the cell centers
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_soil_need_cell_capillarity_pressures(void)
{
  for (int i = 0; i < _n_soils; i++) {

    const cs_gwf_soil_t  *soil = _soils[i];

    if (soil->hydraulic_model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE ||
        soil->hydraulic_model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE) {

      if (soil->model == CS_GWF_SOIL_USER)
        return true;  /* Be sure that all quantities are available even if not
                         needed */

      if (soil->model == CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_TWO_PHASE) {

        cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

        if (!sp->sl_interpolation)
          return true;

      }

    } /* Two-phase flow modelling */

    /* Other modellings do not imply a capillarity pressure */

  } /* Loop on soils */

  return false; /* Default behavior */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build an array storing the dual volume associated to each vertex
 *        taking into account the porosity of the soil
 *        The computed quantity is stored as a static array. Use the function
 *        cs_gwf_soil_get_dual_vol_l()
 *
 * \param[in] cdoq     additional geometrical quantities for CDO schemes
 * \param[in] connect  additional connectivities for CDO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_build_dual_porous_volume(const cs_cdo_quantities_t    *cdoq,
                                     const cs_cdo_connect_t       *connect)
{
  assert(cdoq != NULL && connect != NULL);
  assert(cdoq->dual_vol != NULL);

  const cs_lnum_t  n_vertices = cdoq->n_vertices;

  if (_dual_porous_volume == NULL)
    BFT_MALLOC(_dual_porous_volume, n_vertices, double);
  else
    BFT_REALLOC(_dual_porous_volume, n_vertices, double);

  cs_array_real_fill_zero(n_vertices, _dual_porous_volume);

  if (_n_soils == 1) {

    const cs_gwf_soil_t  *soil = _soils[0];

    cs_array_real_set_wscalar(n_vertices, soil->porosity, cdoq->dual_vol,
                              _dual_porous_volume);

    /* cdoq->dual_vol is already synchronized (parallel sum reduction) */

  }
  else { /* Several soils to handle */

    const cs_adjacency_t  *c2v = connect->c2v;

    for (int s_id = 0; s_id < _n_soils; s_id++) {

      const cs_gwf_soil_t  *soil = _soils[s_id];
      const cs_zone_t  *z = cs_volume_zone_by_id(soil->zone_id);

      assert(z != NULL);

#     pragma omp parallel for if (z->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t  c_id = z->elt_ids[i];
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          _dual_porous_volume[c2v->ids[j]] += soil->porosity * cdoq->pvol_vc[j];

      } /* Loop on cells */

    } /* Loop on soils */

    /* Parallel synchronization */

    if (connect->vtx_ifs != NULL)
      cs_interface_set_sum(connect->vtx_ifs,
                           n_vertices,
                           1, false, /* stride, interlace */
                           CS_REAL_TYPE,
                           _dual_porous_volume);

  } /* Several soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the array storing the dual volume weighted by the soil porosity
 *        Array of size n_vertices
 *
 * \return a pointer to the requested array
 */
/*----------------------------------------------------------------------------*/

const double *
cs_gwf_soil_get_dual_porous_volume(void)
{
  return _dual_porous_volume;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the array storing the associated soil for each cell
 *
 * \return a pointer to the array
 */
/*----------------------------------------------------------------------------*/

const short int *
cs_gwf_soil_get_cell2soil(void)
{
  return _cell2soil_ids;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the array storing the soil state associated to each cell
 *
 * \return a pointer to the array (may be NULL)
 */
/*----------------------------------------------------------------------------*/

const int *
cs_gwf_soil_get_soil_state(void)
{
  return _soil_state_array;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the porosity value for the given soil id
 *
 * \param[in] soil_id      id of the requested soil
 *
 * \return the value of the soil porosity
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_soil_get_porosity(int   soil_id)
{
  cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);

  if (soil == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty soil.\n", __func__);

  return soil->porosity;  /* = saturated moisture or max. liquid saturation */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the saturated moisture for the given soil id
 *
 * \param[in]  soil_id     id of the requested soil
 *
 * \return the value of the saturated moisture
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_soil_get_saturated_moisture(int   soil_id)
{
  /* Avoid a naming which may be disturbing when handling saturated
     single-phase flows */
  return cs_gwf_soil_get_porosity(soil_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the max dim (aniso=9; iso=1) for the absolute permeability
 *         associated to each soil
 *
 * \return the associated max. dimension
 */
/*----------------------------------------------------------------------------*/

int
cs_gwf_soil_get_permeability_max_dim(void)
{
  int dim = 0;

  if (_n_soils < 1)
    return dim;

  if (_soils == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The soil structure is not allocated whereas %d soils"
              " have been added.\n", __func__, _n_soils);

  for (int i = 0; i < _n_soils; i++) {

    cs_gwf_soil_t  *soil = _soils[i];

    dim = CS_MAX(dim, soil->abs_permeability_dim);

  }

  return dim;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a soil defined by a Van Genuchten-Mualem model in the case of
 *        single-phase flow in an (unsaturated) porous media
 *
 *        The (effective) liquid saturation (also called moisture content)
 *        follows the identity
 *        S_l,eff = (S_l - theta_r)/(theta_s - theta_r)
 *                = (1 + |alpha . h|^n)^(-m)
 *
 *        The isotropic relative permeability is defined as:
 *        k_r = S_l,eff^L * (1 - (1 - S_l,eff^(1/m))^m))^2
 *        where m = 1 -  1/n
 *
 * \param[in, out] soil       pointer to a cs_gwf_soil_t structure
 * \param[in]      theta_r    residual moisture
 * \param[in]      alpha      scale parameter (in m^-1)
 * \param[in]      n          shape parameter
 * \param[in]      L          tortuosity parameter
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_vgm_spf_param(cs_gwf_soil_t         *soil,
                              double                 theta_r,
                              double                 alpha,
                              double                 n,
                              double                 L)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  cs_gwf_soil_vgm_spf_param_t  *sp = soil->model_param;

  if (soil->model != CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_SINGLE_PHASE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil model is not Van Genuchten\n", __func__);
  if (sp == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil context not allocated\n", __func__);
  if (n <= FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid value for n = %6.4e (the shape parameter).\n"
              "This value should be > 0.\n", __func__, n);

  sp->residual_moisture = theta_r;

  /* Additional advanced settings */

  sp->n = n;
  sp->m = 1 - 1/sp->n;
  sp->scale = alpha;
  sp->tortuosity = L;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the parameters related to a Van Genuchten-Mualen model to defined
 *        the behavior of a soil in the case of two-phase flow in an porous
 *        media
 *
 *        The (effective) liquid saturation follows the identity
 *        sl_eff = (sl - sl_r)/(sl_s - sl_r)
 *                = (1 + |Pc/Pr_r|^n)^(-m)
 *        where m = 1 -  1/n
 *
 *        The isotropic relative permeability in the liquid and gaz are defined
 *        as:
 *        krl = sl_eff^(1/2) * (1 - (1 - sl_eff^(1/m))^m))^2
 *        krg = (1 - sl_eff)^(1/2) * (1 - sl_eff^(1/m))^(2m)
 *
 * \param[in, out] soil         pointer to a cs_gwf_soil_t structure
 * \param[in]      n            shape parameter
 * \param[in]      pr_r         reference (capillarity) pressure
 * \param[in]      sl_r         residual liquid saturation
 * \param[in]      sl_s         saturated (max.) liquid saturation
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_vgm_tpf_param(cs_gwf_soil_t         *soil,
                              double                 n,
                              double                 pr_r,
                              double                 sl_r,
                              double                 sl_s)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

  if (soil->model != CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_TWO_PHASE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil model is not the one expected\n", __func__);
  if (sp == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil context not allocated\n", __func__);
  if (n - 1 <= FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid value for n = %6.4e (the shape parameter).\n"
              "This value should be > 1.\n", __func__, n);

  assert(pr_r > FLT_MIN);

  /* Main parameters */

  sp->n = n;
  sp->m = 1 - 1/sp->n;
  sp->inv_m = 1 + 1/(sp->n-1);
  sp->pr_r = pr_r;
  sp->inv_pr_r = 1./pr_r;

  /* Parameters to define the effective liquid saturation */

  sp->sl_r = sl_r;
  sp->sl_s = sl_s;
  sp->sl_range = sl_s - sl_r;

  /* Additional advanced settings (default settings) */

  cs_gwf_soil_set_vgm_tpf_advanced_param(soil,
                                         true, /* Sl interpolation */
                                         CS_GWF_SOIL_JOIN_C1_EXPONENTIAL,
                                         CS_GWF_SOIL_JOIN_NOTHING,
                                         0.9999);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set advanced parameter settings related to a Van Genuchten-Mualen
 *        soil model
 *
 * \param[in, out] soil        pointer to a cs_gwf_soil_t structure
 * \param[in]      sl_interp   interpolate Sl rather than Pc
 * \param[in]      sle_jtype   type of joining function for the effective Sl
 * \param[in]      krg_jtype   type of joining function for krg
 * \param[in]      sle_thres   value of the effective liquid saturation above
 *                             which a joining function is used
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_vgm_tpf_advanced_param(cs_gwf_soil_t             *soil,
                                       bool                       sl_interp,
                                       cs_gwf_soil_join_type_t    sle_jtype,
                                       cs_gwf_soil_join_type_t    krg_jtype,
                                       double                     sle_thres)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  cs_gwf_soil_vgm_tpf_param_t  *sp = soil->model_param;

  if (soil->model != CS_GWF_SOIL_VAN_GENUCHTEN_MUALEM_TWO_PHASE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: soil model is not the one expected\n", __func__);
  if (soil->abs_permeability_dim != 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of property for the permeability.\n"
              " Please check your settings.", __func__);

  /* Set the wanted behavior according to the joining type */

  sp->sle_jtype = sle_jtype;
  sp->krg_jtype = krg_jtype;
  sp->sle_thres = sle_thres;
  sp->sl_interpolation = sl_interp;

  if (!sl_interp) {
    if (krg_jtype != CS_GWF_SOIL_JOIN_NOTHING) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf("Joining function for krg not taken into account.");
    }
  }

  switch (soil->abs_permeability_dim) {

  case 1:
    /* Isotropic permeability */
    /* ====================== */
    if (sle_thres > 1 - FLT_MIN) { /* No joining function */

      sp->sle_jtype = CS_GWF_SOIL_JOIN_NOTHING;
      if (sl_interp)
        soil->update_properties = _update_iso_soil_vgm_tpf_sl;
      else
        soil->update_properties = _update_iso_soil_vgm_tpf_pc;

    }
    else {

      switch (sle_jtype) {

      case CS_GWF_SOIL_JOIN_NOTHING:
        if (sl_interp) {
          soil->update_properties = _update_iso_soil_vgm_tpf_sl;
          if (krg_jtype != CS_GWF_SOIL_JOIN_NOTHING) {
            cs_base_warn(__FILE__, __LINE__);
            bft_printf("Joining function for krg not taken into account.");
          }
        }
        else
          soil->update_properties = _update_iso_soil_vgm_tpf_pc;
        break;

      case CS_GWF_SOIL_JOIN_C1_HYPERBOLIC:
        sp->sle_jtype = CS_GWF_SOIL_JOIN_C1_HYPERBOLIC;
        _joining_param_vgm(sp);
        if (sl_interp) {
          if (krg_jtype != CS_GWF_SOIL_JOIN_NOTHING)
            soil->update_properties =
              _update_iso_soil_vgm_tpf_c1_hyperbolic_sl_krg;
          else
            soil->update_properties = _update_iso_soil_vgm_tpf_c1_hyperbolic_sl;
        }
        else
          soil->update_properties = _update_iso_soil_vgm_tpf_c1_hyperbolic_pc;
        break;

      case CS_GWF_SOIL_JOIN_C1_EXPONENTIAL:
        sp->sle_jtype = CS_GWF_SOIL_JOIN_C1_EXPONENTIAL;
        _joining_param_vgm(sp);
        if (sl_interp) {
          if (krg_jtype != CS_GWF_SOIL_JOIN_NOTHING)
            soil->update_properties =
              _update_iso_soil_vgm_tpf_c1_exponential_sl_krg;
          else
            soil->update_properties =
              _update_iso_soil_vgm_tpf_c1_exponential_sl;
        }
        else
          soil->update_properties = _update_iso_soil_vgm_tpf_c1_exponential_pc;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of joining function.\n"
                  " Please check your settings.", __func__);

      } /* switch */

    } /* sle_thres < 1 */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of property for the permeability.\n"
              " Please check your settings.", __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a soil defined by a user-defined model
 *
 * \param[in, out] soil              pointer to a cs_gwf_soil_t structure
 * \param[in]      param             pointer to a structure cast on-the-fly
 * \param[in]      update_func       function pointer to update propoerties
 * \param[in]      free_param_func   function pointer to free the param struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_user_model_param(cs_gwf_soil_t               *soil,
                                 void                        *param,
                                 cs_gwf_soil_update_t        *update_func,
                                 cs_gwf_soil_free_param_t    *free_param_func)
{
  if (soil == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_soil));

  if (soil->model != CS_GWF_SOIL_USER)
    bft_error(__FILE__, __LINE__, 0,
              " %s: soil model is not user-defined.\n", __func__);

  /* Set pointers */

  soil->model_param = param;
  soil->update_properties = update_func;
  soil->free_model_param = free_param_func;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
