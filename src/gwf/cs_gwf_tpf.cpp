/*============================================================================
 * Main functions dedicated to the modelling of two-phase flows in a porous
 * media. This media is always considered as unsaturated. Two sub-models are
 * considered: miscible (MTPF) or immiscible (ITPF)
 *============================================================================*/

/* VERS */

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_array.h"
#include "cs_cdovb_priv.h"
#include "cs_field.h"
#include "cs_gwf_priv.h"
#include "cs_gwf_soil.h"
#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_param_types.h"
#include "cs_physical_constants.h"
#include "cs_post.h"
#include "cs_property.h"
#include "cs_reco.h"
#include "cs_time_plot.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_tpf.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_gwf_tpf.c

  \brief Main functions dedicated to the modelling of two-phase flows in a
         porous media. This media is always considered as unsaturated. Two
         sub-models are considered: miscible (MTPF) or immiscible (ITPF).
*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* Redefined names of function from cs_math to get shorter names */

#define CS_GWF_TPF_DBG 0
#define CS_GWF_TPF_N_OUTPUT_VARS  5

/*============================================================================
 * Local definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for the way to compute the value(s) of the mass
 *        density in the gas phase for the component.
 *
 * \param[in]      id        entity id
 * \param[in]      mh_ov_rt  value of the pre-computed scaling coefficient
 * \param[in]      connect   additional adjacencies for CDO schemes
 * \param[in]      cdoq      additional quantities for CDO schemes
 * \param[in]      pg        values of the gas pressure
 * \param[in, out] tpf       pointer to the model context to update
 * \param[out]     rhog_h    computed value(s)
 */
/*----------------------------------------------------------------------------*/

typedef void
(_compute_rhog_h_t)(cs_lnum_t                    id,
                    double                       mh_ov_rt,
                    const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *cdoq,
                    const cs_real_t             *pg,
                    cs_gwf_tpf_t                *tpf,
                    double                      *rhog_h);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for the way to compute the value(s) of the mass
 *        density in the gas phase for the component.
 *
 * \param[in]      id        entity id
 * \param[in]      h_mh      pre-computed scaling coefficient for rhol_h
 * \param[in]      mh_ov_rt  pre-computed scaling coefficient for rhog_h
 * \param[in]      connect   additional adjacencies for CDO schemes
 * \param[in]      cdoq      additional quantities for CDO schemes
 * \param[in]      pg        values of the gas pressure
 * \param[in, out] tpf       pointer to the model context to update
 * \param[out]     rhog_h    computed value
 * \param[out]     rhol_h    computed value
 */
/*----------------------------------------------------------------------------*/

typedef void
(_compute_rhogl_h_t)(cs_lnum_t                    id,
                     double                       h_mh,
                     double                       mh_ov_rt,
                     const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *cdoq,
                     const cs_real_t             *pg,
                     cs_gwf_tpf_t                *tpf,
                     double                      *rhog_h,
                     double                      *rhol_h);

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_time_plot_t *cs_gwf_tpf_time_plot = nullptr;

static const char _output_varnames[CS_GWF_TPF_N_OUTPUT_VARS][32] = {

  "Pg",
  "Pl",
  "Pc",
  "Sliq",
  "Cliq"

};

static _compute_rhog_h_t  *compute_rhog_h  = nullptr;
static _compute_rhogl_h_t *compute_rhogl_h = nullptr;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the given property with a newly allocated array.
 *
 * \param[in] array_size     size of the array to allocate
 * \param[in] location_flag  where the values are defined
 * \param[in] pty            pointer to a property structure
 *
 */
/*----------------------------------------------------------------------------*/

static inline cs_xdef_t *
_add_pty_array(cs_lnum_t         array_size,
               cs_flag_t         location_flag,
               cs_property_t    *pty)
{
  cs_real_t *array = nullptr;
  BFT_MALLOC(array, array_size, cs_real_t);
  cs_array_real_fill_zero(array_size, array);

  cs_xdef_t *def = cs_property_def_by_array(pty,
                                            nullptr,       /* all cells */
                                            location_flag, /* data location */
                                            array,
                                            true,  /* xdef is owner */
                                            true); /* full length */

  return def;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if one needs a diffusion term in the conservation of the gas
 *        component in the case of (Pc, Pg) coupled solver
 *
 * \param[in] tpf   pointer to the structure managing the TPF models
 *
 * \returns true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
_pcpg_coupled_solver_needs_diffusion(const cs_gwf_tpf_t     *tpf)
{
  assert(tpf != nullptr);
  if (tpf->use_diffusion_view_for_darcy)
    return true;
  else {
    if (tpf->is_miscible)
      return true;
    else
      return false;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Default settings to apply to the given cs_equation_param_t
 *
 * \param[in, out] eqp     pointer to a cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_set_default_eqp_settings(cs_equation_param_t      *eqp)
{
  if (eqp == nullptr)
    return;

  cs_equation_param_set(eqp, CS_EQKEY_BC_ENFORCEMENT, "algebraic");
  cs_equation_param_set(eqp, CS_EQKEY_HODGE_TIME_ALGO, "voronoi");
  cs_equation_param_set(eqp, CS_EQKEY_HODGE_REAC_ALGO, "voronoi");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value(s) of the mass density in the gas phase for the
 *        component. Clip value at mesh vertices.
 *
 * \param[in]      v_id      vertex id
 * \param[in]      mh_ov_rt  value of the pre-computed scaling coefficient
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]      pg        values of the gas pressure
 * \param[in, out] tpf       pointer to the model context to update
 * \param[out]     rhog_h    computed value
 */
/*----------------------------------------------------------------------------*/

static inline void
_rhog_h_vtx(cs_lnum_t                    v_id,
            double                       mh_ov_rt,
            const cs_cdo_connect_t      *connect,
            const cs_cdo_quantities_t   *cdoq,
            const cs_real_t             *pg,
            cs_gwf_tpf_t                *tpf,
            double                      *rhog_h)
{
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(cdoq);
  CS_NO_WARN_IF_UNUSED(tpf);

  const double  pg_vtx = pg[v_id];
  *rhog_h = mh_ov_rt * pg_vtx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for the way to compute the value of the mass density
 *        in the gas phase and in the liquid phase for the component. Clip
 *        value at mesh vertices.
 *
 * \param[in]      v_id      vertex id
 * \param[in]      h_mh      pre-computed scaling coefficient for rhol_h
 * \param[in]      mh_ov_rt  pre-computed scaling coefficient for rhog_h
 * \param[in]      connect   additional adjacencies for CDO schemes
 * \param[in]      cdoq      additional quantities for CDO schemes
 * \param[in]      pg        values of the gas pressure
 * \param[in, out] tpf       pointer to the model context to update
 * \param[out]     rhog_h    computed value
 * \param[out]     rhol_h    computed value
 */
/*----------------------------------------------------------------------------*/

static inline void
_rhogl_h_vtx(cs_lnum_t                    v_id,
             double                       h_mh,
             double                       mh_ov_rt,
             const cs_cdo_connect_t      *connect,
             const cs_cdo_quantities_t   *cdoq,
             const cs_real_t             *pg,
             cs_gwf_tpf_t                *tpf,
             double                      *rhog_h,
             double                      *rhol_h)
{
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(cdoq);
  CS_NO_WARN_IF_UNUSED(tpf);

  const double  p = pg[v_id];

  *rhog_h = mh_ov_rt * p;
  *rhol_h = h_mh * p;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value(s) of the mass density in the gas phase for the
 *        component. Average over all cell vertices. Clip negative gas
 *        pressures.
 *
 * \param[in]      c_id      cell id
 * \param[in]      mh_ov_rt  value of the pre-computed scaling coefficient
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]      pg        values of the gas pressure
 * \param[in, out] tpf       pointer to the model context to update
 * \param[out]     rhog_h    computed value
 */
/*----------------------------------------------------------------------------*/

static void
_rhog_h_cell_mean(cs_lnum_t                    c_id,
                  double                       mh_ov_rt,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *cdoq,
                  const cs_real_t             *pg,
                  cs_gwf_tpf_t                *tpf,
                  double                      *rhog_h)
{
  CS_NO_WARN_IF_UNUSED(tpf);

  const cs_adjacency_t  *c2v = connect->c2v;

  /* Compute the mean value in a cell */

  double  pg_sum = 0;
  for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
    pg_sum += cdoq->pvol_vc[j] * pg[c2v->ids[j]];

  const double  pg_cell = pg_sum/cdoq->cell_vol[c_id];

  *rhog_h = mh_ov_rt*pg_cell;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for the way to compute the value of the mass density
 *        in the gas phase and liquid phase for the component. Clip value at
 *        mesh vertices.
 *
 * \param[in]      c_id      cell id
 * \param[in]      h_mh      pre-computed scaling coefficient for rhol_h
 * \param[in]      mh_ov_rt  pre-computed scaling coefficient for rhog_h
 * \param[in]      connect   additional adjacencies for CDO schemes
 * \param[in]      cdoq      additional quantities for CDO schemes
 * \param[in]      pg        values of the gas pressure
 * \param[in, out] tpf       pointer to the model context to update
 * \param[out]     rhog_h    computed value
 * \param[out]     rhol_h    computed value
 */
/*----------------------------------------------------------------------------*/

static void
_rhogl_h_cell_mean(cs_lnum_t                    c_id,
                   double                       h_mh,
                   double                       mh_ov_rt,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *cdoq,
                   const cs_real_t             *pg,
                   cs_gwf_tpf_t                *tpf,
                   double                      *rhog_h,
                   double                      *rhol_h)
{
  CS_NO_WARN_IF_UNUSED(tpf);

  const cs_adjacency_t  *c2v = connect->c2v;

  /* Compute the mean value in a cell */

  double  pg_sum = 0;
  for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
    pg_sum += cdoq->pvol_vc[j] * pg[c2v->ids[j]];

  const double  pg_cell = pg_sum/cdoq->cell_vol[c_id];

  *rhog_h = mh_ov_rt * pg_cell;
  *rhol_h = h_mh * pg_cell;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value(s) of the mass density in the gas phase for the
 *        component. Enable a portion of upwinding.
 *
 * \param[in]      c_id       cell id
 * \param[in]      mh_ov_rt   value of the pre-computed scaling coefficient
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      pg         values of the gas pressure
 * \param[in, out] tpf        pointer to the model context to update
 * \param[out]     rhog_h     computed value
 */
/*----------------------------------------------------------------------------*/

static void
_rhog_h_cell_upw(cs_lnum_t                    c_id,
                 double                       mh_ov_rt,
                 const cs_cdo_connect_t      *connect,
                 const cs_cdo_quantities_t   *cdoq,
                 const cs_real_t             *pg,
                 cs_gwf_tpf_t                *tpf,
                 double                      *rhog_h)
{
  const cs_adjacency_t  *c2e = connect->c2e;
  const cs_adjacency_t  *e2v = connect->e2v;

  cs_real_3_t  grdpg;
  cs_reco_grad_cell_from_pv(c_id, connect, cdoq, pg, grdpg);

  /* Compute the upwind mean value in a cell */

  double  pg_cen = 0., pg_upw = 0;
  for (cs_lnum_t j = c2e->idx[c_id] + 1; j < c2e->idx[c_id+1]; j++) {

    const cs_lnum_t  e_id = c2e->ids[j];
    const cs_lnum_t  *v_id = e2v->ids + 2*e_id;
    const cs_lnum_t  v0 = (e2v->sgn[2*e_id] < 0) ? v_id[0] : v_id[1];
    const cs_lnum_t  v1 = (v_id[0] == v0) ? v_id[1] : v_id[0];
    const cs_real_t  pg_v0 = pg[v0], pg_v1 = pg[v1];

    pg_cen += cdoq->pvol_ec[j] * (pg_v0 + pg_v1);

    /* Darcy flux is -grad() */

    if (cs_math_3_dot_product(grdpg, cdoq->edge_vector + 3*e_id) < 0)
      pg_upw += cdoq->pvol_ec[j] * pg_v1;
    else
      pg_upw += cdoq->pvol_ec[j] * pg_v0;

  }

  const double  vol_c = cdoq->cell_vol[c_id];
  const double  w = tpf->upwind_weight;

  *rhog_h = mh_ov_rt / vol_c * (w*pg_upw + (1-w)*0.5*pg_cen);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value(s) of the mass density in the gas phase for the
 *        component. Average over all cell vertices. Enable of portion of
 *        upwinding.
 *
 * \param[in]      c_id      cell id
 * \param[in]      mh_ov_rt  value of the pre-computed scaling coefficient
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]      pg        values of the gas pressure
 * \param[in, out] tpf       pointer to the model context to update
 * \param[out]     rhog_h    computed value
 */
/*----------------------------------------------------------------------------*/

static void
_rhogl_h_cell_upw(cs_lnum_t                    c_id,
                  double                       h_mh,
                  double                       mh_ov_rt,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *cdoq,
                  const cs_real_t             *pg,
                  cs_gwf_tpf_t                *tpf,
                  double                      *rhog_h,
                  double                      *rhol_h)
{
  const cs_adjacency_t  *c2e = connect->c2e;
  const cs_adjacency_t  *e2v = connect->e2v;
  const cs_real_t  *pl = tpf->l_pressure->val;

  cs_real_3_t  grdpg, grdpl;
  cs_reco_2grad_cell_from_pv(c_id, connect, cdoq, pg, pl, grdpg, grdpl);

  /* Compute the upwind mean value in a cell */

  double  pg_cen = 0., pg_g_upw = 0;
  for (cs_lnum_t j = c2e->idx[c_id] + 1; j < c2e->idx[c_id+1]; j++) {

    const cs_lnum_t  e_id = c2e->ids[j];
    const cs_lnum_t  *v_id = e2v->ids + 2*e_id;
    const cs_lnum_t  v0 = (e2v->sgn[2*e_id] < 0) ? v_id[0] : v_id[1];
    const cs_lnum_t  v1 = (v_id[0] == v0) ? v_id[1] : v_id[0];
    const cs_real_t  pg_v0 = pg[v0], pg_v1 = pg[v1];

    pg_cen += cdoq->pvol_ec[j] * (pg_v0 + pg_v1);

    /* Darcy flux is -grad() */

    if (cs_math_3_dot_product(grdpg, cdoq->edge_vector + 3*e_id) < 0)
      pg_g_upw += cdoq->pvol_ec[j] * pg_v1;
    else
      pg_g_upw += cdoq->pvol_ec[j] * pg_v0;

  }

  const double  upww = tpf->upwind_weight;
  const double  inv_vol_c = 1./cdoq->cell_vol[c_id];
  const double  p = (upww*pg_g_upw + (1-upww)*0.5*pg_cen);

  *rhog_h = inv_vol_c * mh_ov_rt * p;
  *rhol_h = inv_vol_c * h_mh     * p;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a source term relying on the stiffness matrix of the conservation
          equation for the hydrogen when a segregated solver is used. Case of a
          vertex-based scheme.
 *
 *        Generic function prototype for a hook during the cellwise building
 *        of the linear system.
 *        Fit the cs_equation_build_hook_t prototype. This function may be
 *        called by different OpenMP threads
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context to cast for this discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] context     pointer to a context structure
 * \param[in, out] mass_hodge  pointer to a cs_hodge_t structure (mass matrix)
 * \param[in, out] diff_hodge  pointer to a cs_hodge_t structure (diffusion)
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_build_h_eq_diff_st(const cs_equation_param_t     *eqp,
                    const cs_equation_builder_t   *eqb,
                    const void                    *eq_context,
                    const cs_cell_mesh_t          *cm,
                    void                          *context,
                    cs_hodge_t                    *mass_hodge,
                    cs_hodge_t                    *diff_hodge,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  CS_NO_WARN_IF_UNUSED(eqp);
  CS_NO_WARN_IF_UNUSED(eqb);
  CS_NO_WARN_IF_UNUSED(mass_hodge);

  const cs_cdovb_scaleq_t  *eqc = (const cs_cdovb_scaleq_t *)eq_context;
  assert(csys->n_dofs == cm->n_vc);

  cs_gwf_tpf_t *tpf = (cs_gwf_tpf_t *)context;

  /* Modify the property data associated to the hodge operator related to the
     diffusion term */

  cs_property_data_t  *saved = diff_hodge->pty_data;
  cs_property_data_t  tmp = cs_property_data_define(saved->need_tensor,
                                                    saved->need_eigen,
                                                    tpf->diff_hl_pty);

  diff_hodge->pty_data = &tmp;

  /* Diffusion part of the source term to add */

  cs_hodge_evaluate_property_cw(cm, cb->t_pty_eval, cb->cell_flag, diff_hodge);

  /* Define the local stiffness matrix: local matrix owned by the cellwise
     builder (store in cb->loc) */

  eqc->get_stiffness_matrix(cm, diff_hodge, cb);

  /* Retrieve the values pl^{k,n+1} */

  double  *vec = cb->values, *matvec = cb->values + cm->n_vc;
  for (int v = 0; v < cm->n_vc; v++)
    vec[v] = tpf->l_pressure->val[cm->v_ids[v]];

  cs_sdm_square_matvec(cb->loc, vec, matvec);

  /* Update the rhs */

  for (int v = 0; v < cm->n_vc; v++)
    csys->rhs[v] -= matvec[v];

  /* Set the diffusion property data back to the initial pointer */

  diff_hodge->pty_data = saved;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the advection field/arrays related to the Darcy flux in the
 *        liquid phase when a two-phase flow model is used.
 *        This relies on the case of CDO-Vb schemes and a Darcy flux defined
 *        at dual faces
 *
 *        cell_values could be set to nullptr when the space discretization does
 *        not request these values for the update.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      dof_values    values to consider for the update
 * \param[in]      cell_values   values to consider for the update or nullptr
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      cur2prev      true or false
 * \param[in, out] darcy         pointer to the darcy flux structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_darcy_l(const cs_cdo_connect_t      *connect,
                const cs_cdo_quantities_t   *cdoq,
                const cs_real_t             *dof_values,
                const cs_real_t             *cell_values,
                cs_real_t                    t_eval,
                bool                         cur2prev,
                cs_gwf_darcy_flux_t         *darcy)
{
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(cdoq);

  assert(darcy != nullptr);
  assert(darcy->flux_val != nullptr);
  assert(darcy->update_input != nullptr);

  cs_adv_field_t  *adv = darcy->adv_field;
  assert(adv != nullptr);
  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t  *tpf = (cs_gwf_tpf_t *)darcy->update_input;

  /* Update the array of flux values associated to the advection field.
   *
   * diff_pty for the w_eq --> rho_l * abs_perm * krl / mu_l Thus, one needs to
   * divide by rho_l for the Darcy flux. Thanks to the rescaling, there is no
   * need to define a new property which is very close to the diff_wl_pty
   */

  if (tpf->solver_type != CS_GWF_TPF_SOLVER_PCPG_COUPLED) {

    cs_property_set_scaling_factor(tpf->diff_wl_pty, 1./tpf->l_mass_density);

    cs_equation_compute_diffusive_flux(tpf->w_eq,
                                       nullptr, /* eqp --> default */
                                       nullptr, /* diff_pty --> default */
                                       dof_values,
                                       cell_values,
                                       darcy->flux_location,
                                       t_eval,
                                       darcy->flux_val);

    cs_property_unscale(tpf->diff_wl_pty);

  }
  else { /* CS_GWF_TPF_SOLVER_PCPG_COUPLED */

    cs_property_set_scaling_factor(tpf->diff_wg_pty, 1./tpf->l_mass_density);

    cs_equation_compute_diffusive_flux(tpf->w_eq,
                                       nullptr, /* eqp --> default */
                                       tpf->diff_wg_pty,
                                       dof_values,
                                       cell_values,
                                       darcy->flux_location,
                                       t_eval,
                                       darcy->flux_val);

    cs_property_unscale(tpf->diff_wg_pty);

  }

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != nullptr);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);

  /* Update the Darcy flux at the boundary (take into account the BCs) since
     the liquid pressure is the unknown */

  cs_gwf_darcy_flux_update_on_boundary(t_eval, tpf->w_eq, adv);

  cs_field_t  *bdy_nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  if (bdy_nflx
      != nullptr) { /* Values of the Darcy flux at boundary face exist */

    if (cur2prev)
      cs_field_current_to_previous(bdy_nflx);

    /* Set the new values of the field related to the normal boundary flux */

    cs_advection_field_across_boundary(adv, t_eval, bdy_nflx->val);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the advection field/arrays related to the Darcy flux in the
 *        gas phase when a two-phase flow model is used together with a coupled
 *        solver.
 *        This relies on the case of CDO-Vb schemes and a Darcy flux defined
 *        at dual faces
 *
 *        cell_values could be set to nullptr when the space discretization does
 *        not request these values for the update.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      dof_values    values to consider for the update
 * \param[in]      cell_values   values to consider for the update or nullptr
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      cur2prev      true or false
 * \param[in, out] darcy         pointer to the darcy flux structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_darcy_g_coupled(const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *cdoq,
                        const cs_real_t             *dof_values,
                        const cs_real_t             *cell_values,
                        cs_real_t                    t_eval,
                        bool                         cur2prev,
                        cs_gwf_darcy_flux_t         *darcy)
{
  CS_NO_WARN_IF_UNUSED(cur2prev);

  assert(darcy != nullptr);
  assert(darcy->flux_val != nullptr);
  assert(darcy->update_input != nullptr);
  assert(darcy->adv_field != nullptr);

  cs_adv_field_t  *adv = darcy->adv_field;

  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t   *tpf = (cs_gwf_tpf_t *)darcy->update_input;
  cs_equation_t  *eq = tpf->h_eq;

  /* Update the array of flux values associated to the advection field */

  cs_equation_compute_diffusive_flux(eq,
                                     nullptr, /* eqp --> default */
                                     tpf->diff_g_pty,
                                     dof_values,
                                     cell_values,
                                     darcy->flux_location,
                                     t_eval,
                                     darcy->flux_val);

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != nullptr);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);

  /* Update the Darcy flux at boundary faces if needed */

  cs_gwf_darcy_flux_update_on_boundary_wo_eq(connect, cdoq, vel->val, adv);

  cs_field_t  *bdy_nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  if (bdy_nflx
      != nullptr) { /* Values of the Darcy flux at boundary face exist */

    if (cur2prev)
      cs_field_current_to_previous(bdy_nflx);

    /* Set the new values of the field related to the normal boundary flux */

    cs_advection_field_across_boundary(adv, t_eval, bdy_nflx->val);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the advection field/arrays related to the Darcy flux in the
 *        gas phase when a two-phase flow model is used together with a
 *        segregated solver.
 *        This relies on the case of CDO-Vb schemes and a Darcy flux defined
 *        at dual faces
 *
 *        cell_values could be set to nullptr when the space discretization does
 *        not request these values for the update.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      dof_values    values to consider for the update
 * \param[in]      cell_values   values to consider for the update or nullptr
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      cur2prev      true or false
 * \param[in, out] darcy         pointer to the darcy flux structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_darcy_g_segregated(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *cdoq,
                           const cs_real_t             *dof_values,
                           const cs_real_t             *cell_values,
                           cs_real_t                    t_eval,
                           bool                         cur2prev,
                           cs_gwf_darcy_flux_t         *darcy)
{
  CS_NO_WARN_IF_UNUSED(connect);
  CS_NO_WARN_IF_UNUSED(cdoq);
  CS_NO_WARN_IF_UNUSED(cur2prev);

  assert(darcy != nullptr);
  assert(darcy->flux_val != nullptr);
  assert(darcy->update_input != nullptr);
  assert(darcy->adv_field != nullptr);

  cs_adv_field_t  *adv = darcy->adv_field;

  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t   *tpf = (cs_gwf_tpf_t *)darcy->update_input;
  cs_equation_t  *eq = tpf->h_eq;

  /* Update the array of flux values associated to the advection field */

  cs_equation_compute_diffusive_flux(eq,
                                     nullptr, /* eqp --> default */
                                     nullptr, /* diff_pty --> default */
                                     dof_values,
                                     cell_values,
                                     darcy->flux_location,
                                     t_eval,
                                     darcy->flux_val);

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != nullptr);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);

  /* Update the Darcy flux at boundary faces if needed */

  cs_field_t  *bdy_nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  if (bdy_nflx
      != nullptr) { /* Values of the Darcy flux at boundary face exist */

    if (cur2prev)
      cs_field_current_to_previous(bdy_nflx);

    /* Set the new values of the field related to the normal boundary flux */

    cs_advection_field_across_boundary(adv, t_eval, bdy_nflx->val);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the advection field/arrays related to the Darcy flux in the
 *        gas phase when a two-phase flow model is used.
 *        This relies on the case of CDO-Vb schemes and a Darcy flux defined
 *        at dual faces
 *
 *        cell_values could be set to nullptr when the space discretization does
 *        not request these values for the update.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      dof_values    values to consider for the update
 * \param[in]      cell_values   values to consider for the update or nullptr
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in]      cur2prev      true or false
 * \param[in, out] darcy         pointer to the darcy flux structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_darcy_t(const cs_cdo_connect_t      *connect,
                const cs_cdo_quantities_t   *cdoq,
                const cs_real_t             *dof_values,
                const cs_real_t             *cell_values,
                cs_real_t                    t_eval,
                bool                         cur2prev,
                cs_gwf_darcy_flux_t         *darcy)
{
  CS_NO_WARN_IF_UNUSED(dof_values);
  CS_NO_WARN_IF_UNUSED(cell_values);

  assert(darcy != nullptr);
  assert(darcy->flux_val != nullptr);
  assert(darcy->update_input != nullptr);
  assert(darcy->adv_field != nullptr);

  cs_adv_field_t  *adv = darcy->adv_field;

  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t               *tpf     = (cs_gwf_tpf_t *)darcy->update_input;
  const cs_gwf_darcy_flux_t  *l_darcy = tpf->l_darcy;
  const cs_gwf_darcy_flux_t  *g_darcy = tpf->g_darcy;

  const double  hmh = tpf->h_molar_mass * tpf->henry_constant;
  const double  mh_ov_rt =
    tpf->h_molar_mass / (tpf->ref_temperature * cs_physical_constants_r);
  const cs_adjacency_t  *c2e = connect->c2e;

# pragma omp parallel for if (6*cdoq->n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < c2e->idx[cdoq->n_cells]; i++)
    darcy->flux_val[i] =
      hmh * l_darcy->flux_val[i] + mh_ov_rt * g_darcy->flux_val[i];

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != nullptr);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the Darcy flux structures
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      t_eval       time at which one performs the evaluation
 * \param[in]      cur2prev     true or false
 * \param[in, out] tpf          pointer to a TPF model context
 */
/*----------------------------------------------------------------------------*/

static void
_update_darcy_fluxes(const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *cdoq,
                     cs_real_t                    t_eval,
                     bool                         cur2prev,
                     cs_gwf_tpf_t                *tpf)
{
  cs_real_t *dof_vals = nullptr, *cell_vals = nullptr;

  /* Darcy velocity/flux in the liquid phase. The liquid pressure is the
     unknown for the coupled and the segregated solver */

  cs_gwf_darcy_flux_t  *l_darcy = tpf->l_darcy;

  cs_gwf_get_value_pointers(tpf->w_eq, &dof_vals, &cell_vals);

  l_darcy->update_func(connect,
                       cdoq,
                       dof_vals,
                       cell_vals,
                       t_eval,
                       cur2prev,
                       l_darcy);

  /* Darcy velocity/flux in the gas phase. */

  cs_gwf_darcy_flux_t  *g_darcy = tpf->g_darcy;

  g_darcy->update_func(
    connect, cdoq, tpf->g_pressure->val, nullptr, t_eval, cur2prev, g_darcy);

  if (!tpf->use_diffusion_view_for_darcy) {

    /* Only useful for the computation of the advective term in the
       conservation equation for the hydrogen */

    cs_gwf_darcy_flux_t  *t_darcy = tpf->t_darcy;

    t_darcy->update_func(connect,
                         cdoq,
                         nullptr, /* values are defined inside the call */
                         nullptr, /* values are defined inside the call */
                         t_eval,
                         cur2prev,
                         t_darcy);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the pressure fields and arrays
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      update_flag  metadata associated to type of operation to do
 * \param[in]      cur2prev     true or false
 * \param[in, out] tpf          pointer to a TPF model context
 */
/*----------------------------------------------------------------------------*/

static void
_update_pressures(const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *cdoq,
                  cs_flag_t                    update_flag,
                  bool                         cur2prev,
                  cs_gwf_tpf_t                *tpf)
{
  CS_NO_WARN_IF_UNUSED(connect);

  const cs_lnum_t  n_vertices = cdoq->n_vertices;

  /* Knowing that Pc = Pg - Pl, compute the remaining pressure.  This depends
     on the choice of the solver. Pressure variables which are solved have
     already been updated */

  switch(tpf->solver_type) {

  case CS_GWF_TPF_SOLVER_PCPG_COUPLED:
    {
      /* Main pressure variables: gas and capillarity pressures */

      const cs_real_t  *g_pr = tpf->g_pressure->val;
      const cs_real_t  *c_pr = tpf->c_pressure->val;
      cs_real_t  *l_pr = tpf->l_pressure->val;

      if (cur2prev && !cs_flag_test(update_flag, CS_FLAG_INITIALIZATION))
        cs_field_current_to_previous(tpf->l_pressure);

      /* Compute the new values of the liquid pressure at vertices */

#     pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_vertices; i++)
        l_pr[i] = g_pr[i] - c_pr[i];

      /* Avoid to add an unsteady contribution at the first iteration  */

      if (update_flag & CS_FLAG_INITIALIZATION) {
        cs_array_real_copy(n_vertices, c_pr, tpf->c_pressure->val_pre);
        cs_array_real_copy(n_vertices, g_pr, tpf->g_pressure->val_pre);
        cs_field_current_to_previous(tpf->l_pressure);
      }
    }
    break;

  case CS_GWF_TPF_SOLVER_PLPC_COUPLED:
    {
      /* Main pressure variables: liquid and capillarity pressures */

      const cs_real_t  *l_pr = tpf->l_pressure->val;
      const cs_real_t  *c_pr = tpf->c_pressure->val;
      cs_real_t  *g_pr = tpf->g_pressure->val;

      if (cur2prev && !cs_flag_test(update_flag, CS_FLAG_INITIALIZATION))
        cs_field_current_to_previous(tpf->g_pressure);

      /* Compute the new values of the gas pressure at vertices */

      const double  pg_star = 2e-6; /* 2 times pg_sharp. Rescaling for
                                       positiveness */

#     pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_vertices; i++) {

        const double  pg_tilde = l_pr[i] + c_pr[i];

        if (pg_tilde < pg_star)
          g_pr[i] = pg_star / (1 + (pg_star - pg_tilde)/pg_star);
        else
          g_pr[i] = pg_tilde;

      }

      /* Avoid to add an unsteady contribution at the first iteration  */

      if (update_flag & CS_FLAG_INITIALIZATION) {
        cs_array_real_copy(n_vertices, c_pr, tpf->c_pressure->val_pre);
        cs_array_real_copy(n_vertices, l_pr, tpf->l_pressure->val_pre);
        cs_field_current_to_previous(tpf->g_pressure);
      }
    }
    break;

  case CS_GWF_TPF_SOLVER_PLPG_SEGREGATED:
    {
      /* Main pressure variables: liquid and gas pressures */

      const cs_real_t  *l_pr = tpf->l_pressure->val;
      const cs_real_t  *g_pr = tpf->g_pressure->val;
      cs_real_t  *c_pr = tpf->c_pressure->val;

      if (cur2prev && !cs_flag_test(update_flag, CS_FLAG_INITIALIZATION))
        cs_field_current_to_previous(tpf->c_pressure);

      /* Compute the values of the capillarity pressure at vertices */

#     pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_vertices; i++)
        c_pr[i] = g_pr[i] - l_pr[i];

      /* Avoid to add an unsteady contribution at the first iteration  */

      if (update_flag & CS_FLAG_INITIALIZATION) {
        cs_array_real_copy(n_vertices, c_pr, tpf->c_pressure->val_pre);
        cs_array_real_copy(n_vertices, g_pr, tpf->g_pressure->val_pre);
        cs_field_current_to_previous(tpf->c_pressure);
      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver type", __func__);
    break;

  } /* Switch on the type of solver */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the liquid saturation at cells. Must be done after the update
 *        of the liquid saturation defined following the c2v adjacency
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] tpf        pointer to the model context to update
 */
/*----------------------------------------------------------------------------*/

static void
_update_liquid_saturation_at_cells(const cs_cdo_connect_t      *connect,
                                   const cs_cdo_quantities_t   *cdoq,
                                   cs_gwf_tpf_t                *tpf)
{
  /* Reconstruction of the average liquid saturation in each cell */

  switch (tpf->approx_type) {

  case CS_GWF_TPF_APPROX_VERTEX_SUBCELL:
    cs_reco_scalar_vbyc2c_full(connect->c2v, cdoq,
                               cs_property_get_array(tpf->lsat_pty),
                               tpf->l_saturation->val);
    break;

  default:
    cs_array_real_copy(cdoq->n_cells,
                       cs_property_get_array(tpf->lsat_pty),
                       tpf->l_saturation->val);
    break;

  } /* Switch on the approximation type */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update several arrays associated to the definition of terms involved
 *        in the resolution of a miscible two-phase flow model with an
 *        isotropic absolute permeability.
 *
 *        Numerical options: (Pl, Pc) coupled solver and diffusive viewpoint of
 *        the Darcy terms in the conservation equation of the hydrogen
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] tpf        pointer to the model context to update
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_mtpf_plpc_coupled(const cs_cdo_connect_t      *connect,
                              const cs_cdo_quantities_t   *cdoq,
                              cs_gwf_tpf_t                *tpf)
{
  const cs_adjacency_t  *c2v = connect->c2v;

  const double  hmh = tpf->h_molar_mass * tpf->henry_constant;
  const double  hmh_dhl = hmh * tpf->l_diffusivity_h;
  const double  mh_ov_rt =
    tpf->h_molar_mass / (tpf->ref_temperature * cs_physical_constants_r);

  /* Values of the current pressure fields used for the update stage */

  const cs_real_t  *pg = tpf->g_pressure->val;

  /* Retrieve the arrays storing the property values */

  const cs_real_t  *lsat = cs_property_get_array(tpf->lsat_pty);
  const cs_real_t  *lcap = cs_property_get_array(tpf->lcap_pty);
  const cs_real_t  *krl = cs_property_get_array(tpf->krl_pty);
  const cs_real_t  *krg = cs_property_get_array(tpf->krg_pty);

  /* Retrieve arrays to update */

  cs_real_t  *diff_g_array = cs_property_get_array(tpf->diff_g_pty);
  cs_real_t  *diff_wl_array = cs_property_get_array(tpf->diff_wl_pty);
  cs_real_t  *diff_hc_array = cs_property_get_array(tpf->diff_hc_pty);
  cs_real_t  *diff_hl_array = cs_property_get_array(tpf->diff_hl_pty);
  cs_real_t  *time_wc_array = cs_property_get_array(tpf->time_wc_pty);
  cs_real_t  *time_hc_array = cs_property_get_array(tpf->time_hc_pty);
  cs_real_t  *time_hl_array = cs_property_get_array(tpf->time_hl_pty);

  /* Main loop on soils */

  for (int soil_id = 0; soil_id < cs_gwf_get_n_soils(); soil_id++) {

    cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);

    assert(soil != nullptr);
    assert(soil->hydraulic_model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE);
    assert(soil->abs_permeability_dim == 1);

    const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);
    assert(zone != nullptr);

    /* Soil properties and its derived quantities */

    const double  k_abs = soil->abs_permeability[0][0];
    const double  phi = soil->porosity;
    const double  rhol = tpf->l_mass_density;
    const double  phi_rhol = phi * rhol;
    const double  g_diff_coef = k_abs/tpf->g_viscosity;
    const double  l_diff_coef = k_abs/tpf->l_viscosity;

    /* Loop on cells belonging to this soil */

    switch (tpf->approx_type) {

      /* ================================= */
    case CS_GWF_TPF_APPROX_PC_CELL_AVERAGE:
    case CS_GWF_TPF_APPROX_PC_CELL_VERTEX_AVERAGE:
    case CS_GWF_TPF_APPROX_PC_EDGE_AVERAGE:
    case CS_GWF_TPF_APPROX_PC_VERTEX_AVERAGE:
      /* ================================= */

      for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

        const cs_lnum_t  c_id = zone->elt_ids[i];

        double  rhog_h, rhol_h;
        compute_rhogl_h(c_id, hmh, mh_ov_rt, connect, cdoq, pg, tpf,
                        &rhog_h, &rhol_h);

        const double  krg_coef = krg[c_id] * g_diff_coef;
        const double  krl_coef = krl[c_id] * l_diff_coef;

        /* Update terms for the Darcy flux in the gas phase */

        diff_g_array[c_id] = krg_coef;

        /* Update terms associated to the water conservation equation */

        const double  dsl_dpc = lcap[c_id];

        time_wc_array[c_id] = dsl_dpc * phi_rhol;
        diff_wl_array[c_id] = krl_coef * rhol;

        /* Update terms associated to the hydrogen conservation equation */

        const double  sl = lsat[c_id];
        const double  time_h_coef = mh_ov_rt*(1 - sl) + hmh*sl;

        time_hc_array[c_id] = (time_h_coef + (rhol_h - rhog_h) * dsl_dpc) * phi;
        time_hl_array[c_id] =  time_h_coef * phi;

        const double  diff_h_common_coef = rhog_h * krg_coef + sl * hmh_dhl;

        diff_hc_array[c_id] = diff_h_common_coef;
        diff_hl_array[c_id] = diff_h_common_coef + rhol_h * krl_coef;

      } /* Loop on cells of the zone (= soil) */
      break;

      /* ================================ */
    case CS_GWF_TPF_APPROX_VERTEX_SUBCELL:
      /* ================================ */

      for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

        const cs_lnum_t  c_id = zone->elt_ids[i];

        double  diff_g = 0, time_wc = 0, diff_wl = 0, time_hc = 0, time_hl = 0;
        double  diff_hc = 0, diff_hl = 0;

        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

          const cs_lnum_t  v_id = c2v->ids[j];

          double  rhog_h, rhol_h;
          compute_rhogl_h(v_id, hmh, mh_ov_rt, connect, cdoq, pg, tpf,
                          &rhog_h, &rhol_h);

          const double  pvc = cdoq->pvol_vc[j];
          const double  sl = lsat[j], sg = 1 - sl;
          const double  dsl_dpc = lcap[j];
          const double  krg_coef = krg[j] * g_diff_coef;
          const double  krl_coef = krl[j] * l_diff_coef;

          /* Update terms for the Darcy flux in the gaz phase */

          diff_g += pvc * krg_coef;

          /* Update terms associated to the water conservation equation */

          time_wc += pvc * dsl_dpc;
          diff_wl += pvc * krl_coef;

          /* Update terms associated to the hydrogen conservation equation */

          const double  time_h_coef = mh_ov_rt*sg + hmh*sl;

          time_hc += pvc * (time_h_coef + (rhol_h - rhog_h) * dsl_dpc);
          time_hl += pvc * time_h_coef;

          const double  diff_h_common_coef = rhog_h * krg_coef + sl * hmh_dhl;

          diff_hc += pvc * diff_h_common_coef;
          diff_hl += pvc * (diff_h_common_coef + rhol_h * krl_coef);

        } /* Loop on cell vertices */

        /* Define the cell values (what is stored in properties associated to
           each term of an equation) */

        const double  inv_volc = 1./cdoq->cell_vol[c_id];

        diff_g_array[c_id]  = inv_volc * diff_g;

        time_wc_array[c_id] = inv_volc * time_wc * phi_rhol;
        diff_wl_array[c_id] = inv_volc * diff_wl * rhol;

        time_hc_array[c_id] = inv_volc * time_hc * phi;
        time_hl_array[c_id] = inv_volc * time_hl * phi;
        diff_hc_array[c_id] = inv_volc * diff_hc;
        diff_hl_array[c_id] = inv_volc * diff_hl;

      } /* Loop on cells of the zone (= soil) */
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid way to approximate coefficients.", __func__);
      break;
    }

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update several arrays associated to the definition of terms involved
 *        in the resolution of an immiscible two-phase flow model with an
 *        isotropic absolute permeability.
 *
 *        Numerical options: (Pl, Pc) coupled solver and diffusive viewpoint of
 *        the Darcy terms in the conservation equation of the hydrogen
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] tpf        pointer to the model context to update
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_itpf_plpc_coupled(const cs_cdo_connect_t     *connect,
                              const cs_cdo_quantities_t  *cdoq,
                              cs_gwf_tpf_t               *tpf)
{
  const cs_adjacency_t  *c2v = connect->c2v;

  const double  mh_ov_rt =
    tpf->h_molar_mass / (tpf->ref_temperature * cs_physical_constants_r);

  /* Values of the current gas pressure field used for the update stage */

  const cs_real_t  *pg = tpf->g_pressure->val;

  /* Retrieve the arrays storing the property values */

  const cs_real_t  *lsat = cs_property_get_array(tpf->lsat_pty);
  const cs_real_t  *lcap = cs_property_get_array(tpf->lcap_pty);
  const cs_real_t  *krl = cs_property_get_array(tpf->krl_pty);
  const cs_real_t  *krg = cs_property_get_array(tpf->krg_pty);

  /* Retrieve arrays to update */

  cs_real_t  *diff_g_array = cs_property_get_array(tpf->diff_g_pty);
  cs_real_t  *diff_wl_array = cs_property_get_array(tpf->diff_wl_pty);
  cs_real_t  *diff_hc_array = cs_property_get_array(tpf->diff_hc_pty);
  cs_real_t  *diff_hl_array = cs_property_get_array(tpf->diff_hl_pty);
  cs_real_t  *time_wc_array = cs_property_get_array(tpf->time_wc_pty);
  cs_real_t  *time_hc_array = cs_property_get_array(tpf->time_hc_pty);
  cs_real_t  *time_hl_array = cs_property_get_array(tpf->time_hl_pty);

  /* Main loop on soils */

  for (int soil_id = 0; soil_id < cs_gwf_get_n_soils(); soil_id++) {

    cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);

    assert(soil != nullptr);
    assert(soil->hydraulic_model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE);
    assert(soil->abs_permeability_dim == 1);

    const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);
    assert(zone != nullptr);

    /* Soil properties and its derived quantities */

    const double  k_abs = soil->abs_permeability[0][0];
    const double  phi = soil->porosity;
    const double  phi_rhol = phi * tpf->l_mass_density;

    const double  g_diff_coef = k_abs/tpf->g_viscosity;
    const double  wl_diff_coef = tpf->l_mass_density * k_abs/tpf->l_viscosity;

    /* Loop on cells belonging to this soil */

    switch (tpf->approx_type) {

      /* ================================= */
    case CS_GWF_TPF_APPROX_PC_CELL_AVERAGE:
    case CS_GWF_TPF_APPROX_PC_CELL_VERTEX_AVERAGE:
    case CS_GWF_TPF_APPROX_PC_EDGE_AVERAGE:
    case CS_GWF_TPF_APPROX_PC_VERTEX_AVERAGE:
      /* ================================= */

      for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

        const cs_lnum_t  c_id = zone->elt_ids[i];

        double  rhog_h;
        compute_rhog_h(c_id, mh_ov_rt, connect, cdoq, pg, tpf, &rhog_h);

        /* Update terms for the Darcy flux in the gas phase */

        diff_g_array[c_id] = krg[c_id] * g_diff_coef;

        /* Update terms associated to the water conservation equation */

        time_wc_array[c_id] = lcap[c_id] * phi_rhol;
        diff_wl_array[c_id] = krl[c_id] * wl_diff_coef;

        /* Update terms associated to the hydrogen conservation equation */

        const double  time_h_coef = mh_ov_rt*(1 - lsat[c_id]);

        time_hc_array[c_id] = (time_h_coef - rhog_h * lcap[c_id]) * phi;
        time_hl_array[c_id] = time_h_coef * phi;

        const double  diff_h_coef = krg[c_id] * g_diff_coef * rhog_h;

        diff_hc_array[c_id] = diff_h_coef;
        diff_hl_array[c_id] = diff_h_coef;

      } /* Loop on cells belonging to the (zone) soil */
      break;

      /* ================================ */
    case CS_GWF_TPF_APPROX_VERTEX_SUBCELL:
      /* ================================ */

      for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

        const cs_lnum_t  c_id = zone->elt_ids[i];

        double  time_wc = 0, time_hc = 0, time_hl = 0;
        double  diff_g = 0, diff_wl = 0, diff_hc = 0, diff_hl = 0;

        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

          const cs_lnum_t  v_id = c2v->ids[j];
          double  rhog_h;
          compute_rhog_h(v_id, mh_ov_rt, connect, cdoq, pg, tpf, &rhog_h);

          const double  pvc = cdoq->pvol_vc[j];
          const double  sg = 1 - lsat[j];
          const double  dsldpc = lcap[j];

          /* Update terms for the Darcy flux in the gas phase */

          diff_g += pvc * krg[j] * g_diff_coef;

          /* Update terms associated to the water conservation equation */

          time_wc += pvc * dsldpc; /* (* phi_rhol) */
          diff_wl += pvc * krl[j]; /* (* wl_diff_coef) */

          /* Update terms associated to the hydrogen conservation equation */

          const double  time_h_coef = mh_ov_rt * sg;

          time_hc += pvc * (time_h_coef - rhog_h * dsldpc);
          time_hl += pvc * time_h_coef;

          const double  diff_h_coef = krg[j] * g_diff_coef * rhog_h;

          diff_hc += pvc * diff_h_coef;
          diff_hl += pvc * diff_h_coef;

        }

        const double  inv_volc = 1./cdoq->cell_vol[c_id];

        diff_g_array[c_id] = diff_g * inv_volc;

        time_wc_array[c_id] = time_wc * inv_volc * phi_rhol;
        diff_wl_array[c_id] = diff_wl * inv_volc * wl_diff_coef;

        time_hc_array[c_id] = time_hc * inv_volc * phi;
        time_hl_array[c_id] = time_hl * inv_volc * phi;
        diff_hc_array[c_id] = diff_hc * inv_volc;
        diff_hl_array[c_id] = diff_hl * inv_volc;

      } /* Loop on cells belonging to the (zone) soil */
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid way to approximate coefficients.", __func__);
      break;
    }

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update several arrays associated to the definition of terms involved
 *        in the resolution of an immiscible two-phase flow model with an
 *        isotropic absolute permeability.
 *
 *        Numerical options: (Pc, Pg) coupled solver and diffusive viewpoint of
 *        the Darcy terms in the conservation equation of the hydrogen
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] tpf        pointer to the model context to update
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_itpf_pcpg_coupled_diffview(const cs_cdo_connect_t     *connect,
                                       const cs_cdo_quantities_t  *cdoq,
                                       cs_gwf_tpf_t               *tpf)
{
  const cs_adjacency_t  *c2v = connect->c2v;

  const double  mh_ov_rt =
    tpf->h_molar_mass / (tpf->ref_temperature * cs_physical_constants_r);

  /* Values of the current gas pressure field used for the update stage */

  const cs_real_t  *pg = tpf->g_pressure->val;

  /* Retrieve the arrays storing the property values */

  const cs_real_t  *lsat = cs_property_get_array(tpf->lsat_pty);
  const cs_real_t  *lcap = cs_property_get_array(tpf->lcap_pty);
  const cs_real_t  *krl = cs_property_get_array(tpf->krl_pty);
  const cs_real_t  *krg = cs_property_get_array(tpf->krg_pty);

  /* Retrieve arrays to update */

  cs_real_t  *diff_g_array = cs_property_get_array(tpf->diff_g_pty);

  cs_real_t  *time_wc_array = cs_property_get_array(tpf->time_wc_pty);
  cs_real_t  *diff_wc_array = cs_property_get_array(tpf->diff_wc_pty);
  cs_real_t  *diff_wg_array = cs_property_get_array(tpf->diff_wg_pty);

  cs_real_t  *time_hc_array = cs_property_get_array(tpf->time_hc_pty);
  cs_real_t  *time_hg_array = cs_property_get_array(tpf->time_hg_pty);
  cs_real_t  *diff_hg_array = cs_property_get_array(tpf->diff_hg_pty);

  /* Main loop on soils */

  for (int soil_id = 0; soil_id < cs_gwf_get_n_soils(); soil_id++) {

    cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);

    assert(soil != nullptr);
    assert(soil->hydraulic_model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE);
    assert(soil->abs_permeability_dim == 1);

    const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);
    assert(zone != nullptr);

    /* Soil properties and its derived quantities */

    const double  k_abs = soil->abs_permeability[0][0];
    const double  phi = soil->porosity;
    const double  phi_rhol = phi * tpf->l_mass_density;

    const double  g_diff_coef = k_abs/tpf->g_viscosity;
    const double  rhol_diff_coef = tpf->l_mass_density * k_abs/tpf->l_viscosity;

    /* Loop on cells belonging to this soil */

    for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

      const cs_lnum_t  c_id = zone->elt_ids[i];

      double  diff_g = 0, time_wc = 0, diff_wg = 0, time_hc = 0, time_hg = 0;
      double  diff_hg = 0;

      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

        const cs_lnum_t  v_id = c2v->ids[j];
        const double  pvc = cdoq->pvol_vc[j];

        const double  rhog_h = (pg[v_id] > 0) ? mh_ov_rt * pg[v_id] : 0.;
        const double  sl = lsat[j], sg = 1 - sl;
        const double  dsl_dpc = lcap[j];

        /* Update terms for the Darcy flux in the gas phase */

        const double  diff_g_term = pvc * krg[j] * g_diff_coef;

        diff_g += diff_g_term;

        /* Update terms associated to the water conservation equation */

        time_wc += pvc * dsl_dpc;
        diff_wg += pvc * krl[j];

        /* Update terms associated to the hydrogen conservation equation */

        const double  time_h_coef = mh_ov_rt*sg;

        time_hc += pvc * (-rhog_h) * dsl_dpc;
        time_hg += pvc * time_h_coef;

        diff_hg += diff_g_term * rhog_h;

      } /* Loop on cell vertices */

      /* Define the cell values (what is stored in properties associated to
         each term of an equation) */

      const double  inv_volc = 1./cdoq->cell_vol[c_id];

      diff_g_array[c_id]  = inv_volc * diff_g;

      time_wc_array[c_id] = inv_volc * time_wc * phi_rhol;
      diff_wg_array[c_id] = inv_volc * diff_wg * rhol_diff_coef;
      diff_wc_array[c_id] = -diff_wg_array[c_id];

      time_hc_array[c_id] = inv_volc * time_hc * phi;
      time_hg_array[c_id] = inv_volc * time_hg * phi;
      diff_hg_array[c_id] = inv_volc * diff_hg;

    } /* Loop on cells of the zone (= soil) */

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the non-linear algorithm.
 *
 * \param[in]      pa       1st pressure field
 * \param[in]      pb       2nd pressure field
 * \param[in, out] tpf      pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_non_linear_algo(const cs_field_t     *pa,
                      const cs_field_t     *pb,
                      cs_gwf_tpf_t         *tpf)
{
  if (tpf->nl_algo_type == CS_PARAM_NL_ALGO_NONE)
    return;

  cs_iter_algo_t  *algo = tpf->nl_algo;
  assert(algo != nullptr);

  cs_iter_algo_reset(algo);

  /* No resolution has been done at this stage and also no previous to current
     operation.

     The current values will be the previous state after the current to
     previous operation
  */

  const cs_real_t  *pa_0 = pa->val;
  const cs_real_t  *pb_0 = pb->val;

  /* Set the normalization factor. */

  double  normalization = cs_cdo_blas_square_norm_pvsp(pa_0);
  normalization += cs_cdo_blas_square_norm_pvsp(pb_0);
  normalization = sqrt(normalization);
  if (normalization < cs_math_zero_threshold)
    normalization = 1.0;

  cs_iter_algo_set_normalization(algo, normalization);

  cs_log_printf(CS_LOG_DEFAULT, "%s: Non-linear algo. normalization=%6.4e\n",
                __func__, normalization);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new increment on the capillarity pressure
 *
 * \param[in]      mesh     pointer to a cs_mesh_t structure
 * \param[in]      w_eq     pointer on the water equation structure
 * \param[in]      h_eq     pointer on the hydrogen equation structure
 * \param[in, out] dpc_kp1  values of the increment on the capillarity pressure
 */
/*----------------------------------------------------------------------------*/

static void
_get_capillarity_pressure_increment(const cs_mesh_t       *mesh,
                                    const cs_equation_t   *w_eq,
                                    const cs_equation_t   *h_eq,
                                    cs_real_t             *dpc_kp1)
{
  const cs_equation_param_t  *w_eqp = cs_equation_get_param(w_eq);
  const cs_equation_param_t  *h_eqp = cs_equation_get_param(h_eq);
  const cs_equation_builder_t  *wl_eqb = cs_equation_get_builder(w_eq);
  const cs_equation_builder_t  *hg_eqb = cs_equation_get_builder(h_eq);
  const cs_real_t  *dpl_kp1 = wl_eqb->increment;
  const cs_real_t  *dpg_kp1 = hg_eqb->increment;

  if (w_eqp->incremental_relax_factor < 1 ||
      h_eqp->incremental_relax_factor < 1) {

    assert(w_eqp->incremental_relax_factor > 0);
    assert(h_eqp->incremental_relax_factor > 0);

    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++)
      dpc_kp1[i] = h_eqp->incremental_relax_factor*dpg_kp1[i]
                 - w_eqp->incremental_relax_factor*dpl_kp1[i];

  }
  else {

    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++)
      dpc_kp1[i] = dpg_kp1[i] - dpl_kp1[i];

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT, "%s: dpl=% 6.4e; dpg=% 6.4e\n", __func__,
                sqrt(cs_cdo_blas_square_norm_pvsp(dpl_kp1)),
                sqrt(cs_cdo_blas_square_norm_pvsp(dpg_kp1)));
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the convergence of the non-linear algorithm in case of TPF
 *        model. Work only with the capillarity pressure at vertices
 *
 * \param[in]      dpc_iter    increment values for the capillarity pressure
 * \param[in, out] algo        pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_check_cvg_nl_dpc(const cs_real_t          *dpc_iter,
                  cs_iter_algo_t           *algo)
{
  assert(algo != nullptr);

  double dpc_norm = cs_cdo_blas_square_norm_pvsp(dpc_iter);

  cs_iter_algo_update_residual(algo, sqrt(dpc_norm));

  /* Update the convergence members */

  cs_sles_convergence_state_t
    cvg_status = cs_iter_algo_update_cvg_tol_auto(algo);

  /* Monitoring */

  cs_iter_algo_log_cvg(algo, "# GWF.TPF");

  return cvg_status;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the convergence of the non-linear algorithm in case of TPF
 *        model
 *
 * \param[in]      nl_algo_type  type of non-linear algorithm
 * \param[in]      pa_pre_iter   previous iterate values for the 1st pressure
 * \param[in, out] pa_cur_iter   current  iterate values for the 1st pressure
 * \param[in]      pb_pre_iter   previous iterate values for the 2nd pressure
 * \param[in, out] pb_cur_iter   current  iterate values for the 2nd pressure
 * \param[in, out] algo          pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_check_cvg_nl(cs_param_nl_algo_t        nl_algo_type,
              const cs_real_t          *pa_pre_iter,
              cs_real_t                *pa_cur_iter,
              const cs_real_t          *pb_pre_iter,
              cs_real_t                *pb_cur_iter,
              cs_iter_algo_t           *algo)
{
  if (nl_algo_type == CS_PARAM_NL_ALGO_NONE)
    return CS_SLES_CONVERGED;

  assert(algo != nullptr);

  if (nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON) {

    /* pg_* arrays gather pg and pb (this is done during the solve step) */

    cs_iter_algo_update_anderson(algo,
                                 pa_cur_iter, /* updated during the process */
                                 pa_pre_iter,
                                 cs_cdo_blas_dotprod_2pvsp,
                                 cs_cdo_blas_square_norm_2pvsp);

  } /* Anderson acceleration */

  /* Update the residual (parallel sync. done) */


  double delta_pa = cs_cdo_blas_square_norm_pvsp_diff(pa_pre_iter, pa_cur_iter);
  double delta_pb = cs_cdo_blas_square_norm_pvsp_diff(pb_pre_iter, pb_cur_iter);

  cs_iter_algo_update_residual(algo, sqrt(delta_pa + delta_pb));

  /* Update the convergence members */

  cs_sles_convergence_state_t
    cvg_status = cs_iter_algo_update_cvg_tol_auto(algo);

  /* Monitoring */

  if (algo->verbosity > 1) {

    cs_iter_algo_log_cvg(algo, "# GWF.TPF");

    if (algo->verbosity > 2)
      cs_log_printf(CS_LOG_DEFAULT,
                    "\t||D_Pa||=%10.6e, ||D_Pb||=%10.6e\n",
                    sqrt(delta_pa), sqrt(delta_pb));

  }

  return cvg_status;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new (unsteady) state for the groundwater flows module.
 *        Case of (miscible or immiscible) two-phase flows in porous media
 *        with a non-linear resolution relying on the Picard algorithm and
 *        a coupled algorithm
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step    pointer to a cs_time_step_t structure
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] pa           field pointer for the 1st pressure
 * \param[in, out] pb           field pointer for the 2nd pressure
 * \param[in, out] tpf          pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_compute_coupled_picard(const cs_mesh_t              *mesh,
                        const cs_cdo_connect_t       *connect,
                        const cs_cdo_quantities_t    *cdoq,
                        const cs_time_step_t         *time_step,
                        cs_flag_t                     option_flag,
                        cs_field_t                   *pa,
                        cs_field_t                   *pb,
                        cs_gwf_tpf_t                 *tpf)
{
  bool cur2prev = true;
  cs_flag_t  update_flag = CS_FLAG_CURRENT_TO_PREVIOUS;

  cs_iter_algo_t  *algo = tpf->nl_algo;
  assert(algo != nullptr);
  assert(tpf->nl_algo_type == CS_PARAM_NL_ALGO_PICARD);

  _init_non_linear_algo(pa, pb, tpf);

  cs_real_t *pa_kp1 = nullptr, *pb_kp1 = nullptr; /* at ^{n+1,k+1} */
  cs_real_t *pa_k = nullptr, *pb_k = nullptr;     /* at ^{n+1,k} */

  BFT_MALLOC(pa_k, cdoq->n_vertices, cs_real_t);
  BFT_MALLOC(pb_k, cdoq->n_vertices, cs_real_t);

  pa_kp1 = pa->val; /* val stores always the latest values */
  pb_kp1 = pb->val;

  do {

    /* current values: n+1,k */

    cs_array_real_copy(cdoq->n_vertices, pa->val, pa_k);
    cs_array_real_copy(cdoq->n_vertices, pb->val, pb_k);

    /* pa and pb values are update during the coupled resolution */

    cs_equation_system_solve(cur2prev, tpf->system);

    /* current values: n+1,k+1 now */

    if (cs_iter_algo_get_n_iter(algo) > 0 && tpf->nl_relax_factor < 1.0) {

      const double  relax = tpf->nl_relax_factor;
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
        pa_kp1[i] = pa_kp1[i]*relax + pa_k[i]*(1-relax);
        pb_kp1[i] = pb_kp1[i]*relax + pb_k[i]*(1-relax);
      }

    }

    /* Update the variables related to the groundwater flow system */

    cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                      update_flag, option_flag,
                      tpf);

    /* After the first resolution, no cur2prev anymore */

    cur2prev = false;
    update_flag = 0;

  } while (_check_cvg_nl(tpf->nl_algo_type, pa_k, pa_kp1, pb_k, pb_kp1,
                         algo) == CS_SLES_ITERATING);

  if (algo->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "# GWF.TPF.Picard (exit) n_iter: %4d residual: %9.6e\n",
                  cs_iter_algo_get_n_iter(algo),
                  cs_iter_algo_get_residual(algo));

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_check_warning(__func__,
                             tpf->system->param->name,
                             cs_param_get_nl_algo_label(tpf->nl_algo_type),
                             algo);

  /* Free temporary arrays and structures */

  BFT_FREE(pa_k);
  BFT_FREE(pb_k);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new (unsteady) state for the groundwater flows module.
 *        Case of (miscible or immiscible) two-phase flows in porous media
 *        with a non-linear resolution relying on the Anderson algorithm and
 *        a coupled algorithm
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step    pointer to a cs_time_step_t structure
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] pa           field pointer for the 1st pressure
 * \param[in, out] pb           field pointer for the 2nd pressure
 * \param[in, out] tpf          pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_compute_coupled_anderson(const cs_mesh_t              *mesh,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *cdoq,
                          const cs_time_step_t         *time_step,
                          cs_flag_t                     option_flag,
                          cs_field_t                   *pa,
                          cs_field_t                   *pb,
                          cs_gwf_tpf_t                 *tpf)
{
  bool  cur2prev = false;
  cs_flag_t  update_flag = 0;   /* No current to previous operation */
  cs_iter_algo_t  *algo = tpf->nl_algo;
  assert(algo != nullptr);

  _init_non_linear_algo(pa, pb, tpf);

  /* First resolution */

  cs_equation_system_solve(true, tpf->system); /* cur2prev = true */

  /* Update the variables related to the groundwater flow system */

  cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                    CS_FLAG_CURRENT_TO_PREVIOUS,
                    option_flag,
                    tpf);

  /* A first resolution has been done followed by an update */

  cs_real_t *pa_kp1 = nullptr, *pb_kp1 = nullptr; /* at ^{n+1,k+1} */
  cs_real_t *pa_k = nullptr, *pb_k = nullptr;     /* at ^{n+1,k} */

  BFT_MALLOC(pa_k, 2*cdoq->n_vertices, cs_real_t);
  pb_k = pa_k + cdoq->n_vertices;

  cs_array_real_copy(cdoq->n_vertices, pa->val_pre, pa_k);
  cs_array_real_copy(cdoq->n_vertices, pb->val_pre, pb_k);

  /* One needs only one array gathering the two main pressures */

  BFT_MALLOC(pa_kp1, 2*cdoq->n_vertices, cs_real_t);
  pb_kp1 = pa_kp1 + cdoq->n_vertices;

  cs_array_real_copy(cdoq->n_vertices, pa->val, pa_kp1);
  cs_array_real_copy(cdoq->n_vertices, pb->val, pb_kp1);

  /* Main non-linear loop */

  while (_check_cvg_nl(tpf->nl_algo_type,
                       pa_k, pa_kp1, pb_k, pb_kp1,
                       algo) == CS_SLES_ITERATING) {

    cs_array_real_copy(cdoq->n_vertices, pa_kp1, pa_k);
    cs_array_real_copy(cdoq->n_vertices, pb_kp1, pb_k);

    /* Update the variables related to the groundwater flow system.
     * In case of an Anderson acceleration, pa_kp1 and pb_kp1 may be
     * updated */

    if (cs_iter_algo_get_n_iter(algo) >= tpf->anderson_param.starting_iter) {

      cs_array_real_copy(cdoq->n_vertices, pa_kp1, pa->val);
      cs_array_real_copy(cdoq->n_vertices, pb_kp1, pb->val);

    }

    cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                      update_flag,
                      option_flag,
                      tpf);

    /* Build and solve the linear system related to the coupled system of
       equations. First call: current --> previous and then no operation */

    cs_equation_system_solve(cur2prev, tpf->system);

    cs_array_real_copy(cdoq->n_vertices, pa->val, pa_kp1);
    cs_array_real_copy(cdoq->n_vertices, pb->val, pb_kp1);

  } /* Until convergence */

  /* Get the soil state up to date with the last compute values for the
     liquid and gas pressures */

  cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                    update_flag,
                    option_flag,
                    tpf);

  if (algo->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "# GWF.TPF.AndersonAcc (exit) n_iter: %4d residual: %9.6e\n",
                  cs_iter_algo_get_n_iter(algo),
                  cs_iter_algo_get_residual(algo));

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_check_warning(__func__,
                             tpf->system->param->name,
                             cs_param_get_nl_algo_label(tpf->nl_algo_type),
                             algo);

  /* Free temporary arrays and structures */

  BFT_FREE(pa_k);
  BFT_FREE(pa_kp1);
  cs_iter_algo_release_anderson_arrays((cs_iter_algo_aac_t *)algo->context);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new (unsteady) state for the groundwater flows module.
 *        Case of (miscible or immiscible) two-phase flows in porous media.
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      time_step    pointer to a cs_time_step_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] tpf          pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_compute_segregated(const cs_mesh_t              *mesh,
                    const cs_time_step_t         *time_step,
                    const cs_cdo_connect_t       *connect,
                    const cs_cdo_quantities_t    *cdoq,
                    cs_flag_t                     option_flag,
                    cs_gwf_tpf_t                 *tpf)
{
  assert(tpf->use_incremental_solver);

  /* Copy the state for time t^n into field->val_pre */

  cs_field_current_to_previous(tpf->g_pressure);
  cs_field_current_to_previous(tpf->l_pressure);

  /* Since the cur2prev operation has been done, avoid to do it again */

  bool cur2prev = false;        /* Done just above */
  cs_flag_t  update_flag = 0;   /* No current to previous operation */

  /* Initialize the non-linear algorithm */

  cs_iter_algo_t  *algo = tpf->nl_algo;
  assert(algo != nullptr);

  cs_iter_algo_reset(algo);

  double  normalization = cs_cdo_blas_square_norm_pvsp(tpf->c_pressure->val);

  if (normalization < cs_math_zero_threshold)
    normalization = 1.0;
  else
    normalization = sqrt(normalization);

  cs_iter_algo_set_normalization(algo, normalization);

  /* First solve step:
   * 1. Solve the equation associated to Pl --> Pl^(n+1,1) knowing
   *    Pl^(n+1,0) = Pl^n, Pg^n and thus Pc^n and Sl^n
   *    --> Sl^(n+1,0) - Sl^n = 0 --> no source term
   * 2. Solve the equation associated to Pg --> Pg^(n+1,1) knowing
   *
   */

  cs_equation_solve(cur2prev, mesh, tpf->w_eq);

  cs_equation_solve(cur2prev, mesh, tpf->h_eq);

  /* Update the variables related to the groundwater flow system */

  cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                    CS_FLAG_CURRENT_TO_PREVIOUS, /* Force this operation */
                    option_flag,
                    tpf);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
  char  label[32];
  int  n_iter = cs_iter_algo_get_n_iter(algo);

  sprintf(label, "Pl_iter%02d", n_iter);
  cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                           CS_POST_WRITER_DEFAULT,
                           label,
                           1,
                           false,
                           false,
                           CS_POST_TYPE_cs_real_t,
                           cs_field_by_name("liquid_pressure")->val,
                           time_step);

  sprintf(label, "Pg_iter%02d", n_iter);
  cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                           CS_POST_WRITER_DEFAULT,
                           label,
                           1,
                           false,
                           false,
                           CS_POST_TYPE_cs_real_t,
                           cs_field_by_name("gas_pressure")->val,
                           time_step);

  sprintf(label, "Pc_iter%02d", n_iter);
  cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                           CS_POST_WRITER_DEFAULT,
                           label,
                           1,
                           false,
                           false,
                           CS_POST_TYPE_cs_real_t,
                           cs_field_by_name("capillarity_pressure")->val,
                           time_step);
#endif

  cs_real_t *dpc_kp1 = nullptr;
  BFT_MALLOC(dpc_kp1, mesh->n_vertices, cs_real_t);

  _get_capillarity_pressure_increment(mesh, tpf->w_eq, tpf->h_eq, dpc_kp1);

  while(_check_cvg_nl_dpc(dpc_kp1, algo) == CS_SLES_ITERATING) {

    /* Solve step */

    cs_equation_solve(cur2prev, mesh, tpf->w_eq);

    cs_equation_solve(cur2prev, mesh, tpf->h_eq);

    /* Update the variables related to the groundwater flow system */

    cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                      update_flag,
                      option_flag,
                      tpf);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
    n_iter = cs_iter_algo_get_n_iter(algo);

    sprintf(label, "Pl_iter%02d", n_iter);
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_DEFAULT,
                             label,
                             1,
                             false,
                             false,
                             CS_POST_TYPE_cs_real_t,
                             cs_field_by_name("liquid_pressure")->val,
                             time_step);

    sprintf(label, "Pg_iter%02d", n_iter);
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_DEFAULT,
                             label,
                             1,
                             false,
                             false,
                             CS_POST_TYPE_cs_real_t,
                             cs_field_by_name("gas_pressure")->val,
                             time_step);

    sprintf(label, "Pc_iter%02d", n_iter);
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_DEFAULT,
                             label,
                             1,
                             false,
                             false,
                             CS_POST_TYPE_cs_real_t,
                             cs_field_by_name("capillarity_pressure")->val,
                             time_step);
#endif

    _get_capillarity_pressure_increment(mesh, tpf->w_eq, tpf->h_eq, dpc_kp1);

  } /* while not converged */

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_check_warning(__func__,
                             "Segregated incremental TPF solver",
                             cs_param_get_nl_algo_label(tpf->nl_algo_type),
                             algo);

  BFT_FREE(dpc_kp1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the different blocks building the coupled system
 *
 * \param[in, out] tpf     model context. Point to a cs_gwf_tpf_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_set_coupled_system(cs_gwf_tpf_t    *tpf)
{
  /* Define the coupled system of equations */
  /* -------------------------------------- */

  cs_equation_param_t  *b00_w_eqp = cs_equation_get_param(tpf->w_eq);
  cs_equation_param_t  *b11_h_eqp = cs_equation_get_param(tpf->h_eq);

  _set_default_eqp_settings(b00_w_eqp);
  _set_default_eqp_settings(b11_h_eqp);

  /* Create the (0,1)-block related to the water in the gas phase */

  tpf->b01_w_eqp = cs_equation_param_create("block01_w_eq",
                                           CS_EQUATION_TYPE_GROUNDWATER,
                                           1,
                                           CS_BC_SYMMETRY);

  _set_default_eqp_settings(tpf->b01_w_eqp);

  /* Create the (1,0)-block related to the hydrogen in the liquid phase */

  tpf->b10_h_eqp = cs_equation_param_create("block10_h_eq",
                                           CS_EQUATION_TYPE_GROUNDWATER,
                                           1,
                                           CS_BC_SYMMETRY);

  _set_default_eqp_settings(tpf->b10_h_eqp);

  /* Add a 2x2 system of coupled equations and define each block */

  tpf->system = cs_equation_system_add("two_phase_flow_porous_media",
                                      2,   /* system size */
                                      1);  /* scalar-valued block */

  /* Set all the blocks in the coupled system */

  cs_equation_system_assign_equation(0, tpf->w_eq, tpf->system);  /* (0,0) */
  cs_equation_system_assign_equation(1, tpf->h_eq, tpf->system);  /* (1,1) */
  cs_equation_system_assign_param(0, 1, tpf->b01_w_eqp, tpf->system);
  cs_equation_system_assign_param(1, 0, tpf->b10_h_eqp, tpf->system);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the model context
 *
 *        Case of a two-phase flows model in porous media using a coupled
 *        solver and the couple (Pc, Pg) as main unknowns
 *
 * \param[in, out] tpf         pointer to the model context structure
 * \param[in, out] perm_type   type of permeability to handle
 */
/*----------------------------------------------------------------------------*/

static void
_init_pcpg_coupled_solver(cs_gwf_tpf_t            *tpf,
                          cs_property_type_t       perm_type)
{
  /* Equations and the coupled system
   * ================================
   *
   * Notations are the following :
   * - Two phases: Liquid phase denoted by "l" and gaseous phase denoted by "g"
   * - indice "c" refers to the capillarity pressure
   * - Two components: water denoted by "w" and a gaseous component (let's say
   *   hydrogen) denoted by "h". The gaseous component is present in the two
   *   phases whereas water is only considered in the liquid phase.
   *
   * The resulting linear algebraic system (one applies a linearization) is
   * defined as follows:
   *
   *                              cap.   liq
   * water mass conservation    | b00  | b01 ||P_c|   | b_w |
   *                            |------|-----||---| = |-----|
   * hydrogen mass conservation | b10  | b11 ||P_g|   | b_h |
   *
   * This is a coupled system. Coupling terms are collected inside M_01 and M_10
   */

  tpf->w_eq = cs_equation_add("w_conservation",        /* equation name */
                              "capillarity_pressure",  /* variable name */
                              CS_EQUATION_TYPE_GROUNDWATER,
                              1,
                              CS_BC_SYMMETRY);

  tpf->h_eq = cs_equation_add("h_conservation",   /* equation name */
                              "gas_pressure",     /* variable name */
                              CS_EQUATION_TYPE_GROUNDWATER,
                              1,
                              CS_BC_SYMMETRY);

  cs_equation_param_t  *b00_w_eqp = cs_equation_get_param(tpf->w_eq);
  cs_equation_param_t  *b11_h_eqp = cs_equation_get_param(tpf->h_eq);

  _set_coupled_system(tpf);

  /* Properties
   * ========== */

  /* Conservation of the mass of water
   * --------------------------------- */

  /* (0,0)-block */

  tpf->time_wc_pty = cs_property_add("time_wc_pty", CS_PROPERTY_ISO);
  cs_equation_add_time(b00_w_eqp, tpf->time_wc_pty);

  tpf->diff_wc_pty = cs_property_add("diff_wc_pty", perm_type);
  cs_equation_add_diffusion(b00_w_eqp, tpf->diff_wc_pty);

  /* (0,1)-block */

  tpf->diff_wg_pty = cs_property_add("diff_wg_pty", perm_type);
  cs_equation_add_diffusion(tpf->b01_w_eqp, tpf->diff_wg_pty);

  /* Conservation of the mass of component mainly present in the gas phase
   * --------------------------------------------------------------------- */

  /* (1,0)-block */

  tpf->time_hc_pty = cs_property_add("time_hc_pty", CS_PROPERTY_ISO);
  cs_equation_add_time(tpf->b10_h_eqp, tpf->time_hc_pty);

  if (tpf->is_miscible) {

    tpf->diff_hc_pty = cs_property_add("diff_hc_pty", perm_type);
    cs_equation_add_diffusion(tpf->b10_h_eqp, tpf->diff_hc_pty);

  }

  /* (1,1)-block */

  tpf->time_hg_pty = cs_property_add("time_hg_pty", CS_PROPERTY_ISO);
  cs_equation_add_time(b11_h_eqp, tpf->time_hg_pty);

  if (_pcpg_coupled_solver_needs_diffusion(tpf)) {

    tpf->diff_hg_pty = cs_property_add("diff_hg_pty", perm_type);
    cs_equation_add_diffusion(b11_h_eqp, tpf->diff_hg_pty);

  }

  if (!tpf->use_diffusion_view_for_darcy) /* Advection viewpoint */
    cs_equation_add_advection(b11_h_eqp, tpf->t_darcy->adv_field);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize the setup stage (allocation of arrays for instance) in the
 *        case of a two-phase flows model in porous media using a coupled
 *        solver and the couple (Pc, Pg) as main unknowns
 *
 * \param[in]      connect     set of additional connectivities for CDO
 * \param[in, out] tpf         pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_finalize_setup_pcpg_coupled_solver(const cs_cdo_connect_t  *connect,
                                    cs_gwf_tpf_t            *tpf)
{
  const cs_lnum_t  n_cells = connect->n_cells;

  /* Conservation of the mass of water
   * --------------------------------- */

  /* Unsteady term associated to Pc */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->time_wc_pty);

  /* Diffusion term associated to Pc */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_wc_pty);

  /* Diffusion term associated to Pg */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_wg_pty);

  /* Conservation of the mass of component mainly present in the gas phase
   * --------------------------------------------------------------------- */

  /* Unsteady term associated to Pc */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->time_hc_pty);

  if (tpf->is_miscible) /* Diffusion term associated to Pc */
    _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_hc_pty);

  if (_pcpg_coupled_solver_needs_diffusion(tpf)) {

    /* Diffusion term associated to Pg */

    _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_hg_pty);

  }

  /* Unsteady term associated to Pg */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->time_hg_pty);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the model context
 *
 *        Case of a two-phase flows model in porous media using a coupled
 *        solver and the couple (Pc, Pg) as main unknowns
 *
 * \param[in, out] tpf         pointer to the model context structure
 * \param[in, out] perm_type   type of permeability to handle
 */
/*----------------------------------------------------------------------------*/

static void
_init_plpc_coupled_solver(cs_gwf_tpf_t            *tpf,
                          cs_property_type_t       perm_type)
{
  /* Equations and the coupled system
   * ================================
   *
   * Notations are the following :
   * - Two phases: Liquid phase denoted by "l" and gaseous phase denoted by "g"
   * - indice "c" refers to the capillarity pressure
   * - Two components: water denoted by "w" and a gaseous component (let's say
   *   hydrogen) denoted by "h". The gaseous component is present in the two
   *   phases whereas water is only considered in the liquid phase.
   *
   * The resulting linear algebraic system (one applies a linearization) is
   * defined as follows:
   *
   *                              liq.   cap.
   * water mass conservation    | b00  | b01 ||P_l|   | b_w |
   *                            |------|-----||---| = |-----|
   * hydrogen mass conservation | b10  | b11 ||P_c|   | b_h |
   *
   * This is a coupled system. Coupling terms are collected inside b01 and b10
   */

  tpf->w_eq = cs_equation_add("w_conservation",   /* equation name */
                              "liquid_pressure",  /* variable name */
                              CS_EQUATION_TYPE_GROUNDWATER,
                              1,
                              CS_BC_SYMMETRY);

  tpf->h_eq = cs_equation_add("h_conservation",       /* equation name */
                              "capillarity_pressure", /* variable name */
                              CS_EQUATION_TYPE_GROUNDWATER,
                              1,
                              CS_BC_SYMMETRY);

  cs_equation_param_t  *b00_w_eqp = cs_equation_get_param(tpf->w_eq);
  cs_equation_param_t  *b11_h_eqp = cs_equation_get_param(tpf->h_eq);

  _set_coupled_system(tpf);

  /* Properties
   * ========== */

  /* Conservation of the mass of water
   * --------------------------------- */

  /* (0,0)-block */

  tpf->diff_wl_pty = cs_property_add("diff_wl_pty", perm_type);
  cs_equation_add_diffusion(b00_w_eqp, tpf->diff_wl_pty);

  /* (0,1)-block */

  tpf->time_wc_pty = cs_property_add("time_wc_pty", CS_PROPERTY_ISO);
  cs_equation_add_time(tpf->b01_w_eqp, tpf->time_wc_pty);

  /* Conservation of the mass of component mainly present in the gas phase
   * --------------------------------------------------------------------- */

  /* (1,1)-block */

  tpf->time_hc_pty = cs_property_add("time_hc_pty", CS_PROPERTY_ISO);
  cs_equation_add_time(b11_h_eqp, tpf->time_hc_pty);

  tpf->diff_hc_pty = cs_property_add("diff_hc_pty", perm_type);
  cs_equation_add_diffusion(b11_h_eqp, tpf->diff_hc_pty);

  /* (1,0)-block */

  tpf->time_hl_pty = cs_property_add("time_hl_pty", CS_PROPERTY_ISO);
  cs_equation_add_time(tpf->b10_h_eqp, tpf->time_hl_pty);

  tpf->diff_hl_pty = cs_property_add("diff_hl_pty", perm_type);
  cs_equation_add_diffusion(tpf->b10_h_eqp, tpf->diff_hl_pty);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize the setup stage (allocation of arrays for instance) in the
 *        case of a two-phase flows model in porous media using a coupled
 *        solver and the couple (Pl, Pc) as main unknowns
 *
 * \param[in]      connect     set of additional connectivities for CDO
 * \param[in, out] tpf         pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_finalize_setup_plpc_coupled_solver(const cs_cdo_connect_t  *connect,
                                    cs_gwf_tpf_t            *tpf)
{
  const cs_lnum_t  n_cells = connect->n_cells;

  /* Conservation of the mass of water
   * --------------------------------- */

  /* Unsteady term associated to Pc */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->time_wc_pty);

  /* Diffusion term associated to Pl */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_wl_pty);

  /* Conservation of the mass of component mainly present in the gas phase
   * --------------------------------------------------------------------- */

  /* Unsteady term associated to Pc */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->time_hc_pty);

  /* Diffusion term associated to Pc */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_hc_pty);

  /* Unsteady term associated to Pl */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->time_hl_pty);

  /* Diffusion term associated to Pl */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_hl_pty);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the model context
 *
 *        Case of a two-phase flows model in porous media using a segregated
 *        solver and the couple (Pl, Pg) as main unknowns
 *
 * \param[in, out] tpf         pointer to the model context structure
 * \param[in, out] perm_type   type of permeability to handle
 */
/*----------------------------------------------------------------------------*/

static void
_init_plpg_segregated_solver(cs_gwf_tpf_t            *tpf,
                             cs_property_type_t       perm_type)
{
 /* Equations
  * =========
  *
  * Notations are the following :
  * - Two phases: Liquid phase denoted by "l" and gaseous phase denoted by "g"
  * - indice "c" refers to the capillarity pressure
  * - Two components: water denoted by "w" and a gaseous component (let's say
  *   hydrogen) denoted by "h". The gaseous component is present in the two
  *   phases whereas water is only considered in the liquid phase.
  */

  tpf->w_eq = cs_equation_add("w_conservation",   /* equation name */
                              "liquid_pressure",  /* variable name */
                              CS_EQUATION_TYPE_GROUNDWATER,
                              1,
                              CS_BC_SYMMETRY);

  tpf->h_eq = cs_equation_add("h_conservation",  /* equation name */
                              "gas_pressure",    /* variable name */
                              CS_EQUATION_TYPE_GROUNDWATER,
                              1,
                              CS_BC_SYMMETRY);

  cs_equation_param_t  *w_eqp = cs_equation_get_param(tpf->w_eq);
  cs_equation_param_t  *h_eqp = cs_equation_get_param(tpf->h_eq);

  /* Properties
   * ========== */

  /* Conservation of the mass of water
   * --------------------------------- */

  /* used both for the unsteady term and the computation of the source term.
   * \partial_t Pc = \partial_t Pg - \partial_t Pl */

  tpf->time_wc_pty = cs_property_add("time_wl_pty", CS_PROPERTY_ISO);
  cs_equation_add_time(h_eqp, tpf->time_wc_pty);

  tpf->diff_wl_pty = cs_property_add("diff_wl_pty", perm_type);
  cs_equation_add_diffusion(w_eqp, tpf->diff_wl_pty);

  /* Conservation of the mass of component mainly present in the gas phase
   * --------------------------------------------------------------------- */

  tpf->time_hc_pty = cs_property_add("time_hc_pty", CS_PROPERTY_ISO);
  cs_equation_add_time(h_eqp, tpf->time_hc_pty);

  tpf->diff_hc_pty = cs_property_add("diff_hc_pty", perm_type);
  cs_equation_add_diffusion(h_eqp, tpf->diff_hc_pty);

  tpf->reac_h_pty = cs_property_add("reac_h_pty", CS_PROPERTY_ISO);
  cs_equation_add_reaction(h_eqp, tpf->reac_h_pty);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize the setup stage (allocation of arrays for instance) in the
 *        case of a two-phase flows model in porous media using a segregated
 *        solver and the couple (Pl, Pg) as main unknowns
 *
 * \param[in]      connect     set of additional connectivities for CDO
 * \param[in, out] tpf         pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_finalize_setup_plpg_segregated_solver(const cs_cdo_connect_t  *connect,
                                       cs_gwf_tpf_t            *tpf)
{
  cs_equation_param_t  *w_eqp = cs_equation_get_param(tpf->w_eq);
  cs_equation_param_t  *h_eqp = cs_equation_get_param(tpf->h_eq);

  const cs_lnum_t  n_cells = connect->n_cells;

  /* Conservation of the mass of water
   * --------------------------------- */

  /* Unsteady term */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->time_wc_pty);

  /* Diffusion term */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_wl_pty);

  /* Source term */

  BFT_MALLOC(tpf->srct_w_array, n_cells, cs_real_t);
  cs_array_real_fill_zero(n_cells, tpf->srct_w_array);

  cs_equation_add_source_term_by_array(w_eqp,
                                       nullptr, /* all cells */
                                       cs_flag_primal_cell,
                                       tpf->srct_w_array,
                                       false, /* xdef not owner */
                                       true); /* full length */

  /* Conservation of the mass of component mainly present in the gas phase
   * --------------------------------------------------------------------- */

  /* Unsteady term (associated to Pg in the case of a segregated solver) */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->time_hc_pty);

  /* Diffusion property in the hydrogen equation (associated to Pg in the
     case of a segregated solver) */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_hc_pty);

  /* Diffusion property in the hydrogen equation associated to Pl. This is
     used to define a source term in the case of a segregated solver. */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_hl_pty);

  /* Reaction term */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->reac_h_pty);

  /* Source term (second one; the first is associated to a stiffness term and
     relies on a builder hook function) */

  BFT_MALLOC(tpf->srct_h_array, n_cells, cs_real_t);
  cs_array_real_fill_zero(n_cells, tpf->srct_h_array);

  cs_equation_add_source_term_by_array(h_eqp,
                                       nullptr, /* all cells */
                                       cs_flag_primal_cell,
                                       tpf->srct_h_array,
                                       false, /* xdef not owner */
                                       true); /* full length */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and initialize the model context structure for two-phase
 *        flows in a porous media
 *
 * \param[in] model       type of physical modelling
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tpf_t *
cs_gwf_tpf_create(cs_gwf_model_type_t      model)
{
  assert(model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE ||
         model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE);

  cs_gwf_tpf_t *tpf = nullptr;

  BFT_MALLOC(tpf, 1, cs_gwf_tpf_t);

  /* Set to nullptr since one has to know if a coupled or segregated solver is
     used or what couple of variables is considered */

  tpf->w_eq      = nullptr;
  tpf->h_eq      = nullptr;
  tpf->b01_w_eqp = nullptr;
  tpf->b10_h_eqp = nullptr;
  tpf->system    = nullptr;

  /* Advection fields */
  /* ---------------- */

  tpf->l_darcy = nullptr;
  tpf->g_darcy = nullptr;
  tpf->t_darcy = nullptr;

  /* Properties
   * ---------- */

  tpf->krl_pty  = nullptr;
  tpf->krg_pty  = nullptr;
  tpf->lsat_pty = nullptr;
  tpf->lcap_pty = nullptr;

  /* Note: Adding the properties related to the permeability (diffusion) is
   * postponed since one has to know which type of permeability is considered
   * (iso, ortho, or anisotropic) and this information is a result of the type
   * of soils which have been added.
   */

  /* Properties related to the water equation  */

  tpf->time_wc_pty = nullptr;
  tpf->diff_wl_pty = nullptr;
  tpf->diff_wc_pty = nullptr; /* only if (Pc, Pg) is used */
  tpf->diff_wg_pty = nullptr; /* only if (Pc, Pg) is used */

  /* Properties related to the hydrogen equation  */

  tpf->time_hc_pty = nullptr;
  tpf->time_hl_pty = nullptr; /* only if (Pc, Pl) is used */
  tpf->time_hg_pty = nullptr; /* only if (Pc, Pg) is used */

  tpf->diff_hc_pty = nullptr;
  tpf->diff_hl_pty = nullptr; /* only if (Pc, Pl) is used */
  tpf->diff_hg_pty = nullptr; /* only if (Pc, Pl) is used */

  tpf->reac_h_pty = nullptr;

  /* Fields */
  /* ------ */

  tpf->c_pressure   = nullptr;
  tpf->l_pressure   = nullptr;
  tpf->g_pressure   = nullptr;
  tpf->l_saturation = nullptr;

  /* Arrays */
  /* ------ */

  /* The properties will be defined using arrays.
   * Store these arrays of property values associated to equation terms */

  tpf->srct_w_array = nullptr;
  tpf->srct_h_array = nullptr;

  /* Model parameters (default values) */
  /* ---------------- */

  if (model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE) {

    tpf->is_miscible = true;
    tpf->l_diffusivity_h = 0;      /* immiscible case */
    tpf->henry_constant = 1e-7;    /* default value */

  }
  else { /* immiscible case */

    tpf->is_miscible = false;
    tpf->l_diffusivity_h = 0;
    tpf->henry_constant = 0;

  }

  tpf->l_mass_density = 1000;
  tpf->l_viscosity = 1e-3;
  tpf->g_viscosity = 2e-5;
  tpf->h_molar_mass = 3e-3;
  tpf->ref_temperature = 280;    /* in Kelvin */

  /* Numerical parameters (default values) */

  tpf->approx_type = CS_GWF_TPF_APPROX_PC_CELL_VERTEX_AVERAGE;
  tpf->cell_weight = 0.5;
  tpf->upwind_weight = 0.5;
  tpf->solver_type = CS_GWF_TPF_SOLVER_PLPC_COUPLED;
  tpf->use_coupled_solver = true;
  tpf->use_incremental_solver = false;
  tpf->use_diffusion_view_for_darcy = true;

  tpf->nl_algo_type = CS_PARAM_NL_ALGO_PICARD;

  tpf->nl_relax_factor = 1.0;
  tpf->nl_cvg_param.n_max_iter = 100;
  tpf->nl_cvg_param.rtol = 1e-6;
  tpf->nl_cvg_param.atol = 1e-12;
  tpf->nl_cvg_param.dtol = 1e8; /* Problems with big capillarity steps may
                                   induce locally increment(s) with a strong
                                   variation */

  tpf->anderson_param.n_max_dir = 5;
  tpf->anderson_param.starting_iter = 3;
  tpf->anderson_param.max_cond = -1; /* No test by default */
  tpf->anderson_param.beta = 1.0;    /* No damping by default */
  tpf->anderson_param.dp_type = CS_PARAM_DOTPROD_EUCLIDEAN;

  tpf->nl_algo = nullptr;

  return tpf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to the modelling of two-phase
 *        flows in a porous media
 *
 * \param[in, out] p_tpf   pointer of pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_free(cs_gwf_tpf_t    **p_tpf)
{
  if (cs_gwf_tpf_time_plot != nullptr)
    cs_time_plot_finalize(&cs_gwf_tpf_time_plot);

  if (p_tpf == nullptr)
    return;
  if (*p_tpf == nullptr)
    return;

  cs_gwf_tpf_t  *tpf = *p_tpf;

  /* System of equations are freed elsewhere (just after having freed
     cs_equation_t structures) */

  cs_gwf_darcy_flux_free(&(tpf->l_darcy));
  cs_gwf_darcy_flux_free(&(tpf->g_darcy));
  cs_gwf_darcy_flux_free(&(tpf->t_darcy));

  BFT_FREE(tpf->srct_w_array);
  BFT_FREE(tpf->srct_h_array);

  /* Free the structure handling the convergence of the non-linear algorithm */

  cs_iter_algo_free(&(tpf->nl_algo));

  BFT_FREE(tpf);
  *p_tpf = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the model context of two-phase flows.
 *        Common to the different sub-models relying on two-phase flows.
 *
 * \param[in] tpf      pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_log_setup(cs_gwf_tpf_t          *tpf)
{
  if (tpf == nullptr)
    return;

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Reference temperature: %5.2f K\n",
                tpf->ref_temperature);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | %14s | mass density: %5.3e, viscosity: %5.3e\n",
                "Water", tpf->l_mass_density, tpf->l_viscosity);

  if (tpf->is_miscible) {

    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | %14s | viscosity: %5.3e, saturated diffusivity"
                  " in the liquid phase: %5.3e, molar mass: %5.3e\n",
                  "Gas component", tpf->g_viscosity, tpf->l_diffusivity_h,
                  tpf->h_molar_mass);
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Henry constant: %5.3e\n", tpf->henry_constant);

  }
  else {

    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | %14s | viscosity: %5.3e, molar mass: %5.3e\n",
                  "Gas component", tpf->g_viscosity, tpf->h_molar_mass);

  }

  switch(tpf->approx_type) {

  case CS_GWF_TPF_APPROX_PC_CELL_AVERAGE:
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Approximation: %s\n",
                  "CS_GWF_TPF_APPROX_PC_CELL_AVERAGE");
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Upwind weight: %.2f\n",
                  tpf->upwind_weight);
    break;

  case CS_GWF_TPF_APPROX_PC_CELL_VERTEX_AVERAGE:
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Approximation: V+C average (c:%.3f | v:%.3f)\n",
                  tpf->cell_weight, 1-tpf->cell_weight);
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Upwind weight: %.2f\n",
                  tpf->upwind_weight);
    break;

  case CS_GWF_TPF_APPROX_PC_EDGE_AVERAGE:
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Approximation: %s\n",
                  "CS_GWF_TPF_APPROX_PC_EDGE_AVERAGE");
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Upwind weight: %.2f\n",
                  tpf->upwind_weight);
    break;

  case CS_GWF_TPF_APPROX_PC_VERTEX_AVERAGE:
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Approximation: %s\n",
                  "CS_GWF_TPF_APPROX_PC_VERTEX_AVERAGE");
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Upwind weight: %.2f\n",
                  tpf->upwind_weight);
    break;

  case CS_GWF_TPF_APPROX_VERTEX_SUBCELL:
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Approximation: %s\n",
                  "CS_GWF_TPF_APPROX_VERTEX_SUBCELL");
    break;

  default:
    break;
  }

  switch(tpf->solver_type) {

  case CS_GWF_TPF_SOLVER_PCPG_COUPLED:
    cs_log_printf(CS_LOG_SETUP, "  * GWF | (Pc, Pg) coupled solver\n");
    if (tpf->use_diffusion_view_for_darcy)
      cs_log_printf(CS_LOG_SETUP, "  * GWF | Diffusive view for Darcy terms\n");
    else
      cs_log_printf(CS_LOG_SETUP, "  * GWF | Advective view for Darcy terms\n");
    break;
  case CS_GWF_TPF_SOLVER_PLPC_COUPLED:
    cs_log_printf(CS_LOG_SETUP, "  * GWF | (Pl, Pc) coupled solver\n");
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Diffusive view for Darcy terms\n");
    break;
  case CS_GWF_TPF_SOLVER_PLPG_SEGREGATED:
    cs_log_printf(CS_LOG_SETUP, "  * GWF | (Pl, Pg) segregated solver\n");
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Diffusive view for Darcy terms\n");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid setting", __func__);

  }

  if (tpf->use_incremental_solver)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Incremental solver\n");

  cs_gwf_darcy_flux_log(tpf->l_darcy);
  cs_gwf_darcy_flux_log(tpf->g_darcy);

  if (!tpf->use_diffusion_view_for_darcy)
    cs_gwf_darcy_flux_log(tpf->t_darcy);

  if (tpf->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

    cs_log_printf(CS_LOG_SETUP, "  * GWF | Non-linear algo.: %s\n",
                  cs_param_get_nl_algo_name(tpf->nl_algo_type));

    cs_log_printf(CS_LOG_SETUP, "  * GWF | Tolerances of non-linear algo:"
                  " rtol: %5.3e; atol: %5.3e; dtol: %5.3e; max_iter: %d\n",
                  tpf->nl_cvg_param.rtol, tpf->nl_cvg_param.atol,
                  tpf->nl_cvg_param.dtol, tpf->nl_cvg_param.n_max_iter);

    if (tpf->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON) {

      const cs_iter_algo_param_aac_t  aap = tpf->anderson_param;

      cs_log_printf(CS_LOG_SETUP, "  * GWF | Anderson param: max. dir: %d; "
                    " start: %d; drop. tol: %5.3e; relax: %5.3e\n",
                    aap.n_max_dir, aap.starting_iter, aap.max_cond, aap.beta);
      cs_log_printf(CS_LOG_SETUP, "  * GWF | Anderson param: Dot product: %s\n",
                    cs_param_get_dotprod_type_name(aap.dp_type));

    }

  } /* There is a non-linear algorithm */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the model context according to the settings done inside
 *        the function cs_user_model()
 *        Case of a two-phase flows model in porous media
 *
 * \param[in, out] tpf         pointer to the model context structure
 * \param[in, out] perm_type   type of permeability to handle
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init(cs_gwf_tpf_t            *tpf,
                cs_property_type_t       perm_type)
{
  if (tpf == nullptr)
    return;

  /* Property common to all solvers. This property is used to compute the Darcy
     flux in the gas phase */

  tpf->diff_g_pty = cs_property_add("diff_g_pty", perm_type);

  switch (tpf->approx_type) {

  case CS_GWF_TPF_APPROX_VERTEX_SUBCELL:
    tpf->krl_pty = cs_property_subcell_add("krl_pty", CS_PROPERTY_ISO);
    tpf->krg_pty = cs_property_subcell_add("krg_pty", CS_PROPERTY_ISO);
    tpf->lsat_pty = cs_property_subcell_add("lsat_pty", CS_PROPERTY_ISO);
    tpf->lcap_pty = cs_property_subcell_add("lcap_pty", CS_PROPERTY_ISO);
    break;

  default:
    tpf->krl_pty = cs_property_add("krl_pty", CS_PROPERTY_ISO);
    tpf->krg_pty = cs_property_add("krg_pty", CS_PROPERTY_ISO);
    tpf->lsat_pty = cs_property_add("lsat_pty", CS_PROPERTY_ISO);
    tpf->lcap_pty = cs_property_add("lcap_pty", CS_PROPERTY_ISO);
    break;

  }

  /* Darcy flux (one assumes a CDOVB space discretization) */

  cs_advection_field_status_t  adv_status =
    CS_ADVECTION_FIELD_GWF | CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;

  cs_adv_field_t  *l_adv_field = cs_advection_field_add("l_darcy_field",
                                                        adv_status);
  tpf->l_darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);
  tpf->l_darcy->adv_field = l_adv_field;

  cs_adv_field_t  *g_adv_field = cs_advection_field_add("g_darcy_field",
                                                        adv_status);
  tpf->g_darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);
  tpf->g_darcy->adv_field = g_adv_field;

  if (!tpf->use_diffusion_view_for_darcy) {

    cs_adv_field_t  *t_adv_field = cs_advection_field_add("t_darcy_field",
                                                          adv_status);
    tpf->t_darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);
    tpf->t_darcy->adv_field = t_adv_field;

  }

  /* Define the system of equations */
  /* ------------------------------ */

  switch(tpf->solver_type) {

  case CS_GWF_TPF_SOLVER_PCPG_COUPLED:
    _init_pcpg_coupled_solver(tpf, perm_type);
    break;

  case CS_GWF_TPF_SOLVER_PLPC_COUPLED:
    _init_plpc_coupled_solver(tpf, perm_type);
    break;

  case CS_GWF_TPF_SOLVER_PLPG_SEGREGATED:
    _init_plpg_segregated_solver(tpf, perm_type);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver type", __func__);
  }

  /* Settings for the treatment of the non-linearity */

  if (tpf->use_incremental_solver) {

    cs_equation_param_t  *w_eqp = cs_equation_get_param(tpf->w_eq);
    cs_equation_param_t  *h_eqp = cs_equation_get_param(tpf->h_eq);

    w_eqp->incremental_algo_type = CS_PARAM_NL_ALGO_PICARD;
    h_eqp->incremental_algo_type = CS_PARAM_NL_ALGO_PICARD;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial setup stage for two-phase flows in porous media. At this
 *        stage, all soils have been defined and equation parameters are set.
 *        Case of a miscible or immiscible model.
 *
 * \param[in]      post_flag   optional postprocessing request(s)
 * \param[in, out] tpf         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init_setup(cs_flag_t         post_flag,
                      cs_gwf_tpf_t     *tpf)
{
  if (tpf == nullptr)
    return;

  if (tpf->w_eq == nullptr || tpf->h_eq == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Equations are not defined for this model. Stop execution.\n",
              __func__);

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");

  /* Retrieve the pointers to the structure storing the equation parameters */

  cs_equation_param_t  *w_eqp = cs_equation_get_param(tpf->w_eq);
  cs_equation_param_t  *h_eqp = cs_equation_get_param(tpf->h_eq);

  assert(w_eqp->space_scheme == h_eqp->space_scheme);

  int loc_id = c_loc_id;
  if (w_eqp->space_scheme == CS_SPACE_SCHEME_CDOVB ||
      w_eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB)
    loc_id = v_loc_id;
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme", __func__);

  /* Add the variable fields (Keep always the previous state) */

  cs_equation_predefined_create_field(1, tpf->w_eq);
  cs_equation_predefined_create_field(1, tpf->h_eq);

  /* Set fields related to variables.
   * One has to be consistent with the location of DoFs for the w_eq and h_eq
   */

  switch(tpf->solver_type) {

  case CS_GWF_TPF_SOLVER_PCPG_COUPLED:
    tpf->c_pressure = cs_equation_get_field(tpf->w_eq);
    tpf->g_pressure = cs_equation_get_field(tpf->h_eq);
    tpf->l_pressure = cs_field_create("liquid_pressure",
                                      field_mask,
                                      loc_id,
                                      1,
                                      true); /* has_previous */

    cs_field_set_key_int(tpf->l_pressure, log_key, 1);
    cs_field_set_key_int(tpf->l_pressure, post_key, 1);
    break;

  case CS_GWF_TPF_SOLVER_PLPC_COUPLED:
    tpf->l_pressure = cs_equation_get_field(tpf->w_eq);
    tpf->c_pressure = cs_equation_get_field(tpf->h_eq);
    tpf->g_pressure = cs_field_create("gas_pressure",
                                      field_mask,
                                      loc_id,
                                      1,
                                      true); /* has_previous */

    cs_field_set_key_int(tpf->g_pressure, log_key, 1);
    cs_field_set_key_int(tpf->g_pressure, post_key, 1);
    break;

  case CS_GWF_TPF_SOLVER_PLPG_SEGREGATED:
    tpf->l_pressure = cs_equation_get_field(tpf->w_eq);
    tpf->g_pressure = cs_equation_get_field(tpf->h_eq);
    tpf->c_pressure = cs_field_create("capillarity_pressure",
                                      field_mask,
                                      loc_id,
                                      1,
                                      true); /* has_previous */

    cs_field_set_key_int(tpf->c_pressure, log_key, 1);
    cs_field_set_key_int(tpf->c_pressure, post_key, 1);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver type", __func__);
  }

  /* Create a liquid saturation field attached to cells: S_l (only keep the
     latest computed state) */

  int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;

  tpf->l_saturation = cs_field_create("liquid_saturation",
                                      pty_mask,
                                      c_loc_id,
                                      1,      /* dimension */
                                      false); /* has_previous */

  cs_field_set_key_int(tpf->l_saturation, log_key, 1);
  if (post_flag & CS_GWF_POST_LIQUID_SATURATION)
    cs_field_set_key_int(tpf->l_saturation, post_key, 1);

  /* Set the function for the computation of rhog_h */

  switch (tpf->approx_type) {

  case CS_GWF_TPF_APPROX_VERTEX_SUBCELL:
    compute_rhog_h = _rhog_h_vtx;
    compute_rhogl_h = _rhogl_h_vtx;
    break;

  default:
    if (fabs(tpf->upwind_weight) < FLT_MIN) { /* = 0 */
      compute_rhog_h = _rhog_h_cell_mean;
      compute_rhogl_h = _rhogl_h_cell_mean;
    }
    else {
      tpf->upwind_weight = fmax(0., fmin(1., tpf->upwind_weight));
      compute_rhog_h = _rhog_h_cell_upw;
      compute_rhogl_h = _rhogl_h_cell_upw;
    }
    break;

  }

  /* Handle automatic post-processings */

  if (post_flag & CS_GWF_POST_SOIL_MINMAX) {

    int  n_soils = cs_gwf_get_n_soils();

    /* There are 3 pressures (liquid, gas and capillarity), the liquid
       saturation and the liquid capacity. For each variable, two quantities
       are considered: the min. and the max. */

    int  n_outputs = 2 * CS_GWF_TPF_N_OUTPUT_VARS * n_soils;
    char  **labels;
    BFT_MALLOC(labels, n_outputs, char *);

    for (int i = 0; i < n_soils; i++) {

      const cs_zone_t  *z = cs_gwf_soil_get_zone(i);

      int  min_shift = CS_GWF_TPF_N_OUTPUT_VARS*i;
      int  max_shift = CS_GWF_TPF_N_OUTPUT_VARS*(n_soils + i);

      for (int j = 0; j < CS_GWF_TPF_N_OUTPUT_VARS; j++) {

        const char  *vname = _output_varnames[j];
        int  len = strlen(z->name) + strlen(vname) + 4;

        char *min_name = nullptr;
        BFT_MALLOC(min_name, len + 1, char);
        sprintf(min_name, "%s-%sMin", z->name, vname);
        labels[min_shift + j] = min_name;

        char *max_name = nullptr;
        BFT_MALLOC(max_name, len + 1, char);
        sprintf(max_name, "%s-%sMax", z->name, vname);
        labels[max_shift + j] = max_name;

      } /* Loop on variables */

    } /* Loop on soils */

    cs_gwf_tpf_time_plot = cs_time_plot_init_probe("gwf",
                                                   "",
                                                   CS_TIME_PLOT_CSV,
                                                   false,
                                                   180, /* flush time */
                                                   -1,
                                                   n_outputs,
                                                   nullptr,
                                                   nullptr,
                                                   (const char **)labels);

    for (int i = 0; i < n_outputs; i++)
      BFT_FREE(labels[i]);
    BFT_FREE(labels);

  } /* Min/Max output */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of two-phase flows in a porous media
 *        (miscible or immiscible case)
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]      flag      optional settings for the module
 * \param[in, out] tpf       pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_finalize_setup(const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *cdoq,
                          cs_flag_t                    flag,
                          cs_gwf_tpf_t                *tpf)
{
  CS_NO_WARN_IF_UNUSED(flag);   /* will be useful for gravity effect */

  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_adjacency_t  *c2v = connect->c2v;
  const cs_lnum_t  c2v_size = c2v->idx[n_cells];

  cs_equation_param_t  *w_eqp = cs_equation_get_param(tpf->w_eq);
  cs_equation_param_t  *h_eqp = cs_equation_get_param(tpf->h_eq);

  assert(cs_equation_get_type(tpf->w_eq) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(cs_equation_get_type(tpf->h_eq) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(w_eqp->space_scheme == h_eqp->space_scheme);

  if (w_eqp->space_scheme != CS_SPACE_SCHEME_CDOVB)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme", __func__);

  /* Set the Darcian flux (in the volume and at the boundary) */
  /* -------------------------------------------------------- */

  cs_gwf_darcy_flux_define(connect, cdoq, w_eqp->space_scheme,
                           tpf, _update_darcy_l,
                           tpf->l_darcy);

  if (tpf->use_coupled_solver)
    cs_gwf_darcy_flux_define(connect, cdoq, h_eqp->space_scheme,
                             tpf, _update_darcy_g_coupled,
                             tpf->g_darcy);
  else
    cs_gwf_darcy_flux_define(connect, cdoq, h_eqp->space_scheme,
                             tpf, _update_darcy_g_segregated,
                             tpf->g_darcy);

  if (!tpf->use_diffusion_view_for_darcy)
    cs_gwf_darcy_flux_define(connect, cdoq, w_eqp->space_scheme,
                             tpf, _update_darcy_t,
                             tpf->t_darcy);

  /* Diffusion property for the definition of the Darcy field in the gas
     phase */

  _add_pty_array(n_cells, cs_flag_primal_cell, tpf->diff_g_pty);

  /* Allocate and initialize arrays for the physical properties */
  /* ---------------------------------------------------------- */

  switch (tpf->approx_type) {

  case CS_GWF_TPF_APPROX_VERTEX_SUBCELL:
    {
    cs_xdef_t *d = nullptr;

    /* Relative permeability in the liquid and gas phase */

    d = _add_pty_array(c2v_size, cs_flag_dual_cell_byc, tpf->krl_pty);
    cs_xdef_array_set_adjacency(d, c2v);

    d = _add_pty_array(c2v_size, cs_flag_dual_cell_byc, tpf->krg_pty);
    cs_xdef_array_set_adjacency(d, c2v);

    /* Liquid saturation and liquid capacity */

    _add_pty_array(c2v_size, cs_flag_dual_cell_byc, tpf->lsat_pty);
    cs_xdef_array_set_adjacency(d, c2v);

    _add_pty_array(c2v_size, cs_flag_dual_cell_byc, tpf->lcap_pty);
    cs_xdef_array_set_adjacency(d, c2v);
    }
    break;

  default:

    /* Relative permeability in the liquid and gas phase */

    _add_pty_array(n_cells, cs_flag_primal_cell, tpf->krl_pty);
    _add_pty_array(n_cells, cs_flag_primal_cell, tpf->krg_pty);

    /* Liquid saturation and liquid capacity */

    _add_pty_array(n_cells, cs_flag_primal_cell, tpf->lsat_pty);
    _add_pty_array(n_cells, cs_flag_primal_cell, tpf->lcap_pty);
    break;

  }

  /* Properties associated to terms in equations */
  /* ------------------------------------------- */

  switch(tpf->solver_type) {

  case CS_GWF_TPF_SOLVER_PCPG_COUPLED:
    _finalize_setup_pcpg_coupled_solver(connect, tpf);
    break;

  case CS_GWF_TPF_SOLVER_PLPC_COUPLED:
    _finalize_setup_plpc_coupled_solver(connect, tpf);
    break;

  case CS_GWF_TPF_SOLVER_PLPG_SEGREGATED:
    _finalize_setup_plpg_segregated_solver(connect, tpf);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver type", __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of two-phase flows in a porous media
 *        (miscible or immiscible case)
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in, out] tpf         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init_values(const cs_cdo_connect_t        *connect,
                       const cs_cdo_quantities_t     *cdoq,
                       cs_gwf_tpf_t                  *tpf)
{
  CS_NO_WARN_IF_UNUSED(connect);

  if (tpf == nullptr)
    return;

  /* Set a builder hook function if needed */

  if (tpf->solver_type == CS_GWF_TPF_SOLVER_PLPG_SEGREGATED) {

    /* Add a hook function to define the source term. This operation should
       be done after that the corresponding equation has been initialized
       (a builder structure has to be defined) */

    cs_equation_add_build_hook(tpf->h_eq,
                               tpf,                   /* hook context */
                               _build_h_eq_diff_st); /* hook function */

  }

  /* Define and initialize the non-linear solver */

  if (tpf->use_incremental_solver) {

    /* This model is non-linear so that one needs an iterative algorithm to
       manage this non-linearity */

    if (tpf->nl_algo_type == CS_PARAM_NL_ALGO_NONE)
      tpf->nl_algo_type = CS_PARAM_NL_ALGO_PICARD;

  }

  if (tpf->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

    /* One assumes that the discretization schemes are all set to CDO
       vertex-based schemes */

    assert(tpf->h_eq->param->space_scheme == CS_SPACE_SCHEME_CDOVB);
    assert(tpf->w_eq->param->space_scheme == CS_SPACE_SCHEME_CDOVB);

    if (tpf->nl_algo_type == CS_PARAM_NL_ALGO_PICARD)
      tpf->nl_algo = cs_iter_algo_create_with_settings(CS_ITER_ALGO_DEFAULT,
                                                       tpf->nl_algo_verbosity,
                                                       tpf->nl_cvg_param);

    else if (tpf->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON) {

      tpf->nl_algo = cs_iter_algo_create_with_settings(CS_ITER_ALGO_ANDERSON,
                                                       tpf->nl_algo_verbosity,
                                                       tpf->nl_cvg_param);

      cs_lnum_t  size = cdoq->n_vertices;
      if (tpf->use_coupled_solver)
        size *= 2;

      cs_iter_algo_set_anderson_param(tpf->nl_algo,
                                      tpf->anderson_param,
                                      size);

    } /* Anderson acc. */

  } /* non-linear algo. ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new state for the groundwater flows module.
 *        Case of two-phase flows in porous media.
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step    pointer to a cs_time_step_t structure
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] tpf          pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_compute(const cs_mesh_t               *mesh,
                   const cs_cdo_connect_t        *connect,
                   const cs_cdo_quantities_t     *cdoq,
                   const cs_time_step_t          *time_step,
                   cs_flag_t                      option_flag,
                   cs_gwf_tpf_t                  *tpf)
{
  if (tpf == nullptr)
    return;

  switch(tpf->solver_type) {

  case CS_GWF_TPF_SOLVER_PCPG_COUPLED:
  case CS_GWF_TPF_SOLVER_PLPC_COUPLED:
    {
    cs_field_t *pa = nullptr, *pb = nullptr;

    if (tpf->solver_type == CS_GWF_TPF_SOLVER_PCPG_COUPLED)
      pa = tpf->c_pressure, pb = tpf->g_pressure;
    else
      pa = tpf->l_pressure, pb = tpf->c_pressure;

    switch (tpf->nl_algo_type) {

    case CS_PARAM_NL_ALGO_NONE:                    /* Linear case */
      cs_equation_system_solve(true, tpf->system); /* cur2prev = true */

      /* Update the variables related to the groundwater flow system */

      cs_gwf_tpf_update(mesh,
                        connect,
                        cdoq,
                        time_step,
                        CS_FLAG_CURRENT_TO_PREVIOUS,
                        option_flag,
                        tpf);
      break;

    case CS_PARAM_NL_ALGO_PICARD:
      _compute_coupled_picard(
        mesh, connect, cdoq, time_step, option_flag, pa, pb, tpf);
      break;

    case CS_PARAM_NL_ALGO_ANDERSON: {
      _compute_coupled_anderson(
        mesh, connect, cdoq, time_step, option_flag, pa, pb, tpf);
    } break;

      default:
        break; /* Nothing else to do */

      } /* Switch on the type of algorithm */
    }
    break;

  case CS_GWF_TPF_SOLVER_PLPG_SEGREGATED:
    _compute_segregated(mesh, time_step, connect, cdoq,
                        option_flag,
                        tpf);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver type", __func__);

  } /* Switch on the type of solver */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the update step in the case of a two-phase flow model in
 *        porous media (miscible or immiscible)
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      update_flag  metadata associated to type of operation to do
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] tpf           pointer to a TPF model context
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_tpf_update(const cs_mesh_t             *mesh,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *cdoq,
                  const cs_time_step_t        *ts,
                  cs_flag_t                    update_flag,
                  cs_flag_t                    option_flag,
                  cs_gwf_tpf_t                *tpf)
{
  CS_NO_WARN_IF_UNUSED(option_flag); /* For a later usage (gravity) */

  cs_real_t  time_eval = ts->t_cur;
  bool  cur2prev = false;

  if (update_flag & CS_FLAG_CURRENT_TO_PREVIOUS) {

    time_eval = cs_equation_get_time_eval(ts, tpf->w_eq);
    cur2prev = true;

  }

  if (cs_equation_get_space_scheme(tpf->w_eq) != CS_SPACE_SCHEME_CDOVB)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Incompatible space discretization.", __func__);

  /* Update pressures fields/arrays */

  _update_pressures(connect, cdoq, update_flag, cur2prev, tpf);

  /* Update properties related to soils:
   * - relative permeabilities: krl, krg
   * - liquid_saturation
   * - liquid capacity: \frac{\partial S_l}{\partial P_c}
   *
   * Either a call to a user-defined function or a predefined function if the
   * soil corresponds to a pre-defined model
   */

  cs_gwf_soil_update(time_eval, mesh, connect, cdoq);

  /* Define the liquid saturation in each cell when the liquid saturation has
     been defined on the submesh */

  _update_liquid_saturation_at_cells(connect, cdoq, tpf);

  /* Update arrays associated to each term of the equations */

  switch(tpf->solver_type) {

    /* ============================== */
  case CS_GWF_TPF_SOLVER_PCPG_COUPLED:
    /* ============================== */

    switch (cs_property_get_dim(tpf->diff_wg_pty)) {

    case 1: /* Isotropic case */
      /*       ++++++++++++++ */
      if (tpf->is_miscible)
        bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);
      //_update_iso_mtpf_pcpg_coupled(connect, cdoq, tpf);
      else
        _update_iso_itpf_pcpg_coupled_diffview(connect, cdoq, tpf);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Only the isotropic case is available up to now.",
                __func__);
      break;

    } /* Switch on the permeability type */
    break;

    /* ============================== */
  case CS_GWF_TPF_SOLVER_PLPC_COUPLED:
    /* ============================== */

    switch (cs_property_get_dim(tpf->diff_wl_pty)) {

    case 1: /* Isotropic case */
      /*       ++++++++++++++ */
      if (tpf->is_miscible)
        _update_iso_mtpf_plpc_coupled(connect, cdoq, tpf);
      else
        _update_iso_itpf_plpc_coupled(connect, cdoq, tpf);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Only the isotropic case is available up to now.",
                __func__);
      break;

    } /* Switch on the permeability type */
    break;

    /* ================================= */
  case CS_GWF_TPF_SOLVER_PLPG_SEGREGATED:
    /* ================================= */

    switch (cs_property_get_dim(tpf->diff_wl_pty)) {

    case 1: /* Isotropic case */
      /*       ++++++++++++++ */
      bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Only the isotropic case is available up to now.",
                __func__);
      break;

    } /* Switch on the permeability type */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver type", __func__);
    break;

  } /* Switch on solver type */

  /* Update the Darcy fluxes. This should be done at the end since one needs
     the new state for: the pressure fields, krl, krg, etc. */

  _update_darcy_fluxes(connect, cdoq, time_eval, cur2prev, tpf);

  return time_eval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined extra-operations for the groundwater flow module in case
 *        of miscible or immiscible two-phase flows in porous media
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts         pointer to a cs_time_step_t struct.
 * \param[in]      post_flag  requested quantities to be postprocessed
 * \param[in, out] tpf        pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_extra_op(const cs_cdo_connect_t       *connect,
                    const cs_cdo_quantities_t    *cdoq,
                    const cs_time_step_t         *ts,
                    cs_flag_t                     post_flag,
                    cs_gwf_tpf_t                 *tpf)
{
  assert(tpf != nullptr);

  if (cs_flag_test(post_flag, CS_GWF_POST_DARCY_FLUX_BALANCE)) {

    /* Balance for the Darcy advective flux in the liquid phase */

    cs_gwf_darcy_flux_balance(connect, cdoq,
                              cs_equation_get_param(tpf->w_eq),
                              tpf->l_darcy);

    /* Balance for the Darcy advective flux in the gas phase */

    cs_gwf_darcy_flux_balance(connect, cdoq,
                              cs_equation_get_param(tpf->h_eq),
                              tpf->g_darcy);

  }

  if (cs_flag_test(post_flag, CS_GWF_POST_SOIL_STATE))
    cs_gwf_soil_update_soil_state(cdoq->n_cells, tpf->l_saturation->val);

  if (cs_flag_test(post_flag, CS_GWF_POST_SOIL_MINMAX)) {

    const cs_adjacency_t  *c2v = connect->c2v;

    const cs_real_t  *pc = tpf->c_pressure->val;
    const cs_real_t  *pg = tpf->g_pressure->val;
    const cs_real_t  *pl = tpf->l_pressure->val;
    const cs_real_t  *lsat = tpf->l_saturation->val;
    const cs_real_t  *lcap = cs_property_get_array(tpf->lcap_pty);

    const int  n_soils = cs_gwf_get_n_soils();
    const int  n_min_outputs = n_soils * CS_GWF_TPF_N_OUTPUT_VARS;
    const int  n_outputs = 2 * n_min_outputs;

    cs_real_t *output_values = nullptr;

    BFT_MALLOC(output_values, n_outputs, cs_real_t);

    cs_real_t  *min_outputs = output_values;
    cs_real_t  *max_outputs = output_values + n_min_outputs;

    for (int i = 0; i < n_min_outputs; i++)
      min_outputs[i] =  DBL_MAX; /* min */
    for (int i = 0; i < n_min_outputs; i++)
      max_outputs[i] = -DBL_MAX; /* max */

    for (int id = 0; id < n_soils; id++) {

      cs_real_t  *_minv = min_outputs + CS_GWF_TPF_N_OUTPUT_VARS*id;
      cs_real_t  *_maxv = max_outputs + CS_GWF_TPF_N_OUTPUT_VARS*id;

      const cs_zone_t  *z = cs_gwf_soil_get_zone(id);
      assert(z != nullptr);

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t  c_id = z->elt_ids[i];
        const double  _lsat = lsat[c_id];

        _minv[3] = fmin(_minv[3], _lsat); /* SlMin */
        _maxv[3] = fmax(_maxv[3], _lsat); /* SlMax */

        if (tpf->approx_type != CS_GWF_TPF_APPROX_VERTEX_SUBCELL) {

          const double _lcap = lcap[c_id];
          _minv[4] = fmin(_minv[4], _lcap); /* ClMin */
          _maxv[4] = fmax(_maxv[4], _lcap); /* ClMax */

        }

        for (cs_lnum_t ii = c2v->idx[c_id]; ii < c2v->idx[c_id+1]; ii++) {

          const cs_lnum_t  v_id = c2v->ids[ii];

          _minv[0] = fmin(_minv[0], pg[v_id]); /* PgMin */
          _maxv[0] = fmax(_maxv[0], pg[v_id]); /* PgMax */

          _minv[1] = fmin(_minv[1], pl[v_id]); /* PlMin */
          _maxv[1] = fmax(_maxv[1], pl[v_id]); /* PlMax */

          _minv[2] = fmin(_minv[2], pc[v_id]); /* PcMin */
          _maxv[2] = fmax(_maxv[2], pc[v_id]); /* PcMax */

          if (tpf->approx_type == CS_GWF_TPF_APPROX_VERTEX_SUBCELL) {

            _minv[4] = fmin(_minv[4], lcap[ii]); /* ClMin */
            _maxv[4] = fmax(_maxv[4], lcap[ii]); /* ClMax */

          }

        } /* Loop on cell vertices */

      } /* Loop on zone cells */

    } /* Loop on soils */

    if (cs_glob_n_ranks > 1) {

      cs_parall_min(n_min_outputs, CS_REAL_TYPE, min_outputs);
      cs_parall_max(n_min_outputs, CS_REAL_TYPE, max_outputs);

    }

    if (cs_glob_rank_id < 1 && cs_gwf_tpf_time_plot != nullptr)
      cs_time_plot_vals_write(cs_gwf_tpf_time_plot,
                              ts->nt_cur,
                              ts->t_cur,
                              n_outputs,
                              output_values);

    BFT_FREE(output_values);

  } /* Min./Max. output */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined post-processing output for the groundwater flow module
 *        in case of saturated two-phase flows (tpf) in porous media.
 *
 * \param[in] mesh_id      id of the output mesh for the current call
 * \param[in] n_cells      local number of cells of post_mesh
 * \param[in] cell_ids     list of cells (0 to n-1)
 * \param[in] post_flag    flag gathering quantities to postprocess
 * \param[in] abs_perm     property for the absolute permeability
 * \param[in] tpf          pointer to the model context structure
 * \param[in] connect      pointer to additional connectivities for CDO
 * \param[in] cdoq         pointer to additional mesh quantities for CDO
 * \param[in] time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_extra_post(int                         mesh_id,
                      cs_lnum_t                   n_cells,
                      const cs_lnum_t             cell_ids[],
                      cs_flag_t                   post_flag,
                      const cs_property_t        *abs_perm,
                      const cs_gwf_tpf_t         *tpf,
                      const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *cdoq,
                      const cs_time_step_t       *time_step)
{
  if (mesh_id != CS_POST_MESH_VOLUME)
    return; /* Only postprocessings in the volume are defined */

  assert(tpf != nullptr);

  if (post_flag & CS_GWF_POST_SOIL_CAPACITY) {

    if (tpf->approx_type == CS_GWF_TPF_APPROX_VERTEX_SUBCELL) {

      cs_real_t *lcap_cell = nullptr;
      BFT_MALLOC(lcap_cell, n_cells, cs_real_t);

      cs_reco_scalar_vbyc2c(n_cells, cell_ids,
                            connect->c2v, cdoq,
                            cs_property_get_array(tpf->lcap_pty),
                            true, /* dense output */
                            lcap_cell);

      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_DEFAULT,
                        "l_capacity",
                        1,
                        false, /* interlace */
                        false, /* use parent */
                        CS_POST_TYPE_cs_real_t,
                        lcap_cell,
                        nullptr,
                        nullptr,
                        time_step);

      BFT_FREE(lcap_cell);

    }
    else { /* The array storing the soil capacity is alrady defined at cells */

      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_DEFAULT,
                        "l_capacity",
                        1,
                        false, /* interlace */
                        true,  /* use parent */
                        CS_POST_TYPE_cs_real_t,
                        cs_property_get_array(tpf->lcap_pty),
                        nullptr,
                        nullptr,
                        time_step);
    }

  } /* Postprocess the soil capacity */

  if (post_flag & CS_GWF_POST_PERMEABILITY) {

    /* permeability = krl * abs_permeability */

    cs_real_t *permeability = nullptr;
    int  dim = cs_property_get_dim(abs_perm);
    int  post_dim = (dim == 1) ? 1 : 9;

    BFT_MALLOC(permeability, post_dim*n_cells, cs_real_t);

    if (dim > 1) {

      if (tpf->approx_type == CS_GWF_TPF_APPROX_VERTEX_SUBCELL) {

        const cs_real_t  *krl_c2v = cs_property_get_array(tpf->krl_pty);
        const cs_adjacency_t  *c2v = connect->c2v;

        for (cs_lnum_t i = 0; i < n_cells; i++) {

          cs_lnum_t  c_id = (cell_ids == nullptr || n_cells == cdoq->n_cells)
                              ? i
                              : cell_ids[i];
          cs_real_t  tensor[3][3];

          cs_property_get_cell_tensor(c_id,
                                      time_step->t_cur,
                                      abs_perm,
                                      false, /* inversion */
                                      tensor);

          double  krl_cell = 0;
          for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
            krl_cell = cdoq->pvol_vc[j] * krl_c2v[j];
          krl_cell /= cdoq->cell_vol[c_id];

          cs_real_t  *_cell_perm = permeability + post_dim*i;
          for (int ki = 0; ki < 3; ki++)
            for (int kj = 0; kj < 3; kj++)
              _cell_perm[3*ki+kj] = krl_cell * tensor[ki][kj];

        } /* Loop on selected cells */

      }
      else {

        const cs_real_t  *krl_cell = cs_property_get_array(tpf->krl_pty);

        for (cs_lnum_t i = 0; i < n_cells; i++) {

          cs_lnum_t  c_id = (cell_ids == nullptr || n_cells == cdoq->n_cells)
                              ? i
                              : cell_ids[i];
          cs_real_t  tensor[3][3];

          cs_property_get_cell_tensor(c_id,
                                      time_step->t_cur,
                                      abs_perm,
                                      false, /* inversion */
                                      tensor);

          cs_real_t  *_cell_perm = permeability + post_dim*i;
          for (int ki = 0; ki < 3; ki++)
            for (int kj = 0; kj < 3; kj++)
              _cell_perm[3*ki+kj] = krl_cell[c_id] * tensor[ki][kj];

        } /* Loop on selected cells */
      }

    }
    else { /* dim = 1 */

      /* One uses permeability to store the evaluation of the relative
         permeabilty if needed */

      if (tpf->approx_type == CS_GWF_TPF_APPROX_VERTEX_SUBCELL) {

        cs_reco_scalar_vbyc2c(n_cells, cell_ids,
                              connect->c2v, cdoq,
                              cs_property_get_array(tpf->krl_pty),
                              true, /* dense output */
                              permeability);

        if (cs_property_is_uniform(abs_perm)) {

          const double  abs_perm_value =
            cs_property_get_cell_value(0, time_step->t_cur, abs_perm);

          for (cs_lnum_t i = 0; i < n_cells; i++)
            permeability[i] *= abs_perm_value;

        }
        else {

          for (cs_lnum_t i = 0; i < n_cells; i++) {

            cs_lnum_t c_id = (cell_ids == nullptr || n_cells == cdoq->n_cells)
                               ? i
                               : cell_ids[i];

            double  abs_perm_cell = cs_property_get_cell_value(c_id,
                                                               time_step->t_cur,
                                                               abs_perm);
            permeability[i] *= abs_perm_cell;

          } /* Loop on selected cells */

        } /* abs_perm is uniform ? */

      }
      else {

        const cs_real_t  *krl_values = cs_property_get_array(tpf->krl_pty);

        if (cs_property_is_uniform(abs_perm)) {

          const double  abs_perm_value =
            cs_property_get_cell_value(0, time_step->t_cur, abs_perm);

          for (cs_lnum_t i = 0; i < n_cells; i++) {

            cs_lnum_t c_id = (cell_ids == nullptr || n_cells == cdoq->n_cells)
                               ? i
                               : cell_ids[i];

            permeability[i] = abs_perm_value * krl_values[c_id];

          }

        }
        else {

          for (cs_lnum_t i = 0; i < n_cells; i++) {

            cs_lnum_t c_id = (cell_ids == nullptr || n_cells == cdoq->n_cells)
                               ? i
                               : cell_ids[i];

            double  abs_perm_cell = cs_property_get_cell_value(c_id,
                                                               time_step->t_cur,
                                                               abs_perm);
            permeability[i] = abs_perm_cell * krl_values[c_id];

          }

        } /* abs_perm is uniform ? */

      } /* type of approximation */

    } /* post_dim */

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_DEFAULT,
                      "permeability",
                      post_dim,
                      true,  /* interlace */
                      false, /* use_parent */
                      CS_POST_TYPE_cs_real_t,
                      permeability,
                      nullptr,
                      nullptr,
                      time_step);

    BFT_FREE(permeability);

  } /* Post-processing of the permeability field */

  if (post_flag & CS_GWF_POST_GAS_MASS_DENSITY) {

    cs_real_t *gas_mass_density = nullptr;
    BFT_MALLOC(gas_mass_density, n_cells, cs_real_t);

    cs_reco_scalar_v2c(n_cells, cell_ids,
                       connect->c2v, cdoq,
                       tpf->g_pressure->val, /* at vertices */
                       true, /* dense output */
                       gas_mass_density);

    const double  mh_ov_rt =
      tpf->h_molar_mass / (tpf->ref_temperature * cs_physical_constants_r);

    for (cs_lnum_t c = 0; c < n_cells; c++)
      gas_mass_density[c] *= mh_ov_rt;

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_DEFAULT,
                      "gas_mass_density",
                      1,
                      false, /* interlace */
                      false, /* use parent */
                      CS_POST_TYPE_cs_real_t,
                      gas_mass_density,
                      nullptr,
                      nullptr,
                      time_step);

    BFT_FREE(gas_mass_density);

  } /* Post-processing of the gas mass density */

  if (post_flag & CS_GWF_POST_SOIL_STATE) {

    const int  *soil_state = cs_gwf_soil_get_soil_state();

    if (soil_state != nullptr)
      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_DEFAULT,
                        "soil_state",
                        1,
                        false, /* interlace */
                        true,  /* use_parent */
                        CS_POST_TYPE_int,
                        soil_state,
                        nullptr,
                        nullptr,
                        time_step);

  } /* Post-processing of the soil state */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
