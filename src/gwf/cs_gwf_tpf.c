/*============================================================================
 * Main functions dedicated to the modelling of two-phase flows in a porous
 * media. This media is always considered as unsaturated. Two sub-models are
 * considered: miscible (MTPF) or immiscible (ITPF)
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>

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

/*============================================================================
 * Local definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the modelling of miscible two-phase flows
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_mtpf_log_setup(cs_gwf_tpf_t   *mc)
{
  if (mc == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Water mass density: %5.3e, viscosity: %5.3e,"
                " molar mass: %5.3e\n", mc->l_mass_density, mc->l_viscosity,
                mc->w_molar_mass);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Henry constant: %5.3e\n", mc->henry_constant);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Gas component: viscosity: %5.3e, diffusivity"
                " in the liquid phase: %5.3e, molar mass: %5.3e\n",
                mc->g_viscosity, mc->l_diffusivity_h, mc->h_molar_mass);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the modelling of immiscible two-phase flows
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_itpf_log_setup(cs_gwf_tpf_t   *mc)
{
  if (mc == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Water mass density: %5.3e, viscosity: %5.3e\n",
                mc->l_mass_density, mc->l_viscosity);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Henry constant: %5.3e (Should be very low)\n",
                mc->henry_constant);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Gas component viscosity: %5.3e, molar mass: %5.3e\n",
                mc->g_viscosity, mc->h_molar_mass);
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

  cs_gwf_tpf_t  *mc = context;

  /* Modify the property data associated to the hodge operator related to the
     diffusion term */

  cs_property_data_t  *saved = diff_hodge->pty_data;
  cs_property_data_t  tmp = cs_property_data_define(saved->need_tensor,
                                                    saved->need_eigen,
                                                    mc->diff_hl_pty);

  diff_hodge->pty_data = &tmp;

  /* Diffusion part of the source term to add */

  cs_hodge_evaluate_property_cw(cm, cb->t_pty_eval, cb->cell_flag, diff_hodge);

  /* Define the local stiffness matrix: local matrix owned by the cellwise
     builder (store in cb->loc) */

  eqc->get_stiffness_matrix(cm, diff_hodge, cb);

  /* Retrieve the values pl^{k,n+1} */

  double  *vec = cb->values, *matvec = cb->values + cm->n_vc;
  for (int v = 0; v < cm->n_vc; v++)
    vec[v] = mc->l_pressure->val[cm->v_ids[v]];

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
 *        cell_values could be set to NULL when the space discretization does
 *        not request these values for the update.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      dof_values    values to consider for the update
 * \param[in]      cell_values   values to consider for the update or NULL
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

  assert(darcy != NULL);
  assert(darcy->flux_val != NULL);
  assert(darcy->update_input != NULL);

  cs_adv_field_t  *adv = darcy->adv_field;
  assert(adv != NULL);
  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t  *mc = (cs_gwf_tpf_t *)darcy->update_input;

  /* Update the array of flux values associated to the advection field.
   *
   * diff_pty for the w_eq --> rho_l * abs_perm * krl / mu_l Thus, one needs to
   * divide by rho_l for the Darcy flux. Thanks to the rescaling, there is no
   * need to define a new property which is very close to the diff_wl_pty
   */

  cs_property_set_scaling_factor(mc->diff_wl_pty, 1./mc->l_mass_density);

  cs_equation_compute_diffusive_flux(mc->w_eq,
                                     NULL, /* eqp --> default */
                                     NULL, /* diff_pty --> default */
                                     dof_values,
                                     cell_values,
                                     darcy->flux_location,
                                     t_eval,
                                     darcy->flux_val);

  cs_property_unscale(mc->diff_wl_pty);

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != NULL);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);

  /* Update the Darcy flux at the boundary (take into account the BCs) since
     the liquid pressure is the unknown */

  cs_gwf_darcy_flux_update_on_boundary(t_eval, mc->w_eq, adv);

  cs_field_t  *bdy_nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  if (bdy_nflx != NULL) { /* Values of the Darcy flux at boundary face exist */

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
 *        cell_values could be set to NULL when the space discretization does
 *        not request these values for the update.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      dof_values    values to consider for the update
 * \param[in]      cell_values   values to consider for the update or NULL
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

  assert(darcy != NULL);
  assert(darcy->flux_val != NULL);
  assert(darcy->update_input != NULL);
  assert(darcy->adv_field != NULL);

  cs_adv_field_t  *adv = darcy->adv_field;

  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t  *mc = darcy->update_input;
  cs_equation_t  *eq = mc->h_eq;

  /* Update the array of flux values associated to the advection field */

  cs_equation_compute_diffusive_flux(eq,
                                     NULL, /* eqp --> default */
                                     mc->diff_g_pty,
                                     dof_values,
                                     cell_values,
                                     darcy->flux_location,
                                     t_eval,
                                     darcy->flux_val);

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != NULL);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);

  /* Update the Darcy flux at boundary faces if needed */

  cs_gwf_darcy_flux_update_on_boundary_wo_eq(connect, cdoq, vel->val, adv);

  cs_field_t  *bdy_nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  if (bdy_nflx != NULL) { /* Values of the Darcy flux at boundary face exist */

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
 *        cell_values could be set to NULL when the space discretization does
 *        not request these values for the update.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      dof_values    values to consider for the update
 * \param[in]      cell_values   values to consider for the update or NULL
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

  assert(darcy != NULL);
  assert(darcy->flux_val != NULL);
  assert(darcy->update_input != NULL);
  assert(darcy->adv_field != NULL);

  cs_adv_field_t  *adv = darcy->adv_field;

  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t  *mc = darcy->update_input;
  cs_equation_t  *eq = mc->h_eq;

  /* Update the array of flux values associated to the advection field */

  cs_equation_compute_diffusive_flux(eq,
                                     NULL, /* eqp --> default */
                                     NULL, /* diff_pty --> default */
                                     dof_values,
                                     cell_values,
                                     darcy->flux_location,
                                     t_eval,
                                     darcy->flux_val);

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != NULL);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);

  /* Update the Darcy flux at boundary faces if needed */

  cs_field_t  *bdy_nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  if (bdy_nflx != NULL) { /* Values of the Darcy flux at boundary face exist */

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
 *        cell_values could be set to NULL when the space discretization does
 *        not request these values for the update.
 *
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]      dof_values    values to consider for the update
 * \param[in]      cell_values   values to consider for the update or NULL
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

  assert(darcy != NULL);
  assert(darcy->flux_val != NULL);
  assert(darcy->update_input != NULL);
  assert(darcy->adv_field != NULL);

  cs_adv_field_t  *adv = darcy->adv_field;

  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t  *mc = darcy->update_input;
  const cs_gwf_darcy_flux_t  *l_darcy = mc->l_darcy;
  const cs_gwf_darcy_flux_t  *g_darcy = mc->g_darcy;

  const double  hmh = mc->h_molar_mass * mc->henry_constant;
  const double  mh_ov_rt =
    mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);
  const cs_adjacency_t  *c2e = connect->c2e;

# pragma omp parallel for if (6*cdoq->n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < c2e->idx[cdoq->n_cells]; i++)
    darcy->flux_val[i] =
      hmh * l_darcy->flux_val[i] + mh_ov_rt * g_darcy->flux_val[i];

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != NULL);
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
 * \param[in, out] mc           pointer to a TPF model context
 */
/*----------------------------------------------------------------------------*/

static void
_update_darcy_fluxes(const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *cdoq,
                     cs_real_t                    t_eval,
                     bool                         cur2prev,
                     cs_gwf_tpf_t                *mc)
{
  cs_real_t  *dof_vals = NULL, *cell_vals = NULL;

  /* Darcy velocity/flux in the liquid phase. The liquid pressure is the
     unknown for the coupled and the segregated solver */

  cs_gwf_darcy_flux_t  *l_darcy = mc->l_darcy;

  cs_gwf_get_value_pointers(mc->w_eq, &dof_vals, &cell_vals);

  l_darcy->update_func(connect,
                       cdoq,
                       dof_vals,
                       cell_vals,
                       t_eval,
                       cur2prev,
                       l_darcy);

  /* Darcy velocity/flux in the gas phase. */

  cs_gwf_darcy_flux_t  *g_darcy = mc->g_darcy;

  g_darcy->update_func(connect,
                       cdoq,
                       mc->g_pressure->val,
                       mc->g_pressure_cells,
                       t_eval,
                       cur2prev,
                       g_darcy);

  if (!mc->use_diffusion_view_for_darcy) {

    /* Only useful for the computation of the advective term in the
       conservation equation for the hydrogen */

    cs_gwf_darcy_flux_t  *t_darcy = mc->t_darcy;

    t_darcy->update_func(connect,
                         cdoq,
                         NULL,    /* values are defined inside the call */
                         NULL,    /* values are defined inside the call */
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
 * \param[in, out] mc           pointer to a TPF model context
 */
/*----------------------------------------------------------------------------*/

static void
_update_pressures(const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *cdoq,
                  cs_flag_t                    update_flag,
                  bool                         cur2prev,
                  cs_gwf_tpf_t                *mc)
{
  const cs_lnum_t  n_vertices = cdoq->n_vertices;

  /* Compute either the capillarity pressure or the gas pressure accordind to
     the choice of solver. Pressure variables which are solved have already
     been updated */

  const cs_real_t  *l_pr = mc->l_pressure->val;

  if (mc->use_coupled_solver) {

    /* Main pressure variables: liquid and capillarity pressures */

    const cs_real_t  *c_pr = mc->c_pressure->val;
    cs_real_t  *g_pr = mc->g_pressure->val;

    if (cur2prev && !cs_flag_test(update_flag, CS_FLAG_INITIALIZATION))
      cs_field_current_to_previous(mc->g_pressure);

    /* Compute the new values of the gas pressure at vertices */

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_vertices; i++)
      g_pr[i] = l_pr[i] + c_pr[i];

    /* Avoid to add an unsteady contribution at the first iteration  */

    if (update_flag & CS_FLAG_INITIALIZATION) {
      cs_array_real_copy(n_vertices, c_pr, mc->c_pressure->val_pre);
      cs_array_real_copy(n_vertices, l_pr, mc->l_pressure->val_pre);
      cs_field_current_to_previous(mc->g_pressure);
    }

  }
  else { /* Segregated solver */

    /* Main pressure variables: liquid and gas pressures */

    const cs_real_t  *g_pr = mc->g_pressure->val;
    cs_real_t  *c_pr = mc->c_pressure->val;

    if (cur2prev && !cs_flag_test(update_flag, CS_FLAG_INITIALIZATION))
      cs_field_current_to_previous(mc->c_pressure);

    /* Compute the values of the capillarity pressure at vertices */

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_vertices; i++)
      c_pr[i] = g_pr[i] - l_pr[i];

    /* Avoid to add an unsteady contribution at the first iteration  */

    if (update_flag & CS_FLAG_INITIALIZATION) {
      cs_array_real_copy(n_vertices, c_pr, mc->c_pressure->val_pre);
      cs_array_real_copy(n_vertices, g_pr, mc->g_pressure->val_pre);
      cs_field_current_to_previous(mc->c_pressure);
    }

  }

  /* Interpolate cell values of the capillarity pressure from the vertex
     ones. The capillarity pressure at cells is used to update quantities
     related to a soil model */

  if (mc->c_pressure_cells != NULL)
    cs_reco_pv_at_cell_centers(connect->c2v, cdoq,
                               mc->c_pressure->val,
                               mc->c_pressure_cells);

  /* Interpolate cell values of the gas pressure from the vertex ones */

  if (mc->g_pressure_cells != NULL)
    cs_reco_pv_at_cell_centers(connect->c2v, cdoq,
                               mc->g_pressure->val,
                               mc->g_pressure_cells);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update several arrays associated to the definition of terms involved
 *        in the resolution of an immiscible two-phase flow model with an
 *        isotropic absolute permeability.
 *
 *        Numerical options: coupled solver and diffusive viewpoint of the
 *        Darcy terms in the conservation equation of the hydrogen
 *
 * \param[in, out] mc    pointer to the model context to update
 */
/*----------------------------------------------------------------------------*/

static void
_update_iso_itpf_coupled_diffview_terms(cs_gwf_tpf_t     *mc)
{
  const double  hmh = mc->h_molar_mass * mc->henry_constant;
  const double  mh_ov_rt =
    mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);

  const cs_real_t  *pg_cells = mc->g_pressure_cells;
  const cs_real_t  *l_sat = mc->l_saturation->val;
  const cs_real_t  *l_cap = mc->l_capacity;
  const cs_real_t  *krl =  mc->l_rel_permeability;
  const cs_real_t  *krg =  mc->g_rel_permeability;

  /* Loop on soils */

  for (int soil_id = 0; soil_id < cs_gwf_get_n_soils(); soil_id++) {

    cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);

    assert(soil != NULL);
    assert(soil->hydraulic_model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE);
    assert(soil->abs_permeability_dim == 1);

    const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);
    assert(zone != NULL);

    const double  k_abs = soil->abs_permeability[0][0];
    const double  phi = soil->porosity;
    const double  phi_rhol = phi * mc->l_mass_density;

    const double  l_diff_coef = k_abs/mc->l_viscosity;
    const double  wl_diff_coef = mc->l_mass_density * l_diff_coef;
    const double  hg_diff_coef = k_abs/mc->g_viscosity;

    /* Loop on cells belonging to this soil */

    for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

      const cs_lnum_t  c_id = zone->elt_ids[i];

      const double  rhog_h = mh_ov_rt * pg_cells[c_id];
      const double  rhol_h = hmh * pg_cells[c_id];
      const double  sl = l_sat[c_id], sg = 1 - sl;
      const double  dsl_dpc = l_cap[c_id];

      /* Update terms for the Darcy flux in the gaz phase */

      mc->diff_g_array[c_id] = hg_diff_coef * krg[c_id];

      /* Update terms associated to the water conservation equation */

      mc->time_wc_array[c_id] = phi_rhol * dsl_dpc;
      mc->diff_wl_array[c_id] = wl_diff_coef * krl[c_id];

      /* Update terms associated to the hydrogen conservation equation */

      const double  time_h_coef = phi * ( mh_ov_rt*sg + hmh*sl );
      mc->time_hc_array[c_id] = time_h_coef + phi * (rhol_h - rhog_h)*dsl_dpc;
      mc->time_hl_array[c_id] = time_h_coef;

      const double  diff_hc = rhog_h * mc->diff_g_array[c_id];
      mc->diff_hc_array[c_id] = diff_hc;
      mc->diff_hl_array[c_id] = diff_hc + rhol_h * l_diff_coef * krl[c_id];

#if 0 /* Debugging purpose */
      cs_log_printf(CS_LOG_DEFAULT,
                    "%s: %02d |> di_g: % 8.6e; ti_wc: %8.6e; di_wl: %8.6e;"
                    " ti_hc: % 8.6e; ti_hl: % 8.6e; di_hc: % 8.6e;"
                    " di_hl: % 8.6e\n", __func__, c_id,
                    mc->diff_g_array[c_id], mc->time_wc_array[c_id],
                    mc->diff_wl_array[c_id], mc->time_hc_array[c_id],
                    mc->time_hl_array[c_id], mc->diff_hc_array[c_id],
                    mc->diff_hl_array[c_id]);
#endif
    } /* Loop on cells of the zone (= soil) */

  } /* Loop on soils */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the non-linear algorithm.
 *
 * \param[in, out] mc           pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_non_linear_algo(cs_gwf_tpf_t                  *mc)
{
  if (mc->nl_algo_type == CS_PARAM_NL_ALGO_NONE)
    return;

  cs_iter_algo_t  *algo = mc->nl_algo;
  assert(algo != NULL);

  cs_iter_algo_reset_nl(mc->nl_algo_type, algo);

  /* No resolution has been done at this stage and also no previous to current
     operation.

     The current values will be the previous state after the current to
     previous operation
  */

  const cs_real_t  *pg_0 = mc->g_pressure->val;
  const cs_real_t  *pl_0 = mc->l_pressure->val;

  /* 1. Set the normalization factor. Performed with Pg and Pl even Pc is used
     in the resolution. One prefers to use Pg since it also appears in the
     different terms for the h_eq. */

  algo->normalization = cs_cdo_blas_square_norm_pvsp(pg_0);
  algo->normalization += cs_cdo_blas_square_norm_pvsp(pl_0);
  algo->normalization = sqrt(algo->normalization);
  if (algo->normalization < cs_math_zero_threshold)
    algo->normalization = 1.0;

  cs_log_printf(CS_LOG_DEFAULT, "%s: Non-linear algo. normalization=%6.4e\n",
                __func__, algo->normalization);
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
 * \param[in]      nl_algo_type type of non-linear algorithm
 * \param[in]      dpc_iter     cur. increment values for the capill. pressure
 * \param[in, out] algo         pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_check_cvg_nl_dpc(cs_param_nl_algo_t        nl_algo_type,
                      const cs_real_t          *dpc_iter,
                      cs_iter_algo_t           *algo)
{
  assert(algo != NULL);
  algo->prev_res = algo->res;

  double dpc_norm = cs_cdo_blas_square_norm_pvsp(dpc_iter);

  algo->res = dpc_norm;
  assert(algo->res > -DBL_MIN);
  algo->res = sqrt(algo->res);

  if (algo->n_algo_iter < 1) /* Store the first residual to detect a
                                possible divergence of the algorithm */
    algo->res0 = algo->res;

  /* Update the convergence members */

  cs_iter_algo_update_cvg_default(algo);

  if (algo->verbosity > 0) {

    if (algo->n_algo_iter == 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    "### GWF.TPF %10s.It    Algo.Res   Tolerance\n",
                    cs_param_get_nl_algo_label(nl_algo_type));

    cs_log_printf(CS_LOG_DEFAULT,
                  "### GWF.TPF %10s.It%02d  %5.3e  %6.4e\n",
                  cs_param_get_nl_algo_label(nl_algo_type),
                  algo->n_algo_iter, algo->res, algo->tol);

  } /* verbosity > 0 */

  return algo->cvg_status;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the convergence of the non-linear algorithm in case of TPF
 *        model
 *
 * \param[in]      nl_algo_type type of non-linear algorithm
 * \param[in]      pc_pre_iter  previous iterate values for the gas pressure
 * \param[in, out] pc_cur_iter  current iterate values for the gas pressure
 * \param[in]      pl_pre_iter  previous iterate values for the liquid pressure
 * \param[in, out] pl_cur_iter  current iterate values for the liquid pressure
 * \param[in, out] algo         pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_check_cvg_nl(cs_param_nl_algo_t        nl_algo_type,
              const cs_real_t          *pc_pre_iter,
              cs_real_t                *pc_cur_iter,
              const cs_real_t          *pl_pre_iter,
              cs_real_t                *pl_cur_iter,
              cs_iter_algo_t           *algo)
{
  if (nl_algo_type == CS_PARAM_NL_ALGO_NONE)
    return CS_SLES_CONVERGED;

  assert(algo != NULL);

  if (nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON && algo->n_algo_iter > 0) {

    /* pg_* arrays gather pg and pl (this is done during the solve step) */

    cs_iter_algo_aa_update(algo,
                           pc_cur_iter, /* updated during the process */
                           pc_pre_iter,
                           cs_cdo_blas_dotprod_2pvsp,
                           cs_cdo_blas_square_norm_2pvsp);

  } /* Anderson acceleration */

  algo->prev_res = algo->res;

  double delta_pc = cs_cdo_blas_square_norm_pvsp_diff(pc_pre_iter, pc_cur_iter);
  double delta_pl = cs_cdo_blas_square_norm_pvsp_diff(pl_pre_iter, pl_cur_iter);

  algo->res = delta_pc + delta_pl;
  assert(algo->res > -DBL_MIN);
  algo->res = sqrt(algo->res);

  if (algo->n_algo_iter < 1) { /* Store the first residual to detect a
                                  possible divergence of the algorithm */
    algo->res0 = algo->res;
    algo->prev_res = algo->res;

  }

  /* Update the convergence members */

  cs_iter_algo_update_cvg_default(algo);

  if (algo->verbosity > 0) {

    if (algo->n_algo_iter == 1) {

      if (algo->verbosity > 1)
        cs_log_printf(CS_LOG_DEFAULT,
                      "### GWF.TPF %10s.It    Algo.Res   Tolerance"
                      "  ||D_Pc||  ||D_Pl||\n",
                      cs_param_get_nl_algo_label(nl_algo_type));
      else
        cs_log_printf(CS_LOG_DEFAULT,
                      "### GWF.TPF %10s.It    Algo.Res   Tolerance\n",
                      cs_param_get_nl_algo_label(nl_algo_type));

    }

    if (algo->verbosity > 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    "### GWF.TPF %10s.It%02d  %5.3e  %6.4e %5.3e %5.3e\n",
                    cs_param_get_nl_algo_label(nl_algo_type),
                    algo->n_algo_iter, algo->res, algo->tol,
                    sqrt(delta_pc), sqrt(delta_pl));
    else
      cs_log_printf(CS_LOG_DEFAULT,
                    "### GWF.TPF %10s.It%02d  %5.3e  %6.4e\n",
                    cs_param_get_nl_algo_label(nl_algo_type),
                    algo->n_algo_iter, algo->res, algo->tol);

  } /* verbosity > 0 */

  return algo->cvg_status;
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
 * \param[in, out] mc           pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_compute_coupled_picard(const cs_mesh_t              *mesh,
                        const cs_cdo_connect_t       *connect,
                        const cs_cdo_quantities_t    *cdoq,
                        const cs_time_step_t         *time_step,
                        cs_flag_t                     option_flag,
                        cs_gwf_tpf_t                 *mc)
{
  bool cur2prev = true;
  cs_flag_t  update_flag = CS_FLAG_CURRENT_TO_PREVIOUS;

  cs_iter_algo_t  *algo = mc->nl_algo;
  assert(algo != NULL);
  assert(mc->nl_algo_type == CS_PARAM_NL_ALGO_PICARD);

  _init_non_linear_algo(mc);

  cs_real_t  *pc_kp1 = NULL, *pl_kp1 = NULL;  /* at ^{n+1,k+1} */
  cs_real_t  *pc_k = NULL, *pl_k = NULL;      /* at ^{n+1,k} */

  BFT_MALLOC(pc_k, cdoq->n_vertices, cs_real_t);
  BFT_MALLOC(pl_k, cdoq->n_vertices, cs_real_t);

  pc_kp1 = mc->c_pressure->val; /* val stores always the latest values */
  pl_kp1 = mc->l_pressure->val;

  do {

    /* current values: n+1,k */

    cs_array_real_copy(cdoq->n_vertices, mc->c_pressure->val, pc_k);
    cs_array_real_copy(cdoq->n_vertices, mc->l_pressure->val, pl_k);

    cs_equation_system_solve(cur2prev, mc->system);

    /* current values: n+1,k+1 now */

    if (algo->n_algo_iter > 0 && mc->nl_relax_factor < 1.0) {

      const double  relax = mc->nl_relax_factor;
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
        pc_kp1[i] = pc_kp1[i]*relax + pc_k[i]*(1-relax);
        pl_kp1[i] = pl_kp1[i]*relax + pl_k[i]*(1-relax);
      }

    }

    /* Update the variables related to the groundwater flow system */

    cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                      update_flag, option_flag,
                      mc);

    /* After the first resolution, no cur2prev anymore */

    cur2prev = false;
    update_flag = 0;

  } while (_check_cvg_nl(mc->nl_algo_type,
                         pc_k, pc_kp1, pl_k, pl_kp1,
                         algo) == CS_SLES_ITERATING);

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_post_check(__func__,
                          mc->system->param->name,
                          cs_param_get_nl_algo_label(mc->nl_algo_type),
                          algo);

  /* Free temporary arrays and structures */

  BFT_FREE(pc_k);
  BFT_FREE(pl_k);
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
 * \param[in, out] mc           pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_compute_coupled_anderson(const cs_mesh_t              *mesh,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *cdoq,
                          const cs_time_step_t         *time_step,
                          cs_flag_t                     option_flag,
                          cs_gwf_tpf_t                 *mc)
{
  bool  cur2prev = false;
  cs_flag_t  update_flag = 0;   /* No current to previous operation */
  cs_iter_algo_t  *algo = mc->nl_algo;
  assert(algo != NULL);

  cs_iter_algo_reset_nl(mc->nl_algo_type, algo);

  /* A first resolution has been done followed by an update */

  cs_real_t  *pg_kp1 = NULL, *pl_kp1 = NULL;  /* at ^{n+1,k+1} */
  cs_real_t  *pg_k = NULL, *pl_k = NULL;      /* at ^{n+1,k} */

  BFT_MALLOC(pg_k, 2*cdoq->n_vertices, cs_real_t);
  pl_k = pg_k + cdoq->n_vertices;

  cs_array_real_copy(cdoq->n_vertices, mc->g_pressure->val_pre, pg_k);
  cs_array_real_copy(cdoq->n_vertices, mc->l_pressure->val_pre, pl_k);

  /* One needs only one array gathering the liquid and gas pressures */

  BFT_MALLOC(pg_kp1, 2*cdoq->n_vertices, cs_real_t);
  pl_kp1 = pg_kp1 + cdoq->n_vertices;

  cs_array_real_copy(cdoq->n_vertices, mc->g_pressure->val, pg_kp1);
  cs_array_real_copy(cdoq->n_vertices, mc->l_pressure->val, pl_kp1);

  /* Set the normalization factor */

  algo->normalization = cs_cdo_blas_square_norm_pvsp(pg_k);
  algo->normalization += cs_cdo_blas_square_norm_pvsp(pl_k);
  algo->normalization = sqrt(algo->normalization);
  if (algo->normalization < cs_math_zero_threshold)
    algo->normalization = 1.0;

  /* Main non-linear loop */

  while (_check_cvg_nl(mc->nl_algo_type,
                       pg_k, pg_kp1, pl_k, pl_kp1,
                       algo) == CS_SLES_ITERATING) {

    cs_array_real_copy(cdoq->n_vertices, pg_kp1, pg_k);
    cs_array_real_copy(cdoq->n_vertices, pl_kp1, pl_k);

    /* Update the variables related to the groundwater flow system.
     * In case of an Anderson acceleration, pg_kp1 and pl_kp1 may be
     * updated */

    if (algo->n_algo_iter >= mc->anderson_param.starting_iter) {

      cs_array_real_copy(cdoq->n_vertices, pg_kp1, mc->g_pressure->val);
      cs_array_real_copy(cdoq->n_vertices, pl_kp1, mc->l_pressure->val);
    }

    cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                      update_flag,
                      option_flag,
                      mc);

    /* Build and solve the linear system related to the coupled system of
       equations. First call: current --> previous and then no operation */

    cs_equation_system_solve(cur2prev, mc->system);

    cs_array_real_copy(cdoq->n_vertices, mc->g_pressure->val, pg_kp1);
    cs_array_real_copy(cdoq->n_vertices, mc->l_pressure->val, pl_kp1);

  } /* Until convergence */

  /* Get the soil state up to date with the last compute values for the
     liquid and gas pressures */

  cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                    update_flag,
                    option_flag,
                    mc);

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_post_check(__func__,
                          mc->system->param->name,
                          cs_param_get_nl_algo_label(mc->nl_algo_type),
                          algo);

  /* Free temporary arrays and structures */

  BFT_FREE(pg_k);
  BFT_FREE(pg_kp1);
  cs_iter_algo_aa_free_arrays(algo->context);
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
 * \param[in, out] mc           pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_compute_segregated(const cs_mesh_t              *mesh,
                    const cs_time_step_t         *time_step,
                    const cs_cdo_connect_t       *connect,
                    const cs_cdo_quantities_t    *cdoq,
                    cs_flag_t                     option_flag,
                    cs_gwf_tpf_t                 *mc)
{
  assert(mc->use_incremental_solver);

  /* Copy the state for time t^n into field->val_pre */

  cs_field_current_to_previous(mc->g_pressure);
  cs_field_current_to_previous(mc->l_pressure);

  /* Since the cur2prev operation has been done, avoid to do it again */

  bool cur2prev = false;        /* Done just above */
  cs_flag_t  update_flag = 0;   /* No current to previous operation */

  /* Initialize the non-linear algorithm */

  cs_iter_algo_t  *algo = mc->nl_algo;
  assert(algo != NULL);

  cs_iter_algo_reset_nl(mc->nl_algo_type, algo);

  algo->normalization = cs_cdo_blas_square_norm_pvsp(mc->c_pressure->val);
  if (algo->normalization < cs_math_zero_threshold)
    algo->normalization = 1.0;
  else
    algo->normalization = sqrt(algo->normalization);

  /* First solve step:
   * 1. Solve the equation associated to Pl --> Pl^(n+1,1) knowing
   *    Pl^(n+1,0) = Pl^n, Pg^n and thus Pc^n and Sl^n
   *    --> Sl^(n+1,0) - Sl^n = 0 --> no source term
   * 2. Solve the equation associated to Pg --> Pg^(n+1,1) knowing
   *
   */

  cs_equation_solve(cur2prev, mesh, mc->w_eq);

  cs_equation_solve(cur2prev, mesh, mc->h_eq);

  /* Update the variables related to the groundwater flow system */

  cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                    CS_FLAG_CURRENT_TO_PREVIOUS, /* Force this operation */
                    option_flag,
                    mc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
  char  label[32];

  sprintf(label, "Pl_iter%02d", algo->n_algo_iter);
  cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                           CS_POST_WRITER_DEFAULT,
                           label,
                           1,
                           false,
                           false,
                           CS_POST_TYPE_cs_real_t,
                           cs_field_by_name("liquid_pressure")->val,
                           time_step);

  sprintf(label, "Pg_iter%02d", algo->n_algo_iter);
  cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                           CS_POST_WRITER_DEFAULT,
                           label,
                           1,
                           false,
                           false,
                           CS_POST_TYPE_cs_real_t,
                           cs_field_by_name("gas_pressure")->val,
                           time_step);

  sprintf(label, "Pc_iter%02d", algo->n_algo_iter);
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

  cs_real_t  *dpc_kp1 = NULL;
  BFT_MALLOC(dpc_kp1, mesh->n_vertices, cs_real_t);

  _get_capillarity_pressure_increment(mesh, mc->w_eq, mc->h_eq, dpc_kp1);

  while(_check_cvg_nl_dpc(mc->nl_algo_type,
                          dpc_kp1, algo) == CS_SLES_ITERATING) {

    /* Solve step */

    cs_equation_solve(cur2prev, mesh, mc->w_eq);

    cs_equation_solve(cur2prev, mesh, mc->h_eq);

    /* Update the variables related to the groundwater flow system */

    cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                      update_flag,
                      option_flag,
                      mc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
    sprintf(label, "Pl_iter%02d", algo->n_algo_iter);
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_DEFAULT,
                             label,
                             1,
                             false,
                             false,
                             CS_POST_TYPE_cs_real_t,
                             cs_field_by_name("liquid_pressure")->val,
                             time_step);

    sprintf(label, "Pg_iter%02d", algo->n_algo_iter);
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_DEFAULT,
                             label,
                             1,
                             false,
                             false,
                             CS_POST_TYPE_cs_real_t,
                             cs_field_by_name("gas_pressure")->val,
                             time_step);

    sprintf(label, "Pc_iter%02d", algo->n_algo_iter);
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

    _get_capillarity_pressure_increment(mesh, mc->w_eq, mc->h_eq, dpc_kp1);

  } /* while not converged */

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_post_check(__func__,
                          "Segregated incremental TPF solver",
                          cs_param_get_nl_algo_label(mc->nl_algo_type),
                          algo);

  BFT_FREE(dpc_kp1);
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

  cs_gwf_tpf_t  *mc = NULL;

  BFT_MALLOC(mc, 1, cs_gwf_tpf_t);

  /* Create a new equation associated to the mass conservation of water. The
     unknown is the capillarity pressure. This will stand for the (0,0)-block
     if a coupled system is considered. */

  mc->w_eq = cs_equation_add("w_conservation",   /* equation name */
                             "liquid_pressure",  /* variable name */
                             CS_EQUATION_TYPE_GROUNDWATER,
                             1,
                             CS_PARAM_BC_HMG_NEUMANN);

  cs_equation_param_t  *w_eqp = cs_equation_get_param(mc->w_eq);

  cs_equation_param_set(w_eqp, CS_EQKEY_BC_ENFORCEMENT, "algebraic");
  cs_equation_param_set(w_eqp, CS_EQKEY_HODGE_TIME_ALGO, "voronoi");
  cs_equation_param_set(w_eqp, CS_EQKEY_HODGE_REAC_ALGO, "voronoi");

  /* Set to NULL since one has to know if a coupled or segregated solver is
     used */

  mc->h_eq = NULL;
  mc->b01_w_eqp = NULL;
  mc->b10_h_eqp = NULL;
  mc->system = NULL;

  /* Advection fields */
  /* ---------------- */

  mc->l_darcy = NULL;
  mc->g_darcy = NULL;
  mc->t_darcy = NULL;

  /* Properties
   * ----------
   *
   * Adding the properties related to the permeability (diffusion) is
   * postponed since one has to know which type of permeability is considered
   * (iso, ortho, or anisotropic) and this information is a result of the type
   * of soils which have been added.
   */

  /* Properties related to the water equation  */

  mc->time_wc_pty = NULL;
  mc->diff_wl_pty = NULL;

  /* Properties related to the hydrogen equation  */

  mc->time_hc_pty = NULL;
  mc->diff_hc_pty = NULL;
  mc->time_hl_pty = NULL;
  mc->diff_hl_pty = NULL;
  mc->reac_h_pty = NULL;

  /* Fields */
  /* ------ */

  mc->c_pressure = NULL;
  mc->l_pressure = NULL;
  mc->g_pressure = NULL;
  mc->l_saturation = NULL;

  /* Arrays */
  /* ------ */

  /* The properties will be defined using arrays.
   * Store these arrays of property values associated to equation terms */

  mc->time_wc_array = NULL;
  mc->diff_wl_array = NULL;

  mc->time_hc_array = NULL;
  mc->diff_hc_array = NULL;
  mc->time_hl_array = NULL;
  mc->diff_hl_array = NULL;

  mc->srct_w_array = NULL;
  mc->srct_h_array = NULL;
  mc->reac_h_array = NULL;

  /* Array of additional variable or property values */

  mc->l_rel_permeability = NULL;
  mc->g_rel_permeability = NULL;
  mc->c_pressure_cells = NULL;
  mc->g_pressure_cells = NULL;
  mc->l_capacity = NULL;
  mc->l_saturation_submesh = NULL;

  /* Model parameters (default values) */
  /* ---------------- */

  if (model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE) {

    mc->is_miscible = false;
    mc->l_diffusivity_h = 0;      /* immiscible case */
    mc->henry_constant = 1e-20;   /* nearly immiscible case */

  }
  else {

    mc->is_miscible = false;
    mc->l_diffusivity_h = 0;      /* immiscible case */
    mc->henry_constant = 1e-20;   /* nearly immiscible case */

  }

  mc->l_mass_density = 1000;
  mc->l_viscosity = 1e-3;
  mc->g_viscosity = 2e-5;
  mc->w_molar_mass = 18e-3;
  mc->h_molar_mass = 3e-3;
  mc->ref_temperature = 280;    /* in Kelvin */

  /* Numerical parameters (default values) */

  mc->use_coupled_solver = true;
  mc->use_incremental_solver = false;
  mc->use_definition_on_submesh = false;
  mc->use_diffusion_view_for_darcy = true;

  mc->nl_algo_type = CS_PARAM_NL_ALGO_PICARD;

  mc->nl_relax_factor = 1.0;
  mc->nl_algo_cvg.n_max_iter = 50;
  mc->nl_algo_cvg.rtol = 1e-5;
  mc->nl_algo_cvg.atol = 1e-10;
  mc->nl_algo_cvg.dtol = 1e3;

  mc->anderson_param.n_max_dir = 5;
  mc->anderson_param.starting_iter = 3;
  mc->anderson_param.max_cond = -1; /* No test by default */
  mc->anderson_param.beta = 1.0;    /* No damping by default */
  mc->anderson_param.dp_type = CS_PARAM_DOTPROD_EUCLIDEAN;

  mc->nl_algo = cs_iter_algo_create(0, mc->nl_algo_cvg);

  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to the modelling of two-phase
 *        flows in a porous media
 *
 * \param[in, out] p_mc   pointer of pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_free(cs_gwf_tpf_t    **p_mc)
{
  if (p_mc == NULL)
    return;
  if (*p_mc == NULL)
    return;

  cs_gwf_tpf_t  *mc = *p_mc;

  /* System of equations are freed elsewhere (just after having freed
     cs_equation_t structures) */

  cs_gwf_darcy_flux_free(&(mc->l_darcy));
  cs_gwf_darcy_flux_free(&(mc->g_darcy));
  cs_gwf_darcy_flux_free(&(mc->t_darcy));

  BFT_FREE(mc->time_wc_array);
  BFT_FREE(mc->diff_wl_array);
  BFT_FREE(mc->time_hc_array);
  BFT_FREE(mc->diff_hc_array);
  BFT_FREE(mc->time_hl_array);
  BFT_FREE(mc->diff_hl_array);

  BFT_FREE(mc->srct_w_array);
  BFT_FREE(mc->srct_h_array);
  BFT_FREE(mc->reac_h_array);

  BFT_FREE(mc->diff_g_array);

  BFT_FREE(mc->l_rel_permeability);
  BFT_FREE(mc->g_rel_permeability);
  BFT_FREE(mc->c_pressure_cells);
  BFT_FREE(mc->g_pressure_cells);
  BFT_FREE(mc->l_capacity);
  BFT_FREE(mc->l_saturation_submesh);

  /* Non-linearity: If the context is not NULL, this means that an Anderson
     algorithm has been activated otherwise nothing to do */

  cs_iter_algo_aa_free(mc->nl_algo);

  BFT_FREE(mc->nl_algo);

  BFT_FREE(mc);
  *p_mc = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the model context of two-phase flows.
 *        Common to the different sub-models relying on two-phase flows.
 *
 * \param[in] mc      pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_log_setup(cs_gwf_tpf_t          *mc)
{
  if (mc == NULL)
    return;

  if (mc->is_miscible)
    _mtpf_log_setup(mc);
  else
    _itpf_log_setup(mc);

  cs_log_printf(CS_LOG_SETUP, "  * GWF | Reference temperature: %5.2f K\n",
                mc->ref_temperature);

  cs_gwf_darcy_flux_log(mc->l_darcy);
  cs_gwf_darcy_flux_log(mc->g_darcy);

  if (!mc->use_diffusion_view_for_darcy)
    cs_gwf_darcy_flux_log(mc->t_darcy);

  if (mc->use_coupled_solver)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Coupled solver\n");
  else
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Segregated solver\n");

  if (mc->use_incremental_solver)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Incremental solver\n");

  if (mc->use_definition_on_submesh)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Definition on submesh activated\n");

  if (mc->use_diffusion_view_for_darcy)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Diffusion view for Darcy terms\n");

  if (mc->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

    cs_log_printf(CS_LOG_SETUP, "  * GWF | Non-linear algo.: %s\n",
                  cs_param_get_nl_algo_name(mc->nl_algo_type));

    cs_log_printf(CS_LOG_SETUP, "  * GWF | Tolerances of non-linear algo:"
                  " rtol: %5.3e; atol: %5.3e; dtol: %5.3e; max_iter: %d\n",
                  mc->nl_algo_cvg.rtol, mc->nl_algo_cvg.atol,
                  mc->nl_algo_cvg.dtol, mc->nl_algo_cvg.n_max_iter);

    if (mc->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON) {

      const cs_iter_algo_param_aa_t  aap = mc->anderson_param;

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
 * \param[in, out] mc          pointer to the model context structure
 * \param[in, out] perm_type   type of permeability to handle
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init(cs_gwf_tpf_t            *mc,
                cs_property_type_t       perm_type)
{
  if (mc == NULL)
    return;

  if (mc->w_eq == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Water equation is emty.", __func__);

  /* Property common to all solvers */

  mc->diff_wl_pty = cs_property_add("diff_wl_pty", perm_type);
  mc->diff_g_pty = cs_property_add("diff_g_pty", perm_type);

  /* Darcy flux (one assumes a CDOVB space discretization) */

  cs_advection_field_status_t  adv_status =
    CS_ADVECTION_FIELD_GWF | CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;

  cs_adv_field_t  *l_adv_field = cs_advection_field_add("l_darcy_field",
                                                        adv_status);
  mc->l_darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);
  mc->l_darcy->adv_field = l_adv_field;

  cs_adv_field_t  *g_adv_field = cs_advection_field_add("g_darcy_field",
                                                        adv_status);
  mc->g_darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);
  mc->g_darcy->adv_field = g_adv_field;
  if (!mc->use_diffusion_view_for_darcy) {

    cs_adv_field_t  *t_adv_field = cs_advection_field_add("t_darcy_field",
                                                          adv_status);
    mc->t_darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);
    mc->t_darcy->adv_field = t_adv_field;

  }

  if (mc->use_coupled_solver) {

    /* Define the system of equations */
    /* ------------------------------ */

    /* Create a new equation associated to the mass conservation of gas
       (hydrogen for instance). The unknown is the capillarity pressure. This
       will stand for the (1,1)-block if a coupled system is considered. */

    mc->h_eq = cs_equation_add("h_conservation",       /* equation name */
                               "capillarity_pressure", /* variable name */
                               CS_EQUATION_TYPE_GROUNDWATER,
                               1,
                               CS_PARAM_BC_HMG_NEUMANN);

    cs_equation_param_t  *b00_w_eqp = cs_equation_get_param(mc->w_eq);
    cs_equation_param_t  *b11_h_eqp = cs_equation_get_param(mc->h_eq);

    cs_equation_param_set(b11_h_eqp, CS_EQKEY_BC_ENFORCEMENT, "algebraic");
    cs_equation_param_set(b11_h_eqp, CS_EQKEY_HODGE_TIME_ALGO, "voronoi");
    cs_equation_param_set(b11_h_eqp, CS_EQKEY_HODGE_REAC_ALGO, "voronoi");

    /* Define the coupled system of equations */
    /* -------------------------------------- */

    /* Create the (0,1)-block related to the water in the gas phase */

    mc->b01_w_eqp = cs_equation_param_create("block01_w_eq",
                                             CS_EQUATION_TYPE_GROUNDWATER,
                                             1,
                                             CS_PARAM_BC_HMG_NEUMANN);

    cs_equation_param_set(mc->b01_w_eqp, CS_EQKEY_BC_ENFORCEMENT, "algebraic");
    cs_equation_param_set(mc->b01_w_eqp, CS_EQKEY_HODGE_TIME_ALGO, "voronoi");
    cs_equation_param_set(mc->b01_w_eqp, CS_EQKEY_HODGE_REAC_ALGO, "voronoi");

    /* Create the (1,0)-block related to the hydrogen in the liquid phase */

    mc->b10_h_eqp = cs_equation_param_create("block10_h_eq",
                                             CS_EQUATION_TYPE_GROUNDWATER,
                                             1,
                                             CS_PARAM_BC_HMG_NEUMANN);

    cs_equation_param_set(mc->b10_h_eqp, CS_EQKEY_BC_ENFORCEMENT, "algebraic");
    cs_equation_param_set(mc->b10_h_eqp, CS_EQKEY_HODGE_TIME_ALGO, "voronoi");
    cs_equation_param_set(mc->b10_h_eqp, CS_EQKEY_HODGE_REAC_ALGO, "voronoi");

    /* Add a 2x2 system of coupled equations and define each block */

    mc->system = cs_equation_system_add("two_phase_flow_porous_media",
                                        2,   /* system size */
                                        1);  /* scalar-valued block */

    /* Set all the blocks in the coupled system */

    cs_equation_system_assign_equation(0, mc->w_eq, mc->system);  /* (0,0) */
    cs_equation_system_assign_equation(1, mc->h_eq, mc->system);  /* (1,1) */
    cs_equation_system_assign_param(0, 1, mc->b01_w_eqp, mc->system);
    cs_equation_system_assign_param(1, 0, mc->b10_h_eqp, mc->system);

    /* Properties
     * ----------
     *
     * - diffusion term for water eq. in the liquid phase     (0,0)-block
     * - unsteady term for water eq. in the liquid phase      (0,1) block
     * - unsteady term for hydrogen eq. in the gaseous phase  (1,1)-block
     * - unsteady term for hydrogen eq. in the liquid phase   (1,0)-block
     */

    if (mc->use_definition_on_submesh) {

      mc->time_wc_pty = cs_property_subcell_add("time_wc_pty", CS_PROPERTY_ISO);
      mc->time_hc_pty = cs_property_subcell_add("time_hc_pty", CS_PROPERTY_ISO);
      mc->time_hl_pty = cs_property_subcell_add("time_hl_pty", CS_PROPERTY_ISO);

    }
    else {

      mc->time_wc_pty = cs_property_add("time_wc_pty", CS_PROPERTY_ISO);
      mc->time_hc_pty = cs_property_add("time_hc_pty", CS_PROPERTY_ISO);
      mc->time_hl_pty = cs_property_add("time_hl_pty", CS_PROPERTY_ISO);

    }

    /* Add terms to the water equation */
    /* ------------------------------- */

    cs_equation_add_time(mc->b01_w_eqp, mc->time_wc_pty);

    cs_equation_add_diffusion(b00_w_eqp, mc->diff_wl_pty);

    /* Add terms to the hydrogen equation */
    /* ---------------------------------- */

    cs_equation_add_time(mc->b10_h_eqp, mc->time_hl_pty);
    cs_equation_add_time(b11_h_eqp, mc->time_hc_pty);

    if (mc->use_diffusion_view_for_darcy) {

      mc->diff_hl_pty = cs_property_add("diff_hl_pty", perm_type);
      mc->diff_hc_pty = cs_property_add("diff_hc_pty", perm_type);

      cs_equation_add_diffusion(mc->b10_h_eqp, mc->diff_hl_pty);
      cs_equation_add_diffusion(b11_h_eqp, mc->diff_hc_pty);

    }
    else { /* Advection viewpoint */

      cs_equation_add_advection(mc->b10_h_eqp, mc->t_darcy->adv_field);
      cs_equation_add_advection(b11_h_eqp, mc->t_darcy->adv_field);

    }

    if (mc->use_incremental_solver) {

      b00_w_eqp->incremental_algo_type = CS_PARAM_NL_ALGO_PICARD;
      b11_h_eqp->incremental_algo_type = CS_PARAM_NL_ALGO_PICARD;

    }

  }
  else { /* Segregated solver */

    /* Create a new equation associated to the mass conservation of gas
       (hydrogen for instance). The unknown is the capillarity pressure. This
       will stand for the (1,1)-block if a coupled system is considered. */

    mc->h_eq = cs_equation_add("h_conservation",  /* equation name */
                               "gas_pressure",    /* variable name */
                               CS_EQUATION_TYPE_GROUNDWATER,
                               1,
                               CS_PARAM_BC_HMG_NEUMANN);

    cs_equation_param_t  *w_eqp = cs_equation_get_param(mc->w_eq);
    cs_equation_param_t  *h_eqp = cs_equation_get_param(mc->h_eq);

    /* Segregated solver are always solved by increment and relies on an
     * advection viewpoint of the Darcy terms */

    mc->use_incremental_solver = true;
    mc->use_diffusion_view_for_darcy = false;

    /* Properties */
    /* ---------- */

    if (mc->use_definition_on_submesh)
      mc->time_hc_pty = cs_property_subcell_add("time_hc_pty", CS_PROPERTY_ISO);
    else
      mc->time_hc_pty = cs_property_add("time_hc_pty", CS_PROPERTY_ISO);

    /* Add terms to the water equation */
    /* ------------------------------- */

    cs_equation_add_diffusion(w_eqp, mc->diff_wl_pty);

    /* Add terms to the hydrogen equation */
    /* ---------------------------------- */

    cs_equation_add_time(h_eqp, mc->time_hc_pty);
    cs_equation_add_reaction(h_eqp, mc->reac_h_pty);

    if (mc->use_diffusion_view_for_darcy) {

      mc->diff_hc_pty = cs_property_add("diff_hc_pty", perm_type);

      cs_equation_add_diffusion(h_eqp, mc->diff_hc_pty);

    }
    else
      cs_equation_add_advection(h_eqp, mc->t_darcy->adv_field);

    if (mc->use_incremental_solver) {

      w_eqp->incremental_algo_type = CS_PARAM_NL_ALGO_PICARD;
      h_eqp->incremental_algo_type = CS_PARAM_NL_ALGO_PICARD;

    }

  } /* Segregated or coupled solver */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial setup stage for two-phase flows in porous media. At this
 *        stage, all soils have been defined and equation parameters are set.
 *        Case of a miscible or immiscible model.
 *
 * \param[in]      post_flag   optional postprocessing request(s)
 * \param[in, out] mc          pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init_setup(cs_flag_t         post_flag,
                      cs_gwf_tpf_t     *mc)
{
  if (mc == NULL)
    return;

  if (mc->w_eq == NULL || mc->h_eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Equations are not defined for this model. Stop execution.\n",
              __func__);

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");

  /* Retrieve the pointers to the structure storing the equation parameters */

  cs_equation_param_t  *w_eqp = cs_equation_get_param(mc->w_eq);
  cs_equation_param_t  *h_eqp = cs_equation_get_param(mc->h_eq);

  assert(w_eqp->space_scheme == h_eqp->space_scheme);

  int loc_id = c_loc_id;
  if (w_eqp->space_scheme == CS_SPACE_SCHEME_CDOVB ||
      w_eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB)
    loc_id = v_loc_id;
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme", __func__);

  /* Add the variable fields (Keep always the previous state) */

  cs_equation_predefined_create_field(1, mc->w_eq);
  cs_equation_predefined_create_field(1, mc->h_eq);

  /* Set fields related to variables */

  mc->l_pressure = cs_equation_get_field(mc->w_eq);

  if (mc->use_coupled_solver) {

    mc->c_pressure = cs_equation_get_field(mc->h_eq);

    /* One has to be consistent with the location of DoFs for the w_eq and h_eq
     * which are respectively related to the l_pressure and c_pressure */

    mc->g_pressure = cs_field_create("gas_pressure",
                                     field_mask,
                                     loc_id,
                                     1,
                                     true); /* has_previous */

    cs_field_set_key_int(mc->g_pressure, log_key, 1);
    cs_field_set_key_int(mc->g_pressure, post_key, 1);

  }
  else {

    mc->g_pressure = cs_equation_get_field(mc->h_eq);

    /* One has to be consistent with the location of DoFs for the w_eq and h_eq
     * which are respectively related to the l_pressure and g_pressure */

    mc->c_pressure = cs_field_create("capillarity_pressure",
                                     field_mask,
                                     loc_id,
                                     1,
                                     true); /* has_previous */

    cs_field_set_key_int(mc->c_pressure, log_key, 1);
    cs_field_set_key_int(mc->c_pressure, post_key, 1);

  }

  /* Create a liquid saturation field attached to cells: S_l (only keep the
     latest computed state) */

  int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;

  mc->l_saturation = cs_field_create("liquid_saturation",
                                     pty_mask,
                                     c_loc_id,
                                     1,     /* dimension */
                                     false); /* has_previous */

  cs_field_set_key_int(mc->l_saturation, log_key, 1);
  if (post_flag & CS_GWF_POST_LIQUID_SATURATION)
    cs_field_set_key_int(mc->l_saturation, post_key, 1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of two-phase flows in a porous media
 *        (miscible or immiscible case)
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]      flag      optional settings for the module
 * \param[in, out] mc        pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_finalize_setup(const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *cdoq,
                          cs_flag_t                    flag,
                          cs_gwf_tpf_t                *mc)
{
  CS_NO_WARN_IF_UNUSED(flag);   /* will be useful for gravity effect */

  const cs_adjacency_t  *c2v = connect->c2v;
  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  c2v_size = c2v->idx[n_cells];

  cs_equation_param_t  *w_eqp = cs_equation_get_param(mc->w_eq);
  cs_equation_param_t  *h_eqp = cs_equation_get_param(mc->h_eq);

  assert(cs_equation_get_type(mc->w_eq) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(cs_equation_get_type(mc->h_eq) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(w_eqp->space_scheme == h_eqp->space_scheme);

  if (w_eqp->space_scheme != CS_SPACE_SCHEME_CDOVB)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme", __func__);

  /* Set the Darcian flux (in the volume and at the boundary) */
  /* -------------------------------------------------------- */

  cs_gwf_darcy_flux_define(connect, cdoq, w_eqp->space_scheme,
                           mc, _update_darcy_l,
                           mc->l_darcy);

  if (mc->use_coupled_solver)
    cs_gwf_darcy_flux_define(connect, cdoq, h_eqp->space_scheme,
                             mc, _update_darcy_g_coupled,
                             mc->g_darcy);
  else
    cs_gwf_darcy_flux_define(connect, cdoq, h_eqp->space_scheme,
                             mc, _update_darcy_g_segregated,
                             mc->g_darcy);

  if (!mc->use_diffusion_view_for_darcy)
    cs_gwf_darcy_flux_define(connect, cdoq, w_eqp->space_scheme,
                             mc, _update_darcy_t,
                             mc->t_darcy);

  /* Allocate and initialize arrays for the physical properties */
  /* ---------------------------------------------------------- */

  /* Relative permeability in the liquid and gas phase.
   * One assumes that the medium is saturated by default */

  BFT_MALLOC(mc->l_rel_permeability, n_cells, cs_real_t);
  cs_array_real_set_scalar(n_cells, 1., mc->l_rel_permeability);

  BFT_MALLOC(mc->g_rel_permeability, n_cells, cs_real_t);
  cs_array_real_set_scalar(n_cells, 1., mc->g_rel_permeability);

  /* Interpolation of the gas pressure at cell centers is used to define the
   * properties associated to some terms in hydrogen eq. */

  BFT_MALLOC(mc->g_pressure_cells, n_cells, cs_real_t);
  cs_array_real_fill_zero(n_cells, mc->g_pressure_cells);

  BFT_MALLOC(mc->c_pressure_cells, n_cells, cs_real_t);
  cs_array_real_fill_zero(n_cells, mc->c_pressure_cells);

  /* Handle the possibility to define some properties on a sub-mesh */

  cs_lnum_t  _size = n_cells;
  cs_flag_t  _loc = cs_flag_primal_cell;

  if (mc->use_definition_on_submesh) {

    _size = c2v_size;
    _loc = cs_flag_dual_cell_byc;

    BFT_MALLOC(mc->l_saturation_submesh, _size, cs_real_t);
    cs_array_real_fill_zero(_size, mc->l_saturation_submesh);

  }

  BFT_MALLOC(mc->l_capacity, _size, cs_real_t);
  cs_array_real_fill_zero(_size, mc->l_capacity);

  /* Properties associated to terms in equations */
  /* ------------------------------------------- */

  /* Diffusion property for the water equation */

  BFT_MALLOC(mc->diff_wl_array, n_cells, cs_real_t);
  cs_array_real_fill_zero(n_cells, mc->diff_wl_array);

  cs_property_def_by_array(mc->diff_wl_pty,
                           NULL,                 /* all cells */
                           cs_flag_primal_cell,  /* data location */
                           mc->diff_wl_array,
                           false,                /* xdef not owner */
                           true);                /* full length */

  /* Diffusion property for the definition of the Darcy field in the gas
     phase */

  BFT_MALLOC(mc->diff_g_array, n_cells, cs_real_t);
  cs_array_real_fill_zero(n_cells, mc->diff_g_array);

  cs_property_def_by_array(mc->diff_g_pty,
                           NULL,                 /* all cells */
                           cs_flag_primal_cell,  /* data location */
                           mc->diff_g_array,
                           false,                /* xdef not owner */
                           true);                /* full length */

  if (mc->use_diffusion_view_for_darcy) {

    /* Diffusion property in the hydrogen equation (associated to Pc in the
       case of a coupled solver or to Pg in the case of a segregated solver) */

    BFT_MALLOC(mc->diff_hc_array, n_cells, cs_real_t);
    cs_array_real_fill_zero(n_cells, mc->diff_hc_array);

    cs_property_def_by_array(mc->diff_hc_pty,
                             NULL,                 /* all cells */
                             cs_flag_primal_cell,  /* data location */
                             mc->diff_hc_array,
                             false,                /* xdef not owner */
                             true);                /* full length */

    /* Diffusion property in the hydrogen equation associated to Pl. This is
       used to define a source term in the case of a segregated solver. */

    BFT_MALLOC(mc->diff_hl_array, n_cells, cs_real_t);
    cs_array_real_fill_zero(n_cells, mc->diff_hl_array);

    cs_property_def_by_array(mc->diff_hl_pty,
                             NULL,                /* all cells */
                             cs_flag_primal_cell, /* data location */
                             mc->diff_hl_array,
                             false,               /* xdef not owner */
                             true);               /* full length */

  }

  /* Settings now depending on the choice of the solver */

  if (mc->use_coupled_solver) {

    /* Define the arrays storing the time property for all blocks */

    BFT_MALLOC(mc->time_wc_array, _size, cs_real_t);
    cs_array_real_fill_zero(_size, mc->time_wc_array);

    cs_property_def_by_array(mc->time_wc_pty,
                             NULL,     /* all cells */
                             _loc,     /* data location */
                             mc->time_wc_array,
                             false,    /* xdef not owner */
                             true);    /* full length */

    BFT_MALLOC(mc->time_hl_array, _size, cs_real_t);
    cs_array_real_fill_zero(_size, mc->time_hl_array);

    cs_property_def_by_array(mc->time_hl_pty,
                             NULL,     /* all cells */
                             _loc,     /* data location */
                             mc->time_hl_array,
                             false,    /* xdef not owner */
                             true);    /* full length */

    BFT_MALLOC(mc->time_hc_array, _size, cs_real_t);
    cs_array_real_fill_zero(_size, mc->time_hc_array);

    cs_property_def_by_array(mc->time_hc_pty,
                             NULL,     /* all cells */
                             _loc,     /* data location */
                             mc->time_hc_array,
                             false,    /* xdef not owner */
                             true);    /* full length */

  }
  else { /* Segregated solver */

    /* Treatment of the water mass conservation */
    /* ---------------------------------------- */

    /* Source term */

    BFT_MALLOC(mc->srct_w_array, _size, cs_real_t);
    cs_array_real_fill_zero(_size, mc->srct_w_array);

    cs_xdef_t  *wst_def =
      cs_equation_add_source_term_by_array(w_eqp,
                                           NULL,    /* all cells */
                                           _loc,    /* data location */
                                           mc->srct_w_array,
                                           false,   /* xdef not owner */
                                           true);   /* full length */

    if (mc->use_definition_on_submesh)
      cs_xdef_array_set_adjacency(wst_def, c2v);

    /* Treatment of the hydrogen mass conservation */
    /* ------------------------------------------- */

    /* Source term (second one; the first is associated to a stiffness term and
       relies on a builder hook function) */

    BFT_MALLOC(mc->srct_h_array, _size, cs_real_t);
    cs_array_real_fill_zero(_size, mc->srct_h_array);

    cs_xdef_t  *hst_def =
      cs_equation_add_source_term_by_array(h_eqp,
                                           NULL,   /* all cells */
                                           _loc,   /* data location */
                                           mc->srct_h_array,
                                           false,  /* xdef not owner */
                                           true);  /* full length */

    if (mc->use_definition_on_submesh)
      cs_xdef_array_set_adjacency(hst_def, c2v);

    /* Unsteady term */

    BFT_MALLOC(mc->time_hc_array, _size, cs_real_t);
    cs_array_real_fill_zero(_size, mc->time_hc_array);

    cs_xdef_t  *th_pty_def =
      cs_property_def_by_array(mc->time_hc_pty,
                               NULL,        /* all cells */
                               _loc,        /* data location */
                               mc->time_hc_array,
                               false,       /* xdef not owner */
                               true);       /* full length */

    if (mc->use_definition_on_submesh)
      cs_xdef_array_set_adjacency(th_pty_def, c2v);

    /* Reaction term */

    BFT_MALLOC(mc->reac_h_array, _size, cs_real_t);
    cs_array_real_fill_zero(_size, mc->reac_h_array);

    cs_xdef_t  *rh_pty_def =
      cs_property_def_by_array(mc->reac_h_pty,
                               NULL,     /* all cells */
                               _loc,     /* data location */
                               mc->reac_h_array,
                               false,    /* xdef not owner */
                               true);    /* full length */

    if (mc->use_definition_on_submesh)
      cs_xdef_array_set_adjacency(rh_pty_def, c2v);

  } /* Segregated solver */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of two-phase flows in a porous media
 *        (miscible or immiscible case)
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in, out] mc          pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init_values(const cs_cdo_connect_t        *connect,
                       const cs_cdo_quantities_t     *cdoq,
                       cs_gwf_tpf_t                  *mc)
{
  CS_NO_WARN_IF_UNUSED(connect);

  if (mc == NULL)
    return;

  /* Set a builder hook function if needed */

  if (!mc->use_coupled_solver) {

    /* Add a hook function to define the source term. This operation should
       be done after that the corresponding equation has been initialized
       (a builder structure has to be defined) */

    cs_equation_add_build_hook(mc->h_eq,
                               mc,                   /* hook context */
                               _build_h_eq_diff_st); /* hook function */

  }

  /* Initialize the non-linear algorithm if needed */

  if (mc->use_incremental_solver) {

    /* This model is non-linear so that one needs an iterative algorithm to
       manage this non-linearity */

    if (mc->nl_algo_type == CS_PARAM_NL_ALGO_NONE)
      mc->nl_algo_type = CS_PARAM_NL_ALGO_PICARD;

  }

  if (mc->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

    /* Apply the convergence settings */

    mc->nl_algo->cvg_param.atol = mc->nl_algo_cvg.atol;
    mc->nl_algo->cvg_param.rtol = mc->nl_algo_cvg.rtol;
    mc->nl_algo->cvg_param.dtol = mc->nl_algo_cvg.dtol;
    mc->nl_algo->cvg_param.n_max_iter = mc->nl_algo_cvg.n_max_iter;

    /* One assumes that the discretization schemes are all set to CDO
       vertex-based schemes */

    assert(mc->h_eq->param->space_scheme == CS_SPACE_SCHEME_CDOVB);
    assert(mc->w_eq->param->space_scheme == CS_SPACE_SCHEME_CDOVB);

    cs_lnum_t  size = cdoq->n_vertices;
    if (mc->use_coupled_solver)
      size *= 2;

    if (mc->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
      mc->nl_algo->context = cs_iter_algo_aa_create(mc->anderson_param, size);

  }
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
 * \param[in, out] mc           pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_compute(const cs_mesh_t               *mesh,
                   const cs_cdo_connect_t        *connect,
                   const cs_cdo_quantities_t     *cdoq,
                   const cs_time_step_t          *time_step,
                   cs_flag_t                      option_flag,
                   cs_gwf_tpf_t                  *mc)
{
  if (mc == NULL)
    return;

  if (mc->use_coupled_solver) {

    switch (mc->nl_algo_type) {

    case CS_PARAM_NL_ALGO_NONE: /* Linear case */
      cs_equation_system_solve(true, mc->system); /* cur2prev = true */

      /* Update the variables related to the groundwater flow system */

      cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                        CS_FLAG_CURRENT_TO_PREVIOUS,
                        option_flag,
                        mc);
      break;

    case CS_PARAM_NL_ALGO_PICARD:
      _compute_coupled_picard(mesh, connect, cdoq, time_step,
                              option_flag,
                              mc);
      break;

    case CS_PARAM_NL_ALGO_ANDERSON:
      _compute_coupled_anderson(mesh, connect, cdoq, time_step,
                                option_flag,
                                mc);
      break;

    default:
      break; /* Nothing else to do */
    }

  }
  else
    _compute_segregated(mesh, time_step, connect, cdoq,
                        option_flag,
                        mc);
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
 * \param[in, out] mc           pointer to a TPF model context
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_tpf_update(const cs_mesh_t             *mesh,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *cdoq,
                  const cs_time_step_t        *ts,
                  cs_flag_t                    update_flag,
                  cs_flag_t                    option_flag,
                  cs_gwf_tpf_t                *mc)
{
  CS_NO_WARN_IF_UNUSED(option_flag); /* For a later usage (gravity) */

  cs_real_t  time_eval = ts->t_cur;
  bool  cur2prev = false;

  if (update_flag & CS_FLAG_CURRENT_TO_PREVIOUS) {

    time_eval = cs_equation_get_time_eval(ts, mc->w_eq);
    cur2prev = true;

  }

  if (cs_equation_get_space_scheme(mc->w_eq) !=
      CS_SPACE_SCHEME_CDOVB)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Incompatible space discretization.", __func__);

  /* Update pressures fields/arrays */

  _update_pressures(connect, cdoq, update_flag, cur2prev, mc);

  /* Update properties related to soils:
   * - l_rel_permeability, g_rel_permeability
   * - liquid_saturation
   * - capacity: \frac{\partial S_l}{\partial P_c}
   *
   * Either a call to a user-defined function or a predefined function if the
   * soil corresponds to a known model
   */

  cs_gwf_soil_update(time_eval, mesh, connect, cdoq);

  /* Define the liquid saturation in each cell when the liquid saturation has
     been defined on the submesh */

  if (mc->use_definition_on_submesh) {

    const cs_adjacency_t  *c2v = connect->c2v;
    cs_real_t  *l_sat = mc->l_saturation->val;

    /* Compute the average liquid saturation in each cell */

    for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

      cs_real_t  _sat = 0;
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        _sat += mc->l_saturation_submesh[j] * cdoq->pvol_vc[j];
      l_sat[c_id] = _sat/cdoq->cell_vol[c_id];

    } /* Loop on cells */

  } /* Properties on submesh */

  /* Update arrays associated to each term */

  switch (cs_property_get_dim(mc->diff_wl_pty)) {

  case 1: /* Isotropic case */
    /*       ++++++++++++++ */

    if (mc->is_miscible) {
      /*       ======== */

      bft_error(__FILE__, __LINE__, 0,
                "%s: Only the immiscible case is available up to now.",
                __func__);

      //TODO --> cs_gwf_soil_iso_update_mtpf_terms(mc);

    }
    else { /* Immiscible model
              ================ */

      /* In the immiscible case, l_diffusivity_h should be set to 0 and the
         henry constat should be very low */

      assert(mc->l_diffusivity_h < FLT_MIN);
      assert(mc->henry_constant < 1e-12);

      if (mc->use_coupled_solver) {
        /*    ------------------- */

        if (mc->use_diffusion_view_for_darcy)
          _update_iso_itpf_coupled_diffview_terms(mc);
        else
          bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);

      }
      else { /* segregated solver
                ----------------- */

        if (mc->use_diffusion_view_for_darcy) {

          if (mc->use_incremental_solver) {
            if (mc->use_definition_on_submesh)
              bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);
            else
              bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);
          }
          else {
            if (mc->use_definition_on_submesh)
              bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);
            else
              bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);
          }

        }
        else
          bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);

      }

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Only the isotropic case is available up to now.",
              __func__);
    break;

  } /* Switch on the permeability type */

  /* Update the Darcy fluxes. This should be done at the end since one needs
     the new state for:
     - the pressure fields
     - krl, krg, etc.
  */

  _update_darcy_fluxes(connect, cdoq, time_eval, cur2prev, mc);

  return time_eval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined extra-operations for the groundwater flow module in case
 *        of miscible or immiscible two-phase flows in porous media
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      post_flag  requested quantities to be postprocessed
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_extra_op(const cs_cdo_connect_t          *connect,
                    const cs_cdo_quantities_t       *cdoq,
                    cs_flag_t                        post_flag,
                    cs_gwf_tpf_t                    *mc)
{
  assert(mc != NULL);

  if (cs_flag_test(post_flag, CS_GWF_POST_DARCY_FLUX_BALANCE) == false)
    return; /* Nothing to do */

  /* Balance for the Darcy advective flux in the liquid phase */

  cs_gwf_darcy_flux_balance(connect, cdoq,
                            cs_equation_get_param(mc->w_eq),
                            mc->l_darcy);

  /* Balance for the Darcy advective flux in the gas phase */

  cs_gwf_darcy_flux_balance(connect, cdoq,
                            cs_equation_get_param(mc->h_eq),
                            mc->g_darcy);
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
 * \param[in] mc           pointer to the model context structure
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
                      const cs_gwf_tpf_t         *mc,
                      const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *cdoq,
                      const cs_time_step_t       *time_step)
{
  if (mesh_id != CS_POST_MESH_VOLUME)
    return; /* Only postprocessings in the volume are defined */

  assert(mc != NULL);

  if (post_flag & CS_GWF_POST_SOIL_CAPACITY) {

    if (mc->l_capacity == NULL) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf("%s: Requested postprocessing for capacity but not"
                 " allocated.\n This postprocessing is skipped.\n", __func__);
    }
    else {

      bool  use_parent = true;
      cs_real_t  *l_capacity = NULL;

      if (mc->use_definition_on_submesh) {

        const cs_adjacency_t  *c2v = connect->c2v;

        BFT_MALLOC(l_capacity, n_cells, cs_real_t);
        use_parent = false;

        for (cs_lnum_t i = 0; i < n_cells; i++) {

          cs_lnum_t  c_id =
            (cell_ids == NULL || n_cells == cdoq->n_cells) ? i : cell_ids[i];

          l_capacity[c_id] = 0;
          for (cs_lnum_t j = c2v->idx[c_id]; c2v->idx[c_id+1]; j++)
            l_capacity[c_id] = cdoq->pvol_vc[j] * mc->l_capacity[j];
          l_capacity[c_id] /= cdoq->cell_vol[c_id];

        }

      }
      else
        l_capacity = mc->l_capacity;

      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_DEFAULT,
                        "l_capacity",
                        1,
                        false,        /* interlace */
                        use_parent,
                        CS_POST_TYPE_cs_real_t,
                        l_capacity,
                        NULL,
                        NULL,
                        time_step);

      if (l_capacity != mc->l_capacity)
        BFT_FREE(l_capacity);

    }

  } /* Postprocess the soil capacity */

  if (post_flag & CS_GWF_POST_PERMEABILITY) {

    /* permeability = abs_permeability * l_rel_permeability */

    cs_real_t  *permeability = NULL;
    int  dim = cs_property_get_dim(abs_perm);
    int  post_dim = (dim == 1) ? 1 : 9;

    if (dim > 1) {

      BFT_MALLOC(permeability, post_dim*n_cells, cs_real_t);

      for (cs_lnum_t i = 0; i < n_cells; i++) {

        cs_real_t  tensor[3][3];

        cs_property_get_cell_tensor(cell_ids[i],
                                    time_step->t_cur,
                                    abs_perm,
                                    false, /* inversion */
                                    tensor);

        cs_real_t  *_cell_perm = permeability + post_dim*i;
        for (int ki = 0; ki < 3; ki++)
          for (int kj = 0; kj < 3; kj++)
            _cell_perm[3*ki+kj] = tensor[ki][kj];

      }

    }
    else {

      BFT_MALLOC(permeability, n_cells, cs_real_t);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        permeability[c_id] = cs_property_get_cell_value(cell_ids[c_id],
                                                        time_step->t_cur,
                                                        abs_perm);

    }

    assert(mc->l_rel_permeability != NULL);
    for (cs_lnum_t c = 0; c < n_cells; c++)
      for (int k = 0; k < post_dim; k++)
        permeability[post_dim*c+k] *= mc->l_rel_permeability[cell_ids[c]];

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_DEFAULT,
                      "permeability",
                      post_dim,
                      true,     /* interlace */
                      false,    /* use_parent */
                      CS_POST_TYPE_cs_real_t,
                      permeability,
                      NULL,
                      NULL,
                      time_step);

    BFT_FREE(permeability);

  } /* Post-processing of the permeability field */

  if (post_flag & CS_GWF_POST_GAS_MASS_DENSITY) {

    const double  mh_ov_rt =
      mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);

    cs_real_t  *gas_mass_density = NULL;
    BFT_MALLOC(gas_mass_density, n_cells, cs_real_t);

    for (cs_lnum_t c = 0; c < n_cells; c++)
      gas_mass_density[c] = mh_ov_rt * mc->g_pressure_cells[c];

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_DEFAULT,
                      "gas_mass_density",
                      1,
                      true,
                      false,
                      CS_POST_TYPE_cs_real_t,
                      gas_mass_density,
                      NULL,
                      NULL,
                      time_step);

    BFT_FREE(gas_mass_density);

  } /* Post-processing of the gas mass density */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
