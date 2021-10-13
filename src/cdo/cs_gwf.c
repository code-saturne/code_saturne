/*============================================================================
 * Main functions dedicated to groundwater flows
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_cdovb_scaleq.h"
#include "cs_equation_bc.h"
#include "cs_evaluate.h"
#include "cs_field.h"
#include "cs_gwf_priv.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_param_types.h"
#include "cs_physical_constants.h"
#include "cs_post.h"
#include "cs_reco.h"
#include "cs_zone.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_gwf.c

  \brief Main functions dedicated to groundwater flows when using CDO schemes

*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* Redefined names of function from cs_math to get shorter names */
#define _dp3 cs_math_3_dot_product

#define CS_GWF_DBG 0

/*============================================================================
 * Local definitions
 *============================================================================*/

static const char
cs_gwf_model_name[CS_GWF_N_MODEL_TYPES][CS_BASE_STRING_LEN] =
  { N_("Single-phase Richards equation"),
    N_("Two-phase Richards equation")
  };

/*============================================================================
 * Structure definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_gw[] =
  " Stop execution. The structure related to the groundwater module is empty.\n"
  " Please check your settings.\n";

static cs_gwf_t  *cs_gwf_main_structure = NULL;

/*============================================================================
 * Private static inline function prototypes
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the associated boundary Darcy flux for each vertex of
 *         boundary faces.
 *         Case of a vertex-based discretization and single-phase flows in
 *         porous media.
 *
 * \param[in]      t_eval     time at which one performs the evaluation
 * \param[in, out] mc         pointer to the casted modelling context
 */
/*----------------------------------------------------------------------------*/

static void
_spf_update_darcy_vb_flux_at_boundary(cs_real_t                t_eval,
                                      cs_gwf_single_phase_t   *mc)
{
  cs_adv_field_t  *adv = mc->adv_field;

  if (adv->n_bdy_flux_defs > 1 ||
      adv->bdy_flux_defs[0]->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field at the boundary",
              __func__);

  cs_xdef_t  *def = adv->bdy_flux_defs[0];
  cs_xdef_array_context_t  *actx = (cs_xdef_array_context_t *)def->context;
  cs_real_t  *nflx_val = actx->values;

  if (cs_flag_test(actx->loc, cs_flag_dual_closure_byf) == false)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field at the boundary",
              __func__);

  cs_equation_compute_boundary_diff_flux(t_eval, mc->richards, nflx_val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update head values (pressure head or head values for laws)
 *         Case of single-phase flows in porous media.
 *
 * \param[in, out] mc         pointer to a modelling context structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

static void
_spf_update_head(cs_gwf_single_phase_t       *mc,
                 const cs_cdo_quantities_t   *cdoq,
                 const cs_cdo_connect_t      *connect,
                 bool                         cur2prev)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(mc != NULL && gw != NULL);
  const cs_equation_t  *richards = mc->richards;

  cs_param_space_scheme_t r_scheme = cs_equation_get_space_scheme(richards);
  cs_field_t  *hydraulic_head = cs_equation_get_field(richards);
  cs_field_t  *pressure_head = mc->pressure_head;

  if (gw->flag & CS_GWF_RESCALE_HEAD_TO_ZERO_MEAN_VALUE) {

    cs_real_t  domain_integral = 0.;

    switch (r_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      {
        assert(hydraulic_head->location_id ==
               cs_mesh_location_get_id_by_name("vertices"));

        domain_integral =
          cs_evaluate_scal_domain_integral_by_array(cs_flag_primal_vtx,
                                                    hydraulic_head->val);

        const cs_real_t  mean_value = domain_integral / cdoq->vol_tot;

#       pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++)
          hydraulic_head->val[i] -= mean_value;
      }
      break;

    case CS_SPACE_SCHEME_CDOFB:
      {
        assert(hydraulic_head->location_id ==
               cs_mesh_location_get_id_by_name("cells"));

        domain_integral =
          cs_evaluate_scal_domain_integral_by_array(cs_flag_primal_cell,
                                                    hydraulic_head->val);

        const cs_real_t  mean_value = domain_integral / cdoq->vol_tot;
#       pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < cdoq->n_cells; i++)
          hydraulic_head->val[i] -= mean_value;
      }
      break;

    default:
      break; /* Nothing to do */

    }

  } /* Rescale hydraulic_head */

  if (gw->flag & CS_GWF_GRAVITATION) { /* Update the pressure head */

    cs_physical_constants_t  *phys = cs_get_glob_physical_constants();

    /* Sanity checks */
    if (pressure_head == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " The field related to the pressure head is not allocated.");

    /* Copy current field values to previous values */
    if (cur2prev)
      cs_field_current_to_previous(mc->pressure_head);

    switch (r_scheme) {

    case CS_SPACE_SCHEME_CDOVB:

#     pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
        const cs_real_t  gpot = _dp3(cdoq->vtx_coord + 3*i, phys->gravity);
        pressure_head->val[i] = hydraulic_head->val[i] - gpot;
      }

      /* Update head_in_law */
      cs_reco_pv_at_cell_centers(connect->c2v, cdoq, pressure_head->val,
                                 mc->head_in_law);
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      {
        assert(hydraulic_head->location_id ==
               cs_mesh_location_get_id_by_name("vertices"));

#       pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
          const cs_real_t  gpot = _dp3(cdoq->vtx_coord + 3*i, phys->gravity);
          pressure_head->val[i] = hydraulic_head->val[i] - gpot;
        }

        /* Update head_in_law */
        const cs_real_t  *hydraulic_head_cells =
          cs_equation_get_cell_values(richards, false); /* current values */

#       pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {
          const cs_real_t  gpot = _dp3(cdoq->cell_centers + 3*i, phys->gravity);
          mc->head_in_law[i] = hydraulic_head_cells[i] - gpot;
        }

      }
      break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:

#     pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {
        const cs_real_t  gpot = _dp3(cdoq->cell_centers + 3*i, phys->gravity);
        pressure_head->val[i] = hydraulic_head->val[i] - gpot;
      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");

    } /* Switch on space scheme */

  }
  else { /* No gravity effect is taken into account */

    /* Update head_in_law */
    switch(r_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      cs_reco_pv_at_cell_centers(connect->c2v,
                                 cdoq,
                                 hydraulic_head->val,
                                 mc->head_in_law);
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      {
        const cs_real_t  *hydraulic_head_cells =
          cs_equation_get_cell_values(richards, false); /* current values */

        memcpy(mc->head_in_law, hydraulic_head_cells,
               sizeof(cs_real_t)*cdoq->n_cells);
      }
      break;

    default:
      break; /* Nothing to do for CDO-Fb schemes and HHO schemes */

    }  /* Switch on the space scheme related to the Richards equation */

  } /* Gravity is activated or not */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the advection field related to the Darcean flux.
 *         Case of single phase flows in porous media.
 *
 * \param[in]      t_eval     time at which one performs the evaluation
 * \param[in, out] mc         pointer to the casted modelling context
 * \param[in]      cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

static void
_spf_update_darcy_flux(const cs_real_t              t_eval,
                       cs_gwf_single_phase_t       *mc,
                       bool                         cur2prev)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(mc != NULL && gw != NULL);
  const cs_equation_t  *richards = mc->richards;
  const cs_adv_field_t  *adv = mc->adv_field;

  /* Update the velocity field at cell centers induced by the Darcy flux */
  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);

  assert(vel != NULL); /* Sanity check */
  if (cur2prev)
    cs_field_current_to_previous(vel);

  /* Update arrays related to the Darcy flux:
   * Compute the new darcian flux and darcian velocity inside each cell */
  switch (cs_equation_get_space_scheme(richards)) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:

    /* Update the array mc->darcian_flux associated to the advection field */
    if (cs_flag_test(mc->flux_location, cs_flag_dual_face_byc)) {

      assert(mc->darcian_flux != NULL);
      if (adv->definition->type != CS_XDEF_BY_ARRAY)
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid definition of the advection field", __func__);

      cs_equation_compute_diff_flux_cellwise(richards,
                                             mc->flux_location,
                                             t_eval,
                                             mc->darcian_flux);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 2
      cs_dbg_darray_to_listing("DARCIAN_FLUX_DFbyC",
                               connect->c2e->idx[cdoq->n_cells],
                               gw->darcian_flux, 8);
#endif

      /* Set the new values of the vector field at cell centers */
      cs_advection_field_in_cells(mc->adv_field, t_eval, vel->val);


    }
    else if (cs_flag_test(mc->flux_location, cs_flag_primal_cell))
      cs_equation_compute_diff_flux_cellwise(richards,
                                             mc->flux_location,
                                             t_eval,
                                             vel->val);

    /* Update the Darcy flux at the boundary */
    _spf_update_darcy_vb_flux_at_boundary(t_eval, mc);
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
    bft_error(__FILE__, __LINE__, 0, " TODO.");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");

  } /* End of switch */

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
  cs_dbg_darray_to_listing("DARCIAN_FLUX_CELL", 3*cdoq->n_cells, vel->val, 3);
#endif

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
 * \brief  Predefined extra-operations for the groundwater flow module in case
 *         of single phase flows in porous media
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc        pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_spf_extra_op(const cs_cdo_connect_t      *connect,
              const cs_cdo_quantities_t   *cdoq,
              cs_gwf_single_phase_t       *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  const cs_lnum_t  n_b_faces = cdoq->n_b_faces;
  const cs_adv_field_t  *adv = mc->adv_field;
  assert(adv != NULL);
  const cs_field_t  *nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_t  *flux_val = (nflx == NULL) ? mc->darcian_boundary_flux : nflx->val;

  if (flux_val == NULL && n_b_faces > 0) /* No value on which operates */
    return;

  /* Define the balance by zone (using the splitting arising from the settings
     of the boundary conditions for the Richards equation) */

  const cs_equation_param_t  *eqp = cs_equation_get_param(mc->richards);

  bool  *is_counted = NULL;
  BFT_MALLOC(is_counted, n_b_faces, bool);
# pragma omp parallel for if (n_b_faces > CS_THR_MIN)
  for (int i = 0; i < n_b_faces; i++) is_counted[i] = false;

  cs_real_t  *balances = NULL;
  BFT_MALLOC(balances, eqp->n_bc_defs + 1, cs_real_t);

  for (int ibc = 0; ibc < eqp->n_bc_defs; ibc++) {

    const cs_xdef_t  *def = eqp->bc_defs[ibc];
    const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);

    balances[ibc] = 0;

    if (nflx == NULL) { /* The definition of the boundary flux relies on the
                           bf2v adjacency */

#if defined(DEBUG) && !defined(NDEBUG)
      cs_xdef_t  *_def = adv->bdy_flux_defs[0];
      cs_xdef_array_context_t  *actx = (cs_xdef_array_context_t *)_def->context;

      assert(adv->n_bdy_flux_defs == 1 && _def->type == CS_XDEF_BY_ARRAY);
      assert(cs_flag_test(actx->loc, cs_flag_dual_closure_byf) == true);
#endif

      const cs_adjacency_t  *bf2v = connect->bf2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {
        const cs_lnum_t  bf_id = z->elt_ids[i];
        is_counted[bf_id] = true;
        for (cs_lnum_t j = bf2v->idx[bf_id]; j < bf2v->idx[bf_id+1]; j++)
          balances[ibc] += flux_val[j];
      }

    }
    else {

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {
        const cs_lnum_t  bf_id = z->elt_ids[i];
        is_counted[bf_id] = true;
        balances[ibc] += flux_val[bf_id];
      }

    } /* nflux is NULL ? */

  } /* Loop on BC definitions */

  bool  display = false;
  balances[eqp->n_bc_defs] = 0.;
  for (cs_lnum_t bf_id = 0; bf_id < n_b_faces; bf_id++) {
    if (is_counted[bf_id] == false) {

      display = true;
      if (nflx == NULL) {

        const cs_adjacency_t  *bf2v = connect->bf2v;
        for (cs_lnum_t j = bf2v->idx[bf_id]; j < bf2v->idx[bf_id+1]; j++)
          balances[eqp->n_bc_defs] += flux_val[j];

      }
      else
        balances[eqp->n_bc_defs] += flux_val[bf_id];

    } /* Not already counted */
  } /* Loop on boundary faces */

  int display_flag = display ? 1 : 0;

  cs_parall_max(1, CS_INT_TYPE, &display_flag);
  cs_parall_sum(eqp->n_bc_defs + 1, CS_REAL_TYPE, balances);

  /* Output */
  cs_log_printf(CS_LOG_DEFAULT,
                "-b- Balance of the Darcy flux across the boundary zones:\n");

  for (int ibc = 0; ibc < eqp->n_bc_defs; ibc++) {
    const cs_zone_t  *z = cs_boundary_zone_by_id((eqp->bc_defs[ibc])->z_id);
    cs_log_printf(CS_LOG_DEFAULT, "-b- %-32s: % -5.3e\n",
                  z->name, balances[ibc]);
  }

  if (display_flag > 0)
    cs_log_printf(CS_LOG_DEFAULT, "-b- %-32s: % -5.3e\n",
                  "Remaining part of the boundary", balances[eqp->n_bc_defs]);

  BFT_FREE(is_counted);
  BFT_FREE(balances);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the steady-state of the groundwater flows module.
 *         Nothing is done if all equations are unsteady.
 *         Case of single-phase flows in porous media
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_spf_compute_steady_state(const cs_mesh_t              *mesh,
                          const cs_time_step_t         *time_step,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *cdoq,
                          cs_gwf_single_phase_t        *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);
  assert(gw->model == CS_GWF_MODEL_SINGLE_PHASE_RICHARDS);

  cs_equation_t  *richards = mc->richards;

  /* Sanity check */
  assert(richards != NULL);
  assert(cs_equation_get_type(richards) == CS_EQUATION_TYPE_GROUNDWATER);

  /* Build and solve the linear system related to the Richards equations */
  if (cs_equation_is_steady(richards) ||
      gw->flag & CS_GWF_FORCE_RICHARDS_ITERATIONS) {

    /* Solve the algebraic system */
    cs_equation_solve_steady_state(mesh, richards);

    /* Update the variables related to the groundwater flow system */
    cs_gwf_update(mesh, connect, cdoq, time_step, true);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for the groundwater flows module.
 *         Case of single-phase flows in porous media.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_spf_compute(const cs_mesh_t              *mesh,
             const cs_time_step_t         *time_step,
             const cs_cdo_connect_t       *connect,
             const cs_cdo_quantities_t    *cdoq,
             cs_gwf_single_phase_t        *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);
  assert(gw->model == CS_GWF_MODEL_SINGLE_PHASE_RICHARDS);

  cs_equation_t  *richards = mc->richards;

  /* Sanity check */
  assert(richards != NULL);
  assert(cs_equation_get_type(richards) == CS_EQUATION_TYPE_GROUNDWATER);

  bool cur2prev = true;

  /* Build and solve the linear system related to the Richards equations */
  if (!cs_equation_is_steady(richards) ||
      gw->flag & CS_GWF_FORCE_RICHARDS_ITERATIONS) {

    /* Solve the algebraic system. By default, a current to previous operation
       is performed */
    cs_equation_solve(cur2prev, mesh, richards);

    /* Update the variables related to the groundwater flow system */
    cs_gwf_update(mesh, connect, cdoq, time_step, cur2prev);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the update step in the case of a single-phase flow model in
 *         porous media
 *
 * \param[in]       mesh        pointer to a cs_mesh_t structure
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]       ts          pointer to a cs_time_step_t structure
 * \param[in]       cur2prev    true or false
 * \param[in, out]  mc          pointer to the casted model context
 *
 *\return the value at which one has to evaluate the properties
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_spf_updates(const cs_mesh_t             *mesh,
             const cs_cdo_connect_t      *connect,
             const cs_cdo_quantities_t   *quant,
             const cs_time_step_t        *ts,
             bool                         cur2prev,
             cs_gwf_single_phase_t       *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  const cs_real_t  dt_cur = ts->dt[0];

  cs_real_t  time_eval = ts->t_cur;
  if (cur2prev) {

    /* Define the time at which one evaluates the properties */
    cs_param_time_scheme_t
      rt_scheme = cs_equation_get_time_scheme(mc->richards);
    cs_real_t  theta = -1;

    switch (rt_scheme) {

    case CS_TIME_SCHEME_STEADY:
    case CS_TIME_N_SCHEMES:
      /* Look for tracer equations */
      for (int ieq = 0; ieq < gw->n_tracers; ieq++)
        theta = fmax(theta,
                     cs_equation_get_theta_time_val(gw->tracers[ieq]->eq));
      if (theta > 0)
        time_eval = ts->t_cur + theta*dt_cur;
      else
        time_eval = ts->t_cur;
      break;

    default:
      theta = cs_equation_get_theta_time_val(mc->richards);
      time_eval = ts->t_cur + theta*dt_cur;
      break;

    } /* End of switch on time scheme for the Richards equation */

  } /* cur2prev */

  /* Update head */
  _spf_update_head(mc, quant, connect, cur2prev);

  /* Update the advection field related to the groundwater flow module */
  _spf_update_darcy_flux(time_eval, mc, cur2prev);

  /* Update properties related to soils.
     Handle the moisture content field: do something only if the moisture
     content is set as unsteady, otherwise keep values already set */
  if (gw->flag & CS_GWF_SOIL_ALL_SATURATED) {

    /* Handle only the moisture field if this is the initialization */
    if (cur2prev == false)
      cs_property_eval_at_cells(time_eval,
                                mc->moisture_content,
                                mc->moisture_field->val);

  }
  else {

    /* Handle the permeability, the moisture content and the soil capacity */
    assert(gw->permea_field != NULL);
    if (cur2prev) {

      cs_field_current_to_previous(gw->permea_field);
      cs_field_current_to_previous(mc->moisture_field);
      if (mc->capacity_field != NULL)
        cs_field_current_to_previous(mc->capacity_field);

    }

    const int n_soils = cs_gwf_get_n_soils();
    for (int i = 0; i < n_soils; i++) {

      cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(i);
      const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);

      soil->update_properties(time_eval, mesh, connect, quant,
                              mc->head_in_law,
                              zone,
                              soil->input);

    } /* Loop on soils */

  } /* Not all saturated */

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
  cs_dbg_darray_to_listing("MOISTURE_CONTENT",
                           quant->n_cells,
                           mc->moisture_field->val, 8);
#endif

  return time_eval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the last setup step in the case of a single-phase flow
 *         model in porous media
 *
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       quant      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_spf_last_setup(const cs_cdo_connect_t     *connect,
                const cs_cdo_quantities_t  *quant,
                cs_gwf_single_phase_t      *mc)
{
  size_t  array_size;
  cs_flag_t  array_location;

  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  const cs_field_t  *hydraulic_head = cs_equation_get_field(mc->richards);
  const cs_param_space_scheme_t  richards_scheme =
    cs_equation_get_space_scheme(mc->richards);
  const cs_lnum_t  n_cells = connect->n_cells;

  cs_field_t  *cell_adv_field =
    cs_advection_field_get_field(mc->adv_field, CS_MESH_LOCATION_CELLS);
  assert(cell_adv_field != NULL);

  /* Set the Darcian flux (in the volume and at the boundary) */
  switch (richards_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_adjacency_t  *bf2v = connect->bf2v;

      /* Define the flux of the advection field at the boundary */
      array_size = bf2v->idx[quant->n_b_faces];
      BFT_MALLOC(mc->darcian_boundary_flux, array_size, cs_real_t);
      memset(mc->darcian_boundary_flux, 0, array_size*sizeof(cs_real_t));

      array_location = CS_FLAG_SCALAR | cs_flag_dual_closure_byf;

      /* Do not transfer the ownership */
      cs_advection_field_def_boundary_flux_by_array(mc->adv_field,
                                                    NULL,
                                                    array_location,
                                                    mc->darcian_boundary_flux,
                                                    false,
                                                    bf2v->idx);

      /* Define the advection field in the volume */
      if (cs_flag_test(mc->flux_location, cs_flag_dual_face_byc)) {

        /* Darcian flux settings */
        const cs_adjacency_t  *c2e = connect->c2e;

        array_size = c2e->idx[n_cells];
        BFT_MALLOC(mc->darcian_flux, array_size, cs_real_t);
        memset(mc->darcian_flux, 0, array_size*sizeof(cs_real_t));

        array_location = CS_FLAG_SCALAR | mc->flux_location;

        /* Do not transfer the ownership */
        cs_advection_field_def_by_array(mc->adv_field,
                                        array_location,
                                        mc->darcian_flux,
                                        false, /* transfer ownership */
                                        c2e->idx);

        /* Reset the type of advection field */
        if (mc->adv_field->status & CS_ADVECTION_FIELD_TYPE_VELOCITY_VECTOR)
          mc->adv_field->status -= CS_ADVECTION_FIELD_TYPE_VELOCITY_VECTOR;
        mc->adv_field->status |= CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;

      }
      else if (cs_flag_test(mc->flux_location, cs_flag_primal_cell)) {

        cs_advection_field_def_by_field(mc->adv_field, cell_adv_field);

        /* Reset the type of advection field */
        if (mc->adv_field->status & CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX)
          mc->adv_field->status -= CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;
        mc->adv_field->status |= CS_ADVECTION_FIELD_TYPE_VELOCITY_VECTOR;

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid location for defining the Darcian flux.",
                  __func__);

      /* Allocate a head array defined at cells used to update the soil
         properties */
      BFT_MALLOC(mc->head_in_law, n_cells, cs_real_t);

    }
    break;

  case CS_SPACE_SCHEME_CDOFB:   /* TODO */

    /* Set the head array defined at cells used to update the soil properties */
    if (gw->flag & CS_GWF_GRAVITATION)
      mc->head_in_law = mc->pressure_head->val;
    else
      mc->head_in_law = hydraulic_head->val;

    bft_error(__FILE__, __LINE__, 0,
              " %s: Fb space scheme not fully implemented.", __func__);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);
    break;

  } /* Switch on Richards scheme */

  /* Set permeability, moisture content and soil capacity according to the
     soil settings */
  if (gw->flag & CS_GWF_SOIL_ALL_SATURATED) {

    cs_gwf_soil_set_all_saturated(gw->permeability,
                                  mc->moisture_content,
                                  mc->moisture_field);

    if (gw->permea_field != NULL) /* Fill the values of the permeability field
                                     for post-processing */

      cs_property_eval_at_cells(0, /* Should be a steady-state property */
                                gw->permeability,
                                gw->permea_field->val);

  }
  else
    cs_gwf_soil_set_by_field(gw->permeability,
                             gw->permea_field,
                             mc->moisture_content,
                             mc->moisture_field,
                             mc->soil_capacity,
                             mc->capacity_field);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the initial setup step in the case of a single-phase flow
 *         model in porous media
 *
 * \param[in]       pty_has_previous   unsteady behavior of soil properties
 * \param[in, out]  mc                 pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_spf_init_setup(_Bool                    pty_has_previous,
                cs_gwf_single_phase_t   *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  if (mc->richards == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The Richards equation is not defined. Stop execution.\n",
              __func__);

  const bool has_previous = cs_equation_is_steady(mc->richards) ? false : true;
  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");
  const cs_param_space_scheme_t  space_scheme =
    cs_equation_get_space_scheme(mc->richards);

  /* Handle gravity effects */
  if (gw->flag & CS_GWF_GRAVITATION) {

    switch (space_scheme) {
    case CS_SPACE_SCHEME_CDOVB:
    case CS_SPACE_SCHEME_CDOVCB:
      mc->pressure_head = cs_field_create("pressure_head",
                                          field_mask,
                                          v_loc_id,
                                          1,
                                          has_previous);
      break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
      mc->pressure_head = cs_field_create("pressure_head",
                                          field_mask,
                                          c_loc_id,
                                          1,
                                          has_previous);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);
    }

    cs_field_set_key_int(mc->pressure_head, log_key, 1);
    cs_field_set_key_int(mc->pressure_head, post_key, 1);

  } /* Gravitation effect is activated */


  /* Create a moisture field attached to cells */
  int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;
  mc->moisture_field = cs_field_create("moisture_content",
                                       pty_mask,
                                       c_loc_id,
                                       1,   /* dimension */
                                       pty_has_previous);

  if (pty_has_previous)
    cs_field_set_key_int(mc->moisture_field, log_key, 1);
  if (gw->post_flag & CS_GWF_POST_MOISTURE)
    cs_field_set_key_int(mc->moisture_field, post_key, 1);

  if (!(gw->flag & CS_GWF_SOIL_ALL_SATURATED) ||
      gw->post_flag & CS_GWF_POST_PERMEABILITY) {

    /* Set the dimension of the permeability */
    int  permeability_dim = 0;  /* not set by default */
    if (gw->permeability->type & CS_PROPERTY_ISO)
      permeability_dim = 1;
    else if (gw->permeability->type & CS_PROPERTY_ORTHO)
      permeability_dim = 3;
    else if (gw->permeability->type & CS_PROPERTY_ANISO)
      permeability_dim = 9;
    assert(permeability_dim != 0);

    gw->permea_field = cs_field_create("permeability",
                                       pty_mask,
                                       c_loc_id,
                                       permeability_dim,   /* dimension */
                                       pty_has_previous);

    cs_field_set_key_int(gw->permea_field, log_key, 1);
    if (gw->post_flag & CS_GWF_POST_PERMEABILITY)
      cs_field_set_key_int(gw->permea_field, post_key, 1);

  } /* Need to create a field for the permeability */

  /* Create a capacity field attached to cells */
  if (gw->flag & CS_GWF_RICHARDS_UNSTEADY) {

    mc->capacity_field = cs_field_create("soil_capacity",
                                         pty_mask,
                                         c_loc_id,
                                         1,   /* dimension */
                                         pty_has_previous);

    cs_field_set_key_int(mc->capacity_field, log_key, 1);
    if (gw->post_flag & CS_GWF_POST_CAPACITY)
      cs_field_set_key_int(mc->capacity_field, post_key, 1);

  }

  /* Add default post-processing related to groundwater flow module */
  cs_post_add_time_mesh_dep_output(cs_gwf_extra_post_single_phase, gw);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a modelling context related to single-phase flows
 *
 * \return a pointer to a new allocated cs_gwf_single_phase_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_single_phase_t *
_create_single_phase_context(void)
{
  cs_gwf_single_phase_t  *mc = NULL;

  BFT_MALLOC(mc, 1, cs_gwf_single_phase_t);

  mc->richards = NULL;

  mc->moisture_content = NULL;
  mc->moisture_field = NULL;

  mc->pressure_head = NULL;
  mc->head_in_law = NULL;

  mc->soil_capacity = NULL;
  mc->capacity_field = NULL;

  mc->flux_location = cs_flag_dual_face_byc;
  mc->darcian_flux = NULL;
  mc->darcian_boundary_flux = NULL;
  mc->adv_field = NULL;

  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the groundwater flow module in case
 *         of two phase flows in porous media
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc        pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_extra_op(const cs_cdo_connect_t      *connect,
              const cs_cdo_quantities_t   *cdoq,
              cs_gwf_two_phase_t          *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the steady-state of the groundwater flows module.
 *         Nothing is done if all equations are unsteady.
 *         Case of single-phase flows in porous media
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_compute_steady_state(const cs_mesh_t              *mesh,
                          const cs_time_step_t         *time_step,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *cdoq,
                          cs_gwf_two_phase_t           *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);
  assert(gw->model == CS_GWF_MODEL_TWO_PHASE_RICHARDS);

  cs_equation_t  *water_eq = mc->water_eq;
  cs_equation_t  *gcomp_eq = mc->gcomp_eq;

  /* Sanity check */
  assert(water_eq != NULL && gcomp_eq != NULL);
  assert(cs_equation_get_type(water_eq) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(cs_equation_get_type(gcomp_eq) == CS_EQUATION_TYPE_GROUNDWATER);

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for the groundwater flows module.
 *         Case of single-phase flows in porous media.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_compute(const cs_mesh_t              *mesh,
             const cs_time_step_t         *time_step,
             const cs_cdo_connect_t       *connect,
             const cs_cdo_quantities_t    *cdoq,
             cs_gwf_two_phase_t           *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);
  assert(gw->model == CS_GWF_MODEL_TWO_PHASE_RICHARDS);

  cs_equation_t  *water_eq = mc->water_eq;
  cs_equation_t  *gcomp_eq = mc->gcomp_eq;

  /* Sanity check */
  assert(water_eq != NULL && gcomp_eq != NULL);
  assert(cs_equation_get_type(water_eq) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(cs_equation_get_type(gcomp_eq) == CS_EQUATION_TYPE_GROUNDWATER);

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the update step in the case of a two-phase flow model in
 *         porous media
 *
 * \param[in]       mesh        pointer to a cs_mesh_t structure
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]       ts          pointer to a cs_time_step_t structure
 * \param[in]       cur2prev    true or false
 * \param[in, out]  mc          pointer to the casted model context
 *
 *\return the value at which one has to evaluate the properties
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_tpf_updates(const cs_mesh_t             *mesh,
             const cs_cdo_connect_t      *connect,
             const cs_cdo_quantities_t   *quant,
             const cs_time_step_t        *ts,
             bool                         cur2prev,
             cs_gwf_two_phase_t          *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  /* TODO */
  cs_real_t  time_eval = 0.;

  return time_eval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the initial setup step in the case of a two-phase flow
 *         model in porous media
 *
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       quant      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_last_setup(const cs_cdo_connect_t      *connect,
                const cs_cdo_quantities_t   *quant,
                cs_gwf_two_phase_t          *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the initial setup step in the case of a two-phase flows
 *         model in porous media
 *
 * \param[in]       pty_has_previous   unsteady behavior of soil properties
 * \param[in, out]  mc                 pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_init_setup(_Bool                  pty_has_previous,
                cs_gwf_two_phase_t    *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a modelling context related to two-phase flows
 *
 * \return a pointer to a new allocated cs_gwf_two_phase_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_two_phase_t *
_create_two_phase_context(void)
{
  cs_gwf_two_phase_t  *mc = NULL;

  BFT_MALLOC(mc, 1, cs_gwf_two_phase_t);

  mc->water_eq = NULL;
  mc->gcomp_eq = NULL;

  mc->darcy_l_field = NULL;
  mc->darcy_g_field = NULL;

  mc->flux_location = cs_flag_dual_face_byc;

  mc->darcian_l_flux = NULL;
  mc->darcian_g_flux = NULL;
  mc->darcian_l_boundary_flux = NULL;
  mc->darcian_g_boundary_flux = NULL;

  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a structure dedicated to manage groundwater flows
 *
 * \return a pointer to a new allocated cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_t *
_gwf_create(void)
{
  cs_gwf_t  *gw = NULL;

  BFT_MALLOC(gw, 1, cs_gwf_t);

  /* Default initialization */
  gw->model = CS_GWF_N_MODEL_TYPES;
  gw->flag = 0;
  gw->post_flag = CS_GWF_POST_DARCY_FLUX_BALANCE;

  gw->permeability = NULL;
  gw->permea_field = NULL;

  gw->model_context = NULL;

  /* Tracer part */
  gw->n_tracers = 0;
  gw->tracers = NULL;
  gw->finalize_tracer_setup = NULL;
  gw->add_tracer_terms = NULL;

  return gw;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the groundwater flow module has been activated
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_is_activated(void)
{
  if (cs_gwf_main_structure == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the module dedicated to groundwater flows
 *
 * \param[in]   permeability_type     type of permeability (iso, ortho...)
 * \param[in]   model                 type of physical modelling
 * \param[in]   option_flag           optional flag to specify this module
 *
 * \return a pointer to a new allocated groundwater flow structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_activate(cs_property_type_t           pty_type,
                cs_gwf_model_type_t          model,
                cs_gwf_option_flag_t         option_flag)
{
  cs_gwf_t  *gw = _gwf_create();

  gw->model = model;
  gw->flag = option_flag;

  /* Add a property related to the diffusion term of the Richards eq. */
  gw->permeability = cs_property_add("permeability", pty_type);

  /* Add an advection field related to the darcian flux stemming from the
     Richards equation */
  cs_advection_field_status_t  adv_status =
    CS_ADVECTION_FIELD_GWF | CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;

  switch (model) {

  case CS_GWF_MODEL_SINGLE_PHASE_RICHARDS:
    {
      /* Model context (mc) */
      cs_gwf_single_phase_t  *mc = _create_single_phase_context();

      /* Create a new equation structure for Richards' equation */
      mc->richards = cs_equation_add("Richards",
                                     "hydraulic_head",
                                     CS_EQUATION_TYPE_GROUNDWATER,
                                     1,
                                     CS_PARAM_BC_HMG_NEUMANN);

      cs_equation_param_t  *eqp = cs_equation_get_param(mc->richards);

      /* Associate permeability to the diffusion property of the Richards eq. */
      cs_equation_add_diffusion(eqp, gw->permeability);

      /* Add a property related to the moisture content */
      mc->moisture_content = cs_property_add("moisture_content",
                                             CS_PROPERTY_ISO);

      /* Add a property related to the unsteady term of the Richards eq. */
      if (option_flag & CS_GWF_RICHARDS_UNSTEADY) {

        mc->soil_capacity = cs_property_add("soil_capacity", CS_PROPERTY_ISO);

        /* Associate soil_capacity to the unsteady term of the Richards eq. */
        cs_equation_add_time(eqp, mc->soil_capacity);

      }
      else { /* Steady-state case */

        adv_status |= CS_ADVECTION_FIELD_STEADY;

      }

      /* Add the advection field which is the Darcy flux */
      mc->adv_field = cs_advection_field_add(CS_GWF_ADV_FIELD_NAME, adv_status);

      /* Set the pointer to the modelling context */
      gw->model_context = mc;
    }
    break; /* Single-phase */

  case CS_GWF_MODEL_TWO_PHASE_RICHARDS:
    {
      /* Model context (mc) */
      cs_gwf_two_phase_t  *mc = _create_two_phase_context();

      /* TODO */

      /* Set the pointer to the modelling context */
      gw->model_context = mc;
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);

  }

  /* Store the pointer to the groundawater flow structure */
  cs_gwf_main_structure = gw;

  return gw;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to groundwater flows
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_destroy_all(void)
{
  if (cs_gwf_main_structure == NULL)
    return NULL;

  cs_gwf_t  *gw = cs_gwf_main_structure;

  switch (gw->model) {

  case CS_GWF_MODEL_SINGLE_PHASE_RICHARDS:
    {
      cs_gwf_single_phase_t  *mc = gw->model_context;

      /* darcian_flux and darcian_boundary_flux are allocated only if the
         related advection field is defined by array. At the definition step,
         the GWF module has kept the ownership of the lifecycle of the
         darcian_flux and darcian_boundary_flux arrays. In this case, the
         lifecycle is not managed by the definition */
      BFT_FREE(mc->darcian_boundary_flux);
      BFT_FREE(mc->darcian_flux);
      BFT_FREE(mc->head_in_law);

      /* Free the modelling context */
      BFT_FREE(mc);
    }
    break;

  case CS_GWF_MODEL_TWO_PHASE_RICHARDS:
    {
      cs_gwf_two_phase_t  *mc = gw->model_context;

      BFT_FREE(mc->darcian_l_boundary_flux);
      BFT_FREE(mc->darcian_g_boundary_flux);
      BFT_FREE(mc->darcian_l_flux);
      BFT_FREE(mc->darcian_g_flux);

      /* Free the modelling context */
      BFT_FREE(mc);
    }
    break;

  default:
    /* Nothing to do */
    break;
  }

  /* Free all soils */
  cs_gwf_soil_free_all();

  /* Manage the tracer-related members */
  for (int i = 0; i < gw->n_tracers; i++)
    gw->tracers[i] = cs_gwf_tracer_free(gw->tracers[i]);
  BFT_FREE(gw->tracers);
  BFT_FREE(gw->finalize_tracer_setup);
  BFT_FREE(gw->add_tracer_terms);

  BFT_FREE(gw);

  /* Fields, equations, advection fields and properties are freed elsewhere */
  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_log_setup(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL)
    return;
  assert(gw->model < CS_GWF_N_MODEL_TYPES);

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the groundwater module\n");
  cs_log_printf(CS_LOG_SETUP, "%s", cs_sep_h1);

  /* Tracers */
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Number of tracer equations: %d\n", gw->n_tracers);

  /* Main options */
  if (gw->model == CS_GWF_MODEL_SINGLE_PHASE_RICHARDS) {

    cs_gwf_single_phase_t  *mc = (cs_gwf_single_phase_t *)gw->model_context;

    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Model: %s\n", cs_gwf_model_name[gw->model]);
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Darcy flux location: %s\n",
                  cs_flag_str_location(mc->flux_location));

  }
  else if (gw->model == CS_GWF_MODEL_TWO_PHASE_RICHARDS) {

    cs_gwf_two_phase_t  *mc = (cs_gwf_two_phase_t *)gw->model_context;

    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Model: %s\n", cs_gwf_model_name[gw->model]);
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Darcy flux location: %s\n",
                  cs_flag_str_location(mc->flux_location));

  }

  /* Display information on the general options */
  if (gw->flag & CS_GWF_GRAVITATION)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Gravitation: **True**\n");
  else
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Gravitation: **False**\n");

  if (gw->flag & CS_GWF_ENFORCE_DIVERGENCE_FREE)
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Enforce the divergence-free constraint"
                  " for the Darcy flux\n");
  if (gw->flag & CS_GWF_FORCE_RICHARDS_ITERATIONS)
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Force to solve Richards equation"
                  " at each time step\n");
  if (gw->flag & CS_GWF_RESCALE_HEAD_TO_ZERO_MEAN_VALUE)
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Rescale head w.r.t zero mean value\n");

  /* Display information on the post-processing options */
  bool  post_capacity = (gw->post_flag & CS_GWF_POST_CAPACITY) ? true : false;
  bool  post_moisture = (gw->post_flag & CS_GWF_POST_MOISTURE) ? true : false;
  bool  post_perm = (gw->post_flag & CS_GWF_POST_PERMEABILITY) ? true : false;
  cs_log_printf(CS_LOG_SETUP, "  * GWF | Post: Capacity %s Moisture %s"
                " Permeability %s\n",
                cs_base_strtf(post_capacity), cs_base_strtf(post_moisture),
                cs_base_strtf(post_perm));

  bool  do_balance =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_BALANCE) ? true : false;
  bool  do_divergence =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_DIVERGENCE) ? true : false;
  bool  post_boundary =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY) ? true : false;
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Darcy Flux: Balance %s Divergence %s"
                " At boundary faces: %s\n",
                cs_base_strtf(do_balance), cs_base_strtf(do_divergence),
                cs_base_strtf(post_boundary));

  /* Soils */
  if (gw->flag & CS_GWF_SOIL_ALL_SATURATED)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | All soils are saturated\n");
  if (gw->flag & CS_GWF_SOIL_PROPERTY_UNSTEADY)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Unsteady soil properties\n");

  /* Detailed setup of the soil properties */
  cs_gwf_soil_log_setup();

  /* Detailed setup of the tracer equations */
  for (int i = 0; i < gw->n_tracers; i++)
    cs_gwf_tracer_log_setup(gw->tracers[i]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the flag dedicated to the post-processing of the GWF module
 *
 * \param[in]  post_flag             flag to set
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_post_options(cs_flag_t       post_flag)
{
  if (cs_gwf_main_structure == NULL)
    return;

  cs_gwf_t  *gw = cs_gwf_main_structure;

  gw->post_flag = post_flag;
  if (gw->post_flag & CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY) {

    if (gw->model == CS_GWF_MODEL_SINGLE_PHASE_RICHARDS) {
      cs_gwf_single_phase_t  *mc = (cs_gwf_single_phase_t *)gw->model_context;
      mc->adv_field->status |= CS_ADVECTION_FIELD_DEFINE_AT_BOUNDARY_FACES;
    }

    /* TODO: Two phase flow case */

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a particular type of unsteady advection-diffusion
 *         reaction eq.
 *         Tracer is advected thanks to the darcian velocity and
 *         diffusion/reaction parameters result from a physical modelling.
 *         Terms solved in the equation are activated according to the settings.
 *         The advection field corresponds to that of the liquid phase.
 *
 * \param[in]  tr_model   physical modelling to consider (0 = default settings)
 * \param[in]  eq_name    name of the tracer equation
 * \param[in]  var_name   name of the related variable
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_add_tracer(cs_gwf_tracer_model_t     tr_model,
                  const char               *eq_name,
                  const char               *var_name)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (tr_model & CS_GWF_TRACER_USER)
    bft_error(__FILE__, __LINE__, 0,
              "%s: User-defined is not allowed in this context.\n"
              " Please consider cs_gwf_add_user_tracer() instead.", __func__);

  /* Set the advection field structure */
  cs_adv_field_t  *adv = NULL;
  if (gw->model == CS_GWF_MODEL_SINGLE_PHASE_RICHARDS) {
    cs_gwf_single_phase_t  *mc = (cs_gwf_single_phase_t *)gw->model_context;
    adv = mc->adv_field;
  }
  else if (gw->model == CS_GWF_MODEL_TWO_PHASE_RICHARDS) {
    cs_gwf_two_phase_t  *mc = (cs_gwf_two_phase_t *)gw->model_context;
    adv = mc->darcy_l_field;
  }
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model.\n", __func__);

  int  tr_id = gw->n_tracers;
  cs_gwf_tracer_t  *tracer = cs_gwf_tracer_init(tr_id,
                                                eq_name,
                                                var_name,
                                                adv,
                                                tr_model);

  gw->n_tracers += 1;
  BFT_REALLOC(gw->tracers, gw->n_tracers, cs_gwf_tracer_t *);
  BFT_REALLOC(gw->finalize_tracer_setup,
              gw->n_tracers, cs_gwf_tracer_setup_t *);
  BFT_REALLOC(gw->add_tracer_terms,
              gw->n_tracers, cs_gwf_tracer_add_terms_t *);

  gw->tracers[tr_id] = tracer;
  gw->finalize_tracer_setup[tr_id] = cs_gwf_tracer_setup;
  gw->add_tracer_terms[tr_id] = cs_gwf_tracer_add_terms;

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a particular type of unsteady advection-diffusion
 *         reaction eq.
 *         Tracer is advected thanks to the darcian velocity and
 *         diffusion/reaction parameters result from a physical modelling.
 *         Terms are activated according to the settings.
 *         Modelling of the tracer parameters are left to the user
 *
 * \param[in]   eq_name     name of the tracer equation
 * \param[in]   var_name    name of the related variable
 * \param[in]   setup       function pointer (predefined prototype)
 * \param[in]   add_terms   function pointer (predefined prototype)
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_add_user_tracer(const char                  *eq_name,
                       const char                  *var_name,
                       cs_gwf_tracer_setup_t       *setup,
                       cs_gwf_tracer_add_terms_t   *add_terms)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  /* Set the advection field structure */
  cs_adv_field_t  *adv = NULL;
  if (gw->model == CS_GWF_MODEL_SINGLE_PHASE_RICHARDS) {
    cs_gwf_single_phase_t  *mc = (cs_gwf_single_phase_t *)gw->model_context;
    adv = mc->adv_field;
  }
  else if (gw->model == CS_GWF_MODEL_TWO_PHASE_RICHARDS) {
    cs_gwf_two_phase_t  *mc = (cs_gwf_two_phase_t *)gw->model_context;
    adv = mc->darcy_l_field;
  }
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model.\n", __func__);

  int  tr_id = gw->n_tracers;
  cs_gwf_tracer_t  *tracer = cs_gwf_tracer_init(tr_id,
                                                eq_name,
                                                var_name,
                                                adv,
                                                CS_GWF_TRACER_USER);

  gw->n_tracers += 1;
  BFT_REALLOC(gw->tracers, gw->n_tracers, cs_gwf_tracer_t *);
  BFT_REALLOC(gw->finalize_tracer_setup,
              gw->n_tracers, cs_gwf_tracer_setup_t *);
  BFT_REALLOC(gw->add_tracer_terms,
              gw->n_tracers, cs_gwf_tracer_add_terms_t *);

  gw->tracers[tr_id] = tracer;
  gw->finalize_tracer_setup[tr_id] = setup;
  gw->add_tracer_terms[tr_id] = add_terms;

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to the cs_gwf_tracer_t structure associated to
 *         the name given as parameter
 *
 * \param[in]  eq_name    name of the tracer equation
 *
 * \return the pointer to a cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_tracer_by_name(const char   *eq_name)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (eq_name == NULL)
    return NULL;

  for (int i = 0; i < gw->n_tracers; i++) {
    cs_gwf_tracer_t  *tracer = gw->tracers[i];
    const char *name_to_cmp = cs_equation_get_name(tracer->eq);
    if (strcmp(eq_name, name_to_cmp) == 0)
      return tracer;
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for the Richards equation and the related
 *         equations defining the groundwater flow module.
 *         At this stage, all soils have been defined.
 *         Create new cs_field_t structures according to the setting
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_setup(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  /* Sanity checks */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));
  cs_gwf_soil_check();

  /* 1. Analyze the type of soils to handle.
     Detect if all soils are considered as saturated. Update flags. If this
     not the case, create new fields and properties are time-dependent.
  */
  bool  pty_has_previous;

  if (cs_gwf_soil_all_saturated()) {
    gw->flag |= CS_GWF_SOIL_ALL_SATURATED;
    pty_has_previous = false;
  }
  else {
    gw->flag |= CS_GWF_SOIL_PROPERTY_UNSTEADY;
    pty_has_previous = true;
  }

  switch (gw->model) {

  case CS_GWF_MODEL_SINGLE_PHASE_RICHARDS:
    _spf_init_setup(pty_has_previous, gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE_RICHARDS:
    _tpf_init_setup(pty_has_previous, gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add new terms if needed (such as diffusion or reaction) to tracer
 *         equations according to the settings
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_add_tracer_terms(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  /* Sanity checks */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  /* Loop on tracer equations */
  for (int i = 0; i < gw->n_tracers; i++)
    gw->add_tracer_terms[i](gw->tracers[i]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last initialization step of the groundwater flow module
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_finalize_setup(const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  switch (gw->model) {

  case CS_GWF_MODEL_SINGLE_PHASE_RICHARDS:
    _spf_last_setup(connect, quant, gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE_RICHARDS:
    _tpf_last_setup(connect, quant, gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Store the soil id for each cell */
  cs_gwf_soil_check();
  cs_gwf_build_cell2soil(quant->n_cells);

  /* Loop on tracer equations. Link the advection field to each tracer
     equation */
  for (int i = 0; i < gw->n_tracers; i++)
    gw->finalize_tracer_setup[i](connect, quant, gw->tracers[i]);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the groundwater system (pressure head, head in law, moisture
 *         content, darcian velocity, soil capacity or permeability if needed)
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_update(const cs_mesh_t             *mesh,
              const cs_cdo_connect_t      *connect,
              const cs_cdo_quantities_t   *quant,
              const cs_time_step_t        *ts,
              bool                         cur2prev)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  /* Sanity checks */
  if (gw == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Groundwater module is not allocated.", __func__);

  cs_real_t  time_eval = 0.;
  switch (gw->model) {

  case CS_GWF_MODEL_SINGLE_PHASE_RICHARDS:
    time_eval = _spf_updates(mesh, connect, quant, ts,
                             cur2prev,
                             gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE_RICHARDS:
    time_eval = _tpf_updates(mesh, connect, quant, ts,
                             cur2prev,
                             gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Update the diffusivity tensor associated to each tracer equation since the
     Darcy velocity may have changed */
  for (int i = 0; i < gw->n_tracers; i++) {

    cs_gwf_tracer_t  *tracer = gw->tracers[i];
    if (tracer->update_diff_tensor != NULL)
      tracer->update_diff_tensor(tracer, time_eval, mesh, connect, quant);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the steady-state of the groundwater flows module.
 *         Nothing is done if all equations are unsteady.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_compute_steady_state(const cs_mesh_t              *mesh,
                            const cs_time_step_t         *time_step,
                            const cs_cdo_connect_t       *connect,
                            const cs_cdo_quantities_t    *cdoq)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  switch (gw->model) {

  case CS_GWF_MODEL_SINGLE_PHASE_RICHARDS:
    _spf_compute_steady_state(mesh, time_step, connect, cdoq,
                              gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE_RICHARDS:
    _tpf_compute_steady_state(mesh, time_step, connect, cdoq,
                              gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Compute tracers */
  /* --------------- */

  for (int i = 0; i < gw->n_tracers; i++) {

    cs_gwf_tracer_t  *tracer = gw->tracers[i];

    if (cs_equation_is_steady(tracer->eq)) {

      /* Solve the algebraic system */
      cs_equation_solve_steady_state(mesh, tracer->eq);

      if (tracer->update_precipitation != NULL)
        tracer->update_precipitation(tracer,
                                     time_step->t_cur,
                                     mesh, connect, cdoq);

    } /* Solve this equation which is steady */

  } /* Loop on tracer equations */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the system related to groundwater flows module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_compute(const cs_mesh_t              *mesh,
               const cs_time_step_t         *time_step,
               const cs_cdo_connect_t       *connect,
               const cs_cdo_quantities_t    *cdoq)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  switch (gw->model) {

  case CS_GWF_MODEL_SINGLE_PHASE_RICHARDS:
    _spf_compute(mesh, time_step, connect, cdoq, gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE_RICHARDS:
    _tpf_compute(mesh, time_step, connect, cdoq, gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Compute tracers */
  /* --------------- */
  bool cur2prev = true;

  for (int i = 0; i < gw->n_tracers; i++) {

    cs_gwf_tracer_t  *tracer = gw->tracers[i];

    if (!cs_equation_is_steady(tracer->eq)) { /* unsteady ? */

      /* Solve the algebraic system. By default, a current to previous operation
         is performed */
      cs_equation_solve(cur2prev, mesh, tracer->eq);

      if (tracer->update_precipitation != NULL)
        tracer->update_precipitation(tracer,
                                     time_step->t_cur,
                                     mesh, connect, cdoq);

    } /* Solve this equation which is unsteady */

  } /* Loop on tracer equations */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a given set of cells of the field related
 *         to a tracer equation. This integral turns out to be exact for linear
 *         functions.
 *
 * \param[in]    connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]    cdoq      pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]    tracer    pointer to a \ref cs_gwf_tracer_t structure
 * \param[in]    z_name    name of the volumic zone where the integral is done
 *                         (if NULL or "" all cells are considered)
 *
 * \return the value of the integral
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_integrate_tracer(const cs_cdo_connect_t     *connect,
                        const cs_cdo_quantities_t  *cdoq,
                        const cs_gwf_tracer_t      *tracer,
                        const char                 *z_name)
{
  const int  z_id = cs_get_vol_zone_id(z_name);
  const cs_zone_t  *z = cs_volume_zone_by_id(z_id);
  const short int  *cell2soil = cs_gwf_get_cell2soil();

  const cs_field_t  *moist = cs_field_by_name("moisture_content");
  if (moist == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: \"moisture_content\" not defined",
              __func__);
  assert(tracer != NULL);

  const cs_real_t  *moisture_val = moist->val;
  const cs_equation_param_t  *tr_eqp = cs_equation_get_param(tracer->eq);

  cs_gwf_tracer_input_t  *sti = (cs_gwf_tracer_input_t *)tracer->input;
  cs_real_t  int_value = 0.0;

  switch (tr_eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(tracer->eq,
                                                               false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t  c_id = (z->elt_ids == NULL) ? i : z->elt_ids[i];

        cs_real_t  _int_value = 0.;
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {
          _int_value += cdoq->dcell_vol[j] * v_vals[c2v->ids[j]];
        }

        int_value +=
          (moisture_val[c_id] + sti->rho_kd[cell2soil[c_id]]) * _int_value;

      } /* Loop on selected cells */

    }
    break; /* CS_SPACE_SCHEME_CDOVB */

  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(tracer->eq,
                                                               false);
      const cs_real_t  *c_vals = cs_equation_get_cell_values(tracer->eq,
                                                             false);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t  c_id = (z->elt_ids == NULL) ? i : z->elt_ids[i];

        /* Shares between cell and vertex unknowns:
           - the cell unknown stands for 1/4 of the cell volume
           - the vertex unknown stands for 3/4 of the dual cell volume
         */
        cs_real_t  _int_value = 0.25*cdoq->cell_vol[c_id]*c_vals[c_id];
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {
          _int_value += 0.75 * cdoq->dcell_vol[j] * v_vals[c2v->ids[j]];
        }

        int_value +=
          (moisture_val[c_id] + sti->rho_kd[cell2soil[c_id]]) * _int_value;

      } /* Loop on selected cells */

    }
    break; /* CS_SPACE_SCHEME_CDOVCB */

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme", __func__);
    break;

  } /* End of switch */

  /* Parallel synchronization */
  if (cs_glob_n_ranks > 1)
    cs_parall_sum(1, CS_REAL_TYPE, &int_value);

  return int_value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the groundwater flow module
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_extra_op(const cs_cdo_connect_t      *connect,
                const cs_cdo_quantities_t   *cdoq)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL)
    return;
  if (cs_flag_test(gw->post_flag, CS_GWF_POST_DARCY_FLUX_BALANCE) == false)
    return; /* Nothing to do */

  switch (gw->model) {

  case CS_GWF_MODEL_SINGLE_PHASE_RICHARDS:
    _spf_extra_op(connect, cdoq, gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE_RICHARDS:
    _tpf_extra_op(connect, cdoq, gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the groundwater flow module
 *         in case of single-phase flows in porous media.
 *         Prototype of this function is given since it is a function pointer
 *         defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_extra_post_single_phase(void                      *input,
                               int                        mesh_id,
                               int                        cat_id,
                               int                        ent_flag[5],
                               cs_lnum_t                  n_cells,
                               cs_lnum_t                  n_i_faces,
                               cs_lnum_t                  n_b_faces,
                               const cs_lnum_t            cell_ids[],
                               const cs_lnum_t            i_face_ids[],
                               const cs_lnum_t            b_face_ids[],
                               const cs_time_step_t      *time_step)
{
  CS_UNUSED(cat_id);
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);
  CS_UNUSED(time_step);

  if (input == NULL)
    return;

  const cs_gwf_t  *gw = (const cs_gwf_t *)input;
  const cs_gwf_single_phase_t  *mc =
    (const cs_gwf_single_phase_t  *)gw->model_context;

  if (mesh_id == CS_POST_MESH_VOLUME) {

    if (gw->post_flag & CS_GWF_POST_DARCY_FLUX_DIVERGENCE) {

      /* Only case avalaible up to now */
      if (cs_advection_field_get_deftype(mc->adv_field) == CS_XDEF_BY_ARRAY) {

        cs_real_t  *divergence =
          cs_advection_field_divergence_at_vertices(mc->adv_field,
                                                    time_step->t_cur);

        cs_post_write_vertex_var(mesh_id,
                                 CS_POST_WRITER_DEFAULT,
                                 "darcy_flux_divergence",
                                 1,
                                 false,
                                 false,
                                 CS_POST_TYPE_cs_real_t,
                                 divergence,
                                 time_step);

        BFT_FREE(divergence);
      }

    } /* Post-processing of the divergence is requested */
  } /* volume mesh id */

}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
