/*============================================================================
 * Main functions dedicated to groundwater flows
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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

  \brief Main high-level functions dedicated to groundwater flows when using
         CDO schemes

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
  { N_("Saturated single-phase model"),
    N_("Unsaturated single-phase model"),
    N_("Two-phase model (capillary/gas pressure)")
  };

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_gw[] =
  " Stop execution. The structure related to the groundwater module is empty.\n"
  " Please check your settings.\n";

static cs_gwf_t  *cs_gwf_main_structure = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Estimate the time at which the evaluation of properties has to be
 *         done
 *
 * \param[in]   ts      pointer to a cs_time_step_t structure
 * \param[in]   eq      pointer to an equation structure
 *
 * \return the time value at which one has to perform evaluation
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_get_time_eval(const cs_time_step_t        *ts,
               cs_equation_t               *eq)

{
  cs_real_t  time_eval = ts->t_cur;

  const cs_real_t  dt_cur = ts->dt[0];

  /* Define the time at which one evaluates the properties */

  cs_param_time_scheme_t  time_scheme = cs_equation_get_time_scheme(eq);
  cs_real_t  theta = -1;

  switch (time_scheme) {

  case CS_TIME_SCHEME_STEADY:
  case CS_TIME_N_SCHEMES:

    /* Scan tracer equations */

    cs_gwf_tracer_get_time_theta_max();

    if (theta > 0)
      time_eval = ts->t_cur + theta*dt_cur;
    else
      time_eval = ts->t_cur;
    break;

  default:
    theta = cs_equation_get_theta_time_val(eq);
    time_eval = ts->t_cur + theta*dt_cur;
    break;

  } /* End of switch on the time scheme for the main equation */

  return time_eval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the advection field related to the Darcy flux in the liquid
 *        phase
 *
 * \param[in]  gw     pointer to the main (high-level) GWF structure
 *
 * \return a pointer to a cs_adv_field_t structure or NULL
 */
/*----------------------------------------------------------------------------*/

static cs_adv_field_t *
_get_l_adv_field(const cs_gwf_t   *gw)
{
  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    {
      cs_gwf_saturated_single_phase_t  *mc = gw->model_context;

      return mc->darcy->adv_field;
    }
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    {
      cs_gwf_unsaturated_single_phase_t  *mc = gw->model_context;

      return mc->darcy->adv_field;
    }
    break;

  case CS_GWF_MODEL_TWO_PHASE:
    {
      cs_gwf_miscible_two_phase_t  *mc = gw->model_context;

      if (mc->l_darcy != NULL)
        return mc->l_darcy->adv_field;
    }
    break;

  default:
    break;

  } /* Switch on the model */

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update head values (pressure head or head values for laws)
 *         Case of single-phase flows in porous media (saturated or not).
 *
 * \param[in]      cdoq           pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect        pointer to a cs_cdo_connect_t structure
 * \param[in]      richards       pointer to the Richards equation
 * \param[in, out] pressure_head  pressure head field
 * \param[in, out] head_in_law    values of the head used in law
 * \param[in]      cur2prev       true or false
 */
/*----------------------------------------------------------------------------*/

static void
_spf_update_head(const cs_cdo_quantities_t   *cdoq,
                 const cs_cdo_connect_t      *connect,
                 const cs_equation_t         *richards,
                 cs_field_t                  *pressure_head,
                 cs_real_t                    head_in_law[],
                 bool                         cur2prev)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL);

  cs_param_space_scheme_t r_scheme = cs_equation_get_space_scheme(richards);
  cs_field_t  *hydraulic_head = cs_equation_get_field(richards);

  if (gw->flag & CS_GWF_RESCALE_HEAD_TO_ZERO_MEAN_VALUE) {

    switch (r_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      {
        assert(hydraulic_head->location_id ==
               cs_mesh_location_get_id_by_name("vertices"));

        cs_real_t  domain_integral =
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

        cs_real_t  domain_integral =
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

  if (gw->flag & CS_GWF_GRAVITATION) { /* Update the pressure head (and if
                                          needed head_in_law) */

    cs_physical_constants_t  *phys = cs_get_glob_physical_constants();

    if (pressure_head == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " The field related to the pressure head is not allocated.");

    /* Copy current field values to previous values */

    if (cur2prev)
      cs_field_current_to_previous(pressure_head);

    switch (r_scheme) {

    case CS_SPACE_SCHEME_CDOVB:

#     pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
        const cs_real_t  gpot = _dp3(cdoq->vtx_coord + 3*i, phys->gravity);
        pressure_head->val[i] = hydraulic_head->val[i] - gpot;
      }

      /* Update head_in_law */

      if (head_in_law != NULL)
        cs_reco_pv_at_cell_centers(connect->c2v, cdoq, pressure_head->val,
                                   head_in_law);
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

        if (head_in_law != NULL) {

          const cs_real_t  *hydraulic_head_cells =
            cs_equation_get_cell_values(richards, false); /* current values */

#         pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
          for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {
            const cs_real_t  gpot = _dp3(cdoq->cell_centers + 3*i,
                                         phys->gravity);
            head_in_law[i] = hydraulic_head_cells[i] - gpot;
          }

        } /* head_in_law */

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

      /* Head in law points either to hydraulic_head->val or pressure_head->val */

    default:
      bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");

    } /* Switch on space scheme */

  }
  else { /* No gravity effect is taken into account */

    if (head_in_law == NULL)
      return;

    /* Update head_in_law */

    switch(r_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      cs_reco_pv_at_cell_centers(connect->c2v,
                                 cdoq,
                                 hydraulic_head->val,
                                 head_in_law);
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      {
        const cs_real_t  *hydraulic_head_cells =
          cs_equation_get_cell_values(richards, false); /* current values */

        memcpy(head_in_law, hydraulic_head_cells,
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
 * \brief  Compute the steady-state of the groundwater flows module.
 *         Nothing is done if all equations are unsteady.
 *         Case of unstaturated/saturated single-phase flows in porous media.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] richards   pointer to the Richards equation
 */
/*----------------------------------------------------------------------------*/

static void
_spf_compute_steady_state(const cs_mesh_t                  *mesh,
                          const cs_time_step_t             *time_step,
                          const cs_cdo_connect_t           *connect,
                          const cs_cdo_quantities_t        *cdoq,
                          cs_equation_t                    *richards)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && richards != NULL);

  assert(gw->model == CS_GWF_MODEL_SATURATED_SINGLE_PHASE);
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
 * \brief  Compute the new state for the groundwater flows module.
 *         Case of unstaturated/saturated single-phase flows in porous media.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] richards   pointer to the Richards equation
 */
/*----------------------------------------------------------------------------*/

static void
_spf_compute(const cs_mesh_t                    *mesh,
             const cs_time_step_t               *time_step,
             const cs_cdo_connect_t             *connect,
             const cs_cdo_quantities_t          *cdoq,
             cs_equation_t                      *richards)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  assert(gw != NULL && richards != NULL);

  assert(gw->model == CS_GWF_MODEL_SATURATED_SINGLE_PHASE ||
         gw->model == CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE);
  assert(cs_equation_get_type(richards) == CS_EQUATION_TYPE_GROUNDWATER);

  bool cur2prev = true;

  /* Build and solve the linear system related to the Richards equations */

  if (!cs_equation_is_steady(richards)) {

    /* Solve the algebraic system. By default, a current to previous operation
       is performed */

    cs_equation_solve(cur2prev, mesh, richards);

    /* Update the variables related to the groundwater flow system */

    cs_gwf_update(mesh, connect, cdoq, time_step, cur2prev);

  }
}

/* ==========================================================================
 * Functions related to the model of saturated single phase flows (SSPF)
 * ========================================================================== */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the modelling context for the model of
 *         saturated single-phase flows
 *
 * \return a pointer to a new allocated cs_gwf_saturated_single_phase_t struct
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_saturated_single_phase_t *
_sspf_init_context(void)
{
  cs_gwf_saturated_single_phase_t  *mc = NULL;

  BFT_MALLOC(mc, 1, cs_gwf_saturated_single_phase_t);

  mc->pressure_head = NULL;

  /* Create a new equation structure for Richards' equation */

  mc->richards = cs_equation_add("Richards",       /* equation name */
                                 "hydraulic_head", /* variable name */
                                 CS_EQUATION_TYPE_GROUNDWATER,
                                 1,
                                 CS_PARAM_BC_HMG_NEUMANN);

  /* Define the Darcy flux structure
     Add an advection field related to the darcian flux stemming from the
     Richards equation. This advection field is steady since the head is
     steady in this model.
  */

  mc->darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);

  cs_advection_field_status_t  adv_status = CS_ADVECTION_FIELD_GWF |
    CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX | CS_ADVECTION_FIELD_STEADY;

  mc->darcy->adv_field = cs_advection_field_add("darcy_field",
                                                adv_status);

  /* Add a property related to the moisture content (It should be a constant
     in each soil) */

  mc->moisture_content = cs_property_add("moisture_content", CS_PROPERTY_ISO);

  /* Add the diffusion term to the Richards equation by associating the
     absolute permeability to the diffusion property of the Richards eq. */

  cs_equation_param_t  *eqp = cs_equation_get_param(mc->richards);
  cs_gwf_t  *gw = cs_gwf_main_structure;

  cs_equation_add_diffusion(eqp, gw->abs_permeability);

  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the modelling context for the model of saturated single-phase
 *        flows
 *
 * \param[in, out] p_mc   pointer of pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_sspf_free_context(cs_gwf_saturated_single_phase_t   **p_mc)
{
  if (p_mc == NULL)
    return;
  if ((*p_mc) == NULL)
    return;

  cs_gwf_saturated_single_phase_t  *mc = *p_mc;

  cs_gwf_darcy_flux_free(&(mc->darcy));

  BFT_FREE(mc);
  *p_mc = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the modelling context of saturated
 *        single-phase flows
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_sspf_log_context(cs_gwf_saturated_single_phase_t   *mc)
{
  if (mc == NULL)
    return;

  cs_gwf_darcy_flux_log(mc->darcy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup related to the modelling context of saturated single-phase
 *        flows. At this stage, all soils have been defined and equation
 *        parameters are set.
 *
 * \param[in, out] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_sspf_init_setup(cs_gwf_saturated_single_phase_t   *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (mc == NULL)
    return;
  if (mc->richards == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The Richards equation is not defined. Stop execution.\n",
              __func__);
  assert(cs_equation_is_steady(mc->richards) == true);

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
                                          false); /* has_previous */
      break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
      mc->pressure_head = cs_field_create("pressure_head",
                                          field_mask,
                                          c_loc_id,
                                          1,
                                          false); /* has_previous */
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);
    }

    cs_field_set_key_int(mc->pressure_head, log_key, 1);
    cs_field_set_key_int(mc->pressure_head, post_key, 1);

  } /* Gravitation effect is activated */

  /* Add default post-processing related to groundwater flow module */

  cs_post_add_time_mesh_dep_output(cs_gwf_extra_post_sspf, gw);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the last setup step in the case of a saturated single-phase
 *         flow model in porous media
 *
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       quant      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_sspf_finalize_setup(const cs_cdo_connect_t            *connect,
                     const cs_cdo_quantities_t         *quant,
                     cs_gwf_saturated_single_phase_t   *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  const cs_param_space_scheme_t  richards_scheme =
    cs_equation_get_space_scheme(mc->richards);

  /* Set the Darcian flux (in the volume and at the boundary) */

  cs_gwf_darcy_flux_define(connect, quant, richards_scheme, mc->darcy);

  /* Set the parameter values thanks to the permeability and the moisture
     content */

  cs_gwf_soil_saturated_set_param(gw->abs_permeability, mc->moisture_content);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the update step in the case of a saturated single-phase
 *         flow model in porous media
 *
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
_sspf_updates(const cs_cdo_connect_t            *connect,
              const cs_cdo_quantities_t         *quant,
              const cs_time_step_t              *ts,
              bool                               cur2prev,
              cs_gwf_saturated_single_phase_t   *mc)
{
  cs_real_t  time_eval = ts->t_cur;
  if (cur2prev)
    time_eval = _get_time_eval(ts, mc->richards);

  /* Update head */

  _spf_update_head(quant, connect, mc->richards,
                   mc->pressure_head,
                   NULL,      /* there is no head_in_law for saturated soils */
                   cur2prev);

  /* Update the advection field related to the groundwater flow module */

  cs_gwf_darcy_flux_update(time_eval, mc->richards, cur2prev, mc->darcy);

  /* Properties are constant in saturated soils. Therefore, there is no need fo
     an update of the associated properties (permeability and moisture
     content) */

  return time_eval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the groundwater flow module in case
 *         of saturated single phase flows in porous media
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc        pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_sspf_extra_op(const cs_cdo_connect_t                *connect,
               const cs_cdo_quantities_t             *cdoq,
               cs_gwf_saturated_single_phase_t       *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  if (cs_flag_test(gw->post_flag, CS_GWF_POST_DARCY_FLUX_BALANCE) == false)
    return; /* Nothing to do */

  cs_gwf_darcy_flux_balance(connect, cdoq,
                            cs_equation_get_param(mc->richards),
                            mc->darcy);
}

/* ==========================================================================
 * Functions related to the model of unsaturated single phase flows (USPF)
 * ========================================================================== */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the modelling context for the model of
 *         unsaturated single-phase flows
 *
 * \param[in]   pty_type        type of permeability (iso, ortho...)
 *
 * \return a pointer to a new allocated cs_gwf_unsaturated_single_phase_t struct
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_unsaturated_single_phase_t *
_uspf_init_context(cs_property_type_t           pty_type)
{
  cs_gwf_unsaturated_single_phase_t  *mc = NULL;

  BFT_MALLOC(mc, 1, cs_gwf_unsaturated_single_phase_t);

  mc->permeability_field = NULL;
  mc->moisture_field = NULL;
  mc->capacity_field = NULL;
  mc->pressure_head = NULL;
  mc->head_in_law = NULL;

  /* Create a new equation structure for Richards' equation */

  mc->richards = cs_equation_add("Richards",       /* equation name */
                                 "hydraulic_head", /* variable name */
                                 CS_EQUATION_TYPE_GROUNDWATER,
                                 1,
                                 CS_PARAM_BC_HMG_NEUMANN);

  /* Define the Darcy flux structure
     Add an advection field related to the darcian flux stemming from the
     Richards equation. This advection field is steady since the head is
     steady in this model.
  */

  mc->darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);

  cs_advection_field_status_t  adv_status = CS_ADVECTION_FIELD_GWF |
    CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;

  mc->darcy->adv_field = cs_advection_field_add("darcy_field",
                                                adv_status);

  /* Add several properties:
   * - the moisture content
   * - the soil capacity
   * - the (full) permeability = rel_permeability * abs_permeability
   */

  mc->moisture_content = cs_property_add("moisture_content", CS_PROPERTY_ISO);

  mc->soil_capacity = cs_property_add("soil_capacity", CS_PROPERTY_ISO);

  mc->permeability = cs_property_add("permeability", pty_type);

  /* Add the time + diffusion terms to the Richards equation by associating the
     full permeability to the diffusion property of the Richards eq. and the
     soil capacity to the unsteady term */

  cs_equation_param_t  *eqp = cs_equation_get_param(mc->richards);

  cs_equation_add_diffusion(eqp, mc->permeability);
  cs_equation_add_time(eqp, mc->soil_capacity);

  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the modelling context for the model of unsaturated single-phase
 *        flows
 *
 * \param[in, out] p_mc   pointer of pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_uspf_free_context(cs_gwf_unsaturated_single_phase_t   **p_mc)
{
  if (p_mc == NULL)
    return;
  if (*p_mc == NULL)
    return;

  cs_gwf_unsaturated_single_phase_t  *mc = *p_mc;

  cs_gwf_darcy_flux_free(&(mc->darcy));

  BFT_FREE(mc);
  *p_mc = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the modelling context of unsaturated
 *        single-phase flows
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_uspf_log_context(cs_gwf_unsaturated_single_phase_t   *mc)
{
  if (mc == NULL)
    return;

  cs_gwf_darcy_flux_log(mc->darcy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup related to the modelling context of unsaturated single-phase
 *        flows. At this stage, all soils have been defined and equation
 *        parameters are set.
 *
 * \param[in, out] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_uspf_init_setup(cs_gwf_unsaturated_single_phase_t   *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (mc == NULL)
    return;
  if (mc->richards == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The Richards equation is not defined. Stop execution.\n",
              __func__);
  assert(cs_equation_is_steady(mc->richards) == false);

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
                                          true); /* has_previous */
      break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
      mc->pressure_head = cs_field_create("pressure_head",
                                          field_mask,
                                          c_loc_id,
                                          1,
                                          true); /* has_previous */
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);
    }

    cs_field_set_key_int(mc->pressure_head, log_key, 1);
    cs_field_set_key_int(mc->pressure_head, post_key, 1);

  } /* Gravitation effect is activated */

  /* Create fields at cells for properties */

  int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;
  mc->moisture_field = cs_field_create("liquid_saturation",
                                       pty_mask,
                                       c_loc_id,
                                       1,     /* dimension */
                                       true); /* has_previous */

  cs_field_set_key_int(mc->moisture_field, log_key, 1);
  if (gw->post_flag & CS_GWF_POST_LIQUID_SATURATION)
    cs_field_set_key_int(mc->moisture_field, post_key, 1);

  int  permeability_dim = cs_property_get_dim(gw->abs_permeability);

  mc->permeability_field = cs_field_create("permeability",
                                           pty_mask,
                                           c_loc_id,
                                           permeability_dim,
                                           true); /* has_previous */

  if (gw->post_flag & CS_GWF_POST_PERMEABILITY) {
    cs_field_set_key_int(mc->permeability_field, log_key, 1);
    cs_field_set_key_int(mc->permeability_field, post_key, 1);
  }

  mc->capacity_field = cs_field_create("soil_capacity",
                                       pty_mask,
                                       c_loc_id,
                                       1,   /* dimension */
                                       true);

  cs_field_set_key_int(mc->capacity_field, log_key, 1);
  if (gw->post_flag & CS_GWF_POST_SOIL_CAPACITY)
    cs_field_set_key_int(mc->capacity_field, post_key, 1);

  /* Add default post-processing related to groundwater flow module */

  cs_post_add_time_mesh_dep_output(cs_gwf_extra_post_uspf, gw);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the last setup step in the case of a unsaturated single-phase
 *        flow model in porous media
 *
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       quant      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_uspf_finalize_setup(const cs_cdo_connect_t              *connect,
                     const cs_cdo_quantities_t           *quant,
                     cs_gwf_unsaturated_single_phase_t   *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  const cs_field_t  *hydraulic_head = cs_equation_get_field(mc->richards);
  const cs_param_space_scheme_t  richards_scheme =
    cs_equation_get_space_scheme(mc->richards);
  const cs_lnum_t  n_cells = connect->n_cells;

  /* Set the Darcian flux (in the volume and at the boundary) */

  cs_gwf_darcy_flux_define(connect, quant, richards_scheme, mc->darcy);

  /* Allocate a head array defined at cells and used to update the soil
     properties */

  switch (richards_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    BFT_MALLOC(mc->head_in_law, n_cells, cs_real_t);
#if defined(DEBUG) && !defined(NDEBUG)
    memset(mc->head_in_law, 0, sizeof(cs_real_t)*n_cells);
#endif
    break;

  case CS_SPACE_SCHEME_CDOFB:

    if (gw->flag & CS_GWF_GRAVITATION)
      mc->head_in_law = mc->pressure_head->val;
    else
      mc->head_in_law = hydraulic_head->val;

    bft_error(__FILE__, __LINE__, 0,
              "%s: Fb space scheme not fully implemented.", __func__);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme.", __func__);
    break;

  } /* Switch on Richards scheme */

  /* Define the permeability property using a field */

  cs_property_def_by_field(mc->permeability, mc->permeability_field);

  /* Define the moisture content (liquid saturation) property using a field */

  cs_property_def_by_field(mc->moisture_content, mc->moisture_field);

  /* Define the soil capacity property using a field */

  cs_property_def_by_field(mc->soil_capacity, mc->capacity_field);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the update step in the case of a unsaturated single-phase
 *         flow model in porous media
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
_uspf_updates(const cs_mesh_t                    *mesh,
              const cs_cdo_connect_t             *connect,
              const cs_cdo_quantities_t          *quant,
              const cs_time_step_t               *ts,
              bool                                cur2prev,
              cs_gwf_unsaturated_single_phase_t  *mc)
{
  cs_real_t  time_eval = ts->t_cur;
  if (cur2prev)
    time_eval = _get_time_eval(ts, mc->richards);

  /* Update head */

  _spf_update_head(quant, connect, mc->richards,
                   mc->pressure_head,
                   mc->head_in_law,
                   cur2prev);

  /* Update the advection field related to the groundwater flow module */

  cs_gwf_darcy_flux_update(time_eval, mc->richards, cur2prev, mc->darcy);

  /* Update properties related to soils.
   * Handle the permeability, the moisture content and the soil capacity
   */

  if (cur2prev) {

    cs_field_current_to_previous(mc->permeability_field);
    cs_field_current_to_previous(mc->moisture_field);
    cs_field_current_to_previous(mc->capacity_field);

  }

  /* Update soil properties with the new head values */

  cs_gwf_soil_update(time_eval, mesh, connect, quant);

  return time_eval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the groundwater flow module in case
 *         of unsaturated single phase flows in porous media
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc        pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_uspf_extra_op(const cs_cdo_connect_t                *connect,
               const cs_cdo_quantities_t             *cdoq,
               cs_gwf_unsaturated_single_phase_t     *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  if (cs_flag_test(gw->post_flag, CS_GWF_POST_DARCY_FLUX_BALANCE) == false)
    return; /* Nothing to do */

  cs_gwf_darcy_flux_balance(connect, cdoq,
                            cs_equation_get_param(mc->richards),
                            mc->darcy);
}

/* ==========================================================================
 * Functions related to the model of two phase flows (TPF)
 * ========================================================================== */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the modelling context for the model of
 *         (unsaturated) miscible two-phase flows
 *
 * \param[in]   pty_type        type of permeability (iso, ortho...)
 *
 * \return a pointer to a new allocated cs_gwf_miscible_two_phase_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_miscible_two_phase_t *
_mtpf_init_context(cs_property_type_t           pty_type)
{
  cs_gwf_miscible_two_phase_t  *mc = NULL;

  BFT_MALLOC(mc, 1, cs_gwf_miscible_two_phase_t);

  /* Arrays */

  mc->time_w_eq_array = NULL;
  mc->diff_w_eq_array = NULL;
  mc->time_h_eq_array = NULL;
  mc->diff_hl_eq_array = NULL;
  mc->diff_hg_eq_array = NULL;

  mc->l_rel_permeability = NULL;
  mc->g_rel_permeability = NULL;

  /* Darcy flux (not used up to now) */

  mc->l_darcy = NULL;
  mc->g_darcy = NULL;

  /* Parameters (default values) */

  mc->l_mass_density = 1000;
  mc->l_viscosity = 1e-3;
  mc->g_viscosity = 2e-5;
  mc->l_diffusivity_h = 1e-10;
  mc->w_molar_mass = 18e-3;
  mc->h_molar_mass = 3e-3;
  mc->ref_temperature = 280;    /* in Kelvin */
  mc->henry_constant = 1e-20;   /* immiscible case */

  /* Create a new equation for the water conservation */

  mc->w_eq = cs_equation_add("w_conservation",   /* equation name */
                             "liquid_pressure",  /* variable name */
                             CS_EQUATION_TYPE_GROUNDWATER,
                             1,
                             CS_PARAM_BC_HMG_NEUMANN);

  /* Create a new equation for the hydrogen conservation */

  mc->h_eq = cs_equation_add("h_conservation", /* equation name */
                             "gas_pressure",   /* variable name */
                             CS_EQUATION_TYPE_GROUNDWATER,
                             1,
                             CS_PARAM_BC_HMG_NEUMANN);

  /* Add a 2x2 system of coupled equations and define each block */

  mc->system = cs_equation_system_create("TwoPhaseFlow system", 2);

  /* Set the (0,0)-block */

  cs_equation_system_set_equation(0, mc->w_eq, mc->system);

  /* Set the (1,1)-block */

  cs_equation_system_set_equation(1, mc->h_eq, mc->system);

  /* Create and set the (0,1)-block */

  mc->wh_eqp = cs_equation_param_create("wh_cross_term",
                                        CS_EQUATION_TYPE_GROUNDWATER,
                                        1,
                                        CS_PARAM_BC_HMG_NEUMANN);

  cs_equation_system_set_param(0, 1, mc->wh_eqp, mc->system);

  /* Create and set the (1,0)-block */

  mc->hw_eqp = cs_equation_param_create("hw_cross_term",
                                        CS_EQUATION_TYPE_GROUNDWATER,
                                        1,
                                        CS_PARAM_BC_HMG_NEUMANN);

  cs_equation_system_set_param(1, 0, mc->hw_eqp, mc->system);

  /* Add properties:
   * - unsteady term for w_eq
   * - diffusion term for w_eq
   * - unsteady term for h_eq
   * - diffusion term for h_eq in the gas phase
   * - diffusion term for h_eq in the liquid phase (cross-term)
   */

  mc->time_w_eq_pty = cs_property_add("time_w_eq_pty", CS_PROPERTY_ISO);
  mc->diff_w_eq_pty = cs_property_add("diff_w_eq_pty", pty_type);
  mc->time_h_eq_pty = cs_property_add("time_h_eq_pty", CS_PROPERTY_ISO);
  mc->diff_hg_eq_pty = cs_property_add("diff_hg_eq_pty", pty_type);
  mc->diff_hl_eq_pty = cs_property_add("diff_hl_eq_pty", pty_type);

  /* Associate properties with equation to define new terms in these
     equations */

  cs_equation_param_t  *w_eqp = cs_equation_get_param(mc->w_eq);

  cs_equation_add_time(w_eqp, mc->time_w_eq_pty);
  cs_equation_add_diffusion(w_eqp, mc->diff_w_eq_pty);

  cs_equation_param_t  *h_eqp = cs_equation_get_param(mc->h_eq);

  cs_equation_add_time(h_eqp, mc->time_h_eq_pty);
  cs_equation_add_diffusion(h_eqp, mc->diff_hg_eq_pty);

  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the modelling context for the model of miscible two-phase
 *        flows
 *
 * \param[in, out] p_mc   pointer of pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_mtpf_free_context(cs_gwf_miscible_two_phase_t  **p_mc)
{
  if (p_mc == NULL)
    return;
  if (*p_mc == NULL)
    return;

  cs_gwf_miscible_two_phase_t  *mc = *p_mc;

  /* wh_eqp and hw_eqp are freed inside the next function */

  cs_equation_system_free(&(mc->system));

  cs_gwf_darcy_flux_free(&(mc->l_darcy));
  cs_gwf_darcy_flux_free(&(mc->g_darcy));

  BFT_FREE(mc->time_w_eq_array);
  BFT_FREE(mc->diff_w_eq_array);
  BFT_FREE(mc->time_h_eq_array);
  BFT_FREE(mc->diff_hl_eq_array);
  BFT_FREE(mc->diff_hg_eq_array);

  BFT_FREE(mc->l_rel_permeability);
  BFT_FREE(mc->g_rel_permeability);

  BFT_FREE(mc);
  *p_mc = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the modelling context of miscible
 *        two-phase flows
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_mtpf_log_context(cs_gwf_miscible_two_phase_t   *mc)
{
  if (mc == NULL)
    return;

  cs_gwf_darcy_flux_log(mc->l_darcy);
  cs_gwf_darcy_flux_log(mc->g_darcy);

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Water mass density: %5.3e, viscosity: %5.3e,"
                " molar mass: %5.3e\n", mc->l_mass_density, mc->l_viscosity,
                mc->w_molar_mass);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Main gas component: viscosity: %5.3e, diffusivity"
                " in the liquid phase: %5.3e, molar mass: %5.3e\n",
                mc->g_viscosity, mc->l_diffusivity_h, mc->h_molar_mass);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Reference temperature: %5.3e K\n",
                mc->ref_temperature);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Henry constant: %5.3e\n",
                mc->henry_constant);

  /* Log the system of equations */

  cs_equation_system_log(mc->system);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the initial setup step in the case of a miscible two-phase
 *         flows model in porous media
 *
 * \param[in, out]  mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_mtpf_init_setup(cs_gwf_miscible_two_phase_t    *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  if (mc->w_eq == NULL || mc->h_eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Equations are not defined for this model. Stop execution.\n",
              __func__);

  assert(cs_equation_is_steady(mc->w_eq) == false);
  assert(cs_equation_is_steady(mc->h_eq) == false);

  /* Set fields related to variables */

  mc->l_pressure = cs_equation_get_field(mc->w_eq);
  mc->g_pressure = cs_equation_get_field(mc->h_eq);

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");
  const cs_param_space_scheme_t  space_scheme =
    cs_equation_get_space_scheme(mc->w_eq);
  assert(space_scheme == cs_equation_get_space_scheme(mc->h_eq));

  /* One has to be consistent with the location of DoFs for the w_eq and h_eq
   * which are respectively related to the l_pressure and g_pressure */

  int loc_id = c_loc_id;
  if (space_scheme == CS_SPACE_SCHEME_CDOVB ||
      space_scheme == CS_SPACE_SCHEME_CDOVCB)
    loc_id = v_loc_id;

  mc->c_pressure = cs_field_create("capillarity_pressure",
                                   field_mask,
                                   loc_id,
                                   1,
                                   true); /* has_previous */

  cs_field_set_key_int(mc->c_pressure, log_key, 1);
  cs_field_set_key_int(mc->c_pressure, post_key, 1);

  /* Create a liquid saturation field attached to cells: S_l */

  int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;

  mc->l_saturation = cs_field_create("liquid_saturation",
                                     pty_mask,
                                     c_loc_id,
                                     1,     /* dimension */
                                     true); /* has_previous */

  cs_field_set_key_int(mc->l_saturation, log_key, 1);
  if (gw->post_flag & CS_GWF_POST_LIQUID_SATURATION)
    cs_field_set_key_int(mc->l_saturation, post_key, 1);

  /* Add default post-processing related to groundwater flow module */

  cs_post_add_time_mesh_dep_output(cs_gwf_extra_post_mtpf, gw);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the last setup step in the case of a (miscible) two-phase
 *        flow model in porous media
 *
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       quant      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_mtpf_finalize_setup(const cs_cdo_connect_t          *connect,
                     const cs_cdo_quantities_t       *quant,
                     cs_gwf_miscible_two_phase_t     *mc)
{
  CS_UNUSED(quant);

  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  const cs_lnum_t  n_cells = connect->n_cells;
  const size_t  csize = n_cells*sizeof(cs_real_t);

  /* Set the Darcian flux (in the volume and at the boundary)
   * Not done up to now
   */

  /* Allocate and initialize the relative permeability in the liquid and gas
     phase */

  BFT_MALLOC(mc->l_rel_permeability, n_cells, cs_real_t);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++)
    mc->l_rel_permeability[i] = 1; /* saturated by default */

  BFT_MALLOC(mc->g_rel_permeability, n_cells, cs_real_t);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++)
    mc->g_rel_permeability[i] = 1; /* saturated by default */

  /* Define the array storing the time property for the water eq. */

  BFT_MALLOC(mc->time_w_eq_array, n_cells, cs_real_t);
  memset(mc->time_w_eq_array, 0, csize);

  cs_property_def_by_array(mc->time_w_eq_pty,
                           cs_flag_primal_cell,  /* where data are located */
                           mc->time_w_eq_array,
                           false,                /* not owner of the array */
                           NULL);                /* no index */

  /* Define the array storing the diffusion property for the water eq. */

  BFT_MALLOC(mc->diff_w_eq_array, n_cells, cs_real_t);
  memset(mc->diff_w_eq_array, 0, csize);

  cs_property_def_by_array(mc->diff_w_eq_pty,
                           cs_flag_primal_cell,  /* where data are located */
                           mc->diff_w_eq_array,
                           false,                /* not owner of the array */
                           NULL);                /* no index */

  /* Define the array storing the time property for the hydrogen eq. */

  BFT_MALLOC(mc->time_h_eq_array, n_cells, cs_real_t);
  memset(mc->time_h_eq_array, 0, csize);

  cs_property_def_by_array(mc->time_h_eq_pty,
                           cs_flag_primal_cell,  /* where data are located */
                           mc->time_h_eq_array,
                           false,                /* not owner of the array */
                           NULL);                /* no index */

  /* Define the array storing the diffusion property for the liquid part in the
     hydrogen eq. */

  BFT_MALLOC(mc->diff_hl_eq_array, n_cells, cs_real_t);
  memset(mc->diff_hl_eq_array, 0, csize);

  cs_property_def_by_array(mc->diff_hl_eq_pty,
                           cs_flag_primal_cell,  /* where data are located */
                           mc->diff_hl_eq_array,
                           false,                /* not owner of the array */
                           NULL);                /* no index */

  /* Define the array storing the diffusion property for the gaz part in the
     hydrogen eq. */

  BFT_MALLOC(mc->diff_hg_eq_array, n_cells, cs_real_t);
  memset(mc->diff_hg_eq_array, 0, csize);

  cs_property_def_by_array(mc->diff_hg_eq_pty,
                           cs_flag_primal_cell,  /* where data are located */
                           mc->diff_hg_eq_array,
                           false,                /* not owner of the array */
                           NULL);                /* no index */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the pressure in the liquid phase
 *         Case of two-phase flows in porous media.
 *
 * \param[in, out] mc         pointer to a modelling context structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

/* static void */
/* _tpf_water_pressure(cs_gwf_miscible_two_phase_t          *mc, */
/*                     const cs_cdo_quantities_t   *cdoq, */
/*                     const cs_cdo_connect_t      *connect, */
/*                     bool                         cur2prev) */
/* { */
/*   cs_gwf_t  *gw = cs_gwf_main_structure; */
/*   assert(mc != NULL && gw != NULL); */
/*   const cs_equation_t  *w_eq = mc->w_eq; */
/*   const cs_equation_t  *h_eq = mc->h_eq; */

/*   cs_field_t  *capillary_pressure = cs_equation_get_field(w_eq); */
/*   cs_field_t  *gcomp_pressure = cs_equation_get_field(h_eq); */
/*   cs_field_t  *water_pressure = mc->water_pressure; */

/*   cs_param_space_scheme_t r_scheme = cs_equation_get_space_scheme(w_eq); */

/*   switch (r_scheme) { */

/*   case CS_SPACE_SCHEME_CDOVB: */
/*     { */
/*       assert(capillary_pressure->location_id == */
/*              cs_mesh_location_get_id_by_name("vertices")); */
/*       assert(gcomp_pressure->location_id == */
/*              cs_mesh_location_get_id_by_name("vertices")); */
/*       assert(water_pressure->location_id == */
/*              cs_mesh_location_get_id_by_name("vertices")); */

/* #     pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN) */
/*       for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) */
/*         water_pressure->val[i] = gcomp_pressure->val[i] - capillary_pressure->val[i]; */

/*     } */
/*     break; */

/*   case CS_SPACE_SCHEME_CDOFB: */
/*   default: */
/*     break; /\* Nothing to do *\/ */
/*   } */

/* } */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the update step in the case of a miscible two-phase flow
 *         model in porous media
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
_mtpf_updates(const cs_mesh_t                 *mesh,
              const cs_cdo_connect_t          *connect,
              const cs_cdo_quantities_t       *quant,
              const cs_time_step_t            *ts,
              bool                             cur2prev,
              cs_gwf_miscible_two_phase_t     *mc)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  cs_real_t  time_eval = ts->t_cur;
  if (cur2prev)
    time_eval = _get_time_eval(ts, mc->w_eq);

  cs_param_space_scheme_t space_scheme = cs_equation_get_space_scheme(mc->w_eq);

  const cs_real_t  *l_pr = mc->l_pressure->val;
  const cs_real_t  *g_pr = mc->g_pressure->val;
  cs_real_t  *c_pr = mc->c_pressure->val;

  /* Update the capillary pressure: P_c = P_g - P_l */

  if (cur2prev)
    cs_field_current_to_previous(mc->c_pressure);

  switch (space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    for (cs_lnum_t i = 0; i < quant->n_vertices; i++)
      c_pr[i] = g_pr[i] - l_pr[i];
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme", __func__);

  }

/*   /\* Update properties related to soils (water_saturation, water/gaz permeabilities, */
/*      compressibility coefficients *\/ */

/*   if (cur2prev) { */
/*     cs_field_current_to_previous(gw->permea_field); */
/*     cs_field_current_to_previous(mc->gcomp_permea_field); */
/*     cs_field_current_to_previous(mc->water_pressure); */
/*     cs_field_current_to_previous(mc->l_saturation); */
/*   } */

/*   const cs_equation_t  *w_eq = mc->w_eq; */
/*   const cs_field_t  *capillary_pressure = cs_equation_get_field(w_eq); */

/*   const int n_soils = cs_gwf_get_n_soils(); */
/*   for (int i = 0; i < n_soils; i++) { */

/*     cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(i); */
/*     const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id); */

/*     // Update soil properties in soil struct TO BE EDITED */
/*     soil->update_properties(time_eval, mesh, connect, quant, */
/*                             capillary_pressure->val, */
/*                             zone, */
/*                             soil->input); */

/*   } /\* Loop on soils *\/ */

/* /\* #if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1 *\/ */
/* /\*   cs_dbg_darray_to_listing("MOISTURE_CONTENT", *\/ */
/* /\*                            quant->n_cells, *\/ */
/* /\*                            mc->moisture_field->val, 8); *\/ */
/* /\* #endif *\/ */

  return time_eval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the groundwater flow module in case
 *         of (unsaturated) miscible two-phase flows in porous media
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc        pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_mtpf_extra_op(const cs_cdo_connect_t                *connect,
               const cs_cdo_quantities_t             *cdoq,
               cs_gwf_unsaturated_single_phase_t     *mc)
{
  CS_UNUSED(connect);
  CS_UNUSED(cdoq);

  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  if (cs_flag_test(gw->post_flag, CS_GWF_POST_DARCY_FLUX_BALANCE) == false)
    return; /* Nothing to do */

  /* TO BE DONE (not useful up to now) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for the groundwater flows module.
 *         Case of (miscible) two-phase flows in porous media.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_mtpf_compute(const cs_mesh_t                   *mesh,
              const cs_time_step_t              *time_step,
              const cs_cdo_connect_t            *connect,
              const cs_cdo_quantities_t         *cdoq,
              cs_gwf_miscible_two_phase_t       *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  assert(gw != NULL && mc != NULL);
  assert(gw->model == CS_GWF_MODEL_TWO_PHASE);

  cs_equation_t  *w_eq = mc->w_eq;
  cs_equation_t  *h_eq = mc->h_eq;

  assert(w_eq != NULL && h_eq != NULL);
  assert(cs_equation_get_type(w_eq) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(cs_equation_get_type(h_eq) == CS_EQUATION_TYPE_GROUNDWATER);

  /* TODO */

  bool cur2prev = true;

  /* Build and solve the linear system related to the Richards equations */

  if (!cs_equation_is_steady(w_eq) || !cs_equation_is_steady(h_eq)) {

    /* Solve the algebraic system. By default, a current to previous operation
       is performed */

    /* (@_@) Not clear how to do this: to check with J. Bonelle for this part */

    cs_equation_solve(cur2prev, mesh, w_eq);
    cs_equation_solve(cur2prev, mesh, h_eq);

    /* Update the variables related to the groundwater flow system */

    cs_gwf_update(mesh, connect, cdoq, time_step, cur2prev);

  }
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

  /* Property common to all model */

  gw->abs_permeability = NULL;

  /* Modelling context */

  gw->model_context = NULL;

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
 * \param[in]   pty_type        type of permeability (iso, ortho...)
 * \param[in]   model           type of physical modelling
 * \param[in]   option_flag     optional flag to specify this module
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

  /* Store the pointer to the groundawater flow structure */

  cs_gwf_main_structure = gw;

  gw->model = model;
  gw->flag = option_flag;

  /* Add the absolute permeability property (used in the definition of the
   * diffusion term in the conservation equations and in the definition of the
   * Darcy flux. According to the type of soil model, the full permeability is
   * weigthed by a relative permeability (e.g. in a Van Genuchten-Mualen
   * model). */

  gw->abs_permeability = cs_property_add("absolute_permeability", pty_type);

  /* Allocate and initialize each model context (mc) */

  switch (model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    gw->model_context = _sspf_init_context();
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    gw->model_context = _uspf_init_context(pty_type);
    break;

  case CS_GWF_MODEL_TWO_PHASE:
    gw->post_flag |= CS_GWF_POST_LIQUID_SATURATION;
    gw->model_context = _mtpf_init_context(pty_type);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);

  }

  return gw;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to groundwater flows
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

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    {
      cs_gwf_saturated_single_phase_t  *mc = gw->model_context;

    _sspf_free_context(&(mc));
    }
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    {
      cs_gwf_unsaturated_single_phase_t  *mc = gw->model_context;

      _uspf_free_context(&(mc));
    }
    break;

  case CS_GWF_MODEL_TWO_PHASE:
    {
      cs_gwf_miscible_two_phase_t  *mc = gw->model_context;

      _mtpf_free_context(&(mc));
    }
    break;

  default:
    /* Nothing to do */
    break;
  }

  /* Free all soils */

  cs_gwf_soil_free_all();

  /* Free all tracers */

  cs_gwf_tracer_free_all();

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

  bool  post_capacity =
    (gw->post_flag & CS_GWF_POST_SOIL_CAPACITY) ? true : false;
  bool  post_liquid_saturation =
    (gw->post_flag & CS_GWF_POST_LIQUID_SATURATION) ? true : false;
  bool  post_permeability =
    (gw->post_flag & CS_GWF_POST_PERMEABILITY) ? true : false;

  cs_log_printf(CS_LOG_SETUP, "  * GWF | Post:"
                " Soil capacity %s Liquid saturation %s Permeability %s\n",
                cs_base_strtf(post_capacity),
                cs_base_strtf(post_liquid_saturation),
                cs_base_strtf(post_permeability));

  bool  do_balance =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_BALANCE) ? true : false;
  bool  do_divergence =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_DIVERGENCE) ? true : false;
  bool  post_boundary =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY) ? true : false;

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Darcy Flux:"
                " Balance %s Divergence %s At boundary faces: %s\n",
                cs_base_strtf(do_balance),
                cs_base_strtf(do_divergence),
                cs_base_strtf(post_boundary));

  /* Main options */

  switch(gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Model: %s\n", cs_gwf_model_name[gw->model]);
    _sspf_log_context(gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Model: %s\n", cs_gwf_model_name[gw->model]);
    _uspf_log_context(gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE:
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Model: %s\n", cs_gwf_model_name[gw->model]);
    _mtpf_log_context(gw->model_context);
    break;

  default:
    break;

  }

  /* Soils */

  cs_gwf_soil_log_setup();

  /* Tracers */

  cs_gwf_tracer_log_all();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the parameters defining the two-phase flow model.
 *         Use SI unit if not prescribed otherwise.
 *
 * \param[in] l_mass_density   mass density of the main liquid component
 * \param[in] l_viscosity      viscosity in the liquid phase (Pa.s)
 * \param[in] g_viscosity      viscosity in the gas phase (Pa.s)
 * \param[in] l_diffusivity_h  diffusivity of the main gas component in the
 *                             liquid phase
 * \param[in] w_molar_mass     molar mass of the main liquid component
 * \param[in] h_molar_mass     molar mass of the main gas component
 * \param[in] ref_temperature  reference temperature in Kelvin
 * \param[in] henry_constant   constant in the Henry law
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_two_phase_model(cs_real_t       l_mass_density,
                           cs_real_t       l_viscosity,
                           cs_real_t       g_viscosity,
                           cs_real_t       l_diffusivity_h,
                           cs_real_t       w_molar_mass,
                           cs_real_t       h_molar_mass,
                           cs_real_t       ref_temperature,
                           cs_real_t       henry_constant)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));
  if (gw->model != CS_GWF_MODEL_TWO_PHASE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid model. One expects a two-phase flow model.\n",
              __func__);

  cs_gwf_miscible_two_phase_t  *mc = gw->model_context;

  assert(mc != NULL);
  assert(l_mass_density > 0);
  assert(ref_temperature > 0);  /* In Kelvin */
  assert(w_molar_mass > 0 && h_molar_mass > 0);
  assert(l_viscosity > 0 && g_viscosity > 0);

  /* Set the parameters */

  mc->l_mass_density = l_mass_density;
  mc->l_viscosity = l_viscosity;
  mc->g_viscosity = g_viscosity;
  mc->l_diffusivity_h = l_diffusivity_h;
  mc->w_molar_mass = w_molar_mass;
  mc->h_molar_mass = h_molar_mass;
  mc->ref_temperature = ref_temperature;
  mc->henry_constant = henry_constant;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the flag dedicated to the post-processing of the GWF module
 *
 * \param[in]  post_flag             flag to set
 * \param[in]  reset                 reset post flag before
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_post_options(cs_flag_t       post_flag,
                        bool            reset)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  if (gw == NULL)
    return;

  if (reset)
    gw->post_flag = post_flag;
  else
    gw->post_flag |= post_flag;

  if (gw->post_flag & CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY) {

    cs_adv_field_t  *adv = _get_l_adv_field(gw);
    if (adv != NULL)
      adv->status |= CS_ADVECTION_FIELD_DEFINE_AT_BOUNDARY_FACES;

    if (gw->model == CS_GWF_MODEL_TWO_PHASE) {

      cs_gwf_miscible_two_phase_t  *mc =
        (cs_gwf_miscible_two_phase_t *)gw->model_context;

      if (mc->g_darcy != NULL) {
        adv = mc->g_darcy->adv_field;
        if (adv != NULL)
          adv->status |= CS_ADVECTION_FIELD_DEFINE_AT_BOUNDARY_FACES;
      }

    }

  } /* Darcy flux at boundary */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the advection field related to the Darcy flux in the liquid
 *         phase
 *
 * \return a pointer to a cs_adv_field_t structure or NULL
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_gwf_get_adv_field(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL)
    return NULL;

  return  _get_l_adv_field(gw);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a new cs_gwf_soil_t structure. An initialization by
 *         default of all members is performed.
 *
 * \param[in]  z_name        name of the volume zone corresponding to the soil
 * \param[in]  bulk_density  value of the mass density
 * \param[in]  sat_moisture  value of the saturated moisture content
 * \param[in]  model         type of model for the soil behavior
 *
 * \return a pointer to the new allocated soil structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_add_soil(const char                *z_name,
                double                     bulk_density,
                double                     sat_moisture,
                cs_gwf_soil_model_t        model)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  const cs_zone_t  *zone = cs_volume_zone_by_name_try(z_name);

  if (zone == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Zone %s related to the same soil is not defined.\n"
              " Stop adding a new soil.", z_name);

  assert(bulk_density > 0);
  assert(sat_moisture > 0);

  cs_property_type_t  perm_type = cs_property_get_type(gw->abs_permeability);

  cs_gwf_soil_t  *soil = cs_gwf_soil_create(zone,
                                            gw->model, /* hydraulic model */
                                            model,     /* soil model */
                                            perm_type,
                                            sat_moisture,
                                            bulk_density,
                                            gw->model_context);

  return soil;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module

 *         This equation is a particular type of unsteady advection-diffusion
 *         reaction equation. Tracer is advected thanks to the darcian velocity
 *         and diffusion/reaction parameters result from a physical modelling.
 *         Terms solved in this equation are activated according to predefined
 *         settings. The advection field corresponds to that of the liquid
 *         phase.
 *
 * \param[in]  tr_model   physical modelling to consider (0 = default settings)
 * \param[in]  eq_name    name of the tracer equation
 * \param[in]  var_name   name of the related variable
 *
 * \return a pointer to the new cs_gwf_tracer_t structure
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

  cs_adv_field_t  *adv = _get_l_adv_field(gw);

  /* Set the function pointers */

  cs_gwf_tracer_setup_t  *setup = NULL;

  if (gw->model == CS_GWF_MODEL_SATURATED_SINGLE_PHASE)
    setup = cs_gwf_tracer_saturated_setup;
  else
    setup = cs_gwf_tracer_unsaturated_setup;

  /* Call the main function to add a new tracer */

  cs_gwf_tracer_t  *tracer = cs_gwf_tracer_add(tr_model,
                                               gw->model,
                                               eq_name,
                                               var_name,
                                               adv,
                                               setup,
                                               cs_gwf_tracer_add_default_terms);

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *
 *         This equation is a particular type of unsteady advection-diffusion
 *         reaction equation.  Tracer is advected thanks to the darcian
 *         velocity and diffusion/reaction parameters result from a physical
 *         modelling. Terms are activated according to predefined settings.
 *         Modelling of the tracer parameters are left to the user
 *
 * \param[in]   eq_name     name of the tracer equation
 * \param[in]   var_name    name of the related variable
 * \param[in]   setup       function pointer (predefined prototype)
 * \param[in]   add_terms   function pointer (predefined prototype)
 *
 * \return a pointer to the new cs_gwf_tracer_t structure
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

  cs_adv_field_t  *adv = _get_l_adv_field(gw);

  /* Call the main function to add a new tracer */

  cs_gwf_tracer_t  *tracer = cs_gwf_tracer_add(CS_GWF_TRACER_USER,
                                               gw->model,
                                               eq_name,
                                               var_name,
                                               adv,
                                               setup,
                                               add_terms);

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for the groundwater flow model and its related
 *         equations.
 *         At this stage, all soils have been defined and equation parameters
 *         are set. Create new cs_field_t structures according to the setting.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_setup(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));
  cs_gwf_soil_check();

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    _sspf_init_setup(gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    _uspf_init_setup(gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE:
    _mtpf_init_setup(gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last initialization step of the groundwater flow module. At this
 *        stage, the mesh quantities are defined.
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

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    _sspf_finalize_setup(connect, quant, gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    _uspf_finalize_setup(connect, quant, gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE:
    _mtpf_finalize_setup(connect, quant, gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Store the soil id for each cell */

  cs_gwf_soil_check();
  cs_gwf_build_cell2soil(quant->n_cells);

  /* Finalize the tracer setup */

  cs_gwf_tracer_setup_all(connect, quant);
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

  if (gw == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Groundwater module is not allocated.", __func__);

  cs_real_t  time_eval = ts->t_cur; /* default value */

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    time_eval = _sspf_updates(connect, quant, ts, cur2prev, gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    time_eval = _uspf_updates(mesh, connect, quant, ts, cur2prev,
                              gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE:
    time_eval = _mtpf_updates(mesh, connect, quant, ts, cur2prev,
                              gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Update the diffusivity tensor associated to each tracer equation since the
     Darcy velocity may have changed */

  cs_gwf_tracer_update_diff_tensor(time_eval, mesh, connect, quant);
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

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    {
      cs_gwf_saturated_single_phase_t  *mc = gw->model_context;

      _spf_compute_steady_state(mesh, time_step, connect, cdoq, mc->richards);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for steady-state computations with the"
              " GroundWater Flow module.\n",
              __func__);
  }

  /* Compute tracers */
  /* --------------- */

  cs_gwf_tracer_compute_steady_all(mesh, time_step, connect, cdoq);
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

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    {
      cs_gwf_saturated_single_phase_t  *mc = gw->model_context;

      _spf_compute(mesh, time_step, connect, cdoq, mc->richards);
    }
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    {
      cs_gwf_unsaturated_single_phase_t  *mc = gw->model_context;

      _spf_compute(mesh, time_step, connect, cdoq, mc->richards);
    }
    break;

  case CS_GWF_MODEL_TWO_PHASE:
    _mtpf_compute(mesh, time_step, connect, cdoq, gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Compute tracers */
  /* --------------- */

  cs_gwf_tracer_compute_all(mesh, time_step, connect, cdoq);
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

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    _sspf_extra_op(connect, cdoq, gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    _uspf_extra_op(connect, cdoq, gw->model_context);
    break;

  case CS_GWF_MODEL_TWO_PHASE:
    _mtpf_extra_op(connect, cdoq, gw->model_context);
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
 *         in case of saturated single-phase flows (sspf) in porous media.
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
cs_gwf_extra_post_sspf(void                   *input,
                       int                     mesh_id,
                       int                     cat_id,
                       int                     ent_flag[5],
                       cs_lnum_t               n_cells,
                       cs_lnum_t               n_i_faces,
                       cs_lnum_t               n_b_faces,
                       const cs_lnum_t         cell_ids[],
                       const cs_lnum_t         i_face_ids[],
                       const cs_lnum_t         b_face_ids[],
                       const cs_time_step_t   *time_step)
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

  if (mesh_id == CS_POST_MESH_VOLUME) {

    const cs_gwf_saturated_single_phase_t  *mc = gw->model_context;

    if (gw->post_flag & CS_GWF_POST_DARCY_FLUX_DIVERGENCE) {

      cs_adv_field_t *adv = _get_l_adv_field(gw);

      /* Only case avalaible up to now */

      if (cs_advection_field_get_deftype(adv) == CS_XDEF_BY_ARRAY) {

        cs_real_t  *divergence =
          cs_advection_field_divergence_at_vertices(adv, time_step->t_cur);

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

    if (gw->post_flag & CS_GWF_POST_LIQUID_SATURATION) {

      cs_real_t  *liquid_saturation = NULL;
      BFT_MALLOC(liquid_saturation, n_cells, cs_real_t);

      cs_property_eval_at_cells(time_step->t_cur,
                                mc->moisture_content,
                                liquid_saturation);

      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_DEFAULT,
                        "liquid_saturation",
                        1,
                        false,
                        false,
                        CS_POST_TYPE_cs_real_t,
                        liquid_saturation,
                        NULL,
                        NULL,
                        time_step);

      BFT_FREE(liquid_saturation);

    } /* Post-processing of the liquid saturation */

    if (gw->post_flag & CS_GWF_POST_PERMEABILITY) {

      int  dim = cs_property_get_dim(gw->abs_permeability);
      cs_real_t  *permeability = NULL;
      BFT_MALLOC(permeability, dim*n_cells, cs_real_t);

      cs_property_eval_at_cells(time_step->t_cur,
                                gw->abs_permeability,
                                permeability);

      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_DEFAULT,
                        "permeability",
                        dim,
                        true,
                        false,
                        CS_POST_TYPE_cs_real_t,
                        permeability,
                        NULL,
                        NULL,
                        time_step);

      BFT_FREE(permeability);

    } /* Post-processing of the permeability field */

  } /* volume mesh id */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the groundwater flow module
 *         in case of unsaturated single-phase flows (uspf) in porous media.
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
cs_gwf_extra_post_uspf(void                   *input,
                       int                     mesh_id,
                       int                     cat_id,
                       int                     ent_flag[5],
                       cs_lnum_t               n_cells,
                       cs_lnum_t               n_i_faces,
                       cs_lnum_t               n_b_faces,
                       const cs_lnum_t         cell_ids[],
                       const cs_lnum_t         i_face_ids[],
                       const cs_lnum_t         b_face_ids[],
                       const cs_time_step_t   *time_step)
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

  if (mesh_id == CS_POST_MESH_VOLUME) {
    if (gw->post_flag & CS_GWF_POST_DARCY_FLUX_DIVERGENCE) {

      cs_adv_field_t *adv = _get_l_adv_field(gw);

      /* Only case avalaible up to now */

      if (cs_advection_field_get_deftype(adv) == CS_XDEF_BY_ARRAY) {

        cs_real_t  *divergence =
          cs_advection_field_divergence_at_vertices(adv, time_step->t_cur);

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
/*!
 * \brief  Predefined post-processing output for the groundwater flow module
 *         in case of miscible two-phase flows (mtpf) in porous media.
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
cs_gwf_extra_post_mtpf(void                      *input,
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
  const cs_gwf_miscible_two_phase_t  *mc =
    (const cs_gwf_miscible_two_phase_t  *)gw->model_context;

  if (mesh_id == CS_POST_MESH_VOLUME) {

    if (gw->post_flag & CS_GWF_POST_DARCY_FLUX_DIVERGENCE) {

      cs_adv_field_t *adv = _get_l_adv_field(gw);

      /* Only case avalaible up to now */

      if (cs_advection_field_get_deftype(adv) == CS_XDEF_BY_ARRAY) {

        cs_real_t  *divergence =
          cs_advection_field_divergence_at_vertices(adv, time_step->t_cur);

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

    if (gw->post_flag & CS_GWF_POST_PERMEABILITY) {

      int  dim = cs_property_get_dim(gw->abs_permeability);
      cs_real_t  *permeability = NULL;
      BFT_MALLOC(permeability, dim*n_cells, cs_real_t);

      cs_property_eval_at_cells(time_step->t_cur,
                                gw->abs_permeability,
                                permeability);

      /* permeability = abs_permeability * l_rel_permeability */

      assert(mc->l_rel_permeability != NULL);
      for (cs_lnum_t c = 0; c < n_cells; c++)
        for (int k = 0; k < dim; k++)
          permeability[dim*c+k] *= mc->l_rel_permeability[c];

      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_DEFAULT,
                        "permeability",
                        dim,
                        true,
                        false,
                        CS_POST_TYPE_cs_real_t,
                        permeability,
                        NULL,
                        NULL,
                        time_step);

      BFT_FREE(permeability);

    } /* Post-processing of the permeability field */

  } /* volume mesh id */
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
