/*============================================================================
 * Main functions dedicated to groundwater flows
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

#include "cs_boundary_zone.h"
#include "cs_cdo_blas.h"
#include "cs_cdovb_scaleq.h"
#include "cs_cdovb_priv.h"
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
    N_("Miscible two-phase model (capillary/gas pressure)"),
    N_("Immiscible two-phase model (capillary/gas pressure)")
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

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    {
      cs_gwf_two_phase_t  *mc = gw->model_context;

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

    cs_gwf_update(mesh, connect, cdoq, time_step, CS_FLAG_CURRENT_TO_PREVIOUS);

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

    cs_gwf_update(mesh, connect, cdoq, time_step, CS_FLAG_CURRENT_TO_PREVIOUS);

  }
  else {

    /* Richards is steady but one can force the resolution */

    if (gw->flag & CS_GWF_FORCE_RICHARDS_ITERATIONS) {

      /* Solve the algebraic system */

      cs_equation_solve_steady_state(mesh, richards);

      /* Update the variables related to the groundwater flow system */

      cs_gwf_update(mesh, connect, cdoq, time_step,
                    CS_FLAG_CURRENT_TO_PREVIOUS);

    }

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
_sspf_activate(void)
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

  /* Add a property related to the moisture content (It should be a constant
     in each soil) */

  mc->moisture_content = cs_property_add("moisture_content", CS_PROPERTY_ISO);

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
 * \brief Initialize the model context according to the settings done inside
 *        the function cs_user_model()
 *        Case of a saturated single-phase flows model in porous media
 *
 * \param[in, out] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_sspf_init_model_context(cs_gwf_saturated_single_phase_t   *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (mc == NULL)
    return;

  cs_equation_t  *eq = mc->richards;

  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The Richards equation is not defined. Stop execution.\n",
              __func__);
  assert(cs_equation_is_steady(eq) == true);

  cs_equation_param_t  *eqp = cs_equation_get_param(eq);
  assert(eqp != NULL);

  /* Add the diffusion term to the Richards equation by associating the
     absolute permeability to the diffusion property of the Richards eq. */

  cs_equation_add_diffusion(eqp, gw->abs_permeability);

  /* Add the variable field */

  if (gw->flag & CS_GWF_FORCE_RICHARDS_ITERATIONS)
    cs_equation_predefined_create_field(1, eq); /* Keep two states */
  else
    cs_equation_predefined_create_field(0, eq); /* Keep only one state */
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
  if (mc == NULL)
    return;

  cs_gwf_t  *gw = cs_gwf_main_structure;

  cs_equation_t  *eq = mc->richards;
  assert(eq != NULL);
  assert(cs_equation_is_steady(eq) == true);

  cs_equation_param_t  *eqp = cs_equation_get_param(eq);
  assert(eqp != NULL);

  /* Set the "has_previous" flag */

  bool  has_previous = false;
  if (gw->flag & CS_GWF_FORCE_RICHARDS_ITERATIONS)
    has_previous = true;

  /* Add new fields if needed */

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");

  /* Handle gravity effects */

  if (gw->flag & CS_GWF_GRAVITATION) {

    switch (eqp->space_scheme) {
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
  if (mc == NULL)
    return;

  const cs_param_space_scheme_t  richards_scheme =
    cs_equation_get_space_scheme(mc->richards);

  /* Set the Darcian flux (in the volume and at the boundary) */

  cs_gwf_darcy_flux_define(connect, quant, richards_scheme, mc->darcy);

  /* Set the moisture content from the soil porosity */

  cs_gwf_soil_saturated_set_property(mc->moisture_content);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the update step in the case of a saturated single-phase
 *         flow model in porous media
 *
 * \param[in]       connect      pointer to a cs_cdo_connect_t structure
 * \param[in]       quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]       ts           pointer to a cs_time_step_t structure
 * \param[in]       update_flag  metadata associated to type of operation to do
 * \param[in, out]  mc           pointer to the casted model context
 *
 *\return the value at which one has to evaluate the properties
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_sspf_updates(const cs_cdo_connect_t            *connect,
              const cs_cdo_quantities_t         *quant,
              const cs_time_step_t              *ts,
              cs_flag_t                          update_flag,
              cs_gwf_saturated_single_phase_t   *mc)
{
  cs_real_t  time_eval = ts->t_cur;
  bool  cur2prev = false;

  if (update_flag & CS_FLAG_CURRENT_TO_PREVIOUS) {

    time_eval = _get_time_eval(ts, mc->richards);
    cur2prev = true;

  }

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
 * \return a pointer to a new allocated cs_gwf_unsaturated_single_phase_t struct
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_unsaturated_single_phase_t *
_uspf_activate(void)
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
   * Add an advection field related to the darcian flux stemming from the
   * Richards equation. This advection field is steady since the head is
   * steady in this model.
   */

  mc->darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);

  cs_advection_field_status_t  adv_status = CS_ADVECTION_FIELD_GWF |
    CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;

  mc->darcy->adv_field = cs_advection_field_add("darcy_field", adv_status);

  /* Add several properties:
   * - the moisture content
   * - the soil capacity --> for the unsteady term
   * - the (full) permeability = rel_permeability * abs_permeability
   */

  mc->moisture_content = cs_property_add("moisture_content", CS_PROPERTY_ISO);

  mc->soil_capacity = cs_property_add("soil_capacity", CS_PROPERTY_ISO);

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
 * \brief Initialize the model context according to the settings done inside
 *        the function cs_user_model()
 *        Case of an unsaturated single-phase flows model in porous media
 *
 * \param[in, out] mc          pointer to the model context structure
 * \param[in, out] perm_type   type of permeability to handle
 */
/*----------------------------------------------------------------------------*/

static void
_uspf_init_model_context(cs_gwf_unsaturated_single_phase_t   *mc,
                         cs_property_type_t                   perm_type)
{
  if (mc == NULL)
    return;

  /* Add the property related to the diffusion term */

  mc->permeability = cs_property_add("permeability", perm_type);

  /* Define the Richards equation */

  if (mc->richards == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The Richards equation is not defined. Stop execution.\n",
              __func__);

  cs_equation_param_t  *eqp = cs_equation_get_param(mc->richards);
  assert(eqp != NULL);

  /* Add the diffusion term to the Richards equation by associating the
     full permeability to the diffusion property of the Richards eq. */

  cs_equation_add_diffusion(eqp, mc->permeability);

  /* Add the time term to the Richards equation by associating the the soil
     capacity to the unsteady term */

  assert(mc->soil_capacity != NULL);
  cs_equation_add_time(eqp, mc->soil_capacity);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup related to the modelling context of unsaturated single-phase
 *        flows. At this stage, all soils have been defined and equation
 *        parameters are set.
 *
 * \param[in, out] mc          pointer to the model context structure
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

  cs_equation_param_t  *eqp = cs_equation_get_param(mc->richards);
  assert(eqp != NULL);

  /* Add the variable field (Keep one previous state) */

  cs_equation_predefined_create_field(1, mc->richards);

  /* Set additional fields */

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");

  /* Handle gravity effects */

  if (gw->flag & CS_GWF_GRAVITATION) {

    switch (eqp->space_scheme) {
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
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      update_flag  metadata associated to type of operation to do
 * \param[in, out] mc           pointer to the casted model context
 *
 *\return the value at which one has to evaluate the properties
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_uspf_updates(const cs_mesh_t                    *mesh,
              const cs_cdo_connect_t             *connect,
              const cs_cdo_quantities_t          *quant,
              const cs_time_step_t               *ts,
              cs_flag_t                           update_flag,
              cs_gwf_unsaturated_single_phase_t  *mc)
{
  cs_real_t  time_eval = ts->t_cur;
  bool cur2prev = false;

  if (update_flag & CS_FLAG_CURRENT_TO_PREVIOUS) {

    time_eval = _get_time_eval(ts, mc->richards);
    cur2prev = true;

  }

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
 * Functions related to the model of two phase flows (TPF) miscible or
 * immiscible
 * ========================================================================== */

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add a source term to the hg equation when a segregated approach
 *          is used. Case of a vertex-based scheme.
 *
 *          Generic function prototype for a hook during the cellwise building
 *          of the linear system
 *          Fit the cs_equation_build_hook_t prototype. This function may be
 *          called by different OpenMP threads
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
_build_hg_incr_st(const cs_equation_param_t     *eqp,
                  const cs_equation_builder_t   *eqb,
                  const void                    *eq_context,
                  const cs_cell_mesh_t          *cm,
                  void                          *context,
                  cs_hodge_t                    *mass_hodge,
                  cs_hodge_t                    *diff_hodge,
                  cs_cell_sys_t                 *csys,
                  cs_cell_builder_t             *cb)
{
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);
  CS_UNUSED(mass_hodge);

  const cs_cdovb_scaleq_t  *eqc = (const cs_cdovb_scaleq_t *)eq_context;
  assert(csys->n_dofs == cm->n_vc);

  cs_gwf_two_phase_t  *mc = context;

  /* Modify the property data associated to the hodge operator related to the
     diffusion term */

  cs_property_data_t  *saved = diff_hodge->pty_data;
  cs_property_data_t  tmp = cs_property_data_define(saved->need_tensor,
                                                    saved->need_eigen,
                                                    mc->diff_hl_pty);

  diff_hodge->pty_data = &tmp;

  /* Diffusion part of the source term to add */

  cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, cb->cell_flag, diff_hodge);

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
 * \brief  Allocate and initialize the modelling context for the model of
 *         (unsaturated) miscible two-phase flows
 *
 * \return a pointer to a new allocated cs_gwf_two_phase_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_two_phase_t *
_tpf_activate(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL);

  cs_gwf_two_phase_t  *mc = NULL;

  BFT_MALLOC(mc, 1, cs_gwf_two_phase_t);

  /* Define the system of equations */
  /* ------------------------------ */

  /* Create a new equation for the water conservation associated to the
     pressure in the liquid phase. This will stand for the (0,0)-block */

  mc->wl_eq = cs_equation_add("w_conservation",   /* equation name */
                              "liquid_pressure",  /* variable name */
                              CS_EQUATION_TYPE_GROUNDWATER,
                              1,
                              CS_PARAM_BC_HMG_NEUMANN);

  /* Create a new equation for the hydrogen conservation associated to the
     pressure in the gaseous phase. This will stand for the (1,1)-block */

  mc->hg_eq = cs_equation_add("h_conservation", /* equation name */
                              "gas_pressure",   /* variable name */
                              CS_EQUATION_TYPE_GROUNDWATER,
                              1,
                              CS_PARAM_BC_HMG_NEUMANN);

  /* Set to NULL the remaining pointers. Only a part is set according to the
     numerical settings */

  mc->wg_eqp = NULL;
  mc->hl_eqp = NULL;
  mc->system = NULL;

  /* Advection fields */
  /* ---------------- */

  /* Darcy flux (not used up to now) */

  mc->l_darcy = NULL;
  mc->g_darcy = NULL;

  /* Properties
   * ----------
   *
   * Adding the properties related to the permeability (diffusion) is
   * postponed since one has to know which type of permeability is considered
   * (iso, ortho, or anisotropic) and this information is a result of the type
   * of soils which have been added.
   */

  /* Properties related to the water equation  */

  mc->time_wl_pty = NULL;
  mc->diff_wl_pty = NULL;
  mc->time_wg_pty = NULL;

  /* Properties related to the hydrogen equation  */

  mc->time_hg_pty = NULL;
  mc->diff_hg_pty = NULL;
  mc->reac_hg_pty = NULL;
  mc->time_hl_pty = NULL;
  mc->diff_hl_pty = NULL;

  /* Fields */
  /* ------ */

  mc->l_saturation = NULL;
  mc->c_pressure = NULL;
  mc->l_pressure = NULL;
  mc->g_pressure = NULL;

  /* Arrays */
  /* ------ */

  /* The properties will be defined using arrays.
   * Store these arrays of property values associated to equation terms */

  mc->time_wl_array = NULL;
  mc->diff_wl_array = NULL;
  mc->srct_wl_array = NULL;

  mc->time_wg_array = NULL;

  mc->time_hg_array = NULL;
  mc->diff_hg_array = NULL;
  mc->reac_hg_array = NULL;
  mc->srct_hg_array = NULL;

  mc->time_hl_array = NULL;
  mc->diff_hl_array = NULL;

  /* Array of additional variable or property values */

  mc->l_rel_permeability = NULL;
  mc->g_rel_permeability = NULL;
  mc->l_capacity = NULL;
  mc->capillarity_cell_pressure = NULL;
  mc->g_cell_pressure = NULL;
  mc->l_saturation_submesh = NULL;
  mc->l_saturation_submesh_pre = NULL;

  /* Model parameters (default values) */
  /* ---------------- */

  mc->l_mass_density = 1000;
  mc->l_viscosity = 1e-3;
  mc->g_viscosity = 2e-5;
  mc->l_diffusivity_h = 0;      /* immiscible case */
  mc->w_molar_mass = 18e-3;
  mc->h_molar_mass = 3e-3;
  mc->ref_temperature = 280;    /* in Kelvin */
  mc->henry_constant = 1e-20;   /* nearly immiscible case */

  /* Numerical parameters */

  mc->use_coupled_solver = false;
  mc->use_incremental_solver = true;
  mc->use_properties_on_submesh = true;

  mc->nl_algo_type = CS_PARAM_NL_ALGO_NONE; /* Linear algo. by default */

  mc->nl_algo_param.n_max_algo_iter = 50;
  mc->nl_algo_param.rtol = 1e-5;
  mc->nl_algo_param.atol = 1e-5;
  mc->nl_algo_param.dtol = 1e3;
  mc->nl_algo_param.verbosity = 1;

  mc->anderson_param.n_max_dir = 5;
  mc->anderson_param.starting_iter = 3;
  mc->anderson_param.max_cond = -1; /* No test by default */
  mc->anderson_param.beta = 1.0;    /* No damping by default */
  mc->anderson_param.dp_type = CS_PARAM_DOTPROD_EUCLIDEAN;

  mc->nl_algo = NULL;

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
_tpf_free_context(cs_gwf_two_phase_t  **p_mc)
{
  if (p_mc == NULL)
    return;
  if (*p_mc == NULL)
    return;

  cs_gwf_two_phase_t  *mc = *p_mc;

  /* System of equations are freed elsewhere (just after having freed
     cs_equation_t structures) */

  cs_gwf_darcy_flux_free(&(mc->l_darcy));
  cs_gwf_darcy_flux_free(&(mc->g_darcy));

  BFT_FREE(mc->time_wl_array);
  BFT_FREE(mc->diff_wl_array);
  BFT_FREE(mc->srct_wl_array);
  BFT_FREE(mc->time_wg_array);
  BFT_FREE(mc->time_hg_array);
  BFT_FREE(mc->diff_hg_array);
  BFT_FREE(mc->reac_hg_array);
  BFT_FREE(mc->srct_hg_array);
  BFT_FREE(mc->time_hl_array);
  BFT_FREE(mc->diff_hl_array);

  BFT_FREE(mc->l_rel_permeability);
  BFT_FREE(mc->g_rel_permeability);
  BFT_FREE(mc->l_capacity);
  BFT_FREE(mc->capillarity_cell_pressure);
  BFT_FREE(mc->g_cell_pressure);
  BFT_FREE(mc->l_saturation_submesh);
  BFT_FREE(mc->l_saturation_submesh_pre);

  if (mc->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

    /* Non-linearity: If the context is not NULL, this means that an Anderson
       algorithm has been activated otherwise nothing to do */

    cs_iter_algo_aa_free(mc->nl_algo);

    BFT_FREE(mc->nl_algo);

  }

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
_mtpf_log_context(cs_gwf_two_phase_t   *mc)
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
 * \brief Log the setup related to the modelling context of immiscible
 *        two-phase flows
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_itpf_log_context(cs_gwf_two_phase_t   *mc)
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
 * \brief Log the setup related to the modelling context of two-phase flows
 *        Common to the different models relying on two-phase flows
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_log_context(cs_gwf_two_phase_t   *mc)
{
  if (mc == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "  * GWF | Reference temperature: %5.2f K\n",
                mc->ref_temperature);

  cs_gwf_darcy_flux_log(mc->l_darcy);
  cs_gwf_darcy_flux_log(mc->g_darcy);

  if (mc->use_coupled_solver)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Coupled solver\n");
  else
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Segregated solver\n");

  if (mc->use_incremental_solver)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Incremental solver\n");

  if (mc->use_properties_on_submesh)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Liquid saturation on submesh\n");

  if (mc->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

    cs_log_printf(CS_LOG_SETUP, "  * GWF | Non-linear algo.: %s\n",
                  cs_param_get_nl_algo_name(mc->nl_algo_type));

    cs_log_printf(CS_LOG_SETUP, "  * GWF | Tolerances of non-linear algo:"
                  " rtol: %5.3e; atol: %5.3e; dtol: %5.3e\n",
                  mc->nl_algo_param.rtol, mc->nl_algo_param.atol,
                  mc->nl_algo_param.dtol);
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Max of non-linear iterations: %d\n",
                  mc->nl_algo_param.n_max_algo_iter);

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
 *        Case of a miscible or immiscible two-phase flows model in porous media
 *
 * \param[in, out]  mc         pointer to the casted model context
 * \param[in, out]  perm_type  type of permeability to handle
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_init_model_context(cs_gwf_two_phase_t     *mc,
                        cs_property_type_t      perm_type)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  if (mc->wl_eq == NULL || mc->hg_eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Equations are not defined for this model. Stop execution.\n",
              __func__);

  cs_equation_param_t  *wl_eqp = cs_equation_get_param(mc->wl_eq);
  cs_equation_param_t  *hg_eqp = cs_equation_get_param(mc->hg_eq);

  if (mc->use_coupled_solver) {

    /* Define the coupled system of equations */
    /* -------------------------------------- */

    /* Create the (0,1)-block related to the water in the gas phase */

    mc->wg_eqp = cs_equation_param_create("water_gas_block",
                                          CS_EQUATION_TYPE_GROUNDWATER,
                                          1,
                                          CS_PARAM_BC_HMG_NEUMANN);

    /* Create the (1,0)-block related to the hydrogen in the liquid phase */

    mc->hl_eqp = cs_equation_param_create("h_liquid_block",
                                          CS_EQUATION_TYPE_GROUNDWATER,
                                          1,
                                          CS_PARAM_BC_HMG_NEUMANN);

    /* Add a 2x2 system of coupled equations and define each block */

    mc->system = cs_equation_system_add("PorousTwoPhaseFlow",
                                        2,   /* system size */
                                        1);  /* scalar-valued block */

    /* Set all the blocks in the coupled system */

    cs_equation_system_assign_equation(0, mc->wl_eq, mc->system);  /* (0,0) */
    cs_equation_system_assign_equation(1, mc->hg_eq, mc->system);  /* (1,1) */
    cs_equation_system_assign_param(0, 1, mc->wg_eqp, mc->system); /* (0,1) */
    cs_equation_system_assign_param(1, 0, mc->hl_eqp, mc->system); /* (1,0) */

    /* Properties */
    /* ---------- */

    /* Properties which will be associated to blocks
     * - unsteady term for water eq. in the liquid phase      (0,0) block
     * - diffusion term for water eq. in the liquid phase     (0,0)-block
     * - unsteady term for water eq. in the gaseous phase     (0,1)-block
     * - unsteady term for hydrogen eq. in the gaseous phase  (1,1)-block
     * - diffusion term for hydrogen eq. in the gaseous phase (1,1)-block
     * - unsteady term for hydrogen eq. in the liquid phase   (1,0)-block
     * - diffusion term for hydrogen eq. in the liquid phase  (1,0)-block
     */

    mc->time_wl_pty = cs_property_add("time_wl_pty", CS_PROPERTY_ISO);
    mc->diff_wl_pty = cs_property_add("diff_wl_pty", perm_type);

    mc->time_wg_pty = cs_property_add("time_wg_pty", CS_PROPERTY_ISO);

    mc->diff_hg_pty = cs_property_add("diff_hg_pty", perm_type);
    mc->time_hg_pty = cs_property_add("time_hg_pty", CS_PROPERTY_ISO);

    mc->time_hl_pty = cs_property_add("time_hl_pty", CS_PROPERTY_ISO);
    mc->diff_hl_pty = cs_property_add("diff_hl_pty", perm_type);

    /* Add terms to the water equation */
    /* ------------------------------- */

    cs_equation_add_time(wl_eqp, mc->time_wl_pty);
    cs_equation_add_diffusion(wl_eqp, mc->diff_wl_pty);

    /* Cross-terms for the block (0,1) -- Water equation */

    cs_equation_add_time(mc->wg_eqp, mc->time_wg_pty);

    /* Add terms to the hydrogen equation */
    /* ---------------------------------- */

    cs_equation_add_diffusion(hg_eqp, mc->diff_hg_pty);

    /* Cross-terms for the block (1,0) -- Hydrogen equation */

    cs_equation_add_time(mc->hl_eqp, mc->time_hl_pty);
    cs_equation_add_diffusion(mc->hl_eqp, mc->diff_hl_pty);

  }
  else { /* Segregated solver */

    mc->use_incremental_solver = true; /* Segregated solver are always solved
                                          by increment */

    /* Properties */
    /* ---------- */

    mc->diff_wl_pty = cs_property_add("diff_wl_pty", perm_type);
    mc->diff_hg_pty = cs_property_add("diff_hg_pty", perm_type);
    mc->diff_hl_pty = cs_property_add("diff_hl_pty", perm_type);

    if (mc->use_properties_on_submesh) {

      mc->time_wl_pty = cs_property_subcell_add("time_wl_pty", CS_PROPERTY_ISO);
      mc->time_hg_pty = cs_property_subcell_add("time_hg_pty", CS_PROPERTY_ISO);
      mc->reac_hg_pty = cs_property_subcell_add("reac_hg_pty", CS_PROPERTY_ISO);

    }
    else {

      mc->time_wl_pty = cs_property_add("time_wl_pty", CS_PROPERTY_ISO);
      mc->time_hg_pty = cs_property_add("time_hg_pty", CS_PROPERTY_ISO);
      mc->reac_hg_pty = cs_property_add("reac_hg_pty", CS_PROPERTY_ISO);

    }

    /* Add terms to the water equation */
    /* ------------------------------- */

    cs_equation_add_time(wl_eqp, mc->time_wl_pty);
    cs_equation_add_diffusion(wl_eqp, mc->diff_wl_pty);


    /* Add terms to the hydrogen equation */
    /* ---------------------------------- */

    cs_equation_add_time(hg_eqp, mc->time_hg_pty);
    cs_equation_add_diffusion(hg_eqp, mc->diff_hg_pty);
    cs_equation_add_reaction(hg_eqp, mc->reac_hg_pty);

  } /* Segregated or coupled solver */

  if (mc->use_incremental_solver) {

    wl_eqp->incremental_algo_type = CS_PARAM_NL_ALGO_PICARD;
    hg_eqp->incremental_algo_type = CS_PARAM_NL_ALGO_PICARD;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the initial setup step in the case of a miscible or
 *        immiscible two-phase flows model in porous media
 *
 * \param[in, out]  mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_init_setup(cs_gwf_two_phase_t     *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  if (mc->wl_eq == NULL || mc->hg_eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Equations are not defined for this model. Stop execution.\n",
              __func__);

  /* Add the variable fields (Keep always the previous state) */

  cs_equation_predefined_create_field(1, mc->wl_eq);
  cs_equation_predefined_create_field(1, mc->hg_eq);

  /* Set fields related to variables */

  mc->l_pressure = cs_equation_get_field(mc->wl_eq);
  mc->g_pressure = cs_equation_get_field(mc->hg_eq);

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");

  /* Retrieve the pointers to the structure storing the equation parameters */

  cs_equation_param_t  *wl_eqp = cs_equation_get_param(mc->wl_eq);
  cs_equation_param_t  *hg_eqp = cs_equation_get_param(mc->hg_eq);

  assert(wl_eqp->space_scheme == hg_eqp->space_scheme);

  /* One has to be consistent with the location of DoFs for the w_eq and h_eq
   * which are respectively related to the l_pressure and g_pressure */

  int loc_id = c_loc_id;
  if (wl_eqp->space_scheme == CS_SPACE_SCHEME_CDOVB ||
      wl_eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB)
    loc_id = v_loc_id;
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme", __func__);

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
 * \brief Perform the last setup step in the case of a two-phase flow model in
 *        porous media (miscible or immiscible case)
 *
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       quant      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_finalize_setup(const cs_cdo_connect_t        *connect,
                    const cs_cdo_quantities_t     *quant,
                    cs_gwf_two_phase_t            *mc)
{
  CS_UNUSED(quant);

  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  const cs_lnum_t  n_cells = connect->n_cells;
  const size_t  csize = n_cells*sizeof(cs_real_t);

  /* Set the Darcian flux (in the volume and at the boundary)
   * Not done up to now. */

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

  BFT_MALLOC(mc->capillarity_cell_pressure, n_cells, cs_real_t);
  memset(mc->capillarity_cell_pressure, 0, sizeof(cs_real_t)*n_cells);

  BFT_MALLOC(mc->g_cell_pressure, n_cells, cs_real_t);
  memset(mc->g_cell_pressure, 0, sizeof(cs_real_t)*n_cells);

  if (mc->use_coupled_solver) {

    BFT_MALLOC(mc->l_capacity, n_cells, cs_real_t);
    memset(mc->l_capacity, 0, sizeof(cs_real_t)*n_cells);

    /* Define the array storing the time property for the water eq. */

    BFT_MALLOC(mc->time_wl_array, n_cells, cs_real_t);
    memset(mc->time_wl_array, 0, csize);

    cs_property_def_by_array(mc->time_wl_pty,
                             cs_flag_primal_cell,  /* where data are located */
                             mc->time_wl_array,
                             false,                /* not owner of the array */
                             NULL, NULL);          /* no index, no ids */

    BFT_MALLOC(mc->time_wg_array, n_cells, cs_real_t);
    memset(mc->time_wg_array, 0, csize);

    cs_property_def_by_array(mc->time_wg_pty,
                             cs_flag_primal_cell,  /* where data are located */
                             mc->time_wg_array,
                             false,                /* not owner of the array */
                             NULL, NULL);          /* no index, no ids */

    BFT_MALLOC(mc->time_hl_array, n_cells, cs_real_t);
    memset(mc->time_hl_array, 0, csize);

    cs_property_def_by_array(mc->time_hl_pty,
                             cs_flag_primal_cell,  /* where data are located */
                             mc->time_hl_array,
                             false,                /* not owner of the array */
                             NULL, NULL);          /* no index, no ids */

    /* Define the array storing the time property for the hydrogen eq. */

    BFT_MALLOC(mc->time_hg_array, n_cells, cs_real_t);
    memset(mc->time_hg_array, 0, csize);

    cs_property_def_by_array(mc->time_hg_pty,
                             cs_flag_primal_cell,  /* where data are located */
                             mc->time_hg_array,
                             false,                /* not owner of the array */
                             NULL, NULL);          /* no index, no ids */

  } /* Only defined for a coupled system */

  /* Define the array storing the diffusion property for the water eq. */

  BFT_MALLOC(mc->diff_wl_array, n_cells, cs_real_t);
  memset(mc->diff_wl_array, 0, csize);

  cs_property_def_by_array(mc->diff_wl_pty,
                           cs_flag_primal_cell,  /* where data are located */
                           mc->diff_wl_array,
                           false,                /* not owner of the array */
                           NULL, NULL);          /* no index, no ids */

  if (mc->use_incremental_solver) {

    if (mc->use_properties_on_submesh) {

      const cs_adjacency_t  *c2v = connect->c2v;
      const cs_lnum_t  c2v_size = c2v->idx[n_cells];
      const size_t  c2v_alloc_size = c2v_size*sizeof(cs_real_t);

      BFT_MALLOC(mc->l_saturation_submesh, c2v_size, cs_real_t);
      memset(mc->l_saturation_submesh, 0, c2v_alloc_size);

      BFT_MALLOC(mc->l_saturation_submesh_pre, c2v_size, cs_real_t);
      memset(mc->l_saturation_submesh_pre, 0, c2v_alloc_size);

      BFT_MALLOC(mc->l_capacity, c2v_size, cs_real_t);
      memset(mc->l_capacity, 0, c2v_alloc_size);

      /* Define the array storing the source term values for the water eq. */

      BFT_MALLOC(mc->time_wl_array, c2v_size, cs_real_t);
      memset(mc->time_wl_array, 0, c2v_alloc_size);

      cs_property_def_by_array(mc->time_wl_pty,
                               cs_flag_dual_cell_byc,  /* data location */
                               mc->time_wl_array,
                               false,                  /* not owner */
                               c2v->idx, c2v->ids);

      /* Define the array storing the source term values for the water eq. */

      BFT_MALLOC(mc->srct_wl_array, c2v_size, cs_real_t);
      memset(mc->srct_wl_array, 0, c2v_alloc_size);

      cs_equation_add_source_term_by_array(cs_equation_get_param(mc->wl_eq),
                                           NULL,        /* all cells */
                                           cs_flag_dual_cell_byc,
                                           mc->srct_wl_array,
                                           false,       /* is owner ? */
                                           c2v->idx, c2v->ids);

      BFT_MALLOC(mc->srct_hg_array, c2v_size, cs_real_t);
      memset(mc->srct_hg_array, 0, c2v_alloc_size);

      cs_equation_add_source_term_by_array(cs_equation_get_param(mc->hg_eq),
                                           NULL,   /* all cells */
                                           cs_flag_dual_cell_byc,
                                           mc->srct_hg_array,
                                           false,  /* is owner ? */
                                           c2v->idx, c2v->ids);

      /* Define the array storing the time property for the hydrogen eq. */

      BFT_MALLOC(mc->time_hg_array, c2v_size, cs_real_t);
      memset(mc->time_hg_array, 0, c2v_alloc_size);

      cs_property_def_by_array(mc->time_hg_pty,
                               cs_flag_dual_cell_byc, /* data location */
                               mc->time_hg_array,
                               false,                 /* not owner */
                               c2v->idx, c2v->ids);

      /* Define the array storing the reaction property for the hydrogen eq. */

      BFT_MALLOC(mc->reac_hg_array, c2v_size, cs_real_t);
      memset(mc->reac_hg_array, 0, c2v_alloc_size);

      cs_property_def_by_array(mc->reac_hg_pty,
                               cs_flag_dual_cell_byc, /* data location */
                               mc->reac_hg_array,
                               false,                 /* not owner */
                               c2v->idx, c2v->ids);

    }
    else { /* Liquid saturation only at cells */

      /* Define the array storing the source term values for the water eq. */

      BFT_MALLOC(mc->srct_wl_array, n_cells, cs_real_t);
      memset(mc->srct_wl_array, 0, csize);

      cs_equation_add_source_term_by_array(cs_equation_get_param(mc->wl_eq),
                                           NULL,   /* all cells */
                                           cs_flag_primal_cell,
                                           mc->srct_wl_array,
                                           false,  /* is owner ? */
                                           NULL, NULL);  /* no index/ids */

      /* Define the array storing the time property for the hydrogen eq. */

      BFT_MALLOC(mc->time_hg_array, n_cells, cs_real_t);
      memset(mc->time_hg_array, 0, csize);

      cs_property_def_by_array(mc->time_hg_pty,
                               cs_flag_primal_cell, /* where data are located */
                               mc->time_hg_array,
                               false,               /* not owner of the array */
                               NULL, NULL);         /* no index, no ids */

      /* Define the array storing the reaction property for the hydrogen eq. */

      BFT_MALLOC(mc->reac_hg_array, n_cells, cs_real_t);
      memset(mc->reac_hg_array, 0, csize);

      cs_property_def_by_array(mc->reac_hg_pty,
                               cs_flag_primal_cell, /* where data are located */
                               mc->reac_hg_array,
                               false,               /* not owner of the array */
                               NULL, NULL);         /* no index/ids */

    }

  } /* Incremental solve */

  /* Define the array storing the diffusion property in the hydrogen eq. */

  BFT_MALLOC(mc->diff_hg_array, n_cells, cs_real_t);
  memset(mc->diff_hg_array, 0, csize);

  cs_property_def_by_array(mc->diff_hg_pty,
                           cs_flag_primal_cell,  /* where data are located */
                           mc->diff_hg_array,
                           false,                /* not owner of the array */
                           NULL, NULL);          /* no index/ids */

  BFT_MALLOC(mc->diff_hl_array, n_cells, cs_real_t);
  memset(mc->diff_hl_array, 0, csize);

  cs_property_def_by_array(mc->diff_hl_pty,
                           cs_flag_primal_cell,  /* where data are located */
                           mc->diff_hl_array,
                           false,                /* not owner of the array */
                           NULL, NULL);          /* no index/ids */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the update step in the case of a two-phase flow model in
 *        porous media (miscible or immiscible)
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      update_flag  metadata associated to type of operation to do
 * \param[in, out] mc           pointer to the casted model context
 *
 *\return the value at which one has to evaluate the properties
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_tpf_updates(const cs_mesh_t             *mesh,
             const cs_cdo_connect_t      *connect,
             const cs_cdo_quantities_t   *quant,
             const cs_time_step_t        *ts,
             cs_flag_t                    update_flag,
             cs_gwf_two_phase_t          *mc)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  cs_real_t  time_eval = ts->t_cur;
  bool  cur2prev = false;

  if (update_flag & CS_FLAG_CURRENT_TO_PREVIOUS) {

    time_eval = _get_time_eval(ts, mc->wl_eq);
    cur2prev = true;

  }

  cs_param_space_scheme_t  space_scheme =
    cs_equation_get_space_scheme(mc->wl_eq);

  /* New pressure values for the liquid and the gas have been computed */

  const cs_real_t  *l_pr = mc->l_pressure->val;
  const cs_real_t  *g_pr = mc->g_pressure->val;

  if (update_flag & CS_FLAG_INITIALIZATION)
    memcpy(mc->g_pressure->val_pre, g_pr,          /* dest, src */
           connect->n_vertices*sizeof(cs_real_t)); /* size */

  /* Update the capillary pressure: P_c = P_g - P_l */

  if (cur2prev) {

    cs_field_current_to_previous(mc->c_pressure);
    cs_field_current_to_previous(mc->l_saturation);

    if (mc->use_properties_on_submesh)
      memcpy(mc->l_saturation_submesh_pre,                           /* dest */
             mc->l_saturation_submesh,                               /* src */
             connect->c2v->idx[connect->n_cells]*sizeof(cs_real_t)); /* size */

  }

  cs_real_t  *c_pr = mc->c_pressure->val;

  switch (space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    /* Compute the value at vertices */

    for (cs_lnum_t i = 0; i < quant->n_vertices; i++)
      c_pr[i] = g_pr[i] - l_pr[i];

    /* Now interpolate the vertex values to the cell values. the capillarity
       pressure at cells is the one used to update quantities related to a soil
       model */

    cs_reco_pv_at_cell_centers(connect->c2v, quant, c_pr,
                               mc->capillarity_cell_pressure);

    /* Interpolate the pressure in the gaseous phase at cell centers */

    cs_reco_pv_at_cell_centers(connect->c2v, quant, g_pr,
                               mc->g_cell_pressure);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme", __func__);

  }

  /* Update properties related to soils:
   * - l_rel_permeability, g_rel_permeability
   * - liquid_saturation
   * - capacity: \frac{\partial S_l}{\partial P_c}
   *
   * Either a call to a user-defined function or a predefined function if the
   * soil corresponds to a known model
   */

  cs_gwf_soil_update(time_eval, mesh, connect, quant);

  /* Define the liquid saturation in each cell in some specific situations */

  if (mc->use_properties_on_submesh) {

    const cs_adjacency_t  *c2v = connect->c2v;
    cs_real_t  *l_sat = mc->l_saturation->val;

    /* Compute the average liquid saturation in each cell */

    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      cs_real_t  _sat = 0;
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        _sat += mc->l_saturation_submesh[j] * quant->dcell_vol[j];
      l_sat[c_id] = _sat/quant->cell_vol[c_id];

    } /* Loop on cells */

    if (update_flag & CS_FLAG_INITIALIZATION)
      memcpy(mc->l_saturation_submesh_pre,                           /* dest */
             mc->l_saturation_submesh,                               /* src */
             connect->c2v->idx[connect->n_cells]*sizeof(cs_real_t)); /* size */

  }

  /* TODO: Update the Darcy advection field for the liquid and the gas phase */

  /* Update arrays associated to property terms */

  int  dim = cs_property_get_dim(gw->abs_permeability);

  switch (dim) {

  case 1: /* Isotropic case */
    if (gw->model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE)
      cs_gwf_soil_iso_update_mtpf_terms(mc);

    else {

      if (mc->use_incremental_solver) {

        if (mc->use_properties_on_submesh)
          cs_gwf_soil_iso_update_itpf_terms_incr_submesh(ts, connect, mc);
        else
          cs_gwf_soil_iso_update_itpf_terms_incr(ts, mc);

      }
      else
        cs_gwf_soil_iso_update_itpf_terms(mc);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Only the isotropic case is coded up to now.",
              __func__);
    break;

  } /* Switch on the permeability type */

  return time_eval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the groundwater flow module in case
 *         of miscible or immiscible two-phase flows in porous media
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  mc        pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_extra_op(const cs_cdo_connect_t                *connect,
              const cs_cdo_quantities_t             *cdoq,
              cs_gwf_unsaturated_single_phase_t     *mc)
{
  CS_UNUSED(connect);
  CS_UNUSED(cdoq);
  CS_UNUSED(mc);

  cs_gwf_t  *gw = cs_gwf_main_structure;
  assert(gw != NULL && mc != NULL);

  if (cs_flag_test(gw->post_flag, CS_GWF_POST_DARCY_FLUX_BALANCE) == false)
    return; /* Nothing to do */

  /* TO BE DONE (not useful up to now) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new increment on the capillarity pressure
 *
 * \param[in]      mesh     pointer to a cs_mesh_t structure
 * \param[in]      wl_eq    pointer on tne "wl" equation structure
 * \param[in]      hg_eq    pointer on tne "wl" equation structure
 * \param[in, out] dpc_kp1  values of the increment on the capillarity pressure
 */
/*----------------------------------------------------------------------------*/

static void
_get_capillarity_pressure_increment(const cs_mesh_t       *mesh,
                                    const cs_equation_t   *wl_eq,
                                    const cs_equation_t   *hg_eq,
                                    cs_real_t             *dpc_kp1)
{
  const cs_equation_param_t  *wl_eqp = cs_equation_get_param(wl_eq);
  const cs_equation_param_t  *hg_eqp = cs_equation_get_param(hg_eq);
  const cs_equation_builder_t  *wl_eqb = cs_equation_get_builder(wl_eq);
  const cs_equation_builder_t  *hg_eqb = cs_equation_get_builder(hg_eq);
  const cs_real_t  *dpl_kp1 = wl_eqb->increment;
  const cs_real_t  *dpg_kp1 = hg_eqb->increment;

  if (wl_eqp->incremental_relax_factor < 1 ||
      hg_eqp->incremental_relax_factor < 1) {

    assert(wl_eqp->incremental_relax_factor > 0);
    assert(hg_eqp->incremental_relax_factor > 0);

    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++)
      dpc_kp1[i] = hg_eqp->incremental_relax_factor*dpg_kp1[i]
                 - wl_eqp->incremental_relax_factor*dpl_kp1[i];

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
 * \brief  Check the convergence of the non-linear algorithm in case of TPF
 *         model. Work only with the capillarity pressure at vertices
 *
 * \param[in]       nl_algo_type type of non-linear algorithm
 * \param[in]       dpc_iter     cur. increment values for the capill. pressure
 * \param[in, out]  algo         pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_check_nl_tpf_dpc_cvg(cs_param_nl_algo_t        nl_algo_type,
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

  cs_iter_algo_update_cvg(algo);

  if (algo->param.verbosity > 0) {

    if (algo->n_algo_iter == 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    "### GWF.TPF %10s.It    Algo.Res   Tolerance\n",
                    cs_param_get_nl_algo_label(nl_algo_type));

    cs_log_printf(CS_LOG_DEFAULT,
                  "### GWF.TPF %10s.It%02d  %5.3e  %6.4e\n",
                  cs_param_get_nl_algo_label(nl_algo_type),
                  algo->n_algo_iter, algo->res, algo->tol);

  } /* verbosity > 0 */

  return algo->cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check the convergence of the non-linear algorithm in case of TPF
 *         model
 *
 * \param[in]       nl_algo_type type of non-linear algorithm
 * \param[in]       pg_pre_iter  previous iterate values for the gas pressure
 * \param[in, out]  pg_cur_iter  current iterate values for the gas pressure
 * \param[in]       pl_pre_iter  previous iterate values for the liquid pressure
 * \param[in, out]  pl_cur_iter  current iterate values for the liquid pressure
 * \param[in, out]  algo         pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_check_nl_tpf_cvg(cs_param_nl_algo_t        nl_algo_type,
                  const cs_real_t          *pg_pre_iter,
                  cs_real_t                *pg_cur_iter,
                  const cs_real_t          *pl_pre_iter,
                  cs_real_t                *pl_cur_iter,
                  cs_iter_algo_t           *algo)
{
  assert(algo != NULL);

  if (nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON && algo->n_algo_iter > 0) {

    /* pg_* arrays gather pg and pl (this is done during the solve step) */

    cs_iter_algo_aa_update(algo,
                           pg_cur_iter, /* updated during the process */
                           pg_pre_iter,
                           cs_cdo_blas_dotprod_2pvsp,
                           cs_cdo_blas_square_norm_2pvsp);

  } /* Anderson acceleration */

  algo->prev_res = algo->res;

  double delta_pg = cs_cdo_blas_square_norm_pvsp_diff(pg_pre_iter, pg_cur_iter);
  double delta_pl = cs_cdo_blas_square_norm_pvsp_diff(pl_pre_iter, pl_cur_iter);

  algo->res = delta_pg + delta_pl;
  assert(algo->res > -DBL_MIN);
  algo->res = sqrt(algo->res);

  if (algo->n_algo_iter < 1) /* Store the first residual to detect a
                                possible divergence of the algorithm */
    algo->res0 = algo->res;

  /* Update the convergence members */

  cs_iter_algo_update_cvg(algo);

  if (algo->param.verbosity > 0) {

    if (algo->n_algo_iter == 1) {

      if (algo->param.verbosity > 1)
        cs_log_printf(CS_LOG_DEFAULT,
                      "### GWF.TPF %10s.It    Algo.Res   Tolerance"
                      "  ||D_Pg||  ||D_Pl||\n",
                      cs_param_get_nl_algo_label(nl_algo_type));
      else
        cs_log_printf(CS_LOG_DEFAULT,
                      "### GWF.TPF %10s.It    Algo.Res   Tolerance\n",
                      cs_param_get_nl_algo_label(nl_algo_type));

    }

    if (algo->param.verbosity > 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    "### GWF.TPF %10s.It%02d  %5.3e  %6.4e %5.3e %5.3e\n",
                    cs_param_get_nl_algo_label(nl_algo_type),
                    algo->n_algo_iter, algo->res, algo->tol,
                    sqrt(delta_pg), sqrt(delta_pl));
    else
      cs_log_printf(CS_LOG_DEFAULT,
                    "### GWF.TPF %10s.It%02d  %5.3e  %6.4e\n",
                    cs_param_get_nl_algo_label(nl_algo_type),
                    algo->n_algo_iter, algo->res, algo->tol);

  } /* verbosity > 0 */

  return algo->cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for the groundwater flows module.
 *         Case of (miscible or immiscible) two-phase flows in porous media
 *         with a non-linear resolution relying on the Picard algorithm and
 *         a coupled algorithm
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_coupled_tpf_picard_compute(const cs_mesh_t              *mesh,
                            const cs_cdo_connect_t       *connect,
                            const cs_cdo_quantities_t    *cdoq,
                            const cs_time_step_t         *time_step,
                            cs_gwf_two_phase_t           *mc)
{
  bool cur2prev = false;
  cs_flag_t  update_flag = 0;   /* No current to previous operation */

  cs_iter_algo_t  *algo = mc->nl_algo;
  assert(algo != NULL);

  cs_iter_algo_reset_nl(mc->nl_algo_type, algo);

  /* A first resolution has been done followed by an update */

  cs_real_t  *pg_kp1 = NULL, *pl_kp1 = NULL;  /* at ^{n+1,k+1} */
  cs_real_t  *pg_k = NULL, *pl_k = NULL;      /* at ^{n+1,k} */

  const size_t  vsize = cdoq->n_vertices*sizeof(cs_real_t);

  BFT_MALLOC(pg_k, 2*cdoq->n_vertices, cs_real_t);
  pl_k = pg_k + cdoq->n_vertices;

  memcpy(pg_k, mc->g_pressure->val_pre, vsize);
  memcpy(pl_k, mc->l_pressure->val_pre, vsize);

  pg_kp1 = mc->g_pressure->val;
  pl_kp1 = mc->l_pressure->val;

  /* Set the normalization factor */

  algo->normalization = cs_cdo_blas_square_norm_pvsp(pg_k);
  algo->normalization += cs_cdo_blas_square_norm_pvsp(pl_k);
  algo->normalization = sqrt(algo->normalization);
  if (algo->normalization < cs_math_zero_threshold)
    algo->normalization = 1.0;

  cs_log_printf(CS_LOG_DEFAULT, "%s: normalization=%6.4e\n",
                __func__, algo->normalization);

  /* Main non-linear loop */

  while (_check_nl_tpf_cvg(mc->nl_algo_type,
                           pg_k, pg_kp1, pl_k, pl_kp1,
                           algo) == CS_SLES_ITERATING) {


    memcpy(pg_k, pg_kp1, vsize);
    memcpy(pl_k, pl_kp1, vsize);

    /* Build and solve the linear system related to the coupled system of
       equations. First call: current --> previous and then no operation */

    cs_equation_system_solve(cur2prev, mc->system);

    /* Get the soil state up to date with the last compute values for the
       liquid and gas pressures */

    /* relaxation */

    const double  relax = 0.;
    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++) {

      pg_kp1[i] = (1 - relax)*pg_kp1[i] + relax*pg_k[i];
      pl_kp1[i] = (1 - relax)*pl_kp1[i] + relax*pl_k[i];

    }

    cs_gwf_update(mesh, connect, cdoq, time_step, update_flag);

  } /* Until convergence */

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_post_check(__func__,
                          mc->system->param->name,
                          cs_param_get_nl_algo_label(mc->nl_algo_type),
                          algo);

  /* Free temporary arrays and structures */

  BFT_FREE(pg_k);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for the groundwater flows module.
 *         Case of (miscible or immiscible) two-phase flows in porous media
 *         with a non-linear resolution relying on the Anderson algorithm and
 *         a coupled algorithm
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_coupled_tpf_anderson_compute(const cs_mesh_t              *mesh,
                              const cs_cdo_connect_t       *connect,
                              const cs_cdo_quantities_t    *cdoq,
                              const cs_time_step_t         *time_step,
                              cs_gwf_two_phase_t           *mc)
{
  bool  cur2prev = false;
  cs_flag_t  update_flag = 0;   /* No current to previous operation */
  cs_iter_algo_t  *algo = mc->nl_algo;
  assert(algo != NULL);

  cs_iter_algo_reset_nl(mc->nl_algo_type, algo);

  /* A first resolution has been done followed by an update */

  cs_real_t  *pg_kp1 = NULL, *pl_kp1 = NULL;  /* at ^{n+1,k+1} */
  cs_real_t  *pg_k = NULL, *pl_k = NULL;      /* at ^{n+1,k} */

  const size_t  vsize = cdoq->n_vertices*sizeof(cs_real_t);

  BFT_MALLOC(pg_k, 2*cdoq->n_vertices, cs_real_t);
  pl_k = pg_k + cdoq->n_vertices;

  memcpy(pg_k, mc->g_pressure->val_pre, vsize);
  memcpy(pl_k, mc->l_pressure->val_pre, vsize);

  /* One needs only one array gathering the liquid and gas pressures */

  BFT_MALLOC(pg_kp1, 2*cdoq->n_vertices, cs_real_t);
  pl_kp1 = pg_kp1 + cdoq->n_vertices;

  memcpy(pg_kp1, mc->g_pressure->val, vsize);
  memcpy(pl_kp1, mc->l_pressure->val, vsize);

  /* Set the normalization factor */

  algo->normalization = cs_cdo_blas_square_norm_pvsp(pg_k);
  algo->normalization += cs_cdo_blas_square_norm_pvsp(pl_k);
  algo->normalization = sqrt(algo->normalization);
  if (algo->normalization < cs_math_zero_threshold)
    algo->normalization = 1.0;

  /* Main non-linear loop */

  while (_check_nl_tpf_cvg(mc->nl_algo_type,
                           pg_k, pg_kp1, pl_k, pl_kp1,
                           algo) == CS_SLES_ITERATING) {

    memcpy(pg_k, pg_kp1, vsize);
    memcpy(pl_k, pl_kp1, vsize);

    /* Update the variables related to the groundwater flow system.
     * In case of an Anderson acceleration, pg_kp1 and pl_kp1 may be
     * updated */

    if (algo->n_algo_iter >= mc->anderson_param.starting_iter) {

      memcpy(mc->g_pressure->val, pg_kp1, vsize);
      memcpy(mc->l_pressure->val, pl_kp1, vsize);

    }

    cs_gwf_update(mesh, connect, cdoq, time_step, update_flag);

    /* Build and solve the linear system related to the coupled system of
       equations. First call: current --> previous and then no operation */

    cs_equation_system_solve(cur2prev, mc->system);

    memcpy(pg_kp1, mc->g_pressure->val, vsize);
    memcpy(pl_kp1, mc->l_pressure->val, vsize);

  } /* Until convergence */

  /* Get the soil state up to date with the last compute values for the
     liquid and gas pressures */

  cs_gwf_update(mesh, connect, cdoq, time_step, update_flag);

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
 * \brief  Compute the new (unsteady) state for the groundwater flows module.
 *         Case of (miscible or immiscible) two-phase flows in porous media.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_segregated_tpf_compute(const cs_mesh_t              *mesh,
                        const cs_time_step_t         *time_step,
                        const cs_cdo_connect_t       *connect,
                        const cs_cdo_quantities_t    *cdoq,
                        cs_gwf_two_phase_t           *mc)
{
  assert(mc->use_incremental_solver);

  cs_field_current_to_previous(mc->g_pressure);
  cs_field_current_to_previous(mc->l_pressure);

  bool cur2prev = false;        /* Done just above */
  cs_flag_t  update_flag = 0;   /* No current to previous operation */

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

  cs_equation_solve(cur2prev, mesh, mc->wl_eq);

  cs_equation_solve(cur2prev, mesh, mc->hg_eq);

  /* Update the variables related to the groundwater flow system */

  cs_gwf_update(mesh, connect, cdoq, time_step,
                CS_FLAG_CURRENT_TO_PREVIOUS); /* Force this operation */

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

  _get_capillarity_pressure_increment(mesh, mc->wl_eq, mc->hg_eq, dpc_kp1);

  while(_check_nl_tpf_dpc_cvg(mc->nl_algo_type,
                              dpc_kp1, algo) == CS_SLES_ITERATING) {

    /* Solve step */

    cs_equation_solve(cur2prev, mesh, mc->wl_eq);

    cs_equation_solve(cur2prev, mesh, mc->hg_eq);

    /* Update the variables related to the groundwater flow system */

    cs_gwf_update(mesh, connect, cdoq, time_step, update_flag);

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

    _get_capillarity_pressure_increment(mesh, mc->wl_eq, mc->hg_eq, dpc_kp1);

  } /* while not converged */

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_post_check(__func__,
                          "Segregated incremental TPF solver",
                          cs_param_get_nl_algo_label(mc->nl_algo_type),
                          algo);

  BFT_FREE(dpc_kp1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for the groundwater flows module.
 *         Case of (miscible or immiscible) two-phase flows in porous media.
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
#if defined(DEBUG) && !defined(NDEBUG)
  cs_gwf_t  *gw = cs_gwf_main_structure;

  assert(gw != NULL && mc != NULL);
  assert(gw->model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE ||
         gw->model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE);
  assert(cs_equation_get_type(mc->wl_eq) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(cs_equation_get_type(mc->hg_eq) == CS_EQUATION_TYPE_GROUNDWATER);
#endif

  if (cs_equation_is_steady(mc->wl_eq) && cs_equation_is_steady(mc->hg_eq)) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  "Unsteady computation whereas all equations are steady.\n");
  }

  if (mc->use_coupled_solver) {

    bool cur2prev = true;

    /* Build and solve the linear system related to the coupled system of
     * equations. By default, a current to previous operation is performed so
     * that prev <-- n and cur <-- n+1,k=1 (the first guess for the next time
     * step) */

    cs_equation_system_solve(cur2prev, mc->system);

    /* Update the variables related to the groundwater flow system */

    cs_gwf_update(mesh, connect, cdoq, time_step, CS_FLAG_CURRENT_TO_PREVIOUS);

    switch (mc->nl_algo_type) {

    case CS_PARAM_NL_ALGO_PICARD:
      _coupled_tpf_picard_compute(mesh, connect, cdoq, time_step, mc);
      break;

    case CS_PARAM_NL_ALGO_ANDERSON:
      _coupled_tpf_anderson_compute(mesh, connect, cdoq, time_step, mc);
      break;

    default:
      break; /* Nothing else to do */
    }

  }
  else
    _segregated_tpf_compute(mesh, time_step, connect, cdoq, mc);
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
 * \param[in]   model           type of physical modelling
 * \param[in]   option_flag     optional flag to specify this module
 * \param[in]   post_flag       optional automatic postprocessing
 *
 * \return a pointer to a new allocated groundwater flow structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_activate(cs_gwf_model_type_t      model,
                cs_flag_t                option_flag,
                cs_flag_t                post_flag)
{
  cs_gwf_t  *gw = _gwf_create();

  /* Store the pointer to the groundawater flow structure */

  cs_gwf_main_structure = gw;

  gw->model = model;
  gw->flag = option_flag;

  /* Add the porosity property */

  gw->soil_porosity = cs_property_add("soil_porosity", CS_PROPERTY_ISO);

  /* Allocate and initialize each model context (mc) */

  switch (model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    gw->model_context = _sspf_activate();
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    gw->post_flag |= CS_GWF_POST_LIQUID_SATURATION;
    gw->model_context = _uspf_activate();
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    gw->post_flag |= CS_GWF_POST_LIQUID_SATURATION;
    gw->model_context = _tpf_activate();
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);

  }

  /* Now one can set the post_flag (need the initializtion of the model
     context) */

  cs_gwf_set_post_options(post_flag, false);

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

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    {
      cs_gwf_two_phase_t  *mc = gw->model_context;

      _tpf_free_context(&(mc));
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
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Gravitation: *True*\n");
  else
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Gravitation: *False*\n");

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
  bool  post_gas_density =
    (gw->post_flag & CS_GWF_POST_GAS_MASS_DENSITY) ? true : false;

  cs_log_printf(CS_LOG_SETUP, "  * GWF | Post:"
                " Soil capacity %s Liquid saturation %s Permeability %s\n",
                cs_base_strtf(post_capacity),
                cs_base_strtf(post_liquid_saturation),
                cs_base_strtf(post_permeability));

  if (gw->model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE ||
      gw->model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Post: Gas mass density %s\n",
                  cs_base_strtf(post_gas_density));

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

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Model: **%s**\n", cs_gwf_model_name[gw->model]);

  switch(gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    _sspf_log_context(gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    _uspf_log_context(gw->model_context);
    break;

  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    _itpf_log_context(gw->model_context);
    _tpf_log_context(gw->model_context);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
    _mtpf_log_context(gw->model_context);
    _tpf_log_context(gw->model_context);
    break;

  default:
    break;

  } /* Type of model */

  /* Soils */

  cs_gwf_soil_log_setup();

  /* Tracers */

  cs_gwf_tracer_log_all();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the main structure which manages a two-phase flow model
 *
 * \return a pointer to the structure cs_gwf_two_phase_t
 */
/*----------------------------------------------------------------------------*/

cs_gwf_two_phase_t *
cs_gwf_get_two_phase_model(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (gw->model != CS_GWF_MODEL_MISCIBLE_TWO_PHASE &&
      gw->model != CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid model. One expects a two-phase flow model.\n",
              __func__);

  cs_gwf_two_phase_t  *mc = gw->model_context;

  assert(mc != NULL);
  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the numerical options related to the two phase flow models
 *
 * \param[in] use_coupled_solver         true/false
 * \param[in] use_incremental_solver     true/false
 * \param[in] use_properties_on_submesh  true/false
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_two_phase_numerical_options(bool    use_coupled_solver,
                                       bool    use_incremental_solver,
                                       bool    use_properties_on_submesh)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  cs_gwf_two_phase_t  *mc = gw->model_context;
  assert(mc != NULL);

  mc->use_coupled_solver = use_coupled_solver;
  mc->use_incremental_solver = use_incremental_solver;
  mc->use_properties_on_submesh = use_properties_on_submesh;
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
cs_gwf_set_miscible_two_phase_model(cs_real_t       l_mass_density,
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
  if (gw->model != CS_GWF_MODEL_MISCIBLE_TWO_PHASE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid model. One expects a two-phase flow model.\n",
              __func__);

  cs_gwf_two_phase_t  *mc = gw->model_context;

  assert(mc != NULL);
  assert(l_mass_density > 0);
  assert(ref_temperature > 0);  /* In Kelvin */
  assert(h_molar_mass > 0);
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
 * \brief  Set the parameters defining the immiscible two-phase flow model.
 *         Use SI unit if not prescribed otherwise.
 *
 * \param[in] l_mass_density   mass density of the main liquid component
 * \param[in] l_viscosity      viscosity in the liquid phase (Pa.s)
 * \param[in] g_viscosity      viscosity in the gas phase (Pa.s)
 * \param[in] h_molar_mass     molar mass of the main gas component
 * \param[in] ref_temperature  reference temperature in Kelvin
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_immiscible_two_phase_model(cs_real_t       l_mass_density,
                                      cs_real_t       l_viscosity,
                                      cs_real_t       g_viscosity,
                                      cs_real_t       h_molar_mass,
                                      cs_real_t       ref_temperature)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));
  if (gw->model != CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid model. One expects a two-phase flow model.\n",
              __func__);

  cs_gwf_two_phase_t  *mc = gw->model_context;

  assert(mc != NULL);
  assert(l_mass_density > 0);
  assert(ref_temperature > 0);  /* In Kelvin */
  assert(h_molar_mass > 0);
  assert(l_viscosity > 0 && g_viscosity > 0);

  /* Set the parameters */

  mc->l_mass_density = l_mass_density;
  mc->l_viscosity = l_viscosity;
  mc->g_viscosity = g_viscosity;
  mc->l_diffusivity_h = 0;      /* immiscible case */
  mc->w_molar_mass = 0;         /* immiscible case */
  mc->h_molar_mass = h_molar_mass;
  mc->ref_temperature = ref_temperature;
  mc->henry_constant = 1e-20;  /* immiscible case */
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

    if (gw->model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE ||
        gw->model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE) {

      cs_gwf_two_phase_t  *mc = (cs_gwf_two_phase_t *)gw->model_context;

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
 *         Case of a soil with an isotropic absolute permeability
 *
 * \param[in]  z_name      name of the volume zone corresponding to the soil
 * \param[in]  density     value of the bulk mass density
 * \param[in]  k_abs       absolute (or intrisic) permeability (scalar-valued)
 * \param[in]  porosity    value of the porosity (saturated moisture content)
 * \param[in]  model       type of model for the soil behavior
 *
 * \return a pointer to the new allocated soil structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_add_iso_soil(const char                *z_name,
                    double                     density,
                    double                     k_abs,
                    double                     porosity,
                    cs_gwf_soil_model_t        model)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  const cs_zone_t  *zone = cs_volume_zone_by_name_try(z_name);

  if (zone == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Zone %s related to the same soil is not defined.\n"
              " Stop adding a new soil.", z_name);

  assert(density > 0);
  assert(k_abs > 0);
  assert(porosity > 0);

  double  ktens_abs[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  ktens_abs[0][0] = ktens_abs[1][1] = ktens_abs[2][2] = k_abs;

  cs_gwf_soil_t  *soil = cs_gwf_soil_create(zone,
                                            gw->model, /* hydraulic model */
                                            model,     /* soil model */
                                            CS_PROPERTY_ISO,
                                            ktens_abs,
                                            porosity,
                                            density,
                                            gw->model_context);

  return soil;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and add a new cs_gwf_soil_t structure. An initialization by
 *         default of all members is performed.
 *
 * \param[in]  z_name      name of the volume zone corresponding to the soil
 * \param[in]  density     value of the bulk mass density
 * \param[in]  k_abs       absolute (or intrisic) permeability (tensor-valued)
 * \param[in]  porosity    value of the porosity (saturated moisture content)
 * \param[in]  model       type of model for the soil behavior
 *
 * \return a pointer to the new allocated soil structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_add_aniso_soil(const char                *z_name,
                      double                     density,
                      double                     k_abs[3][3],
                      double                     porosity,
                      cs_gwf_soil_model_t        model)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  const cs_zone_t  *zone = cs_volume_zone_by_name_try(z_name);

  if (zone == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Zone %s related to the same soil is not defined.\n"
              " Stop adding a new soil.", z_name);

  assert(density > 0);
  assert(porosity > 0);
  assert(k_abs[0][0]*k_abs[0][0] +
         k_abs[1][1]*k_abs[1][1] +
         k_abs[2][2]*k_abs[2][2] > 0);

  cs_gwf_soil_t  *soil = cs_gwf_soil_create(zone,
                                            gw->model, /* hydraulic model */
                                            model,     /* soil model */
                                            CS_PROPERTY_ANISO,
                                            k_abs,
                                            porosity,
                                            density,
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

  cs_gwf_tracer_init_setup_t  *init_setup = cs_gwf_tracer_default_init_setup;
  cs_gwf_tracer_finalize_setup_t  *finalize_setup = NULL;

  if (gw->model == CS_GWF_MODEL_SATURATED_SINGLE_PHASE)
    finalize_setup = cs_gwf_tracer_sat_finalize_setup;
  else
    finalize_setup = cs_gwf_tracer_unsat_finalize_setup;

  /* Call the main function to add a new tracer */

  assert(finalize_setup != NULL);
  cs_gwf_tracer_t  *tracer = cs_gwf_tracer_add(tr_model,
                                               gw->model,
                                               eq_name,
                                               var_name,
                                               adv,
                                               init_setup,
                                               finalize_setup);

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
 * \param[in] eq_name         name of the tracer equation
 * \param[in] var_name        name of the related variable
 * \param[in] init_setup      function pointer (predefined prototype)
 * \param[in] finalize_setup  function pointer (predefined prototype)
 *
 * \return a pointer to the new cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_add_user_tracer(const char                       *eq_name,
                       const char                       *var_name,
                       cs_gwf_tracer_init_setup_t       *init_setup,
                       cs_gwf_tracer_finalize_setup_t   *finalize_setup)
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
                                               init_setup,
                                               finalize_setup);

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the context of the model after the activation of the
 *         module and a first settings of the model parameters (physical and
 *         numerical). At this stage, cs_user_parameters() has not been called
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_model_context(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));
  cs_gwf_soil_check();

  int dim = cs_gwf_soil_get_permeability_max_dim();
  cs_property_type_t  perm_type = CS_PROPERTY_ISO;
  if (dim == 9)
    perm_type = CS_PROPERTY_ANISO;

  /* Add the absolute (or intrisic) permeability property (used in the
   * definition of the diffusion term in the conservation equations and in the
   * definition of the Darcy flux).
   */

  gw->abs_permeability = cs_property_add("absolute_permeability", perm_type);

  /* Continue the setup of the model and create new fields */

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    _sspf_init_model_context(gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    _uspf_init_model_context(gw->model_context, perm_type);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    _tpf_init_model_context(gw->model_context, perm_type);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined settings for the groundwater flow model and its related
 *         equations.
 *         At this stage, all soils have been defined and equation parameters
 *         are set (cs_user_parameters() has been called).
 *         Create new cs_field_t structures according to the setting.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_setup(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  /* Continue the setup of the model and create new fields */

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    _sspf_init_setup(gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    _uspf_init_setup(gw->model_context);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    _tpf_init_setup(gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Same step for the tracer equations associated to the GWF module */

  cs_gwf_tracer_init_setup();
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

  /* Set the soil porosity and the absolute permeability from the soil
     definition */

  cs_gwf_soil_set_shared_properties(gw->abs_permeability,
                                    gw->soil_porosity);

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    _sspf_finalize_setup(connect, quant, gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    _uspf_finalize_setup(connect, quant, gw->model_context);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    _tpf_finalize_setup(connect, quant, gw->model_context);
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

  cs_gwf_tracer_finalize_setup(connect, quant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the GWF module (done after all the setup phase and after
 *         the initialization of all equations)
 *         One sets an initial value to all quantities related to this module.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_values(const cs_mesh_t             *mesh,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   const cs_time_step_t        *ts)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  cs_gwf_update(mesh, connect, quant, ts, CS_FLAG_INITIALIZATION);

  /* Further steps dedicated to each model */

  switch (gw->model) {

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    {
      cs_gwf_two_phase_t *mc = gw->model_context;

      if (mc->use_incremental_solver) {

        /* This model is non-linear so that one needs an iterative algorithm to
           manage this non-linearity */

        if (mc->nl_algo_type == CS_PARAM_NL_ALGO_NONE)
          mc->nl_algo_type = CS_PARAM_NL_ALGO_PICARD;

        /* Add a hook function to define the source term. This operation should
           be done after that the corresponding equation has been initialized
           (a builder structure has to be defined) */

        cs_equation_add_build_hook(mc->hg_eq,
                                   mc,                 /* hook context */
                                   _build_hg_incr_st); /* hook function */

        /* Initialise other previous quantities */

        if (mc->use_properties_on_submesh) {

          memcpy(mc->l_saturation_submesh_pre,                      /* dest */
                 mc->l_saturation_submesh,                          /* src */
                 connect->c2v->idx[connect->n_cells]*sizeof(cs_real_t));

        }

      }

      /* Initialize the non-linear algo. if needed */

      if (mc->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

        mc->nl_algo = cs_iter_algo_create(mc->nl_algo_param);

        /* One assumes that the discretization schemes are all set to CDO
           vertex-based schemes */

        cs_lnum_t  size = quant->n_vertices;
        if (mc->use_coupled_solver)
          size *= 2;

        if (mc->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
          mc->nl_algo->context = cs_iter_algo_aa_create(mc->anderson_param,
                                                        size);

      }
    }
    break;

  default:
    break; /* Nothing else to do */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the groundwater system (pressure head, head in law, moisture
 *         content, darcian velocity, soil capacity or permeability if needed)
 *
 * \param[in]  mesh         pointer to a cs_mesh_t structure
 * \param[in]  connect      pointer to a cs_cdo_connect_t structure
 * \param[in]  quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts           pointer to a cs_time_step_t structure
 * \param[in]  update_flag  metadata associated to the status of the update
 *                          step to perform
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_update(const cs_mesh_t             *mesh,
              const cs_cdo_connect_t      *connect,
              const cs_cdo_quantities_t   *quant,
              const cs_time_step_t        *ts,
              cs_flag_t                    update_flag)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Groundwater module is not allocated.", __func__);

  cs_real_t  time_eval = ts->t_cur; /* default value */

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    time_eval = _sspf_updates(connect, quant, ts, update_flag,
                              gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    time_eval = _uspf_updates(mesh, connect, quant, ts, update_flag,
                              gw->model_context);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    time_eval = _tpf_updates(mesh, connect, quant, ts, update_flag,
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

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    break; /* Nothing to do (the w_eq can be steady according to the numerical
              choices but this resolution is performed elsewhere) */
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

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    _tpf_compute(mesh, time_step, connect, cdoq, gw->model_context);
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

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
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
  const cs_gwf_two_phase_t  *mc = gw->model_context;

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

    if (gw->post_flag & CS_GWF_POST_SOIL_CAPACITY) {

      if (mc->l_capacity == NULL) {
        cs_base_warn(__FILE__, __LINE__);
        bft_printf("%s: Requested postprocessing for capacity but not"
                   " allocated.\n This postprocessing is skipped.\n", __func__);
      }
      else {

        if (mc->use_properties_on_submesh) {
          cs_base_warn(__FILE__, __LINE__);
          bft_printf("%s: Requested postprocessing for capacity at cells but"
                     " this functionnality is not yet available\n"
                     "%s: with CS_GWF_LIQUID_SATURATION_ON_SUBMESH\n",
                     __func__, __func__);
        }
        else
          cs_post_write_var(mesh_id,
                            CS_POST_WRITER_DEFAULT,
                            "l_capacity",
                            1,
                            true,
                            false,
                            CS_POST_TYPE_cs_real_t,
                            mc->l_capacity,
                            NULL,
                            NULL,
                            time_step);

      }

    } /* Postprocess the soil capacity */

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

    if (gw->post_flag & CS_GWF_POST_GAS_MASS_DENSITY) {

      const double  mh_ov_rt =
        mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);

      cs_real_t  *gas_mass_density = NULL;
      BFT_MALLOC(gas_mass_density, n_cells, cs_real_t);

      for (cs_lnum_t c = 0; c < n_cells; c++)
        gas_mass_density[c] = mh_ov_rt * mc->g_cell_pressure[c];

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

  } /* volume mesh id */
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
