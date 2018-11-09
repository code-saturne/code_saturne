/*============================================================================
 * Main functions dedicated to groundwater flows
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_param.h"
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
 * Structure definitions
 *============================================================================*/

/* Set of parameters related to the groundwater flow (GWF) module */

struct _gwf_t {

  cs_flag_t          flag;      /* Flag dedicated to general options to handle
                                 * the GWF module*/
  cs_flag_t          post_flag; /* Flag dedicated to the post-processing
                                 * of the GWF module */

  /* Gravity effect */
  cs_real_3_t        gravity;

  /* Set of equations associated to this module */
  cs_equation_t     *richards;  /* "hydraulic_head" is the associated variable
                                   "permeability" is the diffusion property
                                   related to this equation */

  /* Members related to the associated tracer equations */
  int                          n_tracers;
  cs_gwf_tracer_t            **tracers;
  cs_gwf_tracer_setup_t      **finalize_tracer_setup;  /* Function pointers */
  cs_gwf_tracer_add_terms_t  **add_tracer_terms;       /* Function pointers */

  /* Additional heads */
  cs_field_t      *pressure_head;    /* Allocated only if gravitation is active
                                        Location depends on the discretization
                                        scheme used to solve Richards eq.
                                        pressure head is denoted by h
                                        hydraulic head (solved in Richards eq.)
                                        is denoted by H.
                                        h = H - gravity_potential */
  cs_real_t       *head_in_law;      /* Array used as an input in laws */

  /* Moisture content: not constant for unsaturated soils */
  cs_property_t   *moisture_content;
  cs_field_t      *moisture_field;   /* Related cs_field_t structure at cells */

  /* Soil capacity: property attached to the unsteady term in the Richards
     equation */
  cs_property_t   *soil_capacity;
  cs_field_t      *capacity_field;   /* Related cs_field_t structure at cells */

  /* Permeability is the diffusion property related to Richards equation but
     this property plays also a role in the diffusion of tracer equations */
  cs_property_t   *permeability;
  cs_field_t      *permea_field;     /* Related cs_field_t structure at cells */

  /* Settings related to the advection field stemming from the darcian flux */
  cs_flag_t        flux_location;   /* Indicate where the array is defined */

  /* Array defining the advection field (optional) */
  cs_real_t       *darcian_flux;

  /* Array defining the normal flux of the advection field across the domain
     boundary (optional) */
  cs_real_t       *darcian_boundary_flux;

  cs_adv_field_t  *adv_field;

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
 * Private static inline function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate
 *
 * \param[in] gw         pointer to the cs_gwf_t structure
 * \param[in] t_cur      current physical time
 * \param[in] dt_cur     current time step
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_get_time_eval(const cs_gwf_t  *const gw,
               const cs_real_t        t_cur,
               const cs_real_t        dt_cur)
{
  cs_param_time_scheme_t rt_scheme = cs_equation_get_time_scheme(gw->richards);
  cs_real_t  t_eval = 0.;
  cs_real_t  theta = -1;

  switch (rt_scheme) {

  case CS_TIME_SCHEME_STEADY:
  case CS_TIME_N_SCHEMES:
    /* Look for tracer equations */
    for (int ieq = 0; ieq < gw->n_tracers; ieq++)
      theta = fmax(theta, cs_equation_get_theta_time_val(gw->tracers[ieq]->eq));
    if (theta > 0)
      t_eval = t_cur + theta*dt_cur;
    else
      t_eval = t_cur;
    break;

  case CS_TIME_SCHEME_IMPLICIT:
  case CS_TIME_SCHEME_EXPLICIT:
  case CS_TIME_SCHEME_THETA:    /* One assumes that theta != 0 and != 1 */
  case CS_TIME_SCHEME_CRANKNICO:
    theta = cs_equation_get_theta_time_val(gw->richards);
    t_eval = t_cur + theta*dt_cur;
    break;

  } /* End of switch on time scheme for the Richards equation */

  return t_eval;
}

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the associated boundary Darcy flux for each vertex of
 *         boundary faces.
 *         Case of a vertex-based discretization.
 *
 * \param[in]  t_eval     time at which one performs the evaluation
 * \param[in]  gw         pointer to the cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_darcy_vb_flux_at_boundary(cs_real_t                      t_eval,
                                  const cs_gwf_t                *gw)
{
  cs_adv_field_t  *adv = gw->adv_field;

  if (adv->n_bdy_flux_defs > 1 ||
      adv->bdy_flux_defs[0]->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field at the boundary",
              __func__);

  cs_xdef_t  *def = adv->bdy_flux_defs[0];
  cs_xdef_array_input_t  *ai = (cs_xdef_array_input_t *)def->input;
  cs_real_t  *nflx_val = ai->values;

  if (cs_flag_test(ai->loc, cs_flag_dual_closure_byf) == false)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field at the boundary",
              __func__);

  const cs_equation_t  *richards = gw->richards;

  cs_equation_compute_boundary_diff_flux(t_eval, richards, nflx_val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update head values (pressure head or head values for laws)
 *
 * \param[in, out] gw          pointer to a cs_gwf_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      cur2prev    true or false
 */
/*----------------------------------------------------------------------------*/

static void
_update_head(cs_gwf_t                    *gw,
             const cs_cdo_quantities_t   *cdoq,
             const cs_cdo_connect_t      *connect,
             bool                         cur2prev)
{
  const cs_equation_t  *richards = gw->richards;
  assert(richards != NULL);

  cs_param_space_scheme_t r_scheme = cs_equation_get_space_scheme(richards);
  cs_field_t  *hydraulic_head = cs_equation_get_field(richards);
  cs_field_t  *pressure_head = gw->pressure_head;

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

    /* Sanity checks */
    if (pressure_head == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " The field related to the pressure head is not allocated.");

    /* Copy current field values to previous values */
    if (cur2prev)
      cs_field_current_to_previous(gw->pressure_head);

    switch (r_scheme) {

    case CS_SPACE_SCHEME_CDOVB:

#     pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
        const cs_real_t  gpot = _dp3(cdoq->vtx_coord + 3*i, gw->gravity);
        pressure_head->val[i] = hydraulic_head->val[i] - gpot;
      }

      /* Update head_in_law */
      cs_reco_pv_at_cell_centers(connect->c2v, cdoq, pressure_head->val,
                                 gw->head_in_law);
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      {
        assert(hydraulic_head->location_id ==
               cs_mesh_location_get_id_by_name("vertices"));

#       pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
          const cs_real_t  gpot = _dp3(cdoq->vtx_coord + 3*i, gw->gravity);
          pressure_head->val[i] = hydraulic_head->val[i] - gpot;
        }

        /* Update head_in_law */
        const cs_real_t  *hydraulic_head_cells =
          cs_equation_get_cell_values(richards);

#       pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {
          const cs_real_t  gpot = _dp3(cdoq->cell_centers + 3*i, gw->gravity);
          gw->head_in_law[i] = hydraulic_head_cells[i] - gpot;
        }

      }
      break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:

#     pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {
        const cs_real_t  gpot = _dp3(cdoq->cell_centers + 3*i, gw->gravity);
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
                                 gw->head_in_law);
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      {
        const cs_real_t  *hydraulic_head_cells =
          cs_equation_get_cell_values(richards);

        memcpy(gw->head_in_law, hydraulic_head_cells,
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
 * \brief  Update the advection field related to the Darcean flux
 *
 * \param[in]      t_eval     time at which one performs the evaluation
 * \param[in, out] gw         pointer to a cs_gwf_t structure
 * \param[in]      cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

static void
_update_darcy_flux(const cs_real_t              t_eval,
                   cs_gwf_t                    *gw,
                   bool                         cur2prev)
{
  const cs_equation_t  *richards = gw->richards;
  const cs_adv_field_t  *adv = gw->adv_field;

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

    /* Update the array gw->darcian_flux associated to the advection field */
    if (cs_flag_test(gw->flux_location, cs_flag_dual_face_byc)) {

      assert(gw->darcian_flux != NULL);
      if (adv->definition->type != CS_XDEF_BY_ARRAY)
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid definition of the advection field", __func__);

      cs_equation_compute_diff_flux_cellwise(richards,
                                             gw->flux_location,
                                             t_eval,
                                             gw->darcian_flux);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 2
      cs_dbg_darray_to_listing("DARCIAN_FLUX_DFbyC",
                               connect->c2e->idx[cdoq->n_cells],
                               gw->darcian_flux, 8);
#endif

      /* Set the new values of the vector field at cell centers */
      cs_advection_field_in_cells(gw->adv_field, t_eval, vel->val);


    }
    else if (cs_flag_test(gw->flux_location, cs_flag_primal_cell))
      cs_equation_compute_diff_flux_cellwise(richards,
                                             gw->flux_location,
                                             t_eval,
                                             vel->val);

    /* Update the Darcy flux at the boundary */
    _update_darcy_vb_flux_at_boundary(t_eval, gw);
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
  gw->flag = 0;
  gw->post_flag = CS_GWF_POST_DARCY_FLUX_BALANCE;

  gw->richards = NULL;
  gw->n_tracers = 0;
  gw->tracers = NULL;
  gw->finalize_tracer_setup = NULL;
  gw->add_tracer_terms = NULL;

  gw->gravity[0] = 0, gw->gravity[1] = 0, gw->gravity[2] = 0;

  gw->moisture_content = NULL;
  gw->moisture_field = NULL;

  gw->pressure_head = NULL;
  gw->head_in_law = NULL;

  gw->soil_capacity = NULL;
  gw->capacity_field = NULL;

  gw->permeability = NULL;
  gw->permea_field = NULL;

  gw->flux_location = cs_flag_dual_face_byc;
  gw->darcian_flux = NULL;
  gw->darcian_boundary_flux = NULL;
  gw->adv_field = NULL;

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
 * \param[in]      pty_type         type of permeability (iso, ortho...)
 * \param[in]      flag             flag to handle this module
 *
 * \return a pointer to a new allocated groundwater flow structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_activate(cs_property_type_t    pty_type,
                cs_flag_t             flag)
{
  cs_gwf_t  *gw = _gwf_create();

  gw->flag = flag;

  /* Create a new equation structure for Richards' equation */
  cs_equation_t  *richards = cs_equation_add("Richards",
                                             "hydraulic_head",
                                             CS_EQUATION_TYPE_GROUNDWATER,
                                             1,
                                             CS_PARAM_BC_HMG_NEUMANN);

  gw->richards = richards;

  cs_equation_param_t  *eqp = cs_equation_get_param(richards);

  /* Add an advection field related to the darcian flux stemming from the
     Richards equation */
  gw->adv_field = cs_advection_field_add("darcy_velocity",
                                         CS_ADVECTION_FIELD_GWF);

  /* Add a property related to the diffusion term of the Richards eq. */
  gw->permeability = cs_property_add("permeability", pty_type);

  /* Associate permeability to the diffusion property of the Richards eq. */
  cs_equation_add_diffusion(eqp, gw->permeability);

  /* Add a property related to the moisture content */
  gw->moisture_content = cs_property_add("moisture_content", CS_PROPERTY_ISO);

  /* Add a property related to the unsteady term of the Richards eq. */
  if (flag & CS_GWF_RICHARDS_UNSTEADY) {

    gw->soil_capacity = cs_property_add("soil_capacity", CS_PROPERTY_ISO);

    /* Associate soil_capacity to the unsteady term of the Richards eq. */
    cs_equation_add_time(eqp, gw->soil_capacity);

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

  /* darcian_flux and darcian_boundary_flux are allocated only if the related
     advection field is defined by array.
     In this case, the lifecycle is managed by the definition */

  if (gw->head_in_law != NULL)
    BFT_FREE(gw->head_in_law);

  cs_gwf_soil_free_all();

  /* Manage tracer-related members */
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

  cs_log_printf(CS_LOG_SETUP, "\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSummary of the groundwater module\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);

  /* Tracers */
  cs_log_printf(CS_LOG_SETUP,
                " <GWF/Tracer> n_tracer_equations %d\n", gw->n_tracers);

  cs_log_printf(CS_LOG_SETUP, " <GWF/Darcy flux location> %s\n",
                cs_flag_str_location(gw->flux_location));

  /* Display information on the general options */
  if (gw->flag & CS_GWF_GRAVITATION)
    cs_log_printf(CS_LOG_SETUP,
                  " <GWF/Gravitation> true -- Axis = [%.2f %.2f %.2f]\n",
                  gw->gravity[0], gw->gravity[1], gw->gravity[2]);
  else
    cs_log_printf(CS_LOG_SETUP, " <GWF/Gravitation> false\n");

  if (gw->flag & CS_GWF_ENFORCE_DIVERGENCE_FREE)
    cs_log_printf(CS_LOG_SETUP, " <GWF> Enforce the divergence-free constraint"
                  " for the Darcy flux\n");
  if (gw->flag & CS_GWF_FORCE_RICHARDS_ITERATIONS)
    cs_log_printf(CS_LOG_SETUP, " <GWF> Force to solve Richards equation"
                  " at each time step\n");
  if (gw->flag & CS_GWF_RESCALE_HEAD_TO_ZERO_MEAN_VALUE)
    cs_log_printf(CS_LOG_SETUP, " <GWF> Rescale head w.r.t zero mean value\n");

  /* Display information on the post-processing options */
  _Bool  post_capacity = (gw->post_flag & CS_GWF_POST_CAPACITY) ? true : false;
  _Bool  post_moisture = (gw->post_flag & CS_GWF_POST_MOISTURE) ? true : false;
  _Bool  post_perm = (gw->post_flag & CS_GWF_POST_PERMEABILITY) ? true : false;
  cs_log_printf(CS_LOG_SETUP, " <GWF/Post> Capacity: %s; Moisture: %s;"
                " Permeability: %s\n",
                cs_base_strtf(post_capacity), cs_base_strtf(post_moisture),
                cs_base_strtf(post_perm));

  _Bool  do_balance =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_BALANCE) ? true : false;
  _Bool  do_divergence =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_DIVERGENCE) ? true : false;
  _Bool  post_boundary =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY) ? true : false;
  cs_log_printf(CS_LOG_SETUP,
                " <GWF/Darcy_flux_post> Balance: %s; Divergence: %s;"
                " At_boundary_faces: %s\n",
                cs_base_strtf(do_balance), cs_base_strtf(do_divergence),
                cs_base_strtf(post_boundary));

  /* Soils */
  if (gw->flag & CS_GWF_SOIL_ALL_SATURATED)
    cs_log_printf(CS_LOG_SETUP, " <GWF> All soils are saturated\n");
  if (gw->flag & CS_GWF_SOIL_PROPERTY_UNSTEADY)
    cs_log_printf(CS_LOG_SETUP, " <GWF> Soil properties are unsteady\n");

  /* Detailed setup of the soil properties */
  cs_gwf_soil_log_setup();

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
  if (gw->post_flag & CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY)
    cs_advection_field_set_option(gw->adv_field,
                                  CS_ADVKEY_DEFINE_AT_BOUNDARY_FACES);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the gravity and set the gravitaty vector

 * \param[in]       gvec      values of the gravity vector
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_gravity_vector(const cs_real_3_t      gvec)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  gw->flag |= CS_GWF_GRAVITATION;
  gw->gravity[0] = gvec[0];
  gw->gravity[1] = gvec[1];
  gw->gravity[2] = gvec[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Advanced setting: indicate where the darcian flux is stored
 *
 *         cs_flag_dual_face_byc is the default setting for Vb (default space
 *         scheme) whereas cs_flag_primal_cell should be prefered for other
 *         schemes
 *
 * \param[in]  location_flag      where the flux is defined
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_darcian_flux_location(cs_flag_t      location_flag)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  gw->flux_location = location_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a particular type of unsteady advection-diffusion
 *         reaction eq.
 *         Tracer is advected thanks to the darcian velocity and
 *         diffusion/reaction parameters result from a physical modelling.
 *         Terms are activated according to the settings.
 *
 * \param[in]  eq_name    name of the tracer equation
 * \param[in]  var_name   name of the related variable
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_add_tracer(const char               *eq_name,
                  const char               *var_name)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  int  tr_id = gw->n_tracers;
  cs_gwf_tracer_t  *tracer = cs_gwf_tracer_init(tr_id,
                                                eq_name,
                                                var_name,
                                                gw->adv_field,
                                                CS_GWF_TRACER_STANDARD);

  gw->n_tracers += 1;
  BFT_REALLOC(gw->tracers, gw->n_tracers, cs_gwf_tracer_t *);
  BFT_REALLOC(gw->finalize_tracer_setup,
              gw->n_tracers, cs_gwf_tracer_setup_t *);
  BFT_REALLOC(gw->add_tracer_terms,
              gw->n_tracers, cs_gwf_tracer_add_terms_t *);

  gw->tracers[tr_id] = tracer;
  gw->finalize_tracer_setup[tr_id] = cs_gwf_tracer_standard_setup;
  gw->add_tracer_terms[tr_id] = cs_gwf_tracer_standard_add_terms;

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
cs_gwf_add_tracer_user(const char                  *eq_name,
                       const char                  *var_name,
                       cs_gwf_tracer_setup_t       *setup,
                       cs_gwf_tracer_add_terms_t   *add_terms)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  int  tr_id = gw->n_tracers;
  cs_gwf_tracer_t  *tracer = cs_gwf_tracer_init(tr_id,
                                                eq_name,
                                                var_name,
                                                gw->adv_field,
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
 *         equations defining the groundwater flow module
 *         Create new cs_field_t structures according to the setting
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_setup(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  /* Sanity checks */
  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));
  assert(gw->richards != NULL);

  const int  n_soils = cs_gwf_get_n_soils();
  if (n_soils < 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Groundwater module is activated but no soil is defined.",
              __func__);

  const bool has_previous = cs_equation_is_steady(gw->richards) ? false : true;
  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");
  const cs_param_space_scheme_t  space_scheme =
    cs_equation_get_space_scheme(gw->richards);

  /* Handle gravity effects */
  if (gw->flag & CS_GWF_GRAVITATION) {

    switch (space_scheme) {
    case CS_SPACE_SCHEME_CDOVB:
    case CS_SPACE_SCHEME_CDOVCB:
      gw->pressure_head = cs_field_create("pressure_head",
                                          field_mask,
                                          v_loc_id,
                                          1,
                                          has_previous);
      break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
      gw->pressure_head = cs_field_create("pressure_head",
                                          field_mask,
                                          c_loc_id,
                                          1,
                                          has_previous);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);
    }

    cs_field_set_key_int(gw->pressure_head, log_key, 1);
    cs_field_set_key_int(gw->pressure_head, post_key, 1);

  } /* Gravitation is activated */

  /* Detect if all soils are considered as saturated. If this not the case,
     create new fields. Check also if properties are time-dependent. */
  bool  pty_has_previous = false;
  int soil_id = 0;
  for (soil_id = 0; soil_id < n_soils; soil_id++) {

    const cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(soil_id);

    /* Is there a unique model ? */
    if (soil->model != CS_GWF_SOIL_SATURATED) {
      gw->flag |= CS_GWF_SOIL_PROPERTY_UNSTEADY;
      pty_has_previous = true;
      break;
    }

  } /* Loop on soils */

  if (soil_id == n_soils)
    gw->flag |= CS_GWF_SOIL_ALL_SATURATED;

  /* Create a moisture field attached to cells */
  int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  gw->moisture_field = cs_field_create("moisture_content",
                                       pty_mask,
                                       c_loc_id,
                                       1,   /* dimension */
                                       pty_has_previous);

  if (pty_has_previous)
    cs_field_set_key_int(gw->moisture_field, log_key, 1);
  if (gw->post_flag & CS_GWF_POST_MOISTURE)
    cs_field_set_key_int(gw->moisture_field, post_key, 1);

  if (!(gw->flag & CS_GWF_SOIL_ALL_SATURATED) ||
      gw->post_flag & CS_GWF_POST_PERMEABILITY) {

    /* Set the values for the permeability and the moisture content
       and if needed set also the value of the soil capacity */
    int  permeability_dim;
    switch (gw->permeability->type) {

    case CS_PROPERTY_ISO:
      permeability_dim = 1;
      break;
    case CS_PROPERTY_ORTHO:
      permeability_dim = 3;
      break;
    case CS_PROPERTY_ANISO:
      permeability_dim = 9;
      break;

    default:
      permeability_dim = 0;  /* avoid warning */
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of property for %s.",
                __func__, cs_property_get_name(gw->permeability));
      break;

    } /* Switch on property type */

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

    gw->capacity_field = cs_field_create("soil_capacity",
                                         pty_mask,
                                         c_loc_id,
                                         1,   /* dimension */
                                         pty_has_previous);

    cs_field_set_key_int(gw->capacity_field, log_key, 1);
    if (gw->post_flag & CS_GWF_POST_CAPACITY)
      cs_field_set_key_int(gw->capacity_field, post_key, 1);

  }

  /* Add default post-processing related to groundwater flow module */
  cs_post_add_time_mesh_dep_output(cs_gwf_extra_post, gw);
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

  int  n_soils = cs_gwf_get_n_soils();
  if (n_soils < 1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Groundwater module is activated but no soil is defined."));

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
  size_t  array_size;
  cs_flag_t  array_location;

  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  const cs_field_t  *hydraulic_head = cs_equation_get_field(gw->richards);
  const cs_param_space_scheme_t  richards_scheme =
    cs_equation_get_space_scheme(gw->richards);
  const cs_lnum_t  n_cells = connect->n_cells;

  cs_field_t  *cell_adv_field =
    cs_advection_field_get_field(gw->adv_field, CS_MESH_LOCATION_CELLS);
  assert(cell_adv_field != NULL);

  /* Set the Darcian flux (in the volume and at the boundary) */
  switch (richards_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_adjacency_t  *bf2v = connect->bf2v;

      /* Define the flux of the advection field at the boundary */
      array_size = bf2v->idx[quant->n_b_faces];
      BFT_MALLOC(gw->darcian_boundary_flux, array_size, cs_real_t);
      memset(gw->darcian_boundary_flux, 0, array_size*sizeof(cs_real_t));

      array_location = CS_FLAG_SCALAR | cs_flag_dual_closure_byf;
      cs_advection_field_def_boundary_flux_by_array(gw->adv_field,
                                                    NULL,
                                                    array_location,
                                                    gw->darcian_boundary_flux,
                                                    bf2v->idx);

      /* Define the advection field in the volume */
      if (cs_flag_test(gw->flux_location, cs_flag_dual_face_byc)) {

        /* Darcian flux settings */
        const cs_adjacency_t  *c2e = connect->c2e;

        array_size = c2e->idx[n_cells];
        BFT_MALLOC(gw->darcian_flux, array_size, cs_real_t);
        memset(gw->darcian_flux, 0, array_size*sizeof(cs_real_t));

        array_location = CS_FLAG_SCALAR | gw->flux_location;
        cs_advection_field_def_by_array(gw->adv_field,
                                        array_location,
                                        gw->darcian_flux,
                                        c2e->idx);


      }
      else if (cs_flag_test(gw->flux_location, cs_flag_primal_cell)) {

        cs_advection_field_def_by_field(gw->adv_field, cell_adv_field);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid location for defining the Darcian flux.",
                  __func__);

      /* Allocate a head array defined at cells used to update the soil
         properties */
      BFT_MALLOC(gw->head_in_law, n_cells, cs_real_t);

    }
    break;

  case CS_SPACE_SCHEME_CDOFB:   /* TODO */

    /* Set the head array defined at cells used to update the soil properties */
    if (gw->flag & CS_GWF_GRAVITATION)
      gw->head_in_law = gw->pressure_head->val;
    else
      gw->head_in_law = hydraulic_head->val;

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
                                  gw->moisture_content,
                                  gw->moisture_field);

    if (gw->permea_field != NULL) /* Fill the values of the permeability field
                                     for post-processing */

      cs_property_eval_at_cells(0, /* Should be a steady-state property */
                                gw->permeability,
                                gw->permea_field->val);

  }
  else
    cs_gwf_soil_set_by_field(gw->permeability,
                             gw->permea_field,
                             gw->moisture_content,
                             gw->moisture_field,
                             gw->soil_capacity,
                             gw->capacity_field);

  /* Store the soil id for each cell */
  cs_gwf_build_cell2soil(n_cells);

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
 * \param[in]  dt_cur     current value of the time step
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_update(const cs_mesh_t             *mesh,
              const cs_cdo_connect_t      *connect,
              const cs_cdo_quantities_t   *quant,
              const cs_time_step_t        *ts,
              double                       dt_cur,
              bool                         cur2prev)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  /* Sanity checks */
  if (gw == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Groundwater module is not allocated.", __func__);

  /* Update head */
  _update_head(gw, quant, connect, cur2prev);

  cs_real_t  time_eval = cur2prev ?
    _get_time_eval(gw, ts->t_cur, dt_cur) : ts->t_cur;

  /* Update the advection field related to the groundwater flow module */
  _update_darcy_flux(time_eval, gw, cur2prev);

  /* Update properties related to soils.
     Handle the moisture content field: do something only if the moisture
     content is set as unsteady, otherwise keep values already set */
  if (gw->flag & CS_GWF_SOIL_ALL_SATURATED) {

    /* Handle only the moisture field if this is the initialization */
    if (cur2prev == false)
      cs_property_eval_at_cells(time_eval,
                                gw->moisture_content,
                                gw->moisture_field->val);

  }
  else {

    /* Handle the permeability, the moisture content and the soil capacity */
    assert(gw->permea_field != NULL);
    if (cur2prev) {

      cs_field_current_to_previous(gw->permea_field);
      cs_field_current_to_previous(gw->moisture_field);
      if (gw->capacity_field != NULL)
        cs_field_current_to_previous(gw->capacity_field);

    }

    const int n_soils = cs_gwf_get_n_soils();
    for (int i = 0; i < n_soils; i++) {

      cs_gwf_soil_t  *soil = cs_gwf_soil_by_id(i);
      const cs_zone_t  *zone = cs_volume_zone_by_id(soil->zone_id);

      soil->update_properties(time_eval, mesh, connect, quant,
                              gw->head_in_law,
                              zone,
                              soil->input);

    } /* Loop on soils */

  } /* Not all saturated */

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
  cs_dbg_darray_to_listing("MOISTURE_CONTENT",
                           quant->n_cells,
                           gw->moisture_field->val, 8);
#endif

  /* Update the diffusivity associated to each tracer equation if needed */
  for (int i = 0; i < gw->n_tracers; i++) {

    cs_gwf_tracer_t  *tracer = gw->tracers[i];
    if (tracer->update_properties != NULL)
      tracer->update_properties(tracer, mesh, connect, quant, time_eval);

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
  cs_equation_t  *richards = gw->richards;

  double  dt_cur = 0.;  /* Useless in case of steady-state computation */

  /* Sanity check */
  assert(richards != NULL);
  assert(cs_equation_get_type(richards) == CS_EQUATION_TYPE_GROUNDWATER);

  /* Build and solve the linear system related to the Richards equations */
  if (cs_equation_is_steady(richards) ||
      gw->flag & CS_GWF_FORCE_RICHARDS_ITERATIONS) {

    if (cs_equation_uses_new_mechanism(richards))
      cs_equation_solve_steady_state(mesh, richards);

    else { /* Deprecated */

      /* Define the algebraic system */
      cs_equation_build_system(mesh, time_step, dt_cur, richards);

      /* Solve the algebraic system */
      cs_equation_solve_deprecated(richards);

    }

    /* Update the variables related to the groundwater flow system */
    cs_gwf_update(mesh, connect, cdoq, time_step, dt_cur, true);

  }

  for (int i = 0; i < gw->n_tracers; i++) {

    cs_gwf_tracer_t  *tracer = gw->tracers[i];

    if (cs_equation_is_steady(tracer->eq)) {

      if (cs_equation_uses_new_mechanism(tracer->eq))
        cs_equation_solve_steady_state(mesh, tracer->eq);

      else { /* Deprecated */

        /* Define the algebraic system */
        cs_equation_build_system(mesh, time_step, dt_cur, tracer->eq);

        /* Solve the algebraic system */
        cs_equation_solve_deprecated(tracer->eq);

      }

    } /* Solve this equation which is steady */

  } /* Loop on tracer equations */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the system related to groundwater flows module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_compute(const cs_mesh_t              *mesh,
               const cs_time_step_t         *time_step,
               double                        dt_cur,
               const cs_cdo_connect_t       *connect,
               const cs_cdo_quantities_t    *cdoq)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  cs_equation_t  *richards = gw->richards;

  /* Sanity check */
  assert(richards != NULL);
  assert(cs_equation_get_type(richards) == CS_EQUATION_TYPE_GROUNDWATER);

  /* Build and solve the linear system related to the Richards equations */
  if (!cs_equation_is_steady(richards) ||
      gw->flag & CS_GWF_FORCE_RICHARDS_ITERATIONS) {

    if (cs_equation_uses_new_mechanism(richards))
      cs_equation_solve(mesh, dt_cur, richards);

    else { /* Deprecated */

      /* Define the algebraic system */
      cs_equation_build_system(mesh, time_step, dt_cur, richards);

      /* Solve the algebraic system */
      cs_equation_solve_deprecated(richards);

    }

    /* Update the variables related to the groundwater flow system */
    cs_gwf_update(mesh, connect, cdoq, time_step, dt_cur, true);

  }

  for (int i = 0; i < gw->n_tracers; i++) {

    cs_gwf_tracer_t  *tracer = gw->tracers[i];

    if (!cs_equation_is_steady(tracer->eq)) { /* unsteady ? */

      if (cs_equation_uses_new_mechanism(tracer->eq))
        cs_equation_solve(mesh, dt_cur, tracer->eq);

      else { /* Deprecated */

        /* Define the algebraic system */
        cs_equation_build_system(mesh, time_step, dt_cur, tracer->eq);

        /* Solve the algebraic system */
        cs_equation_solve_deprecated(tracer->eq);

      }

    } /* Solve this equation which is steady */

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

  cs_gwf_std_tracer_input_t  *sti = (cs_gwf_std_tracer_input_t *)tracer->input;
  cs_real_t  int_value = 0.0;

  switch (tr_eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    {
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(tracer->eq);
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
      const cs_real_t  *v_vals = cs_equation_get_vertex_values(tracer->eq);
      const cs_real_t  *c_vals = cs_equation_get_cell_values(tracer->eq);
      const cs_adjacency_t  *c2v = connect->c2v;

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t  c_id = (z->elt_ids == NULL) ? i : z->elt_ids[i];

        /* Shares between cell and vertex unknows:
           - the cell unknown stands for 1/4 of the cell volume
           - the vertex unknown stands for 3/4 of the dual cell volume
         */
        cs_real_t  _int_value = 0.25*c_vals[c_id];
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

  /* Parallel synchronisation */
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

  const cs_lnum_t  n_b_faces = cdoq->n_b_faces;
  const cs_adv_field_t  *adv = gw->adv_field;
  assert(adv != NULL);
  const cs_field_t  *nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_t  *flux_val = (nflx == NULL) ? gw->darcian_boundary_flux : nflx->val;

  if (flux_val == NULL && n_b_faces > 0) /* No value on which operates */
    return;

  /* Define the balance by zone (using the splitting arising from the settings
     of the boundary conditions for the Richards equation) */

  const cs_equation_t  *richards = gw->richards;
  const cs_equation_param_t  *eqp = cs_equation_get_param(richards);

  _Bool  *is_counted = NULL;
  BFT_MALLOC(is_counted, n_b_faces, _Bool);
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
      cs_xdef_array_input_t  *ai = (cs_xdef_array_input_t *)_def->input;

      assert(adv->n_bdy_flux_defs == 1 && _def->type == CS_XDEF_BY_ARRAY);
      assert(cs_flag_test(ai->loc, cs_flag_dual_closure_byf) == true);
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

  _Bool  display = false;
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

  if (cs_glob_n_ranks > 1) {
#if defined(HAVE_MPI)
    cs_parall_max(1, CS_INT_TYPE, &display_flag);
    cs_parall_sum(eqp->n_bc_defs + 1, CS_REAL_TYPE, balances);
#endif
  }

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
 * \brief  Predefined post-processing output for the groundwater flow module
 *         prototype of this function is fixed since it is a function pointer
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
cs_gwf_extra_post(void                      *input,
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

  if (mesh_id == CS_POST_MESH_VOLUME) {

    if (gw->post_flag & CS_GWF_POST_DARCY_FLUX_DIVERGENCE) {

      /* Only case avalaible up to now */
      if (cs_advection_field_get_deftype(gw->adv_field) == CS_XDEF_BY_ARRAY) {

        cs_real_t  *divergence =
          cs_advection_field_divergence_at_vertices(gw->adv_field,
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
