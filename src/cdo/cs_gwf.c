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

/* Set of parameters related to the groundwater module */

struct _gwf_t {

  cs_flag_t          flag;

  /* Gravity effect */
  cs_real_3_t        gravity;

  /* Set of equations associated to this module */
  cs_equation_t     *richards;  /* "hydraulic_head" is the associated variable
                                   "permeability" is the diffusion property
                                   related to this equation */

  /* Members related to the associated tracer equations */
  int                          n_tracers;
  cs_gwf_tracer_t            **tracers;
  cs_gwf_tracer_setup_t      **finalize_tracer_setup;  // Function pointers
  cs_gwf_tracer_add_terms_t  **add_tracer_terms;       // Function pointers

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
  cs_field_t      *moisture_field;   // Related cs_field_t structure at cells

  /* Soil capacity: property attached to the unsteady term in the Richards
     equation */
  cs_property_t   *soil_capacity;
  cs_field_t      *capacity_field;   // Related cs_field_t structure at cells

  /* Permeability is the diffusion property related to Richards equation but
     this property plays also a role in the diffusion of tracer equations */
  cs_property_t   *permeability;
  cs_field_t      *permea_field;     // Related cs_field_t structure at cells

  /* Settings related to the advection field stemming from the darcian flux */
  cs_flag_t        flux_location; /* indicate where the array is defined */
  cs_real_t       *darcian_flux;  /* array defining the advection field
                                     (optional) */
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
 * Private function prototypes
 *============================================================================*/

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

  const cs_field_t  *hydraulic_head = cs_equation_get_field(richards);
  cs_field_t  *pressure_head = gw->pressure_head;

  if (gw->flag & CS_GWF_GRAVITATION) { /* Update the pressure head */

    /* Sanity checks */
    if (pressure_head == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " The field related to the pressure head is not allocated.");

    /* Copy current field values to previous values */
    if (cur2prev)
      cs_field_current_to_previous(gw->pressure_head);

    switch (cs_equation_get_space_scheme(richards)) {

    case CS_SPACE_SCHEME_CDOVB:
      assert(hydraulic_head->location_id ==
             cs_mesh_location_get_id_by_name("vertices"));

#     pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {

        const cs_real_t  gpot = _dp3(cdoq->vtx_coord + 3*i,
                                                      gw->gravity);

        pressure_head->val[i] = hydraulic_head->val[i] - gpot;

      }

      /* Update head_in_law */
      cs_reco_pv_at_cell_centers(connect->c2v,
                                 cdoq,
                                 pressure_head->val,
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
      assert(hydraulic_head->location_id ==
             cs_mesh_location_get_id_by_name("cells"));

#     pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {

        const cs_real_t  gpot = _dp3(cdoq->cell_centers + 3*i, gw->gravity);

        pressure_head->val[i] = hydraulic_head->val[i] - gpot;

      }
      break; // Nothing to do (h_head is a pointer to richards field)

    default:
      bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");

    } // Switch on space scheme

  }
  else { // No gravity effect id taken into account

    /* Update head_in_law */
    switch(cs_equation_get_space_scheme(richards)) {

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
      break; // Nothing to do

    } // Switch on the space scheme related to the Richards equation

  } /* Gravity is activated or not */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the advection field related to the Darcean flux
 *
 * \param[in, out] gw          pointer to a cs_gwf_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      ts          pointer to a cs_time_step_t structure
 * \param[in]      cur2prev    true or false
 */
/*----------------------------------------------------------------------------*/

static void
_update_darcy_velocity(cs_gwf_t                    *gw,
                       const cs_cdo_quantities_t   *cdoq,
                       const cs_cdo_connect_t      *connect,
                       const cs_time_step_t        *ts,
                       bool                         cur2prev)
{
  CS_UNUSED(connect);

  const cs_equation_t  *richards = gw->richards;
  const cs_equation_param_t  *r_eqp = cs_equation_get_param(richards);

  cs_field_t  *vel = cs_advection_field_get_field(gw->adv_field,
                                                  CS_MESH_LOCATION_CELLS);
  cs_field_t  *nflx =
    cs_advection_field_get_field(gw->adv_field,
                                 CS_MESH_LOCATION_BOUNDARY_FACES);

  /* Sanity checks */
  assert(vel != NULL && nflx != NULL);
  assert(richards != NULL);

  if (cur2prev) {
    cs_field_current_to_previous(vel);
    cs_field_current_to_previous(nflx);
  }

  /* Compute the darcian flux and the darcian velocity inside each cell */
  switch (cs_equation_get_space_scheme(richards)) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:

    /* Update the array gw->darcian_flux associated to the advection field */
    if (cs_flag_test(gw->flux_location, cs_flag_dual_face_byc)) {

      assert(gw->darcian_flux != NULL);
      cs_equation_compute_diff_flux_cellwise(richards,
                                             gw->flux_location,
                                             ts->t_cur,
                                             gw->darcian_flux);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 2
      if (cs_flag_test(gw->flux_location, cs_flag_dual_face_byc))
        cs_dbg_darray_to_listing("DARCIAN_FLUX_DFbyC",
                                 connect->c2e->idx[cdoq->n_cells],
                                 gw->darcian_flux, 8);
#endif

      /* Set the new values */
      cs_advection_field_in_cells(gw->adv_field, ts->t_cur, vel->val);

    }
    else if (cs_flag_test(gw->flux_location, cs_flag_primal_cell))
      cs_equation_compute_diff_flux_cellwise(richards,
                                             gw->flux_location,
                                             ts->t_cur,
                                             vel->val);

    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
    bft_error(__FILE__, __LINE__, 0, " TODO.");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");

  } // End of switch

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
  cs_dbg_darray_to_listing("DARCIAN_FLUX_CELL", 3*cdoq->n_cells, vel->val, 3);
#endif

  /* Update the value of the normal flux related to Darcean velocity.
     We assume a homogeneous Neumann boundary condition as a default */
  memset(nflx->val, 0, sizeof(cs_real_t)*cdoq->n_b_faces);

  for (int def_id = 0; def_id < r_eqp->n_bc_defs; def_id++) {

    const cs_xdef_t  *def = r_eqp->bc_defs[def_id];
    const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);
    assert(def->support == CS_XDEF_SUPPORT_BOUNDARY);

    if (cs_flag_test(def->meta, CS_CDO_BC_NEUMANN)) {

      switch (def->type) {

      case CS_XDEF_BY_VALUE:
        {
          const cs_real_t  *constant_val = (cs_real_t *)def->input;

          if (z->elt_ids == NULL) {
            assert(z->n_elts == cdoq->n_b_faces);
#           pragma omp parallel for if (cdoq->n_b_faces > CS_THR_MIN)
            for (cs_lnum_t i = 0; i < cdoq->n_b_faces; i++)
              nflx->val[i] = constant_val[0];
          }
          else {
#           pragma omp parallel for if (z->n_elts > CS_THR_MIN)
            for (cs_lnum_t i = 0; i < z->n_elts; i++)
              nflx->val[z->elt_ids[i]] = constant_val[0];
          }

        }
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        {
          cs_xdef_analytic_input_t  *anai =
            (cs_xdef_analytic_input_t *)def->input;

          anai->func(ts->t_cur, z->n_elts, z->elt_ids, cdoq->b_face_center,
                     false,  // compacted output ?
                     anai->input,
                     nflx->val);
        }

        break;

      default:
        bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.", __func__);
      }

    }
    else if (cs_flag_test(def->meta, CS_CDO_BC_HMG_DIRICHLET) ||
             cs_flag_test(def->meta, CS_CDO_BC_DIRICHLET)) {

      cs_nvec3_t  nvec;

      const cs_lnum_t  *bf2c = connect->f2c->ids + 2*cdoq->n_i_faces;

#     pragma omp parallel for if (z->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        const cs_lnum_t  bf_id = (z->elt_ids == NULL) ? i : z->elt_ids[i];
        const cs_lnum_t  c_id = bf2c[bf_id];
        const cs_lnum_t  f_id = cdoq->n_i_faces + bf_id;
        const cs_quant_t  pfq = cs_quant_set_face(f_id, cdoq);

        cs_nvec3(vel->val + 3*c_id, &nvec);
        nflx->val[bf_id] = nvec.meas*pfq.meas * _dp3(pfq.unitv, nvec.unitv);
      }

    }
    else if (cs_flag_test(def->meta, CS_CDO_BC_NEUMANN))
      continue; /* already handled during the initialization */

    else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of boundary condition.", __func__);

  } /* Loop on boundary conditions applied to the Richards equation */

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

  gw->flux_location = cs_flag_primal_cell;
  gw->darcian_flux = NULL;
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

  /* darcian_flux is allocated only if the advection field is defined by array.
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

  if (gw->flag & CS_GWF_GRAVITATION)
    cs_log_printf(CS_LOG_SETUP,
                  "  <GW/Gravitation> true -- Axis = [%.2f %.2f %.2f]\n",
                  gw->gravity[0], gw->gravity[1], gw->gravity[2]);
  else
    cs_log_printf(CS_LOG_SETUP, "  <GW/Gravitation> false\n");

  /* Tracers */
  cs_log_printf(CS_LOG_SETUP,
                "  <GW/Tracer> n_tracer_equations %d\n", gw->n_tracers);

  /* Soils */
  cs_gwf_soil_log_setup();

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
 *         cs_flag_primal_cell is the default setting
 *         cs_flag_dual_face_byc is a valid choice for vertex-based schemes
 *
 * \param[in]       location_flag   where the flux is defined
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

  int  n_soils = cs_gwf_get_n_soils();
  if (n_soils < 1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Groundwater module is activated but no soil is defined."));

  const bool has_previous = cs_equation_is_steady(gw->richards) ? false:true;
  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");

  /* Handle gravity effects */
  if (gw->flag & CS_GWF_GRAVITATION) {

    const cs_param_space_scheme_t  space_scheme =
      cs_equation_get_space_scheme(gw->richards);

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
      bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");
    }

    cs_field_set_key_int(gw->pressure_head, cs_field_key_id("log"), 1);
    cs_field_set_key_int(gw->pressure_head, cs_field_key_id("post_vis"), 1);

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

  } // Loop on soils

  if (soil_id == n_soils)
    gw->flag |= CS_GWF_SOIL_ALL_SATURATED;

  /* Create a moisture field attached to cells */
  int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  gw->moisture_field = cs_field_create("moisture_content",
                                       pty_mask,
                                       c_loc_id,
                                       1,   // dimension
                                       pty_has_previous);

  cs_field_set_key_int(gw->moisture_field, cs_field_key_id("log"), 1);

  if (!(gw->flag & CS_GWF_SOIL_ALL_SATURATED)) {

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
      permeability_dim = 0; // avoid warning
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of property for %s.",
                __func__, cs_property_get_name(gw->permeability));
      break;

    } // Switch on property type

    gw->permea_field = cs_field_create("permeability",
                                       pty_mask,
                                       c_loc_id,
                                       permeability_dim,   // dimension
                                       pty_has_previous);

    cs_field_set_key_int(gw->permea_field, cs_field_key_id("log"), 1);

    /* Add the post-processing of the moisture content */
    cs_field_set_key_int(gw->moisture_field, cs_field_key_id("post_vis"), 1);

    /* Create a capacity field attached to cells */
    if (gw->flag & CS_GWF_RICHARDS_UNSTEADY)
      gw->capacity_field = cs_field_create("soil_capacity",
                                           pty_mask,
                                           c_loc_id,
                                           1,   // dimension
                                           pty_has_previous);

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
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  const cs_lnum_t  n_cells = connect->n_cells;
  cs_field_t  *cell_adv_field =
    cs_advection_field_get_field(gw->adv_field, CS_MESH_LOCATION_CELLS);
  assert(cell_adv_field != NULL);

  /* Set the Darcian flux */
  if (cs_flag_test(gw->flux_location, cs_flag_dual_face_byc)) {

    /* Darcian flux settings */
    const cs_adjacency_t  *c2e = connect->c2e;

    BFT_MALLOC(gw->darcian_flux, c2e->idx[n_cells], cs_real_t);

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < c2e->idx[n_cells]; i++)
      gw->darcian_flux[i] = 0;

    /* Define and then link the advection field to each tracer equations */
    cs_flag_t  array_location = CS_FLAG_SCALAR | gw->flux_location;

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
              " Invalid location for defining the Darcian flux.");

  const cs_field_t  *hydraulic_head = cs_equation_get_field(gw->richards);
  const cs_param_space_scheme_t  richards_scheme =
    cs_equation_get_space_scheme(gw->richards);

  if (richards_scheme == CS_SPACE_SCHEME_CDOFB ||
      richards_scheme == CS_SPACE_SCHEME_HHO_P0 ||
      richards_scheme == CS_SPACE_SCHEME_HHO_P1 ||
      richards_scheme == CS_SPACE_SCHEME_HHO_P2)
    bft_error(__FILE__, __LINE__, 0,
              _(" Richards eq. is only available for vertex-based schemes."));

  /* Up to now Richards equation is only set with vertex-based schemes
     TODO: Face-based schemes */
  switch (richards_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    BFT_MALLOC(gw->head_in_law, n_cells, cs_real_t);
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
    if (gw->flag & CS_GWF_GRAVITATION)
      gw->head_in_law = gw->pressure_head->val;
    else
      gw->head_in_law = hydraulic_head->val;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid space scheme.");
    break;

  }

  /* Set permeability, moisture content and soil capacity according to the
     soil settings */
  if (gw->flag & CS_GWF_SOIL_ALL_SATURATED)
    cs_gwf_soil_set_all_saturated(gw->permeability,
                                  gw->moisture_content,
                                  gw->moisture_field);
  else
    cs_gwf_soil_set_by_field(gw->permeability,
                             gw->permea_field,
                             gw->moisture_content,
                             gw->moisture_field,
                             gw->soil_capacity,
                             gw->capacity_field);

  cs_gwf_build_cell2soil(n_cells);

  /* Loop on tracer equations */
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
              " Groundwater module is not allocated.");

  /* Update head */
  _update_head(gw, quant, connect, cur2prev);

  /* Update properties related to soils.
     Handle the moisture content field: do something only if the moisture
     content is set as unsteady, otherwise keep values already set */
  if (gw->flag & CS_GWF_SOIL_ALL_SATURATED) {

    /* Handle only the moisture field if this is the initialization */
    if (cur2prev == false)
      cs_property_eval_at_cells(ts->t_cur,
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

      soil->update_properties(mesh, connect, quant, ts,
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

  /* Update the advection field related to the groundwater flow module */
  _update_darcy_velocity(gw, quant, connect, ts, cur2prev);

  /* Update the diffusivity associated to each tracer equation if needed */
  for (int i = 0; i < gw->n_tracers; i++) {

    cs_gwf_tracer_t  *tracer = gw->tracers[i];
    if (tracer->update_properties != NULL)
      tracer->update_properties(tracer, mesh, connect, quant, ts->t_cur);

  }

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

  const int  nt_cur = time_step->nt_cur;

  if (nt_cur == 0) {

    /* Build and solve the linear system related to the Richards equations */
    if (cs_equation_is_steady(richards)) {

      /* Define the algebraic system */
      cs_equation_build_system(mesh, time_step, dt_cur, richards);

      /* Solve the algebraic system */
      cs_equation_solve(richards);

      /* Update the variables related to the groundwater flow system */
      cs_gwf_update(mesh, connect, cdoq, time_step, true);

    }

    for (int i = 0; i < gw->n_tracers; i++) {

      cs_gwf_tracer_t  *tracer = gw->tracers[i];

      if (cs_equation_is_steady(tracer->eq)) {

        /* Define the algebraic system */
        cs_equation_build_system(mesh, time_step, dt_cur, tracer->eq);

        /* Solve the algebraic system */
        cs_equation_solve(tracer->eq);

      } /* Solve this equation which is steady */

    } /* Loop on tracer equations */

  }
  else { /* nt_cur > 0 */

    /* Build and solve the linear system related to the Richards equations */
    if (!cs_equation_is_steady(richards)) {

      /* Define the algebraic system */
      if (cs_equation_needs_build(richards)) // unsteady ?
        cs_equation_build_system(mesh, time_step, dt_cur, richards);

      /* Solve the algebraic system */
      cs_equation_solve(richards);

      /* Update the variables related to the groundwater flow system */
      cs_gwf_update(mesh, connect, cdoq, time_step, true);

    }

    for (int i = 0; i < gw->n_tracers; i++) {

      cs_gwf_tracer_t  *tracer = gw->tracers[i];

      if (!cs_equation_is_steady(tracer->eq)) { // unsteady ?

        /* Define the algebraic system */
        if (cs_equation_needs_build(tracer->eq))
          cs_equation_build_system(mesh, time_step, dt_cur, tracer->eq);

        /* Solve the algebraic system */
        cs_equation_solve(tracer->eq);

      } /* Solve this equation which is steady */

    } /* Loop on tracer equations */

  } /* nt_cur > 0 */

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

  if (mesh_id == CS_POST_MESH_BOUNDARY) {

    const cs_field_t  *nflx =
      cs_advection_field_get_field(gw->adv_field,
                                   CS_MESH_LOCATION_BOUNDARY_FACES);

    if (nflx == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Null pointer encounter\n", __func__);

    cs_log_printf(CS_LOG_DEFAULT,
                  " Balance of the Darcy flux across the domain boundary\n");

    cs_real_t balance = 0;
    for (cs_lnum_t  i = 0; i < n_b_faces; i++)
      balance += nflx->val[i];
    cs_real_t  default_balance = balance;

    for (int def_id = 0; def_id < gw->adv_field->n_bdy_flux_defs; def_id++) {

      const cs_xdef_t  *def = gw->adv_field->bdy_flux_defs[def_id];
      const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);
      assert(def->support == CS_XDEF_SUPPORT_BOUNDARY);

      if (z->elt_ids == NULL)
        break; /* Nothing to do (balance = default_balance) */
      else {

        balance = 0;
        for (cs_lnum_t i = 0; i < z->n_elts; i++)
          balance += nflx->val[z->elt_ids[i]];

        cs_log_printf(CS_LOG_DEFAULT, " %32s: % -5.3e\n", z->name, balance);
        default_balance -= balance;

      }
    } /* Loop on boundary definitions for the Darcy flux */

    if (gw->adv_field->n_bdy_flux_defs == 0)
      cs_log_printf(CS_LOG_DEFAULT,
                    " %32s: % -5.3e\n", "Whole boundary", default_balance);

    else
      cs_log_printf(CS_LOG_DEFAULT, " %32s: % -5.3e\n",
                    "Remaining boundary", default_balance);

  } /* boundary mesh_id */

  if (mesh_id == CS_POST_MESH_VOLUME) {

    /* Only case avalaible up to now */
    if (cs_advection_field_get_deftype(gw->adv_field) == CS_XDEF_BY_ARRAY) {

      cs_real_t  *divergence = cs_advection_field_divergence(gw->adv_field,
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

  } /* volume mesh id */

}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
