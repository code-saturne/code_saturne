/*============================================================================
 * Handle the solidification module with CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

#include "cs_navsto_system.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_solidification.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_SOLIDIFICATION_DBG     0

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Set of parameters related to the solidification module */

struct _solidification_t {

  cs_flag_t        model;         /* Modelling for the solidifcation module */
  cs_flag_t        options;       /* Flag dedicated to general options to handle
                                   * the solidification module*/
  cs_flag_t        post_flag;     /* Flag dedicated to the post-processing
                                   * of the solidifcation module */

  /* Equations associated to this module: There are as many equations as the
   * number of tracers to consider (in this case, a tracer is a component of an
   * alloy)
   */

  int              n_alloy_tracers;
  cs_equation_t  **tracer_equations;

  /* Fields associated to this module */
  cs_field_t      *temperature;

  /* Physical parameters to specify the law of variation of the liquid fraction
   * with respect to the temperature
   *
   * gl(T) = 1 if T > t_liquidus and gl(T) = 0 if T < t_solidus
   * Otherwise:
   * gl(T) = (T - t_solidus)/(t_liquidus - t_solidus)
   */
  cs_real_t       t_solidus;
  cs_real_t       t_liquidus;

  cs_field_t      *gl_field;   /* field storing the values of the liquid at
                                  each cell */
  cs_property_t   *gl;         /* liquid fraction property */

  /* A reaction term is introduced in the momentum equation in order to penalize
   * the velocity and tends to forcing_coef when the liquid fraction tends to 0
   * F(u) = forcing_coef * (1- gl)^2/(gl^3 + forcing_eps) * u
   */

  cs_real_t        forcing_eps;       /* Set to 1e-3 by default */
  cs_real_t        forcing_coef;      /* Set to 1600 by default */
  cs_real_t       *forcing_mom_array; /* values of the forcing reaction
                                         coefficient in each cell */
  cs_property_t   *forcing_mom;

  cs_real_t        latent_heat;



};

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_module[] =
  " Stop execution.\n"
  " The structure related to the solidifcation module is empty.\n"
  " Please check your settings.\n";

static cs_solidification_t  *cs_solidification_structure = NULL;

/*============================================================================
 * Private static inline function prototypes
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the structure dedicated to the management of the
 *         solidifcation module
 *
 * \return a pointer to a new allocated cs_solidification_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_solidification_t *
_solidification_create(void)
{
  cs_solidification_t  *solid = NULL;

  BFT_MALLOC(solid, 1, cs_solidification_t);

  /* Default initialization */
  solid->model = 0;
  solid->options = 0;
  solid->post_flag = 0;

  /* Fields related to this module */
  solid->temperature = NULL;

  /* Quantities related to the liquid fraction */
  solid->gl = NULL;
  solid->gl_field = NULL;

  solid->t_solidus = 0.;
  solid->t_liquidus = 0.;

  /* Quantities/structure related to the forcing term treated as a reaction term
     in the momentum equation */
  solid->forcing_mom = NULL;
  solid->forcing_mom_array = NULL;
  solid->forcing_coef = 1600.;
  solid->forcing_eps = 1e-3;

  /* Related alloys equations */
  solid->n_alloy_tracers = 0;
  solid->tracer_equations = NULL;

  /* Other physical parameters */
  solid->latent_heat = 0.;

  return solid;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update/initialize the liquid fraction and its related quantities
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

static void
_update_liquid_fraction(const cs_mesh_t             *mesh,
                        const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_time_step_t        *ts,
                        bool                         cur2prev)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(ts);

  cs_solidification_t  *solid = cs_solidification_structure;
  assert(solid->temperature != NULL);

  if (cur2prev)
    cs_field_current_to_previous(solid->gl_field);

  cs_real_t  *g_l = solid->gl_field->val;
  cs_real_t  *temp = solid->temperature->val;
  assert(temp != NULL);

  /* 1./(t_liquidus - t_solidus) */
  const cs_real_t inv_delta_threshold =
    1./(solid->t_liquidus - solid->t_solidus);

  cs_lnum_t  n_solid_cells = 0;
  cs_lnum_t  *solid_cells = NULL;

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Update the liquid fraction */
    cs_real_t  gl = 0;
    if (temp[c_id] < solid->t_solidus)
      g_l[c_id] = 0;
    else if (temp[c_id] > solid->t_liquidus)
      g_l[c_id] = gl = 1;
    else
      g_l[c_id] = gl = (temp[c_id] - solid->t_solidus)*inv_delta_threshold;

    /* Update the forcing term treated as a reaction term in the momentum eq. */
    const cs_real_t  glm1 = 1 - gl, gl3 = gl*gl*gl;
    solid->forcing_mom_array[c_id] = solid->forcing_coef *
      glm1*glm1/(gl3 + solid->forcing_eps);

    if (gl > 0.99)
      n_solid_cells++;

  } /* Loop on cells */

  if (n_solid_cells > 0) {

    cs_equation_t  *mom_eq = cs_navsto_system_get_momentum_eq();
    cs_equation_param_t  *mom_eqp = cs_equation_get_param(mom_eq);

    BFT_MALLOC(solid_cells, n_solid_cells, cs_lnum_t);

    cs_lnum_t  ii = 0;
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      if (g_l[c_id] > 0.99)
        solid_cells[ii++] = c_id;

    assert(ii == n_solid_cells);
    cs_real_t  zero_velocity[3] = {0, 0, 0};
    cs_equation_enforce_by_cell_selection(mom_eqp,
                                          n_solid_cells, solid_cells,
                                          zero_velocity,
                                          NULL);

    BFT_FREE(solid_cells);

  } /* Set the enforcement of the velocity for solid cells */

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if solidification module is activated
 */
/*----------------------------------------------------------------------------*/

bool
cs_solidification_is_activated(void)
{
  if (cs_solidification_structure == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the solidification module
 *
 * \param[in]  model            type of modelling
 * \param[in]  options          flag to handle optional parameters
 * \param[in]  post_flag        predefined post-processings
 * \param[in]  boundaries       pointer to the domain boundaries
 * \param[in]  algo_coupling    algorithm used for solving the NavSto system
 * \param[in]  ns_option        option flag for the Navier-Stokes system
 * \param[in]  ns_post_flag     predefined post-processings for Navier-Stokes
 *
 * \return a pointer to a new allocated solidification structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_t *
cs_solidification_activate(cs_solidification_model_t      model,
                           cs_flag_t                      options,
                           cs_flag_t                      post_flag,
                           const cs_boundary_t           *boundaries,
                           cs_navsto_param_coupling_t     algo_coupling,
                           cs_flag_t                      ns_option,
                           cs_flag_t                      ns_post_flag)
{
  if (model < 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid modelling. Model = %d\n", __func__, model);

  /* Allocate an empty structure */
  cs_solidification_t  *solid = _solidification_create();

  /* Set members of the structure according to the given settings */
  solid->model = model;
  solid->options = options;
  solid->post_flag = post_flag;

  /* Activate the thermal module */
  cs_flag_t  thm_num = 0, thm_post = 0;
  cs_flag_t  thm_model = CS_THERMAL_MODEL_NAVSTO_VELOCITY;

  cs_thermal_system_t  *thm = cs_thermal_system_activate(thm_model,
                                                         thm_num,
                                                         thm_post);
  cs_navsto_param_model_t  ns_model = 0;
  if (model & CS_SOLIDIFICATION_MODEL_STOKES)
    ns_model |= CS_NAVSTO_MODEL_STOKES;
  else if (model & CS_SOLIDIFICATION_MODEL_NAVIER_STOKES)
    ns_model |= CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES;

  /* Activate the Navier-Stokes module */
  cs_navsto_system_activate(boundaries,
                            ns_model,
                            algo_coupling,
                            ns_option,
                            ns_post_flag);

  /* Add properties related to this module */
  solid->forcing_mom = cs_property_add("forcing_momentum_coef",
                                       CS_PROPERTY_ISO);
  solid->gl = cs_property_add("liquid_fraction", CS_PROPERTY_ISO);

  /* Set the global pointer */
  cs_solidification_structure = solid;

  return solid;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the solidification module
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_solidification_t *
cs_solidification_destroy_all(void)
{
  if (cs_solidification_structure == NULL)
    return NULL;

  cs_solidification_t  *solid = cs_solidification_structure;

  /* The lifecycle of properties, equations and fields is not managed by
   * the current structure.
   * Free only what is owned by this structure */

  if (solid->n_alloy_tracers > 0)
    BFT_FREE(solid->tracer_equations);

  BFT_FREE(solid->forcing_mom_array);
  BFT_FREE(solid);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the main physical parameters dedicated to this module
 *
 * \param[in]  t_solidus      solidus temperature (in K)
 * \param[in]  t_liquidus     liquidus temperatur (in K)
 * \param[in]  latent_heat    latent heat
 * \param[in]  forcing_eps    epsilon used in penalization term to division by
 *                            zero
 * \param[in]  forcing_coef   (< 0) coefficient in the reaction term to reduce
 *                            the velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_parameters(cs_real_t          t_solidus,
                                 cs_real_t          t_liquidus,
                                 cs_real_t          latent_heat,
                                 cs_real_t          forcing_eps,
                                 cs_real_t          forcing_coef)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  /* Sanity checks */
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  solid->t_solidus = t_solidus;
  solid->t_liquidus = t_liquidus;
  solid->latent_heat = latent_heat;
  solid->forcing_coef = forcing_coef;
  solid->forcing_eps = forcing_eps;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup equations/properties related to the Solidification module
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_init_setup(void)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  /* Sanity checks */
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_CDO;
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");

  /* Add a field for the liquid fraction */
  solid->gl_field = cs_field_create("liquid_fraction",
                                    field_mask,
                                    c_loc_id,
                                    1,
                                    true); /* has_previous */

  cs_field_set_key_int(solid->gl_field, log_key, 1);
  cs_field_set_key_int(solid->gl_field, post_key, 1);

  /* Add a reaction term to the momentum equation */
  cs_equation_t  *mom_eq = cs_navsto_system_get_momentum_eq();
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(mom_eq);
  assert(mom_eqp != NULL);

  cs_equation_add_reaction(mom_eqp, solid->forcing_mom);


  /* Add default post-processing related to the Maxwell module */
  cs_post_add_time_mesh_dep_output(cs_solidification_extra_post, solid);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup stage for equations related to the solidification
 *         module
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_finalize_setup(const cs_cdo_connect_t       *connect,
                                 const cs_cdo_quantities_t    *quant)
{
  CS_UNUSED(connect);

  cs_solidification_t  *solid = cs_solidification_structure;

  /* Sanity checks */
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  /* Retrieve the temperature field */
  solid->temperature = cs_field_by_name("temperature");
  assert(solid->temperature != NULL);

  /* Define the liquid fraction */
  cs_property_def_by_field(solid->gl, solid->gl_field);

  /* Define the forcing term acting as a reaction term in the momentum equation
     This term is related to the liquid fraction */
  BFT_MALLOC(solid->forcing_mom_array, quant->n_cells, cs_real_t);
  memset(solid->forcing_mom_array, 0, quant->n_cells*sizeof(cs_real_t));

  cs_property_def_by_array(solid->forcing_mom,
                           cs_flag_primal_cell,
                           solid->forcing_mom_array,
                           false, /* definition is owner ? */
                           NULL); /* no index */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the solidification module in the log file dedicated to
 *         the setup
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_log_setup(void)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the solidification module\n");
  cs_log_printf(CS_LOG_SETUP, "%s\n", h1_sep);

  cs_log_printf(CS_LOG_SETUP, "  * Solidification | Model:");
  cs_log_printf(CS_LOG_SETUP, "\n");
  cs_log_printf(CS_LOG_SETUP, "  * Solidification | Number of alloys: %d",
                solid->n_alloy_tracers);
  if (solid->n_alloy_tracers > 0) {

    cs_log_printf(CS_LOG_SETUP, "  * Solidification | Alloy equations:");
    for (int  i = 0; i < solid->n_alloy_tracers; i++) {
      const char  *eqname = cs_equation_get_name(solid->tracer_equations[i]);
      if (eqname != NULL)
        cs_log_printf(CS_LOG_SETUP, " %s;", eqname);
    }

  } /* n_alloy_tracers > 0 */


}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve equations related to the solidification module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_compute(const cs_mesh_t              *mesh,
                          const cs_time_step_t         *time_step,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  /* Sanity checks */
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  /* Add equations to be solved at each time step */

  /* Update fields and properties which are related to solved variables */
  cs_solidification_update(mesh, connect, quant, time_step,
                           true); /* operate current to previous ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update/initialize the solidification module according to the
 *         settings
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_update(const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const cs_time_step_t        *ts,
                         bool                         cur2prev)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  /* Sanity checks */
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  /* Compute the liquid fraction and the related momentum reaction term */
  _update_liquid_fraction(mesh, connect, quant, ts, cur2prev);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the solidification module
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_extra_op(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL)
    return;

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the solidification module.
 *         Prototype of this function is fixed since it is a function pointer
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
cs_solidification_extra_post(void                      *input,
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
  CS_UNUSED(mesh_id);
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

  cs_solidification_t  *solid = (cs_solidification_t *)input;

  /* TODO */
  CS_UNUSED(solid);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
