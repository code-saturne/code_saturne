/*============================================================================
 * Functions to handle the thermal module with CDO schemes
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
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include "base/cs_mem.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_field_default.h"
#include "pprt/cs_physical_model.h"
#include "cdo/cs_reco.h"
#include "base/cs_syr_coupling.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cdo/cs_thermal_system.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 *  \file cs_thermal_system.h
 *
 *  \brief  Functions to handle the cs_thermal_system_t structure.
 *          This module can be used stand alone or linked with another module
 *          such as Navier-Stokes, groundwater flows or Maxwell...
 *  The temperature field is automatically defined when this module is
 *  activated.
 *
 *  Moreover, one considers according to the modelling
 *  rho    the volumetric mass (mass density)
 *  Cp     the heat capacity
 *  lambda the heat conductivity
 */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_THERMAL_SYSTEM_DBG  0

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_thm[] =
  " Stop execution. The structure related to the thermal system is"
  " empty.\n Please check your settings.\n";

static cs_thermal_system_t *cs_thermal_system = nullptr;

static cs_real_t *_cht_robin_vals = nullptr;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Is the thermal system solved with the temperature as variable
 */
/*----------------------------------------------------------------------------*/

static inline bool
_is_solved_with_temperature(void)
{
  if (cs_thermal_system == nullptr)
    return false;

  else {

    if (cs_thermal_system->model & CS_THERMAL_MODEL_USE_TEMPERATURE)
      return true;
    if (cs_thermal_system->model & CS_THERMAL_MODEL_USE_ENTHALPY)
      return false;
    if (cs_thermal_system->model & CS_THERMAL_MODEL_USE_TOTAL_ENERGY)
      return false;

  }

  return true; /* This is the default choice */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a cs_thermal_system_t structure
 *
 * \return a pointer to a newly allocated \ref cs_thermal_system_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_thermal_system_t *
_init_thermal_system(void)
{
  cs_thermal_system_t *thm = nullptr;

  CS_MALLOC(thm, 1, cs_thermal_system_t);

  /* Flags */

  thm->model = 0;
  thm->numeric = 0;
  thm->post = 0;

  /* Equation */

  thm->thermal_eq = nullptr;

  /* Properties */

  thm->lambda            = nullptr;
  thm->cp                = nullptr;
  thm->rho               = nullptr;
  thm->unsteady_property = nullptr;
  thm->kappa             = nullptr;
  thm->kappa_array       = nullptr;

  /* Fields */

  thm->temperature  = nullptr;
  thm->enthalpy     = nullptr;
  thm->total_energy = nullptr;

  return thm;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the value of the reference temperature associated to a
 *        thermal system.
 *
 * \return the value of the reference temperature
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_thermal_system_get_reference_temperature(void)
{
  if (cs_thermal_system == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              " %s: A reference temperature is requested but no thermal"
              " system is activated.\n"
              " Please check your settings.", __func__);

  return cs_thermal_system->ref_temperature;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the value of the reference temperature associated to the
 *        thermal system.
 *
 * \param[in]  ref     value of the reference temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_set_reference_temperature(cs_real_t    ref)
{
  if (cs_thermal_system == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              " %s: A reference temperature is requested but no thermal"
              " system is activated.\n"
              " Please check your settings.", __func__);

  cs_thermal_system->ref_temperature = ref;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the model flag related to a thermal system
 *
 * \return a flag
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_thermal_system_get_model(void)
{
  if (cs_thermal_system == nullptr)
    return 0;

  return cs_thermal_system->model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Does the thermal system rely on the advection field associated to
 *        the Navier-Stokes equations ?
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_thermal_system_needs_navsto(void)
{
  if (cs_thermal_system == nullptr)
    return false;

  if (cs_thermal_system->model & CS_THERMAL_MODEL_NAVSTO_ADVECTION)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the resolution of the thermal system has been activated
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_thermal_system_is_activated(void)
{
  if (cs_thermal_system == nullptr)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the thermal system
 *
 * \param[in] model     model flag related to the thermal system
 * \param[in] numeric   (optional) numerical flag settings
 * \param[in] post      (optional) post-processing flag settings
 *
 * \return a pointer to a new allocated cs_thermal_system_t structure
 */
/*----------------------------------------------------------------------------*/

cs_thermal_system_t *
cs_thermal_system_activate(cs_thermal_model_type_t    model,
                           cs_flag_t                  numeric,
                           cs_flag_t                  post)
{
  cs_thermal_system_t *thm = nullptr;
  if (cs_thermal_system == nullptr)
    thm = _init_thermal_system();
  else
    thm = cs_thermal_system;    /* Previously allocated when setting the
                                   reference temperature for instance */

  /* Set the physical model type */

  if (cs_glob_physical_model_flag[CS_SOLIDIFICATION] < 1)
    cs_glob_physical_model_flag[CS_HEAT_TRANSFER] = 1;

  /* Set flags */

  thm->model = model;
  thm->numeric = numeric;
  thm->post = post;

  bool has_time = true;
  if (model & CS_THERMAL_MODEL_STEADY)
    has_time = false;

  /* Define or retrieve properties related to this module */
  /* ---------------------------------------------------- */

  /* Mass density */

  thm->rho = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
  if (thm->rho == nullptr)
    thm->rho = cs_property_add(CS_PROPERTY_MASS_DENSITY, CS_PROPERTY_ISO);

  /* Thermal capacity */

  thm->cp = cs_property_by_name(CS_THERMAL_CP_NAME);
  if (thm->cp == nullptr)
    thm->cp = cs_property_add(CS_THERMAL_CP_NAME, CS_PROPERTY_ISO);

  /* Thermal conductivity */

  cs_property_type_t  pty_type = CS_PROPERTY_ISO;
  if (model & CS_THERMAL_MODEL_ANISOTROPIC_CONDUCTIVITY)
    pty_type = CS_PROPERTY_ANISO;
  thm->lambda = cs_property_add(CS_THERMAL_LAMBDA_NAME, pty_type);

  /* Add the associated equation related to this module with respect to
   * the settings */

  cs_equation_t       *eq  = nullptr;
  cs_equation_param_t *eqp = nullptr;

  if (model & CS_THERMAL_MODEL_USE_ENTHALPY) {

    eq = cs_equation_add(CS_THERMAL_EQNAME,
                         "enthalpy",
                         CS_EQUATION_TYPE_THERMAL,
                         1,
                         CS_BC_SYMMETRY);


  }
  else if (model & CS_THERMAL_MODEL_USE_TOTAL_ENERGY) {

    eq = cs_equation_add(CS_THERMAL_EQNAME,
                         "total_energy",
                         CS_EQUATION_TYPE_THERMAL,
                         1,
                         CS_BC_SYMMETRY);

    /* TODO */
    bft_error(__FILE__, __LINE__, 0, " %s: Not yet fully available.\n",
              __func__);

  }
  else { /* Default settings: use temperature as variable */

    thm->model |= CS_THERMAL_MODEL_USE_TEMPERATURE;
    eq = cs_equation_add(CS_THERMAL_EQNAME,
                         "temperature",
                         CS_EQUATION_TYPE_THERMAL,
                         1,
                         CS_BC_SYMMETRY);

    /* Always add a diffusion term */

    eqp = cs_equation_get_param(eq);
    cs_equation_add_diffusion(eqp, thm->lambda);

    if (has_time) {

      thm->unsteady_property = cs_property_add_as_product("rho.cp",
                                                          thm->rho, thm->cp);
      cs_equation_add_time(eqp, thm->unsteady_property);

    }

  } /* Use temperature as variable */

  /* Default numerical settings */
  /* -------------------------- */

  thm->thermal_eq = eq;

  /* Add an advection term  */

  if (thm->model & CS_THERMAL_MODEL_NAVSTO_ADVECTION) {

    cs_equation_add_advection(eqp,
                              cs_advection_field_by_name("mass_flux"));

    /* Set the space discretization by default. One should be consistent with
       the space scheme used for solving the Navier-Stokes system. Either with
       the legacy FV scheme or the CDO-Fb scheme, one needs a CDO Fb scheme. */

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "ocs");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_COEF, "sushi");

    if (thm->model & CS_THERMAL_MODEL_USE_TEMPERATURE) {

      /* The formulation used for the thermal equation with temperature as main
         unknown needs a non-conservative formulation of the advective term */

      cs_equation_add_advection_scaling_property(eqp, thm->cp);
      cs_equation_param_set(eqp, CS_EQKEY_ADV_FORMULATION, "non_conservative");

    }

    cs_equation_param_set(eqp, CS_EQKEY_ADV_SCHEME, "upwind");

  }
  else { /* Stand-alone i.e. not associated to the Navier--Stokes system */

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vb");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "bubble");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_COEF, "frac23");

  }

  /* Linear algebra default settings */

  cs_equation_param_set(eqp, CS_EQKEY_SOLVER_RESNORM_TYPE, "filtered");

  /* Set and return pointer */

  cs_thermal_system = thm;

  return cs_thermal_system;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the thermal system
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_destroy(void)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    return;

  if (thm->kappa_array != nullptr)
    CS_FREE(thm->kappa_array);

  /* Equations, fields and properties related to the thermal system are
   * destroyed elsewhere in a specific stage. The lifecycle of these structures
   * are not managed by cs_thermal_system_t
   */

  CS_FREE(thm);
  cs_thermal_system = nullptr;

  /* Deallocate array used for coupling if used */
  if (_cht_robin_vals != nullptr)
    CS_FREE(_cht_robin_vals);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the main equation related to the thermal system
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_thermal_system_get_equation(void)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    return nullptr;
  else
    return thm->thermal_eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the current temperature at face values
 *
 * \return the pointer to the array of face values.
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_thermal_system_get_face_temperature(void)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    return nullptr;

  else {

    if (_is_solved_with_temperature())
      return cs_equation_get_face_values(thm->thermal_eq, false);

    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: The thermal is not solved with the temperature.\n"
                " Conversion and interpolation are not done automatically.",
                __func__);

  }

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the thermal system
 *         At this stage, numerical settings should be completely determined
 *         but connectivity and geometrical information is not yet available.
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_init_setup(void)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));

  cs_param_space_scheme_t  space_scheme =
    cs_equation_get_space_scheme(thm->thermal_eq);

  int  location_support = CS_MESH_LOCATION_NONE; /* Not set */
  switch (space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    location_support = CS_MESH_LOCATION_VERTICES;
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
  case CS_SPACE_SCHEME_HHO_P1:
  case CS_SPACE_SCHEME_HHO_P2:
    location_support = CS_MESH_LOCATION_CELLS;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space scheme for the thermal system.", __func__);
    break;

  }

  /* Prepare parameters for the creation of fields */

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_CDO;
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");

  bool  has_previous = true;
  if (thm->model & CS_THERMAL_MODEL_STEADY)
    has_previous = false;

  if ((thm->model & CS_THERMAL_MODEL_USE_ENTHALPY) ||
      (thm->model & CS_THERMAL_MODEL_USE_TOTAL_ENERGY)) {

    /* Temperature field is always created. Since the variable solved is either
       the "enthalpy" or the "total_energy", one creates a field for the
       temperature */

    thm->temperature = cs_field_create("temperature",
                                       field_mask,
                                       location_support,
                                       1,
                                       has_previous);

    cs_field_set_key_int(thm->temperature, log_key, 1);
    cs_field_set_key_int(thm->temperature, post_key, 1);

  }

  if (thm->post & CS_THERMAL_POST_ENTHALPY) {

    thm->enthalpy =  cs_field_find_or_create("enthalpy",
                                             field_mask,
                                             location_support,
                                             1,
                                             has_previous);

    cs_field_set_key_int(thm->enthalpy, log_key, 1);
    cs_field_set_key_int(thm->enthalpy, post_key, 1);

  } /* enthalpy */

  /* Add the variable field */

  if (has_previous)
    cs_equation_predefined_create_field(1, thm->thermal_eq);
  else
    cs_equation_predefined_create_field(0, thm->thermal_eq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last step of the setup of the thermal system
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_finalize_setup(const cs_cdo_connect_t     *connect,
                                 const cs_cdo_quantities_t  *quant,
                                 const cs_time_step_t       *time_step)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(time_step);

  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));

  if (thm->temperature == nullptr)
    thm->temperature = cs_field_by_name("temperature");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve a steady-state thermal system
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_compute_steady_state(const cs_mesh_t              *mesh,
                                       const cs_cdo_connect_t       *connect,
                                       const cs_cdo_quantities_t    *quant,
                                       const cs_time_step_t         *time_step)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));
  assert(cs_equation_uses_new_mechanism(thm->thermal_eq));

  cs_equation_solve_steady_state(mesh, thm->thermal_eq);

  /* Update fields and properties which are related to the evolution of the
     variable solved in thermal_eq */

  cs_thermal_system_update(mesh, connect, quant, time_step,
                           true); /* operate current to previous ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the thermal system
 *
 * \param[in]  cur2prev   true="current to previous" operation is performed
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_compute(bool                          cur2prev,
                          const cs_mesh_t              *mesh,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant,
                          const cs_time_step_t         *time_step)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));
  assert(cs_equation_uses_new_mechanism(thm->thermal_eq));

  /* CHT code coupling (if needed) */
  if (cs_syr_coupling_n_couplings() > 0) {

    cs_lnum_t *cpl_ids = nullptr;
    CS_MALLOC(cpl_ids, mesh->n_b_faces, cs_lnum_t);

    cs_real_t *t_bnd = nullptr;
    cs_real_t *hf_bnd = nullptr;
    cs_real_t *tf_bnd = nullptr;
    CS_MALLOC(t_bnd, mesh->n_b_faces, cs_real_t);
    CS_MALLOC(hf_bnd, mesh->n_b_faces, cs_real_t);
    CS_MALLOC(tf_bnd, mesh->n_b_faces, cs_real_t);

    for (int icpl = 0; icpl < cs_syr_coupling_n_couplings(); icpl++) {
      /* Get coupled elements for current coupling */
      const cs_lnum_t n_cpl_elts = cs_syr_coupling_n_elts(icpl, 0);
      cs_syr_coupling_elt_ids(icpl, 0, cpl_ids);

      /* CDO numbering uses global numbering (interior, then boundary) */
      for (cs_lnum_t e_id = 0; e_id < n_cpl_elts; e_id++)
        cpl_ids[e_id] += mesh->n_i_faces;

      /* Project vertex based values to face based values */
      cs_real_t *t_vals = thm->temperature->val;
      cs_reco_scalar_v2f(n_cpl_elts,
                         cpl_ids,
                         connect,
                         quant,
                         t_vals,
                         true,
                         t_bnd);

      for (cs_lnum_t e_id = 0; e_id < n_cpl_elts; e_id++)
        cpl_ids[e_id] -= mesh->n_i_faces;

      // Send boundary temperature
      cs_syr_coupling_send_tsolid(icpl, 0, t_bnd);

      // Recieve Tf,Hf values
      cs_syr_coupling_recv_tf_hf(icpl, 0, nullptr, tf_bnd, hf_bnd);


      // Apply Tf, Hf to boundary condition defintion (xdef)
      for (cs_lnum_t e_id = 0; e_id < n_cpl_elts; e_id++) {
        cs_lnum_t face_id = cpl_ids[e_id];
        cs_real_t *_face_vals = _cht_robin_vals + 3*face_id;
        _face_vals[0] = hf_bnd[e_id];
        _face_vals[1] = tf_bnd[e_id];
        _face_vals[2] = 0.;
      }

    }

    CS_FREE(t_bnd);
    CS_FREE(tf_bnd);
    CS_FREE(hf_bnd);
    CS_FREE(cpl_ids);
  }

  if (!(thm->model & CS_THERMAL_MODEL_STEADY))
    cs_equation_solve(cur2prev, mesh, thm->thermal_eq);

  /* Update fields and properties which are related to the evolution of the
     variable solved in thermal_eq */

  cs_thermal_system_update(mesh, connect, quant, time_step, cur2prev);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set an initial values for all quantities related to this module
 *         This is done after the setup step.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_init_values(const cs_mesh_t             *mesh,
                              const cs_cdo_connect_t      *connect,
                              const cs_cdo_quantities_t   *quant,
                              const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(ts);
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    return;

  /* The value of the variable field is initialize during the initialization
   * of variable field associated to an equation */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the quantities related to the thermal module
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_update(const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const cs_time_step_t        *ts,
                         bool                         cur2prev)
{
  CS_UNUSED(mesh);
  CS_UNUSED(ts);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(cur2prev);

  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    return;

  /* Update the thermal diffusivity if requested */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the Navier-Stokes system
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_extra_op(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *cdoq)
{
  CS_UNUSED(connect);
  CS_UNUSED(cdoq);

  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the thermal system.
 *         The prototype of this function is fixed since it is a function
 *         pointer defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_thermal_system_t structure)
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
cs_thermal_system_extra_post(void                      *input,
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

  if (input == nullptr)
    return;

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main options related to cs_thermal_system_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_log_setup(void)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == nullptr)
    return;

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the thermal module\n");
  cs_log_printf(CS_LOG_SETUP, "%s\n", cs_sep_h1);

  cs_log_printf(CS_LOG_SETUP, "  * Thermal | Model:");
  if (thm->model & CS_THERMAL_MODEL_STEADY)
    cs_log_printf(CS_LOG_SETUP, " Steady-state");
  if (thm->model & CS_THERMAL_MODEL_NAVSTO_ADVECTION)
    cs_log_printf(CS_LOG_SETUP, " + Navsto advection");
  if (thm->model & CS_THERMAL_MODEL_ANISOTROPIC_CONDUCTIVITY)
    cs_log_printf(CS_LOG_SETUP, " + Anistropic conductivity");
  cs_log_printf(CS_LOG_SETUP, "\n");

  cs_log_printf(CS_LOG_SETUP,
                "  * Thermal | Equation solved with the variable");
  if (thm->model & CS_THERMAL_MODEL_USE_ENTHALPY)
    cs_log_printf(CS_LOG_SETUP, " Enthalpy\n");
  else if (thm->model & CS_THERMAL_MODEL_USE_TOTAL_ENERGY)
    cs_log_printf(CS_LOG_SETUP, " Total energy\n");
  else if (thm->model & CS_THERMAL_MODEL_USE_TEMPERATURE)
    cs_log_printf(CS_LOG_SETUP, " Temperature (Kelvin)\n");
  else
    cs_log_printf(CS_LOG_SETUP, " Unknown variable!\n");

  if (thm->post & CS_THERMAL_POST_ENTHALPY)
    cs_log_printf(CS_LOG_SETUP, "  * Thermal | Post: Enthalpy\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize structures and/or boundary conditions used for CHT coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_cht_boundary_conditions_setup(void)
{
  if (cs_syr_coupling_n_couplings() > 0) {

    /* Initialize CHT array */
    if (_cht_robin_vals != nullptr)
      CS_FREE(_cht_robin_vals);
    CS_MALLOC(_cht_robin_vals, 3 * cs_glob_mesh->n_b_faces, cs_real_t);
    cs_array_real_set_scalar(3 * cs_glob_mesh->n_b_faces,
                             1.e20,
                             _cht_robin_vals);

    /* Set boundary condition on zones */
    cs_equation_param_t *eqp = cs_equation_param_by_name(CS_THERMAL_EQNAME);

    for (int z_id = 0; z_id < cs_boundary_zone_n_zones(); z_id++) {
      const cs_zone_t *z = cs_boundary_zone_by_id(z_id);
      if (cs_syr_coupling_is_bnd_zone_coupled(z)) {
        bft_printf("Setting a coupled (CHT) Robin boundary condition for zone %s.\n",
                   z->name);
        cs_equation_add_bc_by_array(eqp,
                                    CS_PARAM_BC_ROBIN,
                                    z->name,
                                    cs_flag_boundary_face,
                                    _cht_robin_vals,
                                    false,  /* not owner */
                                    true); /* not full length -> we use a subset of array */
      }
    }

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
