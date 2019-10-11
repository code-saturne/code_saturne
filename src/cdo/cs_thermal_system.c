/*============================================================================
 * Routines to handle the thermal module with CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_thermal_system.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 *  \file cs_thermal_system.c
 *
 *  \brief  Routines to handle the cs_thermal_system_t structure
 *          This module can be used stand alone or linked with other module
 *          such as Navier-Stokes (for a Boussinesq approximation for instance)
 *          the GroundWater Flow module or the Maxwell module...
 */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_THERMAL_SYSTEM_DBG  0

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Set of parameters related to the thermal module */

struct _thermal_system_t {

  cs_flag_t   model;                  /* Choice of modelling */
  cs_flag_t   numeric;                /* General numerical options */
  cs_flag_t   post;                   /* Post-processing options */

  /* Properties associated to this module */
  cs_property_t  *lambda;        /* Thermal conductivity */
  cs_property_t  *cp;            /* Heat capacity */
  cs_property_t  *rho;           /* Mass density */
  cs_property_t  *rho_cp;        /* Thermal diffusivity */

  /* Equation associated to this module */
  cs_equation_t  *thermal_eq;

  /* Fields associated to this module */
  cs_field_t     *temperature;
  cs_field_t     *enthalpy;

  /* Arrays associated to this module */

};

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_thm[] =
  " Stop execution. The structure related to the thermal system is"
  " empty.\n Please check your settings.\n";

static cs_thermal_system_t  *cs_thermal_system = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/



/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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
  if (cs_thermal_system == NULL)
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
cs_thermal_system_activate(cs_flag_t         model,
                           cs_flag_t         numeric,
                           cs_flag_t         post)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  BFT_MALLOC(thm, 1, cs_thermal_system_t);

  /* Set main flags */
  thm->model = model;
  thm->numeric = numeric;
  thm->post = post;

  /* Set a space discretization by default if not set */
  if (thm->numeric == 0)
    thm->numeric = CS_THERMAL_CDOVB; /* default choice */
  else {

    if (!(numeric & CS_THERMAL_CDOVB)  &&
        !(numeric & CS_THERMAL_CDOVCB) &&
        !(numeric & CS_THERMAL_CDOFB) )
      thm->numeric = CS_THERMAL_CDOVB; /* default choice */

  }

  cs_equation_t  *eq = NULL;
  if (model & CS_THERMAL_MODEL_WITH_ENTHALPY)
    eq = cs_equation_add(CS_THERMAL_EQNAME,
                         "enthalpy",
                         CS_EQUATION_TYPE_THERMAL,
                         1,
                         CS_PARAM_BC_HMG_NEUMANN);
  else
    eq = cs_equation_add(CS_THERMAL_EQNAME,
                         "temperature",
                         CS_EQUATION_TYPE_THERMAL,
                         1,
                         CS_PARAM_BC_HMG_NEUMANN);

  cs_equation_param_t  *eqp = cs_equation_get_param(eq);

  if (thm->numeric & CS_THERMAL_CDOVB) {

    /* Default numerical settings: the linear system should be symmetric */
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "bubble");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_COEF, "frac23");

  }
  else if (thm->numeric & CS_THERMAL_CDOVCB) {

    cs_equation_set_param(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vcb");

  }
  else if (thm->numeric & CS_THERMAL_CDOFB) {

    cs_equation_set_param(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "ocs");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_COEF, "sushi");

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid space discretization for the thermal"
              " module", __func__);

  /* Linear algebra default settings */
  cs_equation_set_param(eqp, CS_EQKEY_SOLVER_FAMILY, "cs");
  cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "amg");
  cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
  cs_equation_set_param(eqp, CS_EQKEY_ITSOL_EPS, "1e-8");
  cs_equation_set_param(eqp, CS_EQKEY_ITSOL_RESNORM_TYPE, "rhs");

  /* Define properties */
  cs_property_type_t  pty_type = CS_PROPERTY_ISO;

  if (thm->model & CS_THERMAL_MODEL_STEADY) {

    thm->cp = NULL;
    thm->rho = NULL;
    thm->rho_cp = NULL;

  }
  else {

    thm->cp = cs_property_add(CS_THERMAL_CP_NAME, pty_type);

    thm->rho = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
    if (thm->rho == NULL)
      thm->rho = cs_property_add(CS_PROPERTY_MASS_DENSITY, pty_type);

    thm->rho_cp = cs_property_add(CS_THERMAL_DIFF_NAME, pty_type);

    /* Link the rho_cp property to the unsteady term of this equation */
    cs_equation_add_time(eqp, thm->rho_cp);
  }

  if (model & CS_THERMAL_MODEL_ANISOTROPIC_CONDUCTIVITY)
    pty_type = CS_PROPERTY_ANISO;

  thm->lambda = cs_property_add(CS_THERMAL_LAMBDA_NAME, pty_type);

  /* Link the conductivity to the diffusion term of this equation */
  cs_equation_add_diffusion(eqp, thm->lambda);

  /* Initialize equation and field pointers */
  thm->temperature = NULL;
  thm->enthalpy = NULL;
  thm->thermal_eq = eq;

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

  if (thm == NULL)
    return;

  /* Equations, fields and properties related to the thermal system are
   * destroyed elsewhere in a specific stage. The lifecycle of these structures
   * are not managed by cs_thermal_system_t
   */

  BFT_FREE(thm);
  cs_thermal_system = NULL;
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

  if (thm == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));

  /* Prepare parameters for the creation of fields */
  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_CDO;
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");

  bool  has_previous = true;
  if (thm->model & CS_THERMAL_MODEL_STEADY)
    has_previous = false;

  int  location_support = CS_MESH_LOCATION_VERTICES;
  if (thm->numeric & CS_THERMAL_CDOFB)
    location_support = CS_MESH_LOCATION_CELLS;

  if (thm->model & CS_THERMAL_MODEL_WITH_ENTHALPY) {

    /* The variable solved is "enthalpy". So, one creates a field for
       the temperature */
    thm->temperature = cs_field_create("temperature",
                                       field_mask,
                                       location_support,
                                       1,
                                       has_previous);

    cs_field_set_key_int(thm->temperature, log_key, 1);
    cs_field_set_key_int(thm->temperature, post_key, 1);

  }
  else { /* The thermal equation is solved using the temperature as
            variable */

    if (thm->post & CS_THERMAL_POST_ENTHALPY) {

      thm->enthalpy =  cs_field_create("enthalpy",
                                       field_mask,
                                       location_support,
                                       1,
                                       has_previous);

      cs_field_set_key_int(thm->enthalpy, log_key, 1);
      cs_field_set_key_int(thm->enthalpy, post_key, 1);

    } /* enthalpy */

  }


}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last step of the setup of the thermal system
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_finalize_setup(const cs_cdo_connect_t     *connect,
                                 const cs_cdo_quantities_t  *quant)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve a steady-state thermal system
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  time_step  pointer to a cs_time_step_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_compute_steady_state(const cs_mesh_t              *mesh,
                                       const cs_time_step_t         *time_step,
                                       const cs_cdo_connect_t       *connect,
                                       const cs_cdo_quantities_t    *quant)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));
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
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  time_step  pointer to a cs_time_step_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_compute(const cs_mesh_t              *mesh,
                          const cs_time_step_t         *time_step,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));
  assert(cs_equation_uses_new_mechanism(thm->thermal_eq));

  if (!(thm->model & CS_THERMAL_MODEL_STEADY))
    cs_equation_solve(mesh, thm->thermal_eq);

  /* Update fields and properties which are related to the evolution of the
     variable solved in thermal_eq */
  cs_thermal_system_update(mesh, connect, quant, time_step,
                           true); /* operate current to previous ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update/initialize the thermal module according to the settings
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

  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));

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

  if (thm == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));

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

  if (input == NULL)
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

  if (thm == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the thermal module\n");
  cs_log_printf(CS_LOG_SETUP, "%s\n", h1_sep);

  cs_log_printf(CS_LOG_SETUP, "  * Thermal | Model:");
  if (thm->model & CS_THERMAL_MODEL_STEADY)
    cs_log_printf(CS_LOG_SETUP, " Steady-state");
  if (thm->model & CS_THERMAL_MODEL_WITH_ENTHALPY)
    cs_log_printf(CS_LOG_SETUP, " + Enthalpy-based");
  if (thm->model & CS_THERMAL_MODEL_ANISOTROPIC_CONDUCTIVITY)
    cs_log_printf(CS_LOG_SETUP, " + Anistropic conductivity");
  cs_log_printf(CS_LOG_SETUP, "\n");

  if (thm->post & CS_THERMAL_POST_ENTHALPY)
    cs_log_printf(CS_LOG_SETUP, "  * Thermal | Post: Enthalpy\n");

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
