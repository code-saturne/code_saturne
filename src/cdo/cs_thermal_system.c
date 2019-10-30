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

  /* Equation associated to this module */
  /* ---------------------------------- */

  cs_equation_t  *thermal_eq;

  /* Properties associated to this module */
  /* ------------------------------------ */

  cs_property_t  *lambda;        /* Thermal conductivity */
  cs_property_t  *cp;            /* Heat capacity */
  cs_property_t  *rho;           /* Mass density */
  cs_property_t  *kappa;         /* Thermal diffusivity */

  /* value of the thermal diffusivity in each cell (allocated if the
   * related property is defined by array) */

  cs_real_t      *kappa_array;

  /* Fields associated to this module */
  /* -------------------------------- */

  cs_field_t     *temperature;
  cs_field_t     *enthalpy;
  cs_field_t     *total_energy;

  /* Additional members */
  /* ------------------ */

  cs_real_t       ref_temperature;
  cs_real_t       thermal_dilatation_coef;

  /* N.B.: Other reference values for properties are stored within each
   * dedicated structure */

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the thermal diffusivity for a given list of cells
 *         Case of an anisotropic thermal conductivity
 *         This function fits the generic prototype of cs_xdef_cell_eval_t
 *
 * \param[in]      n_elts    number of elements to consider
 * \param[in]      elt_ids   list of element ids
 * \param[in]      compact   indirection for output (true or false)
 * \param[in]      mesh      pointer to a cs_mesh_t structure
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in]      input     pointer to an input structure cast on-the_fly
 * \param[in, out] result    array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_eval_aniso_kappa(cs_lnum_t                    n_elts,
                  const cs_lnum_t              elt_ids[],
                  bool                         compact,
                  const cs_mesh_t             *mesh,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *quant,
                  cs_real_t                    t_eval,
                  void                        *input,
                  cs_real_t                   *result)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_thermal_system_t  *thm = (cs_thermal_system_t *)input;
  assert(thm != NULL); /* Sanity check */

  /* rho and cp are assumed to be uniform */
  const cs_real_t  rho = cs_property_get_cell_value(0, t_eval, thm->rho);
  const cs_real_t  cp = cs_property_get_cell_value(0, t_eval, thm->cp);
  const cs_real_t  ov_rhocp = 1./(rho*cp);

  if (elt_ids == NULL) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_real_t  *kappa_c = result + 9*i;
      cs_property_get_cell_tensor(i, t_eval, thm->lambda, false,
                                  (cs_real_t (*)[3])kappa_c);
      for (int k = 0; k < 9; k++)
        kappa_c[k] *= ov_rhocp;

    } /* Loop on all cells */

  }
  else { /* Loop on a selection of cells */

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  c_id = elt_ids[i];
      const cs_lnum_t  id = compact ? i : c_id;

      cs_real_t  *kappa_c = result + 9*id;
      cs_property_get_cell_tensor(c_id, t_eval, thm->lambda, false,
                                  (cs_real_t (*)[3])kappa_c);
      for (int k = 0; k < 9; k++)
        kappa_c[k] *= ov_rhocp;

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the coefficient appearing in time-dependent term of the
 *         simulation of tracer equations
 *         This function fits the generic prototype of cs_xdef_cell_eval_cw_t
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval   time at which one performs the evaluation
 * \param[in]      input    pointer to an input structure cast on-the_fly
 * \param[in, out] result   array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static inline void
_eval_aniso_kappa_cw(const cs_cell_mesh_t    *cm,
                     cs_real_t                t_eval,
                     void                    *input,
                     cs_real_t               *result)
{
  cs_thermal_system_t  *thm = (cs_thermal_system_t *)input;
  assert(thm != NULL); /* Sanity check */

  /* rho and cp are assumed to be uniform */
  const cs_real_t  rho = cs_property_get_cell_value(0, t_eval, thm->rho);
  const cs_real_t  cp = cs_property_get_cell_value(0, t_eval, thm->cp);
  const cs_real_t  ov_rhocp = 1./(rho*cp);

  cs_property_get_cell_tensor(cm->c_id, t_eval, thm->lambda, false,
                              (cs_real_t (*)[3])result);
  for (int k = 0; k < 9; k++)
    result[k] *= ov_rhocp;
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
  cs_thermal_system_t  *thm = NULL;

  BFT_MALLOC(thm, 1, cs_thermal_system_t);

  /* Flags */
  thm->model = 0;
  thm->numeric = 0;
  thm->post = 0;

  /* Equation */
  thm->thermal_eq = NULL;

  /* Properties */
  thm->lambda = NULL;
  thm->cp = NULL;
  thm->rho = NULL;
  thm->kappa = NULL;
  thm->kappa_array = NULL;

  /* Fields */
  thm->temperature = NULL;
  thm->enthalpy = NULL;
  thm->total_energy = NULL;

  /* Other members */
  thm->ref_temperature = 0.;
  thm->thermal_dilatation_coef = 0.;

  return thm;
}

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
cs_thermal_system_activate(cs_thermal_model_type_t    model,
                           cs_flag_t                  numeric,
                           cs_flag_t                  post)
{
  cs_thermal_system_t  *thm = NULL;
  if (cs_thermal_system == NULL)
    thm = _init_thermal_system();
  else
    thm = cs_thermal_system;    /* Previously allocated when setting the
                                   reference temperature for instance */

  /* Set flags */
  thm->model = model;
  thm->numeric = numeric;
  thm->post = post;

  /* Add equation related to this module */
  cs_equation_t  *eq = NULL;
  if (model & CS_THERMAL_MODEL_WITH_ENTHALPY)
    eq = cs_equation_add(CS_THERMAL_EQNAME,
                         "enthalpy",
                         CS_EQUATION_TYPE_THERMAL,
                         1,
                         CS_PARAM_BC_HMG_NEUMANN);
  else if (model & CS_THERMAL_MODEL_WITH_TOTAL_ENERGY)
    eq = cs_equation_add(CS_THERMAL_EQNAME,
                         "total_energy",
                         CS_EQUATION_TYPE_THERMAL,
                         1,
                         CS_PARAM_BC_HMG_NEUMANN);
  else
    eq = cs_equation_add(CS_THERMAL_EQNAME,
                         "temperature",
                         CS_EQUATION_TYPE_THERMAL,
                         1,
                         CS_PARAM_BC_HMG_NEUMANN);

  thm->thermal_eq = eq;

  /* Default numerical settings */
  /* -------------------------- */

  cs_equation_param_t  *eqp = cs_equation_get_param(eq);

  /* Linear algebra default settings */
  cs_equation_set_param(eqp, CS_EQKEY_SOLVER_FAMILY, "cs");
  cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "amg");
  cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
  cs_equation_set_param(eqp, CS_EQKEY_ITSOL_EPS, "1e-8");
  cs_equation_set_param(eqp, CS_EQKEY_ITSOL_RESNORM_TYPE, "rhs");

  /* Set a space discretization by default */
  if (thm->model & CS_THERMAL_MODEL_NAVSTO_VELOCITY) {

    cs_equation_set_param(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "ocs");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_COEF, "sushi");

  }
  else {

    cs_equation_set_param(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vb");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "bubble");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_COEF, "frac23");

  }

  /* Define properties */
  cs_property_type_t  pty_type = CS_PROPERTY_ISO;

  /* Mass density */
  thm->rho = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
  if (thm->rho == NULL)
    thm->rho = cs_property_add(CS_PROPERTY_MASS_DENSITY, pty_type);

  /* Thermal capacity */
  thm->cp = cs_property_add(CS_THERMAL_CP_NAME, pty_type);

  /* Thermal conductivity */
  if (model & CS_THERMAL_MODEL_ANISOTROPIC_CONDUCTIVITY)
    pty_type = CS_PROPERTY_ANISO;
  thm->lambda = cs_property_add(CS_THERMAL_LAMBDA_NAME, pty_type);

  if (thm->model & CS_THERMAL_MODEL_WITH_THERMAL_DIFFUSIVITY) {

    thm->kappa = cs_property_add(CS_THERMAL_DIFF_NAME, pty_type);

    /* Link the conductivity to the diffusion term of this equation */
    cs_equation_add_diffusion(eqp, thm->kappa);

  }
  else /* Link the conductivity to the diffusion term of this equation */
    cs_equation_add_diffusion(eqp, thm->lambda);

  /* Add an advection term  */
  if (thm->model & CS_THERMAL_MODEL_NAVSTO_VELOCITY) {

    cs_equation_add_advection(eqp,
                              cs_advection_field_by_name("velocity_field"));

    cs_equation_set_param(eqp, CS_EQKEY_ADV_FORMULATION, "non_conservative");
    cs_equation_set_param(eqp, CS_EQKEY_ADV_SCHEME, "upwind");

  }

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

  if (thm->kappa_array != NULL)
    BFT_FREE(thm->kappa_array);

  /* Equations, fields and properties related to the thermal system are
   * destroyed elsewhere in a specific stage. The lifecycle of these structures
   * are not managed by cs_thermal_system_t
   */

  BFT_FREE(thm);
  cs_thermal_system = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a reference temperature
 *
 * \param[in]  temp0     reference temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_set_reference_temperature(cs_real_t  temp0)
{
  cs_thermal_system_t  *thm = NULL;
  if (cs_thermal_system == NULL)
    thm = _init_thermal_system();
  else
    thm = cs_thermal_system;    /* Previously allocated when activating
                                   the thermal module for instance */

  thm->ref_temperature = temp0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the thermal dilatation coefficient
 *
 * \param[in]  beta0     reference value of the thermal dilatation coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_system_set_thermal_dilatation_coef(cs_real_t   beta0)
{
  cs_thermal_system_t  *thm = NULL;
  if (cs_thermal_system == NULL)
    thm = _init_thermal_system();
  else
    thm = cs_thermal_system;    /* Previously allocated when activating
                                   the thermal module for instance */

  thm->thermal_dilatation_coef = beta0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a structure to compute the Boussinesq source term
 *
 * \param[in]  gravity    gravity vector
 * \param[in]  rho0       reference value for the mass density
 *
 * \return a pointer to a new allocated \ref cs_source_term_boussinesq_t
 */
/*----------------------------------------------------------------------------*/

cs_source_term_boussinesq_t *
cs_thermal_system_add_boussinesq_source_term(const cs_real_t   *gravity,
                                             cs_real_t          rho0)
{
  cs_thermal_system_t  *thm = cs_thermal_system;

  /* Sanity checks */
  assert(gravity != NULL);
  if (thm == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));
  if (thm->temperature == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: No temperature field allocated.",
              __func__);

  cs_source_term_boussinesq_t  *bq_st = NULL;
  BFT_MALLOC(bq_st, 1, cs_source_term_boussinesq_t);

  bq_st->g[0] = gravity[0];
  bq_st->g[1] = gravity[1];
  bq_st->g[2] = gravity[2];
  bq_st->rho0 = rho0;
  bq_st->beta = thm->thermal_dilatation_coef;
  bq_st->var0 = thm->ref_temperature;
  bq_st->var = thm->temperature->val;

  return bq_st;
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

  if ((thm->model & CS_THERMAL_MODEL_WITH_ENTHALPY) ||
      (thm->model & CS_THERMAL_MODEL_WITH_TOTAL_ENERGY)) {

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

  cs_thermal_system_t  *thm = cs_thermal_system;
  cs_real_t  t_cur = time_step->t_cur;

  if (thm == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));

  if (thm->temperature == NULL)
    thm->temperature = cs_field_by_name("temperature");

  if (thm->model & CS_THERMAL_MODEL_WITH_THERMAL_DIFFUSIVITY) {

    /* The usage of thermal diffusivity relies on a constant mass density and
       on a constant thermal capacity */
    assert(thm->rho->n_definitions == 0 ||
           (thm->rho->n_definitions = 1 &&
            thm->rho->defs[0]->type == CS_XDEF_BY_VALUE));
    assert(thm->cp->n_definitions == 0 ||
           (thm->cp->n_definitions = 1 &&
            thm->cp->defs[0]->type == CS_XDEF_BY_VALUE));

    if (thm->model & CS_THERMAL_MODEL_ANISOTROPIC_CONDUCTIVITY) {

      if (cs_property_is_uniform(thm->lambda)) {

        cs_real_t  tens[3][3];
        /* c_id = 0, do_inversion = false */
        cs_property_get_cell_tensor(0, t_cur, thm->lambda, false, tens);
        cs_real_t  rho = cs_property_get_cell_value(0, t_cur, thm->rho);
        cs_real_t  cp = cs_property_get_cell_value(0, t_cur, thm->cp);

        const cs_real_t  inv_rhocp = 1./(rho*cp);
        for (int ki = 0; ki < 3; ki++)
          for (int kj = 0; kj < 3; kj++)
            tens[ki][kj] *= inv_rhocp;

        cs_property_def_aniso_by_value(thm->kappa, NULL, tens);

      }
      else
        cs_property_def_by_func(thm->kappa,
                                NULL, /* = all cells */
                                thm,
                                _eval_aniso_kappa,
                                _eval_aniso_kappa_cw);

    }
    else {

      /* Choose c_id = 0 since one assumes there is at least one cell */
      cs_real_t  rho = cs_property_get_cell_value(0, t_cur, thm->rho);
      cs_real_t  cp = cs_property_get_cell_value(0, t_cur, thm->cp);
      const cs_real_t  inv_rhocp = 1./(rho*cp);

      /* The thermal conductivity klambda is isotropic */
      if (cs_property_is_uniform(thm->lambda)) {

        /* c_id = 0 */
        cs_real_t  lambda = cs_property_get_cell_value(0, t_cur, thm->lambda);

        cs_property_def_iso_by_value(thm->kappa, NULL, lambda*inv_rhocp);

      }
      else {

        BFT_MALLOC(thm->kappa_array, quant->n_cells, cs_real_t);

        cs_property_eval_at_cells(t_cur, thm->lambda, thm->kappa_array);

#       pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
          thm->kappa_array[c_id] *= inv_rhocp;

        cs_property_def_by_array(thm->kappa,
                                 cs_flag_primal_cell,
                                 thm->kappa_array,
                                 false, /* = not owner */
                                 NULL); /* = no index */

      } /* Lambda uniform ? */

    } /* Anisotropic ? */

  } /* Need to define a thermal diffusivity */

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
  CS_UNUSED(cur2prev);

  cs_thermal_system_t  *thm = cs_thermal_system;

  if (thm == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_thm));

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
  if (thm->model & CS_THERMAL_MODEL_NAVSTO_VELOCITY)
    cs_log_printf(CS_LOG_SETUP, " + Navsto advection");
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
