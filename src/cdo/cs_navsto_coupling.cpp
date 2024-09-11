/*============================================================================
 * Functions and structures to handle the settings of different
 * velocity/pressure coupling algorithms
 *  - Artificial Compressibility algorithm
 *  - Monolithic algorithm
 *  - Projection algorithm
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_array.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_navsto_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_navsto_coupling.c

  \brief  Functions to handle structures used as a context when solving the
          Navier-Stokes equations.
          Structures are cast on-the-fly according to the type of coupling.

          Functions to handle the settings of coupling algorithms are dedicated
          to one of the following algorithm
          - Artificial Compressibility algorithm
          - Monolithic algorithm
          - Projection algorithm
 */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_NAVSTO_COUPLING_DBG 0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the \ref cs_equation_param_t structure related to the
 *         momentum equation according to the type of coupling
 *
 * \param[in]  nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]  context   pointer to a coupling context structure
 *
 * \return a pointer to the corresponding \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_navsto_coupling_get_momentum_eqp(const cs_navsto_param_t    *nsp,
                                    void                       *context)
{
  cs_equation_param_t *mom_eqp = nullptr;

  if (nsp == nullptr)
    return mom_eqp;

  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    {
      cs_navsto_ac_t  *nsc = (cs_navsto_ac_t *)context;
      mom_eqp = cs_equation_get_param(nsc->momentum);
    }
    break;

  case CS_NAVSTO_COUPLING_MONOLITHIC:
    {
      cs_navsto_monolithic_t  *nsc = (cs_navsto_monolithic_t *)context;
      mom_eqp = cs_equation_get_param(nsc->momentum);
    }
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
    {
      cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;
      mom_eqp = cs_equation_get_param(nsc->prediction);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid coupling algorithm\n", __func__);
    break;

  }

  return mom_eqp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using an Artificial Compressibility approach
 *
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_create_context(cs_param_bc_type_t    bc,
                            cs_navsto_param_t    *nsp)
{
  cs_navsto_ac_t *nsc = nullptr;

  BFT_MALLOC(nsc, 1, cs_navsto_ac_t);

  nsc->momentum = cs_equation_add("momentum",
                                  "velocity",
                                  CS_EQUATION_TYPE_NAVSTO,
                                  3,
                                  bc);

  /* Additional property */

  nsc->zeta = cs_property_add("graddiv_coef", CS_PROPERTY_ISO);

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);

  /* Space scheme settings (default) */

  cs_equation_param_set(mom_eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");
  cs_equation_param_set(mom_eqp, CS_EQKEY_HODGE_DIFF_COEF, "sushi");

  /* Set the default solver settings */

  if (nsp->model == CS_NAVSTO_MODEL_STOKES)
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "cg");
  else
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "gcr");

  return nsc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to an Artificial Compressibility
 *         approach
 *
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a nullptr pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_free_context(void     *context)
{
  cs_navsto_ac_t  *nsc = (cs_navsto_ac_t *)context;

  BFT_FREE(nsc);

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when an
 *         Artificial Compressibility algorithm is used to coupled the system.
 *         No mesh information is available at this stage.
 *
 * \param[in]      nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in]      adv_field   pointer to a cs_adv_field_t structure
 * \param[in, out] context     pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_ac_init_setup(const cs_navsto_param_t    *nsp,
                        cs_adv_field_t             *adv_field,
                        void                       *context)
{
  cs_navsto_ac_t  *nsc = (cs_navsto_ac_t *)context;

  assert(nsp != nullptr && nsc != nullptr);

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);

  /* Navier-Stokes parameters induce numerical settings for the related
     equations */

  cs_navsto_param_transfer(nsp, mom_eqp);

  /* Link the time property to the momentum equation */

  if (!cs_navsto_param_is_steady(nsp))
    cs_equation_add_time(mom_eqp, nsp->mass_density);

  /* Add advection term in case of CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES
   * CS_NAVSTO_MODEL_OSEEN: Nothing to do since the Oseen field is set by the
   * user via cs_navsto_add_oseen_field() */

  if (nsp->model & CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES)
    cs_equation_add_advection(mom_eqp, adv_field);

  /* All considered models needs a viscous term */

  cs_equation_add_diffusion(mom_eqp, nsp->tot_viscosity);

  /* Add the variable field (Keep the previous state) */

  cs_equation_predefined_create_field(1, nsc->momentum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when an
 *         Artificial Compressibility algorithm is used to coupled the system.
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_ac_last_setup(const cs_navsto_param_t     *nsp,
                        void                        *context)
{
  cs_navsto_ac_t  *nsc = (cs_navsto_ac_t *)context;

  assert(nsp != nullptr && nsc != nullptr);

  /* Avoid no definition of the zeta coefficient */

  if (nsc->zeta->n_definitions == 0) {

    cs_property_def_iso_by_value(nsc->zeta, nullptr, 1);

    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS,
                  "%s: Add a unitary scaling for the grad-div term.\n"
                  "%s: No definition given.",
                  __func__, __func__);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to the \ref cs_equation_t structure related to
 *         the momentum equation in case of artificial compressibility coupling
 *
 * \param[in] context  pointer to a context structure cast on-the-fly
 *
 * \return a pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_navsto_ac_get_momentum_eq(void       *context)
{
  if (context == nullptr)
    return nullptr;

  cs_navsto_ac_t  *nsc = (cs_navsto_ac_t *)context;

  return nsc->momentum;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using a monolithic approach
 *
 * \param[in]      bc     default \ref cs_param_bc_type_t for the equation
 * \param[in, out] nsp    pointer to a \ref cs_navsto_param_t structure
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_monolithic_create_context(cs_param_bc_type_t    bc,
                                    cs_navsto_param_t    *nsp)
{
  cs_navsto_monolithic_t *nsc = nullptr;

  BFT_MALLOC(nsc, 1, cs_navsto_monolithic_t);

  /* Add an equation for the momentum conservation */

  nsc->momentum = cs_equation_add("momentum",
                                  "velocity",
                                  CS_EQUATION_TYPE_NAVSTO,
                                  3,
                                  bc);

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);

  /* Space scheme settings (default) */

  cs_equation_param_set(mom_eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");
  cs_equation_param_set(mom_eqp, CS_EQKEY_HODGE_DIFF_COEF, "gcr");

  /* Forcing steady state in order to avoid inconsistencies */

  if (nsp->model_flag & CS_NAVSTO_MODEL_STEADY)
    mom_eqp->time_scheme = CS_TIME_SCHEME_STEADY;

  /* Set the name of the saddle-point problem to solve */

  cs_param_saddle_set_name("NavSto", mom_eqp->saddle_param);

  /* Associate the cs_param_sles_t structure related to the (1,1)-block */

  cs_param_saddle_set_block11_sles_param(mom_eqp->saddle_param,
                                         mom_eqp->sles_param);

  /* Solver settings: Only the linear algebra settings related to the momentum
   * equation. */

  if (nsp->model ==  CS_NAVSTO_MODEL_STOKES) {

    /* The linear system to solve should be symmetric */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "gcr");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_PRECOND, "sgs");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SCHUR_APPROX, "mass_scaled");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER_RESTART, "40");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_MAX_ITER, "100");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_RTOL, "1e-8");

    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "fcg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_MAX_ITER, "50");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_RTOL, "1e-1");

#if defined(HAVE_PETSC)
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_FAMILY, "petsc");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND_BLOCK_TYPE, "diag");

    if (cs_param_sles_hypre_from_petsc())
      cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "boomer");
    else
      cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "gamg");
#else   /* Not HAVE_PETSC */
#if defined(HAVE_HYPRE)
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "boomer");
#else  /* Neither HYPRE or PETSc */
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "k_cycle");
#endif  /* HAVE_HYPRE */
#endif  /* HAVE_PETSC */

  }
  else {

#if defined(HAVE_MUMPS)
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "alu");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_AUGMENT_SCALING, "100");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "mumps");
#else
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "gcr");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_PRECOND, "sgs");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SCHUR_APPROX, "lumped_inv");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER_RESTART, "40");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_MAX_ITER, "100");

    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "jacobi");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "gcr");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_MAX_ITER, "50");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_RTOL, "1e-1");
#endif

  }

  return nsc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to a monolithic approach
 *
 * \param[in, out] context     pointer to a context structure cast on-the-fly
 *
 * \return a nullptr pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_monolithic_free_context(void             *context)
{
  cs_navsto_monolithic_t  *nsc = (cs_navsto_monolithic_t *)context;

  BFT_FREE(nsc);

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a monolithic
 *         algorithm is used to coupled the system.
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp        pointer to a \ref cs_navsto_param_t structure
 * \param[in]      adv_field  pointer to a cs_adv_field_t structure
 * \param[in, out] context    pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_monolithic_init_setup(const cs_navsto_param_t    *nsp,
                                cs_adv_field_t             *adv_field,
                                void                       *context)
{
  cs_navsto_monolithic_t  *nsc = (cs_navsto_monolithic_t *)context;

  assert(nsp != nullptr && nsc != nullptr);

  /* Handle the momentum equation */

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);

  /* Navier-Stokes parameters induce numerical settings for the related
     equations */

  cs_navsto_param_transfer(nsp, mom_eqp);

  /* Link the time property to the momentum equation */

  if (!cs_navsto_param_is_steady(nsp))
    cs_equation_add_time(mom_eqp, nsp->mass_density);

  /* Add advection term in case of CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES
   * CS_NAVSTO_MODEL_OSEEN: Nothing to do since the Oseen field is set by the
   * user via cs_navsto_add_oseen_field() */

  if (nsp->model & CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES)
    cs_equation_add_advection(mom_eqp, adv_field);

  /* The matrix is symmetric for Stokes */
  if (nsp->model == CS_NAVSTO_MODEL_STOKES) {
    mom_eqp->sles_param->mat_is_sym = true;
  }

  /* All considered models needs a viscous term */

  cs_equation_add_diffusion(mom_eqp, nsp->tot_viscosity);

  /* Add the variable field (Always keep a previous state) */

  cs_equation_predefined_create_field(1, nsc->momentum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when a monolithic
 *         algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_monolithic_last_setup(const cs_navsto_param_t     *nsp,
                                void                        *context)
{
  cs_navsto_monolithic_t  *nsc = (cs_navsto_monolithic_t *)context;

  assert(nsp != nullptr && nsc != nullptr);

  return; /* Nothing to do up to now */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to the \ref cs_equation_t structure related to
 *         the momentum equation in case of a monolithic coupling
 *
 * \param[in] context  pointer to a context structure cast on-the-fly
 *
 * \return a pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_navsto_monolithic_get_momentum_eq(void       *context)
{
  if (context == nullptr)
    return nullptr;

  cs_navsto_monolithic_t  *nsc = (cs_navsto_monolithic_t *)context;

  return nsc->momentum;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using an incremental Projection approach in the
 *         the rotational form (see Minev & Guermond, 2006, JCP)
 *
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_projection_create_context(cs_param_bc_type_t    bc,
                                    cs_navsto_param_t    *nsp)
{
  cs_navsto_projection_t *nsc = nullptr;

  BFT_MALLOC(nsc, 1, cs_navsto_projection_t);

  nsc->prediction = cs_equation_add("velocity_prediction",
                                    "predicted_velocity",
                                    CS_EQUATION_TYPE_NAVSTO,
                                    3,
                                    bc);

  /* Set the default solver settings */
  {
    cs_equation_param_t  *eqp = cs_equation_get_param(nsc->prediction);

    /* Space scheme settings (default) */

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_COEF, "sushi");

    /* Solver settings */

    if (nsp->model == CS_NAVSTO_MODEL_STOKES)
      cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "cg");
    else
      cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "gcr");
  }

  /* The default boundary condition on the pressure field is always a
     homogeneous Neumann */

  nsc->correction = cs_equation_add("pressure_correction",
                                    "phi",
                                    CS_EQUATION_TYPE_NAVSTO,
                                    1,
                                    CS_BC_SYMMETRY);

  /* Set the default solver settings */
  {
    cs_equation_param_t  *eqp = cs_equation_get_param(nsc->correction);

    /* Default sles parameters are set following the CDO-Cb scheme setting */

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_cb");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_COEF, "gcr");

  }

  nsc->div_st             = nullptr;
  nsc->bdy_pressure_incr  = nullptr;
  nsc->predicted_velocity = nullptr;
  nsc->phi                = nullptr;

  return nsc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to a Projection approach
 *
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a nullptr pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_projection_free_context(void           *context)
{
  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  BFT_FREE(nsc->div_st);
  BFT_FREE(nsc->bdy_pressure_incr);

  BFT_FREE(nsc);

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Start setting-up the Navier-Stokes equations when a projection
 *        algorithm is used to coupled the system.
 *        No mesh information is available at this stage.
 *
 * \param[in]      nsp           pointer to a \ref cs_navsto_param_t structure
 * \param[in]      adv_field     pointer to a cs_adv_field_t structure
 * \param[in]      has_previous  values at different time steps (true/false)
 * \param[in, out] context       pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_projection_init_setup(const cs_navsto_param_t *nsp,
                                cs_adv_field_t           *adv_field,
                                bool                      has_previous,
                                void                     *context)
{
  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  assert(nsp != nullptr && nsc != nullptr);

  /* Prediction step: Approximate the velocity */
  /* ----------------------------------------- */

  cs_equation_param_t *u_eqp = cs_equation_get_param(nsc->prediction);

  cs_navsto_param_transfer(nsp, u_eqp);

  /* There is always a time derivative with projection algorithm */

  cs_equation_add_time(u_eqp, nsp->mass_density);

  /* All considered models needs a viscous term */

  cs_equation_add_diffusion(u_eqp, nsp->tot_viscosity);

  /* Add an advection term in case of
   * CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES
   *
   * If the model is set to CS_NAVSTO_MODEL_OSEEN, then there is nothing to do
   * since the Oseen field is set by the user via cs_navsto_add_oseen_field()
   */

  if (nsp->model & CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES)
    cs_equation_add_advection(u_eqp, adv_field);

  /* Correction step: Approximate the pressure */
  /* ----------------------------------------- */

  cs_equation_param_t *p_eqp = cs_equation_get_param(nsc->correction);

  if (nsp->dof_reduction_mode != p_eqp->dof_reduction)
    p_eqp->dof_reduction = nsp->dof_reduction_mode;

  cs_equation_add_diffusion(p_eqp, cs_property_by_name("unity"));

  /* Add the variable field */

  cs_equation_predefined_create_field((has_previous ? 1 : 0), nsc->prediction);
  cs_equation_predefined_create_field(0, nsc->correction);

  nsc->predicted_velocity = cs_equation_get_field(nsc->prediction);
  nsc->phi = cs_equation_get_field(nsc->correction);

  if (nsp->verbosity > 1)
    nsc->pressure_incr_gradient =
      cs_field_find_or_create("pressure_increment_gradient",
                              CS_FIELD_INTENSIVE,
                              CS_MESH_LOCATION_CELLS,
                              3,
                              false);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when a
 *         projection algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_projection_last_setup(const cs_cdo_quantities_t  *quant,
                                const cs_navsto_param_t    *nsp,
                                void                       *context)
{
  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  assert(nsp != nullptr && nsc != nullptr);

  /* Source term in the correction step stems from the divergence of the
     predicted velocity */

  BFT_MALLOC(nsc->div_st, quant->n_cells, cs_real_t);
  cs_array_real_fill_zero(quant->n_cells, nsc->div_st);

  cs_equation_t  *corr_eq = nsc->correction;
  cs_equation_param_t  *corr_eqp = cs_equation_get_param(corr_eq);

  cs_equation_add_source_term_by_array(corr_eqp,
                                       nullptr, /* all cells */
                                       cs_flag_primal_cell,
                                       nsc->div_st,
                                       false, /* xdef is not owner */
                                       true); /* full length */

  /* Defined BC for the pressure increment in the correction step */

  BFT_MALLOC(nsc->bdy_pressure_incr, quant->n_b_faces, cs_real_t);
  cs_array_real_fill_zero(quant->n_b_faces, nsc->bdy_pressure_incr);

  for (int id = 0; id < nsp->n_pressure_bc_defs; id++) {
    cs_xdef_t *pdef = nsp->pressure_bc_defs[id];
    if (pdef->meta & CS_CDO_BC_DIRICHLET) {
      const cs_zone_t *z = cs_boundary_zone_by_id(pdef->z_id);
      cs_equation_add_bc_by_array(corr_eqp,
                                  CS_BC_DIRICHLET,
                                  z->name,
                                  cs_flag_primal_face,
                                  nsc->bdy_pressure_incr,
                                  false, /* is not owner */
                                  true); /* full list */
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to the \ref cs_equation_t structure related to
 *         the momentum equation in case of a projection coupling
 *
 * \param[in] context  pointer to a context structure cast on-the-fly
 *
 * \return a pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_navsto_projection_get_momentum_eq(void       *context)
{
  if (context == nullptr)
    return nullptr;

  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  return nsc->prediction;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
