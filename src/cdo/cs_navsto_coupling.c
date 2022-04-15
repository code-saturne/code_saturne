/*============================================================================
 * Functions to handle the settings of coupling algorithms
 *  - Artificial Compressibility algorithm
 *  - Its variant VVP (Vector Projection Penalty) algorithm
 *  - Projection algorithm
 *============================================================================*/

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
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

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

#define CS_NAVSTO_COUPLING_DBG  0

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
  cs_navsto_ac_t  *nsc = NULL;

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
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "cg");
  else
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "gcr");

  return nsc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to an Artificial Compressibility
 *         approach
 *
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_free_context(void     *context)
{
  cs_navsto_ac_t  *nsc = (cs_navsto_ac_t *)context;

  BFT_FREE(nsc);

  return NULL;
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

  assert(nsp != NULL && nsc != NULL);

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

  assert(nsp != NULL && nsc != NULL);

  /* Avoid no definition of the zeta coefficient */

  if (nsc->zeta->n_definitions == 0)
    cs_property_def_iso_by_value(nsc->zeta, NULL, nsp->gd_scale_coef);

  /* Set the quadrature level for BCs, if needed */

  const cs_equation_param_t *eqp = cs_equation_get_param(nsc->momentum);

  for (short int i = 0; i < eqp->n_bc_defs; i++) {

    cs_xdef_t *def = eqp->bc_defs[i];
    if (def->type == CS_XDEF_BY_ANALYTIC_FUNCTION) /* Otherwise not useful */
      cs_xdef_set_quadrature(def, nsp->qtype);

  } /* Loop on BC definitions */
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
  if (context == NULL)
    return NULL;

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
  cs_navsto_monolithic_t  *nsc = NULL;

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
  cs_equation_param_set(mom_eqp, CS_EQKEY_HODGE_DIFF_COEF, "sushi");

  /* Solver settings: Only the linear algebra settings related to the momentum
  *  equation.  The strategy is set in _navsto_param_sles_create() */

  if (nsp->model ==  CS_NAVSTO_MODEL_STOKES) {
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "cg");
  }
  else {
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "poly1");
    cs_equation_param_set(mom_eqp, CS_EQKEY_ITSOL, "gcr");
  }

  return nsc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to a monolithic approach
 *
 * \param[in, out] context     pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_monolithic_free_context(void             *context)
{
  cs_navsto_monolithic_t  *nsc = (cs_navsto_monolithic_t *)context;

  BFT_FREE(nsc);

  return NULL;
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

  assert(nsp != NULL && nsc != NULL);

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

  assert(nsp != NULL && nsc != NULL);

  /* Set the quadrature level for BCs, if needed */

  const cs_equation_param_t *eqp = cs_equation_get_param(nsc->momentum);

  for (short int i = 0; i < eqp->n_bc_defs; i++) {

    cs_xdef_t *def = eqp->bc_defs[i];
    if (def->type == CS_XDEF_BY_ANALYTIC_FUNCTION) /* Otherwise not useful */
      cs_xdef_set_quadrature(def, nsp->qtype);

  } /* Loop on BC definitions */
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
  if (context == NULL)
    return NULL;

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
  cs_navsto_projection_t  *nsc = NULL;

  BFT_MALLOC(nsc, 1, cs_navsto_projection_t);

  nsc->prediction = cs_equation_add("velocity_prediction",
                                    "velocity",
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
      cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "cg");
    else
      cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "gcr");
  }

  /* The default boundary condition on the pressure field is always a
     homogeneous Neumann */

  nsc->correction = cs_equation_add("pressure_correction",
                                    "phi",
                                    CS_EQUATION_TYPE_NAVSTO,
                                    1,
                                    CS_PARAM_BC_HMG_NEUMANN);

  /* Set the default solver settings */
  {
    cs_equation_param_t  *eqp = cs_equation_get_param(nsc->correction);

    /* Space scheme settings (default) */

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_COEF, "sushi");

    /* Solver settings */

    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "cg");
  }

  nsc->div_st = NULL;
  nsc->bdy_pressure_incr = NULL;
  nsc->predicted_velocity = NULL;

  return nsc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to a Projection approach
 *
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_projection_free_context(void           *context)
{
  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  BFT_FREE(nsc->div_st);
  BFT_FREE(nsc->bdy_pressure_incr);

  BFT_FREE(nsc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a projection
 *         algorithm is used to coupled the system.
 *         No mesh information is available at this stage.
 *
 * \param[in]      nsp           pointer to a \ref cs_navsto_param_t structure
 * \param[in]      adv_field     pointer to a cs_adv_field_t structure
 * \param[in]      loc_id        id related to a mesh location
 * \param[in]      has_previous  values at different time steps (true/false)
 * \param[in, out] context       pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_projection_init_setup(const cs_navsto_param_t    *nsp,
                                cs_adv_field_t             *adv_field,
                                int                         loc_id,
                                bool                        has_previous,
                                void                       *context)
{
  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  assert(nsp != NULL && nsc != NULL);

  /* Prediction step: Approximate the velocity */
  /* ----------------------------------------- */

  cs_equation_param_t *u_eqp = cs_equation_get_param(nsc->prediction);

  cs_navsto_param_transfer(nsp, u_eqp);

  /* There is always a time derivative with projection algorithm */

  cs_equation_add_time(u_eqp, nsp->mass_density);

  /* All considered models needs a viscous term */

  cs_equation_add_diffusion(u_eqp, nsp->tot_viscosity);

  /* Add advection term in case of CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES
   * CS_NAVSTO_MODEL_OSEEN: Nothing to do since the Oseen field is set by the
   * user via cs_navsto_add_oseen_field() */

  if (nsp->model & CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES)
    cs_equation_add_advection(u_eqp, adv_field);

  /* Add the variable field (Always keep a previous state) */

  cs_equation_predefined_create_field(1, nsc->prediction);

  /* Correction step: Approximate the pressure */
  /* ----------------------------------------- */

  cs_equation_param_t *p_eqp = cs_equation_get_param(nsc->correction);

  cs_navsto_param_transfer(nsp, p_eqp);

  cs_equation_add_diffusion(p_eqp, cs_property_by_name("time_step"));

  /* Add the predicted velocity field */

  nsc->predicted_velocity = cs_field_create("predicted_velocity",
                                            CS_FIELD_INTENSIVE,
                                            loc_id,
                                            3,
                                            has_previous);

  /* Add the variable field */

  cs_equation_predefined_create_field((has_previous ? 1 : 0), nsc->prediction);
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

  assert(nsp != NULL && nsc != NULL);

  /* Source term in the correction step stems from the divergence of the
     predicted velocity */

  BFT_MALLOC(nsc->div_st, quant->n_cells, cs_real_t);
  memset(nsc->div_st, 0, quant->n_cells*sizeof(cs_real_t));

  cs_equation_t  *corr_eq = nsc->correction;
  cs_equation_param_t  *corr_eqp = cs_equation_get_param(corr_eq);

  cs_equation_add_source_term_by_array(corr_eqp,
                                       NULL,
                                       cs_flag_primal_cell,
                                       nsc->div_st,
                                       false,       /* xdef is not owner */
                                       NULL, NULL); /* no index/ids */

  /* Defined BC for the pressure increment in the correction step */

  BFT_MALLOC(nsc->bdy_pressure_incr, quant->n_b_faces, cs_real_t);
  memset(nsc->bdy_pressure_incr, 0, quant->n_b_faces*sizeof(cs_real_t));

  for (int i = 0; i < nsp->n_pressure_bc_defs; i++) {

    const cs_xdef_t  *pdef = nsp->pressure_bc_defs[i];
    const cs_zone_t  *z = cs_boundary_zone_by_id(pdef->z_id);

    cs_equation_add_bc_by_array(corr_eqp,
                                CS_PARAM_BC_DIRICHLET,
                                z->name,
                                cs_flag_primal_face,
                                nsc->bdy_pressure_incr,
                                false, /* xdef is not owner */
                                NULL, NULL); /* no index, no ids */

  } /* Loop on pressure definitions */
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
  if (context == NULL)
    return NULL;

  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  return nsc->prediction;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
