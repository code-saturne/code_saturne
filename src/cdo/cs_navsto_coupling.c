/*============================================================================
 * Routines to handle the settings of coupling algorithms
 *  - Uzawa-Augmented Lagrangian algorithm
 *  - Artificial Compressibility algorithm
 *  - Its variant VVP (Vector Projection Penalty) algorithm
 *  - Projection algorithm
 *============================================================================*/

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
 *  \file cs_navsto_coupling.c
 *
 * \brief Routines to handle the settings of coupling algorithms
 *  - Uzawa-Augmented Lagrangian algorithm
 *  - Artificial Compressibility algorithm
 *  - Its variant VVP (Vector Projection Penalty) algorithm
 *  - Projection algorithm
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
 *         system is coupled using an Uzawa-Augmented Lagrangian approach
 *
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_uzawa_create_context(cs_navsto_param_t    *nsp,
                               cs_param_bc_type_t    bc)
{
  assert(nsp != NULL);
  CS_UNUSED(nsp); /* Avoid a warning when compiling */

  cs_navsto_uzawa_t  *nsc = NULL;

  BFT_MALLOC(nsc, 1, cs_navsto_uzawa_t);

  /* Add an equation for the momentum conservation */
  nsc->momentum = cs_equation_add("Momentum",
                                  "velocity",
                                  CS_EQUATION_TYPE_PREDEFINED,
                                  3,
                                  bc);

  /* Set the default settings for the momentum equation */
  {
    cs_equation_param_t  *eqp = cs_equation_get_param(nsc->momentum);

    /* Solver settings */
    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "jacobi");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "bicg");
  }

  nsc->energy = NULL;   /* Not used up to now */

  nsc->zeta = cs_property_add("graddiv_coef", CS_PROPERTY_ISO);
  nsc->relax  = 1.0;

  return nsc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to an Uzawa-Augmented Lagrangian
 *         approach
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_uzawa_free_context(const cs_navsto_param_t    *nsp,
                             void                       *context)
{
  assert(nsp != NULL);
  CS_UNUSED(nsp); /* Avoid a warning when compiling */

  cs_navsto_uzawa_t  *nsc = (cs_navsto_uzawa_t *)context;

  BFT_FREE(nsc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a Uzawa
 *         Augmented Lagrangian algorithm is used to coupled the system
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_uzawa_init_setup(const cs_navsto_param_t    *nsp,
                           void                       *context)
{
  cs_navsto_uzawa_t  *nsc = (cs_navsto_uzawa_t *)context;

  assert(nsp != NULL && nsc != NULL);

  /* Handle the momentum equation */
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);

  cs_navsto_param_transfer(nsp, mom_eqp);

  /* Link the time property to the momentum equation */
  switch (nsp->time_state) {

  case CS_NAVSTO_TIME_STATE_UNSTEADY:
  case CS_NAVSTO_TIME_STATE_LIMIT_STEADY:
    cs_equation_add_time(mom_eqp, cs_property_by_name("unity"));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid choice for the time state", __func__);
  }

  /* All considered models needs a viscous term */
  cs_equation_add_diffusion(mom_eqp, nsp->lami_viscosity);

  /* Handle the energy equation */
  if (nsc->energy != NULL)
    cs_navsto_param_transfer(nsp, cs_equation_get_param(nsc->energy));

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when an Uzawa
 *         Augmented Lagrangian algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      connect  pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_uzawa_last_setup(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_navsto_param_t     *nsp,
                           void                        *context)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_navsto_uzawa_t  *nsc = (cs_navsto_uzawa_t *)context;

  assert(nsp != NULL && nsc != NULL);

  /* Avoid no definition of the zeta coefficient */
  if (nsc->zeta->n_definitions == 0)
    cs_property_def_iso_by_value(nsc->zeta, NULL, nsp->gd_scale_coef);

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using an Artificial Compressibility approach
 *
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_create_context(cs_navsto_param_t    *nsp,
                            cs_param_bc_type_t    bc)
{
  assert(nsp != NULL);
  CS_UNUSED(nsp); /* Avoid warning when compiling */

  cs_navsto_ac_t  *nsc = NULL;

  BFT_MALLOC(nsc, 1, cs_navsto_ac_t);

  nsc->momentum = cs_equation_add("Momentum",
                                  "velocity",
                                  CS_EQUATION_TYPE_PREDEFINED,
                                  3,
                                  bc);

  /* Set the default solver settings */
  {
    cs_equation_param_t  *eqp = cs_equation_get_param(nsc->momentum);

    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "jacobi");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "bicg");
  }

  /* Additional property */
  nsc->zeta = cs_property_add("graddiv_coef", CS_PROPERTY_ISO);

  return nsc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to an Artificial Compressibility
 *         approach
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_free_context(const cs_navsto_param_t    *nsp,
                          void                       *context)
{
  assert(nsp != NULL);
  CS_UNUSED(nsp); /* Avoid warning when compiling with optimizations */

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
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_ac_init_setup(const cs_navsto_param_t    *nsp,
                        void                       *context)
{
  cs_navsto_ac_t  *nsc = (cs_navsto_ac_t *)context;

  assert(nsp != NULL && nsc != NULL);

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);

  /* Navier-Stokes parameters induce numerical settings for the related
     equations */
  cs_navsto_param_transfer(nsp, mom_eqp);

  /* Link the time property to the momentum equation */
  switch (nsp->time_state) {

  case CS_NAVSTO_TIME_STATE_UNSTEADY:
  case CS_NAVSTO_TIME_STATE_LIMIT_STEADY:
    cs_equation_add_time(mom_eqp, cs_property_by_name("unity"));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid choice for the time state", __func__);
  }

  /* All considered models needs a viscous term */
  cs_equation_add_diffusion(mom_eqp, nsp->lami_viscosity);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when an
 *         Artificial Compressibility algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      connect  pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_ac_last_setup(const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_navsto_param_t     *nsp,
                        void                        *context)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_navsto_ac_t  *nsc = (cs_navsto_ac_t *)context;

  assert(nsp != NULL && nsc != NULL);

  /* Avoid no definition of the zeta coefficient */
  if (nsc->zeta->n_definitions == 0)
    cs_property_def_iso_by_value(nsc->zeta, NULL, nsp->gd_scale_coef);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using an Artificial Compressibility - VPP approach
 *
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_vpp_create_context(cs_navsto_param_t    *nsp,
                                 cs_param_bc_type_t    bc)
{
  assert(nsp != NULL);
  CS_UNUSED(nsp); /* Avoid warning when compiling */

  cs_navsto_ac_vpp_t  *nsc = NULL;

  BFT_MALLOC(nsc, 1, cs_navsto_ac_vpp_t);

  nsc->momentum = cs_equation_add("Momentum",
                                  "Utilda",
                                  CS_EQUATION_TYPE_PREDEFINED,
                                  3,
                                  bc);

  /* Set the default solver settings for "Momentum" */
  {
    cs_equation_param_t  *eqp = cs_equation_get_param(nsc->momentum);

    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "jacobi");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
  }

  /* The grad-div equation is usually always with homogeneous Dirichlet */
  nsc->graddiv = cs_equation_add("Graddiv",
                                 "Uhat",
                                 CS_EQUATION_TYPE_PREDEFINED,
                                 3,
                                 CS_PARAM_BC_HMG_DIRICHLET);

  /* Set the default solver settings for "Graddiv" */
  {
    cs_equation_param_t  *eqp = cs_equation_get_param(nsc->graddiv);

    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "jacobi");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
  }

  nsc->zeta = cs_property_add("graddiv_coef", CS_PROPERTY_ISO);

  return nsc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to an Artificial Compressibility
 *         with the VPP approach
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_ac_vpp_free_context(const cs_navsto_param_t    *nsp,
                              void                       *context)
{
  assert(nsp != NULL);
  CS_UNUSED(nsp); /* Avoid warning when compiling with optimizations */

  cs_navsto_ac_vpp_t  *nsc = (cs_navsto_ac_vpp_t *)context;

  BFT_FREE(nsc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when an Artificial
 *         Compressibility with VPP algorithm is used to coupled the system.
 *         No mesh information is available at this stage.
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_ac_vpp_init_setup(const cs_navsto_param_t    *nsp,
                            void                       *context)
{
  cs_navsto_ac_vpp_t *nsc = (cs_navsto_ac_vpp_t *)context;

  assert(nsp != NULL && nsc != NULL);

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);
  cs_equation_param_t  *gd_eqp = cs_equation_get_param(nsc->graddiv);

  /* Navier-Stokes parameters induce numerical settings for the related
     equations */
  cs_navsto_param_transfer(nsp, mom_eqp);

  /* TODO: Is this good? Should we force BC or alikes? */
  cs_navsto_param_transfer(nsp, gd_eqp);

  /* Link the time property to the momentum equation */
  switch (nsp->time_state) {

  case CS_NAVSTO_TIME_STATE_UNSTEADY:
  case CS_NAVSTO_TIME_STATE_LIMIT_STEADY:
    cs_equation_add_time(mom_eqp, cs_property_by_name("unity"));
    cs_equation_add_time(gd_eqp, cs_property_by_name("unity"));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid choice for the time state", __func__);
  }

  /* All considered models needs a viscous term */
  cs_equation_add_diffusion(mom_eqp, nsp->lami_viscosity);
  cs_equation_add_diffusion(gd_eqp, nsp->lami_viscosity);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when an
 *         Artificial Compressibility algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      connect  pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_ac_vpp_last_setup(const cs_cdo_connect_t      *connect,
                            const cs_cdo_quantities_t   *quant,
                            const cs_navsto_param_t     *nsp,
                            void                        *context)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_navsto_ac_vpp_t  *nsc = (cs_navsto_ac_vpp_t *)context;

  assert(nsp != NULL && nsc != NULL);

  /* Avoid no definition of the zeta coefficient */
  if (nsc->zeta->n_definitions == 0)
    cs_property_def_iso_by_value(nsc->zeta, NULL, nsp->gd_scale_coef);

  /* TODO: Setting quadrature for the source terms */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using an incremental Projection approach in the
 *         the rotational form (see Minev & Guermond, 2006, JCP)
 *
 * \param[in]  nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in]  bc     default \ref cs_param_bc_type_t for the equation
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_projection_create_context(cs_navsto_param_t    *nsp,
                                    cs_param_bc_type_t    bc)
{
  assert(nsp != NULL);
  CS_UNUSED(nsp); /* Avoid warning when compiling */

  cs_navsto_projection_t  *nsc = NULL;

  BFT_MALLOC(nsc, 1, cs_navsto_projection_t);

  nsc->prediction = cs_equation_add("Velocity_Prediction",
                                    "velocity",
                                    CS_EQUATION_TYPE_PREDEFINED,
                                    3,
                                    bc);

  /* Set the default solver settings */
  {
    cs_equation_param_t  *eqp = cs_equation_get_param(nsc->prediction);

    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "jacobi");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "bicg");
  }

  nsc->correction = cs_equation_add("Pressure_Correction",
                                    "phi",
                                    CS_EQUATION_TYPE_PREDEFINED,
                                    1,
                                    CS_PARAM_BC_HMG_NEUMANN);

  /* Set the default solver settings */
  {
    cs_equation_param_t  *eqp = cs_equation_get_param(nsc->correction);

    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
  }

  return nsc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the context structure related to a Projection approach
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_navsto_projection_free_context(const cs_navsto_param_t    *nsp,
                                  void                       *context)
{
  assert(nsp != NULL);
  CS_UNUSED(nsp); /* Avoid warning when compiling */

  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  BFT_FREE(nsc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a projection
 *         algorithm is used to coupled the system.
 *         No mesh information is available at this stage.
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_projection_init_setup(const cs_navsto_param_t    *nsp,
                                void                       *context)
{
  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  assert(nsp != NULL && nsc != NULL);

  /* Prediction step: Approximate the velocity */
  cs_equation_param_t *p_eqp = cs_equation_get_param(nsc->prediction);

  cs_navsto_param_transfer(nsp, p_eqp);

  cs_equation_add_time(p_eqp, cs_property_by_name("unity"));

  /* Correction step: Approximate the pressure */
  cs_navsto_param_transfer(nsp, cs_equation_get_param(nsc->correction));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when a
 *         projection algorithm is used to coupled the system.
 *         Connectivity and geometric quantities are available at this stage.
 *
 * \param[in]      connect  pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_projection_last_setup(const cs_cdo_connect_t     *connect,
                                const cs_cdo_quantities_t  *quant,
                                const cs_navsto_param_t    *nsp,
                                void                       *context)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_navsto_projection_t  *nsc = (cs_navsto_projection_t *)context;

  assert(nsp != NULL && nsc != NULL);

  /* TODO */
  CS_UNUSED(nsp);
  CS_UNUSED(nsc);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
