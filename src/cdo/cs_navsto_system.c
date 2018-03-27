/*============================================================================
 * Routines to handle cs_navsto_system_t structure
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

#include "cs_cdofb_navsto.h"
#include "cs_hho_stokes.h"
#include "cs_equation.h"
#include "cs_log.h"
#include "cs_post.h"
#include "cs_flag.h"
#include "cs_volume_zone.h"
#include "cs_evaluate.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_navsto_coupling.h"
#include "cs_navsto_system.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 *  \file cs_navsto_system.c
 *
 *  \brief  Routines to handle the cs_navsto_system_t structure
 */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_NAVSTO_SYSTEM_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_ns[] =
  " Stop execution. The structure related to the Navier-Stokes system is"
  " empty.\n Please check your settings.\n";

static const char _err_invalid_coupling[] =
  " %s: Invalid case for the coupling algorithm.\n";

static const char
_space_scheme_key[CS_SPACE_N_SCHEMES][CS_BASE_STRING_LEN] =
  { "fv",
    "cdo_vb",
    "cdo_vcb",
    "cdo_fb",
    "hho_p0",
    "hho_p1",
    "hho_p2"
  };

static const char
_time_scheme_key[CS_TIME_N_SCHEMES][CS_BASE_STRING_LEN] =
  { "implicit",
    "explicit",
    "crank_nicolson",
    "theta_scheme"
  };

static const char
_dof_reduction_key[CS_PARAM_N_REDUCTIONS][CS_BASE_STRING_LEN] =
  { "derham",
    "average"
  };

static const char
_quad_type_key[CS_QUADRATURE_N_TYPES][CS_BASE_STRING_LEN] =
  { "none",
    "bary",
    "bary_subdiv",
    "higher",
    "highest"
  };

static cs_navsto_system_t  *cs_navsto_system = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate an empty Navier-Stokes system
 *
 * \return a pointer to a new allocated groundwater flow structure
 */
/*----------------------------------------------------------------------------*/

static cs_navsto_system_t *
_allocate_navsto_system(void)
{
  cs_navsto_system_t  *navsto = NULL;

  BFT_MALLOC(navsto, 1, cs_navsto_system_t);

  navsto->param = NULL;

  navsto->adv_field = NULL;

  /* Main set of variables */
  navsto->velocity = NULL;
  navsto->pressure = NULL;
  navsto->temperature = NULL;

  /* Main set of properties */
  navsto->density = NULL;
  navsto->lami_viscosity = NULL;

  /* Additional data fitting the choice of model */
  navsto->context = NULL;

  /* Function pointers */
  navsto->init = NULL;
  navsto->free = NULL;
  navsto->compute = NULL;

  /* Initial conditions */
  navsto->n_velocity_ic_defs = 0;
  navsto->n_pressure_ic_defs = 0;
  navsto->velocity_ic_defs = NULL;
  navsto->pressure_ic_defs = NULL;

  return navsto;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the numerical settings defined for the Navier-Stokes system
 *         to an equation related to this system.
 *
 * \param[in]       nsp    pointer to a cs_navsto_param_t structure
 * \param[in, out]  eqp    pointer to a cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_apply_param(const cs_navsto_param_t    *nsp,
             cs_equation_param_t        *eqp)
{
  assert(nsp != NULL && eqp != NULL);

  /*  Set the space discretization scheme */
  const char  *ss_key = _space_scheme_key[nsp->space_scheme];

  cs_equation_set_param(eqp, CS_EQKEY_SPACE_SCHEME, ss_key);

  /*  Set the time discretization scheme */
  const char  *ts_key = _time_scheme_key[nsp->time_scheme];

  cs_equation_set_param(eqp, CS_EQKEY_TIME_SCHEME, ts_key);
  if (nsp->time_scheme == CS_TIME_SCHEME_THETA) {
    char  cvalue[36]; /* include '\0' */
    snprintf(cvalue, 35*sizeof(char), "%g", nsp->theta);
    cs_equation_set_param(eqp, CS_EQKEY_TIME_THETA, cvalue);
  }

  /*  Set the way DoFs are defined */
  const char  *dof_key = _dof_reduction_key[nsp->dof_reduction_mode];

  cs_equation_set_param(eqp, CS_EQKEY_DOF_REDUCTION, dof_key);

  /*  Set quadratures type */
  const char  *quad_key = _quad_type_key[nsp->qtype];

  cs_equation_set_param(eqp, CS_EQKEY_BC_QUADRATURE, quad_key);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a context structure when the Navier-Stokes
 *         system is coupled using an Uzawa-Augmented Lagrangian approach
 *
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

static void *
_create_uzawa_context(cs_navsto_param_t    *nsp)
{
  assert(nsp != NULL);

  cs_navsto_coupling_uzawa_t  *nsc = NULL;

  BFT_MALLOC(nsc, 1, cs_navsto_coupling_uzawa_t);

  /* Add an equation for the momentum conservation */
  nsc->momentum = cs_equation_add("Momentum",
                                  "velocity",
                                  CS_EQUATION_TYPE_PREDEFINED,
                                  3,
                                  CS_PARAM_BC_HMG_DIRICHLET);

  /* Set the default settings */
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
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static void *
_free_uzawa_context(const cs_navsto_param_t    *nsp,
                    void                       *context)
{
  assert(nsp != NULL);

  cs_navsto_coupling_uzawa_t  *nsc = (cs_navsto_coupling_uzawa_t *)context;

  BFT_FREE(nsc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when an Uzawa
 *         Augmented Lagrangian algorithm is used to coupled the system
 *         No mesh information is available
 *
 * \param[in, out] ns       pointer to a cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_uzawa_init_setup(cs_navsto_system_t          *ns)
{
  assert(ns != NULL);

  cs_navsto_param_t  *nsp = ns->param;
  cs_navsto_coupling_uzawa_t  *nsc = (cs_navsto_coupling_uzawa_t *)ns->context;

  assert(nsp != NULL && nsc != NULL);

  /* Handle the momentum equation */
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);

  _apply_param(nsp, mom_eqp);

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
  cs_equation_add_diffusion(mom_eqp, ns->lami_viscosity);

  /* Handle the energy equation */
  if (nsc->energy != NULL)
    _apply_param(nsp, cs_equation_get_param(nsc->energy));

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when an Uzawa
 *         Augmented Lagrangian algorithm is used to coupled the system
 *
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a cs_cdo_quantities_t structure
 * \param[in, out] ns       pointer to a cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_uzawa_last_setup(const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *quant,
                  cs_navsto_system_t          *ns)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  assert(ns != NULL);

  cs_navsto_param_t  *nsp = ns->param;
  cs_navsto_coupling_uzawa_t  *nsc = (cs_navsto_coupling_uzawa_t *)ns->context;

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
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

static void *
_create_ac_context(cs_navsto_param_t    *nsp)
{
  assert(nsp != NULL);

  cs_navsto_coupling_ac_t  *nsc = NULL;

  BFT_MALLOC(nsc, 1, cs_navsto_coupling_ac_t);

  nsc->momentum = cs_equation_add("Momentum",
                                  "velocity",
                                  CS_EQUATION_TYPE_PREDEFINED,
                                  3,
                                  CS_PARAM_BC_HMG_DIRICHLET);

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
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static void *
_free_ac_context(const cs_navsto_param_t    *nsp,
                 void                       *context)
{
  assert(nsp != NULL);

  cs_navsto_coupling_ac_t  *nsc = (cs_navsto_coupling_ac_t *)context;

  BFT_FREE(nsc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when an
 *         Artificial Compressibility algorithm is used to coupled the system
 *         No mesh information is available
 *
 * \param[in, out] ns       pointer to a cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_ac_init_setup(cs_navsto_system_t     *ns)
{
  assert(ns != NULL);

  cs_navsto_param_t  *nsp = ns->param;
  cs_navsto_coupling_ac_t  *nsc = (cs_navsto_coupling_ac_t *)ns->context;

  assert(nsp != NULL && nsc != NULL);

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);

  /* Navier-Stokes parameters induce numerical settings for the related
   equations */
  _apply_param(nsp, mom_eqp);

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
  cs_equation_add_diffusion(mom_eqp, ns->lami_viscosity);


}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when an
 *         Artificial Compressibility algorithm is used to coupled the system
 *
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a cs_cdo_quantities_t structure
 * \param[in, out] ns       pointer to a cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_ac_last_setup(const cs_cdo_connect_t      *connect,
               const cs_cdo_quantities_t   *quant,
               cs_navsto_system_t          *ns)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  assert(ns != NULL);

  cs_navsto_param_t  *nsp = ns->param;
  cs_navsto_coupling_ac_t  *nsc = (cs_navsto_coupling_ac_t *)ns->context;

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
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

static void *
_create_ac_vpp_context(cs_navsto_param_t    *nsp)
{
  assert(nsp != NULL);

  cs_navsto_coupling_ac_vpp_t  *nsc = NULL;

  BFT_MALLOC(nsc, 1, cs_navsto_coupling_ac_vpp_t);

  nsc->momentum = cs_equation_add("Momentum",
                                  "Utilde",
                                  CS_EQUATION_TYPE_PREDEFINED,
                                  3,
                                  CS_PARAM_BC_HMG_DIRICHLET);
  nsc->graddiv = cs_equation_add("Graddiv",
                                  "Uhat",
                                  CS_EQUATION_TYPE_PREDEFINED,
                                  3,
                                  CS_PARAM_BC_HMG_DIRICHLET);

  /* Set the default solver settings for "Momentum" */
  {
    cs_equation_param_t  *eqp = cs_equation_get_param(nsc->momentum);

    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "jacobi");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
  }

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
 *          - VPP approach
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static void *
_free_ac_vpp_context(const cs_navsto_param_t    *nsp,
                     void                       *context)
{
  assert(nsp != NULL);

  cs_navsto_coupling_ac_vpp_t  *nsc = (cs_navsto_coupling_ac_vpp_t *)context;

  BFT_FREE(nsc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when an
 *         Artificial Compressibility with VPP algorithm is used to coupled the
 *         system
 *         No mesh information is available
 *
 * \param[in, out] ns       pointer to a cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_ac_vpp_init_setup(cs_navsto_system_t     *ns)
{
  assert(ns != NULL);

  cs_navsto_param_t  *nsp = ns->param;
  cs_navsto_coupling_ac_vpp_t *nsc = (cs_navsto_coupling_ac_vpp_t *)ns->context;

  assert(nsp != NULL && nsc != NULL);

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);
  cs_equation_param_t  *gd_eqp = cs_equation_get_param(nsc->graddiv);

  /* Navier-Stokes parameters induce numerical settings for the related
     equations */
  _apply_param(nsp, mom_eqp);

  // TODO: Is this good? Should we force BC or alikes?
  _apply_param(nsp, gd_eqp);

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
  cs_equation_add_diffusion(mom_eqp, ns->lami_viscosity);
  cs_equation_add_diffusion(gd_eqp, ns->lami_viscosity);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when an
 *         Artificial Compressibility algorithm is used to coupled the system
 *
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a cs_cdo_quantities_t structure
 * \param[in, out] ns       pointer to a cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_ac_vpp_last_setup(const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   cs_navsto_system_t          *ns)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  assert(ns != NULL);

  cs_navsto_param_t  *nsp = ns->param;
  cs_navsto_coupling_ac_vpp_t *nsc = (cs_navsto_coupling_ac_vpp_t *)ns->context;

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
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 *
 * \return a pointer to the context structure
 */
/*----------------------------------------------------------------------------*/

static void *
_create_projection_context(cs_navsto_param_t    *nsp)
{
  assert(nsp != NULL);

  cs_navsto_coupling_projection_t  *nsc = NULL;

  BFT_MALLOC(nsc, 1, cs_navsto_coupling_projection_t);

  nsc->prediction = cs_equation_add("Velocity_Prediction",
                                    "velocity",
                                    CS_EQUATION_TYPE_PREDEFINED,
                                    3,
                                    CS_PARAM_BC_HMG_DIRICHLET);

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
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static void *
_free_projection_context(const cs_navsto_param_t    *nsp,
                         void                       *context)
{
  assert(nsp != NULL);

  cs_navsto_coupling_projection_t  *nsc =
    (cs_navsto_coupling_projection_t *)context;

  BFT_FREE(nsc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a
 *         projection algorithm is used to coupled the system
 *         No mesh information is available
 *
 * \param[in, out] ns       pointer to a cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_projection_init_setup(cs_navsto_system_t         *ns)
{
  assert(ns != NULL);

  cs_navsto_param_t  *nsp = ns->param;
  cs_navsto_coupling_projection_t  *nsc =
    (cs_navsto_coupling_projection_t *)ns->context;

  assert(nsp != NULL && nsc != NULL);

  /* Prediction step: Approximate the velocity */
  cs_equation_param_t *p_eqp = cs_equation_get_param(nsc->prediction);

  _apply_param(nsp, p_eqp);

  cs_equation_add_time(p_eqp, cs_property_by_name("unity"));

  /* Correction step: Approximate the pressure */
  _apply_param(nsp, cs_equation_get_param(nsc->correction));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup for the Navier-Stokes equations when a
 *         projection algorithm is used to coupled the system
 *
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a cs_cdo_quantities_t structure
 * \param[in, out] ns       pointer to a cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_projection_last_setup(const cs_cdo_connect_t     *connect,
                       const cs_cdo_quantities_t  *quant,
                       cs_navsto_system_t         *ns)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  assert(ns != NULL);

  cs_navsto_param_t  *nsp = ns->param;
  cs_navsto_coupling_projection_t  *nsc =
    (cs_navsto_coupling_projection_t *)ns->context;

  assert(nsp != NULL && nsc != NULL);

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the face values for the velocity unknows (in case of
 *        CDO Face-based scheme) in order to have a initial guess in accordance
 *        with the user requirements when solving the momentum equation
 */
/*----------------------------------------------------------------------------*/

static inline void
_init_face_velocity_values(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  const cs_navsto_param_t *nsp = navsto->param;
  const cs_flag_t f_dof_flag = CS_FLAG_VECTOR | cs_flag_primal_face;
  const cs_param_dof_reduction_t red = nsp->dof_reduction_mode;

  cs_equation_t *first_eq = NULL;

  /* Switching on coupling in order to set the values on the faces which will
   * be used as initial guess in the first solve */
  switch (nsp->coupling) {
  case CS_NAVSTO_COUPLING_UZAWA:
    first_eq = ((cs_navsto_coupling_uzawa_t*)navsto->context)->momentum;
    break;

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    first_eq = ((cs_navsto_coupling_ac_t*)navsto->context)->momentum;
    break;

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
    first_eq = ((cs_navsto_coupling_ac_vpp_t*)navsto->context)->momentum;
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
    first_eq =
      ((cs_navsto_coupling_projection_t*)navsto->context)->prediction;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;

  } /* Switch on coupling algorithm */

  cs_real_t  *f_values = NULL;
  f_values = cs_equation_get_face_values(first_eq);

  for (int def_id = 0; def_id < navsto->n_velocity_ic_defs; def_id++) {

    /* Get and then set the definition of the initial condition */
    cs_xdef_t  *def = navsto->velocity_ic_defs[def_id];

    /* Forcing definition on all faces
       TODO: enable multiple definitions */
    const cs_flag_t meta_cpy = def->meta;
    def->meta = CS_FLAG_FULL_LOC;

    /* Initialize face-based array */
    switch(def->type) {

    case CS_XDEF_BY_VALUE:
      cs_evaluate_potential_by_value(f_dof_flag, def, f_values);
      break;

    case CS_XDEF_BY_ANALYTIC_FUNCTION:

      /* An evaluation at the barycenter will be enough */
      cs_xdef_set_quadrature(def, CS_QUADRATURE_BARY);

      switch (red) {
      case CS_PARAM_REDUCTION_DERHAM:
        cs_evaluate_potential_by_analytic(f_dof_flag, def, f_values);
        break;
      case CS_PARAM_REDUCTION_AVERAGE:
        cs_evaluate_average_on_faces_by_analytic(f_dof_flag, def, f_values);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: Incompatible reduction.\n"), __func__);

      } /* Switch on possible reduction types */

      /* Switch to the initial setting */
      cs_xdef_set_quadrature(def, nsp->qtype);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Incompatible initialization.\n"), __func__);

    } // Switch on possible type of definition

    /* Resetting */
    def->meta = meta_cpy;

  } // Loop on definitions
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the resolution of the Navier-Stokes system has been
 *        activated
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_navsto_system_is_activated(void)
{
  if (cs_navsto_system == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the Navier-Stokes (NS) system
 *
 * \param[in] model          type of model related to the NS system
 * \param[in] time_state     state of the time for the NS equations
 * \param[in] algo_coupling  algorithm used for solving the NS system
 *
 * \return a pointer to a new allocated cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_system_t *
cs_navsto_system_activate(cs_navsto_param_model_t        model,
                          cs_navsto_param_time_state_t   time_state,
                          cs_navsto_param_coupling_t     algo_coupling)
{
  /* Sanity checks */
  if (model == CS_NAVSTO_N_MODELS)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model for Navier-Stokes.\n",
              __func__);

  /* Allocate an empty structure */
  cs_navsto_system_t  *navsto = _allocate_navsto_system();

  /* Initialize the set of parameters */
  navsto->param = cs_navsto_param_create(model, time_state, algo_coupling);

  /* Main set of properties */
  navsto->density = cs_property_add("density", CS_PROPERTY_ISO);
  navsto->lami_viscosity = cs_property_add("laminar_viscosity",
                                           CS_PROPERTY_ISO);

  /* Advection field related to the resolved velocity */
  navsto->adv_field = cs_advection_field_add("velocity_field");

  /* Additional initialization fitting the choice of model */
  switch (navsto->param->coupling) {

  case CS_NAVSTO_COUPLING_UZAWA:
    navsto->context = _create_uzawa_context(navsto->param);
    break;

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    navsto->context = _create_ac_context(navsto->param);
    break;

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
    navsto->context = _create_ac_vpp_context(navsto->param);
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
    navsto->context = _create_projection_context(navsto->param);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    return NULL;
  }

  /* Set the static variable */
  cs_navsto_system = navsto;

  return navsto;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the Navier-Stokes system
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_destroy(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL)
    return;

  if (navsto->n_velocity_ic_defs > 0) {

    for (short int i = 0; i < navsto->n_velocity_ic_defs; i++)
      navsto->velocity_ic_defs[i] = cs_xdef_free(navsto->velocity_ic_defs[i]);
    BFT_FREE(navsto->velocity_ic_defs);
    navsto->velocity_ic_defs = NULL;

  } // Velocity IC

  if (navsto->n_pressure_ic_defs > 0) {

    for (short int i = 0; i < navsto->n_pressure_ic_defs; i++)
      navsto->pressure_ic_defs[i] = cs_xdef_free(navsto->pressure_ic_defs[i]);
    BFT_FREE(navsto->pressure_ic_defs);
    navsto->pressure_ic_defs = NULL;

  } // Pressure IC

  /*
    Properties, advection fields, equations and fields are all destroyed
    respectively inside cs_property_destroy_all(),
    cs_advection_field_destroy_all(), cs_equation_destroy_all() and
    cs_field_destroy_all()
  */

  cs_navsto_param_t  *nsp = navsto->param;

  /* Free the context according to the model choice */
  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_UZAWA:
    navsto->context = _free_uzawa_context(nsp, navsto->context);
    break;

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    navsto->context = _free_ac_context(nsp, navsto->context);
    break;

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
    navsto->context = _free_ac_vpp_context(nsp, navsto->context);
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
    navsto->context = _free_projection_context(nsp, navsto->context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;
  }

  /* Destroy the context related to the discretization scheme */
  navsto->free(nsp);

  /* Set of numerical parameters */
  navsto->param = cs_navsto_param_free(nsp);

  BFT_FREE(navsto);
  cs_navsto_system = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Recover the structure storing the parameters for the Navier--Stokes
 *         system
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_system_get_param(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL)
    return NULL;

  return navsto->param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes system
 *         At this stage, numerical settings should be completely determined
 *         but connectivity and geometrical information is not yet available.
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_init_setup(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  cs_navsto_param_t  *nsp = navsto->param;

  /* Create if needed velocity and pressure fields */

  int  location_id = -1; // initialize values to avoid a warning
  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;

  const bool has_previous = cs_navsto_param_is_steady(nsp) ? false:true;

  if (!has_previous)
    field_mask |= CS_FIELD_STEADY;

  switch (nsp->space_scheme) {

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
  case CS_SPACE_SCHEME_HHO_P1:
  case CS_SPACE_SCHEME_HHO_P2:
    location_id = cs_mesh_location_get_id_by_name("cells");
    break; /* Face-based scheme family */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space discretization scheme.", __func__);
  }

  /* Nothing to do if the velocity or pressure fields have been already
     created */
  navsto->velocity = cs_field_find_or_create("velocity",
                                             field_mask,
                                             location_id,
                                             3, /* dimension */
                                             has_previous);

  navsto->pressure = cs_field_find_or_create("pressure",
                                             field_mask,
                                             location_id,
                                             1, /* dimension */
                                             has_previous);

  /* TODO: temperature for the energy equation */

  /* Setup data according to the type of coupling */
  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_UZAWA:
    _uzawa_init_setup(navsto);
    break;

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    _ac_init_setup(navsto);
    break;

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
    _ac_vpp_init_setup(navsto);
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
    _projection_init_setup(navsto);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last step of the setup of the Navier-Stokes system
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_finalize_setup(const cs_cdo_connect_t     *connect,
                                const cs_cdo_quantities_t  *quant)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  assert(connect != NULL && quant != NULL);
  if (navsto == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  cs_navsto_param_t  *nsp = navsto->param;

  /* Avoid an error if no definition is given for the mandatory physical
     properties */
  cs_real_t  one = 1.0;
  if (navsto->density->n_definitions == 0) /* Not set by the user */
    cs_property_def_iso_by_value(navsto->density,
                                 NULL, /* all cells */
                                 one);

  if (navsto->lami_viscosity->n_definitions == 0) /* Not set by the user */
    cs_property_def_iso_by_value(navsto->lami_viscosity,
                                 NULL, /* all cells */
                                 one);

  /* Set functions according to the discretization scheme */
  switch (nsp->space_scheme) {

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
    {
      /* Setup data according to the type of coupling */
      switch (nsp->coupling) {

      case CS_NAVSTO_COUPLING_UZAWA:
        navsto->init = cs_cdofb_navsto_init_uzawa_context;
        navsto->compute = cs_cdofb_navsto_uzawa_compute;

        _uzawa_last_setup(connect, quant, navsto);
        break;

      case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
        navsto->init = cs_cdofb_navsto_init_ac_context;
        navsto->compute = cs_cdofb_navsto_ac_compute;

        _ac_last_setup(connect, quant, navsto);
        break;

      case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
        navsto->init = cs_cdofb_navsto_init_ac_vpp_context;
        navsto->compute = cs_cdofb_navsto_ac_vpp_compute;

        _ac_vpp_last_setup(connect, quant, navsto);
        break;

      case CS_NAVSTO_COUPLING_PROJECTION:
        navsto->init = cs_cdofb_navsto_init_proj_context;
        navsto->compute = cs_cdofb_navsto_proj_compute;

        _projection_last_setup(connect, quant, navsto);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
        break;

      }

      navsto->free = cs_cdofb_navsto_free_context;

    }
    break; /* Lowest-order face-based schemes */

  case CS_SPACE_SCHEME_HHO_P1:
  case CS_SPACE_SCHEME_HHO_P2:
    {
      //TODO: set function pointers

      /* Setup data according to the type of coupling */
      switch (nsp->coupling) {

      case CS_NAVSTO_COUPLING_UZAWA:
        _uzawa_last_setup(connect, quant, navsto);
        break;

      case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
        _ac_last_setup(connect, quant, navsto);
        break;

      case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
        _ac_vpp_last_setup(connect, quant, navsto);
        break;

      case CS_NAVSTO_COUPLING_PROJECTION:
        _projection_last_setup(connect, quant, navsto);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
        break;

      }

    }
    break; /* HHO schemes */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space discretization scheme.", __func__);
  }

  /* Add default post-processing related to the Navier-Stokes system */
  cs_post_add_time_mesh_dep_output(cs_navsto_system_extra_post, navsto);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_velocity_ic_by_value(const char    *z_name,
                                   cs_real_t     *val)
{
  assert(cs_navsto_system != NULL);
  cs_navsto_system_t *nss = cs_navsto_system;

  /* Add a new cs_xdef_t structure */
  int z_id = 0;
  if (z_name != NULL && strlen(z_name) > 0)
    z_id = (cs_volume_zone_by_name(z_name))->id;

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        3, // dim
                                        z_id,
                                        CS_FLAG_STATE_UNIFORM, // state flag
                                        meta_flag,
                                        val);

  int  new_id = nss->n_velocity_ic_defs;
  nss->n_velocity_ic_defs += 1;
  BFT_REALLOC(nss->velocity_ic_defs, nss->n_velocity_ic_defs, cs_xdef_t *);
  nss->velocity_ic_defs[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the pressure unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_pressure_ic_by_value(const char    *z_name,
                                   cs_real_t     *val)
{
  assert(cs_navsto_system != NULL);
  cs_navsto_system_t *nss = cs_navsto_system;

  /* Add a new cs_xdef_t structure */
  int z_id = 0;
  if (z_name != NULL && strlen(z_name) > 0)
    z_id = (cs_volume_zone_by_name(z_name))->id;

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        1, // dim
                                        z_id,
                                        CS_FLAG_STATE_UNIFORM, // state flag
                                        meta_flag,
                                        val);

  int  new_id = nss->n_pressure_ic_defs;
  nss->n_pressure_ic_defs += 1;
  BFT_REALLOC(nss->pressure_ic_defs, nss->n_pressure_ic_defs, cs_xdef_t *);
  nss->pressure_ic_defs[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unkowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_velocity_ic_by_analytic(const char             *z_name,
                                      cs_analytic_func_t     *analytic,
                                      void                   *input)
{
  assert(cs_navsto_system != NULL);
  cs_navsto_system_t *nss = cs_navsto_system;

  /* Add a new cs_xdef_t structure */
  int z_id = 0;
  if (z_name != NULL && strlen(z_name) > 0)
    z_id = (cs_volume_zone_by_name(z_name))->id;

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_input_t  anai = {.func = analytic,
                                    .input = input };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        3, // dim
                                        z_id,
                                        0, // state flag
                                        meta_flag,
                                        &anai);

  int  new_id = nss->n_velocity_ic_defs;
  nss->n_velocity_ic_defs += 1;
  BFT_REALLOC(nss->velocity_ic_defs, nss->n_velocity_ic_defs, cs_xdef_t *);
  nss->velocity_ic_defs[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the pressure unkowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_pressure_ic_by_analytic(const char             *z_name,
                                      cs_analytic_func_t     *analytic,
                                      void                   *input)
{
  assert(cs_navsto_system != NULL);
  cs_navsto_system_t *nss = cs_navsto_system;

  /* Add a new cs_xdef_t structure */
  int z_id = 0;
  if (z_name != NULL && strlen(z_name) > 0)
    z_id = (cs_volume_zone_by_name(z_name))->id;

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_input_t  anai = {.func = analytic,
                                    .input = input };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        1, // dim
                                        z_id,
                                        0, // state flag
                                        meta_flag,
                                        &anai);

  int  new_id = nss->n_pressure_ic_defs;
  nss->n_pressure_ic_defs += 1;
  BFT_REALLOC(nss->pressure_ic_defs, nss->n_pressure_ic_defs, cs_xdef_t *);
  nss->pressure_ic_defs[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the context structure used to build the algebraic system
 *         This is done after the setup step.
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_initialize(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;
  const cs_navsto_param_t *nsp = navsto->param;

  if (navsto == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));
  assert(nsp != NULL);
  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space discretization scheme.", __func__);

  navsto->init(nsp, navsto->context);

  const cs_param_dof_reduction_t  red = nsp->dof_reduction_mode;

  /* Initial conditions for the velocity */
  if (navsto->n_velocity_ic_defs > 0) {

    assert(navsto->velocity_ic_defs != NULL);

    /* Initialize cell-based array */
    const cs_flag_t c_dof_flag = CS_FLAG_VECTOR | cs_flag_primal_cell;
    cs_field_t *field  = cs_field_by_name("velocity");
    cs_real_t  *c_values = field->val;

    for (int def_id = 0; def_id < navsto->n_velocity_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      cs_xdef_t  *def = navsto->velocity_ic_defs[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(c_dof_flag, def, c_values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_xdef_set_quadrature(def, nsp->qtype);
        switch (red) {
        case CS_PARAM_REDUCTION_DERHAM:
          cs_evaluate_potential_by_analytic(c_dof_flag, def, c_values);
          break;
        case CS_PARAM_REDUCTION_AVERAGE:
          cs_evaluate_average_on_cells_by_analytic(c_dof_flag, def, c_values);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    _(" Incompatible reduction for the field %s.\n"),
                    field->name);

        } // Switch on possible reduction types
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Incompatible way to initialize the field %s.\n"),
                  field->name);

      } // Switch on possible type of definition

    } // Loop on definitions

    /* Initialize face-based array */
    if (true) /* DEBUG */
      _init_face_velocity_values();

  } // If velocity IC

  /* Initial conditions for the pressure */
  if (navsto->n_pressure_ic_defs > 0) {

    assert(navsto->pressure_ic_defs != NULL);
    cs_field_t *field  = cs_field_by_name("pressure");
    cs_real_t  *values = field->val;
    const cs_flag_t   dof_flag = CS_FLAG_SCALAR | cs_flag_primal_cell;

    for (int def_id = 0; def_id < navsto->n_pressure_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      cs_xdef_t  *def = navsto->pressure_ic_defs[def_id];

      /* Initialize face-based array */
      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(dof_flag, def, values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_xdef_set_quadrature(def, nsp->qtype);
        switch (red) {
        case CS_PARAM_REDUCTION_DERHAM:
          cs_evaluate_potential_by_analytic(dof_flag, def, values);
          break;
        case CS_PARAM_REDUCTION_AVERAGE:
          cs_evaluate_average_on_cells_by_analytic(dof_flag, def, values);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    _(" Incompatible reduction for the field %s.\n"),
                    field->name);

        } // Switch on possible reduction types
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Incompatible way to initialize the field %s.\n"),
                  field->name);
        break;

      } // Switch on possible type of definition

    } // Loop on definitions

  } // If pressure IC

  /* TODO: Set the initial condition for variables not directly related
     to an equation */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build, solve and update the Navier-Stokes system
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      dt_cur     current value of the time step
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_compute(const cs_mesh_t              *mesh,
                         double                        dt_cur)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  /* Build and solve the Navier-Stokes system */
  navsto->compute(mesh, dt_cur, navsto->param, navsto->context);

  /* TODO: Update the variable states */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the Navier-Stokes system.
 *         The prototype of this function is fixed since it is a function
 *         pointer defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
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
cs_navsto_system_extra_post(void                      *input,
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

  cs_navsto_system_t  *navsto = (cs_navsto_system_t *)input;

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_log_setup(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSummary of the Navier-Stokes system\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);

  /* Main set of numerical parameters */
  cs_navsto_param_log(navsto->param);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
