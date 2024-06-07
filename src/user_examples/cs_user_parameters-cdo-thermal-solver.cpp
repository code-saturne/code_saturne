/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* code_saturne version 8.1-alpha */

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
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*!
 * \brief Set the initial temperature in the computational domain
 *
 * For the calling function, elt_ids is optional. If not NULL, the coords
 * array should be accessed with an indirection. The same indirection can
 * be applied to fill retval if dense_output is set to false.
 *
 * \param[in]      time          when ?
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids (in coords and retval)
 * \param[in]      coords        where ?
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*--------------------------------------------------------------------------*/

/*! [param_cdo_initial_temperature_function] */
static void
_initial_temperature(cs_real_t            time,
                     cs_lnum_t            n_elts,
                     const cs_lnum_t     *elt_ids,
                     const cs_real_t     *coords,
                     bool                 dense_output,
                     void                *input,
                     cs_real_t           *retval)
{
  CS_NO_WARN_IF_UNUSED(time);
  CS_NO_WARN_IF_UNUSED(input);

  if (elt_ids == NULL) { /* No indirection */

    for (cs_lnum_t i = 0; i < n_elts; i++)
      retval[i] = -1 + 2*coords[3*i];            /* T(t=0) = -1 + 2*x */

  }
  else {

    if (dense_output)
      for (cs_lnum_t i = 0; i < n_elts; i++)
        retval[i] = -1 + 2*coords[3*elt_ids[i]]; /* T(t=0) = -1 + 2*x */

    else {

      for (cs_lnum_t i = 0; i < n_elts; i++) {

        const cs_lnum_t  elt_id = elt_ids[i];
        retval[elt_id] = -1 + 2*coords[3*elt_id]; /* T(t=0) = -1 + 2*x */

      }

    } /* dense_output ? */

  } /* elt_ids == NULL ? */
}
/*! [param_cdo_initial_temperature_function] */

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{
  /* Activate the CDO module so that the main additional structures are
     built */

  cs_param_cdo_mode_set(CS_PARAM_CDO_MODE_ONLY);

  /*! [param_cdo_activate_thermal_solver] */
  {
    cs_thermal_system_activate(0,  /* model flag (0=default) */
                               0,  /* numeric flag (0=default) */
                               0); /* post flag (0=no automatic post-process) */
  }
  /*! [param_cdo_activate_thermal_solver] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t   *domain)
{
  /*! [param_cdo_set_time_step] */
  {
    /* Time step management
     * ====================
     * If there is an inconsistency between the max. number of iteration in
     * time and the final physical time, the first condition encountered stops
     * the calculation.
     */

    cs_domain_set_time_param(domain,
                             250,     /* nt_max or -1 (automatic) */
                             -1.);    /* t_max or < 0. (automatic) */

    /* Define the value of the time step. Two functions are available to do
       this:
       1. cs_domain_def_time_step_by_value(domain, dt_val);
       2. cs_domain_def_time_step_by_func(domain, dt_func);

       The second way enables more complex definitions of the time step.
    */

    cs_domain_def_time_step_by_value(domain, 5.0);
  }
  /*! [param_cdo_set_time_step] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify output user parameters.
 *
 * For CDO schemes, this function concludes the setup of properties, equations,
 * source terms...
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /*! [param_cdo_define_thermal_properties] */
  {
    /* Define the predefined properties */
    /* -------------------------------- */

    cs_property_t  *lambda = cs_property_by_name(CS_THERMAL_LAMBDA_NAME);

    cs_property_def_iso_by_value(lambda, "D1", 1.0);
    cs_property_def_iso_by_value(lambda, "D2", 2.5);

    cs_property_t  *cp = cs_property_by_name(CS_THERMAL_CP_NAME);
    cs_property_def_iso_by_value(cp, NULL, 1); /* NULL means all cells */

    cs_property_t  *rho = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
    cs_property_def_iso_by_value(rho, NULL, 1979);

  }
  /*! [param_cdo_define_thermal_properties] */

  /*! [param_cdo_define_thermal_bc] */
  {
    /* Define the boundary conditions */
    /* ------------------------------ */

    cs_equation_param_t  *eqp = cs_equation_param_by_name(CS_THERMAL_EQNAME);

    double  Tleft = -1;           /* Dirichlet BC value */

    cs_equation_add_bc_by_value(eqp,
                                CS_BC_DIRICHLET,
                                "left",
                                &Tleft);

    double  WallFlux = 2;        /* Neumann BC value */

    cs_equation_add_bc_by_value(eqp,
                                CS_BC_NEUMANN,
                                "wall",
                                &WallFlux);

    /* Robin BC values are defined such that:
     *  normal_flux = alpha * (T - T0) + phi0 */

    double  RobinBCs[3] = {10,    /* alpha */
                            1,    /* T0 */
                           -2};   /* phi0 */

    cs_equation_add_bc_by_value(eqp,
                                CS_BC_ROBIN,
                                "right",
                                RobinBCs);

    /* The remaining part of the boundary faces are set to a homogeneous Neumann
       BCs (the default BC) */
  }
  /*! [param_cdo_define_thermal_bc] */

  /*! [param_cdo_define_thermal_ic] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name(CS_THERMAL_EQNAME);

    /* Define the initial conditions */
    /* ----------------------------- */

    /* If nothing is done, then a zero value is set by default */

    cs_equation_add_ic_by_analytic(eqp,
                                   NULL,                 /* all cells */
                                   _initial_temperature, /* function def. */
                                   NULL);                /* no context */
  }
  /*! [param_cdo_define_thermal_ic] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
