/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
 * \file cs_user_parameters-cdo-gwf.c
 *
 * \brief User functions for setting a calculation using the groundwater flow
 *        module with CDO schemes
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/* Permeability in each subdomain */

static const double k1 = 1e5, k2 = 1;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the initial condition for the tracer attached to the GWF problem
 *         Generic function pointer for an evaluation relying on an analytic
 *         function.
 *         pt_ids is optional. If not NULL, it enables to access to the xyz
 *         array with an indirection. The same indirection can be applied to
 *         fill retval if dense_output is set to false.
 *
 * \param[in]      time          when ?
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      pt_ids        list of elements ids (in coords and retval)
 * \param[in]      xyz           where ? Coordinates array
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static inline void
get_tracer_ic(cs_real_t          time,
              cs_lnum_t          n_elts,
              const cs_lnum_t   *pt_ids,
              const cs_real_t   *xyz,
              bool               dense_output,
              void              *input,
              cs_real_t         *retval)
{
  CS_UNUSED(input);

  /* Physical parameters */

  const double  magnitude = 2*k1/(k1 + k2);
  const double  x_front = magnitude * time;

  /* Fill retval array */

  for (cs_lnum_t  i = 0; i < n_elts; i++) {

    const cs_lnum_t  id = (pt_ids == NULL) ? i : pt_ids[i];
    const cs_lnum_t  r_id = (dense_output) ? i : id;
    const double  x = xyz[3*id];

    if (x <= x_front)
      retval[r_id] = 1;
    else
      retval[r_id] = 0;

  } /* Loop on selected elements */

}

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
  cs_domain_t  *domain = cs_glob_domain;

  /* ======================
     Boundary of the domain
     ====================== */

  cs_boundary_t  *bdy = domain->boundaries;

  /* Choose a boundary by default */

  cs_boundary_set_default(bdy, CS_BOUNDARY_SYMMETRY);

  /* Add new boundaries */

  cs_boundary_add(bdy, CS_BOUNDARY_INLET, "left");
  cs_boundary_add(bdy, CS_BOUNDARY_OUTLET, "right");

  /* Activate CDO/HHO module so that main additional structure are built */

  cs_domain_set_cdo_mode(domain, CS_DOMAIN_CDO_MODE_ONLY);

  /* =========================
     Generic output management
     ========================= */

  cs_domain_set_output_param(domain,
                             -1,      // restart frequency
                             100,     // output log frequency
                             3);      // verbosity (-1: no, 0, ...)

  /* ====================
     Time step management
     ==================== */

  /*
    If there is an inconsistency between the max. number of iteration in
    time and the final physical time, the first condition encountered stops
    the calculation.
  */

  cs_domain_set_time_param(domain,
                           500,     // nt_max or -1 (automatic)
                           -1.);    // t_max or < 0. (automatic)

  /* Define the value of the time step
     >> cs_domain_def_time_step_by_value(domain, dt_val);
     >> cs_domain_def_time_step_by_func(domain, dt_func, context);
     context may be NULL.
     This second way to define the time step enable complex definitions.
  */

  cs_domain_def_time_step_by_value(domain, 1.0);

  /* =================================
     Define the groundwater flow model
     ================================= */

  /* 1. Activate groundwater flow module */
  /* ----------------------------------- */

  /*! [param_cdo_activate_gwf] */
  {
    /* For the groundwater flow module:
       cs_gwf_activate(permeability_type,
                       model_type,
                       option_flag);

     - permeability_type is one of the following keywords:
       CS_PROPERTY_ISO, CS_PROPERTY_ORTHO or CS_PROPERTY_ANISO

     - model_type is one among
       CS_GWF_MODEL_SATURATED_SINGLE_PHASE
       CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE

     - option_flag can be defined from the following flags
       CS_GWF_GRAVITATION
       CS_GWF_RICHARDS_UNSTEADY
       CS_GWF_SOIL_PROPERTY_UNSTEADY
       CS_GWF_SOIL_ALL_SATURATED
       or 0 if there is no option to set
    */

    cs_flag_t  option_flag = 0;

    cs_gwf_activate(CS_PROPERTY_ISO,
                    CS_GWF_MODEL_SATURATED_SINGLE_PHASE,
                    option_flag);
  }
  /*! [param_cdo_activate_gwf] */

  /*! [param_cdo_activate_gwf_b] */
  {
    /* Take into account the gravity effect */

    cs_flag_t  option_flag = CS_GWF_GRAVITATION;

    /* In this case, the gravity vector has to be defined (either using the GUI
       or in cs_user_parameters() function */

    cs_gwf_activate(CS_PROPERTY_ISO,
                    CS_GWF_MODEL_SATURATED_SINGLE_PHASE,
                    option_flag);
  }
  /*! [param_cdo_activate_gwf_b] */

  /* 2. Add and define soils (must be done before adding tracers) */
  /* ----------------------- */

  /*! [param_cdo_gwf_add_define_simple_soil] */
  {
    /* Example 1: two "saturated" soils (simple case) */

    const cs_real_t  theta_s = 1;
    const cs_real_t  bulk_density = 1.0; /* useless if no tracer is
                                            considered */

    /* Two (volume) zones have be defined called "soil1" and "soil2". */

    cs_gwf_soil_t  *s1 = cs_gwf_add_soil("soil1",
                                         CS_GWF_SOIL_SATURATED,
                                         theta_s,
                                         bulk_density);

    cs_gwf_soil_set_iso_saturated(s1, k1);  /* saturated permeability */

    /* For "simple" soil definitions, definition can be made here. For more
       complex definition, please use \ref cs_user_gwf_setup */

    cs_gwf_soil_t  *s2 = cs_gwf_add_soil("soil2",
                                         CS_GWF_SOIL_SATURATED,
                                         theta_s,
                                         bulk_density);

    cs_gwf_soil_set_iso_saturated(s2, k2);  /* saturated permeability */
  }
  /*! [param_cdo_gwf_add_define_simple_soil] */

  /*! [param_cdo_gwf_add_user_soil] */
  {
    /* Example 2: Add a new user-defined soil for all the cells */

    const cs_real_t  theta_s = 0.9;         /* max. liquid saturation */
    const cs_real_t  bulk_density = 1800.0; /* useless if no tracer is
                                               considered */

    cs_gwf_add_soil("cells", CS_GWF_SOIL_USER, theta_s, bulk_density);

    /* The explicit definition of this soil is done in \ref cs_user_gwf_setup */
  }
  /*! [param_cdo_gwf_add_user_soil] */

  /*! [param_cdo_gwf_get_soil] */
  {
    /* If a soil structure has been added but not already defined, one can
     * retrieve it thanks to \ref cs_gwf_soil_by_name
     *
     * The name of the soil is the same as the name of the volume zone used at
     * the creation of the soil
     */

    cs_gwf_soil_t  *s1 = cs_gwf_soil_by_name("soil1");
  }
  /*! [param_cdo_gwf_get_soil] */

  /* 3. Add and define tracer equations */
  /* ---------------------------------- */

  /*! [param_cdo_gwf_add_tracer] */
  {
    /*
      Add a tracer equation which is unsteady and convected by the darcean flux
      This implies the creation of a new equation called eqname along with a
      new field called varname.

      For standard tracer:
       cs_gwf_add_tracer(eqname, varname);
       For user-defined tracer
       cs_gwf_add_user_tracer(eqname, varname, setup_func);
    */

    cs_gwf_tracer_model_t  model = 0; /* default model */

    cs_gwf_tracer_t  *tr = cs_gwf_add_tracer(model,       /* tracer model */
                                             "Tracer_01", /* eq. name */
                                             "C1");       /* var. name */

    /* For "simple" tracer definitions, definition can be made here. For more
     * complex definition, please use \ref cs_user_gwf_setup
     *
     * The parameters defining the tracer behavior can be set soil by soil
     * (give the soil name as the second argument) or to all soils in one call
     * (give NULL as the second argument)
     */

    cs_gwf_set_main_tracer_param(tr,
                                 NULL,     /* soil name or NULL for all */
                                 0.,       /* water molecular diffusivity */
                                 0., 0.,   /* alpha (longi. and transvesal) */
                                 0.,       /* distribution coef. */
                                 0.);      /* 1st order decay coef. */
  }
  /*! [param_cdo_gwf_add_tracer] */

  /*! [param_cdo_gwf_get_tracer] */
  {
    /* If the tracer structure has been added but not already defined, one can
     * retrieve it thanks to \ref cs_gwf_tracer_by_name
     *
     * The name of the tracer is the same as the name of the associated
     * equation given at the creation of the tracer
     */

    cs_gwf_tracer_t *tr = cs_gwf_tracer_by_name("Tracer_01");
  }
  /*! [param_cdo_gwf_get_tracer] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify the elements such as properties, advection fields,
 *           user-defined equations and modules which have been previously
 *           added.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* 1. Set the Richards equation */
  /* ---------------------------- */

  cs_equation_param_t  *r_eqp = cs_equation_param_by_name("Richards");

  /* Define the boundary conditions  */

  cs_real_t  right_val = 0.0, left_val = 1.0;

  cs_equation_add_bc_by_value(r_eqp,
                              CS_PARAM_BC_DIRICHLET,
                              "left",       /* boundary zone name */
                              &left_val);   /* value to set */

  cs_equation_add_bc_by_value(r_eqp,
                              CS_PARAM_BC_DIRICHLET,
                              "right",       /* boundary zone name */
                              &right_val);   /* value to set */

  /* 2. Set the Richards equation */
  /* ---------------------------- */

  cs_equation_param_t  *t_eqp = cs_equation_param_by_name("Tracer_01");

  cs_equation_add_bc_by_value(t_eqp,
                              CS_PARAM_BC_DIRICHLET,
                              "left",      /* boundary zone name */
                              &left_val);  /* value to set */

  /* Define the initial condition with an analytic function (if nothing is
     done, the default initialization is zero) */

  cs_equation_add_ic_by_analytic(t_eqp, "cells", get_tracer_ic, NULL);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
