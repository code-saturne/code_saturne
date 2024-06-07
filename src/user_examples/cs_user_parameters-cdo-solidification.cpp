/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

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
 * \file cs_user_parameters-cdo-solidification.c
 *
 * \brief User functions for setting a calculation using the solidification
 *        module with CDO schemes
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

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

  cs_boundary_add(bdy, CS_BOUNDARY_WALL, "left");
  cs_boundary_add(bdy, CS_BOUNDARY_WALL, "right");
  cs_boundary_add(bdy, CS_BOUNDARY_WALL, "top");
  cs_boundary_add(bdy, CS_BOUNDARY_WALL, "bottom");

  /* Activate CDO/HHO module so that main additional structure are built */

  cs_param_cdo_mode_set(CS_PARAM_CDO_MODE_ONLY);

  /* ===============================
     Define the solidification model
     =============================== */

  /* 1. Activate the solidification module */
  /* ------------------------------------- */

  /*! [param_cdo_activate_solidification_binary] */
  {
    /* For the solidification module:
       cs_solidification_activate(solidification_model_type,
                                  solid_option_flag,
                                  solid_post_flag,
                                  boundaries,
                                  navsto_model,
                                  navsto_model_flag,
                                  navsto_coupling_type,
                                  navsto_post_flag);

       If a flag is set to 0, then there is no option to add.
       To add options to a flag:
       flag = option1 | option2 | option3 | ...
    */

    cs_flag_t  solid_option_flag = 0;

    cs_flag_t  solid_post_flag = 0 |
      CS_SOLIDIFICATION_POST_SOLIDIFICATION_RATE |
      CS_SOLIDIFICATION_POST_SEGREGATION_INDEX;

    cs_flag_t  navsto_model_flag = 0;

    cs_flag_t  navsto_post_flag = 0     |
      CS_NAVSTO_POST_VELOCITY_DIVERGENCE |
      CS_NAVSTO_POST_MASS_DENSITY;

    /* Activate the solidification module with a binary alloy model (the
       Navier-Stokes and the thermal modules are also activated in back-end) */

    cs_solidification_activate(/* Main solidification model */
                               CS_SOLIDIFICATION_MODEL_BINARY_ALLOY,
                               /* Solidification options */
                               solid_option_flag,
                               /* Solidification automatic post options */
                               solid_post_flag,
                               /* NavSto parameters */
                               domain->boundaries,
                               CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES,
                               navsto_model_flag,
                               CS_NAVSTO_COUPLING_MONOLITHIC,
                               navsto_post_flag);

  }
  /*! [param_cdo_activate_solidification_binary] */

  /*! [param_cdo_solidification_set_binary_alloy] */
  {
    /* Physical data for the settings a binay alloy model. The main
       solidification model is CS_SOLIDIFICATION_MODEL_BINARY_ALLOY */

    cs_real_t  T0 = 0.5, beta_t = 0.01;
    cs_real_t  conc0 = 1.0, beta_c = 0.01;
    cs_real_t  ml = -0.1, kp = 0.1;
    cs_real_t  t_eutec = -0.1, t_melt = 0.2;
    cs_real_t  diff_val = 0;
    cs_real_t  latent_heat = 5;

    /* Set the parameters for the binary alloy model */

    cs_solidification_set_binary_alloy_model("C_solute",
                                             "C_bulk",
                                             /* Boussinesq approximation */
                                             beta_t,
                                             T0,
                                             beta_c,
                                             conc0,
                                             /* Phase diagram settings */
                                             kp,
                                             ml,
                                             t_eutec,
                                             t_melt,
                                             /* Solute transport equation */
                                             diff_val,
                                             /* Physical constants */
                                             latent_heat);

    /* Set the parameters related to the Kozeny-Carman relation */

    cs_solidification_set_kozeny_carman_parameters(4.55945,  // Kozeny constant
                                                   2.0,      // Tortuosity
                                                   0.33541); // s_das
  }
  /*! [param_cdo_solidification_set_binary_alloy] */

  /*! [param_cdo_activate_solidification_voller] */
  {
    cs_flag_t  solid_option_flag = 0;
    cs_flag_t  solid_post_flag = CS_SOLIDIFICATION_POST_CELL_STATE;
    cs_flag_t  navsto_model_flag = 0;
    cs_flag_t  navsto_post_flag = CS_NAVSTO_POST_VELOCITY_DIVERGENCE;

    /* Activate the solidification module with a binary alloy model (the
       Navier-Stokes and the thermal modules are also activated in back-end) */

    cs_solidification_activate(/* Main solidification model */
                               CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87,
                               /* Solidification options */
                               solid_option_flag,
                               /* Solidification automatic post options */
                               solid_post_flag,
                               /* NavSto parameters */
                               domain->boundaries,
                               CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES,
                               navsto_model_flag,
                               CS_NAVSTO_COUPLING_MONOLITHIC,
                               navsto_post_flag);

  }
  /*! [param_cdo_activate_solidification_voller] */

  /*! [param_cdo_solidification_set_voller] */
  {
    /* Physical data for the settings a Voller & Prakash model (i.e. without
       segregation, only solidification is taken into account). The main model
       is set to CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87 */

    cs_real_t  T0 = 0.5, beta_t = 0.01;
    cs_real_t  t_solidus = -0.1, t_liquidus = 0.1;
    cs_real_t  latent_heat = 5;

    /* Set the parameters for the Voller & Prakash model */

    cs_solidification_set_voller_model(/* Boussinesq approximation */
                                       beta_t,
                                       T0,
                                       /* Phase diagram settings */
                                       t_solidus,
                                       t_liquidus,
                                       /* Physical constants */
                                       latent_heat);

    /* Set the parameters related to the Kozeny-Carman relation */

    cs_solidification_set_kozeny_carman_parameters(4.55945,  // Kozeny constant
                                                   2.0,      // Tortuosity
                                                   0.33541); // s_das
  }
  /*! [param_cdo_solidification_set_voller] */

  /*! [param_cdo_activate_solidification_voller_no_vel] */
  {
    cs_flag_t  solid_option_flag = CS_SOLIDIFICATION_NO_VELOCITY_FIELD;
    cs_flag_t  solid_post_flag = CS_SOLIDIFICATION_POST_ENTHALPY;
    cs_flag_t  navsto_model_flag = 0;
    cs_flag_t  navsto_post_flag = 0;

    /* Activate the solidification module with a binary alloy model (the
       Navier-Stokes and the thermal modules are also activated in back-end) */

    cs_solidification_activate(/* Main solidification model */
                               CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87,
                               /* Solidification options */
                               solid_option_flag,
                               /* Solidification automatic post options */
                               solid_post_flag,
                               /* NavSto parameters */
                               domain->boundaries,
                               CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES,
                               navsto_model_flag,
                               CS_NAVSTO_COUPLING_MONOLITHIC,
                               navsto_post_flag);
  }
  /*! [param_cdo_activate_solidification_voller_no_vel] */

  /*! [param_cdo_solidification_set_voller_no_vel] */
  {
    /* Physical data for the settings a Voller & Prakash model (i.e. without
       segregation, only solidification is taken into account). The main model
       is set to CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87 but the optional
       flag CS_SOLIDIFICATION_NO_VELOCITY_FIELD has been set. Then, one uses a
       simplified version of the previous function. */

    cs_real_t  t_solidus = -0.1, t_liquidus = 0.1;
    cs_real_t  latent_heat = 5;

    /* Set the parameters for the Voller & Prakash model when there is no
       resolution of the velocity */

    cs_solidification_set_voller_model_no_velocity(/* Phase diagram settings*/
                                                   t_solidus,
                                                   t_liquidus,
                                                   /* Physical constants */
                                                   latent_heat);
  }
  /*! [param_cdo_solidification_set_voller_no_vel] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t    *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);

  /*! [param_cdo_solidification_set_strategy] */
  {
    cs_solidification_set_strategy(CS_SOLIDIFICATION_STRATEGY_PATH);
  }
  /*! [param_cdo_solidification_set_strategy] */

  /*! [param_cdo_solidification_nl_voller_advanced] */
  {
    /* If a non-linear Voller model has been activated */

    cs_solidification_voller_t
      *model_struct = cs_solidification_get_voller_struct();

    /* Drive the convergence of the non-linear algorithm to update the thermal
     * source term. */

    model_struct->nl_algo->cvg_param.n_max_iter = 20;
    model_struct->nl_algo->cvg_param.rtol = 1e-3; /* relative tolerance */
  }
  /*! [param_cdo_solidification_nl_voller_advanced] */

  /*! [param_cdo_solidification_binary_advanced] */
  {
    /* If a "binary alloy" model has been activated */

    cs_solidification_binary_alloy_t
      *model_struct = cs_solidification_get_binary_alloy_struct();

    int  n_iter_max = 10;
    double  delta_eps = 1e-3;

    /* Drive the convergence of the coupled system (solute transport and thermal
     * equation) with respect to the following criteria (taken from Voller and
     * Swaminathan'91)
     *   max_{c\in C} |Temp^(k+1) - Temp^(k)| < delta_tolerance
     *   max_{c\in C} |Cbulk^(k+1) - Cbulk*^(k)| < delta_tolerance
     *   n_iter < n_iter_max
     */

    model_struct->n_iter_max = n_iter_max;
    model_struct->delta_tolerance = delta_eps;
  }
  /*! [param_cdo_solidification_binary_advanced] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify the elements such as properties, advection fields,
 *         user-defined equations and modules which have been previously added.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);

  /*! [param_cdo_solidification_properties] */
  {
    /* All the following properties are isotropic and set on all the mesh cells
       (this implies NULL for the second argument). If the property is
       piecewise constant, then replace NULL by the name of a volum zone. */

    /* Mass density (kg.m^-3) */

    cs_property_t  *rho = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
    cs_real_t  rho0 = 1;
    cs_property_def_iso_by_value(rho, NULL, rho0);

    /* Laminar dynamic viscosity (Pa.s) */

    cs_property_t  *mu = cs_property_by_name(CS_NAVSTO_LAM_VISCOSITY);
    cs_real_t  mu0 = 1;
    cs_property_def_iso_by_value(mu, NULL, mu0);

    /* Thermal heat capacity */

    cs_property_t  *cp = cs_property_by_name(CS_THERMAL_CP_NAME);
    cs_real_t  cp0 = 1.;
    cs_property_def_iso_by_value(cp, NULL, cp0);

    /* Thermal conductivity */

    cs_property_t  *lambda = cs_property_by_name(CS_THERMAL_LAMBDA_NAME);
    cs_real_t  lambda0 = 0.001;
    cs_property_def_iso_by_value(lambda, NULL, lambda0);
  }
  /*! [param_cdo_solidification_properties] */


  /*! [param_cdo_solidification_thermal_eq] */
  {
    cs_real_t  t_ref = 0.5;

    cs_equation_param_t  *th_eqp = cs_equation_param_by_name(CS_THERMAL_EQNAME);

    /* Set the initial value for the temperature */

    cs_equation_add_ic_by_value(th_eqp, NULL, &t_ref);

    /* Set the value of the boundary conditions.
     *
     * The other boundary zones are associated to the default boundary (by
     * default, the thermal equation is associated to a no flux (Homogeneous
     * Neumann boundary condition) */

    cs_real_t  Th = t_ref, Tc = -t_ref;
    cs_equation_add_bc_by_value(th_eqp, CS_BC_DIRICHLET, "left", &Tc);
    cs_equation_add_bc_by_value(th_eqp, CS_BC_DIRICHLET, "right", &Th);
  }
  /*! [param_cdo_solidification_thermal_eq] */

  /*! [param_cdo_solidification_solute_eq] */
  {
    cs_real_t  c_ref = 0.5;

    /* One assumes that an equation called "C_solute" has been added */

    cs_equation_param_t  *c_eqp = cs_equation_param_by_name("C_solute");

    /* Set the initial value for the solute concentration to the reference
       value */

    cs_equation_add_ic_by_value(c_eqp, NULL, &c_ref);

  }
  /*! [param_cdo_solidification_solute_eq] */

}
/*----------------------------------------------------------------------------*/

END_C_DECLS
