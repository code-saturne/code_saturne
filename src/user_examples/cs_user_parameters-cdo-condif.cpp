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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-cdo-condif.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Give the explicit definition of the advection field.
 *         Generic function pointer for an evaluation relying on an analytic
 *         function
 *         pt_ids is optional. If not NULL, it enables to access to the coords
 *         array with an indirection. The same indirection can be applied to
 *         fill retval if dense_output is set to false.
 *
 * \param[in]      time          when ?
 * \param[in]      n_pts         number of elements to consider
 * \param[in]      pt_ids        list of elements ids (in coords and retval)
 * \param[in]      xyz           where ? Coordinates array
 * \param[in]      dense_output  perform an indirection in res or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] res           resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_define_adv_field(cs_real_t           time,
                  cs_lnum_t           n_pts,
                  const cs_lnum_t    *pt_ids,
                  const cs_real_t    *xyz,
                  bool                dense_output,
                  void               *input,
                  cs_real_t          *res)
{
  CS_NO_WARN_IF_UNUSED(time);
  CS_NO_WARN_IF_UNUSED(input);

  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const cs_lnum_t  id = (pt_ids == NULL) ? p : pt_ids[p];
    const cs_lnum_t  ii = dense_output ? p : id;
    const cs_real_t  *pxyz = xyz + 3*id;
    cs_real_t  *pres = res + 3*ii;

    pres[0] = pxyz[1] - 0.5;
    pres[1] = 0.5 - pxyz[0];
    pres[2] = pxyz[2] + 1;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Give the explicit definition of the dirichlet boundary conditions
 *         Generic function pointer for an evaluation relying on an analytic
 *         function
 *         pt_ids is optional. If not NULL, it enables to access to the coords
 *         array with an indirection. The same indirection can be applied to
 *         fill retval if dense_output is set to false.
 *
 * \param[in]      time          when ?
 * \param[in]      n_pts         number of elements to consider
 * \param[in]      pt_ids        list of elements ids (in coords and retval)
 * \param[in]      xyz           where ? Coordinates array
 * \param[in]      dense_output  perform an indirection in res or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] res           resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_define_bcs(cs_real_t           time,
            cs_lnum_t           n_pts,
            const cs_lnum_t    *pt_ids,
            const cs_real_t    *xyz,
            bool                dense_output,
            void               *input,
            cs_real_t          *res)
{
  CS_NO_WARN_IF_UNUSED(time);
  CS_NO_WARN_IF_UNUSED(input);

  const double  pi = cs_math_pi;
  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const cs_lnum_t  id = (pt_ids == NULL) ? p : pt_ids[p];
    const cs_lnum_t  ii = dense_output ? p : id;
    const cs_real_t  *_xyz = xyz + 3*id;
    const double  x = _xyz[0], y = _xyz[1], z = _xyz[2];

    res[ii] = 1 + sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+0.25));

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Give the explicit definition of the source term
 *        pt_ids is optional. If not NULL, it enables to access in xyz
 *        at the right location and to fill the resulting array if a dense
 *        output is not requested.
 *        Rely on the generic function pointer for an analytic function
 *
 * \param[in]      time         when ?
 * \param[in]      n_elts       number of elements to consider
 * \param[in]      pt_ids       list of elements ids (to access coords and fill)
 * \param[in]      coords       where ?
 * \param[in]      dense_output true:no indirection, false:apply pt_ids
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] values       result of the function
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_condif_analytic_st] */
static void
_my_source_term(cs_real_t           time,
                cs_lnum_t           n_pts,
                const cs_lnum_t    *pt_ids,
                const cs_real_t    *xyz,
                bool                dense_output,
                void               *input,
                cs_real_t          *values)
{
  CS_NO_WARN_IF_UNUSED(time);

  const double  *pcoefs = (double *)input;
  const double  mu = pcoefs[0];
  const double  pi = pcoefs[1];

  for (cs_lnum_t i = 0; i < n_pts; i++) {

    const cs_lnum_t  id = (pt_ids == NULL) ? i : pt_ids[i];
    const cs_lnum_t  ii = dense_output ? i : id;
    const double  x = xyz[3*id], y = xyz[3*id+1], z = xyz[3*id+2];
    const double  px = pi*x, cpx = cos(px), spx = sin(px);
    const double  py = 2*pi*y, cpy = cos(py), spy = sin(py);
    const double  pz = pi*z, cpz = cos(pz), spz = sin(pz);
    const cs_real_t  gx = cpx*spy*spz;
    const cs_real_t  gy = 2*spx*cpy*spz;
    const cs_real_t  gz = spx*spy*cpz;

    values[ii] = pi*( (y-0.5)*gx + (0.5-x)*gy + (z+1.0)*gz ) + mu*spx*spy*spz;

  }
}
/*! [param_cdo_condif_analytic_st] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy an input data structure.
 *        Complex data structure can be used when a \ref cs_xdef_t structure
 *        is defined by an analytic function, a DoF function or a time
 *        function. Please refer to \ref cs_xdef_analytic_context_t,
 *        \ref cs_xdef_time_func_context_t or \ref cs_xdef_dof_context_t
 *
 * \param[in, out] input    pointer to an input structure associated to a
 *                          context structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_condif_free_input] */
static void *
_free_input(void   *input)
{
  double *_input = (double *)input;

  BFT_FREE(_input);
  input = NULL;

  return NULL;
}
/*! [param_cdo_condif_free_input] */

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
  /* Activate CDO/HHO module so that the main additional structures are built */

  /*! [param_cdo_activation] */
  {
    cs_param_cdo_mode_set(CS_PARAM_CDO_MODE_ONLY);
  }
  /*! [param_cdo_activation] */

  /* ======================
     Boundary of the domain
     ====================== */

  /*! [param_cdo_domain_boundary] */
  {
    cs_domain_t  *domain = cs_glob_domain;
    cs_boundary_t  *bdy = domain->boundaries;

    /* Choose a boundary by default */
    cs_boundary_set_default(bdy, CS_BOUNDARY_WALL);

    /* Add a new boundary
     *  >> cs_domain_add_boundary(boundary, type_of_boundary, zone_name);
     *
     * zone_name is either a predefined one or user-defined one
     * type_of_boundary is one of the following keywords:
     *   CS_BOUNDARY_WALL,
     *   CS_BOUNDARY_INLET,
     *   CS_BOUNDARY_OUTLET,
     *   CS_BOUNDARY_SYMMETRY
     */

    cs_boundary_add(bdy, CS_BOUNDARY_SYMMETRY, "sym");
  }
  /*! [param_cdo_domain_boundary] */

  /*! [param_cdo_wall_distance] */
  {
    /*  Activate predefined module as the computation of the wall distance */
    cs_walldistance_activate();
  }
  /*! [param_cdo_wall_distance] */

  /*! [param_cdo_add_user_equation] */
  {
    cs_equation_add_user("my_equation", /* name of the equation */
                         "my_variable", /* name of the associated variable */
                         1,             /* dimension of the variable */
                         CS_BC_SYMMETRY); /* default BC */

    /* Add a new user equation.
     *   The default boundary condition has to be chosen among:
     *    CS_BC_HMG_DIRICHLET
     *    CS_BC_SYMMETRY
     */

    cs_equation_add_user("AdvDiff.Upw", // equation name
                         "Pot.Upw",     // associated variable field name
                         1,             // dimension of the unknown
                         CS_BC_HMG_DIRICHLET); // default boundary

    cs_equation_add_user("AdvDiff.SG",  // equation name
                         "Pot.SG",      // associated variable field name
                         1,             // dimension of the unknown
                         CS_BC_HMG_DIRICHLET); // default boundary
  }
  /*! [param_cdo_add_user_equation] */

  /* ========================================
     Add material properties
     ======================================== */

  /*! [param_cdo_add_user_properties] */
  {
    /* For an anistropic property (tensor-valued) */

    cs_property_add("conductivity",      /* property name */
                    CS_PROPERTY_ANISO);  /* type of material property */

    /* For an isotropic property (scalar-valued) */

    cs_property_add("rho.cp",            /* property name */
                    CS_PROPERTY_ISO);    /* type of material property */

  }
  /*! [param_cdo_add_user_properties] */

  /*! [param_cdo_add_user_properties_opt] */
  {
    /* Retrieve a property named "conductivity"  */
    cs_property_t  *pty = cs_property_by_name("conductivity");

    /* Activate the computation of the Fourier number for this property */
    cs_property_set_option(pty, CS_PTYKEY_POST_FOURIER);
  }
  /*! [param_cdo_add_user_properties_opt] */

  /* ========================================
     Add advection fields
     ======================================== */

  /*! [param_cdo_add_user_adv_field] */
  {
    /* Add a user-defined advection field named "adv_field"  */

    cs_adv_field_t  *adv = cs_advection_field_add_user("adv_field");

    assert(adv != NULL);
  }
  /*! [param_cdo_add_user_adv_field] */

  /*! [param_cdo_add_adv_field] */
  {
    /* Add a user-defined advection field named "adv_field"  */
    cs_advection_field_status_t  adv_status =
      CS_ADVECTION_FIELD_USER                 | /* = user-defined */
      CS_ADVECTION_FIELD_TYPE_VELOCITY_VECTOR | /* = define by a vector field */
      CS_ADVECTION_FIELD_DEFINE_AT_VERTICES   | /* = add a field at vertices */
      CS_ADVECTION_FIELD_DEFINE_AT_BOUNDARY_FACES;  /* = add boundary fluxes */

    cs_adv_field_t  *adv = cs_advection_field_add("adv_field", adv_status);

    CS_NO_WARN_IF_UNUSED(adv); /* adv can be used to set options */
  }
  /*! [param_cdo_add_adv_field] */

  /*! [param_cdo_add_user_adv_field_post] */
  {
    /* Retrieve an advection field named "adv_field"  */
    cs_adv_field_t  *adv = cs_advection_field_by_name("adv_field");

    /* Compute the Courant number (if unsteady simulation) */

    cs_advection_field_set_postprocess(adv, CS_ADVECTION_FIELD_POST_COURANT);
  }
  /*! [param_cdo_add_user_adv_field_post] */


  /*! [param_cdo_add_user_tracer] */
  {
    /* Add the user-defined advection field  */

    cs_adv_field_t  *adv_field = cs_advection_field_add_user("adv_field");

    /* Add the user-defined diffusion property */

    cs_property_t  *diff_pty = cs_property_add("diff_pty", CS_PROPERTY_ISO);

    /* Add the user-defined time property */

    cs_property_t  *time_pty = cs_property_add("time_pty", CS_PROPERTY_ISO);

    /* Add a new user-defined equation and associate this equation with some
       properties to get a scalar-valued unsteady convection/diffusion
       equation */

    cs_equation_add_user_tracer("MyTracerEq",  /* Eq. name */
                                "MyTracerVar", /* Variable name */
                                1,             /* Variable dim. */
                                CS_BC_SYMMETRY,
                                time_pty,
                                adv_field,
                                diff_pty);
  }
  /*! [param_cdo_add_user_tracer] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t   *domain)
{
  /*! [param_cdo_domain_output] */
  {
    /* Manage the output logging/checkpoint
     * ====================================
     *   - log frequency in the run_solver.log
     *   - level of verbosity (-1: no, 0, 1, higher is for debugging purpose)
     *   - checkpoint/restart frequency
     */

    cs_domain_set_output_param(domain,
                               CS_RESTART_INTERVAL_ONLY_AT_END,
                               10,  /* log frequency */
                               2);  /* verbosity  */
  }
  /*! [param_cdo_domain_output] */

  /*! [param_cdo_time_step] */
  {
    /* Time step management
     * ====================
     * If there is an inconsistency between the max. number of iteration in
     * time and the final physical time, the first condition encountered stops
     * the calculation.
     */

    cs_domain_set_time_param(domain,
                             100,     /* nt_max or -1 (automatic) */
                             -1.);    /* t_max or < 0. (automatic) */

    /* Define the value of the time step. Two functions are available to do
       this:
       1. cs_domain_def_time_step_by_value(domain, dt_val);
       2. cs_domain_def_time_step_by_func(domain, dt_func);

       The second way enables more complex definitions of the time step.
    */

    cs_domain_def_time_step_by_value(domain, 1.0);
  }
  /*! [param_cdo_time_step] */

  /* ─────────────────────────────────────────────────────────────────────────
   *  Modify the setting of an equation using a generic process
   *  cf. the DOXYGEN documentation (website or in your installation path)
   *
   *         cs_equation_set_param(eqp, CS_EQKEY_*, "key_val")
   *
   * ─────────────────────────────────────────────────────────────────────────*/

  /*! [param_cdo_post_equation] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("AdvDiff.Upw");

    cs_equation_param_set(eqp, CS_EQKEY_VERBOSITY, "2");
    cs_equation_param_set(eqp, CS_EQKEY_EXTRA_OP, "peclet");
  }
  /*! [param_cdo_post_equation] */

  /*! [param_cdo_numerics] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEq");

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");
  }
  /*! [param_cdo_numerics] */

  /*! [param_cdo_conv_numerics] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("AdvDiff.Upw");

    /* Set the advection scheme */

    cs_equation_param_set(eqp, CS_EQKEY_ADV_SCHEME, "upwind");

    /* Set the advection formulation
       - "u.grad(Y)" for a "non-conservative" or gradient formulation
       - "div(u.Y)" for a "conservtive" or divergence formulation
    */

    cs_equation_param_set(eqp, CS_EQKEY_ADV_FORMULATION, "non_conservative");
  }
  /*! [param_cdo_conv_numerics] */

  /*! [param_cdo_conv_schemes] */
  {
    cs_equation_param_t  *eqp1 = cs_equation_param_by_name("Eq1");

    /* Set the advection scheme */

    cs_equation_param_set(eqp1, CS_EQKEY_ADV_SCHEME, "mix_centered_upwind");

    /* Set the portion of upwinding to add to a centered scheme. The same
       portion is added in the whole domain (contrary to the "sg" or
       "samarskii" advection scheme) */

    cs_equation_param_set(eqp1, CS_EQKEY_ADV_UPWIND_PORTION, "0.45");

    /* It's possible to set automatically the upwinding portion using a
       Scharfetter-Gummel scheme or the Samarskii scheme. These two advection
       schemes differ on the weighting function used to compute the portion of
       upwinding. In both cases, one uses an estimation of the local Péclet
       number to evaluate the needed portion of upwinding. */

    cs_equation_param_t  *eqp_sg = cs_equation_param_by_name("Eq2_SG");

    cs_equation_param_set(eqp_sg, CS_EQKEY_ADV_SCHEME, "sg");

    cs_equation_param_t  *eqp_sa = cs_equation_param_by_name("Eq2_SA");

    cs_equation_param_set(eqp_sa, CS_EQKEY_ADV_SCHEME, "samarskii");

    /* In case of a CDO-VCb schemes, one sets a CIP scheme and then modify the
       scaling coefficient in front of the stabilization term */

    cs_equation_param_t  *eqp_vcb = cs_equation_param_by_name("Eq3_CIP");

    cs_equation_param_set(eqp_vcb, CS_EQKEY_ADV_SCHEME, "cip");

    cs_equation_param_set(eqp_vcb, CS_EQKEY_ADV_CIP_COEF, "1");
  }
  /*! [param_cdo_conv_schemes] */

  /*! [param_cdo_diff_numerics] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEq");

    /* Modify other parameters than the space discretization */

    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "cost");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_COEF, "dga");
  }
  /*! [param_cdo_diff_numerics] */

  /*! [param_cdo_time_schemes] */
  {
    cs_equation_param_t  *eqp1 = cs_equation_param_by_name("Eq1");

    /* Set the time scheme */

    cs_equation_param_set(eqp1, CS_EQKEY_TIME_SCHEME, "theta_scheme");

    /* When a theta-scheme is used. It is possible to specify the value of the
       theta weighting (value should be between 0 and 1 */

    cs_equation_param_set(eqp1, CS_EQKEY_TIME_THETA, "0.57");
  }
  /*! [param_cdo_time_schemes] */

  /*! [param_cdo_time_hodge] */
  {
    cs_equation_param_t  *eqp_fe = cs_equation_param_by_name("Eq_FE_like");

    /* Set the algorithm to build the Hodge operator related to the unsteady
       term. This is the same as defining the mass matrix in Finite Element
       schemes */

    cs_equation_param_set(eqp_fe, CS_EQKEY_HODGE_TIME_ALGO, "wbs");

    /* The other choice */

    cs_equation_param_t  *eqp_fv = cs_equation_param_by_name("Eq_FV_like");

    cs_equation_param_set(eqp_fv, CS_EQKEY_HODGE_TIME_ALGO, "voronoi");
  }
  /*! [param_cdo_time_hodge] */

  /*! [param_cdo_sles_settings1] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEq");

    /* Linear algebra settings */

#if defined(HAVE_PETSC)
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_FAMILY, "petsc");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "fcg");
    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
#else
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_FAMILY, "cs");
    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "jacobi");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "cg");
#endif

    /* Manage the stopping criteria */

    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_MAX_ITER, "2500");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_RTOL, "1e-12");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_ATOL, "1e-15");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_DTOL, "1e2");
  }

  /*! [param_cdo_sles_settings1] */

  /*! [param_cdo_sles_amg_type] */
  {
    /* First example with a K-cycle settings */
    /* ------------------------------------- */

    cs_equation_param_t  *eqp0 = cs_equation_param_by_name("Eq_Kcycle");

    cs_equation_param_set(eqp0, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(eqp0, CS_EQKEY_AMG_TYPE, "k_cycle");

    /* Second example with boomer AMG from the HYPRE's library */
    /* ------------------------------------------------------- */

    cs_equation_param_t  *eqp1 = cs_equation_param_by_name("Eq_Boomer");

    cs_equation_param_set(eqp1, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(eqp1, CS_EQKEY_AMG_TYPE, "bamg");

    /* If one wants to access BoomerAMG through the PETSc library, add this
       line. Warning: PETSc should be installed with the HYPRE library */

    cs_equation_param_set(eqp1, CS_EQKEY_SOLVER_FAMILY, "petsc");
  }
  /*! [param_cdo_sles_amg_type] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify the elements such as properties, advection fields,
 *         user-defined equations and modules which have been previously
 *         added.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);

  /* =======================
     User-defined properties
     ======================= */

  /*! [param_cdo_setup_property] */
  {
    cs_property_t  *conductivity = cs_property_by_name("conductivity");
    cs_real_33_t  tensor = {{1.0,  0.5, 0.0}, {0.5, 1.0, 0.5}, {0.0, 0.5, 1.0}};

    cs_property_def_aniso_by_value(conductivity, // property structure
                                   "cells",      // name of the volume zone
                                   tensor);      // values of the property

    cs_property_t  *rhocp = cs_property_by_name("rho.cp");
    cs_real_t  iso_val = 2.0;

    cs_property_def_iso_by_value(rhocp,    // property structure
                                 "cells",  // name of the volume zone
                                 iso_val); // value of the property

  }
  /*! [param_cdo_setup_property] */

  /* =============================
     User-defined advection fields
     ============================= */

  /*! [param_cdo_setup_advfield] */
  {
    cs_adv_field_t  *adv = cs_advection_field_by_name("adv_field");

    cs_advection_field_def_by_analytic(adv, _define_adv_field, NULL);
  }
  /*! [param_cdo_setup_advfield] */

  /* ======================
     User-defined equations
     ====================== */

  /* Define the boundary conditions
     >> cs_equation_add_bc_by_analytic(eqp,
                                       bc_type,
                                       "zone_name",
                                       analytic_function);

     -> eq is the structure related to the equation to set
     -> type of boundary condition:
        CS_BC_HMG_DIRICHLET, CS_BC_DIRICHLET, CS_BC_SYMMETRY,
        CS_BC_NEUMANN, CS_BC_NEUMANN_FULL, CS_BC_ROBIN,
        CS_BC_CIRCULATION, CS_BC_WALL_MODELLED:

     >> cs_equation_add_bc_by_value(eqp,
                                    bc_type,
                                    "mesh_location_name",
                                    values); // pointer
  */

  /*! [param_cdo_setup_bcs] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("AdvDiff.Upw");

    cs_equation_add_bc_by_analytic(eqp,
                                   CS_BC_DIRICHLET,
                                   "boundary_faces",  // zone name
                                   _define_bcs,       // pointer to the function
                                   NULL);             // input structure
  }
  /*! [param_cdo_setup_bcs] */

  /*! [param_cdo_add_terms] */
  {
    /* Retrieve the set of equation parameters of the equation to modify */

    cs_equation_param_t  *eqp = cs_equation_param_by_name("AdvDiff.Upw");

    /* Add an unsteady term (draw a link between this equation and a
       property) */

    cs_property_t  *rhocp = cs_property_by_name("rho.cp");

    cs_equation_add_time(eqp, rhocp);

    /* Add a diffusion term (draw a link between this equation and a
       property) */

    cs_property_t  *conductivity = cs_property_by_name("conductivity");

    cs_equation_add_diffusion(eqp, conductivity);

    /* Add an advection term (draw a link between this equation and an
       advection field) */

    cs_adv_field_t  *adv = cs_advection_field_by_name("adv_field");

    cs_equation_add_advection(eqp, adv);
  }
  /*! [param_cdo_add_terms] */

    /*! [param_cdo_add_terms_2] */
  {
    /* Retrieve the set of equation parameters of the equation to modify */

    cs_equation_param_t  *eqp = cs_equation_param_by_name("MyEq");

    /* Add a curl-curl term (draw a link between this equation and a
       property) */

    cs_property_t  *mu = cs_property_by_name("mu");

    cs_equation_add_curlcurl(eqp,
                             mu,
                             false); /* reciprocal of the property */

    /* Add a grad-div term (draw a link between this equation and a
       property) */

    cs_property_t  *gamma = cs_property_by_name("gamma");

    cs_equation_add_graddiv(eqp, gamma);
  }
  /*! [param_cdo_add_terms_2] */

  /*! [param_cdo_add_simple_source_terms] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("AdvDiff.Upw");

    /* Simple definition with cs_equation_add_source_term_by_val
       where the value of the source term is given by m^3
    */

    cs_real_t  st_val = -0.1;
    cs_equation_add_source_term_by_val(eqp, "cells", &st_val);
  }
  /*! [param_cdo_add_simple_source_terms] */

  /*! [param_cdo_add_source_terms] */
  {
    /* More advanced definition relying on an analytic function and a context
       structure */

    cs_equation_param_t  *eqp2 = cs_equation_param_by_name("MyEq");

    double  *input = NULL;
    BFT_MALLOC(input, 2, double);

    input[0] = 0.5;            /* Value of the reaction coefficient */
    input[1] = 4.0*atan(1.0);  /* Value of pi (computed only once) */

    cs_xdef_t  *st =
      cs_equation_add_source_term_by_analytic(eqp2,
                                              "x_leq_0", /* zone name */
                                              _my_source_term,
                                              input);    /* context structure */

    /* Specify how the cs_xdef_t structure can free its input structure */

    cs_xdef_set_free_input_function(st, _free_input);

    /* st can be used for advanced settings like quadrature rules. By default,
       a barycentric quadrature rule is used (one evaluation at the cell center
       by mesh cell) */

    cs_xdef_set_quadrature(st, CS_QUADRATURE_HIGHER);

    /* Each cell is subdivided into an implicit tetrahedral sub-mesh relying on
       the vertices, the middle of each edge, the face barycenter and the cell
       center. A quadrature exact up to 3rd order polynomial is set. */
  }
  /*! [param_cdo_add_source_terms] */

  /*! [param_cdo_copy_settings] */
  {
    /* Copy the settings from AdvDiff.Upw to AdvDiff.SG */
    cs_equation_param_t  *eqp_ref = cs_equation_param_by_name("AdvDiff.Upw");

    /* Another example to retrieve the cs_equation_param structure */
    cs_equation_t  *eq = cs_equation_by_name("AdvDiff.SG");
    cs_equation_param_t  *eqp = cs_equation_get_param(eq);

    /* Copy the settings to start from a common state as in the reference
     * equation (in this case, do not copy the associated field id since these
     * are two different equations with their own variable field) */
    bool  copy_field_id = false;
    cs_equation_copy_param_from(eqp_ref, eqp, copy_field_id);

    /* Keep all the settings from "AdvDiff.Upw and then only change the
       advection scheme for the second equation */
    cs_equation_param_set(eqp, CS_EQKEY_ADV_SCHEME, "sg");

    /* Call this function to be sure that the linear solver is set to what
       one wants (if there is no modification of the SLES parameters, this
       step can be skipped) */
    cs_equation_param_set_sles(eqp);

  }
  /*! [param_cdo_copy_settings] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
