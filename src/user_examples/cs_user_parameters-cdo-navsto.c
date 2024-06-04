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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
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
 * See \subpage parameters for examples.
 */
/*----------------------------------------------------------------------------*/

#define _PI   3.14159265358979323846

static const cs_real_t _2pi  = 2*_PI, _pis  = _PI*_PI;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Velocity at (x, y, z)
 *    u =  0.5 * sin(2pi x) cos(2pi y) cos(2pi z)
 *    v =  0.5 * cos(2pi x) sin(2pi y) cos(2pi z)
 *    w = -      cos(2pi x) cos(2pi y) sin(2pi z)
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_navsto_vel_function] */
static inline void
__vel(const cs_real_3_t pxyz,
      cs_real_3_t       res)
{
  const cs_real_t x  = pxyz[0], y = pxyz[1], z = pxyz[2];
  const cs_real_t sx = sin(_2pi*x), cx = cos(_2pi*x),
                  sy = sin(_2pi*y), cy = cos(_2pi*y),
                  sz = sin(_2pi*z), cz = cos(_2pi*z);

  res[0] = .5 * sx * cy * cz;
  res[1] = .5 * cx * sy * cz;
  res[2] = -    cx * cy * sz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function used to define the velocity value everywhere. This prototype
 *         follows the one defined for all analytic functions.
 *         elt_ids is optional. If not NULL, it enables to get an access in
 *         coords at the right location and the same thing to fill retval if
 *         compact is set to false
 *
 * \param[in]      time     when ?
 * \param[in]      n_pts    number of elements to consider
 * \param[in]      pt_ids   list of elements ids (to access coords and fill)
 * \param[in]      coords   where ?
 * \param[in]      compact  true:no indirection, false:indirection for filling
 * \param[in]      input    pointer to a structure cast on-the-fly (may be NULL)
 * \param[in, out] res      result of the function
 */
/*----------------------------------------------------------------------------*/

static void
_vel_def(cs_real_t           time,
         cs_lnum_t           n_pts,
         const cs_lnum_t    *pt_ids,
         const cs_real_t    *xyz,
         bool                dense_output,
         void               *input,
         cs_real_t          *res)
{
  CS_UNUSED(input);
  CS_UNUSED(time);

  if (pt_ids != NULL && !dense_output) {

#   pragma omp parallel for if(n_pts > CS_THR_MIN)
    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_lnum_t id = pt_ids[p];
      __vel(xyz + 3*id, res + 3*id);
    }

  }
  else if (pt_ids != NULL && dense_output) {

#   pragma omp parallel for if(n_pts > CS_THR_MIN)
    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_lnum_t id = pt_ids[p];
      __vel(xyz + 3*id, res + 3*p);
    }

  }
  else {

#   pragma omp parallel for if(n_pts > CS_THR_MIN)
    for (cs_lnum_t p = 0; p < n_pts; p++)
      __vel(xyz + 3*p, res + 3*p);

  }
}
/*! [param_cdo_navsto_vel_function] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief vector-valued source term definition at (x, y, z)
 *
 * sx =  6 sin(2pi x)cos(2pi y)cos(2pi z) + 2pi cos(2pi x)sin(2pi y)sin(2pi z)
 * sy =  6 cos(2pi x)sin(2pi y)cos(2pi z) + 2pi sin(2pi x)cos(2pi y)sin(2pi z)
 * sz =-12 cos(2pi x)cos(2pi y)sin(2pi z) + 2pi sin(2pi x)sin(2pi y)cos(2pi z)
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_navsto_st_function] */
static inline void
__st(const cs_real_3_t pxyz,
     cs_real_3_t       res)
{
  const cs_real_t x  = pxyz[0], y = pxyz[1], z = pxyz[2];
  const cs_real_t sx = sin(_2pi*x), cx = cos(_2pi*x),
                  sy = sin(_2pi*y), cy = cos(_2pi*y),
                  sz = sin(_2pi*z), cz = cos(_2pi*z);
  res[0] =   6.*_pis * sx * cy * cz + _2pi * cx * sy * sz;
  res[1] =   6.*_pis * cx * sy * cz + _2pi * sx * cy * sz;
  res[2] = -12.*_pis * cx * cy * sz + _2pi * sx * sy * cz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the source term definition (steady case)
 *        This prototype follows the one defined for all analytic functions.
 *        pt_ids is optional. If not NULL, it enables to access the xyz array
 *        at the right location and the same thing to fill the res array if
 *        compact is set to false
 *
 * \param[in]      time     when ?
 * \param[in]      n_pts    number of elements to consider
 * \param[in]      pt_ids   list of elements ids (to access coords and fill)
 * \param[in]      xyz      where ?
 * \param[in]      compact  true:no indirection, false:indirection for filling
 * \param[in]      input    pointer to a structure cast on-the-fly (may be NULL)
 * \param[in, out] res      result of the function
 */
/*----------------------------------------------------------------------------*/

static void
_src_def(cs_real_t           time,
         cs_lnum_t           n_pts,
         const cs_lnum_t    *pt_ids,
         const cs_real_t    *xyz,
         bool                dense_output,
         void               *input,
         cs_real_t          *res)
{
  CS_UNUSED(input);
  CS_UNUSED(time);

  if (pt_ids != NULL && !dense_output) {

#   pragma omp parallel for if(n_pts > CS_THR_MIN)
    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_lnum_t id = pt_ids[p];
      __st(xyz + 3*id, res + 3*id);
    }

  }
  else if (pt_ids != NULL && dense_output) {

#   pragma omp parallel for if(n_pts > CS_THR_MIN)
    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_lnum_t id = pt_ids[p];
      __st(xyz + 3*id, res + 3*p);
    }

  }
  else {

#   pragma omp parallel for if(n_pts > CS_THR_MIN)
    for (cs_lnum_t p = 0; p < n_pts; p++)
      __st(xyz + 3*p, res + 3*p);

  }
}
/*! [param_cdo_navsto_st_function] */

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

  /*! [param_cdo_navsto_boundary] */
  {
    cs_domain_t  *domain = cs_glob_domain;

    cs_boundary_set_default(domain->boundaries, CS_BOUNDARY_WALL);

    /* Define an "inlet" boundary */

    cs_boundary_add(domain->boundaries,
                    CS_BOUNDARY_INLET | CS_BOUNDARY_IMPOSED_VEL,
                    "inlet_faces"); /* boundary zone name */

    /* Define an "outlet" boundary */

    cs_boundary_add(domain->boundaries,
                    CS_BOUNDARY_OUTLET,
                    "outlet_faces"); /* boundary zone name */
  }
  /*! [param_cdo_navsto_boundary] */

  /*! [param_cdo_navsto_activate] */
  {
    cs_domain_t  *domain = cs_glob_domain;

    cs_navsto_system_activate(domain->boundaries,
                              /* Model: type + additional flag */
                              CS_NAVSTO_MODEL_STOKES,
                              CS_NAVSTO_MODEL_STEADY,
                              /* Velocity-pressure coupling */
                              CS_NAVSTO_COUPLING_MONOLITHIC,
                              /* Post flag */
                              CS_NAVSTO_POST_VELOCITY_DIVERGENCE   |
                              0);
  }
  /*! [param_cdo_navsto_activate] */

  /*! [param_cdo_navsto_activate2] */
  {
    cs_domain_t  *domain = cs_glob_domain;
    cs_flag_t  post_flag =
      CS_NAVSTO_POST_VELOCITY_DIVERGENCE |
      CS_NAVSTO_POST_VELOCITY_GRADIENT   |
      CS_NAVSTO_POST_MASS_DENSITY        |
      CS_NAVSTO_POST_PRESSURE_GRADIENT;

    cs_navsto_system_activate(domain->boundaries,
                              /* Model: type + additional flag */
                              CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES,
                              0,
                              /* Velocity-pressure coupling */
                              CS_NAVSTO_COUPLING_MONOLITHIC,
                              /* Post flag */
                              post_flag);
  }
  /*! [param_cdo_navsto_activate2] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t    *domain)
{
  /*! [param_cdo_navsto_sles_alu] */
  {
    /* Parameters related to the momentum equation */

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings for the saddle-point system */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "alu");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_AUGMENT_SCALING, "100");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_MAX_ITER, "1000");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_RTOL, "1e-8");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_ATOL, "1e-14");

#if defined(HAVE_MUMPS)
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "mumps");
#else
    bft_error(__FILE__, __LINE__, 0, "%s: MUMPS is not available\n", __func__);
#endif
  }
  /*! [param_cdo_navsto_sles_alu] */

  /*! [param_cdo_navsto_sles_mumps] */
  {
    /* Parameters related to the momentum equation */

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");
    cs_param_sles_t  *slesp = cs_equation_param_get_sles_param(mom_eqp);

    cs_param_sles_mumps(slesp,
                        false,                       /* single-precision ? */
                        CS_PARAM_MUMPS_FACTO_LDLT_SPD);  /* type of facto. */

    cs_param_sles_mumps_advanced(slesp,
                                 CS_PARAM_MUMPS_ANALYSIS_PTSCOTCH,
                                 3,  /* size of the block for analysis */
                                 -1, /* pct memory increase < 0 --> not used */
                                 0,  /* BLR compression: 0 --> not used */
                                 0,  /* iterative refinement steps */
                                 CS_PARAM_MUMPS_MEMORY_CONSTRAINED,
                                 true); /* advanced optimizations */
  }
  /*! [param_cdo_navsto_sles_mumps] */

  /*! [param_cdo_navsto_sles_notay] */
  {
    /* Parameters related to the momentum equation */

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    /* Linear algebra settings for the saddle-point system.
     *  Notay's algebraic transformation is an experimental feature
     *  The main solver is a FGMRES on the full saddle-point problem.
     */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_SOLVER, "notay");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_MAX_ITER, "1000");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_RTOL, "1e-8");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SADDLE_ATOL, "1e-14");

    /* Set the solver for the (1,1)-block preconditioner --> the velocity block
       in the case of the Navier-Stokes system */

    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER, "fgmres");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(mom_eqp, CS_EQKEY_PRECOND_BLOCK_TYPE, "upper");
    cs_equation_param_set(mom_eqp, CS_EQKEY_AMG_TYPE, "boomer");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_RTOL, "1e-2");
    cs_equation_param_set(mom_eqp, CS_EQKEY_SOLVER_MAX_ITER, "100");

  }
  /*! [param_cdo_navsto_sles_notay] */

  /*! [param_cdo_navsto_schur_mass_scaled] */
  {
    /* Parameters related to the momentum equation */

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    cs_equation_param_set(mom_eqp,
                          CS_EQKEY_SADDLE_SCHUR_APPROX, "mass_scaled");

    /* Nothing else to add for the approximation of the Schur complement */
  }
  /*! [param_cdo_navsto_schur_mass_scaled] */

  /*! [param_cdo_navsto_schur_mass_scaled_diag_inv] */
  {
    /* Parameters related to the momentum equation */

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    cs_equation_param_set(mom_eqp,
                          CS_EQKEY_SADDLE_SCHUR_APPROX, "mass_scaled_diag_inv");

    cs_param_saddle_t  *saddlep = cs_equation_param_get_saddle_param(mom_eqp);

    /* Retrieve the set of parameters to handle the system associated to the
       Schur complement approximation */

    cs_param_sles_t  *schur_slesp =
      cs_param_saddle_get_schur_sles_param(saddlep);

    /* Set the solver, its preconditionner and the stopping criteria */

    int ierr = cs_param_sles_set_solver("fgmres", schur_slesp);
    assert(ierr == 0);

    ierr = cs_param_sles_set_precond("amg", schur_slesp);
    assert(ierr == 0);

    ierr = cs_param_sles_set_amg_type("boomer", schur_slesp);
    assert(ierr == 0);

    /* One uses CS_CDO_KEEP_DEFAULT as parameter to keep unchanged the default
       settings */

    cs_param_sles_set_cvg_param(schur_slesp,
                                1e-2,                 /* rtol */
                                CS_CDO_KEEP_DEFAULT,  /* atol */
                                CS_CDO_KEEP_DEFAULT,  /* dtol */
                                CS_CDO_KEEP_DEFAULT); /* max. iter. */
  }
  /*! [param_cdo_navsto_schur_mass_scaled_diag_inv] */

  /*! [param_cdo_navsto_schur_lumped_inv] */
  {
    /* Parameters related to the momentum equation */

    cs_equation_param_t  *mom_eqp = cs_equation_param_by_name("momentum");

    cs_equation_param_set(mom_eqp,
                          CS_EQKEY_SADDLE_SCHUR_APPROX, "lumped_inv");

    cs_param_saddle_t  *saddlep = cs_equation_param_get_saddle_param(mom_eqp);

    /* Retrieve the set of parameters to handle the system associated to the
       Schur complement approximation */

    cs_param_sles_t  *schur_slesp =
      cs_param_saddle_get_schur_sles_param(saddlep);

    /* Set the solver, its preconditionner and the stopping criteria */

    int ierr = cs_param_sles_set_solver("fcg", schur_slesp);
    assert(ierr == 0);

    ierr = cs_param_sles_set_precond("amg", schur_slesp);
    assert(ierr == 0);

    /* One uses CS_CDO_KEEP_DEFAULT as parameter to keep unchanged the default
       settings */

    cs_param_sles_set_cvg_param(schur_slesp,
                                1e-2,                 /* rtol */
                                CS_CDO_KEEP_DEFAULT,  /* atol */
                                CS_CDO_KEEP_DEFAULT,  /* dtol */
                                15);                  /* max. iter. */


    ierr = cs_param_sles_set_amg_type("gamg", schur_slesp);
    assert(ierr == 0);

    /* Retrieve the additional system (the one used for the lumped inverse") */

    cs_param_sles_t  *xtra_slesp = cs_param_saddle_get_xtra_sles_param(saddlep);

    ierr = cs_param_sles_set_solver("fgmres", xtra_slesp);
    assert(ierr == 0);

    ierr = cs_param_sles_set_precond("amg", xtra_slesp);
    assert(ierr == 0);

    ierr = cs_param_sles_set_amg_type("boomer", xtra_slesp);
    assert(ierr == 0);

    ierr = cs_param_sles_set_precond_block_type("diag", xtra_slesp);
    assert(ierr == 0);

    cs_param_sles_set_cvg_param(xtra_slesp,
                                1e-2,                 /* rtol */
                                CS_CDO_KEEP_DEFAULT,  /* atol */
                                CS_CDO_KEEP_DEFAULT,  /* dtol */
                                20);                  /* max. iter. */
  }
  /*! [param_cdo_navsto_schur_lumped_inv] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  - Specify the elements such as properties, advection fields,
 *           user-defined equations and modules which have been previously
 *           added.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);

  /*! [param_cdo_navsto_bc1] */
  {
    cs_navsto_param_t  *nsp = cs_navsto_system_get_param();

    /* Define the boundary conditions from a function called "_vel_def" */

    cs_navsto_set_velocity_inlet_by_analytic(nsp,
                                             "boundary_faces",
                                             _vel_def,
                                             NULL);
  }
  /*! [param_cdo_navsto_bc1] */

  /*! [param_cdo_navsto_st1] */
  {
    cs_navsto_param_t  *nsp = cs_navsto_system_get_param();

    /* Define the source term */

    cs_xdef_t *st_def =
      cs_navsto_add_source_term_by_analytic(nsp, "cells", _src_def, NULL);

    /* Optionally, the parameters of this definition can be updated */

    cs_xdef_set_quadrature(st_def, CS_QUADRATURE_BARY_SUBDIV);

  }
  /*! [param_cdo_navsto_st1] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
