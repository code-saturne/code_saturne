/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_advection_field.h"
#include "cs_base.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_param.h"
#include "cs_property.h"
#include "cs_prototypes.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_walldistance.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-base.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \subpage parameters for examples.
 */
/*----------------------------------------------------------------------------*/

static const double  one_third = 1./3.;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Give the explicit definition of the advection field.
 *         pt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *         Rely on a generic function pointer for an analytic function
 *
 * \param[in]      time       when ?
 * \param[in]      n_elts     number of elements to consider
 * \param[in]      pt_ids     list of elements ids (to access coords and fill)
 * \param[in]      coords     where ?
 * \param[in]      compact    true:no indirection, false:indirection for filling
 * \param[in]      input      NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval     result of the function
 */
/*----------------------------------------------------------------------------*/

static void
_define_adv_field(cs_real_t           time,
                  cs_lnum_t           n_pts,
                  const cs_lnum_t    *pt_ids,
                  const cs_real_t    *xyz,
                  bool                compact,
                  void               *input,
                  cs_real_t          *res)
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  if (pt_ids != NULL && !compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {

      const cs_lnum_t  id = pt_ids[p];
      const cs_real_t  *pxyz = xyz + 3*id;
      cs_real_t  *pres = res + 3*id;

      pres[0] = pxyz[1] - 0.5;
      pres[1] = 0.5 - pxyz[0];
      pres[2] = pxyz[2];

    }

  }
  else if (pt_ids != NULL && compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {

      const cs_real_t  *pxyz = xyz + 3*pt_ids[p];
      cs_real_t  *pres = res + 3*p;

      pres[0] = pxyz[1] - 0.5;
      pres[1] = 0.5 - pxyz[0];
      pres[2] = pxyz[2];

    }

  }
  else {

    assert(pt_ids == NULL);
    for (cs_lnum_t p = 0; p < n_pts; p++) {

      const cs_real_t  *pxyz = xyz + 3*p;
      cs_real_t  *pres = res + 3*p;

      pres[0] = pxyz[1] - 0.5;
      pres[1] = 0.5 - pxyz[0];
      pres[2] = pxyz[2];

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Give the explicit definition of the dirichlet boundary conditions
 *         pt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *         Rely on a generic function pointer for an analytic function
 *
 * \param[in]      time      when ?
 * \param[in]      n_elts    number of elements to consider
 * \param[in]      pt_ids    list of elements ids (to access coords and fill)
 * \param[in]      coords    where ?
 * \param[in]      compact   true:no indirection, false:indirection for filling
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 * \param[in, out] res       result of the function
 */
/*----------------------------------------------------------------------------*/

static void
_define_bcs(cs_real_t           time,
            cs_lnum_t           n_pts,
            const cs_lnum_t    *pt_ids,
            const cs_real_t    *xyz,
            bool                compact,
            void               *input,
            cs_real_t          *res)
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  const double  pi = 4.0*atan(1.0);
  if (pt_ids != NULL && !compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {

      const cs_lnum_t  id = pt_ids[p];
      const cs_real_t  *_xyz = xyz + 3*id;
      const double  x = _xyz[0], y = _xyz[1], z = _xyz[2];

      res[id] = 1 + sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+cs_math_onethird));

    }

  }
  else if (pt_ids != NULL && compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_real_t  *_xyz = xyz + 3*pt_ids[p];
      const double  x = _xyz[0], y = _xyz[1], z = _xyz[2];

      res[p] = 1 + sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+cs_math_onethird));
    }

  }
  else {

    assert(pt_ids == NULL);
    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_real_t  *_xyz = xyz + 3*p;
      const double  x = _xyz[0], y = _xyz[1], z = _xyz[2];

      res[p] = 1 + sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+cs_math_onethird));
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Give the explicit definition of the source term
 *         pt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *         Rely on a generic function pointer for an analytic function
 *
 * \param[in]      time      when ?
 * \param[in]      n_elts    number of elements to consider
 * \param[in]      pt_ids    list of elements ids (to access coords and fill)
 * \param[in]      coords    where ?
 * \param[in]      compact   true:no indirection, false:indirection for filling
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 * \param[in, out] res       result of the function
 */
/*----------------------------------------------------------------------------*/

static void
_define_source(cs_real_t           time,
               cs_lnum_t           n_pts,
               const cs_lnum_t     pt_ids[],
               const cs_real_t    *xyz,
               bool                compact,
               void               *input,
               cs_real_t          *res)
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  const double  pi = 4.0*atan(1.0), pi2 = pi*pi;
  if (pt_ids != NULL && !compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {

      const cs_lnum_t  id = pt_ids[p];
      const double  x = xyz[3*id], y = xyz[3*id+1], z = xyz[3*id+2];
      const double  cpx = cos(pi*x), spx = sin(pi*x);
      const double  cpy = cos(pi*(y+0.5)), spy = sin(pi*(y+0.5));
      const double  cpz = cos(pi*(z+one_third)), spz = sin(pi*(z+one_third));

      /* first derivatives */
      cs_real_t gx = pi*cpx*spy*spz;
      cs_real_t gy = pi*spx*cpy*spz;
      cs_real_t gz = pi*spx*spy*cpz;

      /* second derivatives */
      cs_real_t gxx, gyy, gzz, gxy, gxz, gyz;
      gxx = gyy = gzz = -pi2*spx*spy*spz;
      gxy = pi2*cpx*cpy*spz, gxz = pi2*cpx*spy*cpz, gyz = pi2*spx*cpy*cpz;

      /* Material property */
      cs_real_33_t  cond;
      cond[0][0] = 1.0, cond[0][1] = 0.5, cond[0][2] = 0.0;
      cond[1][0] = 0.5, cond[1][1] = 1.0, cond[1][2] = 0.5;
      cond[2][0] = 0.0, cond[2][1] = 0.5, cond[2][2] = 1.0;

      /* Contribution of the diffusive part */
      res[id] = cond[0][0]*gxx + cond[1][1]*gyy + cond[2][2]*gzz +
        2*( cond[0][1]*gxy + cond[0][2]*gxz + cond[1][2]*gyz);
      res[id] *= -1;

      /* Contribution of the advection term */
      res[id] += (y - 0.5)*gx + (0.5 - x)*gy + z*gz + 1 + spx*spy*spz;

    } // Loop on evaluation points

  }
  else if (pt_ids != NULL && compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {

      const cs_lnum_t  id = pt_ids[p];
      const double  x = xyz[3*id], y = xyz[3*id+1], z = xyz[3*id+2];
      const double  cpx = cos(pi*x), spx = sin(pi*x);
      const double  cpy = cos(pi*(y+0.5)), spy = sin(pi*(y+0.5));
      const double  cpz = cos(pi*(z+one_third)), spz = sin(pi*(z+one_third));

      /* first derivatives */
      cs_real_t gx = pi*cpx*spy*spz;
      cs_real_t gy = pi*spx*cpy*spz;
      cs_real_t gz = pi*spx*spy*cpz;

      /* second derivatives */
      cs_real_t gxx, gyy, gzz, gxy, gxz, gyz;
      gxx = gyy = gzz = -pi2*spx*spy*spz;
      gxy = pi2*cpx*cpy*spz, gxz = pi2*cpx*spy*cpz, gyz = pi2*spx*cpy*cpz;

      /* Material property */
      cs_real_33_t  cond;
      cond[0][0] = 1.0, cond[0][1] = 0.5, cond[0][2] = 0.0;
      cond[1][0] = 0.5, cond[1][1] = 1.0, cond[1][2] = 0.5;
      cond[2][0] = 0.0, cond[2][1] = 0.5, cond[2][2] = 1.0;

      /* Contribution of the diffusive part */
      res[p] = cond[0][0]*gxx + cond[1][1]*gyy + cond[2][2]*gzz +
        2*( cond[0][1]*gxy + cond[0][2]*gxz + cond[1][2]*gyz);
      res[p] *= -1;

      /* Contribution of the advection term */
      res[p] += (y - 0.5)*gx + (0.5 - x)*gy + z*gz + 1 + spx*spy*spz;

    } // Loop on evaluation points

  }
  else {

    assert(pt_ids == NULL);
    for (cs_lnum_t p = 0; p < n_pts; p++) {

      const double  x = xyz[3*p], y = xyz[3*p+1], z = xyz[3*p+2];
      const double  cpx = cos(pi*x), spx = sin(pi*x);
      const double  cpy = cos(pi*(y+0.5)), spy = sin(pi*(y+0.5));
      const double  cpz = cos(pi*(z+one_third)), spz = sin(pi*(z+one_third));

      /* first derivatives */
      cs_real_t gx = pi*cpx*spy*spz;
      cs_real_t gy = pi*spx*cpy*spz;
      cs_real_t gz = pi*spx*spy*cpz;

      /* second derivatives */
      cs_real_t gxx, gyy, gzz, gxy, gxz, gyz;
      gxx = gyy = gzz = -pi2*spx*spy*spz;
      gxy = pi2*cpx*cpy*spz, gxz = pi2*cpx*spy*cpz, gyz = pi2*spx*cpy*cpz;

      /* Material property */
      cs_real_33_t  cond;
      cond[0][0] = 1.0, cond[0][1] = 0.5, cond[0][2] = 0.0;
      cond[1][0] = 0.5, cond[1][1] = 1.0, cond[1][2] = 0.5;
      cond[2][0] = 0.0, cond[2][1] = 0.5, cond[2][2] = 1.0;

      /* Contribution of the diffusive part */
      res[p] = cond[0][0]*gxx + cond[1][1]*gyy + cond[2][2]*gzz +
        2*( cond[0][1]*gxy + cond[0][2]*gxz + cond[1][2]*gyz);
      res[p] *= -1;

      /* Contribution of the advection term */
      res[p] += (y - 0.5)*gx + (0.5 - x)*gy + z*gz + 1 + spx*spy*spz;

    } // Loop on evaluation points

  }

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

  /*! [param_cdo_activation] */
  {
    /* Activate CDO/HHO module so that main additional structure are built */

    cs_domain_t  *domain = cs_glob_domain;

    cs_domain_set_cdo_mode(domain, CS_DOMAIN_CDO_MODE_ONLY);
  }
  /*! [param_cdo_activation] */

  /*! [param_cdo_domain_boundary] */
  {
    /* ======================
       Boundary of the domain
       ====================== */

    cs_domain_t  *domain = cs_glob_domain;

    /* Choose a boundary by default */
    cs_domain_set_default_boundary(domain, CS_DOMAIN_BOUNDARY_WALL);

    /* Add a new boundary
       >> cs_domain_add_boundary(domain, type_of_boundary, zone_name);

       * zone_name is either a predefined one or user-defined one
       * type_of_boundary is one of the following keywords:
         CS_DOMAIN_BOUNDARY_WALL,
         CS_DOMAIN_BOUNDARY_INLET,
         CS_DOMAIN_BOUNDARY_OUTLET,
         CS_DOMAIN_BOUNDARY_SYMMETRY
    */

    cs_domain_add_boundary(domain, CS_DOMAIN_BOUNDARY_INLET, "in");
    cs_domain_add_boundary(domain, CS_DOMAIN_BOUNDARY_OUTLET, "out");

  }
  /*! [param_cdo_domain_boundary] */


  /*! [param_cdo_domain_ouput] */
  {
    cs_domain_t  *domain = cs_glob_domain;

    /* =========================
       Generic output management
       ========================= */

    cs_domain_set_output_param(domain,
                               10,     // output log frequency
                               2);     // verbosity (-1: no, 0, ...)

  }
  /*! [param_cdo_domain_ouput] */

  /*! [param_cdo_time_step] */
  {
    cs_domain_t  *domain = cs_glob_domain;

    /* ====================
       Time step management
       ====================

       If there is an inconsistency between the max. number of iteration in
       time and the final physical time, the first condition encountered stops
       the calculation.
    */

    cs_domain_set_time_param(domain,
                             100,     // nt_max or -1 (automatic)
                             -1.);    // t_max or < 0. (automatic)

    /* Define the value of the time step
       >> cs_domain_def_time_step_by_value(domain, dt_val);
       >> cs_domain_def_time_step_by_func(domain, dt_func);

       The second way to define the time step enable complex definitions.
       dt_func must have the following prototype:

       double dt_func(int  nt_cur, double  time_cur)
    */

    cs_domain_def_time_step_by_value(domain, 1.0);

  }
  /*! [param_cdo_time_step] */

  /*! [param_cdo_wall_distance] */
  {
    cs_walldistance_activate();
  }
  /*! [param_cdo_wall_distance] */

  /*! [param_cdo_add_user_equation] */
  {
    /* Add user equation:
       -- default boundary condition is among:
       CS_PARAM_BC_HMG_DIRICHLET or CS_PARAM_BC_HMG_NEUMANN

       By default, initial values are set to zero (or the value given by the
       restart file in case of restart).
    */

    cs_equation_add_user("AdvDiff",     // equation name
                         "Potential",   // associated variable field name
                         1,             // dimension of the unknown
                         CS_PARAM_BC_HMG_DIRICHLET); // default boundary
  }
  /*! [param_cdo_add_user_equation] */

  /*! [param_cdo_add_user_properties] */
  {

    /* ========================================
       Add material properties/advection fields
       ======================================== */

    cs_property_add("conductivity",      // property name
                    CS_PROPERTY_ANISO);  // type of material property
    cs_property_add("rho.cp",            // property name
                    CS_PROPERTY_ISO);    // type of material property

    cs_advection_field_add("adv_field"); // name of the new advection field

  }
  /*! [param_cdo_add_user_properties] */

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
cs_user_parameters(void)
{
  /*! [param_cdo_numerics] */
  {

    /* Modify the setting of an equation
       =================================

       Available keywords and related values for keywords are described in
       the DOXYGEN documentation.
    */

    cs_equation_param_t  *eqp = cs_equation_param_by_name("FVCA6.1");

    /* The modification of the space discretization should be apply first */
    cs_equation_set_param(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vb");

    /* Modify other parameters than the space discretization */
    cs_equation_set_param(eqp, CS_EQKEY_VERBOSITY, "2");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "cost");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_COEF, "dga");

    /* Linear algebra settings */
#if defined(HAVE_PETSC)
    cs_equation_set_param(eqp, CS_EQKEY_SOLVER_FAMILY, "petsc");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "amg");
#else
    cs_equation_set_param(eqp, CS_EQKEY_SOLVER_FAMILY, "cs");
    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "jacobi");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
#endif
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL_MAX_ITER, "2500");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL_EPS, "1e-12");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL_RESNORM, "false");

  }
  /*! [param_cdo_numerics] */
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
cs_user_cdo_finalize_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

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

    /* Enable also the defintion of the advection field at mesh vertices */
    cs_advection_field_set_option(adv, CS_ADVKEY_DEFINE_AT_VERTICES);

    /* Activate the post-processing of the related Courant number */
    cs_advection_field_set_option(adv, CS_ADVKEY_POST_COURANT);
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
        CS_PARAM_BC_DIRICHLET, CS_PARAM_BC_HMG_DIRICHLET,
        CS_PARAM_BC_NEUMANN, CS_PARAM_BC_HMG_NEUMANN, CS_PARAM_BC_ROBIN

     >> cs_equation_add_bc_by_value(eqp,
                                    bc_type,
                                    "mesh_location_name",
                                    values); // pointer
  */

  /*! [param_cdo_setup_bcs] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("AdvDiff");

    cs_equation_add_bc_by_analytic(eqp,
                                   CS_PARAM_BC_DIRICHLET,
                                   "boundary_faces",  // zone name
                                   _define_bcs,       // pointer to the function
                                   NULL);             // input structure

  }
  /*! [param_cdo_setup_bcs] */

  /*! [param_cdo_add_terms] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("AdvDiff");
    cs_property_t  *rhocp = cs_property_by_name("rho.cp");
    cs_property_t  *conductivity = cs_property_by_name("conductivity");
    cs_adv_field_t  *adv = cs_advection_field_by_name("adv_field");

    /* Activate the unsteady term */
    cs_equation_add_time(eqp, rhocp);

    /* Activate the diffusion term */
    cs_equation_add_diffusion(eqp, conductivity);

    /* Activate advection effect */
    cs_equation_add_advection(eqp, adv);

    /* Add a source term: _define_source is a user-defined function with a
       a specific list of arguments.
       Simpler definition can be used: cs_equation_add_source_term_by_val
       where the value of the source term is given by m^3
    */
    cs_xdef_t  *st = cs_equation_add_source_term_by_analytic(eqp,
                                                             "cells",
                                                             _define_source,
                                                             NULL);

    /* Optional: specify the quadrature used for computing the source term */
    cs_xdef_set_quadrature(st, CS_QUADRATURE_BARY);
  }
  /*! [param_cdo_add_terms] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
