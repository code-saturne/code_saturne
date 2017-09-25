/*============================================================================
 * Set main parameters for the current simulation when the CDO kernel is used
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

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_boundary_zone.h"
#include "cs_mesh_location.h"
#include "cs_sdm.h"
#include "cs_property.h"
#include "cs_advection_field.h"
#include "cs_walldistance.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_cdo-condif.c
 *
 * \brief Main user function for setting up a calculation with CDO.
 *
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! \endcond (end ignore by Doxygen) */

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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate or not the CDO module
 */
/*----------------------------------------------------------------------------*/

int
cs_user_cdo_activated(void)
{
  /* CS_CDO_OFF     = -1 --> CDO schemes are not used (no activation)
     CS_CDO_WITH_FV =  0 --> CDO schemes are used as well as finite volume
     CS_CDO_ONLY    =  1 --> CDO schemes are exclusively used */

  return  CS_CDO_ONLY;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify additional mesh locations
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_add_mesh_locations(void)
{
  cs_boundary_zone_define("in", "x < 1e-5", 0);
  cs_boundary_zone_define("out", "x > 0.9999", 0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify for the computational domain:
 *         -- which type of boundaries closed the computational domain
 *         -- the settings for the time step
 *         -- activate predefined equations or modules
 *         -- add user-defined properties and/or advection fields
 *         -- add user-defined equations
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_init_setup(cs_domain_t   *domain)
{
  /* ======================
     Boundary of the domain
     ====================== */

  /* Choose a boundary by default.
     Valid choice is CS_PARAM_BOUNDARY_WALL or CS_PARAM_BOUNDARY_SYMMETRY */
  cs_domain_set_default_boundary(domain, CS_PARAM_BOUNDARY_WALL);

  /* Add a new boundary
     >> cs_domain_add_boundary(domain,
                               type_of_boundary,
                               mesh_location_name);

     * mesh_location_name is either a predefined mesh location or one defined
     by the user
     * type_of_boundary is one of the following keyword
        CS_PARAM_BOUNDARY_WALL,
        CS_PARAM_BOUNDARY_INLET,
        CS_PARAM_BOUNDARY_OUTLET,
        CS_PARAM_BOUNDARY_SYMMETRY
  */

  cs_domain_add_boundary(domain, CS_PARAM_BOUNDARY_INLET, "in");
  cs_domain_add_boundary(domain, CS_PARAM_BOUNDARY_OUTLET, "out");

  /* =========================
     Generic output management
     ========================= */

  cs_domain_set_output_param(domain,
                             10,     // output log frequency
                             2);     // verbosity (-1: no, 0, ...)

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

  /* ============================
     Add equations/model to solve
     ============================

     Activate predefined equations/modules

     Add user equation:
       default boundary condition is among:
       CS_PARAM_BC_HMG_DIRICHLET or CS_PARAM_BC_HMG_NEUMANN

     By default, initial values are set to zero (or the value given by the
     restart file in case of restart).
  */

  cs_walldistance_activate();

  cs_equation_add_user("AdvDiff",     // equation name
                       "Potential",   // associated variable field name
                       1,             // dimension of the unknown
                       CS_PARAM_BC_HMG_DIRICHLET); // default boundary

  /* ========================================
     Add material properties/advection fields
     ======================================== */

  cs_property_add("conductivity",      // property name
                  CS_PROPERTY_ANISO);  // type of material property
  cs_property_add("rho.cp",            // property name
                  CS_PROPERTY_ISO);    // type of material property

  cs_advection_field_add("adv_field"); // name of the new advection field
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

  cs_property_t  *conductivity = cs_property_by_name("conductivity");
  cs_real_33_t  tensor = {{1.0,  0.5, 0.0}, {0.5, 1.0, 0.5}, {0.0, 0.5, 1.0}};

  cs_property_def_aniso_by_value(conductivity, // property structure
                                 "cells",      // name of the volume zone
                                 tensor);      // values of the property

  cs_property_t  *rhocp = cs_property_by_name("rho.cp");
  cs_real_t  iso_val = 2.0;

  cs_property_def_iso_by_value(rhocp,    // property structure
                               "cells",  // name of the volume zone
                               iso_val);   // value of the property

  /* =============================
     User-defined advection fields
     ============================= */

  cs_adv_field_t  *adv = cs_advection_field_by_name("adv_field");

  cs_advection_field_def_by_analytic(adv, _define_adv_field, NULL);

  /* Enable also the defintion of the advection field at mesh vertices */
  cs_advection_field_set_option(adv, CS_ADVKEY_DEFINE_AT_VERTICES);

  /* Activate the post-processing of the related Courant number */
  cs_advection_field_set_option(adv, CS_ADVKEY_POST_COURANT);

  /* ======================
     User-defined equations
     ====================== */

  /* Retrieve the equation to set
     >> cs_equation_t  *eq = cs_equation_by_name("eq_name");  */

  cs_equation_t  *eq = cs_equation_by_name("AdvDiff");

  /* Define the boundary conditions
     >> cs_equation_add_bc_by_analytic(eq,
                                       bc_type,
                                       "zone_name",
                                       analytic_function);

     -> eq is the structure related to the equation to set
     -> type of boundary condition:
        CS_PARAM_BC_DIRICHLET, CS_PARAM_BC_HMG_DIRICHLET,
        CS_PARAM_BC_NEUMANN, CS_PARAM_BC_HMG_NEUMANN, CS_PARAM_BC_ROBIN

     >> cs_equation_add_bc_by_value(eq,
                                    bc_type,
                                    "mesh_location_name",
                                    values); // pointer
  */

  cs_equation_add_bc_by_analytic(eq,
                                 CS_PARAM_BC_DIRICHLET,
                                 "boundary_faces",  // zone name
                                 _define_bcs,       // pointer to the function
                                 NULL);             // input structure

  /* Link properties to different terms of this equation
     >> cs_equation_link(eq,
                         "term_keyword",
                         structure_to_link);

     - eq is the structure related to the equation to set
     - Keyword related to the term to set is a choice among:
       >> "diffusion", "time" or "advection"
     - If keyword is "time" or "diffusion", the structure to link is a
       property.
       If keyword is "advection", the structure to link is an advection field
  */

  /* Activate unsteady effect */
  cs_equation_link(eq, "time", rhocp);
  /* Activate diffusion effect */
  cs_equation_link(eq, "diffusion", conductivity);
  /* Activate advection effect */
  cs_equation_link(eq, "advection", adv);

  /* Add a source term:
     >> cs_equation_add_source_term_by_val(eq,
                                           volume_zone__name,
                                           val);
     or
     >> cs_equation_add_source_term_by_analytic(eq,
                                                volume_zone_name,
                                                analytic_func,
                                                input_pointer);


     where val is the value of the source term by m^3
     or where analytic_func is the name of the analytical function
   */

  cs_xdef_t  *st
    = cs_equation_add_source_term_by_analytic(eq,
                                              "cells",
                                              _define_source,
                                              NULL);

  /* Optional: specify the quadrature used for computing a source term

     >> key takes value among
     CS_QUADRATURE_BARY     barycenter approximation
     CS_QUADRATURE_HIGHER   4 Gauss points for approximating the integral
     CS_QUADRATURE_HIGHEST  5 Gauss points for approximating the integral
  */

  cs_xdef_set_quadrature(st, CS_QUADRATURE_BARY);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
