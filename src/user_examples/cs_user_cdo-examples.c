/*============================================================================
 * Set main parameters for the current simulation when the CDO kernel is used
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include "cs_cdo_toolbox.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_user_cdo.c

  \brief  Set main parameters for the current simulation when the CDO kernel
          is used

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! \endcond (end ignore by Doxygen) */

static const double  one_third = 1./3.;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/* ---------------------------------------------------------------------------
 * FVAC6: TEST 1 -- Boundary conditions
 * Tbc =  1 + sin(pi*x_fac) * sin(pi*(y_fac+0.5)) * sin(pi*(z_fac+1/3))
 * -------------------------------------------------------------------------- */

static void
_define_bcs(cs_real_t     time,
            cs_real_3_t   xyz,
            cs_get_t     *get)
{
  double  bcval = 0.0;

  const double  x = xyz[0], y = xyz[1], z = xyz[2];
  const double  pi = 4.0*atan(1.0);

  bcval = 1 + sin(pi*x) * sin(pi*(y + 0.5)) * sin(pi*(z + one_third));

  (*get).val = bcval;
}

/* ---------------------------------------------------------------------------
 * FVAC6: TEST 1 -- Source term
 * material property is constant
 * -------------------------------------------------------------------------- */

static void
_define_source_term(cs_real_t     time,
                    cs_real_3_t   xyz,
                    cs_get_t     *get)
{
  cs_real_33_t  cond;

  double  stval = 0.0;
  double  sx = DBL_MAX, sy = DBL_MAX, sz = DBL_MAX;
  double  sxx = DBL_MAX, syy = DBL_MAX, szz = DBL_MAX;
  double  sxy = DBL_MAX, sxz = DBL_MAX, syz = DBL_MAX;

  const double  x = xyz[0], y = xyz[1], z = xyz[2];
  const double  pi = 4.0*atan(1.0), pi2 = pi*pi;
  const double  cpx = cos(pi*x), spx = sin(pi*x);
  const double  cpy = cos(pi*(y+0.5)), spy = sin(pi*(y+0.5));
  const double  cpz = cos(pi*(z+one_third)), spz = sin(pi*(z+one_third));

  cond[0][0] = 1.0, cond[0][1] = 0.5, cond[0][2] = 0.0;
  cond[1][0] = 0.5, cond[1][1] = 1.0, cond[1][2] = 0.5;
  cond[2][0] = 0.0, cond[2][1] = 0.5, cond[2][2] = 1.0;

  /* first derivatives */
  sx = pi*cpx*spy*spz, sy = pi*spx*cpy*spz, sz = pi*spx*spy*cpz;

  /* second derivatives */
  sxx = syy = szz = -pi2*spx*spy*spz;
  sxy = pi2*cpx*cpy*spz, sxz = pi2*cpx*spy*cpz, syz = pi2*spx*cpy*cpz;

  stval = cond[0][0]*sxx + cond[1][1]*syy + cond[2][2]*szz +
    2*(cond[0][1]*sxy + cond[0][2]*sxz + cond[1][2]*syz);
  stval *= -1;

  (*get).val = stval;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify additional mesh locations
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_add_mesh_locations(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* ===========================
     Define mesh locations
     ===========================

     By default several mesh locations are predefined
     >> "cells"
     >> "interior faces"
     >> "boundary faces"
     >> "vertices"

 */

  //  cs_mesh_location_add("Bottom", CS_MESH_LOCATION_BOUNDARY_FACES, "x < 1e-6");

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add user-defined material properties and/or advection fields
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_add_properties(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* ================================
     User-defined material properties
     ================================

     By default, several material properties are defined:
     >> "unity"
     >> "mass_density"
     >> "laminar_viscosity"

     Users can also define additional material properties
     1) Set a name specific to this material property
     2) Set a type among the following choices:
        >> "isotropic", "orthotropic" or "anisotropic"
     3) Set the frequency of post-processing:
        >> -1: no post-processing
        >>  0: at the beginning of the simulation
        >>  n: at each n iteration(s)
  */

  cs_param_pty_add("conductivity",  // property name
                   "anisotropic",   // type of material property
                   -1);             // frequency of post-processing

  /* =============================
     User-defined advection fields
     =============================

     Users can also define additional material properties
     1) Set a name specific to this convection field
     3) Set the frequency of post-processing:
        >> -1: no post-processing
        >>  0: at the beginning of the simulation
        >>  n: at each n iteration(s)
  */

  cs_param_adv_field_add("adv_field",  // property name
                         0);           // frequency of post-processing

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify the definition of additional material properties and/or
 *         advection fields
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_set_properties(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* In a second step, please set the value of the material property.

     1) Give the name of the property
     2) Set the type of definition among the following choices:
        >> "value", "analytic", "user"
     3) Set the value
        >> "1.0"     for an isotropic property set by value
        >> _my_func  name of the function for a property set by analytic
  */

  cs_param_pty_set("conductivity", // property name
                   "value",
                   "1.0  0.5  0.0\n" // value of the material property
                   "0.5  1.0  0.5\n"
                   "0.0  0.5  1.0\n");

  /* In a second step, please set the value of the advection field.

     1) Give the name of the property
     2) Set the type of definition among the following choices:
        >> "value", "analytic", "user"
     3) Set the value
        >> "1.0"     for an isotropic property set by value
        >> _my_func  name of the function for a property set by analytic
  */

  cs_param_adv_field_set("adv_field",        // property name
                         "value",            // type of definition
                         "1.0  0.0  0.0\n"); // value of the advection field

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify which type of boundaries closed the computational domain
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_setup_domain_boundary(cs_domain_t   *domain)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* Choose a boundary by default
     boundary keyword is one of the following keyword
     >> wall or symmetry
  */

  cs_domain_set_default_boundary(domain, "wall");

  /* Add a boundary
     >>   cs_domain_add_boundary(domain,
                                 mesh location name,
                                 boundary keyword)

     mesh location name is either a predefined mesh location or one defined
     by user

     boundary keyword is one of the following keyword
     >> wall, inlet, outlet, symmetry
  */

  //  cs_domain_add_boundary(domain, "Top", "wall");

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify which are the equatinos to be solved in this computational
 *         domain
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_add_domain_equations(cs_domain_t   *domain)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* =========================================
     Activate predefined equations
     =========================================

     For the wall distance: cs_domain_activate_wall_distance(domain);

  */

  cs_domain_activate_wall_distance(domain);

  /* =========================================
     Define additional user equations to solve
     =========================================

     cs_domain_add_user_equation
     >> arguements: domain,
                    equation name,
                    associated field name,
                    type of equation: "scalar", "vector" or "tensor"
                    is steady ?       true or false
                    do_convection ?   true or false
                    do_diffusion ?    true or false
                    default_bc:       "zero_value" or "zero_flux"
  */

  cs_domain_add_user_equation(domain,
                              "FVCA6.1",     // equation name
                              "Potential",   // associated field name
                              "scalar",      // type of equation
                              true,         // steady ?
                              false,         // convection ?
                              true,          // diffusion ?
                              "zero_value"); // default boundary condition

  cs_domain_add_user_equation(domain,
                              "CODITS.FVCA", // equation name
                              "Pot2",        // associated field name
                              "scalar",      // type of equation
                              true,          // steady ?
                              true,          // convection ?
                              true,          // diffusion ?
                              "zero_value"); // default boundary condition

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate material property to user-defined equations and specify
 *         boundary conditions, source terms for thes additional equations
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_setup_equations(cs_domain_t   *domain)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  cs_equation_t  *eq = NULL;

  eq = cs_domain_get_equation(domain, "FVCA6.1");

  if (eq != NULL) { // Equation has been found

    /* Define the boundary conditions (BCs) for equation eq

       boundary conditions are among the following choices:
       >> "dirichlet", "neumann" or "robin"

       type of definition is among the following choices:
       >> "value", "analytic", "user"

     */

    cs_equation_add_bc(eq,                // equation
                       "boundary_faces",  // name of the mesh location
                       "dirichlet",       // BC type
                       "analytic",        // type of definition
                       _define_bcs);      // pointer to the analytic function

    /* Associate properties with terms at play in the equation */

    cs_equation_link(eq,              // equation
                     "diffusion",     // for the diffusion term
                     "conductivity"); // property name

    /* Add a source term: There are several types of source terms

       label of the source term is optional (i.e. NULL is possible)

       type of the source is among the following choices:
       >> "implicit", "explicit", "imex"
          >> "implicit" means that the source term is added to the diagonal of
          the linear system
          >> "explicit" means that the source term is added to the right-hand
          side of the linear system
          >> "imex" means that there are two contributions: the first one is
          implicit and the second one is explicit (Not implemented yet)

       type of definition is among the following choices:
       >> "value", "analytic", "user"

     */

    cs_equation_add_source_term(eq,
                                "SourceTerm",    // label of the source term
                                "cells",         // name of the mesh location
                                "explicit",      // type of source term
                                "analytic",      // type of definition
                                _define_source_term); // analytic function

  } /* eq != NULL */

  eq = cs_domain_get_equation(domain, "CODITS.FVCA");

  if (eq != NULL) { // Equation has been found

    /* Define the boundary conditions (BCs) for equation eq

       boundary conditions are among the following choices:
       >> "dirichlet", "neumann" or "robin"

       type of definition is among the following choices:
       >> "value", "analytic", "user"

     */

    cs_equation_add_bc(eq,                // equation
                       "boundary_faces",  // name of the mesh location
                       "dirichlet",       // BC type
                       "analytic",        // type of definition
                       _define_bcs);      // pointer to the analytic function

    /* Associate properties with terms at play in the equation */

    cs_equation_link(eq,              // equation
                     "diffusion",     // for the diffusion term
                     "conductivity"); // property name

    cs_equation_link(eq,
                     "advection",
                     "adv_field");

    /* Add a source term: There are several types of source terms

       label of the source term is optional (i.e. NULL is possible)

       type of the source is among the following choices:
       >> "implicit", "explicit", "imex"
          >> "implicit" means that the source term is added to the diagonal of
          the linear system
          >> "explicit" means that the source term is added to the right-hand
          side of the linear system
          >> "imex" means that there are two contributions: the first one is
          implicit and the second one is explicit (Not implemented yet)

       type of definition is among the following choices:
       >> "value", "analytic", "user"

     */

    cs_equation_add_source_term(eq,
                                "SourceTerm",    // label of the source term
                                "cells",         // name of the mesh location
                                "explicit",      // type of source term
                                "analytic",      // type of definition
                                _define_source_term); // analytic function

  } /* eq != NULL */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
