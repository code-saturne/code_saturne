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

#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_order.h"
#include "cs_cdo.h"
#include "cs_cdo_toolbox.h"
#include "cs_param_eq.h"

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

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the specifications for the current simulation
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_setup(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* ===========================
     Define material properties
     ===========================

     By default, several material properties are defined:
     --> "Unity"
     --> "MassDensity"
     --> "LaminarViscosity"

     Users can define additional material properties
     1) Set a name specific to this material property
     2) Set a type of material property
     --> CS_PARAM_PTY_ISO
     --> CS_PARAM_PTY_ORTHO
     --> CS_PARAM_PTY_ANISO
     3) Set the frequency of post-processing:
     --> -1: no post-processing
     -->  0: at the beginning of the simulation
     -->  n: at each n iteration(s)

     In a second step, please set the value of the material property.
     There exist several ways to set this value:
     --> By value (the simplest): direclty set the value. Possible when the
         material property is uniform and steady.
     --> By user function
     --> By field
     ...

  */

  cs_param_pty_add("Conductivity",     // property name
                   CS_PARAM_PTY_ISO,   // type of material property
                   0);                 // frequency of post-processing

  cs_get_t  conductivity;
  conductivity.val = 1.0;
  cs_param_pty_set_by_val("Conductivity", // property name
                          conductivity);  // value of the material property

  /* ===========================
     Define mesh locations
     =========================== */

  /* Define two new mesh locations */
  cs_mesh_location_add("Left",
                       CS_MESH_LOCATION_BOUNDARY_FACES,
                       "x < 1e-6");

  cs_mesh_location_add("Right",
                       CS_MESH_LOCATION_BOUNDARY_FACES,
                       "x > 0.99999");

  /* ===========================
     Define equations to solve
     ===========================

     Specified the behaviour of each additional equation
     1) Indicate a name specific to each equation
     2) Indicate a name specific for the variable related to this equation
     3) Set a type of equation:
     --> CS_PARAM_EQ_TYPE_SCAL
     --> CS_PARAM_EQ_TYPE_VECT
     --> CS_PARAM_EQ_TYPE_TENS (TODO)

     4) Choose which term plays a role in this equation
     --> Unsteady term:   true/false
     --> Diffusion term:  true/false
     --> Convection term: true/false

     5) Set a (basic) boundary condition (BC) by default on all the boundary
     --> CS_PARAM_BC_BASIC_NONE (no BC defined, useful for checkings)
     --> CS_PARAM_BC_HMG_DIRICHLET (set the unknow to 0)
     --> CS_PARAM_BC_HMG_NEUMANN (set the flux to 0)

  */

  cs_param_eq_add("Laplace",                // equation
                  "Potential",              // variable name
                  CS_PARAM_EQ_TYPE_SCAL,    // type of equation
                  true,                     // steady ?
                  false,                    // convection term ?
                  true,                     // diffusion term ?
                  CS_PARAM_BC_HMG_NEUMANN); // default BC

  /* Define the boundary conditions */
  cs_param_eq_add_scalbc_by_val("Laplace", "Left",  CS_PARAM_BC_DIRICHLET, 0.0);
  cs_param_eq_add_scalbc_by_val("Laplace", "Right", CS_PARAM_BC_DIRICHLET, 1.0);

  /* Associate the property related to the diffusion term to an equation */
  cs_param_eq_set_diffusion_pty("Laplace",        // equation
                                "Conductivity");  // material pty

  /* Add a source term: There are several types of source terms
     --> CS_PARAM_SOURCE_TERM_BASIC (explicit = in the right-hand side)
     --> CS_PARAM_SOURCE_TERM_IMPLICIT (in the matrix to invert)
     --> CS_PARAM_SOURCE_TERM_IMEX (two contributions: implicit and explicit)
     --> CS_PARAM_SOURCE_TERM_MASS
     --> CS_PARAM_SOURCE_TERM_HEADLOSS
  */

  // Example 1
  cs_get_t  imp_part, exp_part;

  imp_part.val = 0.0;
  exp_part.val = 0.0;
  cs_param_eq_add_source_term_by_val("Laplace",                // equation
                                     "SourceTerm",             // label
                                     "cells",                  // mesh location
                                     CS_PARAM_SOURCE_TERM_BASIC, // type
                                     imp_part,                 // implicit part
                                     exp_part);                // explicit part

  // Example 2
  /* cs_param_eq_add_source_term_by_analytic("Laplace",              // equation */
  /*                                         "SourceTerm",           // label */
  /*                                         "cells",                // location */
  /*                                         true,                   // post ? */
  /*                                         CS_PARAM_SOURCE_TERM_STD, // type */
  /*                                         CS_QUADRATURE_BARY,     // quadrature */
  /*                                         laplace_st);            // function */

}


/*----------------------------------------------------------------------------*/

END_C_DECLS
