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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_domain.h"

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
 * \brief  Specify additional mesh locations
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_add_mesh_locations(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify additional material properties
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_setup_properties(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */
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
}


/*----------------------------------------------------------------------------*/

END_C_DECLS
