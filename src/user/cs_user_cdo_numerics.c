/*============================================================================
 * Set advanced numerical parameters for the current simulation
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_quantities.h"
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

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_cdo_numerics.c
 *
 * \brief Set advanced parameters about the numerical schemes for each
 *        equation to solve.
 *        Useful to change the default behaviour.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the way geometric quantities
 *         are built
 *
 * \return the type of computation to evaluate the cell center
 */
/*----------------------------------------------------------------------------*/

cs_cdo_cell_center_algo_t
cs_user_cdo_geometric_settings(void)
{
  /* Algorithm for computing cell centers */
  /* ==================================== */

  /* Choice between:
     CS_CDO_CCENTER_MEANV: Cell center is computed as the mean of cell vertices
     CS_CDO_CCENTER_BARYC: Cell center is computed as the real cell barycenter
     CS_CDO_CCENTER_SATURNE: Cell center is given by Code_Saturne
   */

  return CS_CDO_CCENTER_SATURNE;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the numerical parameters
 *         of the equation resolved during the computation
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_numeric_settings(void)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
