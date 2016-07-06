/*============================================================================
 * Additional post-processing functions defined by user related to CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_post.h"
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
 * \file cs_user_cdo_extra_op.c
 *
 * \brief Additional user-defined post-processing and analysis functions
*/
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initial step for user-defined operations on results provided by the
 *         CDO kernel.
 *
 * \param[in]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_start_extra_op(const cs_domain_t     *domain)
{
  CS_UNUSED(domain);

  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Additional user-defined operations on results provided by the CDO
 *         kernel. Define advanced post-processing and analysis for example.
 *
 * \param[in]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_extra_op(const cs_domain_t           *domain)
{
  CS_UNUSED(domain);

  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Final step for user-defined operations on results provided by the
 *         CDO kernel.
 *
 * \param[in]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_end_extra_op(const cs_domain_t     *domain)
{
  CS_UNUSED(domain);

  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
