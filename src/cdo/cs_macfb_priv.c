/*============================================================================
 * Common functions between scalar-valued and vector-valued CDO face-based
 * schemes
 *============================================================================*/

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
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_macfb_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_MACFB_PRIV_DBG 0

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the advection-related parameters in the context structure of
 *         MAC face-based schemes
 *
 * \param[in]      eqp    pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb    pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] eqc    pointer to a cs_macfb_priv_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_set_advection_function(const cs_equation_param_t *eqp,
                                cs_equation_builder_t     *eqb,
                                cs_macfb_priv_t           *eqc)
{
  if (eqc == NULL || eqb == NULL)
    return;

  assert(eqp != NULL);

  if (cs_equation_param_has_convection(eqp) == false)
    return;

  bft_error(__FILE__,
            __LINE__,
            0,
            "%s: Advection not implemented for MAC-fb sheme",
            __func__);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
