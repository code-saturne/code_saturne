/*============================================================================
 * Additional post-processing functions defined by user related to CDO schemes
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_post.h"
#include "cs_cdo.h"
#include "cs_quadrature.h"
#include "cs_param.h"
#include "cs_param_eq.h"
#include "cs_cdovb_codits.h"
#include "cs_cdofb_codits.h"

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
  \file cs_user_cdo_post.c

  \brief Additional user-defined post-processing and analysis functions
*/

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
 * \brief  Additional operations on results provided by CDO schemes.
 *         Define advanced post-processing and/or analysis for instance.
 *
 * \param[in]  m        pointer to a cs_mesh_t struct.
 * \param[in]  mq       pointer to a cs_quantities_t struct.
 * \param[in]  connect  pointer to a cs_cdo_connect_t structure
 * \param[in]  quant    pointer to a cs_cdo_quantities_t structure
 * \param[in]  tcur     current physical time of the simulation
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_extra_op(const cs_mesh_t             *m,
                     const cs_mesh_quantities_t  *mq,
                     const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     double                       tcur)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
