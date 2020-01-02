/*============================================================================
 * User head loss definitions.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_head_losses.c
 *
 * \brief User head loss definitions.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define head losses for a given volume zone.
 *
 * Head loss tensor coefficients for each cell are organized as follows:
 * ck11, ck22, ck33, ck12, ck13, ck23.
 *
 * Coefficients are set to zero (then computed based on definitions provided
 * through the GUI if this is the case) before calling this function, so
 * setting values to zero is usually not necessary, unless we want to fully
 * overwrite a GUI-based definition.
 *
 * Diagonal coefficients must be positive; the calculation may diverge
 * if this is not the case.
 *
 * \param[in]       zone  pointer to zone structure
 * \param[in, out]  cku   head loss coefficients
 */
/*----------------------------------------------------------------------------*/

void
cs_user_head_losses(const  cs_zone_t  *zone,
                    cs_real_t          cku[][6])
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
