/*============================================================================
 * General functions or variables for the CDO module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include <float.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"

#include "cs_defs.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Global variables
 *============================================================================*/

/* Activation of the CDO/HHO module */
_Bool  cs_cdo_is_activated = false;

/* Separation lines: long, medium, short */
const char lsepline[80] =
  "# =======================================================================\n";
const char msepline[60] =
  "# =========================================\n";
const char ssepline[40] =
  "# =================\n";

/* Default locations */
const cs_flag_t  cs_cdo_primal_vtx  = CS_FLAG_PRIMAL | CS_FLAG_VERTEX;
const cs_flag_t  cs_cdo_primal_face = CS_FLAG_PRIMAL | CS_FLAG_FACE;
const cs_flag_t  cs_cdo_primal_cell = CS_FLAG_PRIMAL | CS_FLAG_CELL;
const cs_flag_t  cs_cdo_dual_vtx  = CS_FLAG_DUAL | CS_FLAG_VERTEX;
const cs_flag_t  cs_cdo_dual_face = CS_FLAG_DUAL | CS_FLAG_FACE;
const cs_flag_t  cs_cdo_dual_cell = CS_FLAG_DUAL | CS_FLAG_CELL;
const cs_flag_t  cs_cdo_dual_face_byc =
  CS_FLAG_DUAL | CS_FLAG_FACE | CS_FLAG_BY_CELL;

/*=============================================================================
 * Local static variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a string "true" or "false" according to the boolean
 *
 * \param[in]  boolean  bool type
 *
 * \return a string "true" or "false"
 */
/*----------------------------------------------------------------------------*/

const char *
cs_base_strtf(bool  boolean)
{
  if (boolean)
    return "true";
  else
    return "false";
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_nvec3_t structure from a cs_real_3_t
 *
 * \param[in]  v     vector of size 3
 * \param[out] qv    pointer to a cs_nvec3_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_nvec3(const cs_real_3_t    v,
         cs_nvec3_t          *qv)
{
  cs_real_t  magnitude = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  qv->meas = magnitude;
  if (fabs(magnitude) > cs_math_zero_threshold) {

    cs_real_t  inv = 1/magnitude;
    for (int k = 0; k < 3; k++)
      qv->unitv[k] = inv * v[k];

  }
  else
    for (int k = 0; k < 3; k++)
      qv->unitv[k] = 0;

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
