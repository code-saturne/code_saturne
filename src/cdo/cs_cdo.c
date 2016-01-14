/*============================================================================
 * General functions or variables for the INNOV module
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
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Global variables
 *============================================================================*/

/* Separation lines: long, medium, short */
const char lsepline[] =
  "# =======================================================================\n";
const char msepline[] =
  "# =========================================\n";
const char ssepline[] =
  "# =================\n";

/*=============================================================================
 * Local static variables
 *============================================================================*/

static double  cs_base_eps_machine;
static double  cs_base_zthreshold = FLT_MIN;

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
 * \brief  Compute epsilon which is the machine precision
 */
/*----------------------------------------------------------------------------*/

void
cs_set_eps_machine(void)
{
  double  y;

  double  eps = 5e-16;

  y = 1.0 + eps;

  while (y > 1.0) {
    eps /= 2.0;
    y = 1.0 + eps;
  }
  eps *= 2.0;
  printf("# Machine precision: %12.8e\n", eps);

  cs_base_eps_machine = eps;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the machine precision
 */
/*----------------------------------------------------------------------------*/

double
cs_get_eps_machine(void)
{
  return cs_base_eps_machine;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the threshold under which one considers it's zero
 */
/*----------------------------------------------------------------------------*/

double
cs_get_zero_threshold(void)
{
  return cs_base_zthreshold;
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
  if (fabs(magnitude) > cs_base_zthreshold) {

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
