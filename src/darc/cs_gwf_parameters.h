#ifndef __CS_GWF_PARAMETERS_H__
#define __CS_GWF_PARAMETERS_H__

/*============================================================================
 * General parameters management for groundwater flow module.
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Parameter check behavior when an error is detected */

/*----------------------------------------------------------------------------
 * Structure defining the sorption model
 *----------------------------------------------------------------------------*/

typedef struct {

  int     kinetic;       /* 0 : sorption at equilibirum
                            1 : sorption with kinetic   */
  int     ikd;           /* id for kd */
  int     idel;          /* id for delay */
  int     ikp;           /* id for kplus */
  int     ikm;           /* id for kminus */
  int     imxsol;        /* id for solubility index */
  int     anai;          /* 0: solve EK model in an explicit way
                            1: use analytical solution for EK model */

} cs_gwf_soilwater_partition_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define field keys for the ground water flow module.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_parameters_define_field_keys(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_PARAMETERS_H__ */
