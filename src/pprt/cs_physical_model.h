#ifndef __CS_PHYSICAL_MODEL_H__
#define __CS_PHYSICAL_MODEL_H__

/*============================================================================
 * General parameters management.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*! Specific physical model types */
/*--------------------------------*/

typedef enum {

  CS_PHYSICAL_MODEL_FLAG,          /*!< global specific physics flag */

  CS_COMBUSTION_3PT,               /*!< 3-point combustion model */
  CS_COMBUSTION_SLFM,              /*!< Steady laminar flamelet model */
  CS_COMBUSTION_EBU,               /*!< EBU combustion model */
  CS_COMBUSTION_LW,                /*!< Libby-Williams combustion model */
  CS_COMBUSTION_PCLC,              /*!< pulverized coal combustion with
                                        Lagrangian coupling */
  CS_COMBUSTION_COAL,              /*!< coal combustion model */
  CS_COMBUSTION_FUEL,              /*!< fuel combustion */
  CS_JOULE_EFFECT,                 /*!< Joule effect */
  CS_ELECTRIC_ARCS,                /*!< electric arcs */
  CS_COMPRESSIBLE,                 /*!< compressible model */
  CS_ATMOSPHERIC,                  /*!< atmospheric model */
  CS_COOLING_TOWERS,               /*!< cooling towers */
  CS_GAS_MIX,                      /*!< gas mix model */
  CS_GROUNDWATER,                  /*!< Richards model */

  CS_N_PHYSICAL_MODEL_TYPES        /*!< Number of physical model types */

} cs_physical_model_type_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Names of specific physical models */

extern int cs_glob_physical_model_flag[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PHYSICAL_MODEL_H__ */
