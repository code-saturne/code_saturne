#ifndef __CS_THERMAL_MODEL_H__
#define __CS_THERMAL_MODEL_H__

/*============================================================================
 * Base thermal model data.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Thermal model type
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_THERMAL_MODEL_NONE,
  CS_THERMAL_MODEL_TEMPERATURE,
  CS_THERMAL_MODEL_ENTHALPY,
  CS_THERMAL_MODEL_TOTAL_ENERGY,
  CS_THERMAL_MODEL_N_TYPES

} cs_thermal_model_variable_t;

typedef enum {

  CS_TEMPERATURE_SCALE_NONE,
  CS_TEMPERATURE_SCALE_KELVIN,
  CS_TEMPERATURE_SCALE_CELSIUS

} cs_temperature_scale_t;

/* thermal model descriptor */
/*--------------------------*/

typedef struct {

  int           itherm;      /* thermal model
                                - 0: no thermal model
                                - 1: temperature
                                - 2: enthalpy
                                - 3: total energy (only for compressible
                                     module) */
  int           itpscl;      /* temperature scale
                                - 0: none
                                - 1: Kelvin
                                - 2: Celsius */
  int           iscalt;      /* index of the thermal scalar (temperature,
                                energy or enthalpy), the index of the
                                corresponding variable is isca(iscalt) */

} cs_thermal_model_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to thermal model structure */

extern const cs_thermal_model_t  *cs_glob_thermal_model;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_thermal_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_thermal_model_t *
cs_get_glob_thermal_model(void);

/*----------------------------------------------------------------------------
 * Return thermal field (temperature, enthalpy, total energy according to
 * thermal model).
 *
 * returns:
 *   pointer to thermal field
 *----------------------------------------------------------------------------*/

cs_field_t *
cs_thermal_model_field(void);

/*----------------------------------------------------------------------------
 * Print the thermal model structure to setup.log.
 *----------------------------------------------------------------------------*/

void
cs_thermal_model_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_THERMAL_MODEL_H__ */
