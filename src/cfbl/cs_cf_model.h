#ifndef __CS_CF_MODEL_H__
#define __CS_CF_MODEL_H__

/*============================================================================
 * Thermodynamic laws for the compressible module
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/* compressible model general options descriptor */
/*-----------------------------------------------*/

typedef struct {

  int           ieos;         /* indicator of equation of state */

  int           ithvar;       /* indicator for thermodynamic
                                 variables initialization */

  int           icfgrp;       /* indicator for hydrostatic balance
                                 in boundary conditions */

  double        psginf;       /* stiffened gas limit pressure (zero in
                                 perfect gas) (Pa) for single phase model */

  double        gammasg;      /* stiffened gas polytropic coefficient,
                                 (dimensionless) for single phase model */

  int           hgn_relax_eq_st;  /* source term step:
                                     - -1 disabled
                                     -  0 enabled */

} cs_cf_model_t;

typedef enum {

  CS_EOS_NONE                  = -1,
  CS_EOS_IDEAL_GAS             = 1,
  CS_EOS_STIFFENED_GAS         = 2,
  CS_EOS_GAS_MIX               = 3,
  CS_EOS_HOMOGENEOUS_TWO_PHASE = 4,
  CS_EOS_MOIST_AIR             = 5

} cs_cf_model_eos_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* pointer to main compressible model descriptor structure */

extern const cs_cf_model_t         *cs_glob_cf_model;

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to compressible model global structure cs_glob_cf_model
 */
/*----------------------------------------------------------------------------*/

cs_cf_model_t *
cs_get_glob_cf_model(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to boundary face indicator array of convection flux
 *        - 0 upwind scheme
 *        - 1 imposed flux
 */
/*----------------------------------------------------------------------------*/

int *
cs_cf_get_icvfli(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to imposed thermal flux indicator at the boundary
 *        (some boundary contributions of the total energy eq. have to be
 *         cancelled)
 */
/*----------------------------------------------------------------------------*/

int *
cs_cf_get_ifbet(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Property field definitions for the compressible module,
 *        according to calculation type selected by the user.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_add_property_fields(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the compressible module options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_model_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CF_MODEL_H__ */
