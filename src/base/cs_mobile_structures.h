#ifndef __CS_MOBILE_STRUCTURES_H__
#define __CS_MOBILE_STRUCTURES_H__

/*============================================================================
 * Mobile structures management.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_restart.h"
#include "cs_time_control.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Maximum number of implicitation iterations of the structure displacement */
extern int cs_glob_mobile_structures_i_max;

/*! Relative precision of implicitation of the structure displacement */
extern double cs_glob_mobile_structures_i_eps;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize mobile structures with ALE for internal coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize mobile structures with ALE for internal coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize mobile structures with ALE for internal coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log structures and coupling information
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query number of internal mobile structures defined.
 *
 * \return  number of internal mobile structures
 */
/*----------------------------------------------------------------------------*/

int
cs_mobile_structures_get_n_structures(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add internal mobile structures.
 *
 * This function may be called multiple time to change the number of
 * mobile structures.
 *
 * \param[in]   n_structures  number of internal mobile structures
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_add_n_structures(int  n_structures);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add external mobile structures.
 *
 * This function may be called multiple time to change the number of
 * mobile structures.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_add_external_structures(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Newmark coefficients for internal mobile structures.
 *
 * \param[in]   alpha  alpha coefficient for Newmark algorithm
 * \param[in]   beta   beta coefficient for Newmark algorithm
 * \param[in]   gamma  gamma coefficient for Newmark algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_set_newmark_coefficients(cs_real_t  alpha,
                                              cs_real_t  beta,
                                              cs_real_t  gamma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predict displacement of mobile structures with ALE.
 *
 * \param[in]   itrale   ALE iteration number
 * \param[in]   italim   implicit coupling iteration number
 * \param[in]   ineefl   indicate whether fluxes should be saved
 * \param[out]  impale   imposed displacement indicator
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_prediction(int  itrale,
                                int  italim,
                                int  ineefl,
                                int  impale[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Displacement of mobile structures with ALE for internal coupling.
 *
 * \param[in]       itrale   ALE iteration number
 * \param[in]       italim   implicit coupling iteration number
 * \param[in, out]  itrfin   indicator for last iteration of implicit coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_displacement(int   itrale,
                                  int   italim,
                                  int  *itrfin);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read mobile structures data to checkpoint.
 *
 * \param[in, out]  r   associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_restart_read(cs_restart_t  *r);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write mobile structures data to checkpoint.
 *
 * \param[in, out]  r   associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_restart_write(cs_restart_t  *r);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MOBILE_STRUCTURES_H__ */
