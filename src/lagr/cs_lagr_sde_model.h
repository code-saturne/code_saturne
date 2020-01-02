#ifndef __CS_LAGR_SDE_MODEL_H__
#define __CS_LAGR_SDE_MODEL_H__

/*============================================================================
 * Lagrangian module physical model
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of particle stochastic differential equations
 *        for specific physical models.
 *
 * - fluid temperature seen by particles,
 * - particle temperature,
 * - particle diameter
 * - particle mass
 * - variables related to coal grains (Temp, MCH, MCK)
 * - additional user parameters
 *
 * \param[in]  tempct       thermal characteristic time
 * \param[out] cpgd1,cpgd2  devolatilisation terms 1 and 2
 * \param[out] cpght        heterogeneos combusion terms (coal with thermal
 *                          return coupling)
 */
/*------------------------------------------------------------------------- */

void
cs_lagr_sde_model(const cs_real_t  tempct[],
                  cs_real_t        cpgd1[],
                  cs_real_t        cpgd2[],
                  cs_real_t        cpght[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_SDE_MODEL_H__ */
