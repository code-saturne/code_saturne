#ifndef __CS_WALLDISTANCE_H__
#define __CS_WALLDISTANCE_H__

/*============================================================================
 * Compute the wall distance using the CDO framework
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if the computation of the wall distance is activated
 */
/*----------------------------------------------------------------------------*/

bool
cs_walldistance_is_activated(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the future computation of the wall distance
 */
/*----------------------------------------------------------------------------*/

void
cs_walldistance_activate(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the equation related to the wall distance
 */
/*----------------------------------------------------------------------------*/

void
cs_walldistance_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup stage for the equation related to the wall
 *         distance. Only useful for Hamilton-Jacobi equation
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_walldistance_finalize_setup(const cs_cdo_connect_t       *connect,
                               const cs_cdo_quantities_t    *cdoq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the wall distance
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_walldistance_compute(const cs_mesh_t              *mesh,
                        const cs_time_step_t         *time_step,
                        const cs_cdo_connect_t       *connect,
                        const cs_cdo_quantities_t    *cdoq);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_WALLDISTANCE_H__ */
