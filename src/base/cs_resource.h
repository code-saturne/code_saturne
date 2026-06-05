#ifndef CS_RESOURCE_H
#define CS_RESOURCE_H

/*============================================================================
 * Resource allocation management (available time).
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/


/*=============================================================================
 * Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*--------------------------------------------------------------------------*/
/*!
 * \brief Set wall time limit check to false or true.
 */
/*--------------------------------------------------------------------------*/

void
cs_resource_wt_limit_check_set_status
(
  const bool status /*!<[in] status to set (activation or deactivation) */
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get current wall-clock time limit.
 *
 * \return current wall-time limit (in seconds), or -1
 */
/*----------------------------------------------------------------------------*/

double
cs_resource_get_wt_limit(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set wall-clock time limit.
 *
 * \param[in]  wt  wall-time limit (in seconds), or -1
 */
/*----------------------------------------------------------------------------*/

void
cs_resource_set_wt_limit(double  wt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Limit number of remaining time steps if the remaining allocated
 *        time is too small to attain the requested number of steps.
 *
 * \param[in]       ts_cur  current time step number
 * \param[in, out]  ts_max  maximum time step number
 */
/*----------------------------------------------------------------------------*/

void
cs_resource_get_max_timestep(int   ts_cur,
                             int  *ts_max);

/*----------------------------------------------------------------------------*/

#endif /* CS_RESOURCE_H */