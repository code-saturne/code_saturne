#ifndef __CS_INTPRF_H__
#define __CS_INTPRF_H__

/*============================================================================
 * Temporal and z-axis interpolation for meteorological profiles.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Temporal and z-axis interpolation for meteorological profiles
 *
 * An optimized linear interpolation is used.
 *
 * \param[in]   nprofz      total number of measure points
 * \param[in]   nproft      total number of time values
 * \param[in]   profz       z coordinates of measure points
 * \param[in]   proft       physical times of dataset acquisition
 * \param[in]   profv       measured values
 * \param[in]   xz          interpolation elevation
 * \param[in]   t           interpolation time
 *
 * \return  interpolated value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_intprf(int              nprofz,
          int              nproft,
          const cs_real_t  profz[],
          const cs_real_t  proft[],
          const cs_real_t  profv[],
          cs_real_t        xz,
          cs_real_t        t);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Z-axis interpolation for meteorological profiles
 *
 * An optimized linear interpolation is used.
 *
 * \param[in]   nprofz      total number of measure points
 * \param[in]   profz       z coordinates of measure points
 * \param[in]   profv       measured values
 * \param[in]   xz          interpolation elevation
 * \param[out]  z_lv        index of lower and upel level (optional, size 2),
 *                          or NULL
 * \param[out]  var         interpolated value
 */
/*----------------------------------------------------------------------------*/

void
cs_intprz(int               nprofz,
          const cs_real_t   profz[],
          const cs_real_t   profv[],
          cs_real_t         xz,
          int              *z_lv,
          cs_real_t        *var);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_INTPRF_H__ */
