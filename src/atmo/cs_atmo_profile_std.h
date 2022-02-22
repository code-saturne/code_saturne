#ifndef __CS_ATMO_PROFILE_STD_H__
#define __CS_ATMO_PROFILE_STD_H__

/*============================================================================
 * Compute standard atmospheric profile
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
 * \brief compute standard atmospheric profile (Holton p 374)
 *
 * \param[in]       z          absolute altitude in m
 * \param[out]      p          pressure in pa
 * \param[out]      t          temperature in k
 * \param[out]      r          density in kg/m3
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_profile_std(cs_real_t   z,
                    cs_real_t  *p,
                    cs_real_t  *t,
                    cs_real_t  *r);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ATMO_PROFILE_STD_H__ */
