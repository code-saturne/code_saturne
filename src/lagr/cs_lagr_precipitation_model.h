#ifndef __CS_LAGR_PRECIPITATION_H__
#define __CS_LAGR_PRECIPITATION_H__

/*============================================================================
 * Functions and types for precipitation model
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*============================================================================
 * Public function prototypes
 *============================================================================*/

void
CS_PROCF(precst, PRECST)(cs_real_t *dtref,
                         cs_real_t *crom,
                         cs_real_t *cvar_scal,
                         cs_real_t  crvexp[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Management of the injection of particles formed by precipitation.
 *
 * \param[in]   vela       pointer to fluid velocity array (per cell)
 * \param[out]  val        number of particles to inject (with weight)
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_precipitation_injection(cs_real_t   *vela,
                                cs_real_t   *val);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR__H_PRECIPITATION_ */
