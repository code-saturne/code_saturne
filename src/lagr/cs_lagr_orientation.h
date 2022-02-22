#ifndef __CS_LAGR_ORIENTATION_H__
#define __CS_LAGR_ORIENTATION_H__

/*============================================================================
 * Functions and types for orientation
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of the equations of motion for particle orientation:
 *
 * Details in publication by M. Bossy and L. Campana (INRIA)
 *
 * \param[in] iprev     time step indicator for fields
 *                        0: use fields at current time step
 *                        1: use fields at previous time step
 * \param[in] dt_p      lagrangian time step
 * \param[in] gradvf    fluid velocity gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_orientation_dyn_spheroids(int              iprev,
                                  cs_real_t        dt_p,
                                  const cs_real_t  gradvf[][3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of the Jeffey equations in DNS mode
 *
 * \param[in] dt_p      lagrangian time step
 * \param[in] gradvf    fluid velocity gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_orientation_dyn_jeffery(cs_real_t        dt_p,
                                const cs_real_t  gradvf[][3][3]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_ORIENTATION_H__ */
