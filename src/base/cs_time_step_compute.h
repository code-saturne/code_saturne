#ifndef __CS_TIME_STEP_COMPUTE__
#define __CS_TIME_STEP_COMPUTE__

/*============================================================================
 * Compute the local time step and the local Courant and Fourier number
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \file cs_time_step_compute.c
 *
 * \brief Compute the local time step and the local Courant and Fourier number.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
 * \param[in]     itrale        ALE iteration number
 */
/*----------------------------------------------------------------------------*/

void
cs_local_time_step_compute(int itrale);

/*----------------------------------------------------------------------------*/

void
cs_courant_fourier_compute(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TIME_STEP_COMPUTE__ */
