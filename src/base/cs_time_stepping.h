#ifndef __CS_TIME_STEPPING_H__
#define __CS_TIME_STEPPING_H__

/*============================================================================
 * Main time loop.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \file cs_time_stepping.cpp
 *
 * \brief Main time loop.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Main time loop.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_stepping(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output a checkpoint.
 *
 * If needed, the mesh is also output in the checkpoint directory,
 * exect if this function is called for checkpoint serialized in memory
 * (which is a special case for FMI exchange).
 *
 * \param[in]  checkpoint_mesh  also save mesh in checkpoint directory
 */
/*----------------------------------------------------------------------------*/

void
cs_time_stepping_write_checkpoint(bool  checkpoint_mesh);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TIME_STEPPING_H__ */
