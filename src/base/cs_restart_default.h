#ifndef __CS_RESTART_DEFAULT_H__
#define __CS_RESTART_DEFAULT_H__

/*============================================================================
 * Checkpoint/restart handling for default application.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
#include "cs_field.h"
#include "cs_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read boundary condition coefficients for all fields from checkpoint.
 *
 * parameters:
 *   restart <-> associated restart file pointer
 *----------------------------------------------------------------------------*/

void
cs_restart_read_bc_coeffs(cs_restart_t  *restart);

/*----------------------------------------------------------------------------
 * Write boundary condition coefficients for all fields to checkpoint.
 *
 * parameters:
 *   restart <-> associated restart file pointer
 *----------------------------------------------------------------------------*/

void
cs_restart_write_bc_coeffs(cs_restart_t  *restart);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RESTART_DEFAULT_H__ */
