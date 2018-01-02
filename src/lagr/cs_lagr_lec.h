#ifndef __CS_LAGR_LEC_H__
#define __CS_LAGR_LEC_H__

/*============================================================================
 * Functions and types for lagrangian specific prints
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 * Fortran wrapper function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Fortran wrapper for restart files readings
 */
/*----------------------------------------------------------------------------*/

void
CS_PROCF(laglec, LAGLEC)(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fortran wrapper for restart files output.
 */
/*----------------------------------------------------------------------------*/

void
CS_PROCF(lagout, LAGOUT)(void);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 *\brief  Read Lagrangian restart files.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_lagrangian_checkpoint_read(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read Lagrangian particle and statistics restart files.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_restart_read_p(void);

/*--------------------------------------------------------------------*/
/*! \brief Output Lagrangian restart files.
 */
/*--------------------------------------------------------------------*/

void
cs_restart_lagrangian_checkpoint_write(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_LEC_H__ */
