#ifndef __CS_SADDLE_SOLVER_SETUP_H__
#define __CS_SADDLE_SOLVER_SETUP_H__

/*============================================================================
 * Routines to handle the SLES settings for a saddle-point problem
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
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_saddle_solver_setup.h

  \brief Routines to handle the setup of SLES relying on a \ref
         cs_param_saddle_t structure
*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the SLES structures (potentially several) related to a
 *        saddle-point system
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_setup_sles(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SADDLE_SOLVER_SETUP_H__ */
