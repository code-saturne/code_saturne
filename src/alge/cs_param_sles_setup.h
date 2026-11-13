#ifndef CS_PARAM_SLES_SETUP_H
#define CS_PARAM_SLES_SETUP_H

/*============================================================================
 * Routines to handle the SLES settings relying on a cs_param_sles_t
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "alge/cs_param_sles.h"

/*----------------------------------------------------------------------------*/

/*!
  \file cs_param_sles_setup.h

  \brief Routines to handle the setup of SLES relying on a \ref cs_param_sles_t
         structure
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
/*
 * Define cs_sles_t structure in accordance with the settings of a
 * cs_param_sles_t structure (SLES = Sparse Linear Equation Solver).
 */
/*----------------------------------------------------------------------------*/

int
cs_param_sles_setup(bool             use_field_id,
                    cs_param_sles_t *slesp);

/*----------------------------------------------------------------------------*/
/*
 * Update the settings associated to a cs_sles_t structure and apply those
 * defined in the given cs_param_sles_t structure.
 * This function is used only when a first setup has been performed.
 *
 * One modifies only some specific options like the max. number of iterations
 * or the relative tolerance
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_cvg_param(bool                   use_field_id,
                              const cs_param_sles_t *slesp);

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*
 * Set the command line option for PETSc.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_petsc_cmd(bool        use_prefix,
                              const char *prefix,
                              const char *keyword,
                              const char *keyval);

/*----------------------------------------------------------------------------*/
/*
 * Set a KSP structure in PETSc. This is equivalent to set a solver and its
 * related preconditioner.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_petsc_ksp(const char      *label,
                              cs_param_sles_t *slesp,
                              void            *p_ksp);

/*----------------------------------------------------------------------------*/
/*
 * Set an AMG preconditioner in PETSc.
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_petsc_pc_amg(const char      *prefix,
                                 cs_param_sles_t *slesp,
                                 void            *p_pc);
#endif

/*----------------------------------------------------------------------------*/

#endif /* CS_PARAM_SLES_SETUP_H */
