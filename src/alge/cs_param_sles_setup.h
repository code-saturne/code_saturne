#ifndef __CS_PARAM_SLES_SETUP_H__
#define __CS_PARAM_SLES_SETUP_H__

/*============================================================================
 * Routines to handle the SLES settings relying on a cs_param_sles_t
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "alge/cs_param_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
/*!
 * \brief Define cs_sles_t structure in accordance with the settings of a
 *        cs_param_sles_t structure (SLES = Sparse Linear Equation Solver)
 *
 * \param[in]      use_field_id  if false use system name to define a SLES
 * \param[in, out] slesp         pointer to a cs_param_sles_t structure
 *
 * \return an error code (-1 if a problem is encountered, 0 otherwise)
 */
/*----------------------------------------------------------------------------*/

int
cs_param_sles_setup(bool              use_field_id,
                    cs_param_sles_t  *slesp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the settings associated to a cs_sles_t structure and apply
 *        those defined in the given cs_param_sles_t structure.
 *        This function is used only when a first setup has been performed.
 *
 *        One modifies only some specific options like the max. number of
 *        iterations or the relative tolerance
 *
 * \param[in] use_field_id  if false use a name to retrieve the cs_sles_t struc.
 * \param[in] slesp         pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_cvg_param(bool                    use_field_id,
                              const cs_param_sles_t  *slesp);

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the command line option for PETSc
 *
 * \param[in]      use_prefix    need a prefix
 * \param[in]      prefix        optional prefix
 * \param[in]      keyword       command keyword
 * \param[in]      keyval        command value
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_petsc_cmd(bool         use_prefix,
                              const char  *prefix,
                              const char  *keyword,
                              const char  *keyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a KSP structure in PETSc. This is equivalent to set a solver and
 *        its related preconditioner
 *
 * \param[in]      label  label to identify this (part of) system
 * \param[in, out] slesp  pointer to a set of SLES parameters
 * \param[in, out] p_ksp  solver structure for PETSc
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_petsc_ksp(const char       *label,
                              cs_param_sles_t  *slesp,
                              void             *p_ksp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set an AMG preconditioner in PETSc
 *
 * \param[in]      prefix  label to identify this (part of) system
 * \param[in]      slesp   pointer to a set of SLES parameters
 * \param[in, out] p_pc    preconditioner structure for PETsc
 */
/*----------------------------------------------------------------------------*/

void
cs_param_sles_setup_petsc_pc_amg(const char       *prefix,
                                 cs_param_sles_t  *slesp,
                                 void             *p_pc);
#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_SLES_SETUP_H__ */
