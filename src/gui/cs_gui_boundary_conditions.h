#ifndef __CS_GUI_BOUNDARY_CONDITION_H__
#define __CS_GUI_BOUNDARY_CONDITION_H__

/*============================================================================
 * Management of the GUI parameters file: boundary conditions
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Boundary conditions treatment
 *
 * Remember: rdoccl[k][j][i] = rcodcl[ k * dim1 *dim2 + j *dim1 + i]
 *
 * Fortran Interface:
 *
 * subroutine uiclim
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclim, UICLIM)(const int  *nozppm,
                               int        *iqimp,
                               int        *icalke,
                               int        *ientat,
                               int        *ientcp,
                               int        *inmoxy,
                               int        *ientox,
                               int        *ientfu,
                               int        *ientgf,
                               int        *ientgb,
                               int        *iprofm,
                               int        *iautom,
                               int        *itypfb,
                               int        *izfppp,
                               int        *icodcl,
                               double     *qimp,
                               double     *qimpat,
                               double     *qimpcp,
                               double     *dh,
                               double     *xintur,
                               double     *timpat,
                               double     *timpcp,
                               double     *tkent,
                               double     *fment,
                               double     *distch,
                               int        *nvar,
                               double     *rcodcl);

/*----------------------------------------------------------------------------
 * Boundary conditions input verification
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICLVE
 * *****************
 *
 * integer          nozppm  <-- max number of boundary conditions zone
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclve, UICLVE)(const int  *nozppm);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define boundary conditions based on setup file.
 *
 * \param[in, out]  bdy   boundaries structure to update
 *                        (if NULL, default to cs_glob_domain->boundaries)
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_define(cs_boundary_t  *bdy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free GUI boundary condition structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_free_memory(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_BOUNDARY_CONDITION_H__ */
