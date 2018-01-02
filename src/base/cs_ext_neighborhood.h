#ifndef __CS_EXT_NEIGHBOR_H__
#define __CS_EXT_NEIGHBOR_H__

/*============================================================================
 * Fortran interfaces of functions needing a synchronization of the extended
 * neighborhood.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Reduce the "cell -> cells" connectivity for the  extended neighborhood
 * using a non-orthogonality criterion.
 *
 * Note: Only cells sharing only a vertex or vertices (not a face)
 *       belong to the "cell -> cells" connectivity.
 *
 * Fortran Interface :
 *
 * SUBROUTINE REDVSE
 * *****************
 *    & ( ANOMAX )
 *
 * parameters:
 *   anomax  -->  non-orthogonality angle (rad) above which cells
 *                are selected for the extended neighborhood
 *----------------------------------------------------------------------------*/

void
CS_PROCF (redvse, REDVSE) (const cs_real_t  *anomax);

/*----------------------------------------------------------------------------
 * Reduce the "cell -> cells" connectivity for the  extended neighborhood
 * using a non-orthogonality criterion.
 *
 * Note: Only cells sharing only a vertex or vertices (not a face)
 *       belong to the "cell -> cells" connectivity.
 *
 * parameters:
 *   mesh            <-> pointer to mesh structure
 *   mesh_quantities <-> associated mesh quantities
 *   non_ortho_max   <-- non-orthogonality angle (rad) above which cells
 *                       are selected for the extended neighborhood
 *----------------------------------------------------------------------------*/

void
cs_ext_neighborhood_reduce(cs_mesh_t             *mesh,
                           cs_mesh_quantities_t  *mesh_quantities,
                           double                 non_ortho_max);

/*----------------------------------------------------------------------------
 * Create the  "cell -> cells" connectivity.
 *
 * parameters:
 *   mesh <-> pointer to a mesh structure
 *---------------------------------------------------------------------------*/

void
cs_ext_neighborhood_define(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EXT_NEIGHBOR_H__ */
