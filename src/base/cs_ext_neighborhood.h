#ifndef __CS_EXT_NEIGHBOR_H__
#define __CS_EXT_NEIGHBOR_H__

/*============================================================================
 * Fortran interfaces of functions needing a synchronization of the extended
 * neighborhood.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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
 * Define a new "cell -> cells" connectivity for the  extended neighborhood
 * in case of computation of gradient whith the least square algorithm
 * (imrgra = 3).
 * The "cell -> cells" connectivity is clipped by a non-orthogonality
 * criterion.
 *
 * Warning   :  Only cells sharing a vertex or vertices
 *              (not a face => mesh->face_cells) belong to the
 *              "cell -> cells" connectivity.
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
 * Compute filters for dynamic models. This function deals with the standard
 * or extended neighborhood.
 *
 * Fortran Interface :
 *
 * SUBROUTINE CFILTR (VAR, F_VAR, WBUF1, WBUF2)
 * *****************
 *
 * DOUBLE PRECISION(*) var[]      --> array of variables to filter
 * DOUBLE PRECISION(*) f_var[]    --> filtered variable array
 * DOUBLE PRECISION(*) wbuf1[]    --> working buffer
 * DOUBLE PRECISION(*) wbuf2[]    --> working buffer
 *----------------------------------------------------------------------------*/

void
CS_PROCF (cfiltr, CFILTR)(cs_real_t         var[],
                          cs_real_t         f_var[],
                          cs_real_t         wbuf1[],
                          cs_real_t         wbuf2[]);

/*----------------------------------------------------------------------------
 * Create the  "cell -> cells" connectivity
 *
 * parameters:
 *   mesh           <->  pointer to a mesh structure.
 *---------------------------------------------------------------------------*/

void
cs_ext_neighborhood_define(cs_mesh_t   *mesh);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EXT_NEIGHBOR_H__ */
