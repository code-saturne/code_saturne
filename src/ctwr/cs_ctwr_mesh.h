#ifndef __CS_CTWR_MESH_H__
#define __CS_CTWR_MESH_H__

/*============================================================================
 * Specific mesh functions for cooling towers modelling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create nodal coupled mesh.
 * Send vertices's coordinates and connectivity of coupled mesh.
 *
 * Fortran Interface:
 *
 * SUBROUTINE GEOct
 * *****************
 *
 * INTEGER          n_ct     : <-- : number of exchange area
 *----------------------------------------------------------------------------*/

void CS_PROCF(geoct, GEOCT) (void);

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Construction du maillage eau
 *----------------------------------------------------------------------------*/

void cs_ctwr_maille
(
  const cs_mesh_t             *mesh,      /* <-- structure maillage associee  */
  const cs_mesh_quantities_t  *mesh_quantities   /* <-- grandeurs du maillage */
);

/*----------------------------------------------------------------------------
 * Interpolation AIR -> EAU
 *----------------------------------------------------------------------------*/

void
cs_ctwr_adeau(const cs_mesh_t             *mesh,
              const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Interpolation EAU -> AIR
 *----------------------------------------------------------------------------*/

void cs_ctwr_adair (void);

/*----------------------------------------------------------------------------*
 * Chaining of the exchange area                                              *
 *----------------------------------------------------------------------------*/

void
cs_ctwr_stacking (void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CTWR_MESH_H__ */
