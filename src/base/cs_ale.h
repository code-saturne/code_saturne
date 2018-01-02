#ifndef __CS_ALE_H__
#define __CS_ALE_H__

/*============================================================================
 * Functions associated to ALE formulation
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute cell and face center of gravity, cell volume.
 *
 * Fortran Interface
 *
 * subroutine algrma
 * *****************
 *
 * min_vol           : --> : Minimum cell volume
 * max_vol           : --> : Maximum cell volume
 * tot_vol           : --> : Total mesh volume
 *----------------------------------------------------------------------------*/

void
CS_PROCF (algrma, ALGRMA)(cs_real_t  *min_vol,
                          cs_real_t  *max_vol,
                          cs_real_t  *tot_vol);

/*----------------------------------------------------------------------------
 * Projection on mesh vertices of the displacement (computed on cell center)
 *
 * Fortran Interface
 *
 * subroutine aledis
 * *****************
 *
 * ialtyb            : <-- : Type of boundary for ALE
 * meshv             : <-- : Mesh velocity
 * gradm             : <-- : Mesh velocity gradient (du_i/dx_j : gradv[][i][j])
 * claale            : <-- : Boundary conditions A
 * clbale            : <-- : Boundary conditions B
 * dt                : <-- : Time step
 * disp_proj         : --> : Displacement projected on vertices
 *----------------------------------------------------------------------------*/

void
CS_PROCF (aledis, ALEDIS)(const cs_int_t      ialtyb[],
                          const cs_real_3_t  *meshv,
                          const cs_real_33_t  gradm[],
                          const cs_real_3_t  *claale,
                          const cs_real_33_t *clbale,
                          const cs_real_t    *dt,
                          cs_real_3_t        *disp_proj);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ALE_H__ */

