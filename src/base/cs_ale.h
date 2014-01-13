#ifndef __CS_ALE_H__
#define __CS_ALE_H__

/*============================================================================
 * Functions associated to ALE formulation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 * Compute cell and face centre of gravity, cell volume.
 *
 * Fortran Interface
 *
 * SUBROUTINE ALGRMA
 * *****************
 *
 *----------------------------------------------------------------------------*/

void
CS_PROCF (algrma, ALGRMA)(void);

/*----------------------------------------------------------------------------
 * Projection on mesh vertices of the displacement (computed on cell center)
 *
 * Fortran Interface
 *
 * SUBROUTINE ALDEPL
 * *****************
 *
 * INTEGER         IFACEL(2,NFAC)  : --> : Interior faces -> cells connectivity
 * INTEGER         IFABOR(NFABOR)  : --> : Border faces -> cells connectivity
 * INTEGER         IPNFAC(NFAC+1)  : --> : Interior faces -> vertices index
 * INTEGER         NODFAC(LNDFAC)  : --> : Interior faces -> vertices list
 * INTEGER         IPNFBR(NFABOR+1): --> : Border faces -> vertices index
 * INTEGER         NODFBR(LNDFBR)  : --> : Border faces -> vertices list
 * DOUBLE PRECISION UMA(NCELET)    : --> : Mesh velocity along X
 * DOUBLE PRECISION VMA(NCELET)    : --> : Mesh velocity along Y
 * DOUBLE PRECISION WMA(NCELET)    : --> : Mesh velocity along Z
 * DOUBLE PRECISION COEFAU(NCELET) : --> : Boundary condition A for UMA
 * DOUBLE PRECISION COEFAV(NCELET) : --> : Boundary condition A pour VMA
 * DOUBLE PRECISION COEFAW(NCELET) : --> : Boundary condition A pour WMA
 * DOUBLE PRECISION COEFBU(NCELET) : --> : Boundary condition B pour UMA
 * DOUBLE PRECISION COEFBV(NCELET) : --> : Boundary condition B pour VMA
 * DOUBLE PRECISION COEFBW(NCELET) : --> : Boundary condition B pour WMA
 * DOUBLE PRECISION DT(NCELET)     : --> : Time step
 * DOUBLE PRECISION DEPROJ(NNOD,3)): <-- : Displacement projected on vertices
 *----------------------------------------------------------------------------*/

void
CS_PROCF (aldepl, ALDEPL)(const cs_int_t    i_face_cells[],
                          const cs_int_t    b_face_cells[],
                          const cs_int_t    i_face_vtx_idx[],
                          const cs_int_t    i_face_vtx_lst[],
                          const cs_int_t    b_face_vtx_idx[],
                          const cs_int_t    b_face_vtx_lst[],
                          cs_real_t        *uma,
                          cs_real_t        *vma,
                          cs_real_t        *wma,
                          cs_real_t        *coefau,
                          cs_real_t        *coefav,
                          cs_real_t        *coefaw,
                          cs_real_t        *coefbu,
                          cs_real_t        *coefbv,
                          cs_real_t        *coefbw,
                          cs_real_t        *dt,
                          cs_real_t        *disp_proj);

/*----------------------------------------------------------------------------
 * Projection on mesh vertices of the displacement (computed on cell center)
 *
 * Fortran Interface
 *
 * subroutine aledis
 * *****************
 *
 * ifacel            : <-- : Interior faces -> cells connectivity
 * ifabor            : <-- : Border faces -> cells connectivity
 * ipnfac            : <-- : Interior faces -> vertices index
 * nodfac            : <-- : Interior faces -> vertices list
 * ipnfbr            : <-- : Border faces -> vertices index
 * nodfbr            : <-- : Border faces -> vertices list
 * ialtyb            : <-- : Type of boundary for ALE
 * meshv             : <-- : Mesh velocity
 * gradm             : <-- : Mesh velocity gradient (du_i/dx_j : gradv[][i][j])
 * claale            : <-- : Boundary conditions A
 * clbale            : <-- : Boundary conditions B
 * dt                : <-- : Time step
 * deproj            : --> : Displacement projected on vertices
 *----------------------------------------------------------------------------*/

void
CS_PROCF (aledis, ALEDIS)(const cs_int_t      i_face_cells[],
                          const cs_int_t      b_face_cells[],
                          const cs_int_t      i_face_vtx_idx[],
                          const cs_int_t      i_face_vtx_lst[],
                          const cs_int_t      b_face_vtx_idx[],
                          const cs_int_t      b_face_vtx_lst[],
                          const cs_int_t      ialtyb[],
                          const cs_real_t    *meshv,
                          const cs_real_33_t  gradm[],
                          const cs_real_t    *claale,
                          const cs_real_t    *clbale,
                          const cs_real_t    *dt,
                          cs_real_t          *disp_proj);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ALE_H__ */

