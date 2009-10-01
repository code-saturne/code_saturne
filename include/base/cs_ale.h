/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_ALE_H__
#define __CS_ALE_H__

/*============================================================================
 * Functions associated to ALE formulation
 *============================================================================*/

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
 * Public function definitions
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
 * Destroy the associated fvm_interface_set_t structure if necessary
 *
 * Fortran Interface
 *
 * SUBROUTINE LBRALE
 * *****************
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lbrale, LBRALE)(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ALE_H__ */

