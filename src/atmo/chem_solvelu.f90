!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!> \file chem_solvelu.f90
!> \brief Solver of AX=B with LU decomposition of A for atmospheric chemical
!>        systems
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     kindlu        flag for LU decomposition
!> \param[in]     dla           Matrix A in AX=B
!> \param[in,out] dlalu         Matrix A in AX=B after LU decomposition
!> \param[out]    dlx           Vector of unknowns X in AX=B
!> \param[in]     dlb           Vector B in AX=B
!______________________________________________________________________________

subroutine cs_solvlin (kindlu,dla,dlalu,dlx,dlb)

!===============================================================================
! Module files
!===============================================================================

use atchem

implicit none

procedure() :: lu_decompose_1, lu_decompose_2, lu_decompose_3
procedure() :: lu_solve_1, lu_solve_2, lu_solve_3

! Arguments

integer kindlu
double precision dla(nespg,nespg)
double precision dlalu(nespg,nespg)
double precision dlx(nespg), dlb(nespg)

! Local variables

integer ji, jj

!------------------------------------------------------------------------
!*    0. Setup

do ji = 1, nespg
  dlx(ji) = dlb(ji)
enddo

!------------------------------------------------------------------------
!*    1. Compute DLx

if (kindlu .eq. 0) then
  do jj = 1, nespg
    do ji = 1, nespg
      dlalu(ji,jj) = dla(ji,jj)
    enddo
  enddo

  if (ichemistry.eq.1) then
    call lu_decompose_1(nespg,dlalu)
    call lu_solve_1(nespg,dlalu,dlx)
  else if (ichemistry.eq.2) then
    call lu_decompose_2(nespg,dlalu)
    call lu_solve_2(nespg,dlalu,dlx)
  else if (ichemistry.eq.3) then
    call lu_decompose_3(nespg,dlalu)
    call lu_solve_3(nespg,dlalu,dlx)
  else if (ichemistry.eq.4) then
    call ssh_lu_decompose(nespg,dlalu)
    call ssh_lu_solve(nespg,dlalu,dlx)
  endif
else
  if (ichemistry.eq.1) then
    call lu_solve_1(nespg,dlalu,dlx)
  else if (ichemistry.eq.2) then
    call lu_solve_2(nespg,dlalu,dlx)
  else if (ichemistry.eq.3) then
    call lu_solve_3(nespg,dlalu,dlx)
  else if (ichemistry.eq.4) then
    call ssh_lu_solve(nespg,dlalu,dlx)
  endif
endif

end subroutine cs_solvlin
