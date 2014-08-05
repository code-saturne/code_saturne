!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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
!===============================================================================
! Function:
! ---------

!> \file condensation_source_terms.f90
!> \brief Explicit sources terms from sources condensation computation.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     iterns        iteration number on Navier-Stokes
!> \param[in]     ifbpcd        index of cells with condensation source terms
!> \param[in]     itypcd        type of condensation source terms for each ivar
!> \param[in]     pvara         variable value at time step beginning
!> \param[in]     spcondp       value of the variable associated
!>                              to condensation source term
!> \param[in]     gam_s         flow condensation rate value
!> \param[in,out] tsexp         explicit source term part linear in the variable
!> \param[in,out] tsimp         associated value withr \ref tsexp
!>                              to be stored in the matrix
!______________________________________________________________________________

subroutine condensation_source_terms &
  (ncelet , ncel  , nfbpcd  , ifbpcd , itypcd ,    &
   pvara  , spcondp , gam_s  ,                     &
   tsexp  , tsimp )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use mesh, only: ifabor, surfbn

!===============================================================================

implicit none

! Variables

integer          ncelet, ncel , nfbpcd
integer          ifbpcd(nfbpcd), itypcd(nfbpcd)

double precision pvara (ncelet)
double precision spcondp(nfbpcd), gam_s(nfbpcd)

double precision tsexp(ncelet), tsimp(ncelet)

! Local variables

integer ii, ifac, iel

!===============================================================================

!--- Compute the explicit term of the condensation model
do ii = 1, nfbpcd
  ifac= ifbpcd(ii)
  iel = ifabor(ifac)
  tsexp(iel) = tsexp(iel) - surfbn(ifac) * gam_s(ii) * pvara(iel)
  if (itypcd(ii).eq.1) then
    tsexp(iel) = tsexp(iel) + surfbn(ifac) *gam_s(ii) * spcondp(ii)
  endif
enddo

!--- Compute the implicit term of the condensation model
!--- here we just use actually a explicit form no implicit term
do ii = 1, nfbpcd
  ifac= ifbpcd(ii)
  iel = ifabor(ifac)
  if (gam_s(ii).gt.0.d0) then
    tsimp(iel) = tsimp(iel) + surfbn(ifac) *gam_s(ii)
  endif
enddo

return
end subroutine condensation_source_terms
