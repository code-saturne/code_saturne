!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
!> \param[in]     iscal         scalar number
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ifbpcd        index of faces with condensation source terms
!> \param[in]     itypcd        type of condensation source terms for each ivar
!> \param[in]     ncmast        number of cells with metal mass condensation
!> \param[in]     ltmast       index of cells with  metal mass condensation
!> \param[in]     itypst        type of metal mass condensation source terms
!> \param[in]     spcondp       value of the variable associated
!>                              to surface condensation source term
!> \param[in]     gam_s         surface condensation flow rate value
!> \param[in]     svcondp       value of the variable associated
!>                              to metal mass condensation source term
!> \param[in]     gam_ms        metal mass condensation flow rate value
!> \param[in]     fluxv_ms      metal mass condensation heat transfer flux
!> \param[in]     pvara         variable value at time step beginning
!> \param[in,out] tsexp         explicit source term part linear in the variable
!> \param[in,out] tsimp         associated value with \c tsexp
!>                              to be stored in the matrix
!______________________________________________________________________________

subroutine condensation_source_terms &
  (ncelet  , ncel    ,                                   &
   iscal   ,                                             &
   nfbpcd  , ifbpcd  , itypcd ,                          &
   ncmast  , ltmast  , itypst ,                          &
   spcondp , gam_s   ,                                   &
   svcondp , gam_ms  , fluxv_ms ,                        &
   pvara  ,                                              &
   tsexp  , tsimp )

!===============================================================================
! Module files
!===============================================================================

use optcal  , only:iscalt, itherm
use cstphy  , only:voltot
use ppincl  , only:icondb, icondv
use mesh    , only:ifabor, surfbn, volume
use cs_tagms, only:s_metal
!===============================================================================

implicit none

! Variables

integer          ncelet, ncel
integer          nfbpcd, ncmast
integer          iscal
integer          ifbpcd(nfbpcd), itypcd(nfbpcd)
integer          ltmast(ncelet), itypst(ncelet)

double precision spcondp(nfbpcd), gam_s(nfbpcd)
double precision svcondp(ncelet), gam_ms(ncelet)
double precision fluxv_ms(ncelet)
double precision pvara (ncelet)

double precision tsexp(ncelet), tsimp(ncelet)

! Local variables

integer ii, ifac, iel

double precision, allocatable, dimension(:) :: surfbm
!===============================================================================

if (icondb.eq.0) then
 !-----------------------------------------------------------------
 !--- Compute condensation source terms associated to surface zones
 !-----------------------------------------------------------------

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
endif

if (icondv.eq.0) then
  !-----------------------------------------------------------------
  !--- Compute condensation source terms associated to volume zones
  !--- with the metal mass structures modelling
  !-----------------------------------------------------------------
  allocate(surfbm(ncelet))
  surfbm(:) = 0.d0

  !--- Compute the explicit term of the condensation model
  do ii = 1, ncmast
    iel= ltmast(ii)
    surfbm(iel) = s_metal*volume(iel)/voltot

    tsexp(iel)  = tsexp(iel) - surfbm(iel) * gam_ms(iel) * pvara(iel)
    if (itypst(iel).eq.1) then
      if (iscal.eq.iscalt.and.itherm.eq.2) then
        tsexp(iel) = tsexp(iel) + surfbm(iel) *gam_ms(iel) * svcondp(iel) &
                                - fluxv_ms(iel)
      else
        tsexp(iel) = tsexp(iel) + surfbm(iel) *gam_ms(iel) * svcondp(iel)
      endif
    endif
  enddo

  !--- Compute the implicit term of the condensation model
  !--- here we just use actually a explicit form no implicit term
  do ii = 1, ncmast
    iel= ltmast(ii)
    surfbm(iel) = s_metal*volume(iel)/voltot

    if (gam_ms(iel).gt.0.d0) then
      tsimp(iel) = tsimp(iel) + surfbm(ifac) *gam_ms(iel)
    endif
  enddo
  deallocate(surfbm)

endif

return
end subroutine condensation_source_terms
