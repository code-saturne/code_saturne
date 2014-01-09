!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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
! Purpose:
! -------

!> \file cs_user_extra_operations-print_statistical_moment.f90
!> This is an example of cs_user_extra_operations.f90 which
!> prints statistical moment

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     , rtpa   , rtp    , propce )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use lagpar
use lagran
use lagdim
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)

! Local variables

!< [loc_var_dec]
integer          iel
integer          imom   , ipcmom , idtcm
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!===============================================================================
! Example: print first calculated statistical moment
!===============================================================================

!< [example_1]
if (nbmomt.gt.0) then

  imom = 1 ! Moment number

  ! Position in 'propce' of the array of temporal accumulation for moments,
  ! propce(iel,ipcmom)

  ipcmom = ipproc(icmome(imom))

  ! The temporal accumulation for moments must be divided by the accumulated
  ! time, which id an array of size ncel or a single real number:
  ! - array of size ncel if idtmom(imom) > 0 : propce(iel, idtcm)
  ! - or simple real     if idtmom(imom) < 0 : dtcmom(idtcm)

  if (idtmom(imom).gt.0) then
    idtcm = ipproc(icdtmo(idtmom(imom)))
    do iel = 1, ncel
      write(nfecra, 4000)  &
           iel, propce(iel, ipcmom)/max(propce(iel, idtcm), epzero)
    enddo
  elseif (idtmom(imom).lt.0) then
    idtcm = -idtmom(imom)
    do iel = 1, ncel
      write(nfecra,4000)  &
           iel, propce(iel, ipcmom)/max(dtcmom(idtcm), epzero)
    enddo
  endif

endif

!< [example_1]
4000 format(' Cell ',i10,'   First moment ',e14.5)

return
end subroutine cs_f_user_extra_operations
