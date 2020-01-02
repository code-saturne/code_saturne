!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
! --------
!> \file cs_fuel_thfieldconv1.f90
!>
!> \brief Calculation of the gas temperature
!>        Function with the gas enthalpy and concentrations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     location_id   mesh location id (cells or boundary faces)
!> \param[in]     eh            gas enthalpy
!>                              (\f$ j . kg \f$ of gaseous mixture)
!> \param[in,out] tp            gas temperature (in kelvin)
!______________________________________________________________________________!

subroutine cs_fuel_thfieldconv1 &
 ( location_id     ,                                              &
   eh     ,                                                       &
   tp     )

!==============================================================================
! Module files
!==============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use ppcpfu
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          location_id

double precision eh(*)
double precision tp(*)

! Local variables

integer          ii, icel, ifac
double precision eh0, eh1

double precision, dimension(:), pointer :: fuel1, fuel2, fuel3, fuel4, fuel5
double precision, dimension(:), pointer :: fuel6, fuel7, oxyd
double precision, dimension(:), pointer :: prod1, prod2, prod3 , xiner

!===============================================================================

call field_get_val_s(iym1(ifo0), fuel1)
call field_get_val_s(iym1(ifov), fuel2)
call field_get_val_s(iym1(ico ), fuel3)
call field_get_val_s(iym1(ih2s), fuel4)
call field_get_val_s(iym1(ihy ), fuel5)
call field_get_val_s(iym1(ihcn), fuel6)
call field_get_val_s(iym1(inh3), fuel7)
call field_get_val_s(iym1(io2 ), oxyd)
call field_get_val_s(iym1(ico2), prod1)
call field_get_val_s(iym1(ih2o), prod2)
call field_get_val_s(iym1(iso2), prod3)
call field_get_val_s(iym1(in2 ), xiner)

if (location_id .eq. MESH_LOCATION_CELLS) then

  ii = npo-1
  do icel = 1, ncel

    ! --- Eventual clipping of temperature at TH(NPO) if EH > EH1

    eh1 = fuel1(icel)*ehgaze(ifo0,ii+1)                       &
         +fuel2(icel)*ehgaze(ifov,ii+1)                       &
         +fuel3(icel)*ehgaze(ico ,ii+1)                       &
         +fuel4(icel)*ehgaze(ih2s,ii+1)                       &
         +fuel5(icel)*ehgaze(ihy ,ii+1)                       &
         +fuel6(icel)*ehgaze(ihcn,ii+1)                       &
         +fuel7(icel)*ehgaze(inh3,ii+1)                       &
         +oxyd(icel) *ehgaze(io2 ,ii+1)                       &
         +prod1(icel)*ehgaze(ico2,ii+1)                       &
         +prod2(icel)*ehgaze(ih2o,ii+1)                       &
         +prod3(icel)*ehgaze(iso2,ii+1)                       &
         +xiner(icel)*ehgaze(in2 ,ii+1)
    if (eh(icel) .ge. eh1) tp(icel) = th(ii+1)

  enddo

  ii = 1
  do icel = 1, ncel

    ! --- Eventual clipping of temperature at TH(1) if EH < EH0

    eh0 = fuel1(icel)*ehgaze(ifo0,ii)                         &
         +fuel2(icel)*ehgaze(ifov,ii)                         &
         +fuel3(icel)*ehgaze(ico ,ii)                         &
         +fuel4(icel)*ehgaze(ih2s,ii)                         &
         +fuel5(icel)*ehgaze(ihy ,ii)                         &
         +fuel6(icel)*ehgaze(ihcn,ii)                         &
         +fuel7(icel)*ehgaze(inh3,ii)                         &
         +oxyd(icel) *ehgaze(io2 ,ii)                         &
         +prod1(icel)*ehgaze(ico2,ii)                         &
         +prod2(icel)*ehgaze(ih2o,ii)                         &
         +prod3(icel)*ehgaze(iso2,ii)                         &
         +xiner(icel)*ehgaze(in2 ,ii)

    if (eh(icel) .le. eh0) then
      tp(icel) = th(1)
    endif

  enddo

  do ii = 1, npo-1

    do icel = 1, ncel

      eh0 = fuel1(icel)*ehgaze(ifo0,ii)                       &
           +fuel2(icel)*ehgaze(ifov,ii)                       &
           +fuel3(icel)*ehgaze(ico ,ii)                       &
           +fuel4(icel)*ehgaze(ih2s,ii)                       &
           +fuel5(icel)*ehgaze(ihy ,ii)                       &
           +fuel6(icel)*ehgaze(ihcn,ii)                       &
           +fuel7(icel)*ehgaze(inh3,ii)                       &
           +oxyd(icel) *ehgaze(io2 ,ii)                       &
           +prod1(icel)*ehgaze(ico2,ii)                       &
           +prod2(icel)*ehgaze(ih2o,ii)                       &
           +prod3(icel)*ehgaze(iso2,ii)                       &
           +xiner(icel)*ehgaze(in2 ,ii)

      eh1 = fuel1(icel)*ehgaze(ifo0,ii+1)                     &
           +fuel2(icel)*ehgaze(ifov,ii+1)                     &
           +fuel3(icel)*ehgaze(ico ,ii+1)                     &
           +fuel4(icel)*ehgaze(ih2s,ii+1)                     &
           +fuel5(icel)*ehgaze(ihy ,ii+1)                     &
           +fuel6(icel)*ehgaze(ihcn,ii+1)                     &
           +fuel7(icel)*ehgaze(inh3,ii+1)                     &
           +oxyd(icel) *ehgaze(io2 ,ii+1)                     &
           +prod1(icel)*ehgaze(ico2,ii+1)                     &
           +prod2(icel)*ehgaze(ih2o,ii+1)                     &
           +prod3(icel)*ehgaze(iso2,ii+1)                     &
           +xiner(icel)*ehgaze(in2 ,ii+1)

      if (eh(icel) .ge. eh0 .and. eh(icel) .le. eh1) then
        tp(icel) = th(ii) + (eh(icel)-eh0) * (th(ii+1)-th(ii)) / (eh1-eh0)
      endif

    enddo

  enddo

else if (location_id .eq. MESH_LOCATION_BOUNDARY_FACES) then

  ii = npo-1
  do ifac = 1, nfabor

    icel = ifabor(ifac)

    ! --- Eventual clipping of temperature at TH(NPO) if EH > EH1

    eh1 = fuel1(icel)*ehgaze(ifo0,ii+1)                       &
         +fuel2(icel)*ehgaze(ifov,ii+1)                       &
         +fuel3(icel)*ehgaze(ico ,ii+1)                       &
         +fuel4(icel)*ehgaze(ih2s,ii+1)                       &
         +fuel5(icel)*ehgaze(ihy ,ii+1)                       &
         +fuel6(icel)*ehgaze(ihcn,ii+1)                       &
         +fuel7(icel)*ehgaze(inh3,ii+1)                       &
         +oxyd(icel) *ehgaze(io2 ,ii+1)                       &
         +prod1(icel)*ehgaze(ico2,ii+1)                       &
         +prod2(icel)*ehgaze(ih2o,ii+1)                       &
         +prod3(icel)*ehgaze(iso2,ii+1)                       &
         +xiner(icel)*ehgaze(in2 ,ii+1)
    if (eh(ifac) .ge. eh1) tp(ifac) = th(ii+1)

  enddo

  ii = 1
  do ifac = 1, nfabor

    icel = ifabor(ifac)

    ! --- Eventual clipping of temperature at TH(1) if EH < EH0

    eh0 = fuel1(icel)*ehgaze(ifo0,ii)                         &
         +fuel2(icel)*ehgaze(ifov,ii)                         &
         +fuel3(icel)*ehgaze(ico ,ii)                         &
         +fuel4(icel)*ehgaze(ih2s,ii)                         &
         +fuel5(icel)*ehgaze(ihy ,ii)                         &
         +fuel6(icel)*ehgaze(ihcn,ii)                         &
         +fuel7(icel)*ehgaze(inh3,ii)                         &
         +oxyd(icel) *ehgaze(io2 ,ii)                         &
         +prod1(icel)*ehgaze(ico2,ii)                         &
         +prod2(icel)*ehgaze(ih2o,ii)                         &
         +prod3(icel)*ehgaze(iso2,ii)                         &
         +xiner(icel)*ehgaze(in2 ,ii)

    if (eh(ifac) .le. eh0) then
      tp(ifac) = th(1)
    endif

  enddo

  do ii = 1, npo-1

    do ifac = 1, nfabor

      icel = ifabor(ifac)

      eh0 = fuel1(icel)*ehgaze(ifo0,ii)                       &
           +fuel2(icel)*ehgaze(ifov,ii)                       &
           +fuel3(icel)*ehgaze(ico ,ii)                       &
           +fuel4(icel)*ehgaze(ih2s,ii)                       &
           +fuel5(icel)*ehgaze(ihy ,ii)                       &
           +fuel6(icel)*ehgaze(ihcn,ii)                       &
           +fuel7(icel)*ehgaze(inh3,ii)                       &
           +oxyd(icel) *ehgaze(io2 ,ii)                       &
           +prod1(icel)*ehgaze(ico2,ii)                       &
           +prod2(icel)*ehgaze(ih2o,ii)                       &
           +prod3(icel)*ehgaze(iso2,ii)                       &
           +xiner(icel)*ehgaze(in2 ,ii)

      eh1 = fuel1(icel)*ehgaze(ifo0,ii+1)                     &
           +fuel2(icel)*ehgaze(ifov,ii+1)                     &
           +fuel3(icel)*ehgaze(ico ,ii+1)                     &
           +fuel4(icel)*ehgaze(ih2s,ii+1)                     &
           +fuel5(icel)*ehgaze(ihy ,ii+1)                     &
           +fuel6(icel)*ehgaze(ihcn,ii+1)                     &
           +fuel7(icel)*ehgaze(inh3,ii+1)                     &
           +oxyd(icel) *ehgaze(io2 ,ii+1)                     &
           +prod1(icel)*ehgaze(ico2,ii+1)                     &
           +prod2(icel)*ehgaze(ih2o,ii+1)                     &
           +prod3(icel)*ehgaze(iso2,ii+1)                     &
           +xiner(icel)*ehgaze(in2 ,ii+1)

      if (eh(ifac) .ge. eh0 .and. eh(ifac) .le. eh1) then
        tp(ifac) = th(ii) + (eh(ifac)-eh0) * (th(ii+1)-th(ii)) / (eh1-eh0)
      endif

    enddo

  enddo

endif

!----
! End
!----

return
end subroutine
