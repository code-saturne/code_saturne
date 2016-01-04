!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
!> \file cs_fuel_thfieldconv2.f90
!> \brief Calculation of the particles temperature
!>  Fonction with the fuel enthalpy and concentrations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in,out] propce        physical properties at cell centers
!______________________________________________________________________________!

subroutine cs_fuel_thfieldconv2 &
 ( ncelet , ncel   ,            &
   propce )

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
use ppincl
use field
use cs_fuel_incl

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel
double precision propce(ncelet,*)

! Local variables

integer          iel    , icla
integer          ipcte1 , ipcte2
integer          mode

double precision eh2
double precision xsolid(2), mkfini , diamgt
double precision masgut  , mfgout , mkgout , rhofol

double precision, dimension(:), pointer :: cvar_yfolcl, cvar_h2cl

!===============================================================================
! 1. Preliminary calculations
!===============================================================================
! --- Initialization of T2 as T1

ipcte1 = ipproc(itemp1)
do icla = 1, nclafu
  ipcte2 = ipproc(itemp2(icla))
  do iel = 1, ncel
    propce(iel,ipcte2) = propce(iel,ipcte1)
  enddo
enddo

!===============================================================================
! 2. Calculation of the particles temperature
!===============================================================================

do icla = 1, nclafu

  call field_get_val_s(ivarfl(isca(iyfol(icla))), cvar_yfolcl)
  call field_get_val_s(ivarfl(isca(ih2(icla))), cvar_h2cl)
  ipcte2 = ipproc(itemp2(icla))

  mkfini = rho0fl*pi/6.d0*dinikf(icla)**3

  do iel = 1, ncel

    rhofol = propce(iel,ipproc(irom2 (icla)))
    diamgt = propce(iel,ipproc(idiam2(icla)))
    masgut = rho0fl*pi/6.d0*(diamgt**3.d0)
    if (diamgt.le.dinikf(icla)) then
      mkgout = masgut
    else
      mkgout = mkfini
    endif
    mfgout = masgut - mkgout
    xsolid(1) = 1.d0-fkc
    xsolid(2) = fkc
    if(masgut.gt.zero) then
      xsolid(1) = mfgout / masgut
      xsolid(2) = mkgout / masgut
    endif
    xsolid(1) = min(1.d0,max(0.d0,xsolid(1)))
    xsolid(2) = min(1.d0,max(0.d0,xsolid(2)))

    if ( cvar_yfolcl(iel) .gt. (3.d3*epsifl) ) then

      eh2 =  cvar_h2cl(iel)/cvar_yfolcl(iel)

      mode = 1
      call cs_fuel_htconvers2(mode, eh2, xsolid, propce(iel,ipcte2))

    endif

  enddo
!
enddo

!----
! End
!----
return
end subroutine
