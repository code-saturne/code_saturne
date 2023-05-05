!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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
!> \file cs_tagmri.f90
!>
!> \brief The 1D thermal model to compute the temperature to impose
!> at the cold wall. This one is used by the COPAIN model to estimate
!> the heat flux at the wall where the condensation occurs.
!>
!> This subroutine is used to compute at each face the
!> \f$ \mbox{tmur} \f$ at cold wall.
!-------------------------------------------------------------------------------

subroutine cs_tagmri()  &
  bind(C, name='cs_f_tagmri')

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use cstphy
use entsor
use pointe
use field
use mesh, only:ifabor
use cs_nz_condensation, only: nfbpcd, izzftcd, iztag1d, ztpar, ifbpcd
use cs_nz_tagmr, only: ztmur
use cs_c_bindings

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

! Arguments


! Local variables


integer          ii, iz, ivar
integer          ifac, iel
integer          icldef
double precision temper, enthal

double precision, dimension(:), pointer :: cpro_cp

integer, pointer, dimension(:,:) :: icodcl
double precision, pointer, dimension(:,:,:) :: rcodcl

!===============================================================================

! Sans specification, une face couplee est une face de type paroi

icldef = 5

ivar = isca(iscalt)

call field_build_bc_codes_all(icodcl, rcodcl) ! Get map

do ii = 1, nfbpcd

  ifac = ifbpcd(ii) + 1
  iz = izzftcd(ii) + 1
  if(iztag1d(iz).eq.1) then
    icodcl(ifac,ivar)   = 1
    rcodcl(ifac,ivar,1) = ztmur(ii,1)
    rcodcl(ifac,ivar,2) = rinfin
    rcodcl(ifac,ivar,3) = 0.d0
  else
    icodcl(ifac,ivar)   = 1
    rcodcl(ifac,ivar,1) = ztpar(iz)
    rcodcl(ifac,ivar,2) = rinfin
    rcodcl(ifac,ivar,3) = 0.d0
  endif
enddo

! Conversion eventuelle temperature -> enthalpie

if (itherm.eq.2) then

  ! --- Specific heat
  if (icp.ge.0) then
    call field_get_val_s(icp, cpro_cp)
    ! --- Stop if Cp is not variable
  else
    write(nfecra,1000) icp
    call csexit (1)
  endif

  do ii = 1, nfbpcd

    ifac = ifbpcd(ii) + 1
    iel = ifabor(ifac)

    temper = rcodcl(ifac,ivar,1)

    enthal = (temper + tkelvi) * cpro_cp(iel)

    rcodcl(ifac,ivar,1) = enthal

  enddo

endif

!--------
! Formats
!--------

 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@      usipsu specifies that the specific heat is uniform',/,    &
'@        icp = ',i10   ,' while',/,                              &
'@      cs_user_physical_properties prescribes a variable specific heat.',/,&
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or cs_user_physical_properties.',/,           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

!----
! End
!----

return
end subroutine
