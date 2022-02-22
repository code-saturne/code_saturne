!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
! Function :
! --------

!> \file cs_fuel_bt2h.f90
!>
!> \brief Convert temperature to enthalpy at boundary for fuel combustion.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n_faces       number of faces in list
!> \param[in]     face_ids      list of boundary faces at which conversion
!>                              is requested (0-based numbering)
!> \param[in]     t_b           temperature at boundary
!> \param[out]    h_b           enthalpy at boundary
!_______________________________________________________________________________


subroutine cs_fuel_bt2h(n_faces, face_ids, t_b, h_b)  &
 bind(C, name='cs_fuel_bt2h')

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe
use ppppar
use ppthch
use cpincl
use cs_fuel_incl
use ppincl
use radiat
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer(c_int), value :: n_faces
integer(c_int), dimension(*), intent(in) :: face_ids
real(kind=c_double), dimension(*), intent(in) :: t_b
real(kind=c_double), dimension(*), intent(out) :: h_b

! Local variables

integer          iel, ilst, ifac, icla, ige
integer          mode

double precision coefe(ngazem)
double precision x2t, h2, x2h2, hf, xsolid(nsolim), tbl
double precision diamgt,masgut,mkgout,mfgout,mkfini,rhofol

type(pmapper_double_r1), dimension(:), allocatable :: cvara_yfol
type(pmapper_double_r1), dimension(:), allocatable :: cpro_ym1
type(pmapper_double_r1), dimension(:), allocatable :: cpro_rom2, cpro_diam2

!===============================================================================

mode = -1

! Arrays of pointers containing the fields values for each class
! (loop on cells outside loop on classes)

allocate(cvara_yfol(nclafu))
allocate(cpro_rom2(nclafu))
allocate(cpro_diam2(nclafu))

allocate(cpro_ym1(ngaze))

do icla = 1, nclafu
  call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol(icla)%p)
  call field_get_val_prev_s(irom2(icla), cpro_rom2(icla)%p)
  call field_get_val_prev_s(idiam2(icla), cpro_diam2(icla)%p)
enddo
do ige = 1, ngaze
  call field_get_val_s(iym1(ige), cpro_ym1(ige)%p)
enddo

! Now loop on faces

do ilst = 1, n_faces

  ifac = face_ids(ilst) + 1
  iel = ifabor(ifac)

  tbl = t_b(ifac)

  ! Coal combustion

  x2t  = zero
  x2h2 = zero

  do icla = 1, nclafu

    x2t = x2t+cvara_yfol(icla)%p(iel)

    mkfini = rho0fl*pi/6.d0*dinikf(icla)**3
    rhofol = cpro_rom2(icla)%p(iel)
    diamgt = cpro_diam2(icla)%p(iel)
    masgut = rhofol*pi/6.d0*diamgt**3
    if (diamgt.le.dinikf(icla)) then
      mkgout = masgut
    else
      mkgout = mkfini
    endif
    mfgout = masgut - mkgout
    xsolid(1) = 1.d0-fkc
    xsolid(2) = fkc
    if (masgut.gt.epzero) then
      xsolid(1) = mfgout / masgut
      xsolid(2) = mkgout / masgut
    endif
    xsolid(1) = min(1.d0,max(0.d0,xsolid(1)))
    xsolid(2) = min(1.d0,max(0.d0,xsolid(2)))

    call cs_fuel_htconvers2(mode, h2, xsolid, tbl)

    x2h2 = x2h2 + cvara_yfol(icla)%p(iel)*h2

  enddo

  do ige = 1,ngaze
    coefe(ige) = cpro_ym1(ige)%p(iel)
  enddo
  call cs_fuel_htconvers1(mode, hf, coefe, tbl)
  h_b(ifac) = (1.d0-x2t)*hf+x2h2

enddo

deallocate(cvara_yfol)
deallocate(cpro_rom2, cpro_diam2, cpro_ym1)

return
end subroutine cs_fuel_bt2h
