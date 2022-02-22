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

!> \file cs_coal_bt2h.f90
!>
!> \brief Convert temperature to enthalpy at boundary for coal combustion.
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


subroutine cs_coal_bt2h(n_faces, face_ids, t_b, h_b)  &
 bind(C, name='cs_coal_bt2h')

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

integer          iel , ilst, ifac , icla , icha , isol , ige
integer          ipcx2c, iii
integer          mode

double precision coefe(ngazem)
double precision f1mc(ncharm), f2mc(ncharm)
double precision x2t, h2, x2h2, hf, xsolid(nsolim), t1, tbl

type(pmapper_double_r1), dimension(:), allocatable :: cvar_xch, cvar_xck
type(pmapper_double_r1), dimension(:), allocatable :: cvar_xnp, cvar_xwt
type(pmapper_double_r1), dimension(:), allocatable :: cvar_f1m, cvar_f2m
type(pmapper_double_r1), dimension(:), allocatable :: cpro_x2, cpro_ym1

!===============================================================================

mode = -1

! Arrays of pointers containing the fields values for each class
! (loop on cells outside loop on classes)

allocate(cvar_xch(nclacp))
allocate(cvar_xck(nclacp))
allocate(cvar_xnp(nclacp))
allocate(cvar_xwt(nclacp))
allocate(cpro_x2(nclacp))
allocate(cpro_ym1(ngazg))

do icla = 1, nclacp
  call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xch(icla)%p)
  call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xck(icla)%p)
  call field_get_val_s(ivarfl(isca(inp(icla))), cvar_xnp(icla)%p)
  if (ippmod(iccoal).eq.1) then
    call field_get_val_s(ivarfl(isca(ixwt(icla))), cvar_xwt(icla)%p)
  endif
  call field_get_val_s(ix2(icla), cpro_x2(icla)%p)
enddo

allocate(cvar_f1m(ncharb))
allocate(cvar_f2m(ncharb))

do icha = 1, ncharb
  call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m(icha)%p)
  call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m(icha)%p)
enddo

do ige = 1, ngazg
  call field_get_val_s(iym1(ige), cpro_ym1(ige)%p)
enddo

! Now loop on faces

do ilst = 1, n_faces

  ifac = face_ids(ilst) + 1
  iel = ifabor(ifac)

  tbl = t_b(ifac)

  x2t  = zero
  x2h2 = zero
  do icla = 1, nclacp
    icha   = ichcor(icla)
    ipcx2c = ix2(icla)
    x2t = x2t + cpro_x2(icla)%p(iel)
    h2 = zero
    do isol = 1, nsolim
      xsolid(isol) = zero
    enddo
    if (cpro_x2(icla)%p(iel).gt.epsicp) then
      xsolid(ich(icha)) = cvar_xch(icla)%p(iel) / cpro_x2(icla)%p(iel)
      xsolid(ick(icha)) = cvar_xck(icla)%p(iel) / cpro_x2(icla)%p(iel)
      xsolid(iash(icha)) = cvar_xnp(icla)%p(iel)*xmash(icla) &
                           / cpro_x2(icla)%p(iel)
      if (ippmod(iccoal).eq.1) then
        xsolid(iwat(icha)) = cvar_xwt(icla)%p(iel) / cpro_x2(icla)%p(iel)
      endif
      iii = icla
      t1 = tbl

      call cs_coal_htconvers2(mode, iii, h2, xsolid, tbl, t1)

    endif
    x2h2 = x2h2 + cpro_x2(icla)%p(iel)*h2
  enddo

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo
  do icha = 1, ncharb
    f1mc(icha) = cvar_f1m(icha)%p(iel) / (1.d0-x2t)
    f2mc(icha) = cvar_f2m(icha)%p(iel) / (1.d0-x2t)
  enddo
  do icha = (ncharb+1), ncharm
    f1mc(icha) = zero
    f2mc(icha) = zero
  enddo
  do ige = 1, ngazg
    coefe(ige) = cpro_ym1(ige)%p(iel)
  enddo

  call cs_coal_htconvers1(mode, hf, coefe, f1mc, f2mc, tbl)

  h_b(ifac) = (1.d0-x2t)*hf+x2h2

enddo

deallocate(cvar_xch, cvar_xck, cvar_xnp)
if (ippmod(iccoal).eq.1) then
  deallocate(cvar_xwt)
endif
deallocate(cvar_f1m, cvar_f2m)
deallocate(cpro_ym1, cpro_x2)

return
end subroutine cs_coal_bt2h
