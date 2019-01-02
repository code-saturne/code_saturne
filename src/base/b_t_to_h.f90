!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file b_t_to_h.f90
!>
!> \brief Convert temperature to enthalpy at boundary
!>
!> This handles both user and model enthalpy conversions, so can be used
!> safely whenever conversion is needed.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nlst          number of faces in list
!> \param[in]     lstfac        list of boundary faces at which conversion
!>                              is requested
!> \param[in]     t_b           temperature at boundary
!> \param[out]    h_b           enthalpy at boundary
!_______________________________________________________________________________


subroutine b_t_to_h            &
 ( nlst, lstfac, t_b , h_b )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use radiat
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nlst
integer          lstfac(nlst)

double precision, dimension(nfabor),intent(in) :: t_b
double precision, dimension(nfabor), intent(out), target :: h_b

! Local variables

integer          iel , ilst, ifac , icla , icha , isol , ige
integer          ipcx2c , igg, iii
integer          iesp, mode

double precision coefg(ngazgm), coefe(ngazem)
double precision f1mc(ncharm), f2mc(ncharm)
double precision x2t, h2, x2h2, hf, xsolid(nsolim), t1, tbl
double precision ym(ngazgm)
double precision diamgt,masgut,mkgout,mfgout,mkfini,rhofol
character(len=80) :: f_name

double precision, dimension(:), pointer :: bym1, bym2, bym3

type(pmapper_double_r1), dimension(:), allocatable :: cvar_xch, cvar_xck
type(pmapper_double_r1), dimension(:), allocatable :: cvar_xnp, cvar_xwt
type(pmapper_double_r1), dimension(:), allocatable :: cvar_f1m, cvar_f2m
type(pmapper_double_r1), dimension(:), allocatable :: cvara_yfol
type(pmapper_double_r1), dimension(:), allocatable :: cvar_ycoel
type(pmapper_double_r1), dimension(:), allocatable :: cpro_x2, cpro_ym1
type(pmapper_double_r1), dimension(:), allocatable :: cpro_rom2, cpro_diam2

!===============================================================================

mode = -1

! Non-specific physics

if (ippmod(iphpar).le.1) then

  do ilst = 1, nlst
    ifac = lstfac(ilst)
    tbl = t_b(ifac)
    call usthht(mode, h_b(ifac), tbl)
  enddo

  return

endif

! Mappings for specific physics

call field_get_val_s(ibym(1), bym1)
call field_get_val_s(ibym(2), bym2)
call field_get_val_s(ibym(3), bym3)

! Arrays of pointers containing the fields values for each class
! (loop on cells outside loop on classes)

if (ippmod(iccoal).ge.0) then

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

else if (ippmod(icfuel).ge.0) then

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

else if (ippmod(ielarc).ge.1) then

  allocate(cvar_ycoel(ngazg-1))

  do iesp = 1, ngazg-1
    write(f_name,'(a13,i2.2)') 'esl_fraction_',iesp
    call field_get_val_prev_s_by_name(trim(f_name), cvar_ycoel(iesp)%p)
  enddo

endif

! Now loop on faces

do ilst = 1, nlst

  ifac = lstfac(ilst)
  iel = ifabor(ifac)

  tbl = t_b(ifac)

  ! Gas combustion: premix or diffusion flame

  if (ippmod(icoebu).ge.0 .or. ippmod(icod3p).ge.0) then

    do igg = 1, ngazgm
      coefg(igg) = zero
    enddo
    coefg(1) = bym1(ifac)
    coefg(2) = bym2(ifac)
    coefg(3) = bym3(ifac)
    call cothht(mode, ngazg, ngazgm, coefg, npo, npot, th, ehgazg,   &
                h_b(ifac), tbl)

  ! Pulverized coal combustion

  else if (ippmod(iccoal).ge.0) then

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

  ! Fuel combustion

  else if (ippmod(icfuel).ge.0) then

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

  ! Electric arcs

  else if (ippmod(ieljou).ge.1) then

    call usthht(mode, h_b(ifac), tbl)

  else if (ippmod(ielarc).ge.1) then

    if (ngazg .eq. 1) then
      ym(1) = 1.d0
      call elthht(mode, ym, h_b(ifac), tbl)
    else
      ym(ngazg) = 1.d0
      do iesp = 1, ngazg-1
        ym(iesp) = cvar_ycoel(iesp)%p(iel)
        ym(ngazg) = ym(ngazg) - ym(iesp)
      enddo
      call elthht(mode, ym, h_b(ifac), tbl)
    endif

  endif

enddo

if (ippmod(iccoal).ge.0) then
  deallocate(cvar_xch, cvar_xck, cvar_xnp)
  if (ippmod(iccoal).eq.1) then
    deallocate(cvar_xwt)
  endif
  deallocate(cvar_f1m, cvar_f2m)
  deallocate(cpro_ym1, cpro_x2)
else if (ippmod(icfuel).ge.0) then
  deallocate(cvara_yfol)
  deallocate(cpro_rom2, cpro_diam2, cpro_ym1)
else if (ippmod(ielarc).ge.1) then
  deallocate(cvar_ycoel)
endif

return
end subroutine b_t_to_h
