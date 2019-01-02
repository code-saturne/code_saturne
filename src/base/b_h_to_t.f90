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
!> \brief Convert enthalpy to temperature at boundary
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
!> \param[in]     h_b           enthalpy at boundary
!> \param[in,out] t_b           temperature at boundary
!_______________________________________________________________________________


subroutine b_h_to_t            &
 ( h_b , t_b )

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
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

double precision, dimension(nfabor),intent(in) :: h_b
double precision, dimension(nfabor), intent(out), target :: t_b

! Local variables

integer          iel , ifac
integer          igg
integer          iesp, mode

double precision coefg(ngazgm)
double precision hbl
double precision ym(ngazgm)
character(len=80) :: f_name

double precision, dimension(:), pointer :: bym1, bym2, bym3

type(pmapper_double_r1), dimension(:), allocatable :: cvar_ycoel

!===============================================================================

mode = 1

! Non-specific physics

if (ippmod(iphpar).le.1) then

  do ifac = 1, nfabor
    hbl = h_b(ifac)
    call usthht(mode, hbl, t_b(ifac))
  enddo

! Gas combustion: premix or diffusion flame

else if (ippmod(icoebu).ge.0 .or. ippmod(icod3p).ge.0) then

  call field_get_val_s(ibym(1), bym1)
  call field_get_val_s(ibym(2), bym2)
  call field_get_val_s(ibym(3), bym3)

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    hbl = h_b(ifac)
    do igg = 1, ngazgm
      coefg(igg) = zero
    enddo
    coefg(1) = bym1(ifac)
    coefg(2) = bym2(ifac)
    coefg(3) = bym3(ifac)
    call cothht(mode, ngazg, ngazgm, coefg, npo, npot, th, ehgazg,   &
                hbl, t_b(ifac))
  enddo

else if (ippmod(iccoal).ge.0) then

  ! Use mixture temperature
  call cs_coal_thfieldconv1(MESH_LOCATION_BOUNDARY_FACES, h_b, t_b)

else if (ippmod(icfuel).ge.0) then

  ! Use mixture temperature
  call cs_fuel_thfieldconv1(MESH_LOCATION_BOUNDARY_FACES, h_b, t_b)

! Electric arcs

else if (ippmod(ieljou).ge.1) then

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    hbl = h_b(ifac)
    call usthht(mode, hbl, t_b(ifac))
  enddo

else if (ippmod(ielarc).ge.1) then

  if (ngazge .gt. 1) then
    allocate(cvar_ycoel(ngazge-1))
    do iesp = 1, ngazge-1
      write(f_name,'(a13,i2.2)') 'esl_fraction_',iesp
      call field_get_val_prev_s_by_name(trim(f_name), cvar_ycoel(iesp)%p)
    enddo
  endif

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    hbl = h_b(ifac)
    if (ngazge .eq. 1) then
      ym(1) = 1.d0
      call elthht(mode, ym, hbl, t_b(ifac))
    else
      ym(ngazge) = 1.d0
      do iesp = 1, ngazge-1
        ym(iesp) = cvar_ycoel(iesp)%p(iel)
        ym(ngazge) = ym(ngazge) - ym(iesp)
      enddo
      call elthht(mode, ym, hbl, t_b(ifac))
    endif
  enddo

  if (ngazge .gt. 1) then
    deallocate(cvar_ycoel)
  endif

endif

return
end subroutine b_h_to_t
