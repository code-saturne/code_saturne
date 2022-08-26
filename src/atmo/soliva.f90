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
!> \file soliva.f90
!> \brief Atmospheric soil module - soil variables initialisation

!> \brief  Initialize soil model variables
!>  NB : data structures are define in module atsoil.f90
!-------------------------------------------------------------------------------

subroutine soliva

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use atincl
use atsoil
use field
use ctincl, only: cp_a, cp_v
use mesh

implicit none

!===============================================================================

! Local variables

integer          isol,iphysi
double precision rscp
double precision esaini,qsaini,huini,psini
double precision cpvcpa

double precision, pointer, dimension(:)   :: bvar_temp_sol
double precision, pointer, dimension(:)   :: bvar_tempp
double precision, pointer, dimension(:)   :: bvar_total_water
double precision, pointer, dimension(:)   :: bvar_w1
double precision, pointer, dimension(:)   :: bvar_w2

!===============================================================================

!     1 - initialisation du tableau solva
!     -----------------------------------

cpvcpa = cp_v / cp_a

call field_get_val_s_by_name("soil_temperature", bvar_temp_sol)
call field_get_val_s_by_name("soil_pot_temperature", bvar_tempp)
call field_get_val_s_by_name("soil_total_water", bvar_total_water)
call field_get_val_s_by_name("soil_w1", bvar_w1)
call field_get_val_s_by_name("soil_w2", bvar_w2)

!  initialisation de t et qv en z0

if (qvsini.gt.1.d0) then
  !  si qvsini>1 qvsini represente l'humidite relative en %
  esaini = 610.78d0*exp(17.2694d0*tsini/                           &
    (tsini + tkelvi - 35.86d0))
  qsaini = esaini/(rvsra*p0 + esaini*(1.d0 - rvsra))
  qvsini = qvsini*qsaini/100.d0
endif

!==============================================================

if (ippmod(iatmos).eq.1) iphysi = 0
if (ippmod(iatmos).eq.2) iphysi = 3

do isol = 1, nfmodsol

  psini = p0
  if (iphysi.gt.0) then
    rscp = (rair/cp0)*(1.d0 + (rvsra - cpvcpa)*qvsini)
  else
    rscp = (rair/cp0)
  endif

  bvar_temp_sol(isol) = tsini

  bvar_tempp(isol) = (tsini+tkelvi)*((ps/psini)**rscp)

  if (iphysi.eq.0) then
    bvar_total_water(isol) = 0.d0
    bvar_w1(isol) = 0.d0
    bvar_w2(isol) = 0.d0
  else
    bvar_total_water(isol) = qvsini

    bvar_w1(isol) = w1ini
    bvar_w2(isol) = w2ini
    ! If not initialized, we compute an approximation of the initial
    ! water content of the top reservoir
    if (bvar_w1(isol).lt.1.d-20) then
      esaini = 610.78d0*exp(17.2694d0*tsini/                           &
        (tsini + tkelvi - 35.86d0))
      qsaini = esaini/(rvsra*psini + esaini*(1.d0 - rvsra))
      huini  = qvsini/qsaini
      huini  = min(huini,1.d0)
      bvar_w1(isol) = acos(1.d0-2.d0*huini)/acos(-1.d0)
    endif

    ! If the deep reservoir is uninitialized,
    ! we set the ratio w1/w2 to 1 (layer equilibrium)
    if (bvar_w2(isol).lt.1.d-20) then
      bvar_w2(isol) = bvar_w1(isol)
    endif

  endif

enddo

return
end subroutine
