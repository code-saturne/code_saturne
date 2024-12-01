!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!-------------------------------------------------------------------------------
!> \file ppiniv.f90
!> \brief Initialisation of specific physic variables
!>
!>  Physical properties SHOULD not be modified here
!>  First stage BEFOR user input (GUI and user subroutines).
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!-------------------------------------------------------------------------------

subroutine ppiniv0() &
  bind(C, name='cs_f_ppiniv0')

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
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_atmo_fields_init0()  &
    bind(C, name='cs_atmo_fields_init0')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_atmo_fields_init0

end interface

!===============================================================================

! ---> Combustion gaz
!      Flamme de diffusion : chimie 3 points
if (ippmod(icod3p).ge.0) then
  call d3pini
endif

if (ippmod(islfm).ge.0) then
  call cs_steady_laminar_flamelet_init
endif
! ---> Combustion gaz
!      Flamme de premelange : modele EBU

if (ippmod(icoebu).ge.0) then
  call ebuini
endif

! ---> Combustion gaz
!      Flamme de premelange : modele LWC

if (ippmod(icolwc).ge.0) then
  call lwcini
endif

! ---> Combustion charbon pulverise

if (ippmod(iccoal).ge.0) then
  call cs_coal_varini
endif

! Atmospheric flows, first stage
if (ippmod(iatmos).ge.0) then
  call atiniv0
  call cs_atmo_fields_init0
endif

!----
! Formats
!----

!----
! End
!----

return
end subroutine
