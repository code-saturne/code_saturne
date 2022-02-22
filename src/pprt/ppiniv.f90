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

subroutine ppiniv0

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

!===============================================================================

implicit none

! Local variables

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

! ---> Combustion charbon pulverise couples Lagrangien

if (ippmod(icpl3c).ge.0) then
  call cplini
endif

! ---> Combustion fuel

if (ippmod(icfuel).ge.0) then
  call cs_fuel_varini
endif

! Atmospheric flows, first stage
if (ippmod(iatmos).ge.0) then
  call atiniv0
endif

! ---> Cooling towers

if (ippmod(iaeros).ge.0) then
  call ctiniv0
endif

! Electric arcs, Joule effect or ionic conduction

if (ippmod(ieljou).ge.1 .or. ippmod(ielarc).ge.1) then
  call eliniv(isuite)
endif

!----
! Formats
!----

!----
! End
!----

return
end subroutine

!-------------------------------------------------------------------------------
!> \file ppiniv.f90
!> \brief Initialisation of specific physic variables
!>
!>  Physical properties SHOULD not be modified here
!>  Second stage after the GUI and user modifications.
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!-------------------------------------------------------------------------------


subroutine ppiniv1

!===============================================================================

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

!===============================================================================

implicit none

! Local variables

!===============================================================================

! Atmospheric flows, second stage
if (ippmod(iatmos).ge.0) then
  call atiniv1
endif

! ---> Cooling towers
if (ippmod(iaeros).ge.0) then
  call ctiniv1
endif

! Gas mixture modelling in presence of noncondensable gases and
! condensable gas as stream.
if (ippmod(igmix).ge.0) then
  call cs_gas_mix_initialization
endif

! Compressible
! Has to be called AFTER the gas mix initialization because the
! mixture composition is taken into account in the thermodynamic
! law, if gas mix specific physics is enabled.
if (ippmod(icompf).ge.0) then
  call cfiniv
endif

!----
! Formats
!----

!----
! End
!----

return
end subroutine
