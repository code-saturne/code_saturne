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
!> \param[in]   nvar        total number of variables
!> \param[in]   nscal       total number of scalars
!> \param[in]   dt          time step value
!-------------------------------------------------------------------------------

subroutine ppiniv0 &
 ( nvar   , nscal  ,                                              &
   dt     )

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

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!===============================================================================

! ---> Combustion gaz
!      Flamme de diffusion : chimie 3 points
if ( ippmod(icod3p).ge.0 ) then
  call d3pini                                                     &
    ( nvar   , nscal  ,                                              &
    dt     )
endif

! ---> Combustion gaz
!      Flamme de premelange : modele EBU

if (ippmod(icoebu).ge.0) then
  call ebuini(nvar, nscal, dt)
endif

! ---> Combustion gaz
!      Flamme de premelange : modele LWC

if (ippmod(icolwc).ge.0) then
  call lwcini(nvar, nscal, dt)
endif

! ---> Combustion charbon pulverise

if (ippmod(iccoal).ge.0) then
  call cs_coal_varini(nvar, nscal, dt)
endif

! ---> Combustion charbon pulverise couples Lagrangien

if (ippmod(icpl3c).ge.0) then
  call cplini
endif

! ---> Combustion fuel

if  (ippmod(icfuel).ge.0) then
  call cs_fuel_varini(nvar, nscal, dt)
endif

! ---> Version electrique
!      Effet Joule
!      Conduction ionique

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1       ) then

  call eliniv(isuite)
endif

! Atmospheric flows, first stage
if (ippmod(iatmos).ge.0) then

  call atiniv0                                                     &
 ( nvar   , nscal  ,                                              &
   dt     )

endif

! ---> Cooling towers

if (ippmod(iaeros).ge.0) then
  call ctiniv0(nvar, nscal, dt)
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
!> \param[in]   nvar        total number of variables
!> \param[in]   nscal       total number of scalars
!> \param[in]   dt          time step value
!-------------------------------------------------------------------------------


subroutine ppiniv1 &
 ( nvar   , nscal  ,                                              &
   dt     )

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

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!===============================================================================

! Atmospheric flows, second stage
if (ippmod(iatmos).ge.0) then

  call atiniv1                                                     &
 ( nvar   , nscal  ,                                              &
   dt     )

endif

! ---> Cooling towers
if (ippmod(iaeros).ge.0) then
  call ctiniv1(nvar, nscal, dt)
endif

! Gas mixture modelling in presence of noncondensable gases and
! condensable gas as stream.
if (ippmod(igmix).ge.0) then
  call cs_gas_mix_initialization &
  ( nvar   , nscal  ,                                            &
    dt     )
endif

! Compressible
! Has to be called AFTER the gas mix initialization because the
! mixture composition is taken into account in the thermodynamic
! law, if gas mix specific physics is enabled.
if (ippmod(icompf).ge.0) then
  call cfiniv(nvar, nscal, dt)
endif

!----
! Formats
!----

!----
! End
!----

return
end subroutine
