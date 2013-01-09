!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine ppini1
!================


!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES OPTIONS DES VARIABLES SELON
!   LE TYPE DE PHYSIQUE PARTICULIERE
!   EN COMPLEMENT DE CE QUI A ETTE DEJA FAIT DANS USINI1

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

!===============================================================================

!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! ---> Combustion gaz : Flamme de diffusion  (Chimie 3 points)
!                       Flamme de premelange (Modele EBU)
!                       Flamme de premelange (Modele LWC)

if ( ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0                 &
                         .or. ippmod(icolwc).ge.0 ) then
  call coini1
  !==========
endif

! ---> Physique particuliere : Combustion Charbon Pulverise

if ( ippmod(iccoal).ge.0 ) then
  call cs_coal_param
  !=================
endif

! ---> Physique particuliere : Combustion Eulerienne Charbon Pulverise
!                              Couplee Transport Lagrangien

if ( ippmod(icpl3c).ge.0 ) then
  call cplin1
  !==========
endif

! ---> Physique particuliere : Combustion fuel

if ( ippmod(icfuel).ge.0 ) then
  call cs_fuel_param
  !=================
endif

! ---> Physique particuliere : Compressible

if ( ippmod(icompf).ge.0) then
  call cfini1
  !==========
endif

! ---> Physique particuliere : Versions electriques

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then
  call elini1
  !==========
endif

! ---> Physique particuliere : Ecoulements atmospheriques

if ( ippmod(iatmos).ge.0 ) then
  call atini1
  !==========
endif

! ---> Physique particuliere : Aerorefrigerants

if ( ippmod(iaeros).ge.0) then
  call ctini1
  !==========
endif

return
end subroutine
