!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine pplecd
!================

!===============================================================================
!  FONCTION  :
!  ---------

! LECTURE DU FICHIER DE DONNEES PHYSIQUE PARTICULIERE

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
use pointe
use entsor
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

! Arguments


! Local variables


!===============================================================================

!===============================================================================
! 1. AIGUILLAGE VERS LE MODELE ADEQUAT
!===============================================================================

! ---> Flamme de diffusion  - Chimie 3 points
!      Flamme de premelange - Modele EBU
!      Flamme de premelange - Modele LWC

if ( ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0                 &
                         .or. ippmod(icolwc).ge.0  ) then
  call colecd
  !==========
endif


! ---> Flamme charbon pulverise

if ( ippmod(icp3pl).ge.0 ) then
  call cplecd
  !==========
endif

if ( ippmod(iccoal).ge.0 ) then
  call cs_coal_readata
  !==================
endif
! ---> Combustion charbon pulverise couple transport Lagrangien
!      des particules de charbon

if ( ippmod(icpl3c).ge.0 ) then
  call cplecd
  !==========
endif

! ---> Flamme fuel

if ( ippmod(icfuel).ge.0 ) then
  call cs_fuel_readata
  !==================
endif

! ---> Version Electrique : Effet Joule, Arc Electrique,
!                           Conduction Ionique

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1      ) then
  call ellecd
  !==========
endif

!----
! FIN
!----

return
end subroutine

