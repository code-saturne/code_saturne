!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

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

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "dimens.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "cstnum.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"

!===============================================================================


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

if ( ippmod(icp3pl).ge.0 ) then
  call cpini1
  !==========
endif

! ---> Physique particuliere : Combustion Eulerienne Charbon Pulverise
!                              Couplee Transport Lagrangien

if ( ippmod(icpl3c).ge.0 ) then
  call cplin1
  !==========
endif

! ---> Physique particuliere : Combustion fuel

if ( ippmod(icfuel).ge.0 ) then
  call fuini1
  !==========
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
