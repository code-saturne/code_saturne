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

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[out]    nmodpp        number of activated paricle physic models
!______________________________________________________________________________

subroutine ppvarp &
( nmodpp )
!================


!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES POSITIONS DES VARIABLES SELON
!   LE TYPE DE PHYSIQUE PARTICULIERE
! REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP

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
use elincl
use atincl
use ihmpre

!===============================================================================

implicit none

! Arguments

integer       nmodpp

!===============================================================================

! ---> Physique particuliere : Combustion Gaz

if (itherm.ne.0 .and. nmodpp.eq.0) then
  if (iihmpr.eq.1) then
    call uithsc(iscalt)
  endif

endif

if (  ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0                &
                          .or. ippmod(icolwc).ge.0 ) then
  call covarp
  !==========
endif

if ( ippmod(iccoal).ge.0 ) then
  call cs_coal_varpos
  !==================
endif

! ---> Physique particuliere :  Combustion Charbon Pulverise
!      Couplee Transport Lagrangien des particules de charbon

if ( ippmod(icpl3c).ge.0 ) then
  call cplvar
  !==========
endif

! ---> Physique particuliere :  Combustion Fuel

if ( ippmod(icfuel).ge.0 ) then
  call cs_fuel_varpos
  !==================
endif

! ---> Physique particuliere : Compressible

if ( ippmod(icompf).ge.0) then
  call cfvarp
  !==========
endif

! ---> Physique particuliere : Versions Electriques

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then
  call elvarp
  !==========
endif

! ---> Physique particuliere : Version Atmospherique

if ( ippmod(iatmos).ge.1 ) then
  call atvarp
  !==========
endif

! ---> Physique particuliere : Aerorefrigerants

if ( ippmod(iaeros).ge.0 ) then
  call ctvarp
  !==========
endif

return
end subroutine
