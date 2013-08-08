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

subroutine ppprop &
!================

 ( ipropp , ipppst )

!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES POSITIONS DES VARIABLES D'ETAT SELON
!   LE TYPE DE PHYSIQUE PARTICULIERE
!   (DANS VECTEURS PROPCE, PROPFB)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ipropp           ! e  ! <-- ! numero de la derniere propriete                !
!                  !    !     !  (les proprietes sont dans propce ou propfb)   !
! ipppst           ! e  ! <-- ! pointeur indiquant le rang de la               !
!                  !    !     !  derniere grandeur definie aux                 !
!                  !    !     !  cellules (rtp,propce...) pour le              !
!                  !    !     !  post traitement                               !
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
use atincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer       ipropp, ipppst

!===============================================================================

! ---> Physique particuliere : Combustion Gaz

if ( ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0                 &
                         .or. ippmod(icolwc).ge.0 ) then
  call coprop(ipropp,ipppst)
  !==========
endif

! ---> Physique particuliere :  Combustion Charbon Pulverise

if ( ippmod(iccoal).ge.0 ) then
  call cs_coal_prop(ipropp,ipppst)
  !================
endif

! ---> Physique particuliere :  Combustion Charbon Pulverise
!      Couplee Transport Lagrangien des particules de charbon

if ( ippmod(icpl3c).ge.0 ) then
  call cplpro (ipropp,ipppst)
  !==========
endif

! ---> Physique particuliere : Combustion Fuel

if ( ippmod(icfuel).ge.0 ) then
  call cs_fuel_prop(ipropp,ipppst)
  !================
endif

! ---> Physique particuliere : Compressible

if ( ippmod(icompf).ge.0 ) then
  call cfprop(ipropp,ipppst)
  !==========
endif

! ---> Physique particuliere : Versions electriques

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then
  call elprop(ipropp,ipppst)
  !==========
endif

! ---> Physique particuliere : Atmospherique

if ( ippmod(iatmos).ge.1 ) then
  call atprop(ipropp,ipppst)
  !==========
endif

end subroutine
