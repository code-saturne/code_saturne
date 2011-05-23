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

subroutine cregeo &
!================

 ( idbia0 , idbra0 ,                                              &
   ia     ,                                                       &
   ra     )

!===============================================================================
!  FONCTION
!  ---------

!   CREATION DES ENTITES GEOMETRIQUES, PAR  CALCUL

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ia(*)            ! ia ! --- ! main integer work array                        !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use entsor
use optcal
use cstphy
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0

integer          ia(*)

double precision ra(*)

! Local variables

integer          idebia , idebra
integer          ilcel  , ilfaci , ilfacb
integer          ifinia , ifinra , nbrsyr , nbzech

character        ficsui*32

!===============================================================================
! 1. INITIALISATION DE LA MEMOIRE EN LOCAL POUR LES TAB DE TRAV
!===============================================================================

idebia = idbia0
idebra = idbra0

!     MEMOIRE DE TRAVAIL POUR LA DEFINITION DE COUPES

ilcel  = idebia
ilfaci = ilcel  + ncelet
ilfacb = ilfaci + nfac
ifinia = ilfacb + nfabor

!     VERIFICATION DE LA DISPONIBILITE DE LA MEMOIRE

call iasize ('cregeo', ifinia)
!==========

!===============================================================================
! 2. DEFINITION DE MAILLAGES ET FORMATS DE POST TRAITEMENT UTILISATEUR
!===============================================================================

call usdpst                                                       &
!==========
 ( ia(ilcel) , ia(ilfaci) , ia(ilfacb)  ,                         &
   ia     ,                                                       &
   ra     )


!===============================================================================
! 3. CREATION DU MAILLAGE EXTRAIT COUPLE AVEC SYRTHES
!    ENVOI DES DONNEES GEOMETRIQUES A SYRTHES SI NECESSAIRE
!===============================================================================

!     NOMBRE DE COUPLAGES SYRTHES DEFINIS

call nbcsyr(nbrsyr)
!==========

if (nbrsyr .gt. 0) then
  call geosyr
  !==========
endif


!===============================================================================
! 3. CREATION DU MAILLAGE EXTRUDE POUR LES ZONES D'ECHANGE AERO
!===============================================================================

if (ippmod(iaeros).ge.0) then


  call usctdz                                                     &
  !==========
 ( ia     ,                                                       &
   ra     )

  call nbzect(nbzech)
  !==========

  if (nbzech .gt. 0) then
    call geoct
    !=========
    if (ichrze.gt.0) then
      call pstict
      !==========
    endif
  endif

  if (ippmod(iaeros).ge.0.and.isuict.eq.1) then
     ficsui = 'cooling_towers'
     call lecctw ( ficsui , len(ficsui) )
     !==========
  endif

endif


!===============================================================================
! 4. ECRITURE DES MAILLAGES DE POST TRAITEMENT INDEPENDANTS DU TEMPS
!===============================================================================

call pstema (ntcabs, ttcabs)
!==========


!===============================================================================
! 5. FILTRAGE DU VOISINAGE ETENDU POUR LE GRADIENT PAR MOINDRES CARRES
!===============================================================================

if (imrgra.eq.3) then
  call redvse (anomax)
  !==========
endif


return
end subroutine
