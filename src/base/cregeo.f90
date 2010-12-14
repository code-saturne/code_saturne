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
   nideve , nrdeve , nituse , nrtuse ,                            &
   idevel , ituser , ia     ,                                     &
   rdevel , rtuser , ra     )

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
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
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
integer          nideve , nrdeve , nituse , nrtuse

integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia , idebra
integer          ilcel  , ilfaci , ilfacb
integer          maxelt , ils
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
 ( ifinia , ifinra ,                                              &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ia(ilcel) , ia(ilfaci) , ia(ilfacb)  ,                         &
   idevel , ituser , ia     ,                                     &
   rdevel , rtuser , ra     )


!===============================================================================
! 3. CREATION DU MAILLAGE EXTRAIT COUPLE AVEC SYRTHES
!    ENVOI DES DONNEES GEOMETRIQUES A SYRTHES SI NECESSAIRE
!===============================================================================

!     NOMBRE DE COUPLAGES SYRTHES DEFINIS

call nbcsyr(nbrsyr)
!==========

if (nbrsyr .gt. 0) then
  call geosyr(ichrsy)
  !==========
endif


!===============================================================================
! 3. CREATION DU MAILLAGE EXTRUDE POUR LES ZONES D'ECHANGE AERO
!===============================================================================

if (ippmod(iaeros).ge.0) then

  maxelt = max(ncelet,nfac,nfabor)

  ils    = idebia
  ifinia = ils + maxelt
  call iasize('cregeo',ifinia)
  !==========

  call usctdz                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   nideve , nrdeve , nituse , nrtuse ,                            &
   maxelt , ia(ils),                                              &
   idevel , ituser , ia     ,                                     &
   rdevel , rtuser , ra     )

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
! 5. CALCUL DES TABLEAUX COMPLEMENTAIRES
!===============================================================================

call calgeo                                                       &
!==========
 ( idebia , idebra ,                                              &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idevel , ituser , ia     ,                                     &
   volmin , volmax , voltot ,                                     &
   rdevel , rtuser , ra     )

return
end subroutine
