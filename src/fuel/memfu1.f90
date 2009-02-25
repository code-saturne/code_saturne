!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

                  subroutine memfu1                               &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , ncelet , ncel   , nfac   , nfabor ,                   &
   ntbfui , ifuwi  ,                                              &
   ntbfur , ifuwr  ,                                              &
   ntbwoi , iwori  ,                                              &
   ntbwor , iworr  ,                                              &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  PHYSIQUE PARTICULIERE : FLAMME FUEL
!    GESTION MEMOIRE CALCUL DES PROPRIETES PHYSIQUES PHASE GAZ

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0/idbra0    ! e  ! <-- ! pointeur de la premiere cas libre des          !
!                  !    !     !  tableaux ia/ra                                !
! nvar             ! e  ! <-- ! nombre de variables                            !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! ntbcpi           ! e  ! <-- ! taille du macro tableau cp entiers             !
! icpwi            ! e  ! --> ! pointeur  macro tableau cp entiers             !
! ntbcpr           ! e  ! <-- ! taille du macro tableau cp reels               !
! icpwr            ! e  ! --> ! pointeur  macro tableau cp reels               !
! ntbmci           ! e  ! <-- ! taille du macro tableau mc entiers             !
! imcwi            ! e  ! --> ! pointeur  macro tableau mc entiers             !
! ntbmcr           ! e  ! <-- ! taille du macro tableau mc reels               !
! imcwr            ! e  ! --> ! pointeur  macro tableau mc reels               !
! ntbwoi           ! e  ! <-- ! taille du macro tableau work entiers           !
! iwori            ! e  ! --> ! pointeur  macro tableau work entiers           !
! ntbwor           ! e  ! <-- ! taille du macro tableau work reels             !
! iworr            ! e  ! --> ! pointeur  macro tableau work reels             !
! ifinia           ! e  ! --> ! pointeur de la premiere cas libre              !
!                  !    !     !  dans ia en sortie                             !
! ifinra           ! e  ! --> ! pointeur de la premiere cas libre              !
!                  !    !     !  dans ra en sortie                             !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "optcal.h"

!===============================================================================

! Arguments

integer          idbia0 ,idbra0
integer          nvar
integer          ncelet , ncel   , nfac   , nfabor
integer          ntbfui , ifuwi  , ntbfur , ifuwr
integer          ntbwoi , iwori  , ntbwor , iworr
integer          ifinia , ifinra

! VARIABLES LOCALES

integer          idebia , idebra

!===============================================================================

!---> INITIALISATION

idebia = idbia0
idebra = idbra0

!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA

ifuwi   =       idebia
iwori   =       ifuwi   + ntbfui*ncelet
ifinia =        iwori   + ntbwoi*ncelet

ifuwr   =       idebra
iworr   =       ifuwr   + ntbfur*ncelet
ifinra =        iworr   + ntbwor*ncelet

!---> VERIFICATION

CALL IASIZE('MEMFU1',IFINIA)
!==========

CALL RASIZE('MEMFU1',IFINRA)
!==========

return
end
