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

subroutine memcbr &
!================

 ( iicelb , ncelet , ncel   , nfabor ,                            &
   ncelbr , ifinia ,                                              &
   ifabor ,                                                       &
   ia     )

!===============================================================================

!  FONCTION  :
!  --------
!            REPERAGE DES ELEMENTS AYANT AU MOINS UNE FACE DE BORD
!            POUR LA RECONSTRUCTION DES GRADIENTS (RECALCUL DE
!            COCG UNIQUEMENT AUX BORDS)

!          ON CALCULE NCELBR ET ON REMPLIT IA(IICELB)



!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iicelb           ! e  ! <-- ! pointeur de la premiere cas libre du           !
!                  !    !     !  tableau  ia (et pointeur sur icelbr)          !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! ncelbr           ! e  ! --> ! nombre d'elements ayant au moins une           !
! ifinia           ! e  ! --> ! iicelb+ncelbr debut de zone libre              !
!                  !    !     ! dans ia en sortie                              !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ia               ! te ! --- ! tableau de travail entier                      !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

!===============================================================================

! Arguments

integer          iicelb , ncelet , ncel   , nfabor
integer          ncelbr , ifinia
integer          ifabor(nfabor)
integer          ia(*)

! Local variables

integer          ifac, ii, iel , iiasse

!===============================================================================


!===============================================================================
! 1. MEMOIRE
!===============================================================================

! On surdimensionne ICELBR (NFABOR)
!     et un tableau d'assemblage NCELET

iiasse = iicelb + nfabor
ifinia = iiasse + ncelet

CALL IASIZE('MEMCBR',IFINIA)
!==========

!===============================================================================
! 1. INITIALISATION
!===============================================================================

do iel = 1, ncelet
  ia(iiasse+iel-1) = 0
enddo

!===============================================================================
! 2. REPERAGE DES ELEMENTS DE BORD
!===============================================================================

do ifac = 1, nfabor
  iel = ifabor(ifac)
  ia(iiasse+iel-1) = ia(iiasse+iel-1) + 1
enddo

!===============================================================================
! 3. REMPLISSAGE DE ICELBR
!===============================================================================

ii = 0
do iel = 1, ncel
  if (ia(iiasse+iel-1).gt.0) then
    ii = ii + 1
    ia(iicelb+ii-1) = iel
  endif
enddo

ncelbr = ii

!===============================================================================
! 4. DESALLOCATION DE LA MEMOIRE EN TROP ET DU TABLEAU DE TRAVAIL IASSE
!===============================================================================

ifinia = iicelb + ncelbr

!===============================================================================
! 5. FIN
!===============================================================================

return
end subroutine
