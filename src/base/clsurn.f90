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

subroutine clsurn &
!================

 ( idbia0 , idbra0 ,                                              &
   nfac   , nfabor ,                                              &
   surfac , surfbo ,                                              &
   surfan , surfbn ,                                              &
   ia     , ra     )

!===============================================================================

!  FONCTION  :
!  ----------

!            CALCUL DES SURFACES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0/idbra0    ! e  ! <-- ! pointeur de la premiere cas libre des          !
!                  !    !     !  tableaux ia/ra                                !
! nfac  /nfabor    ! e  ! <-- ! nombre total de faces internes/de brd          !
! surfac           ! tr ! <-- ! coords du vecteur surface des nfaglo           !
! (3,nfac  )       !    !     !  faces internes ; dirige du vois 1             !
!                  !    !     !  vers le voisin 2 (ifacel)                     !
!                  !    !     !  non unitaire                                  !
! surfbo           ! tr ! <-- ! coords du vecteur surface des nfagbr           !
! (3,nfabor)       !    !     !  faces de bord ; dirige vers l'                !
!                  !    !     !  exterieur du domaine ; non unitaire           !
! surfan           ! tr ! <-- ! surface des faces internes                     !
! (nfac    )       !    !     !     (norme de surfac)                          !
! surfbn           ! tr ! <-- ! surface des faces de bord                      !
! (nfabor  )       !    !     !     (norme de surfbo)                          !
! ia               ! te ! --- ! tableau de travail entier                      !
! ra               ! tr ! --- ! tableau de travail reel                        !
!__________________.______________!_____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

!===============================================================================

! Arguments

integer idbia0,idbra0
integer nfac,nfabor
integer ia(*)
double precision surfac(3,nfac  ),surfbo(3,nfabor)
double precision surfan(nfac    ),surfbn(nfabor  )
double precision ra(*)

integer ifac
double precision surfx,surfy,surfz

!-------------------------------------------------------------------------------
! 1. FACES INTERNES
!-------------------------------------------------------------------------------

do ifac = 1, nfac

   surfx = surfac(1,ifac)
   surfy = surfac(2,ifac)
   surfz = surfac(3,ifac)
   surfan(ifac) = sqrt(surfx**2 + surfy**2 + surfz**2)

enddo

!-------------------------------------------------------------------------------
! 2. FACES DE BORD
!-------------------------------------------------------------------------------

do ifac = 1, nfabor

   surfx = surfbo(1,ifac)
   surfy = surfbo(2,ifac)
   surfz = surfbo(3,ifac)
   surfbn(ifac) = sqrt(surfx**2 + surfy**2 + surfz**2)

enddo

!-------------------------------------------------------------------------------
! 3. FIN
!-------------------------------------------------------------------------------

return
end
