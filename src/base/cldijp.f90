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

subroutine cldijp &
!================

 ( idbia0 , idbra0 ,                                              &
   nfac   , nfabor , ncelet , ncel   ,                            &
   ifacel , ifabor ,                                              &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   surfan , surfbn ,                                              &
   pond   ,                                                       &
   dijpf  , diipb  , dofij  ,                                     &
   ia     , ra     )

!===============================================================================

!  FONCTION  :
!  --------

!    CALCUL DE VECTEURS POUR TERMES DE NON ORTHOGONALITE :

!       SOIT UNE FACE A,B,C ET I,J LES CENTRES DES VOLUMES VOISINS.
!         (SEUL I EST DEFINI POUR UNE FACE DE BORD)

!       LA FACE EST ORIENTEE DE I VERS J, DE NORMALE NIJ.
!         (LES FACES DE BORD SONT ORIENTEES VERS L'EXTERIEUR)
!       LA NORME DE NIJ VAUT 1.
!       SIJ EST LA SURFACE DE LA FACE ABC.
!       ON DEFINIT I' ET J' LES PROJETES ORTHOGONAUX DE
!         I ET J SUR LA DROITE ORTHOGONALE A ABC PASSANT PAR SON CDG F
!         (SEUL I' EST DEFINI POUR LES FACES DE BORD)

!       ON CALCULE ICI LE VECTEUR I'J' AUX FACES INTERNES (DIJPF)
!                      LE VECTEUR II'  AUX FACES DE BORD  (DIIPB)
!                      LE VECTEUR OF   AUX FACES INTERNES (DOFIJ)

!       NOTER LES RELATIONS
!                 II' = IG - (IG.NIJ)NIJ
!                 JJ' = JG - (JG.NIJ)NIJ

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0/idbra0    ! e  ! <-- ! pointeur de la premiere cas libre des          !
!                  !    !     !  tableaux ia/ra                                !
! nfac             ! e  ! <-- ! nombre total de faces internes                 !
! nfabor           ! e  ! <-- ! nombre total de faces         de bord          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! ifacel           ! te ! <-- ! ifacel(iel,ifac) num glob de l'elemnt          !
! (2,nfac)         !    !     !  vois iel (1 ou 2) de la fac int ifac          !
! ifabor           ! te ! <-- ! ifabor(ifac    ) num glob de l'elt             !
! nfabor  )        !    !     !  voisin iel (1) de la fac de brd ifac          !
! xyzcen           ! tr ! <-- ! coords du "centre" des nelem elements          !
! (3,ncelet        !    !     !                                                !
! surfac           ! tr ! <-- ! coords du vecteur surface des nfac             !
! (3,nfac  )       !    !     ! faces internes; dirige du vois 1 vers          !
!                  !    !     !  le vois 2 (ifacel) ; non unitaire             !
! surfbo           ! tr ! <-- ! coords du vecteur surface des nfabor           !
! (3,nfabor)       !    !     !  faces de bord ; dirige vers                   !
!                  !    !     !  l'exterieur du domaine; non unitaire          !
! surfan           ! tr ! <-- ! norme de surfac (surface des faces             !
! (nfac    )       !    !     ! internes)                                      !
! surfbn           ! tr ! <-- ! norme de surfbo (surface des faces             !
! (nfabor  )       !    !     !  de bord)                                      !
! cdgfac           ! tr ! <-- ! coords du centre de gravite des faces          !
! 3,nfac  )        !    !     !           internes                             !
! cdgfbo           ! tr ! <-- ! coords du centre de gravite des faces          !
! 3,nfabor)        !    !     !           de  bord                             !
! pond             ! tr ! --> ! ponderation pour face interne                  !
! (nfac  )         !    !     !  = d2/(d1+d2) ou d1 et d2 sont les             !
!                  !    !     !  projetes sur la normale a la face             !
!                  !    !     !  des vecteurs definis resp. par :              !
!                  !    !     !d1(orig: voisin 1, extremite: cdg face          !
!                  !    !     !d2(orig: cdg face, extremite: voisin 2          !
! dijpf            ! tr ! --> ! vecteur i'j' pour les faces internes           !
! (ndim,nfac  )    !    !     !                                                !
! diipb            ! tr ! --> ! vecteur ii' pour les faces de bord             !
! (ndim,nfabor)    !    !     !                                                !
! dofij            ! tr ! --> ! vecteur of pour les faces internes             !
! (ndim,nfac  )    !    !     ! o : intersection de ij et la face              !
!                  !    !     ! f : centre de la face                          !
! ia               ! te ! --- ! tableau de travail entier                      !
! ra               ! tr ! --- ! tableau de travail reel                        !
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

!===============================================================================

! Arguments

integer idbia0,idbra0
integer nfac,nfabor,ncelet,ncel
integer ifacel(2,nfac)
integer ifabor(nfabor)
integer ia(*)
double precision xyzcen(3,ncelet)
double precision surfac(3,nfac),surfbo(3,nfabor)
double precision surfan(nfac),surfbn(nfabor)
double precision cdgfac(3,nfac),cdgfbo(3,nfabor)
double precision pond(nfac)
double precision dijpf(3,nfac)
double precision diipb(3,nfabor)
double precision dofij(3,nfac)
double precision ra(*)

integer ifac,ivois1,ivois2
double precision surfnx,surfny,surfnz
double precision vecigx,vecigy,vecigz
double precision vecijx,vecijy,vecijz
double precision dipjp
double precision psi

!===============================================================================
! 1. FACES INTERNES
!===============================================================================

do ifac = 1, nfac

!--->  NUMERO DES VOISINS

   ivois1 = ifacel(1,ifac)
   ivois2 = ifacel(2,ifac)

!---> NORMALE NORMEE

   surfnx = surfac(1,ifac)/surfan(ifac)
   surfny = surfac(2,ifac)/surfan(ifac)
   surfnz = surfac(3,ifac)/surfan(ifac)

!---> IJ

   vecijx = xyzcen(1,ivois2)-xyzcen(1,ivois1)
   vecijy = xyzcen(2,ivois2)-xyzcen(2,ivois1)
   vecijz = xyzcen(3,ivois2)-xyzcen(3,ivois1)

!---> DIJPP = IJ.NIJ

   dipjp  = vecijx*surfnx + vecijy*surfny +                       &
            vecijz*surfnz

!---> DIJPF = (IJ.NIJ).NIJ

   dijpf(1,ifac) = dipjp * surfnx
   dijpf(2,ifac) = dipjp * surfny
   dijpf(3,ifac) = dipjp * surfnz

!---> DOFIJ = OF

  dofij(1,ifac) = cdgfac(1,ifac) -                                &
                  (pond(ifac) * xyzcen(1,ivois1) +                &
                  (1.d0-pond(ifac))*xyzcen(1,ivois2))
  dofij(2,ifac) = cdgfac(2,ifac) -                                &
                  (pond(ifac)*xyzcen(2,ivois1) +                  &
                  (1.d0-pond(ifac))*xyzcen(2,ivois2))
  dofij(3,ifac) = cdgfac(3,ifac) -                                &
                  (pond(ifac)*xyzcen(3,ivois1) +                  &
                  (1.d0-pond(ifac))*xyzcen(3,ivois2))

enddo

!===============================================================================
! 2. FACES DE BORD
!===============================================================================

do ifac = 1, nfabor

!--->  NUMERO DU VOISIN

   ivois1 = ifabor(ifac)

!---> NORMALE NORMEE

   surfnx = surfbo(1,ifac)/surfbn(ifac)
   surfny = surfbo(2,ifac)/surfbn(ifac)
   surfnz = surfbo(3,ifac)/surfbn(ifac)

!---> IG

   vecigx = cdgfbo(1,ifac)-xyzcen(1,ivois1)
   vecigy = cdgfbo(2,ifac)-xyzcen(2,ivois1)
   vecigz = cdgfbo(3,ifac)-xyzcen(3,ivois1)

!---> PSI = IG.NIJ

   psi = vecigx*surfnx+vecigy*surfny+vecigz*surfnz

!---> DIIPB = IG - (IG.NIJ)NIJ

   diipb(1,ifac) = vecigx - psi*surfnx
   diipb(2,ifac) = vecigy - psi*surfny
   diipb(3,ifac) = vecigz - psi*surfnz

enddo

!===============================================================================
! 3. FIN
!===============================================================================

return
end subroutine
