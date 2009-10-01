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

subroutine cldipo &
!================

 ( idbia0 , idbra0 ,                                              &
   nfac   , nfabor , ncelet , ncel   ,                            &
   ifacel , ifabor ,                                              &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , surfan , surfbn , &
   dist   , distbo , pond   ,                                     &
   ia     , ra     )

!===============================================================================

!  FONCTION  :
!  --------
!         CALCUL DES DISTANCES RELATIVES AUX FACES
!                     ET DES PONDERATIONS ASSOCIEES


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0/idbra0    ! e  ! <-- ! pointeur de la premiere cas libre des          !
!                  !    !     !  tableaux ia/ra                                !
! nfac  /nfabor    ! e  ! <-- ! nombre total de faces internes/de bor          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! ifacel           ! te ! <-- ! ifacel(iel,ifac) num glob de l'elemnt          !
! (2,nfac)         !    !     !  vois iel (1 ou 2) de la fac int ifac          !
! ifabor           ! te ! <-- ! ifabor(ifac) num glob de l'element             !
! nfac  ,2)        !    !     !  voisin iel (1) de la fac de brd ifac          !
! xyzcen           ! tr ! <-- ! coords du "centre" des nelem elements          !
! (3,ncelet )      !    !     !                                                !
! surfac           ! tr ! <-- ! coords du vecteur surface des nfaglo           !
! (3,nfac  )       !    !     ! faces internes; dirige du vois 1 vers          !
!                  !    !     !  le voisin 2 (ifacel) ; non unitaire           !
! surfbo           ! tr ! <-- ! coords du vecteur surface des nfagbr           !
! (3,nfabor)       !    !     !  faces de brd;dirige vers l'exterieur          !
!                  !    !     !  du domaine ; non unitaire                     !
! cdgfac           ! tr ! <-- ! coords du centre de gravite des faces          !
! (ndim,nfac  )    !    !     !                              internes          !
! cdgfbo           ! tr ! <-- ! coords du centre de gravite des faces          !
! (ndim,nfabor)    !    !     !                              de bord           !
! dist             ! tr ! --> ! prod.scal. de la normale normee a une          !
! (nfac  )         !    !     !  face interne et du vecteur defini             !
!                  !    !     !  par les voisins (voisin 1 : origine,          !
!                  !    !     !  voisin 2 : extremite)                         !
! pond             ! tr ! --> ! ponderation pour face interne                  !
! (nfac  )         !    !     !  = d2/(d1+d2) ou d1 et d2 sont les             !
!                  !    !     !  projetes sur la normale a la face             !
!                  !    !     !  des vecteurs definis resp. par :              !
!                  !    !     !d1(orig: voisin 1, extremite: cdg face          !
!                  !    !     !d2(orig: cdg face, extremite: voisin 2          !
! distbo           ! tr ! --> ! prod. scal. de la normale normee               !
! (nfabor)         !    !     !  une face de bord et du vecteur                !
!                  !    !     !  voisin cdg face ((voisin : origine,           !
!                  !    !     !  cdg face : extremite)                         !
! ia               ! te ! --- ! tableau de travail entier                      !
! ra               ! tr ! --- ! tableau de travail reel                        !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

include "paramx.h"
include "cstnum.h"
include "entsor.h"

integer           idbia0, idbra0
integer           nfac  , nfabor, ncelet, ncel
integer           ifacel(2,nfac), ifabor(nfabor)
integer           ia(*)
double precision  xyzcen(3,ncelet)
double precision  surfac(3,nfac  ),surfbo(3,nfabor)
double precision  cdgfac(3,nfac  ),cdgfbo(3,nfabor)
double precision  surfan(nfac    ),surfbn(nfabor)
double precision  dist(nfac)
double precision  distbo(nfabor)
double precision  pond(nfac)
double precision  ra(*)

integer           ifac, ivois1,ivois2,ipond0
double precision  surfx,surfy,surfz, surfn,dist2f
double precision  xvn,yvn,zvn,xvv,yvv,zvv

!===============================================================================
! 1.  TRAITEMENT DES FACES INTERNES
!===============================================================================

ipond0 = 0
do ifac = 1, nfac

   surfx = surfac(1,ifac)
   surfy = surfac(2,ifac)
   surfz = surfac(3,ifac)
   surfn = surfan(ifac)

   ivois1 = ifacel(1,ifac)
   ivois2 = ifacel(2,ifac)

! (CDG FACE,VOISIN2)
   xvn = xyzcen(1,ivois2) - cdgfac(1,ifac)
   yvn = xyzcen(2,ivois2) - cdgfac(2,ifac)
   zvn = xyzcen(3,ivois2) - cdgfac(3,ifac)

! (PRODUIT SCALAIRE AVEC LA NORMALE NORMEE)
   dist2f       =    (xvn*surfx + yvn*surfy + zvn*surfz)/surfn

! (VOISIN1, VOISIN2)
   xvv = xyzcen(1,ivois2) - xyzcen(1,ivois1)
   yvv = xyzcen(2,ivois2) - xyzcen(2,ivois1)
   zvv = xyzcen(3,ivois2) - xyzcen(3,ivois1)

! (PRODUIT SCALAIRE AVEC LA NORMALE NORMEE)
   dist(ifac) =    (xvv*surfx + yvv*surfy + zvv*surfz)/surfn

   if (abs(dist(ifac)).ge.epzero) then
      pond(ifac) = dist2f / dist(ifac)
   else
      ipond0 = ipond0 + 1
      pond(ifac) = 0.5d0
   endif

enddo

!===============================================================================
! 2.  TRAITEMENT DES FACES DE BORD
!===============================================================================

do ifac = 1, nfabor

   surfx = surfbo(1,ifac)
   surfy = surfbo(2,ifac)
   surfz = surfbo(3,ifac)
   surfn = surfbn(ifac)

   ivois1 = ifabor(ifac)

! (VOISIN 1, CDG FACE)
   xvn = cdgfbo(1,ifac) - xyzcen(1,ivois1)
   yvn = cdgfbo(2,ifac) - xyzcen(2,ivois1)
   zvn = cdgfbo(3,ifac) - xyzcen(3,ivois1)

   distbo(ifac) =    (xvn*surfx + yvn*surfy + zvn*surfz)/surfn

enddo

!===============================================================================
! 3.  FIN
!===============================================================================

if(ipond0.ne.0) then
  write(nfecra,1000) ipond0
endif

#if defined(_CS_LANG_FR)

 1000 format(' CLDIPO : ',I10,' DISTANCES NULLES ENTRE CENTRES ',/,     &
       '          ON MET POND A 0.5 ')

#else

 1000 format(' CLDIPO : ',I10,' NULL  DISTANCES BETWEEN CENTRES ',/,    &
       '          POND IS SET TO 0.5 ')

#endif

return
end subroutine
