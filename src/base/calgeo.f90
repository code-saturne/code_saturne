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

subroutine calgeo &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   volmin , volmax , voltot ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
!  FONCTION  :
!  ---------

!  CALCUL DES ENTITES GEOMETRIQUES DEDUITES
!     DU JEU DE DONNEES MINIMAL

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nnod             ! e  ! <-- ! nombre de sommets                              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr (optionnel          !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! nfml  ,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds (optionnel)              !
! (ndim,nnod)      !    !     !                                                !
! volume           ! tr ! --> ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! volmin           ! r  ! --> ! volume de controle minimal                     !
! volmax           ! r  ! --> ! volume de controle maximal                     !
! voltot           ! r  ! --> ! volume total du domaine                        !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

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
include "pointe.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision volmin, volmax, voltot
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer idebia, idebra

!===============================================================================
! 1. ON SAUVEGARDE LA POSITION DE LA MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2. ON CALCULE LE VOLUME MIN et TOTAL DES ELEMENTS
!===============================================================================

call clvolc                                                       &
!==========
     ( ncelet , ncel   ,                                          &
       volmin , volmax , voltot , volume )

!===============================================================================
! 3. ON CALCULE LES SURFACES DES FACES
!===============================================================================

call clsurn                                                       &
!==========
     ( idebia , idebra ,                                          &
       nfac   , nfabor ,                                          &
       surfac , surfbo ,                                          &
       ra(isrfan) , ra(isrfbn) ,                                  &
       ia     , ra     )


!===============================================================================
! 4. ON CALCULE LE PRODUIT SCALAIRE DE LA NORMALE NORMEE A UNE FACE ET
!        DU VECTEUR DEFINI PAR LES VOISINS (VOISIN 1 : ORIGINE,
!        VOISIN 2 : EXTREMITE)
!               LA PONDERATION RESULTANTE   POND  = D2/(D1+D2)
!        OU D1 ET D2 SONT LES PROJETES SUR LA NORMALE A LA FACE DES
!        VECTEURS DEFINIS RESPECTIVEMENT PAR
!    D1 : (ORIGINE : VOISIN 1, EXTREMITE : CENTRE DE GRAVITE DE LA FACE)
!    D2 : (ORIGINE : CENTRE DE GRAVITE DE LA FACE, EXTREMITE : VOISIN 2)
!===============================================================================

call cldipo                                                       &
!==========
 ( idebia , idebra ,                                              &
   nfac   , nfabor , ncelet , ncel   ,                            &
   ifacel , ifabor ,                                              &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   ra(isrfan) , ra(isrfbn) ,                                      &
   ra(idist)  , ra(idistb) , ra(ipond) ,                          &
   ia     , ra     )

!===============================================================================
! 5. ON CALCULE LES VECTEURS IIP ET JJP POUR LES RECONSTRUCTIONS
!===============================================================================

call cldijp                                                       &
!==========
 ( idebia , idebra ,                                              &
   nfac   , nfabor , ncelet , ncel   ,                            &
   ifacel , ifabor ,                                              &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   ra(isrfan) , ra(isrfbn) ,                                      &
   ra(ipond)  ,                                                   &
   ra(idijpf) , ra(idiipb)  , ra(idofij) ,                        &
   ia     , ra     )


!===============================================================================
! 6. FILTRAGE DU VOISINAGE ETENDU POUR LE GRADIENT PAR MOINDRES CARRES
!===============================================================================

if (imrgra.eq.3) then

  call redvse (anomax)
  !==========

endif


!===============================================================================
! 8. FIN
!===============================================================================


return
end
