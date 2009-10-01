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

subroutine visort &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse , imvisf ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w2     , w3     ,                                     &
   viscf  , viscb  ,                                              &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DE LA VITESSE DE DIFFUSION "ORTHOTROPE"
! VISCF,B = VISCOSITE*SURFACE/DISTANCE, HOMOGENE A UN DEBIT EN KG/S

!         =
! (NX**2*VISC11_MOY_FACE
! +NY**2*VISC22_MOY_FACE+NZ**2*VISC33_MOY_FACE)*SURFACE/DISTANCE

! LA VISCOSITE EST DONNE PAR W1, W2, W3

! RQE : A PRIORI, PAS BESOIN DE TECHNIQUE DE RECONSTRUCTION
!  ( A AMELIORER SI NECESSAIRE )

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
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! imvisf           ! e  ! <-- ! methode de calcul de la visc face              !
!                  !    !     !  = 0 arithmetique                              !
!                  !    !     !  = 1 harmonique                                !
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
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! w1,2,3(ncelet    ! tr ! <-- ! valeurs de la viscosite                        !
! viscf(nfac)      ! tr ! --> ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --> ! visc*surface/dist aux faces de bord            !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "pointe.h"
include "period.h"
include "parall.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nideve , nrdeve , nituse , nrtuse , imvisf

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          ifac, ii, jj
integer          idimte, itenso
double precision viscxi, viscxj, viscyi, viscyj, visczi, visczj
double precision sx2, sy2, sz2, dist, pond, surfn

!===============================================================================

! ---> TRAITEMENT DU PARALLELISME

if(irangp.ge.0) then
  call parcom (w1)
  !==========
  call parcom (w2)
  !==========
  call parcom (w3)
  !==========
endif

! ---> TRAITEMENT DE LA PERIODICITE

if(iperio.eq.1) then
  idimte = 21
  itenso = 0
  call percom                                                     &
  !==========
  ( idimte , itenso ,                                             &
    w1     , w1     , w1    ,                                     &
    w2     , w2     , w2    ,                                     &
    w3     , w3     , w3    )
endif


if( imvisf.eq.0 ) then

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    surfn = ra(isrfan-1+ifac)
    dist  = ra(idist -1+ifac)

    viscxi = w1(ii)
    viscxj = w1(jj)
    viscyi = w2(ii)
    viscyj = w2(jj)
    visczi = w3(ii)
    visczj = w3(jj)

    sx2    = surfac(1,ifac)**2
    sy2    = surfac(2,ifac)**2
    sz2    = surfac(3,ifac)**2

    viscf(ifac) = 0.5d0*(                                         &
       (viscxi+viscxj)*sx2                                        &
     + (viscyi+viscyj)*sy2                                        &
     + (visczi+visczj)*sz2 ) / (surfn*dist)

  enddo

else

  do ifac = 1,nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    surfn = ra(isrfan-1+ifac)
    dist  = ra(idist -1+ifac)
    pond  = ra(ipond -1+ifac)

    viscxi = w1(ii)
    viscxj = w1(jj)
    viscyi = w2(ii)
    viscyj = w2(jj)
    visczi = w3(ii)
    visczj = w3(jj)

    sx2    = surfac(1,ifac)**2
    sy2    = surfac(2,ifac)**2
    sz2    = surfac(3,ifac)**2

    viscf(ifac) =                                                 &
      ( viscxi*viscxj*sx2                                         &
              /(pond*viscxi+(1.d0-pond)*viscxj)                   &
      + viscyi*viscyj*sy2                                         &
              /(pond*viscyi+(1.d0-pond)*viscyj)                   &
      + visczi*visczj*sz2                                         &
              /(pond*visczi+(1.d0-pond)*visczj)                   &
       ) /(surfn*dist)
  enddo

endif

do ifac=1,nfabor

  ii = ifabor(ifac)

  surfn = ra(isrfbn-1+ifac)
  dist  = ra(idistb-1+ifac)

  viscxi = w1(ii)
  viscyi = w2(ii)
  visczi = w3(ii)

  sx2    = surfbo(1,ifac)**2
  sy2    = surfbo(2,ifac)**2
  sz2    = surfbo(3,ifac)**2

  viscb(ifac) =                                                   &
    (viscxi*sx2+viscyi*sy2+visczi*sz2)/(surfn*dist)

enddo

!----
! FIN
!----

return

end subroutine
