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

subroutine attycl &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbmetd , nbmett , nbmetm ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , itrifb , itypfb , izfppp , iprofm ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   tmprom , ztprom , zdprom , xmet   , ymet   , pmer   ,          &
   ttprom , qvprom , uprom  , vprom  , ekprom , epprom ,          &
   rprom  , tpprom , phprom ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!    CONDITIONS AUX LIMITES AUTOMATIQUES

!           ECOULEMENTS ATMOSPHERIQUES


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
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
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
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! itrifb(nfabor    ! te ! <-- ! indirection pour tri des faces de brd          !
!  nphas      )    !    !     !                                                !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! w1,2,3,4,5,6     ! tr ! --- ! tableaux de travail                            !
!  (ncelet         !    !     !  (calcul du gradient de pression)              !
! coefu            ! tr ! --- ! tab de trav                                    !
!  (nfabor,3)      !    !     !  (calcul du gradient de pression)              !
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

! Arguments

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbmetd , nbmett , nbmetm
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          izfppp(nfabor), iprofm(nozppm)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision tmprom(nbmetm)
double precision ztprom(nbmett) , zdprom(nbmetd)
double precision xmet(nbmetm)   , ymet(nbmetm)  , pmer(nbmetm)
double precision ttprom(nbmett,nbmetm) , qvprom(nbmett,nbmetm)
double precision uprom(nbmetd,nbmetm)  , vprom(nbmetd,nbmetm)
double precision ekprom(nbmetd,nbmetm) , epprom(nbmetd,nbmetm)
double precision rprom(nbmett,nbmetm)  , tpprom(nbmett,nbmetm)
double precision phprom(nbmett,nbmetm)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,ndim)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ifac, izone, iphas
double precision d2s3, zent, vs, xuent, xvent
double precision xkent, xeent, tpent

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

d2s3 = 2.d0/3.d0

xuent = 0.d0
xvent = 0.d0
xkent = 0.d0
xeent = 0.d0
tpent = 0.d0

!===============================================================================
! 2.  SI IPROFM = 1 : CHOIX ENTREE/SORTIE SUIVANT LE PROFIL METEO SI
!                       ITYPFB N'A PAS ETE MODIFIE
!                     VARIABLES TIREES DU PROFIL METEO SI
!                       RCODCL(IFAC,IVAR,1) N'A PAS ETE MODIFIE

!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!===============================================================================


do ifac = 1, nfabor

  izone = izfppp(ifac)

  if (iprofm(izone).eq.1) then

!     On recupere les valeurs du profil et on met a jour RCODCL s'il n'a pas
!       ete modifie. Il servira si la face est une face d'entree ou si c'est une
!       face de sortie (si le flux est rentrant).
    zent=cdgfbo(3,ifac)

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, uprom , zent  , ttcabs, xuent )

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, vprom , zent  , ttcabs, xvent )

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, ekprom, zent  , ttcabs, xkent )

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, epprom, zent  , ttcabs, xeent )

    call intprf                                                   &
    !==========
   (nbmett, nbmetm,                                               &
    ztprom, tmprom, tpprom, zent  , ttcabs, tpent )
!
    vs = xuent*surfbo(1,ifac) + xvent*surfbo(2,ifac)

    do iphas = 1, nphas


!     On met a jour le type de face de bord s'il n'a pas ete specifie
!       par l'utilisateur.
!     Pour une entree, on remplit la condition de Dirichlet si elle n'a pas
!     ete  specifiee par utilisateur.

      if (vs.gt.0) then
        if (itypfb(ifac,iphas).eq.0) itypfb(ifac,iphas) = isolib
      else
        if (itypfb(ifac,iphas).eq.0) itypfb(ifac,iphas) = ientre


        if (rcodcl(ifac,iu(iphas),1).gt.rinfin*0.5d0)             &
           rcodcl(ifac,iu(iphas),1) = xuent
        if (rcodcl(ifac,iv(iphas),1).gt.rinfin*0.5d0)             &
           rcodcl(ifac,iv(iphas),1) = xvent
        if (rcodcl(ifac,iw(iphas),1).gt.rinfin*0.5d0)             &
           rcodcl(ifac,iw(iphas),1) = 0.d0

        if    (itytur(iphas).eq.2) then

          if (rcodcl(ifac,ik(iphas),1).gt.rinfin*0.5d0)           &
             rcodcl(ifac,ik(iphas),1) = xkent
          if (rcodcl(ifac,iep(iphas),1).gt.rinfin*0.5d0)          &
             rcodcl(ifac,iep(iphas),1) = xeent

        elseif(itytur(iphas).eq.3) then

          if (rcodcl(ifac,ir11(iphas),1).gt.rinfin*0.5d0)         &
             rcodcl(ifac,ir11(iphas),1) = d2s3*xkent
          if (rcodcl(ifac,ir22(iphas),1).gt.rinfin*0.5d0)         &
             rcodcl(ifac,ir22(iphas),1) = d2s3*xkent
          if (rcodcl(ifac,ir33(iphas),1).gt.rinfin*0.5d0)         &
             rcodcl(ifac,ir33(iphas),1) = d2s3*xkent
          if (rcodcl(ifac,ir12(iphas),1).gt.rinfin*0.5d0)         &
             rcodcl(ifac,ir12(iphas),1) = 0.d0
          if (rcodcl(ifac,ir13(iphas),1).gt.rinfin*0.5d0)         &
             rcodcl(ifac,ir13(iphas),1) = 0.d0
          if (rcodcl(ifac,ir23(iphas),1).gt.rinfin*0.5d0)         &
             rcodcl(ifac,ir23(iphas),1) = 0.d0
          if (rcodcl(ifac,iep(iphas),1).gt.rinfin*0.5d0)          &
             rcodcl(ifac,iep(iphas),1) = xeent

        elseif(iturb(iphas).eq.50) then

          if (rcodcl(ifac,ik(iphas),1).gt.rinfin*0.5d0)           &
             rcodcl(ifac,ik(iphas),1) = xkent
          if (rcodcl(ifac,iep(iphas),1).gt.rinfin*0.5d0)          &
             rcodcl(ifac,iep(iphas),1) = xeent
          if (rcodcl(ifac,iphi(iphas),1).gt.rinfin*0.5d0)         &
             rcodcl(ifac,iphi(iphas),1) = d2s3
          if (rcodcl(ifac,ifb(iphas),1).gt.rinfin*0.5d0)          &
             rcodcl(ifac,ifb(iphas),1) = 0.d0

        elseif(iturb(iphas).eq.60) then

          if (rcodcl(ifac,ik(iphas),1).gt.rinfin*0.5d0)           &
             rcodcl(ifac,ik(iphas),1) = xkent
          if (rcodcl(ifac,iomg(iphas),1).gt.rinfin*0.5d0)         &
             rcodcl(ifac,iomg(iphas),1) = xeent/cmu/xkent

        endif

        if (iscalt(iphas).ne.-1) then

          if (rcodcl(ifac,isca(iscalt(iphas)),1).gt.rinfin*0.5d0) &
           rcodcl(ifac,isca(iscalt(iphas)),1) = tpent

        endif

      endif

    enddo

  endif

enddo

!----
! FORMATS
!----


!----
! FIN
!----

return
end
