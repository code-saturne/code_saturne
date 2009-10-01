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

subroutine strpre &
!================

 ( idbia0 , idbra0 , itrale , italim , ineefl ,                   &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   impale ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   flmalf , flmalb , xprale , cofale , depale , rdevel , rtuser , &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! PREDICTION DU DEPLACEMENT DES STRUCTURES EN ALE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! itrale           ! e  ! <-- ! numero d'iteration pour l'ale                  !
! italim           ! e  ! <-- ! numero d'iteration couplage implicite          !
! ineedf           ! e  ! <-- ! indicateur de sauvegarde des flux              !
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
!   (nfac+1)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (lndfac)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfabor+1)     !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (lndfbr)       !    !     !  (optionnel)                                   !
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
! flmalf(nfac)     ! tr ! --> ! sauvegarde du flux de masse faces int          !
! flmalb(nfabor    ! tr ! --> ! sauvegarde du flux de masse faces brd          !
! cofale           ! tr ! --> ! sauvegarde des cl de p et u                    !
!    (nfabor,8)    !    !     !                                                !
! xprale(ncelet    ! tr ! --> ! sauvegarde de la pression, si nterup           !
!                  !    !     !    est >1                                      !
! depale(nnod,3    ! tr ! <-- ! deplacement aux noeuds                         !
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

include "dimfbr.h"
include "paramx.h"
include "optcal.h"
include "numvar.h"
include "pointe.h"
include "albase.h"
include "alstru.h"
include "alaste.h"
include "period.h"
include "parall.h"
include "entsor.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0 , itrale , italim , ineefl
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          impale(nnod)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision flmalf(nfac), flmalb(nfabor), xprale(ncelet)
double precision cofale(nfabor,8)
double precision depale(nnod,3)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra, ifinia
integer          istr, ii, ifac, inod, iel, indast
integer          iflmas, iflmab,iclp,iclu,iclv,iclw
integer          ilstfa

!===============================================================================


!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

iflmas = ipprof(ifluma(iu(1)))
iflmab = ipprob(ifluma(iu(1)))
iclp = iclrtp(ipr(1),icoef)
iclu = iclrtp(iu(1),icoef)
iclv = iclrtp(iv(1),icoef)
iclw = iclrtp(iw(1),icoef)

!===============================================================================
! 2. PREDICTION DU DEPLACEMENT DES STRUCTURES
!===============================================================================

! 2.1 STRUCTURES INTERNES :
! -----------------------

! Lors de la phase d'initialisation de l'ALE (ITRALE=0), XSTP contient
!    - la valeur du deplacement initial des structures si l'utilisateur
!        les a touchees (suite ou pas)
!    - 0 si debut de calcul avec les structures
!    - le deplacement utilise pour le calcul precedent si suite sans
!        modification utilisateur
!   Il faut cependant transferer sa valeur dans XSTR (qui est utilise
!   par Newmark)
! Lors des iterations suivantes (ITRALE>0) on utilise les schemas standard
!   de calcul de XSTP

if (nbstru.gt.0) then

  if (itrale.eq.0) then

    do istr = 1, nbstru
      do ii = 1, ndim
        xstr(ii,istr) = xstp(ii,istr)
      enddo
    enddo

  else

! 2.1.1 : SCHEMA DE COUPLAGE EXPLICITE
!---------------------------------------------
    if (nalimx.eq.1) then
      do istr = 1, nbstru
        do ii = 1, 3
          xstp(ii,istr) = xstr(ii,istr)                           &
             + aexxst*dtstr(istr)*xpstr(ii,istr)                  &
             + bexxst*dtstr(istr)*(xpstr(ii,istr)-xpsta(ii,istr))
        enddo
      enddo

! 2.1.2 : SCHEMA DE COUPLAGE IMPLICITE
!---------------------------------------------
    else
      do istr = 1, nbstru
        do ii = 1, 3
          xstp(ii,istr)  = xstr(ii,istr)
        enddo
      enddo
    endif

  endif

  do ifac = 1, nfabor
    istr = ia(iidfst+ifac-1)
    if (istr.gt.0) then
      do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
        inod = nodfbr(ii)
        impale(inod) = 1
        depale(inod,1) = xstp(1,istr)
        depale(inod,2) = xstp(2,istr)
        depale(inod,3) = xstp(3,istr)
      enddo
    endif
  enddo

endif


! 2.2 STRUCTURES EXTERNES (COUPLAGE CODE_ASTER) :
! -----------------------


if (nbaste.gt.0) then

  do ifac = 1, nfabor
    istr = ia(iidfst+ifac-1)
    if (istr.lt.0) then
      do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
        inod = nodfbr(ii)
        impale(inod) = 1
      enddo
    endif
  enddo

! Si ITRALE = 0, on ne fait rien pour l'instant, mais il faudrait
! prevoir la reception de deplacement initiaux venant de Code_Aster
  if (itrale.gt.0) then

    ntcast = ntcast + 1

! Reception des deplacements predits et remplissage de depale

    ilstfa = idebia
    ifinia = ilstfa + nbfast
    CALL IASIZE('STRPRE',IFINIA)

    indast = 0
    do ifac = 1, nfabor
      istr = ia(iidfst+ifac-1)
      if (istr.lt.0) then
        indast = indast + 1
        ia(ilstfa + indast-1) = ifac
      endif
    enddo

    call astcin(ntcast, nbfast, ia(ilstfa), depale)
    !==========

  endif

endif

!===============================================================================
! 3. DEPLACEMENT AU PAS DE TEMPS PRECEDENT ET SAUVEGARDE FLUX ET DE P
!===============================================================================

if (italim.eq.1) then
  do istr = 1, nbstru
    do ii = 1, 3
      xsta(ii,istr)   = xstr(ii,istr)
      xpsta(ii,istr)  = xpstr(ii,istr)
      xppsta(ii,istr) = xppstr(ii,istr)
    enddo
  enddo
  if (ineefl.eq.1) then
    do ifac = 1, nfac
      flmalf(ifac) = propfa(ifac,iflmas)
    enddo
    do ifac = 1, nfabor
      flmalb(ifac) = propfb(ifac,iflmab)
      cofale(ifac,1) = coefa(ifac,iclp)
      cofale(ifac,2) = coefa(ifac,iclu)
      cofale(ifac,3) = coefa(ifac,iclv)
      cofale(ifac,4) = coefa(ifac,iclw)
      cofale(ifac,5) = coefb(ifac,iclp)
      cofale(ifac,6) = coefb(ifac,iclu)
      cofale(ifac,7) = coefb(ifac,iclv)
      cofale(ifac,8) = coefb(ifac,iclw)
    enddo
    if (nterup.gt.1) then
      do iel = 1, ncelet
        xprale(iel) = rtpa(iel,ipr(1))
      enddo
    endif
  endif
endif

!----
! FORMATS
!----



!----
! FIN
!----

end subroutine
