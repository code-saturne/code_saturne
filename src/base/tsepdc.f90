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

subroutine tsepdc &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp ,                            &
   nideve , nrdeve , nituse , nrtuse , iphas  , idiaex ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , trav   ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DES TERMES DE PERTE DE CHARGE
!        POUR LE BILAN EXPLICITE
!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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
! ncepdp           ! e  ! <-- ! nombre de cellules avec pdc                    !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! iphas            ! e  ! <-- ! numero de phase courante                       !
! idiaex           ! e  ! <-- ! indicateur de traitement de la                 !
!                  !    !     ! diagonale (=1) ou extradiagonale (=2)          !
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
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
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
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant prec)                     !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
! trav(ncelet,3    ! tr ! <-- ! tableau des second membres                     !
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
include "numvar.h"
include "optcal.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp
integer          nideve , nrdeve , nituse , nrtuse , iphas
integer          idiaex

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision trav(ncelet,3)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          iel   , ielpdc
integer          iuiph , iviph , iwiph , ipcrom, ipcroo
double precision romvom, vit1  , vit2  , vit3
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23

!===============================================================================

idebia = idbia0
idebra = idbra0

iuiph  = iu(iphas)
iviph  = iv(iphas)
iwiph  = iw(iphas)
ipcrom = ipproc(irom  (iphas))

ipcroo = ipcrom
if(iroext(iphas).gt.0.and.isno2t(iphas).gt.0) then
  ipcroo = ipproc(iroma (iphas))
endif

!     La diagonale est toujours "implicite"

if(idiaex.eq.1) then

  do ielpdc = 1, ncepdp

    iel    = icepdc(ielpdc)
    romvom =-propce(iel,ipcrom)*volume(iel)
    cpdc11 = ckupdc(ielpdc,1)
    cpdc22 = ckupdc(ielpdc,2)
    cpdc33 = ckupdc(ielpdc,3)
    vit1   = rtpa(iel,iuiph)
    vit2   = rtpa(iel,iviph)
    vit3   = rtpa(iel,iwiph)

    trav(iel,1) = trav(iel,1) +                                   &
         romvom * ( cpdc11*vit1                             )
    trav(iel,2) = trav(iel,2) +                                   &
         romvom * (               cpdc22*vit2               )
    trav(iel,3) = trav(iel,3) +                                   &
         romvom * (                             cpdc33*vit3 )

  enddo

endif

!     L'extradiagonale est explicite mais peut etre extrapolee

if(idiaex.eq.2) then

  do ielpdc = 1, ncepdp

    iel    = icepdc(ielpdc)
    romvom =-propce(iel,ipcroo)*volume(iel)
    cpdc12 = ckupdc(ielpdc,4)
    cpdc13 = ckupdc(ielpdc,5)
    cpdc23 = ckupdc(ielpdc,6)
    vit1   = rtpa(iel,iuiph)
    vit2   = rtpa(iel,iviph)
    vit3   = rtpa(iel,iwiph)

    trav(iel,1) = trav(iel,1) +                                   &
         romvom * (               cpdc12*vit2 + cpdc13*vit3 )
    trav(iel,2) = trav(iel,2) +                                   &
         romvom * ( cpdc12*vit1               + cpdc23*vit3 )
    trav(iel,3) = trav(iel,3) +                                   &
         romvom * ( cpdc13*vit1 + cpdc23*vit2               )

  enddo

endif


return

end subroutine
