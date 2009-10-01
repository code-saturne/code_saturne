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

subroutine lwcphy &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ibrom  , izfppp ,          &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE PREMELANGE MODELE LWC
! Calcul de RHO adiabatique ou permeatique (transport de H)


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
! nphmx            ! e  ! <-- ! nphsmx                                         !
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
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!   (nphmx   )     !    !     !                                                !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
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
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
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
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse , nphmx

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr), ibrom(nphmx)
integer          izfppp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          igg, iphas, iel, ipcrom
integer          izone , ifac, ipbrom
double precision coefg(ngazgm)
double precision nbmol , temsmm
double precision masmg

integer       ipass
data          ipass /0/
save          ipass

!===============================================================================

!===============================================================================
! 0. ON COMPTE LES PASSAGES
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0

! ---> Initialisation

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

! ---> Positions des variables, coefficients

iphas = 1
ipcrom = ipproc(irom(iphas))
ipbrom = ipprob(irom(iphas))

!===============================================================================
! 2. DETERMINATION DES GRANDEURS THERMOCHIMIQUES MOYENNES
!===============================================================================


if ( (ippmod(icolwc).eq.0) .or. (ippmod(icolwc).eq.1) ) then

  call pdflwc                                                     &
  !==========
   ( ncelet        , ncel          ,                              &
     rtp(1,isca(ifm))    , rtp(1,isca(ifp2m))  ,                  &
     rtp(1,isca(iyfm))   , rtp(1,isca(iyfp2m)) ,                  &
     propce   )

endif

 if ( (ippmod(icolwc).eq.2) .or. (ippmod(icolwc).eq.3) ) then

  call pdfpp3                                                     &
  !==========
   ( ncelet        , ncel          ,                              &
     rtp(1,isca(ifm))    , rtp(1,isca(ifp2m))  ,                  &
     rtp(1,isca(iyfm))   , rtp(1,isca(iyfp2m)) ,                  &
     rtp(1,isca(icoyfp)) ,                                        &
     propce   )

 endif

 if ( (ippmod(icolwc).eq.4).or.(ippmod(icolwc).eq.5) ) then

  call pdfpp4                                                     &
  !==========
   ( ncelet        , ncel          ,                              &
     rtp(1,isca(ifm))    , rtp(1,isca(ifp2m))  ,                  &
     rtp(1,isca(iyfm))   , rtp(1,isca(iyfp2m)) ,                  &
     rtp(1,isca(icoyfp)) ,                                        &
     propce   )

 endif

!===============================================================================
! 3. CALCUL DE RHO ET DES FRACTIONS MASSIQUES DES ESPECES GLOBALES
!    SUR LES BORDS
!===============================================================================

! --> Masse volumique au bord

iphas = 1
ibrom(iphas) = 1

! ---- Masse volumique au bord pour toutes les facettes
!      Les facettes d'entree seront recalculees.

do ifac = 1, nfabor
iel = ifabor(ifac)
  propfb(ifac,ipbrom) = propce(iel,ipcrom)
enddo

! ---- Masse volumique au bord pour les facettes d'entree UNIQUEMENT
!      Le test sur IZONE sert pour les reprises de calcul

if ( ipass.gt.1 .or. isuite.eq.1 ) then
  do ifac = 1, nfabor
    izone = izfppp(ifac)
    if(izone.gt.0) then
      if ( ientgb(izone).eq.1 .or. ientgf(izone).eq. 1) then
        coefg(1) = fment(izone)
        coefg(2) = 1.d0-fment(izone)
        coefg(3) = zero
        if ( ientgb(izone).eq.1 ) then
          coefg(1) = max(zero,(fment(izone)-fs(1))/(1.d0-fs(1)))
          coefg(3) = (fment(izone)-coefg(1))/fs(1)
          coefg(2) = 1.d0 - coefg(1) - coefg(3)
        endif
        nbmol = 0.d0
        do igg = 1, ngazg
          nbmol = nbmol + coefg(igg)/wmolg(igg)
        enddo
       masmg = 1.d0/nbmol
       temsmm = tkent(izone)/masmg
       propfb(ifac,ipbrom) = p0(iphas)/(rr*temsmm)
      endif
    endif
  enddo
endif

!----
! FIN
!----

return
end subroutine
