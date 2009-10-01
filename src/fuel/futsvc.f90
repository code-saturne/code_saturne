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

subroutine futsvc &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iscal  , iscala ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   smbrs  , rovsdt ,                                              &
   wfb    ,                                                       &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME FUEL
!   TERMES SOURCES DE PRODUCTION ET DE DISSIPATION POUR
!   LA VARIANCE (BILANS EXPLICITE ET IMPLICITE)

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
! ncesmp           ! e  ! <-- ! nombre de cellules a source de masse           !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! iscal            ! e  ! <-- ! numero du scalaire                             !
! iscala           ! e  ! <-- ! numero du scalaire associe                     !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
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
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
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
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
! wfb(nfabor)      ! tr ! --- ! tableau de travail    faces de bord            !
! w1..8(ncelet)    ! tr ! --- ! tableau de travail    cellules                 !
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
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "parall.h"
include "period.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "fuincl.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse
integer          iscal  , iscala

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml) , itypfb(nfabor,nphas)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          idevel(nideve)
integer          ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision smbrs(ncelet), rovsdt(ncelet)
double precision wfb(nfabor)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia , idebra
integer          ivar   , ivarsc , ivarut , ivar0
integer          iel    , iphas  , ifac   , icla
integer          ipcrom , ipcvst
integer          ikiph  , ieiph , iomgip , iphydp
integer          ir11ip , ir22ip, ir33ip
integer          inc    , iccocg , nswrgp , imligp , iwarnp
integer          ifinra , icoefa , icoefb
integer          idimte , itenso

double precision xk     , xe     , rhovst
double precision epsrgp , climgp , extrap

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar   = isca(iscal)

! --- Numero du scalaire eventuel associe dans le cas fluctuation
!         ISCALA et numero de variable de calcul
if (iscala.gt.0) then
  ivarsc = isca(iscala)
else
  ivarsc = 0
endif

! --- Numero de phase associee au scalaire ISCAL
iphas = iphsca(iscal)

! --- Numero des variables de calcul
if ( itytur(iphas).eq.2 .or. iturb(iphas).eq.50 ) then
  ikiph  = ik  (iphas)
  ieiph  = iep (iphas)
elseif ( itytur(iphas).eq.3 ) then
  ir11ip = ir11(iphas)
  ir22ip = ir22(iphas)
  ir33ip = ir33(iphas)
  ieiph  = iep (iphas)
elseif ( iturb(iphas).eq.60 ) then
  ikiph  = ik  (iphas)
  iomgip = iomg(iphas)
endif

! --- Numero des grandeurs physiques
ipcrom = ipproc(irom(iphas))
ipcvst = ipproc(ivisct(iphas))


!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES DE PRODUCTION PAR LES GRADIENTS
!    ET DE DISSIPATION
!===============================================================================

if ( itytur(iphas).eq.2 .or. itytur(iphas).eq.3                   &
     .or. iturb(iphas).eq.50 .or. iturb(iphas).eq.60) then

  inc = 1
  iccocg = 1
  if (ivarsc.gt.0) then
    ivarut = ivarsc
  else
! A defaut de savoir pour F4M on prend comme pour IFVAP
    ivarut = isca(ifvap)
  endif
  nswrgp = nswrgr(ivarut)
  imligp = imligr(ivarut)
  iwarnp = iwarni(ivarut)
  epsrgp = epsrgr(ivarut)
  climgp = climgr(ivarut)
  extrap = extrag(ivarut)

! --> Calcul de FIM et X2 dans W7 et W8

  do iel = 1, ncel
    w1(iel) = zero
    w2(iel) = zero
    w3(iel) = zero
    w4(iel) = zero
    w5(iel) = zero
    w6(iel) = zero
    w7(iel) = zero
    w8(iel) = 1.d0
  enddo

! ---- W8 = X1

  do iel = 1, ncel
    w8(iel) = 1.d0
  enddo
  do icla = 1, nclafu
    do iel = 1, ncel
      w8(iel) = w8(iel) - rtp(iel,isca(iyfol(icla)))
    enddo
  enddo

! ---- W7 = FJM (kg/kg du melange gazeux)

  if (ivarsc.eq.0) then
    do iel = 1, ncel
      w7(iel) = 1.d0                                              &
               -( rtp(iel,isca(ifhtf))                            &
                 +rtp(iel,isca(ifvap)) ) / w8(iel)

    enddo
  else
    do iel = 1, ncel
      w7(iel) = rtp(iel,ivarsc) / w8(iel)
    enddo
  endif

! --> Calcul des COEFA et COEFB de FIM afin d'en calculer son gradient
!     On alloue localement 2 tableaux de NFABOR pour le calcul
!       de COEFA et COEFB de FIM

  icoefa = idebra
  icoefb = icoefa + nfabor
  ifinra = icoefb + nfabor
  CALL RASIZE ('FUTSVC',IFINRA)
  !==========

  do ifac = 1, nfabor
    ra(icoefa+ifac-1) = zero
    ra(icoefb+ifac-1) = 1.d0
    if ( itypfb(ifac,iphas).eq.ientre ) then
      ra(icoefa+ifac-1) = zero
      ra(icoefb+ifac-1) = zero
      if (ivarsc.eq.0) ra(icoefa+ifac-1) = 1.d0
    endif
  enddo

! En periodique et parallele, echange avant calcul du gradient

!    Parallele
  if(irangp.ge.0) then
    call parcom(w7)
    !==========
  endif

!    Periodique
  if(iperio.eq.1) then
    idimte = 0
    itenso = 0
    call percom                                                   &
    !==========
  ( idimte , itenso ,                                             &
    w7     , w7     , w7     ,                                    &
    w7     , w7     , w7     ,                                    &
    w7     , w7     , w7     )
  endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
  ivar0  = 0
  iphydp = 0
  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w7     , w7     , w7     ,                                     &
   w7     , ra(icoefa) , ra(icoefb)  ,                            &
!          FIM      COEFA        COEFB
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

  do iel = 1, ncel
    if ( itytur(iphas).eq.2 .or. iturb(iphas).eq.50 ) then
      xk = rtpa(iel,ikiph)
      xe = rtpa(iel,ieiph)
    elseif ( itytur(iphas).eq.3 ) then
      xk =                                                        &
       0.5d0*(rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))
      xe = rtpa(iel,ieiph)
    elseif ( iturb(iphas).eq.60 ) then
      xk = rtpa(iel,ikiph)
      xe = cmu*xk*rtpa(iel,iomgip)
    endif

    rhovst = propce(iel,ipcrom)*xe/                               &
             (xk * rvarfl(iscal))*volume(iel)
    rovsdt(iel) = rovsdt(iel) + max(zero,rhovst)
    smbrs (iel) = smbrs (iel) +                                   &
                2.d0*propce(iel,ipcvst)*volume(iel)/sigmas(iscal) &
                * (w1(iel)**2 + w2(iel)**2 + w3(iel)**2) * w8(iel)&
                - rhovst*rtpa(iel,ivar)

  enddo

endif

!--------
! FORMATS
!--------



!----
! FIN
!----

return

end subroutine
