!-------------------------------------------------------------------------------

!VERS


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

subroutine uselen &
!================

 ( idbia0 , idbra0 , nummai ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncelps , nfacps , nfbrps ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   lstcel , lstfac , lstfbr ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     ,                                              &
   tracel , trafac , trafbr , rdevel , rtuser , ra     )

!===============================================================================
! Purpose :
! --------

! For post-processing in electric module

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! nummai           ! ec ! <-- ! numero du maillage post                        !
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
! ncelps           ! e  ! <-- ! nombre de cellules du maillage post            !
! nfacps           ! e  ! <-- ! nombre de faces interieur post                 !
! nfbrps           ! e  ! <-- ! nombre de faces de bord post                   !
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
! lstcel(ncelps    ! te ! <-- ! liste des cellules du maillage post            !
! lstfac(nfacps    ! te ! <-- ! liste des faces interieures post               !
! lstfbr(nfbrps    ! te ! <-- ! liste des faces de bord post                   !
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
! tracel(*)        ! tr ! <-- ! tab reel valeurs cellules post                 !
! trafac(*)        ! tr ! <-- ! tab reel valeurs faces int. post               !
! trafbr(*)        ! tr ! <-- ! tab reel valeurs faces bord post               !
! w1-w2            ! tr ! --- ! tab reel pour calcul gradient                  !
! (ncelet,3)       !    !     !                                                !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use elincl

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nummai
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncelps , nfacps , nfbrps
integer          nideve , nrdeve , nituse , nrtuse
integer          idimt

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision tracel(ncelps*3)
double precision trafac(nfacps*3), trafbr(nfbrps*3)
double precision w1(ncelet,3), w2(ncelet,3)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

character*32     namevr
integer          idebia, idebra, iel   , iloc
integer          ivar  , ivar0 , inc   , iccocg
integer          iphydp, nswrgp, imligp, iwarnp, iclimv
integer          ipcsii
integer          ientla, ivarpr
double precision epsrgp, climgp, extrap
double precision rbid(1)

!===============================================================================
!===============================================================================
! 0.  PAR DEFAUT, ON CONSIDERE QUE LE SOUS PROGRAMME CI-DESSOUS CONVIENT
!       A L'UTILISATEUR, C'EST-A-DIRE QUE LA MISE EN OEUVRE DU MODULE
!       ELECTRIQUE DECLENCHE LA PRODUCTION DE CHAMPS STANDARD DANS LE
!       POST-TRAITEMENT.
!     L'UTILISATEUR N'A PAS A MODIFIER LE PRESENT SOUS-PROGRAMME DANS
!       LES CONDITIONS D'UTILISATION STANDARD.
!     DANS LE CAS OU IL SOUHAITE PRODUIRE DES VARIABLES SUPPLEMENTAIRES
!       IL PEUT LES AJOUTER A LA FIN, VOIR LA DOCUMENTATION DE USEEVO
!===============================================================================


idebia = idbia0
idebra = idbra0

if(nummai.eq.-1) then

!===============================================================================
! 1.   Graident of the real potential
!===============================================================================

  idimt  = 3
  NAMEVR = 'Gr_PotR'

  ivar = isca(ipotr)
  iclimv = iclrtp(ivar,icoef)

  inc = 1
  iccocg = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  ivar0 = 0
  iphydp = 0
!
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
   ra     , ra     , ra     ,                                     &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
!       POTR
   w1(1,1) , w1(1,2) , w1(1,3) ,                                  &
!       d POTR /dx   d POTR /dy   d POTR /dz
   w2(1,1) , w2(1,2) , w2(1,3) ,                                  &
   rdevel , rtuser , ra     )

!
  ientla = 0
  ivarpr = 1

  call psteva(nummai, namevr, idimt, ientla, ivarpr,              &
  !==========
              ntcabs, ttcabs, w1, rbid, rbid)

!===============================================================================
! 2.   For Joule Heating by direct conduction :
!                           gradient of the imaginary component of the potential
!===============================================================================

  if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then

    idimt  = 3
    NAMEVR = 'Gr_PotI'

    ivar = isca(ipoti)
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
!
    ivar0 = 0
    iphydp = 0
!
    call grdcel                                                   &
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
   ra     , ra     , ra     ,                                     &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
!       POTI
   w1(1,1) , w1(1,2) , w1(1,3) ,                                  &
!       d POTI /dx   d POTI /dy   d POTI /dz
   w2(1,1) , w2(1,2) , w2(1,3) ,                                  &
   rdevel , rtuser , ra     )

!
    ientla = 0
    ivarpr = 1

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, w1, rbid, rbid)

  endif

!===============================================================================
! 3.  For Joule heating by direct conduction :
!                                     imaginary component of the current density
!===============================================================================

  if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4 ) then

    idimt  = 3
    NAMEVR = 'Cour_Im'

    ivar = isca(ipoti)
    iclimv = iclrtp(ivar,icoef)

!    As in elflux
    ipcsii = ipproc(ivisls(ipoti))

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
!
    ivar0 = 0
    iphydp = 0

    call grdcel                                                   &
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
   ra     , ra     , ra     ,                                     &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
!       POTI
   w1(1,1) , w1(1,2) , w1(1,3) ,                                  &
!       d POTI /dx   d POTI /dy   d POTI /dz
   w2(1,1) , w2(1,2) , w2(1,3) ,                                  &
   rdevel , rtuser , ra     )

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          = -propce(iel,ipcsii)*w1(iel,1)
      tracel(iloc+ncelps)   = -propce(iel,ipcsii)*w1(iel,2)
      tracel(iloc+2*ncelps) = -propce(iel,ipcsii)*w1(iel,3)
    enddo
!
    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

  endif

!==========================================================
! 5.   For electric arc : electromagnetic field calculation
!==========================================================

  if( ippmod(ielarc).ge.2 ) then

    idimt  = 3
    NAMEVR = 'Ch_Mag'

!   Ax Component

    ivar = isca(ipotva(1))
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
!
    ivar0 = 0
    iphydp = 0
!
    call grdcel                                                   &
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
   ra     , ra     , ra     ,                                     &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
   w1(1,1) , w1(1,2) , w1(1,3) ,                                  &
!       d Ax /dx   d Ax /dy   d Ax /dz
   w2(1,1) , w2(1,2) , w2(1,3) ,                                  &
   rdevel , rtuser , ra     )

!       B = rot A ( B = curl A)

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          =  zero
      tracel(iloc+ncelps)   =  w1(iel,3)
      tracel(iloc+2*ncelps) = -w1(iel,2)
    enddo

!    Ay component

    ivar = isca(ipotva(2))
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
!
    ivar0 = 0
    iphydp = 0
!
    call grdcel                                                   &
    !==========
  ( idbia0 , idbra0 ,                                             &
    ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml ,&
    nnod   , lndfac , lndfbr , ncelbr , nphas  ,                  &
    nideve , nrdeve , nituse , nrtuse ,                           &
    ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp ,&
    iwarnp , nfecra , epsrgp , climgp , extrap ,                  &
    ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                  &
    ipnfac , nodfac , ipnfbr , nodfbr ,                           &
    idevel , ituser , ia     ,                                    &
    xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume ,&
    ra     , ra     , ra     ,                                    &
    rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv) ,              &
    w1(1,1) , w1(1,2) , w1(1,3) ,                                 &
!       d Ay /dx   d Ay /dy   d Ay /dz
    w2(1,1) , w2(1,2) , w2(1,3) ,                                 &
    rdevel , rtuser , ra     )

!       B = rot A (B = curl A)

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          = tracel(iloc)          - w1(iel,3)
      tracel(iloc+ncelps)   = tracel(iloc + ncelps) + zero
      tracel(iloc+2*ncelps) = tracel(iloc+2*ncelps) + w1(iel,1)
    enddo

!    Az component

    ivar = isca(ipotva(3))
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)
!
    ivar0 = 0
   iphydp = 0
!
    call grdcel                                                   &
    !==========
  ( idbia0 , idbra0 ,                                             &
    ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml ,&
    nnod   , lndfac , lndfbr , ncelbr , nphas  ,                  &
    nideve , nrdeve , nituse , nrtuse ,                           &
    ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp ,&
    iwarnp , nfecra , epsrgp , climgp , extrap ,                  &
    ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                  &
    ipnfac , nodfac , ipnfbr , nodfbr ,                           &
    idevel , ituser , ia     ,                                    &
    xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume ,&
    ra     , ra     , ra     ,                                    &
    rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv) ,              &
    w1(1,1) , w1(1,2) , w1(1,3) ,                                 &
!       d Az /dx   d Az /dy   d Az /dz
    w2(1,1) , w2(1,2) , w2(1,3) ,                                 &
    rdevel , rtuser , ra     )

!       B = rot A (B = curl A)

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          = tracel(iloc)          + w1(iel,2)
      tracel(iloc+ncelps)   = tracel(iloc+ncelps)   - w1(iel,1)
      tracel(iloc+2*ncelps) = tracel(iloc+2*ncelps) + zero
    enddo
!
    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

  endif


!===============================================================================
! 4.   Calculation of Module and Argument of the complex potential if IELJOU = 4
!===============================================================================

  if (ippmod(ieljou).eq.4) then

    idimt  = 1
    NAMEVR = 'ModPot'

    ivar = isca(ipotr)

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc) =                                              &
        sqrt( rtp(iel,isca(ipotr))*rtp(iel,isca(ipotr))           &
             +rtp(iel,isca(ipoti))*rtp(iel,isca(ipoti)) )
    enddo

    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

    idimt  = 1
    NAMEVR = 'ArgPot'

    ivar = isca(ipotr)

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      if ( rtp(iel,isca(ipotr)) .ne. 0.d0 ) then
        if ( rtp(iel,isca(ipotr)) .ge. 0.d0 ) then
          tracel(iloc) =                                          &
           atan( rtp(iel,isca(ipoti))/rtp(iel,isca(ipotr)))
        else
          if ( rtp(iel,isca(ipoti)) .gt. 0.d0 ) then
            tracel(iloc) =                                        &
              4.d0*atan(1.d0)                                     &
             +atan( rtp(iel,isca(ipoti))                          &
                   /rtp(iel,isca(ipotr)))
          else
            tracel(iloc) =                                        &
             -4.d0*atan(1.d0)                                     &
             +atan( rtp(iel,isca(ipoti))                          &
                   /rtp(iel,isca(ipotr)))
          endif
        endif
      else
        tracel(iloc) = 2.d0*atan(1.d0)
      endif

      if (tracel(iloc) .lt. 0.d0) then
        tracel(iloc) = tracel(iloc) + 8.d0**atan(1.d0)
      endif

    enddo

    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

  endif

endif

return

end subroutine
