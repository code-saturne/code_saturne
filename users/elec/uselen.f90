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
! FONCTION :
! --------

! POUR LA SORTIE POST-TRAITEMENT MODULE ELECTRIQUE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nummai           ! ec ! <-- ! numero du maillage post                        !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncelps           ! e  ! <-- ! nombre de cellules du maillage post            !
! nfacps           ! e  ! <-- ! nombre de faces interieur post                 !
! nfbrps           ! e  ! <-- ! nombre de faces de bord post                   !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! lstcel(ncelps    ! te ! <-- ! liste des cellules du maillage post            !
! lstfac(nfacps    ! te ! <-- ! liste des faces interieures post               !
! lstfbr(nfbrps    ! te ! <-- ! liste des faces de bord post                   !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! tracel(*)        ! tr ! <-- ! tab reel valeurs cellules post                 !
! trafac(*)        ! tr ! <-- ! tab reel valeurs faces int. post               !
! trafbr(*)        ! tr ! <-- ! tab reel valeurs faces bord post               !
! w1-w2            ! tr ! --- ! tab reel pour calcul gradient                  !
! (ncelet,3)       !    !     !                                                !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "dimfbr.h"
include "paramx.h"
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "period.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "elincl.h"

!===============================================================================

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
!       ILPEUT LES AJOUTER A LA FIN, VOIR LA DOCUMENTATION DE USEEVO
!===============================================================================


idebia = idbia0
idebra = idbra0

if(nummai.eq.-1) then

!===============================================================================
! 1.   Gradient du potentiel reel
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

! En periodique et parallele, echange avant calcul du gradient
!     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
  ivar0 = 0

!    Sans prise en compte de la pression hydrostatique

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
   ra     , ra     , ra     ,                                     &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
!       POTR
   w1(1,1) , w1(1,2) , w1(1,3) ,                                  &
!       d POTR /dx   d POTR /dy   d POTR /dz
   w2(1,1) , w2(1,2) , w2(1,3) ,                                  &
   rdevel , rtuser , ra     )

!       Le gradient est defini sur le maillage principal tout entier ;
!       inutile de le recopier, on utilise l'indirection (IVARPR = 1),
!       et les valeurs sont non entrelacees (IENTLA = 0)
  ientla = 0
  ivarpr = 1

  call psteva(nummai, namevr, idimt, ientla, ivarpr,              &
  !==========
              ntcabs, ttcabs, w1, rbid, rbid)

!===============================================================================
! 2.   Gradient du potentiel imaginaire si Joule
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

! En periodique et parallele, echange avant calcul du gradient
!     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
    ivar0 = 0

!    Sans prise en compte de la pression hydrostatique

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

!         Le gradient est defini sur le maillage principal tout entier ;
!         inutile de le recopier, on utilise l'indirection (IVARPR = 1),
!         et les valeurs sont non entrelacees (IENTLA = 0)
    ientla = 0
    ivarpr = 1

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, w1, rbid, rbid)

  endif

!===============================================================================
! 3.   Courant imaginaire si Joule
!===============================================================================

  if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4 ) then

    idimt  = 3
    NAMEVR = 'Cour_Im'

    ivar = isca(ipoti)
    iclimv = iclrtp(ivar,icoef)

!     Commme dans elflux
    ipcsii = ipproc(ivisls(ipoti))

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

! En periodique et parallele, echange avant calcul du gradient
!     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
    ivar0 = 0

!    Sans prise en compte de la pression hydrostatique

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

!         La variable est définie sur le tableau de travail. On
!         a deja utilise l'indirection via LSTCEL (donc IVARPR = 0),
!         et les valeurs sont non entrelacees (IENTLA = 0)
    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

  endif

!===============================================================================
! 5.   Champ Magnetique si Arc
!===============================================================================

  if( ippmod(ielarc).ge.2 ) then

    idimt  = 3
    NAMEVR = 'Ch_Mag'

!    Sur Ax

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

! En periodique et parallele, echange avant calcul du gradient
!     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)

    ivar0 = 0

!    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE

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
   w1(1,1) , w1(1,2) , w1(1,3) ,                                  &
!       d Ax /dx   d Ax /dy   d Ax /dz
   w2(1,1) , w2(1,2) , w2(1,3) ,                                  &
   rdevel , rtuser , ra     )

!       B = rot A

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          =  zero
      tracel(iloc+ncelps)   =  w1(iel,3)
      tracel(iloc+2*ncelps) = -w1(iel,2)
    enddo

!    Sur Ay

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

! En periodique et parallele, echange avant calcul du gradient
!     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)

    ivar0 = 0

!    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE

    iphydp = 0

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

!       B = rot A

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          = tracel(iloc)          - w1(iel,3)
      tracel(iloc+ncelps)   = tracel(iloc + ncelps) + zero
      tracel(iloc+2*ncelps) = tracel(iloc+2*ncelps) + w1(iel,1)
    enddo

!    Sur Az

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

! En periodique et parallele, echange avant calcul du gradient
!     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)

    ivar0 = 0

!    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE

    iphydp = 0

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

!       B = rot A

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc)          = tracel(iloc)          + w1(iel,2)
      tracel(iloc+ncelps)   = tracel(iloc+ncelps)   - w1(iel,1)
      tracel(iloc+2*ncelps) = tracel(iloc+2*ncelps) + zero
    enddo

!         La variable est définie sur le tableau de travail. On
!         a deja utilise l'indirection via LSTCEL (donc IVARPR = 0),
!         et les valeurs sont non entrelacees (IENTLA = 0)
    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

  endif


!===============================================================================
! 4.   Module et Argument du potentiel si IELJOU = 4
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
