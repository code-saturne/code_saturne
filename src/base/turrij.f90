!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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

subroutine turrij &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iphas  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   tslagr ,                                                       &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  , coefax ,                                     &
   dam    , xam    , drtp   ,                                     &
   smbr   , rovsdt , grdvit , produc , grarox , graroy , graroz , &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     ,                   &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS Rij-EPS 1 PHASE INCOMPRESSIBLE OU
! RHO VARIABLE SUR UN PAS DE TEMPS

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
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
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iphas            ! i  ! <-- ! phase number                                   !
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
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
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
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! coefax(nfabor    ! tr ! --- ! tab de trav pour cond.lim. paroi               !
!                  ! tr ! --- !   attention : uniquement avec echo             !
!                  ! tr ! --- !   de paroi et abs(icdpar) = 1                  !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr (ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! drtp(ncelet)     ! tr ! --- ! tableau de travail pour increment              !
! smbr?(ncelet)    ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! grdvit           ! tr ! --- ! tableau de travail pour terme grad             !
!  (ncelet,3,3)    !    !     !    de vitesse     uniqt pour iturb=31          !
! produc           ! tr ! <-- ! tableau de travail pour production             !
!  (6,ncelet)      !    !     ! (sans rho volume) uniqt pour iturb=30          !
! grarox,y,z       ! tr ! --- ! tableau de travail pour grad rom               !
!  (ncelet)        !    !     !                                                !
! w?(ncelet)       ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "dimfbr.h"
include "paramx.h"
include "numvar.h"
include "entsor.h"
include "cstphy.h"
include "optcal.h"
include "lagpar.h"
include "lagran.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse , iphas

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp)
integer          itypsm(*) ! mapped to (ncesmp,nvar)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6)
double precision viscf(nfac), viscb(nfabor), coefax(nfabor)
double precision smacel(*) ! mapped to (ncesmp,nvar)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbr(ncelet), rovsdt(ncelet)
double precision grdvit(ncelet,3,3), produc(6,ncelet)
double precision grarox(ncelet), graroy(ncelet), graroz(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ifac  , iel   , ivar  , isou  , ii
integer          inc   , iccocg
integer          ipp   , iwarnp, iclip
integer          ipriph, iuiph , iviph , iwiph
integer          ir11ip, ir22ip, ir33ip, ir12ip, ir13ip, ir23ip
integer          ieiph
integer          icliup, iclivp, icliwp
integer          nswrgp, imligp, iphydp
integer          ipcrom, ipbrom, ipcroo, ipbroo, iivar
integer          iitsla, iitsli, iicsmp, iigamm
double precision epsrgp, climgp, extrap

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

ipriph = ipr (iphas)
iuiph  = iu  (iphas)
iviph  = iv  (iphas)
iwiph  = iw  (iphas)
ir11ip = ir11(iphas)
ir22ip = ir22(iphas)
ir33ip = ir33(iphas)
ir12ip = ir12(iphas)
ir13ip = ir13(iphas)
ir23ip = ir23(iphas)
ieiph  = iep (iphas)

icliup = iclrtp(iuiph,icoef)
iclivp = iclrtp(iviph,icoef)
icliwp = iclrtp(iwiph,icoef)

ipcrom = ipproc(irom  (iphas))
ipbrom = ipprob(irom  (iphas))

if(iwarni(ieiph).ge.1) then
  if (iturb(iphas).eq.30) then
    write(nfecra,1000) iphas
  else
    write(nfecra,1001) iphas
  endif
endif

iitsla = 1
iitsli = 1

!     SI ITURB=30 (RIJ STD) ON STOCKE DIRECTEMENT LA PRODUCTION DANS
!     LE TABLEAU PRODUC
!     SI ITURB=31 (SSG) ON STOCKE LE GRADIENT DE VITESSE DANS GRDVIT

!===============================================================================
! 2.a CALCUL DU TENSEUR DE PRODUCTION POUR LE RIJ STANDARD
!     W7 = P11 , W8 = P22 , W9 = P33
!     W10 = P12 , W11 = P13 , W9 = P23
!===============================================================================

if (iturb(iphas).eq.30) then
! INITIALISATIONS DE W7 ... W12

  do ii = 1 , 6
    do iel = 1, ncel
      produc(ii,iel) = 0.0d0
    enddo
  enddo

! CALCUL DU GRADIENT DES 3 COMPOSANTES DE LA VITESSE

  iccocg = 1
  inc    = 1

! GRADIENT SUIVANT X

  nswrgp = nswrgr(iuiph)
  imligp = imligr(iuiph)
  iwarnp = iwarni(iuiph)
  epsrgp = epsrgr(iuiph)
  climgp = climgr(iuiph)
  extrap = extrag(iuiph)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iuiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iuiph)   , coefa(1,icliup) , coefb(1,icliup) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )


  do iel = 1 , ncel

    produc(1,iel) = produc(1,iel)                                 &
         - 2.0d0*(rtpa(iel,ir11ip)*w1(iel) +                      &
         rtpa(iel,ir12ip)*w2(iel) +                               &
         rtpa(iel,ir13ip)*w3(iel) )

    produc(4,iel) = produc(4,iel)                                 &
         - (rtpa(iel,ir12ip)*w1(iel) +                            &
         rtpa(iel,ir22ip)*w2(iel) +                               &
         rtpa(iel,ir23ip)*w3(iel) )

    produc(5,iel) = produc(5,iel)                                 &
         - (rtpa(iel,ir13ip)*w1(iel) +                            &
         rtpa(iel,ir23ip)*w2(iel) +                               &
         rtpa(iel,ir33ip)*w3(iel) )

  enddo

! Gradient suivant Y

  nswrgp = nswrgr(iviph)
  imligp = imligr(iviph)
  iwarnp = iwarni(iviph)
  epsrgp = epsrgr(iviph)
  climgp = climgr(iviph)
  extrap = extrag(iviph)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iviph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iviph)   , coefa(1,iclivp) , coefb(1,iclivp) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

  do iel = 1 , ncel

    produc(2,iel) = produc(2,iel)                                 &
         - 2.0d0*(rtpa(iel,ir12ip)*w1(iel) +                      &
         rtpa(iel,ir22ip)*w2(iel) +                               &
         rtpa(iel,ir23ip)*w3(iel) )

    produc(4,iel) = produc(4,iel)                                 &
         - (rtpa(iel,ir11ip)*w1(iel) +                            &
         rtpa(iel,ir12ip)*w2(iel) +                               &
         rtpa(iel,ir13ip)*w3(iel) )

    produc(6,iel) = produc(6,iel)                                 &
         - (rtpa(iel,ir13ip)*w1(iel) +                            &
         rtpa(iel,ir23ip)*w2(iel) +                               &
         rtpa(iel,ir33ip)*w3(iel) )

  enddo

! Gradient suivant Z

  nswrgp = nswrgr(iwiph)
  imligp = imligr(iwiph)
  iwarnp = iwarni(iwiph)
  epsrgp = epsrgr(iwiph)
  climgp = climgr(iwiph)
  extrap = extrag(iwiph)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iwiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iwiph)   , coefa(1,icliwp) , coefb(1,icliwp) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

  do iel = 1 , ncel

    produc(3,iel) = produc(3,iel)                                 &
         - 2.0d0*(rtpa(iel,ir13ip)*w1(iel) +                      &
         rtpa(iel,ir23ip)*w2(iel) +                               &
         rtpa(iel,ir33ip)*w3(iel) )

    produc(5,iel) = produc(5,iel)                                 &
         - (rtpa(iel,ir11ip)*w1(iel) +                            &
         rtpa(iel,ir12ip)*w2(iel) +                               &
         rtpa(iel,ir13ip)*w3(iel) )

    produc(6,iel) = produc(6,iel)                                 &
         - (rtpa(iel,ir12ip)*w1(iel) +                            &
         rtpa(iel,ir22ip)*w2(iel) +                               &
         rtpa(iel,ir23ip)*w3(iel) )

  enddo

else

!===============================================================================
! 2.b CALCUL DU GRADIENT DE VITESSE POUR LE RIJ SSG
!     GRDVIT(IEL,I,J) = dUi/dxj(IEL)
!===============================================================================

! CALCUL DU GRADIENT DES 3 COMPOSANTES DE LA VITESSE

  iccocg = 1
  inc    = 1

! GRADIENT SUIVANT X

  nswrgp = nswrgr(iuiph)
  imligp = imligr(iuiph)
  iwarnp = iwarni(iuiph)
  epsrgp = epsrgr(iuiph)
  climgp = climgr(iuiph)
  extrap = extrag(iuiph)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iuiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iuiph)   , coefa(1,icliup) , coefb(1,icliup) ,          &
   grdvit(1,1,1)   , grdvit(1,1,2)   , grdvit(1,1,3)   ,          &
!        -------------     -------------     -------------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )


! Gradient suivant Y

  nswrgp = nswrgr(iviph)
  imligp = imligr(iviph)
  iwarnp = iwarni(iviph)
  epsrgp = epsrgr(iviph)
  climgp = climgr(iviph)
  extrap = extrag(iviph)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iviph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iviph)   , coefa(1,iclivp) , coefb(1,iclivp) ,          &
   grdvit(1,2,1)   , grdvit(1,2,2)   , grdvit(1,2,3)   ,          &
!        -------------     -------------     -------------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )


! Gradient suivant Z

  nswrgp = nswrgr(iwiph)
  imligp = imligr(iwiph)
  iwarnp = iwarni(iwiph)
  epsrgp = epsrgr(iwiph)
  climgp = climgr(iwiph)
  extrap = extrag(iwiph)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iwiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iwiph)   , coefa(1,icliwp) , coefb(1,icliwp) ,          &
   grdvit(1,3,1)   , grdvit(1,3,2)   , grdvit(1,3,3)   ,          &
!        -------------     -------------     -------------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

endif


!===============================================================================
! 3.  CALCUL DU GRADIENT DE ROM POUR LES TERMES DE GRAVITE
!===============================================================================

if(igrari(iphas).eq.1) then

! Conditions aux limites : Dirichlet ROMB
!   On utilise VISCB pour stocker le coefb relatif a ROM
!   On impose en Dirichlet (COEFA) la valeur ROMB

  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

! Le choix ci dessous a l'avantage d'etre simple

  nswrgp = nswrgr(ir11ip)
  imligp = imligr(ir11ip)
  iwarnp = iwarni(ir11ip)
  epsrgp = epsrgr(ir11ip)
  climgp = climgr(ir11ip)
  extrap = extrag(ir11ip)
  iphydp = 0

  iivar = 0

!     Si on extrapole les termes sources et rho, on utilise cpdt rho^n
  ipcroo = ipcrom
  ipbroo = ipbrom
  if(isto2t(iphas).gt.0.and.iroext(iphas).gt.0) then
    ipcroo = ipproc(iroma(iphas))
    ipbroo = ipprob(iroma(iphas))
  endif

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iivar  , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   propce(1,ipcroo), propfb(1,ipbroo), viscb           ,          &
   grarox , graroy , graroz ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

endif


!===============================================================================
! 4.  Boucle sur les variables Rij (6 variables)
!     L'ordre est R11 R22 R33 R12 R13 R23 (La place de ces variables
!     est IR11.    ..
!     On resout les equation dans une routine semblable a covofi.F
!===============================================================================


do isou = 1, 6
  if    (isou.eq.1) then
    ivar   = ir11ip
  elseif(isou.eq.2) then
    ivar   = ir22ip
  elseif(isou.eq.3) then
    ivar   = ir33ip
  elseif(isou.eq.4) then
    ivar   = ir12ip
  elseif(isou.eq.5) then
    ivar   = ir13ip
  elseif(isou.eq.6) then
    ivar   = ir23ip
  endif
  ipp    = ipprtp(ivar)

  if (iilagr.eq.2 .and. iphas.eq.1) then
    iitsla = itsr11 + (isou-1)
    iitsli = itsli
  endif

  iicsmp = 1 + ncesmp*(ivar-1)
  iigamm = 1 + ncesmp*(ipriph-1)

  !     Rij-epsilon standard (LRR)
  if (iturb(iphas).eq.30) then
    call resrij                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , ivar   , isou   , ipp    ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm(iicsmp)  ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , produc , grarox , graroy , graroz ,          &
   ckupdc , smacel(iicsmp)  , smacel(iigamm)  ,                   &
   viscf  , viscb  , coefax ,                                     &
   tslagr(1,iitsla) , tslagr(1,iitsli) ,                          &
   dam    , xam    , drtp   , smbr   , rovsdt ,                   &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     ,                   &
   rdevel , rtuser , ra     )

  else
    !     Rij-epsilon SSG
    call resssg                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , ivar   , isou   , ipp    ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm(iicsmp)  ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , grdvit , grarox , graroy , graroz ,          &
   ckupdc , smacel(iicsmp)  , smacel(iigamm)  ,                   &
   viscf  , viscb  , coefax ,                                     &
   tslagr(1,iitsla) , tslagr(1,iitsli) ,                          &
   dam    , xam    , drtp   , smbr   , rovsdt ,                   &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     ,                   &
   rdevel , rtuser , ra     )
  endif

enddo

!===============================================================================
! 5.  RESOLUTION DE EPSILON
!===============================================================================

ivar   = ieiph
ipp    = ipprtp(ivar)
isou   = 7

iicsmp = 1 + ncesmp*(ivar-1)
iigamm = 1 + ncesmp*(ipriph-1)

call reseps                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , ivar   , isou   , ipp    ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm(iicsmp)  ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , grdvit , produc ,grarox , graroy , graroz ,  &
   ckupdc , smacel(iicsmp)  , smacel(iigamm) ,                    &
   viscf  , viscb  ,                                              &
   tslagr ,                                                       &
   dam    , xam    , drtp   , smbr   , rovsdt ,                   &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     ,                   &
   rdevel , rtuser , ra     )

!===============================================================================
! 6. CLIPPING
!===============================================================================

iclip  = 2
call clprij                                                       &
!==========
 ( ncelet , ncel   , nvar   , nphas  ,                            &
   iphas  , iclip  ,                                              &
   propce , rtpa   , rtp    )


!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** PHASE ',I4,' RESOLUTION DU Rij-EPSILON LRR             ',/,&
'      ------------------------------------------             ',/)
 1001 format(/,                                                   &
'   ** PHASE ',I4,' RESOLUTION DU Rij-EPSILON SSG             ',/,&
'      ------------------------------------------             ',/)

#else

 1000 format(/,                                                   &
'   ** PHASE ',I4,' SOLVING Rij-EPSILON LRR'                   ,/,&
'      ------------------------------------'                   ,/)
 1001 format(/,                                                   &
'   ** PHASE ',I4,' SOLVING Rij-EPSILON SSG'                   ,/,&
'      ------------------------------------'                   ,/)

#endif

!----
! FIN
!----

return

end subroutine
