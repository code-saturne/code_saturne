!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine phyvar &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , w10    , w11    , w12    ,          &
   xmij   ,                                                       &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! REMPLISSAGE DES GRANDEURS PHYSIQUES VARIABLES EN TEMPS
!    ESSENTIELLEMENT LA VISCOSITE TURBULENTE.
!    ON APPELLE UN SOUS PROGRAMME UTILISATEUR QUI PERMET DE
!    SPECIFIER ROM, ROMB, VISCL, VISCLS ...

! Arguments
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
! w1..12(ncelet    ! tr ! --- ! tableau de travail                             !
! xmij(ncelet,6    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "pointe.h"
include "albase.h"
include "period.h"
include "parall.h"
include "ihmpre.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "matiss.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

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
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision w9(ncelet),w10(ncelet),w11(ncelet),w12(ncelet)
double precision xmij(ncelet,6)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          ivar  , iel   , ifac  , iscal , iphas , iph
integer          ii    , iok   , iok1  , iok2  , iisct
integer          nphmx , nn
integer          ibrom(nphsmx) , ipcrom, ipbrom, ipcvst
integer          ikiph , ieiph , ir11ip, ir22ip, ir33ip
integer          ipccp , ipcvis, iphiph, ipcvma
integer          iclipc
double precision xk, xe, xnu, xrom, vismax(nscamx), vismin(nscamx)
double precision varmn(4), varmx(4), tt, ttmin, ttke, vistot

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

ipass = ipass + 1

!===============================================================================
! 2.  PREPARATION DE LA PERIODICITE DE ROTATION
!       CALCUL DE DUDXYZ ET DRDXYZ (gradients sur les halos avec prise
!       en compte des periodicites pour exploitation dans pering, inimas)
!===============================================================================

if(iperot.gt.0) then

  do iphas = 1,nphas

    iph = iphas
    call perinu                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , iph    ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser ,                                              &
   ra(idudxy) ,                                                   &
   ra     )

  enddo


  do iphas = 1, nphas

    if(itytur(iphas).eq.3) then

      iph = iphas
      call perinr                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , iph    ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser ,                                              &
   ra(idrdxy) ,                                                   &
   ra     )

    endif

  enddo

endif

!===============================================================================
! 4.  ON REND LA MAIN A L'UTILISATEUR POUR LA PROGRAMMATION DES
!      GRANDEURS PHYSIQUES VARIABLES QUI LUI SONT PROPRES
!===============================================================================

nphmx = nphsmx
do iphas = 1, nphsmx
  ibrom(iphas) = 0
enddo


if (ippmod(iphpar).ge.1) then
  call ppphyv                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ibrom  ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

endif


if (imatis.eq.1) then

!     Matisse
  call mtphyv                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ibrom  ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   ,                                     &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

else

  ! - Interface Code_Saturne
  !   ======================

  if (iihmpr.eq.1) then
    call uiphyv                                                    &
    !===========
  ( ncel, ncelet, nscaus,                                         &
    irom, iviscl, icp,    ivisls, irovar, ivivar,                 &
    isca, iscalt, iscavr, ipproc,                                 &
    ro0,  cp0,    rtp,    propce)
  endif

!     Utilisateur standard
  call usphyv                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ibrom  ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   ,                                     &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

endif

!  ROMB SUR LES BORDS : VALEUR PAR DEFAUT (CELLE DE LA CELLULE VOISINE)

do iphas = 1, nphas
  if (ibrom(iphas).eq.0) then
    ipcrom = ipproc(irom(iphas))
    ipbrom = ipprob(irom(iphas))
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      propfb(ifac,ipbrom) = propce(iel ,ipcrom)
    enddo
  endif
enddo

!  Au premier pas de temps du calcul
!     Si on a indique que rho (visc) etait constant
!       et qu'on l'a modifie dans usphyv, ca ne va pas
!     On se sert de irovar (ivivar) pour ecrire et lire
!       rho (visc) dans le fichier suite

if(ntcabs.eq.ntpabs+1) then

!     Masse volumique aux cellules et aux faces de bord
  iok1 = 0
  do iphas = 1, nphas
    if(irovar(iphas).eq.0) then
      ipcrom = ipproc(irom(iphas))
      ipbrom = ipprob(irom(iphas))
      do iel = 1, ncel
        if( abs(propce(iel ,ipcrom)-ro0   (iphas)).gt.epzero) then
          iok1 = 1
        endif
      enddo
      do ifac = 1, nfabor
        if( abs(propfb(ifac,ipbrom)-ro0   (iphas)).gt.epzero) then
          iok1 = 1
        endif
      enddo
    endif
  enddo
  if(iok1.ne.0) then
    write(nfecra,9001)
  endif

!     Viscosite moleculaire aux cellules
  iok2 = 0
  do iphas = 1, nphas
    if(ivivar(iphas).eq.0) then
      ipcvis = ipproc(iviscl(iphas))
      do iel = 1, ncel
        if( abs(propce(iel ,ipcvis)-viscl0(iphas)).gt.epzero) then
          iok2 = 1
        endif
      enddo
    endif
  enddo
  if(iok2.ne.0) then
    if ( ippmod(icompf) .ge. 0 ) then
      write(nfecra,9003)
    else
      write(nfecra,9002)
    endif
  endif

  if(iok1.ne.0.or.iok2.ne.0) then
    call csexit(1)
  endif

endif

!===============================================================================
! 3.  CALCUL DE LA VISCOSITE TURBULENTE
!===============================================================================

! --- Boucle sur les phases : Debut
do iphas = 1, nphas

  if     (iturb(iphas).eq. 0) then

! 3.1 LAMINAIRE
! ==============

    ipcvst = ipproc(ivisct(iphas))

    do iel = 1, ncel
      propce(iel,ipcvst) = 0.d0
    enddo

  elseif (iturb(iphas).eq.10) then

! 3.2 LONGUEUR DE MELANGE
! ========================

    iph = iphas
    call vislmg                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdc(iphas) , ncetsm(iphas) ,                                &
   nideve , nrdeve , nituse , nrtuse , iph    ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) , ia(iicesm(iphas)) , ia(iitpsm(iphas)) ,    &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ra(ickupd(iphas)), ra(ismace(iphas)),        &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

  elseif (itytur(iphas).eq.2) then

! 3.3 K-EPSILON
! ==============

    ikiph  = ik(iphas)
    ieiph  = iep(iphas)
    ipcvst = ipproc(ivisct(iphas))
    ipcrom = ipproc(irom  (iphas))

    do iel = 1, ncel
      xk = rtp(iel,ikiph)
      xe = rtp(iel,ieiph)
      propce(iel,ipcvst) = propce(iel,ipcrom)*cmu*xk**2/xe
    enddo

  elseif (itytur(iphas).eq.3) then

! 3.4 Rij-EPSILON
! ================

    ir11ip = ir11(iphas)
    ir22ip = ir22(iphas)
    ir33ip = ir33(iphas)
    ieiph  = iep(iphas)
    ipcvst = ipproc(ivisct(iphas))
    ipcrom = ipproc(irom  (iphas))

    do iel = 1, ncel
      xk = 0.5d0*(rtp(iel,ir11ip)+rtp(iel,ir22ip)+rtp(iel,ir33ip))
      xe = rtp(iel,ieiph)
      propce(iel,ipcvst) = propce(iel,ipcrom)*cmu*xk**2/xe
    enddo

  elseif (iturb(iphas).eq.40) then

! 3.5 LES Smagorinsky
! ===================

    iph = iphas

    call vissma                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdc(iphas) , ncetsm(iphas) ,                                &
   nideve , nrdeve , nituse , nrtuse , iph    ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) ,                                            &
   ia(iicesm(iphas)) , ia(iitpsm(iphas)),                         &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ra(ickupd(iphas)), ra(ismace(iphas)),        &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

  elseif(iturb(iphas).eq.41) then

! 3.6 LES dynamique
! =================

    iph = iphas

    call visdyn                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdc(iphas) , ncetsm(iphas) ,                                &
   nideve , nrdeve , nituse , nrtuse , iph    ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) ,                                            &
   ia(iicesm(iphas)) , ia(iitpsm(iphas)),                         &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ra(ickupd(iphas)), ra(ismace(iphas)),        &
   propce(1,ipproc(ismago(iphas))) ,                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     , w10    , xmij   , &
   rdevel , rtuser , ra     )

  elseif (iturb(iphas).eq.42) then

! 3.7 LES WALE
! ============

    iph = iphas

    call viswal                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdc(iphas) , ncetsm(iphas) ,                                &
   nideve , nrdeve , nituse , nrtuse , iph    ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) ,                                            &
   ia(iicesm(iphas)) , ia(iitpsm(iphas)),                         &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ra(ickupd(iphas)), ra(ismace(iphas)),        &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   w9     , w10    , w11    , w12    ,                            &
   rdevel , rtuser , ra     )

  elseif (iturb(iphas).eq.50) then

! 3.8 v2f phi-model
! =================
    ikiph  = ik(iphas)
    ieiph  = iep(iphas)
    iphiph = iphi(iphas)
    ipcvis = ipproc(iviscl(iphas))
    ipcvst = ipproc(ivisct(iphas))
    ipcrom = ipproc(irom  (iphas))

    do iel = 1, ncel
      xk = rtp(iel,ikiph)
      xe = rtp(iel,ieiph)
      xrom = propce(iel,ipcrom)
      xnu = propce(iel,ipcvis)/xrom
      ttke = xk / xe
      ttmin = cv2fct*sqrt(xnu/xe)
      tt = max(ttke,ttmin)
      propce(iel,ipcvst) = cv2fmu*xrom*tt*rtp(iel,iphiph)         &
           *rtp(iel,ikiph)
    enddo

  elseif (iturb(iphas).eq.60) then

! 3.9 K-OMEGA SST
! ===============

    iph = iphas
    call vissst                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdc(iphas) , ncetsm(iphas) ,                                &
   nideve , nrdeve , nituse , nrtuse , iph    ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) , ia(iicesm(iphas)) , ia(iitpsm(iphas)) ,    &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ra(ickupd(iphas)), ra(ismace(iphas)),        &
   ra(is2kw(iphas)), ra(idvukw(iphas)),                           &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

  endif

enddo
! --- Boucle sur les phases : fin

!===============================================================================
! 4.  MODIFICATION UTILISATEUR DE LA VISCOSITE TURBULENTE
!===============================================================================

do iphas = 1, nphas

  iph = iphas
  call usvist                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdc(iphas)   , ncetsm(iphas)   ,                            &
   nideve , nrdeve , nituse , nrtuse , iph    ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) , ia(iicesm(iphas)) , ia(iitpsm(iphas)) ,    &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ra(ickupd(iphas)), ra(ismace(iphas)),        &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

enddo

!===============================================================================
! 5.  CLIPPING DE LA VISCOSITE TURBULENTE EN LES DYNAMIQUE
!===============================================================================

do iphas = 1, nphas
! Pour la LES en modele dynamique on clippe la viscosite turbulente de maniere
! a ce que mu+mu_t soit positif, .e. on autorise mu_t legerement negatif
! La diffusivite turbulente des scalaires (mu_t/sigma), elle, sera clippee a 0
! dans covofi

  if (iturb(iphas).eq.41) then
    ipcvis = ipproc(iviscl(iphas))
    ipcvst = ipproc(ivisct(iphas))
    iclipc = 0
    do iel = 1, ncel
      vistot = propce(iel,ipcvis) + propce(iel,ipcvst)
      if(vistot.lt.0.d0) then
        propce(iel,ipcvst) = 0.d0
        iclipc = iclipc + 1
      endif
    enddo
    if(iwarni(iu(iphas)).ge.1) then
      if(irangp.ge.0) then
        call parcpt(iclipc)
        !==========
      endif
      write(nfecra,1000) iclipc
    endif
  endif

enddo

!===============================================================================
! 6.  MODIFICATION UTILISATEUR DE LA VISCOSITE DE MAILLAGE EN ALE
!===============================================================================

if (iale.eq.1.and.ntcabs.eq.0) then

  ! - Interface Code_Saturne
  !   ======================

  if (iihmpr.eq.1) then

    call uivima                         &
    !==========
  ( ncel,                             &
    propce(1,ipproc(ivisma(1))),      &
    propce(1,ipproc(ivisma(2))),      &
    propce(1,ipproc(ivisma(3))),      &
    xyzcen, dtref, ttcabs, ntcabs )

  endif

  call usvima                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  ,                                              &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , propce(1,ipproc(ivisma(1))) ,                &
   propce(1,ipproc(ivisma(2))) , propce(1,ipproc(ivisma(3))) ,    &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

endif

!===============================================================================
! 7.  IMPRESSIONS DE CONTROLE DES VALEURS ENTREES PAR L'UTILISATEUR
!===============================================================================


! ---> Calcul des bornes des variables par phase et impressions

! Indicateur d'erreur
iok = 0

do iphas = 1, nphas

! Rang des variables dans PROPCE
  ipcrom = ipproc(irom(iphas))
  ipcvis = ipproc(iviscl(iphas))
  ipcvst = ipproc(ivisct(iphas))
  if(icp(iphas).gt.0) then
    ipccp  = ipproc(icp   (iphas))
    nn     = 4
  else
    ipccp = 0
    nn    = 3
  endif

! Rang des variables dans PROPFB
  ipbrom = ipprob(irom(iphas))

! Min et max sur les cellules
  do ii = 1, nn
    ivar = 0
    if(ii.eq.1) ivar = ipcrom
    if(ii.eq.2) ivar = ipcvis
    if(ii.eq.3) ivar = ipcvst
    if(ii.eq.4) ivar = ipccp
    if(ivar.gt.0) then
      varmx(ii) = propce(1,ivar)
      varmn(ii) = propce(1,ivar)
      do iel = 2, ncel
        varmx(ii) = max(varmx(ii),propce(iel,ivar))
        varmn(ii) = min(varmn(ii),propce(iel,ivar))
      enddo
      if (irangp.ge.0) then
        call parmax (varmx(ii))
        !==========
        call parmin (varmn(ii))
        !==========
      endif
    endif
  enddo

! Min et max sur les faces de bord (masse volumique uniquement)
  ii   = 1
  ivar = ipbrom
  do ifac = 1, nfabor
    varmx(ii) = max(varmx(ii),propfb(ifac,ivar))
    varmn(ii) = min(varmn(ii),propfb(ifac,ivar))
  enddo
  if (irangp.ge.0) then
    call parmax (varmx(ii))
    !==========
    call parmin (varmn(ii))
    !==========
  endif

! Impressions
  iok1 = 0
  do ii = 1, nn
    if(ii.eq.1) chaine = nomvar(ipppro(ipproc(irom  (iphas))))
    if(ii.eq.2) chaine = nomvar(ipppro(ipproc(iviscl(iphas))))
    if(ii.eq.3) chaine = nomvar(ipppro(ipproc(ivisct(iphas))))
    if(ii.eq.4) chaine = nomvar(ipppro(ipproc(icp   (iphas))))
    if(iwarni(iu(iphas)).ge.1.or.ipass.eq.1.or.                   &
                                          varmn(ii).lt.0.d0)then
      if(iok1.eq.0) then
        write(nfecra,3010)iphas
        iok1 = 1
      endif
      if ((ii.ne.3).or.(iturb(iphas).ne.0))                       &
           write(nfecra,3011)chaine(1:8),varmn(ii),varmx(ii)
    endif
  enddo
  if(iok1.eq.1) write(nfecra,3012)

! Verifications de valeur physique

! Masse volumique definie
  ii = 1
  chaine = nomvar(ipppro(ipproc(irom  (iphas))))
  if (varmn(ii).lt.0.d0) then
    write(nfecra,9011)chaine(1:8),varmn(ii)
    iok = iok + 1
  endif

! Viscosite moleculaire definie
  ii = 2
  chaine = nomvar(ipppro(ipproc(iviscl(iphas))))
  if (varmn(ii).lt.0.d0) then
    write(nfecra,9011)chaine(1:8),varmn(ii)
    iok = iok + 1
  endif

! Viscosite turbulente definie
! on ne clippe pas mu_t en modele LES dynamique, car on a fait
! un clipping sur la viscosite totale
  ii = 3
  chaine = nomvar(ipppro(ipproc(ivisct(iphas))))
  if (varmn(ii).lt.0.d0.and.iturb(iphas).ne.41) then
    write(nfecra,9012)iphas,varmn(ii)
    iok = iok + 1
  endif

! Chaleur specifique definie
  if(icp(iphas).gt.0) then
    ii = 4
    chaine = nomvar(ipppro(ipproc(icp   (iphas))))
    if (varmn(ii).lt.0.d0) then
      iisct = 0
      do iscal = 1, nscal
        if (iscsth(iscal).ne.0) then
          iisct = 1
        endif
      enddo
      if (iisct.eq.1) then
        write(nfecra,9011)chaine(1:8),varmn(ii)
        iok = iok + 1
      endif
    endif
  endif

enddo



! ---> Calcul des bornes des scalaires et impressions

if(nscal.ge.1) then

  iok1 = 0
  do iscal = 1, nscal

    if(ivisls(iscal).gt.0) then
      ipcvis = ipproc(ivisls(iscal))
    else
      ipcvis = 0
    endif

    vismax(iscal) = -grand
    vismin(iscal) =  grand
    if(ipcvis.gt.0) then
      do iel = 1, ncel
        vismax(iscal) = max(vismax(iscal),propce(iel,ipcvis))
        vismin(iscal) = min(vismin(iscal),propce(iel,ipcvis))
      enddo
      if (irangp.ge.0) then
        call parmax (vismax(iscal))
        !==========
        call parmin (vismin(iscal))
        !==========
      endif
    else
      vismax(iscal) = visls0(iscal)
      vismin(iscal) = visls0(iscal)
    endif

    ivar = isca(iscal)
    if(iwarni(ivar).ge.1.or.ipass.eq.1.or.                        &
                                        vismin(iscal).le.0.d0)then
      chaine = nomvar(ipprtp(ivar))
      if(iok1.eq.0) then
        write(nfecra,3110)
        iok1 = 1
      endif
      write(nfecra,3111)chaine(1:8),iscal,                        &
                        vismin(iscal),vismax(iscal)
    endif

  enddo
  if(iok1.eq.1) write(nfecra,3112)

! Verifications de valeur physique

! IOK a deja ete initialise pour les valeurs liees aux phases

  do iscal = 1, nscal

    ivar = isca(iscal)

    if (vismin(iscal).lt.0.d0) then
      chaine = nomvar(ipprtp(ivar))
      write(nfecra,9111)chaine(1:8),iscal,vismin(iscal)
      iok = iok + 1
    endif

  enddo

endif

! ---> Calcul des bornes de viscosite de maillage en ALE

if (iale.eq.1.and.ntcabs.eq.0) then

  iok1 = 0
  nn = 1
  if (iortvm.eq.1) nn = 3
  do ii = 1, nn
    ipcvma = ipproc(ivisma(ii))

! Min et max sur les cellules
    varmx(1) = propce(1,ipcvma)
    varmn(1) = propce(1,ipcvma)
    do iel = 2, ncel
      varmx(1) = max(varmx(1),propce(iel,ipcvma))
      varmn(1) = min(varmn(1),propce(iel,ipcvma))
    enddo
    if (irangp.ge.0) then
      call parmax (varmx(1))
      !==========
      call parmin (varmn(1))
      !==========
    endif

! Impressions
    chaine = nomvar(ipppro(ipcvma))
    if(iwarni(iuma).ge.1.or.ipass.eq.1.or.                        &
         varmn(1).lt.0.d0)then
      if (iok1.eq.0) then
        write(nfecra,3210)
        iok1 = 1
      endif
      write(nfecra,3211)chaine(1:8),varmn(1),varmx(1)
    endif

! Verifications de valeur physique

! Viscosite de maillage definie
    chaine = nomvar(ipppro(ipcvma))
    if (varmn(1).le.0.d0) then
      write(nfecra,9211) varmn(1)
      iok = iok + 1
    endif

  enddo

  if (iok1.eq.1) write(nfecra,3212)

endif

! --->  arret eventuel

if(iok.ne.0) then
  write(nfecra,9999)iok
  call csexit (1)
endif

!===============================================================================
! 8.  ECHANGES
!===============================================================================

! Pour navsto et vissec on a besoin de ROM dans le halo


do iphas = 1, nphas

  ipcrom = ipproc(irom(iphas))

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(propce(1,ipcrom))
    !==========
  endif

enddo

!----
! FORMATS
!----

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
' Nb de clippings de la viscosite totale (mu+mu_t>0) :',I10,/)
 3010 format(                                                           &
' --- Phase : ',I10                                            ,/,&
' ---------------------------------                           ',/,&
' Propriete  Valeur min  Valeur max                           ',/,&
' ---------------------------------                           '  )
 3011 format(                                                           &
 2x,    a8,      e12.4,      e12.4                               )
 3012 format(                                                           &
' ---------------------------------                           ',/)
 3110 format(                                                           &
' --- Diffusivite :                                           ',/,&
' ---------------------------------------                     ',/,&
' Scalaire Numero  Valeur min  Valeur max                     ',/,&
' ---------------------------------------                     '  )
 3111 format(                                                           &
 1x,    a8,   i7,      e12.4,      e12.4                         )
 3112 format(                                                           &
' ---------------------------------------                     ',/)
 3210 format(                                                           &
' --- Viscosite de maillage (methode ALE)                     ',/,&
' ---------------------------------                           ',/,&
' Propriete  Valeur min  Valeur max                           ',/,&
' ---------------------------------                           '  )
 3211 format(                                                           &
 2x,    a8,      e12.4,      e12.4                               )
 3212 format(                                                           &
' ---------------------------------                           ',/)

 9001  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    INCOHERENCE ENTRE USPHYV ET USINI1 POUR                 ',/,&
'@                                          LA MASSE VOLUMIQUE',/,&
'@                                                            ',/,&
'@  On a indique dans usini1 que la masse volumique etait     ',/,&
'@     constante (IROVAR(IPHAS)=0) mais on a modifie ses      ',/,&
'@     valeurs aux cellules ou aux faces de bord dans usphyv. ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1 et usphyv.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9002  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    INCOHERENCE ENTRE USPHYV ET USINI1 POUR                 ',/,&
'@                                    LA VISCOSITE MOLECULAIRE',/,&
'@                                                            ',/,&
'@  On a indique dans usini1 que la viscosite moleculaire     ',/,&
'@     etait constante (IVIVAR(IPHAS)=0) mais on a modifie ses',/,&
'@     valeurs aux cellules dans usphyv.                      ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1 et usphyv.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9003  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    MODULE COMPRESSIBLE                                     ',/,&
'@    INCOHERENCE ENTRE USCFPV ET USCFX1 POUR                 ',/,&
'@                                    LA VISCOSITE MOLECULAIRE',/,&
'@                                                            ',/,&
'@  En compressible la viscosite moleculaire est constante par',/,&
'@     defaut (IVIVAR(IPHAS)=0) et la valeur de IVIVAR n''a   ',/,&
'@     pas ete modifiee dans uscfx1. Pourtant, on a modifie   ',/,&
'@     les valeurs de la viscosite moleculaire dans uscfpv.   ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfx1 et uscfpv.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9011  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    LA PROPRIETE PHYSIQUE ',A8  ,' N A PAS ETE              ',/,&
'@                                       CORRECTEMENT DEFINIE.',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  La propriete physique identifiee ci-dessus est variable et',/,&
'@    le minimum atteint est ',E12.4                           ,/,&
'@  Verifier que cette propriete a ete definie dans usphyv et ',/,&
'@    que la loi adoptee conduit a des valeurs correctes.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9012  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    LA VISCOSITE TURBULENTE DE LA PHASE ',I10                ,/,&
'@                           N A PAS ETE CORRECTEMENT DEFINIE.',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Le minimum atteint est ',E12.4                             ,/,&
'@  Verifier le cas echeant la definition de la masse         ',/,&
'@    volumique dans usphyv et la modification de la viscosite',/,&
'@    turbulente dans usvist.                                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9111  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    LA DIFFUSIVITE DU SCALAIRE ',A8                          ,/,&
'@       (SCALAIRE NUMERO ',I10   ,') N A PAS ETE             ',/,&
'@                                       CORRECTEMENT DEFINIE.',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  La propriete physique identifiee ci-dessus est variable et',/,&
'@    le minimum atteint est ',E12.4                           ,/,&
'@  Verifier que cette propriete a ete definie dans usphyv et ',/,&
'@    que la loi adoptee conduit a des valeurs correctes.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9211  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    LA VISCOSITE DE MAILLAGE N A PAS ETE                    ',/,&
'@                                       CORRECTEMENT DEFINIE.',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Le minimum atteint est ',E12.4                             ,/,&
'@  Verifier le cas echeant la modification de la viscosite   ',/,&
'@    dans usvima ou dans l interface graphique.              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DES GRANDEURS PHYSIQUES ONT DES VALEURS INCORRECTES     ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier usphyv.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                           &
' Nb of clippings for the effective viscosity (mu+mu_t>0):',I10,/)
 3010 format(                                                           &
' --- Phase: ',I10                                             ,/,&
' ---------------------------------                           ',/,&
' Property   Min. value  Max. value                           ',/,&
' ---------------------------------                           '  )
 3011 format(                                                           &
 2x,    a8,      e12.4,      e12.4                               )
 3012 format(                                                           &
' ---------------------------------                           ',/)
 3110 format(                                                           &
' --- Diffusivity:                                            ',/,&
' ---------------------------------------                     ',/,&
' Scalar   Number  Min. value  Max. value                     ',/,&
' ---------------------------------------                     '  )
 3111 format(                                                           &
 1x,    a8,   i7,      e12.4,      e12.4                         )
 3112 format(                                                           &
' ---------------------------------------                     ',/)
 3210 format(                                                           &
' --- Mesh viscosity (ALE method)                             ',/,&
' ---------------------------------                           ',/,&
' Property   Min. value  Max. value                           ',/,&
' ---------------------------------                           '  )
 3211 format(                                                           &
 2x,    a8,      e12.4,      e12.4                               )
 3212 format(                                                           &
' ---------------------------------                           ',/)

 9001  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION   ',/,&
'@    ========                                                ',/,&
'@    INCOHERENCY BETWEEN USPHYV AND USINI1 FOR THE DENSITY   ',/,&
'@                                                            ',/,&
'@  The density has been declared constant in usini1          ',/,&
'@     (IROVAR(IPHAS)=0) but its value has been modified      ',/,&
'@     in cells of at boundary faces in usphyv.               ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usini1 and usphyv.                                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9002  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION   ',/,&
'@    ========                                                ',/,&
'@    INCOHERENCY BETWEEN USPHYV AND USINI1 FOR               ',/,&
'@                                     THE MOLECULAR VISCOSITY',/,&
'@                                                            ',/,&
'@  The molecular viscosity has been declared constant in     ',/,&
'@     usini1 (IVIVAR(IPHAS)=0) but its value has been        ',/,&
'@     modified in cells or at boundary faces in usphyv.      ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usini1 and usphyv.                                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9003  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION   ',/,&
'@    ========                                                ',/,&
'@    INCOHERENCY BETWEEN USCFPV AND USCFX1 FOR               ',/,&
'@                                     THE MOLECULAR VISCOSITY',/,&
'@                                                            ',/,&
'@  In the compressible module, the molecular viscosity is    ',/,&
'@     constant by default (IVIVAR(IPHAS)=0) and the value    ',/,&
'@     of IVIVAR  has not been modified in uscfx1. Yet, its   ',/,&
'@     value has been modified in uscfpv.                     ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify uscfx1 and uscfpv.                                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9011  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION   ',/,&
'@    ========                                                ',/,&
'@    THE PHYSICAL PROPERTY ',A8  ,' HAS NOT BEEN             ',/,&
'@                                          CORRECTLY DEFINED.',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  The physical property identified is variable and the      ',/,&
'@    minimum reached is ',E12.4                               ,/,&
'@  Verify that this property has been defined in usphyv and  ',/,&
'@    that the chosen law leads to correct values.            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9012  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION   ',/,&
'@    ========                                                ',/,&
'@    THE TURBULENT VISCOSITY OF PHASE ',I10 ,' HAS NOT BEEN  ',/,&
'@                                          CORRECTLY DEFINED.',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  The  minimum reached is ',E12.4                            ,/,&
'@  Verify the density definition in usphyv and a turbulent   ',/,&
'@    viscosity modification in usvist (if any).              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9111  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION   ',/,&
'@    ========                                                ',/,&
'@    THE DIFFUSIVITY OF THE SCALAR ',A8                       ,/,&
'@       (SCALAR NUMBER ',I10   ,') HAS NOT BEEN              ',/,&
'@                                          CORRECTLY DEFINED.',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  The physical property identified is variable and the      ',/,&
'@    minimum reached is ',E12.4                               ,/,&
'@  Verify that this property has been defined in usphyv and  ',/,&
'@    that the chosen law leads to correct values.            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9211  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION   ',/,&
'@    ========                                                ',/,&
'@    THE MESH VISCOSITY HAS NOT BEEN CORRECTLY DEFINED.      ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  The  minimum reached is ',E12.4                            ,/,&
'@  Verify that this property has been defined in usvima and  ',/,&
'@    that the chosen law leads to correct values.            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION   ',/,&
'@    ========                                                ',/,&
'@    SOME PHYSICAL QUANTITIES HAVE INCORRECT VALUES          ',/,&
'@                                                            ',/,&
'@  The calculation will not be run (',I10,' errors).         ',/,&
'@                                                            ',/,&
'@  Refer to previous warnings for further information.       ',/,&
'@  Verify usphyv.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! FIN
!----

return
end subroutine
