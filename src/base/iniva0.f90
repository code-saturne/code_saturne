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

subroutine iniva0 &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncofab ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , propce , propfa , propfb ,                   &
   coefa  , coefb  , frcxt  ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL, DU PAS DE TEMPS
!  ET DU TABLEAU INDICATEUR DU CALCUL DE LA DISTANCE A LA PAROI
!  AUX VALEURS PAR DEFAUT
!                AVANT LECTURE EVENTUELLE DU FICHIER SUITE ET
!                AVANT DE PASSER LA MAIN A L'UTILISATEUR
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
! ncofab           ! e  ! <-- ! nombre de couples coefa/b pour les cl          !
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
!                  !    !     !  a la face interne                             !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa coefb      ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! frcxt(ncelet,    ! tr ! <-- ! force exterieure generant la pression          !
!   3,nphas)       !    !     !  hydrostatique                                 !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use numvar
use optcal
use cstphy
use cstnum
use pointe
use entsor
use albase
use parall
use period
use ppppar
use ppthch
use ppincl
use cplsat

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , ncofab
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,ncofab), coefb(nfabor,ncofab)
double precision frcxt(ncelet,3,nphas)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          iis   , ivar  , iphas , iphass, iscal , imom
integer          iel   , ifac
integer          iclip , ii    , jj    , idim
integer          iiflum, iiflua
integer          iirom , iiromb, iiroma
integer          iivisl, iivist, iivisa, iivism
integer          iicp  , iicpa
integer          iiviss, iiptot
integer          iptsna, iptsta, iptsca
integer          ikiph , ieiph , iphiph, ifbiph, iomgip
integer          ir11ip, ir22ip, ir33ip, ir12ip, ir13ip, ir23ip
integer          nn
double precision ro0iph, visiph
double precision xxk, xcmu, trii

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Initialize variables to avoid compiler warnings

jj = 0

! Memoire

idebia = idbia0
idebra = idbra0

! En compressible, ISYMPA initialise (= 1) car utile dans le calcul
!     du pas de temps variable avant passage dans les C.L.

if ( ippmod(icompf).ge.0 ) then
  do ifac = 1, nfabor*nphas
    ia(iisymp+ifac-1) = 1
  enddo
endif

!===============================================================================
! 2. PAS DE TEMPS
!===============================================================================

do iel = 1, ncel
  dt (iel) = dtref
enddo

!===============================================================================
! 3.  INITIALISATION DES PROPRIETES PHYSIQUES
!===============================================================================

do iphas = 1, nphas

!     Masse volumique
  iirom  = ipproc(irom  (iphas))
  iiromb = ipprob(irom  (iphas))
  ro0iph = ro0   (iphas)

!     Masse volumique aux cellules (et au pdt precedent si ordre2 ou icalhy)
  do iel = 1, ncel
    propce(iel,iirom)  = ro0iph
  enddo
  if(iroext(iphas).gt.0.or.icalhy.eq.1) then
    iiroma = ipproc(iroma (iphas))
    do iel = 1, ncel
      propce(iel,iiroma) = propce(iel,iirom)
    enddo
  endif
!     Masse volumique aux faces de bord (et au pdt precedent si ordre2)
  do ifac = 1, nfabor
    propfb(ifac,iiromb) = ro0iph
  enddo
  if(iroext(iphas).gt.0) then
    iiroma = ipprob(iroma (iphas))
    do ifac = 1, nfabor
      propfb(ifac,iiroma) = propfb(ifac,iiromb)
    enddo
  endif

!     Viscosite moleculaire
  iivisl = ipproc(iviscl(iphas))
  iivist = ipproc(ivisct(iphas))
  visiph = viscl0(iphas)

!     Viscosite moleculaire aux cellules (et au pdt precedent si ordre2)
  do iel = 1, ncel
    propce(iel,iivisl) = visiph
  enddo
  if(iviext(iphas).gt.0) then
    iivisa = ipproc(ivisla(iphas))
    do iel = 1, ncel
      propce(iel,iivisa) = propce(iel,iivisl)
    enddo
  endif
!     Viscosite turbulente aux cellules (et au pdt precedent si ordre2)
  do iel = 1, ncel
    propce(iel,iivist) = 0.d0
  enddo
  if(iviext(iphas).gt.0) then
    iivisa = ipproc(ivista(iphas))
    do iel = 1, ncel
      propce(iel,iivisa) = propce(iel,iivist)
    enddo
  endif

!     Chaleur massique aux cellules (et au pdt precedent si ordre2)
  if(icp(iphas).gt.0) then
    iicp = ipproc(icp(iphas))
    do iel = 1, ncel
      propce(iel,iicp) = cp0(iphas)
    enddo
    if(icpext(iphas).gt.0) then
      iicpa  = ipproc(icpa(iphas))
      do iel = 1, ncel
        propce(iel,iicpa ) = propce(iel,iicp)
      enddo
    endif
  endif

! La pression totale sera initialisee a P0 + rho.g.r dans INIVAR
!  si l'utilisateur n'a pas fait d'initialisation personnelle
! Non valable en compressible
  if (ippmod(icompf).lt.0) then
    iiptot = ipproc(iprtot(iphas))
    do iel = 1, ncel
      propce(iel,iiptot) = - rinfin
    enddo
  endif

enddo


!     Diffusivite des scalaires
do iscal = 1, nscal
  if(ivisls(iscal).gt.0) then
    iiviss = ipproc(ivisls(iscal))
!     Diffusivite aux cellules (et au pdt precedent si ordre2)
    do iel = 1, ncel
      propce(iel,iiviss) = visls0(iscal)
    enddo
    if(ivsext(iphsca(iscal)).gt.0) then
      iivisa = ipproc(ivissa(iscal))
      do iel = 1, ncel
        propce(iel,iivisa) = propce(iel,iiviss)
      enddo
    endif
  endif
enddo


!     Viscosite de maillage en ALE
if (iale.eq.1) then
  nn = 1
  if (iortvm.eq.1) nn = 3
  do ii = 1, nn
    iivism = ipproc(ivisma(ii))
    do iel = 1, ncel
      propce(iel,iivism) = 1.d0
    enddo
  enddo
endif

!===============================================================================
! 4. INITIALISATION STANDARD DES VARIABLES DE CALCUL
!     On complete ensuite pour les variables turbulentes et les scalaires
!===============================================================================

!     Toutes les variables a 0
do ivar = 1, nvar
  do iel = 1, ncel
    rtp(iel,ivar) = 0.d0
  enddo
enddo

!     On met la pression P* a PRED0
do iphas = 1, nphas
  do iel = 1, ncel
    rtp(iel,ipr(iphas)) = pred0(iphas)
  enddo
enddo

!     Couplage U-P
if(ipucou.eq.1) then
  do iel = 1, ncel
    ra(itpuco+iel-1         ) = 0.d0
    ra(itpuco+iel-1+  ncelet) = 0.d0
    ra(itpuco+iel-1+2*ncelet) = 0.d0
  enddo
endif

!===============================================================================
! 5. INITIALISATION DE K, RIJ ET EPS
!===============================================================================

!  Si UREF n'a pas ete donnee par l'utilisateur ou a ete mal initialisee
!    (valeur negative), on met les valeurs de k, Rij, eps et omega a
!    -10*GRAND. On testera ensuite si l'utilisateur les a modifiees dans
!    usiniv ou en lisant un fichier suite.

do iphas = 1, nphas

  if(itytur(iphas).eq.2 .or. iturb(iphas).eq.50) then

    ikiph  = ik (iphas)
    ieiph  = iep(iphas)

    xcmu = cmu
    if (iturb(iphas).eq.50) xcmu = cv2fmu


    if (uref(iphas).ge.0.d0) then
      do iel = 1, ncel
        rtp(iel,ikiph) = 1.5d0*(0.02d0*uref(iphas))**2
        rtp(iel,ieiph) = rtp(iel,ikiph)**1.5d0*xcmu/almax(iphas)
      enddo

      iclip = 1
      iphass = iphas
      call clipke(ncelet , ncel   , nvar    , nphas  ,            &
                  iphass , iclip  , iwarni(ikiph),                &
                  propce , rtp    )

    else
      do iel = 1, ncel
        rtp(iel,ikiph) = -grand
        rtp(iel,ieiph) = -grand
      enddo
    endif

    if (iturb(iphas).eq.50) then
      iphiph = iphi(iphas)
      ifbiph = ifb (iphas)
      do iel = 1, ncel
        rtp(iel,iphiph) = 2.d0/3.d0
        rtp(iel,ifbiph) = 0.d0
      enddo
    endif

  elseif(itytur(iphas).eq.3) then

    ir11ip = ir11(iphas)
    ir22ip = ir22(iphas)
    ir33ip = ir33(iphas)
    ir12ip = ir12(iphas)
    ir13ip = ir13(iphas)
    ir23ip = ir23(iphas)
    ieiph  = iep (iphas)

    if (uref(iphas).ge.0.d0) then

      trii   = (0.02d0*uref(iphas))**2

      do iel = 1, ncel
        rtp(iel,ir11ip) = trii
        rtp(iel,ir22ip) = trii
        rtp(iel,ir33ip) = trii
        rtp(iel,ir12ip) = 0.d0
        rtp(iel,ir13ip) = 0.d0
        rtp(iel,ir23ip) = 0.d0
        xxk = 0.5d0*(rtp(iel,ir11ip)+                             &
                     rtp(iel,ir22ip)+rtp(iel,ir33ip))
        rtp(iel,ieiph) = xxk**1.5d0*cmu/almax(iphas)
      enddo
      iclip = 1
      iphass = iphas
      call clprij(ncelet , ncel   , nvar    , nphas  ,            &
                  iphass , iclip  ,                               &
                  propce , rtp    , rtp    )

    else

      do iel = 1, ncel
        rtp(iel,ir11ip) = -grand
        rtp(iel,ir22ip) = -grand
        rtp(iel,ir33ip) = -grand
        rtp(iel,ir12ip) = -grand
        rtp(iel,ir13ip) = -grand
        rtp(iel,ir23ip) = -grand
        rtp(iel,ieiph)  = -grand
      enddo

    endif

  elseif(iturb(iphas).eq.60) then

    ikiph   = ik  (iphas)
    iomgip  = iomg(iphas)

    if (uref(iphas).ge.0.d0) then

      do iel = 1, ncel
        rtp(iel,ikiph ) = 1.5d0*(0.02d0*uref(iphas))**2
!     on utilise la formule classique eps=k**1.5/Cmu/ALMAX et omega=eps/Cmu/k
        rtp(iel,iomgip) = rtp(iel,ikiph)**0.5d0/almax(iphas)
      enddo
!     pas la peine de clipper, les valeurs sont forcement positives

    else

      do iel = 1, ncel
        rtp(iel,ikiph ) = -grand
        rtp(iel,iomgip) = -grand
      enddo

    endif

  endif

enddo

!===============================================================================
! 6.  CLIPPING DES GRANDEURS SCALAIRES (SF K-EPS VOIR CI DESSUS)
!===============================================================================

if (nscal.gt.0) then

!    Clipping des scalaires non variance
  do iis = 1, nscal
    if(iscavr(iis).eq.0) then
      iscal = iis
      call clpsca                                                 &
      !==========
     ( ncelet , ncel   , nvar   , nscal  , iscal  ,               &
       propce , ra     , rtp    )
    endif
  enddo

!     Clipping des variances qui sont clippees sans recours au scalaire
!        associe
  do iis = 1, nscal
    if(iscavr(iis).ne.0.and.iclvfl(iis).ne.1) then
      iscal = iis
      call clpsca                                                 &
      !==========
     ( ncelet , ncel   , nvar   , nscal  , iscal  ,               &
       propce , ra     , rtp    )
    endif
  enddo

!     Clipping des variances qui sont clippees avec recours au scalaire
!        associe s'il est connu
  do iis = 1, nscal
    if( iscavr(iis).le.nscal.and.iscavr(iis).ge.1.and.            &
        iclvfl(iis).eq.1 ) then
      iscal = iis
      call clpsca                                                 &
      !==========
     ( ncelet , ncel   , nvar   , nscal  , iscal  ,               &
       propce , rtp(1,isca(iscavr(iis))) , rtp      )
    endif
  enddo

endif

!===============================================================================
! 7.  INITIALISATION DE CONDITIONS AUX LIMITES ET FLUX DE MASSE
!      NOTER QUE LES CONDITIONS AUX LIMITES PEUVENT ETRE UTILISEES DANS
!      PHYVAR, PRECLI
!===============================================================================

! Conditions aux limites
do ii = 1, ncofab
  do ifac = 1, nfabor
    coefa(ifac,ii) = 0.d0
    coefb(ifac,ii) = 1.d0
  enddo
enddo

do iphas = 1, nphas
  do ifac = 1, nfabor
    ia(iitypf-1+ifac+nfabor*(iphas-1)) = 0
    ia(iitrif-1+ifac+nfabor*(iphas-1)) = 0
  enddo
enddo

! Type symétrie : on en a besoin dans le cas du calcul des gradients
!     par moindres carrés étendu avec extrapolation du gradient au bord
!     La valeur 0 permet de ne pas extrapoler le gradient sur les faces.
!     Habituellement, on évite l'extrapolation sur les faces de symétries
!     pour ne pas tomber sur une indétermination et une matrice 3*3 non
!     inversible dans les configurations 2D).
do iphas = 1, nphas
  do ifac = 1, nfabor
    ia(iisymp-1+ifac+nfabor*(iphas-1)) = 0
  enddo
enddo

! Flux de masse (on essaye de ne pas trop faire les choses 2 fois,
!  sans toutefois faire des tests trop compliques)
iiflum = 0
do ivar = 1, nvar
  if(ifluma(ivar).gt.0.and.ifluma(ivar).ne.iiflum) then
    iiflum = ifluma(ivar)
    do ifac = 1, nfac
      propfa(ifac,ipprof(iiflum)) = 0.d0
    enddo
    do ifac = 1, nfabor
      propfb(ifac,ipprob(iiflum)) = 0.d0
    enddo
  endif
enddo

! Flux de masse "ancien" : on utilise le fait  que IFLUMAA = -1 si
!     le flux ancien n'est pas defini (voir le test dans varpos).

iiflua = 0
do ivar = 1, nvar
  if(ifluaa(ivar).gt.0.and.ifluaa(ivar).ne.iiflua) then
    iiflua = ifluaa(ivar)
    iiflum = ifluma(ivar)
    do ifac = 1, nfac
      propfa(ifac,ipprof(iiflua)) = propfa(ifac,ipprof(iiflum))
    enddo
    do ifac = 1, nfabor
      propfb(ifac,ipprob(iiflua)) = propfb(ifac,ipprob(iiflum))
    enddo
  endif
enddo


!===============================================================================
! 8.  INITIALISATION DES TERMES SOURCES SI EXTRAPOLES
!===============================================================================

do iphas = 1, nphas

!     les termes sources de Navier Stokes
  if(isno2t(iphas).gt.0) then
    iptsna = ipproc(itsnsa(iphas))
    do ii = 1, ndim
      do iel = 1, ncel
        propce(iel,iptsna+ii-1) = 0.d0
      enddo
    enddo
  endif

!     les termes sources turbulents
  if(isto2t(iphas).gt.0) then
    if(itytur(iphas).eq.2) jj = 2
    if(itytur(iphas).eq.3) jj = 7
    if(iturb(iphas).eq.50) jj = 4
    if(iturb(iphas).eq.60) jj = 2
    iptsta = ipproc(itstua(iphas))
    do ii = 1, jj
      do iel = 1, ncel
        propce(iel,iptsta+ii-1) = 0.d0
      enddo
    enddo
  endif

enddo

!     les termes sources des scalaires
do iis = 1, nscal
  if(isso2t(iis).gt.0) then
    iptsca = ipproc(itssca(iis))
    do iel = 1, ncel
      propce(iel,iptsca) = 0.d0
    enddo
  endif
enddo

!===============================================================================
! 9.  INITIALISATION DES MOYENNES
!===============================================================================

do imom = 1, nbmomt
  do iel = 1, ncel
    propce(iel,ipproc(icmome(imom))) = 0.d0
  enddo
enddo
do ii = 1,  nbdtcm
  do iel = 1, ncel
    propce(iel,ipproc(icdtmo(ii))) = 0.d0
  enddo
enddo
do ii = 1,  nbmomx
  dtcmom(ii) = 0.d0
enddo

!===============================================================================
! 10.  INITIALISATION CONSTANTE DE SMAGORINSKY EN MODELE DYNAMIQUE
!===============================================================================

do iphas = 1, nphas
  if(iturb(iphas).eq.41) then
    do iel = 1, ncel
      propce(iel,ipproc(ismago(iphas))) = 0.d0
    enddo
  endif
enddo

!===============================================================================
! 11.  INITIALISATION DU NUMERO DE LA FACE DE PAROI 5 LA PLUS PROCHE
!===============================================================================

!     Si IFAPA existe,
!     on suppose qu'il faut le (re)calculer : on init le tab a -1.

do iphas = 1, nphas
  if(iifapa(iphas).gt.0) then
    do iel = 1, ncel
      ia(iifapa(iphas)-1+ iel) = -1
    enddo
  endif
enddo

!===============================================================================
! 12.  INITIALISATION DE LA FORCE EXTERIEURE QUAND IPHYDR=1
!===============================================================================

if(iphydr.eq.1) then
  do iphas = 1, nphas
    do iel = 1, ncel
      frcxt(iel,1,iphas) = 0.d0
      frcxt(iel,2,iphas) = 0.d0
      frcxt(iel,3,iphas) = 0.d0
    enddo
  enddo
endif

!===============================================================================
! 13.  INITIALISATIONS EN ALE OU MAILLAGE MOBILE
!===============================================================================

if (iale.eq.1) then
  do ii = 1, nnod
    ia(iimpal+ii-1) = 0
    do idim = 1, 3
      ra(idepal+(idim-1)*nnod+ii-1) = 0.d0
    enddo
  enddo
endif

if (iale.eq.1.or.imobil.eq.1) then
  do ii = 1, nnod
    do idim = 1, 3
      ra(ixyzn0+(ii-1)*ndim+idim-1) = xyznod(idim,ii)
    enddo
  enddo
endif

!----
! FIN
!----

return
end subroutine
