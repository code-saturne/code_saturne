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

subroutine covofi &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  , itspdv ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb , tslagr , &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt ,                                              &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR UN SCALAIRE SUR UN PAS DE TEMPS

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! itspdv           ! e  ! <-- ! calcul termes sources prod et dissip           !
!                  !    !     !  (0 : non , 1 : oui)                           !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
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
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
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
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use fuincl
use ppincl
use lagpar
use lagran
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal  , itspdv

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          ivar
integer          ifac  , iel
integer          init  , inc   , iccocg, isqrt, iii, iiun, ibcl
integer          ivarsc, iscala
integer          iiscav, iicp
integer          iclvar, iclvaf
integer          ipcrom, ipcvst, ipcvsl, iflmas, iflmab
integer          ippvar, ipp   , iptsca, ipcvso
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          idbia1

double precision epsrgp, climgp, extrap, relaxp, blencp, epsilp
double precision epsrsp
double precision rhovst, xk    , xe    , sclnor
double precision thetv , thets , thetap, thetp1
double precision smbexp

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:,:) :: grad

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(w1(ncelet))

! Initialize variables to avoid compiler warnings

xe = 0.d0
xk = 0.d0

! --- Memoire

idebia = idbia0
idebra = idbra0

! --- Numero de variable de calcul et de post associe au scalaire traite
ivar   = isca(iscal)
ippvar = ipprtp(ivar)

! --- Numero du scalaire eventuel associe dans le cas fluctuation
!         et numero de variable de calcul
iiscav = iscavr(iscal)
if(iiscav.gt.0.and.iiscav.le.nscal) then
  ivarsc = isca(iiscav)
else
  ivarsc = 0
endif

! --- Numero des conditions aux limites
iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

! --- Numero des grandeurs physiques
ipcrom = ipproc(irom  )
ipcvst = ipproc(ivisct)
iflmas = ipprof(ifluma(ivar ))
iflmab = ipprob(ifluma(ivar ))
if(ivisls(iscal).gt.0) then
  ipcvsl = ipproc(ivisls(iscal))
else
  ipcvsl = 0
endif

! --- Numero du terme source dans PROPCE si extrapolation
if(isso2t(iscal).gt.0) then
  iptsca = ipproc(itssca(iscal))
else
  iptsca = 0
endif

!     S pour Source, V pour Variable
thets  = thetss(iscal)
thetv  = thetav(ivar )

chaine = nomvar(ippvar)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!===============================================================================
! 2. TERMES SOURCES
!===============================================================================

! --> Initialisation


do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

do iel = 1, ncel
  smbrs(iel) = 0.d0
enddo

iscala = iscal

call ustssc &
!==========
( nvar   , nscal  , ncepdp , ncesmp ,                            &
  iscala ,                                                       &
  icepdc , icetsm , itypsm ,                                     &
  ia     ,                                                       &
  dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
  coefa  , coefb  , ckupdc , smacel ,                            &
  smbrs  , rovsdt ,                                              &
!        ------   ------
  ra     )

!     Si on extrapole les TS :
!       SMBRS recoit -theta PROPCE du pas de temps precedent
!         (on aurait pu le faire avant ustssc, mais avec le risque que
!          l'utilisateur l'ecrase)
!       SMBRS recoit la partie du terme source qui depend de la variable
!       A l'ordre 2, on suppose que le ROVSDT fourni par l'utilisateur est <0
!         on implicite le terme (donc ROVSDT*RTPA va dans SMBRS)
!       En std, on adapte le traitement au signe de ROVSDT, mais ROVSDT*RTPA va
!         quand meme dans SMBRS (pas d'autre choix)
if(isso2t(iscal).gt.0) then
  do iel = 1, ncel
!        Stockage temporaire pour economiser un tableau
    smbexp = propce(iel,iptsca)
!        Terme source utilisateur explicite
    propce(iel,iptsca) = smbrs(iel)
!        Terme source du pas de temps precedent et
!        On suppose -ROVSDT > 0 : on implicite
!           le terme source utilisateur (le reste)
    smbrs(iel) = rovsdt(iel)*rtpa(iel,ivar) - thets*smbexp
!        Diagonale
    rovsdt(iel) = - thetv*rovsdt(iel)
  enddo
!     Si on n'extrapole pas les TS :
else
  do iel = 1, ncel
!        Terme source utilisateur
    smbrs(iel) = smbrs(iel) + rovsdt(iel)*rtpa(iel,ivar)
!        Diagonale
    rovsdt(iel) = max(-rovsdt(iel),zero)
  enddo
endif

! --> Physique particulieres
!     Ordre 2 non pris en compte

if (ippmod(iphpar).ge.1) then
  call pptssc                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscala ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt , tslagr ,                                     &
!        ------   ------
   ra     )
endif

! --> Rayonnement
!     Ordre 2 non pris en compte

if ( iirayo.ge.1 ) then

  if (iscal.eq.iscalt) then
    call raysca                    &
    !==========
  ( iscalt,ncelet,ncel,     &
    smbrs, rovsdt,volume,propce )
  endif

  !-> Charbon pulverise
  !   Ordre 2 non pris en compte

  if ( ippmod(icp3pl) .gt. 0 ) then
    if ( isca(iscal).ge.isca(ih2(1)) .and.       &
         isca(iscal).le.isca(ih2(nclacp)) ) then

      call cprays                  &
      !==========
    ( ivar  ,ncelet, ncel  ,       &
      volume,propce,smbrs,rovsdt)

    endif
  endif

  ! -> Fuel
  !    Ordre 2 non pris en compte
  !    Pour l'instant rayonnement non compatible avec Fuel

  if ( ippmod(icfuel) .ge. 0 ) then
    if ( isca(iscal).ge.isca(ihlf(1)) .and.       &
         isca(iscal).le.isca(ihlf(nclafu)) ) then

      call furays                  &
      !==========
    ( ivar  ,ncelet, ncel  ,       &
      volume,rtpa  , propce,smbrs,rovsdt)

    endif
  endif

endif

! --> Lagrangien (couplage retour thermique)
!     Ordre 2 non pris en compte

if (iilagr.eq.2 .and. ltsthe.eq.1)  then

  if (iscsth(iscal).eq.2) then

!    --> Enthalpie

    do iel = 1,ncel
      smbrs (iel) = smbrs(iel)  + tslagr(iel,itste)
      rovsdt(iel) = rovsdt(iel) + max(tslagr(iel,itsti),zero)
    enddo

  else if (iscsth(iscal).eq.1 .or. iscsth(iscal).eq.-1 ) then

!    --> Temperature  :

    if (icp.eq.0) then

!       --> Cp constant

      do iel = 1,ncel
        smbrs (iel) = smbrs(iel)  + tslagr(iel,itste)/cp0
        rovsdt(iel) = rovsdt(iel) + max(tslagr(iel,itsti),zero)
      enddo

    else if (icp.gt.0) then

!       --> Cp variable

      iicp = ipproc(icp)
      do iel = 1,ncel
        smbrs (iel) = smbrs(iel) + tslagr(iel,itste)              &
                                 / propce(iel,iicp)
        rovsdt(iel) = rovsdt(iel) + max(tslagr(iel,itsti),zero)
      enddo
    endif

  endif

endif


!     TERMES DE SOURCE DE MASSE

if (ncesmp.gt.0) then

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

!       On incremente SMBRS par -Gamma RTPA et ROVSDT par Gamma (*theta)
  call catsma                                                     &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   , isso2t(iscal) , thetv  ,   &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipr), &
   smbrs , rovsdt , w1)

!       Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if(isso2t(iscal).gt.0) then
    do iel = 1, ncel
      propce(iel,iptsca) = propce(iel,iptsca) + w1(iel)
    enddo
!       Sinon on le met directement dans SMBRS
  else
    do iel = 1, ncel
      smbrs(iel) = smbrs(iel) + w1(iel)
    enddo
  endif

endif


!     TERME D'ACCUMULATION DE MASSE -(dRO/dt)*VOLUME

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w1)

!     Extrapolation ou non, le terme d'accumulation de masse va dans SMBRS
do iel = 1, ncel
  smbrs(iel) = smbrs(iel)                                         &
              + iconv(ivar)*w1(iel)*rtpa(iel,ivar)
enddo

!     Extrapolation ou non le terme d'accumulation de masse a la meme forme
!       par coherence avec bilsc2
do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel) - iconv(ivar)*w1(iel)*thetv
enddo


!     TERME INSTATIONNAIRE

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                       &
           + istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel)
enddo



!    SI ON CALCULE LA VARIANCE DES FLUCTUATIONS D'UN SCALAIRE,
!      ON RAJOUTE LES TERMES DE PRODUCTION ET DE DISSIPATION

if (itspdv.eq.1) then

  if(itytur.eq.2.or.itytur.eq.3                     &
       .or.iturb.eq.50 .or. iturb.eq.60) then

    ! Allocate a temporary array for the gradient reconstruction
    allocate(grad(ncelet,3))

    ! Remarque : on a prevu la possibilite de scalaire associe non
    !  variable de calcul, mais des adaptations sont requises

    if(ivarsc.gt.0) then
      iii = ivarsc
    else
      write(nfecra,9000)ivarsc
      call csexit(1)
    endif

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(iii)
    imligp = imligr(iii)
    iwarnp = iwarni(iii)
    epsrgp = epsrgr(iii)
    climgp = climgr(iii)
    extrap = extrag(iii)

    call grdcel                                                   &
    !==========
 ( iii    , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   rtpa(1,iii) , coefa(1,iclrtp(iii,icoef)) ,                     &
                 coefb(1,iclrtp(iii,icoef)) ,                     &
   grad   ,                                                       &
   ra     )

!     Traitement de la production
!     On utilise MAX(PROPCE,ZERO) car en LES dynamique on fait un clipping
!     tel que (mu + mu_t)>0, donc mu_t peut etre negatif et donc
!     potentiellement (lambda/Cp + mu_t/sigma) aussi
!     Ceci ne pose probleme que quand on resout une equation de variance
!     de scalaire avec un modele LES ... ce qui serait curieux mais n'est
!     pas interdit par le code.
!       Si extrapolation : dans PROPCE
    if (isso2t(iscal).gt.0) then
!         On prend la viscosite a l'instant n, meme si elle est extrapolee
      ipcvso = ipcvst
      if(iviext.gt.0) ipcvso = ipproc(ivista)
      do iel = 1, ncel
        propce(iel,iptsca) = propce(iel,iptsca)                   &
             + 2.d0*max(propce(iel,ipcvso),zero)                  &
             *volume(iel)/sigmas(iscal)                           &
             *(grad(iel,1)**2 + grad(iel,2)**2 + grad(iel,3)**2)
      enddo
!       Sinon : dans SMBRS
    else
      ipcvso = ipcvst
      do iel = 1, ncel
        smbrs(iel) = smbrs(iel)                                   &
             + 2.d0*max(propce(iel,ipcvso),zero)                  &
             *volume(iel)/sigmas(iscal)                           &
             *(grad(iel,1)**2 + grad(iel,2)**2 + grad(iel,3)**2)
      enddo
    endif

    ! Free memory
    deallocate(grad)

!     Traitement de la dissipation
    if (isso2t(iscal).gt.0) then
      thetap = thetv
    else
      thetap = 1.d0
    endif
    do iel = 1, ncel
      if(itytur.eq.2 .or. iturb.eq.50) then
        xk     = rtpa(iel,ik)
        xe     = rtpa(iel,iep)
      elseif(itytur.eq.3) then
        xk     =                                                  &
        0.5d0*(rtpa(iel,ir11)+rtpa(iel,ir22)+rtpa(iel,ir33))
        xe     = rtpa(iel,iep)
      elseif(iturb.eq.60) then
        xk     = rtpa(iel,ik)
        xe     = cmu*xk*rtpa(iel,iomg)
      endif
      rhovst = propce(iel,ipcrom)*xe/                             &
                 (xk * rvarfl(iscal))*volume(iel)
!     La diagonale recoit eps/Rk, (*theta eventuellement)
      rovsdt(iel) = rovsdt(iel) + rhovst*thetap
!     SMBRS recoit la dissipation
      smbrs(iel) = smbrs(iel) - rhovst*rtpa(iel,ivar)
    enddo

  endif

endif


if(isso2t(iscal).gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + thetp1 * propce(iel,iptsca)
  enddo
endif

!     "VITESSE" DE DIFFUSION FACETTE

!     On prend le MAX(mu_t,0) car en LES dynamique mu_t peut etre negatif
!     (clipping sur (mu + mu_t)). On aurait pu prendre
!     MAX(K + K_t,0) mais cela autoriserait des K_t negatif, ce qui est
!     considere ici comme non physique.
if( idiff(ivar).ge. 1 ) then
  if(ipcvsl.eq.0)then
    do iel = 1, ncel
      w1(iel) = visls0(iscal)                                     &
         + idifft(ivar)*max(propce(iel,ipcvst),zero)/sigmas(iscal)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvsl)                                &
         + idifft(ivar)*max(propce(iel,ipcvst),zero)/sigmas(iscal)
    enddo
  endif
  call viscfa                                                     &
  !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  ,                                              &
   ra     )

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif


!===============================================================================
! 3. RESOLUTION
!===============================================================================

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
ndircp = ndircl(ivar)
nitmap = nitmax(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
iescap = 0
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
ipp    = ippvar
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbrs  , rtp(1,ivar)     ,                            &
   rvoid  ,                                                       &
   ra     )

!===============================================================================
! 4. IMPRESSIONS ET CLIPPINGS
!===============================================================================

if(ivarsc.gt.0) then
  iii = ivarsc
else
! Valeur bidon
  iii = 1
endif

call clpsca                                                       &
!==========
 ( ncelet , ncel   , nvar   , nscal  , iscal  ,                   &
   propce , rtp(1,iii)      , rtp    )


! BILAN EXPLICITE (VOIR CODITS : ON ENLEVE L'INCREMENT)
! Ceci devrait etre valable avec le theta schema sur les Termes source

if (iwarni(ivar).ge.2) then
  if(nswrsm(ivar).gt.1) then
    ibcl = 1
  else
    ibcl = 0
  endif
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                       &
            - istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel)&
                *(rtp(iel,ivar)-rtpa(iel,ivar))*ibcl
  enddo
  isqrt = 1
  call prodsc(ncelet,ncel,isqrt,smbrs,smbrs,sclnor)
  write(nfecra,1200)chaine(1:8) ,sclnor
endif

! Free memory
deallocate(w1)

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A8                        ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)
 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS COVOFI                          ',/,&
'@    =========                                               ',/,&
'@    IVARSC DOIT ETRE UN ENTIER POSITIF STRICTEMENT          ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(/,                                                   &
'   ** SOLVING VARIABLE ',A8                                   ,/,&
'      ----------------'                                       ,/)
 1200 format(1X,A8,' : EXPLICIT BALANCE = ',E14.5)
 9000 format(                                                           &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ERROR IN COVOFI'                                ,/,&
'@    ========'                                                ,/,&
'@    IVARSC MUST BE A STRICTLY POSITIVE INTEGER'              ,/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run.'                          ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! FIN
!----

return

end subroutine
