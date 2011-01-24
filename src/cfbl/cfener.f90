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

subroutine cfener &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR L'ENERGIE TOTALE SUR UN PAS DE TEMPS

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
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
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
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
use entsor
use optcal
use cstphy
use cstnum
use pointe
use parall
use period
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbrs(ncelet)
double precision rovsdt(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          ifinia, ifinra
integer          ivar  , iphas
integer          ifac  , iel
integer          init  , isqrt , iii
integer          iclvar, iclvaf
integer          ipcrom, ipcvst, ipcvsl, iflmas, iflmab
integer          ippvar, ipp
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
double precision epsrgp, climgp, extrap, blencp, epsilp
double precision sclnor, thetap, epsrsp

integer          iwb    , inc    , iccocg , icoefa , icoefb
integer          ivar0  , iphydp , iij , ii , jj
integer          iccfth , imodif
integer          iel1  , iel2, iifru, iifbe
integer          iterns
integer          maxelt, ils, idbia1
double precision flux
double precision dijpfx, dijpfy, dijpfz, pond  , pip   , pjp
double precision diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz
!      DOUBLE PRECISION FLUI  , FLUJ

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! ---  Numero de phase associee au scalaire traite
iphas  = iphsca(iscal)

! --- Numero de variable de calcul et de post associe au scalaire traite
ivar   = isca(iscal)
ippvar = ipprtp(ivar)

! --- Numero des conditions aux limites
iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

! --- Numero des grandeurs physiques
ipcrom = ipproc(irom  (iphas))
ipcvst = ipproc(ivisct(iphas))
iflmas = ipprof(ifluma(ivar ))
iflmab = ipprob(ifluma(ivar ))
if(ivisls(iscal).gt.0) then
  ipcvsl = ipproc(ivisls(iscal))
else
  ipcvsl = 0
endif

! --- Indicateur flux de bord Rusanov
if(iifbru.gt.0) then
  iifru = iifbru+(iphas-1)*nfabor
else
  iifru = 1
endif

! --- Indicateur flux conductif de bord imposé
if(iifbet.gt.0) then
  iifbe = iifbet+(iphas-1)*nfabor
else
  iifbe = 1
endif

! --- Impressions
chaine = nomvar(ippvar)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

! --- Reservation de la memoire

call memcfe                                                       &
!==========
 ( idebia , idebra ,                                              &
   iwb    ,                                                       &
   ifinia , ifinra )

idebia = ifinia
idebra = ifinra

!===============================================================================
! 2. TERMES SOURCES
!===============================================================================

! --> Theta-schema de resolution

! Pour l'instant on prend THETA=1 et on ne code pas le theta-schema

! --> Initialisation

do iel = 1, ncel
  smbrs(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo


!     TERME SOURCE VOLUMIQUE DE CHALEUR : RHO*PHI *VOLUME
!     =================================          v

maxelt = max(ncelet,nfac,nfabor)
ils = idebia
idbia1 = ils + maxelt
call iasize('cfener',idbia1)

call ustssc                                                       &
!==========
 ( idbia1 , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iscal  ,                                                       &
   maxelt , ia(ils),                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt ,                                              &
!        ------   ------
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , drtp   , dam    ,          &
   ra     )

do iel = 1, ncel
  smbrs(iel) = smbrs(iel) + rovsdt(iel)*rtp(iel,ivar)
  rovsdt(iel) = max(-rovsdt(iel),zero)
enddo


!     TERMES DE SOURCE DE MASSE
!     =========================

!     GAMMA(IEL) = SMACEL(IEL,IPR(IPHAS))

!     Terme implicite : GAMMA*VOLUME
!                                                        n
!     Terme explicite : GAMMA*VOLUME*e   - GAMMA*VOLUME*e
!                                     inj
if (ncesmp.gt.0) then
  iterns = 1
  call catsma ( ncelet , ncel , ncesmp , iterns ,                 &
                isno2t(iphas), thetav(ivar),                      &
                icetsm , itypsm(1,ivar) ,                         &
                volume , rtpa(1,ivar) , smacel(1,ivar) ,          &
                smacel(1,ipr(iphas)) , smbrs , rovsdt , w1)
endif

!                                     __    n+1
!     TERME D'ACCUMULATION DE MASSE : >  (Q    .n)  *S
!     =============================   --   pr     ij  ij

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w1)

!                                      __    n+1             n
!     TERME INSTATIONNAIRE EXPLICITE : >  (Q    .n)  *S   * e
!     ==============================   --   pr     ij  ij

do iel = 1, ncel
  smbrs(iel) = smbrs(iel) + iconv(ivar)*w1(iel)*rtpa(iel,ivar)
enddo

!                                      RHO*VOLUME   __    n+1
!     TERME INSTATIONNAIRE IMPLICITE : ---------- - >  (Q    .n)  *S
!     ==============================       DT       --   pr     ij  ij

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel) - iconv(ivar)*w1(iel)                 &
           + istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel)
enddo

!                                       __        v
!     TERME DE DISSIPATION VISQUEUSE  : >  ((SIGMA *U).n)  *S
!     ==============================    --               ij  ij

if( idiff(iu(iphas)).ge. 1 ) then
!                             ^^^
  call cfdivs                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rtp(1,iu(iphas)), rtp(1,iv(iphas)), rtp(1,iw(iphas)), &
!        ------
   w9     ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

endif


!                                         __   P        n+1
!     TERME DE TRANSPORT DE PRESSION  : - >  (---)  *(Q    .n)  *S
!     ==============================      --  RHO ij   pr     ij  ij


if(igrdpp(iphas).gt.0) then
  do iel = 1, ncel
    w9(iel) = rtp(iel,isca(irho(iphas)))
  enddo
else
  do iel = 1, ncel
    w9(iel) = rtpa(iel,isca(irho(iphas)))
  enddo
endif

!     Avec Reconstruction : ca pose probleme pour l'instant


!   Calcul du gradient de P/RHO

!      do iel = 1, ncel
!        w7(iel) = rtp(iel,ipr(iphas))/w9(iel)
!      enddo

! Rq : A defaut de connaitre les parametres pour P/RHO on prend ceux de P

!      iii = ipr(iphas)
!      inc = 1
!      iccocg = 1
!      nswrgp = nswrgr(iii)
!      imligp = imligr(iii)
!      iwarnp = iwarni(iii)
!      epsrgp = epsrgr(iii)
!      climgp = climgr(iii)
!      extrap = extrag(iii)

!       On alloue localement 2 tableaux de NFABOR pour le calcul
!       de COEFA et COEFB de P/RHO

!      icoefa = idebra
!      icoefb = icoefa + nfabor
!      ifinra = icoefb + nfabor
!      call rasize ('cfener',ifinra)
!      !==========

!      do ifac = 1, nfabor
!        ra(icoefa+ifac-1) = zero
!        ra(icoefb+ifac-1) = 1.d0
!      enddo

! En periodique et parallele, echange avant calcul du gradient
!      if (irangp.ge.0.or.iperio.eq.1) then
!        call synsca(w7)
!        !==========
!      endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
!      ivar0 = 0
!      iphydp = 0
!      call grdcel
!      !==========
!     & ( idebia , ifinra ,
!     &   nphas  ,
!     &   ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp ,
!     &   iwarnp , nfecra , epsrgp , climgp , extrap ,
!     &   ia     ,
!     &   w7     , w7     , w7     ,
!     &   w7     , ra(icoefa) , ra(icoefb)  ,
!     &   w1     , w2     , w3     ,
!     &   w4     , w5     , w6     ,
!     &   ra     )

!       On libere la place dans RA

!      ifinra = idebra

!     Faces internes
!      do ifac = 1, nfac

!        ii = ifacel(1,ifac)
!        jj = ifacel(2,ifac)

!        iij = idijpf-1+3*(ifac-1)
!        dijpfx = ra(iij+1)
!        dijpfy = ra(iij+2)
!        dijpfz = ra(iij+3)

!        pond   = ra(ipond-1+ifac)

!        Calcul II' et JJ'

!        diipfx = cdgfac(1,ifac) - (xyzcen(1,ii)+
!     &           (1.d0-pond) * dijpfx)
!        diipfy = cdgfac(2,ifac) - (xyzcen(2,ii)+
!     &           (1.d0-pond) * dijpfy)
!        diipfz = cdgfac(3,ifac) - (xyzcen(3,ii)+
!     &           (1.d0-pond) * dijpfz)
!        djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj)+
!     &               pond  * dijpfx
!        djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj)+
!     &               pond  * dijpfy
!        djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj)+
!     &               pond  * dijpfz

!        pip = w7(ii)
!     &       +w1(ii)*diipfx+w2(ii)*diipfy+w3(ii)*diipfz

!        pjp = w7(jj)
!     &       +w1(jj)*djjpfx+w2(jj)*djjpfy+w3(jj)*djjpfz

!        flui = (propfa(ifac,iflmas)+abs(propfa(ifac,iflmas)))
!        fluj = (propfa(ifac,iflmas)-abs(propfa(ifac,iflmas)))

!        viscf(ifac) = -(pond*pip*flui+pond*pjp*fluj)

!      enddo

!     Sans Reconstruction

!     En periodique et parallele, echange avant utilisation
!       des valeurs aux faces
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w9)
  !==========
endif

!     Faces internes
do ifac = 1, nfac
  iel1 = ifacel(1,ifac)
  iel2 = ifacel(2,ifac)
  viscf(ifac) =                                                   &
     - rtp(iel1,ipr(iphas))/w9(iel1)                              &
     *0.5d0*( propfa(ifac,iflmas) +abs(propfa(ifac,iflmas)) )     &
     - rtp(iel2,ipr(iphas))/w9(iel2)                              &
     *0.5d0*( propfa(ifac,iflmas) -abs(propfa(ifac,iflmas)) )
enddo

!     Faces de bord : pour les faces ou on a calcule un flux de Rusanov,
!       on remplace la contribution standard par le flux de Rusanov qui
!       contient tous les flux convectifs (et il faudra donc eliminer le
!       flux convectif dans cfbsc2)

if(iifbru.gt.0) then

  do ifac = 1, nfabor
    if(ia(iifru+ifac-1).eq.0) then

      iel = ifabor(ifac)
      viscb(ifac) = - propfb(ifac,iflmab)                         &
    * ( coefa(ifac,iclrtp(ipr(iphas),icoef))                      &
      + coefb(ifac,iclrtp(ipr(iphas),icoef))*rtp(iel,ipr(iphas)) )&
    / ( coefa(ifac,iclrtp(isca(irho(iphas)),icoef))               &
      + coefb(ifac,iclrtp(isca(irho(iphas)),icoef))*w9(iel) )

    else
      viscb(ifac) = - propfb(ifac,ipprob(ifbene(iphas)))
    endif
  enddo

else
  do ifac = 1, nfabor

    iel = ifabor(ifac)
    viscb(ifac) = - propfb(ifac,iflmab)                           &
    * ( coefa(ifac,iclrtp(ipr(iphas),icoef))                      &
      + coefb(ifac,iclrtp(ipr(iphas),icoef))*rtp(iel,ipr(iphas)) )&
    / ( coefa(ifac,iclrtp(isca(irho(iphas)),icoef))               &
      + coefb(ifac,iclrtp(isca(irho(iphas)),icoef))*w9(iel) )

  enddo
endif

!     Divergence
init = 0
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
              ifacel,ifabor,viscf,viscb,smbrs)


!     TERME DE FORCES DE PESANTEUR : RHO*g.U *VOLUME
!     ============================

do iel = 1, ncel
  smbrs(iel) = smbrs(iel) + w9(iel)*volume(iel)                   &
                           *( gx*rtp(iel,iu(iphas))               &
                            + gy*rtp(iel,iv(iphas))               &
                            + gz*rtp(iel,iw(iphas)) )
enddo

!                                       Kij*Sij           LAMBDA   Cp   MUT
!     "VITESSE" DE DIFFUSION FACETTE : --------- avec K = ------ + -- .------
!     ==============================    IJ.nij              Cv     Cv  SIGMAS

if( idiff(ivar).ge. 1 ) then

!     MUT/SIGMAS
  do iel = 1, ncel
    w1(iel) = propce(iel,ipcvst)/sigmas(iscal)
  enddo
!     CP*MUT/SIGMAS
  if(icp(iphas).gt.0) then
    do iel = 1, ncel
      w1(iel) = w1(iel)*propce(iel,ipproc(icp(iphas)))
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel)*cp0(iphas)
    enddo
  endif
!     (CP/CV)*MUT/SIGMAS
  if(icv(iphas).gt.0) then
    do iel = 1, ncel
      w1(iel) = w1(iel)/propce(iel,ipproc(icv(iphas)))
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel)/cv0(iphas)
    enddo
  endif
!     (CP/CV)*MUT/SIGMAS+LAMBDA/CV
  if(ipcvsl.eq.0)then
    do iel = 1, ncel
      w1(iel) = w1(iel) + visls0(iscal)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel) + propce(iel,ipcvsl)
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


!     TERME DIFFUSIF COMPLEMENTAIRE : - div( K grad ( epsilon - Cv.T ) )
!     =============================                   1  2
!                                     - div( K grad ( -.u  ) )
!                                                     2

!     Terme complementaire au centre des cellules
  iccfth = 7
  imodif = 0
  call uscfth                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iccfth , imodif , iphas  ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w9     , ra(iwb), w8     , w1     ,                            &
!        ------   -------
   ra     )

!     Calcul de la divergence avec reconstruction


!   Calcul du gradient de (0.5*u*u+EPSILONsup)


do iel = 1, ncel
  w7(iel) =0.5d0*( rtp(iel,iu(iphas))**2                          &
                  +rtp(iel,iv(iphas))**2                          &
                  +rtp(iel,iw(iphas))**2 ) + w9(iel)
enddo

! Rq : A defaut de connaitre les parametres, on prend ceux de la Vitesse

iii = iu(iphas)
inc = 1
iccocg = 1
nswrgp = nswrgr(iii)
imligp = imligr(iii)
iwarnp = iwarni(iii)
epsrgp = epsrgr(iii)
climgp = climgr(iii)
extrap = extrag(iii)

!       On alloue localement 2 tableaux de NFABOR pour le calcul
!       de COEFA et COEFB

icoefa = idebra
icoefb = icoefa + nfabor
ifinra = icoefb + nfabor
call rasize ('cfener',ifinra)
!==========

do ifac = 1, nfabor
  ra(icoefa+ifac-1) = zero
  ra(icoefb+ifac-1) = 1.d0
enddo

! En periodique et parallele, echange avant calcul du gradient
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w7)
  !==========
endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
ivar0 = 0
iphydp = 0
call grdcel                                                       &
!==========
 ( idebia , ifinra ,                                              &
   nphas  ,                                                       &
   ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w7     , w7     , w7     ,                                     &
   w7     , ra(icoefa) , ra(icoefb)  ,                            &
   w1     , w2     , w3     ,                                     &
   w4     , w5     , w6     ,                                     &
   ra     )

!       On libere la place dans RA

ifinra = idebra

!     Faces internes

do ifac = 1, nfac

  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)

  iij = idijpf-1+3*(ifac-1)
  dijpfx = ra(iij+1)
  dijpfy = ra(iij+2)
  dijpfz = ra(iij+3)

  pond   = ra(ipond-1+ifac)

!        Calcul II' et JJ'

  diipfx = cdgfac(1,ifac) - (xyzcen(1,ii)+                        &
           (1.d0-pond) * dijpfx)
  diipfy = cdgfac(2,ifac) - (xyzcen(2,ii)+                        &
           (1.d0-pond) * dijpfy)
  diipfz = cdgfac(3,ifac) - (xyzcen(3,ii)+                        &
           (1.d0-pond) * dijpfz)
  djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj)+                        &
               pond  * dijpfx
  djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj)+                        &
               pond  * dijpfy
  djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj)+                        &
               pond  * dijpfz

  pip = w7(ii)                                                    &
       +w1(ii)*diipfx+w2(ii)*diipfy+w3(ii)*diipfz

  pjp = w7(jj)                                                    &
       +w1(jj)*djjpfx+w2(jj)*djjpfy+w3(jj)*djjpfz

  flux = viscf(ifac)*(pip-pjp)

  smbrs(ii) = smbrs(ii) - flux
  smbrs(jj) = smbrs(jj) + flux

enddo


!       Assemblage a partir des facettes de bord
!     Pour les faces à flux imposé ou temperature imposée, tout est
!       pris par le terme de diffusion de l'energie. On ne doit donc
!       pas prendre en compte la contribution des termes en u2 et e-CvT
!       quand IA(IIFBE+IFAC-1).NE.0

  if(iifbet.gt.0) then
    do ifac = 1, nfabor

      if(ia(iifbe+ifac-1).eq.0) then
        iel = ifabor(ifac)

        flux = viscb(ifac)*( w9(iel) - ra(iwb +ifac-1)            &
           + 0.5d0*( rtp(iel,iu(iphas))**2 -                      &
   ( coefa(ifac,iclrtp(iu(iphas),icoef))                          &
   + coefb(ifac,iclrtp(iu(iphas),icoef))*rtp(iel,iu(iphas)) )**2  &
                   + rtp(iel,iv(iphas))**2 -                      &
   ( coefa(ifac,iclrtp(iv(iphas),icoef))                          &
   + coefb(ifac,iclrtp(iv(iphas),icoef))*rtp(iel,iv(iphas)) )**2  &
                   + rtp(iel,iw(iphas))**2 -                      &
   ( coefa(ifac,iclrtp(iw(iphas),icoef))                          &
   + coefb(ifac,iclrtp(iw(iphas),icoef))*rtp(iel,iw(iphas)) )**2))

        smbrs(iel) = smbrs(iel) - flux

      endif

    enddo

!     Sinon : meme code, mais sans le test
  else

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      flux = viscb(ifac)*( w9(iel) - ra(iwb +ifac-1)              &
           + 0.5d0*( rtp(iel,iu(iphas))**2 -                      &
   ( coefa(ifac,iclrtp(iu(iphas),icoef))                          &
   + coefb(ifac,iclrtp(iu(iphas),icoef))*rtp(iel,iu(iphas)) )**2  &
                   + rtp(iel,iv(iphas))**2 -                      &
   ( coefa(ifac,iclrtp(iv(iphas),icoef))                          &
   + coefb(ifac,iclrtp(iv(iphas),icoef))*rtp(iel,iv(iphas)) )**2  &
                   + rtp(iel,iw(iphas))**2 -                      &
   ( coefa(ifac,iclrtp(iw(iphas),icoef))                          &
   + coefb(ifac,iclrtp(iw(iphas),icoef))*rtp(iel,iw(iphas)) )**2))

      smbrs(iel) = smbrs(iel) - flux

    enddo
  endif

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

!===============================================================================
! 4. RESOLUTION
!===============================================================================

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
nitmap = nitmax(ivar)
ndircp = ndircl(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
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
thetap = thetav(ivar)
iescap = 0

call cfcdts                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iscal  , iconvp , idiffp , ireslp , ndircp , nitmap ,          &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , iifbru ,                            &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap , thetap , &
   ia(iifru) ,                                                    &
   ia     ,                                                       &
   rtpa(1,ivar)    , coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbrs  , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )

!===============================================================================
! 5. IMPRESSIONS ET CLIPPINGS
!===============================================================================

! Valeur bidon
  iii = 1

call clpsca                                                       &
!==========
 ( ncelet , ncel   , nvar   , nscal  , iscal  ,                   &
   propce , rtp(1,iii)      , rtp    )

! --- Traitement utilisateur pour gestion plus fine des bornes
!       et actions correctives éventuelles.
  iccfth = -4
  imodif = 0
  call uscfth                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iccfth , imodif , iphas  ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w6     , w7     , w8     , w9     ,                            &
   ra     )


! --- Bilan explicite (voir codits : on enleve l'increment)

if (iwarni(ivar).ge.2) then
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                       &
            - istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel)&
                *(rtp(iel,ivar)-rtpa(iel,ivar))                   &
                * max(0,min(nswrsm(ivar)-2,1))
  enddo
  isqrt = 1
  call prodsc(ncelet,ncel,isqrt,smbrs,smbrs,sclnor)
  write(nfecra,1200)chaine(1:8) ,sclnor
endif

!===============================================================================
! 6. ACTUALISATION FINALE DE LA PRESSION (et calcul de la température)
!===============================================================================
!                               n+1      n+1  n+1
! On utilise l'equation d'etat P   =P(RHO   ,H   )

! --- Calcul de P et T au centre des cellules
  iccfth = 24
  imodif = 0
  call uscfth                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iccfth , imodif , iphas  ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   rtp(1,ipr(iphas)) , rtp(1,isca(itempk(iphas))) , w8     , w9 , &
!        -----------------   --------------------------
   ra     )

!===============================================================================
! 7. COMMUNICATION DE LA PRESSION, DE L'ENERGIE ET DE LA TEMPERATURE
!===============================================================================

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtp(1,ipr(iphas)))
  !==========
  call synsca(rtp(1,ivar))
  !==========
  call synsca(rtp(1,isca(itempk(iphas))))
  !==========
endif


!--------
! FORMATS
!--------

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A8                        ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)

!----
! FIN
!----

return

end subroutine
