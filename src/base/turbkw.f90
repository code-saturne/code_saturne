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

subroutine turbkw &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   tslagr , coefa  , coefb  , ckupdc , smacel ,                   &
   s2kw   , divukw , viscf  , viscb  ,                            &
   dam    , xam    ,                                              &
   drtp   , smbrk  , smbrw  , tinstk , tinstw ,                   &
   xf1    , w1     , w2     , w3     , w4     ,                   &
   w5     , w6     , w7     , w8     , w9     ,                   &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS K-OMEGA SST 1 PHASE INCOMPRESSIBLE OU
! RHO VARIABLE SUR UN PAS DE TEMPS

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
! iphas            ! i  ! <-- ! phase number                                   !
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
! s2kw(ncelet)     ! tr ! <-- ! tableau contenant s2=2.sijsij                  !
! divukw(ncelet    ! tr ! <-- ! divergence de u (calculee par grdcel)          !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr.(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! tinst.(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! xf1(ncelet)      ! tr ! --- ! tableau de travail pour coef f1                !
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
use dimens, only: ndimfb
use numvar
use entsor
use cstnum
use cstphy
use optcal
use lagpar
use lagran
use pointe
use parall
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          iphas

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision s2kw(ncelet), divukw(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet),   smbrk(ncelet),  smbrw(ncelet)
double precision tinstk(ncelet), tinstw(ncelet), xf1(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra, ifinia
integer          iel   , ifac  , init  , inc   , iccocg, ivar
integer          ii, iivar , iiun  , ifacpt
integer          iclipk, iclipw, isqrt
integer          nswrgp, imligp, iphydp
integer          ipriph, ikiph , iomgip
integer          iclikp, iclomg
integer          iclvar, iclvaf
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          ipcrom, ipbrom, ipcvst, ipcvis, iflmas, iflmab
integer          iwarnp, ipp
integer          iptsta
integer          ipcroo, ipbroo, ipcvto, ipcvlo
integer          maxelt, ils
double precision rnorm , d2s3, divp23, epz2
double precision deltk , deltw, a11, a12, a22, a21
double precision unsdet, romvsd
double precision prdtur, xk, xw, xeps, xnu
double precision visct , rom
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision thetp1, thetak, thetaw, thets, thetap, epsrsp
double precision tuexpk, tuexpw
double precision cdkw, xarg1, xxf1, xgamma, xbeta, sigma, produc
double precision var, vrmin, vrmax

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

epz2 = epzero**2

ipriph = ipr (iphas)
ikiph  = ik  (iphas)
iomgip = iomg(iphas)

iclikp = iclrtp(ikiph ,icoef)
iclomg = iclrtp(iomgip,icoef)

ipcrom = ipproc(irom  (iphas))
ipcvst = ipproc(ivisct(iphas))
ipcvis = ipproc(iviscl(iphas))
iflmas = ipprof(ifluma(ikiph))
iflmab = ipprob(ifluma(ikiph))
ipbrom = ipprob(irom  (iphas))

thets  = thetst(iphas)

ipcroo = ipcrom
ipbroo = ipbrom
ipcvto = ipcvst
ipcvlo = ipcvis
if(isto2t(iphas).gt.0) then
  if (iroext(iphas).gt.0) then
    ipcroo = ipproc(iroma(iphas))
    ipbroo = ipprob(iroma(iphas))
  endif
  if(iviext(iphas).gt.0) then
    ipcvto = ipproc(ivista(iphas))
    ipcvlo = ipproc(ivisla(iphas))
  endif
endif

if(isto2t(iphas).gt.0) then
  iptsta = ipproc(itstua(iphas))
else
  iptsta = 0
endif

if(iwarni(ikiph).ge.1) then
  write(nfecra,1000)iphas
endif


!===============================================================================
! 2. CALCUL DE dk/dxj.dw/dxj
!      Le terme est stocke dans         W1
!      En sortie de l'etape on conserve W1
!===============================================================================

iccocg = 1
inc = 1


nswrgp = nswrgr(ikiph)
imligp = imligr(ikiph)
iwarnp = iwarni(ikiph)
epsrgp = epsrgr(ikiph)
climgp = climgr(ikiph)
extrap = extrag(ikiph)
iphydp = 0

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   ikiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w7     , w7     , w7     ,                                     &
   rtpa(1,ikiph)   , coefa(1,iclikp) , coefb(1,iclikp) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w7     , w8     , xf1    ,                                     &
   ra     )


nswrgp = nswrgr(iomgip)
imligp = imligr(iomgip)
iwarnp = iwarni(iomgip)
epsrgp = epsrgr(iomgip)
climgp = climgr(iomgip)
extrap = extrag(iomgip)
iphydp = 0

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   iomgip , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w7     , w7     , w7     ,                                     &
   rtpa(1,iomgip)  , coefa(1,iclomg) , coefb(1,iclomg) ,          &
   w4     , w5     , w6     ,                                     &
!        ------   ------   ------
   w7     , w8     , xf1    ,                                     &
   ra     )

do iel = 1, ncel
  w1(iel) = w1(iel)*w4(iel)+w2(iel)*w5(iel)+w3(iel)*w6(iel)
enddo

!====================================================
! 3. CALCUL DU COEFFICIENT DE PONDERATION F1
!      Le terme est stocke dans         XF1
!      En sortie de l'etape on conserve W1,XF1
!====================================================


if(abs(icdpar).eq.2) then
  do iel = 1, ncel
    ifacpt = ia(iifapa(iphas)-1+iel)
    w2(iel) =                                                     &
         (cdgfbo(1,ifacpt)-xyzcen(1,iel))**2                      &
         +(cdgfbo(2,ifacpt)-xyzcen(2,iel))**2                     &
         +(cdgfbo(3,ifacpt)-xyzcen(3,iel))**2
    w2(iel) = sqrt(w2(iel))
  enddo
else
  do iel = 1, ncel
    w2(iel) =  max(ra(idipar+iel-1),epzero)
  enddo
endif

!     En cas d'ordre 2 on utilise les valeurs en n car le terme en (1-F1)*W1
!     sera dans PROPCE. Du coup, on aura quand meme certaines "constantes"
!     intervenant dans des termes en n+1/2 (ex sigma_k pour la diffusion) calcules
!     a partir de F1 en n -> mais l'effet sur les "constantes" est faible
!     -> a garder en tete si on fait vraiment de l'ordre 2 en temps en k-omega
do iel = 1, ncel
  rom = propce(iel,ipcroo)
  xnu = propce(iel,ipcvlo)/rom
  xk = rtpa(iel,ikiph)
  xw  = rtpa(iel,iomgip)
  cdkw = 2*rom/ckwsw2/xw*w1(iel)
  cdkw = max(cdkw,1.d-20)
  xarg1 = max( sqrt(xk)/cmu/xw/w2(iel),                           &
       500.d0*xnu/xw/w2(iel)**2 )
  xarg1 = min(xarg1,                                              &
       4.d0*rom*xk/ckwsw2/cdkw/w2(iel)**2)
  xf1(iel) = tanh(xarg1**4)
enddo

!===============================================================================
! 4. CALCUL DU TERME DE PRODUCTION
!      Les termes sont stockes dans     TINSTK,TINSTW
!      En sortie de l'etape on conserve W1,XF1,TINSTK,TINSTW
!===============================================================================

d2s3 = 2.d0/3.d0
do iel = 1, ncel
  xk   = rtpa(iel,ikiph)
  xeps = cmu*rtpa(iel,iomgip)*xk
  tinstk(iel) = s2kw(iel) - d2s3*divukw(iel)*divukw(iel)
  tinstw(iel) = propce(iel,ipcvto)*tinstk(iel)                    &
       -d2s3*propce(iel,ipcroo)*xk*divukw(iel)
  tinstk(iel) = min(tinstw(iel),ckwc1*propce(iel,ipcroo)*xeps)
enddo

!===============================================================================
! 4. CALCUL DU TERME DE GRAVITE
!      Les termes sont stockes dans     TINSTK,TINSTW,W2
!      En sortie de l'etape on conserve W1,W2,XF1,TINSTK,TINSTW
!===============================================================================

if(igrake(iphas).eq.1) then

! --- Terme de gravite G = BETA*G*GRAD(SCA)/PRDTUR/RHO
!     Ici on calcule   G =-G*GRAD(RHO)/PRDTUR/RHO

  iccocg = 1
  inc = 1

!     Le choix ci dessous a l'avantage d'etre simple

  nswrgp = nswrgr(ikiph)
  epsrgp = epsrgr(ikiph)
  imligp = imligr(ikiph)
  iwarnp = iwarni(ikiph)
  climgp = climgr(ikiph)
  extrap = extrag(ikiph)

!     Conditions aux limites sur ROM : Dirichlet ROMB
!       On utilise VISCB pour stocker le COEFB relatif a ROM
!       On impose en Dirichlet (COEFA) la valeur ROMB

  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

  iivar = 0

  iphydp = 0
  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   iivar  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w5     , w5     , w5     ,                                     &
   propce(1,ipcroo), propfb(1,ipbroo), viscb  ,                   &
   w2     , w3     , w4     ,                                     &
!        ------   ------   ------
   w5     , w6     , w7     ,                                     &
   ra     )


!      Production et terme de gravite
!        TINSTK=MIN(P,C1*EPS)+G et TINSTW=P+(1-CE3)*G
!        On conserve G dans W2 pour la phase de couplage des termes sources

  if(iscalt(iphas).gt.0.and.nscal.ge.iscalt(iphas)) then
    prdtur = sigmas(iscalt(iphas))
  else
    prdtur = 1.d0
  endif

  do iel = 1, ncel
    w2(iel) = -(w2(iel)*gx+w3(iel)*gy+w4(iel)*gz)/                &
         (propce(iel,ipcroo)*prdtur)
    tinstw(iel)=tinstw(iel)+propce(iel,ipcvto)*max(w2(iel),zero)
    tinstk(iel)=tinstk(iel)+propce(iel,ipcvto)*w2(iel)
  enddo

endif


!===============================================================================
! 5. PRISE EN COMPTE DES TERMES SOURCES UTILISATEURS
!      Les termes sont stockes dans     SMBRK,SMBRW,DAM,W3
!      En sortie de l'etape on conserve W1-3,XF1,TINSTK,TINSTW,SMBRK,SMBRW,DAM
!===============================================================================

do iel = 1, ncel
  smbrk(iel) = 0.d0
  smbrw(iel) = 0.d0
  dam  (iel) = 0.d0
  w3   (iel) = 0.d0
enddo

maxelt = max(ncelet, nfac, nfabor)
ils    = idebia
ifinia = ils + maxelt
call iasize('turbkw',ifinia)

call ustskw                                                       &
!==========
 ( ifinia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   maxelt , ia(ils),                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , s2kw   , divukw ,          &
   w1     , w2     , xf1    ,                                     &
   smbrk  , smbrw  , dam    , w3     ,                            &
!        ------   ------   ------   ------
   viscf  , viscb  , xam    , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , drtp   ,                   &
   ra     )

!===============================================================================
! 6. TERME D'ACCUMULATION DE MASSE -(dRO/dt)*Volume

!      Le terme est stocke dans         W4
!      En sortie de l'etape on conserve W1-4,XF1,TINSTK,TINSTW,SMBRK,SMBRW,DAM
!===============================================================================

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w4)

!===============================================================================
! 7. ON FINALISE LE CALCUL DES TERMES SOURCES

!      Les termes sont stockes dans     SMBRK, SMBRW
!      En sortie de l'etape on conserve SMBRK,SMBRW,DAM,W1-4
!===============================================================================

do iel = 1, ncel

  visct  = propce(iel,ipcvto)
  rom    = propce(iel,ipcroo)
  xk     = rtpa(iel,ikiph)
  xw     = rtpa(iel,iomgip)
  xxf1   = xf1(iel)
  xgamma = xxf1*ckwgm1 + (1.d0-xxf1)*ckwgm2
  xbeta  = xxf1*ckwbt1 + (1.d0-xxf1)*ckwbt2

  smbrk(iel) = smbrk(iel) + volume(iel)*(                         &
       tinstk(iel)                                                &
       -cmu*rom*xw*xk )

  smbrw(iel) = smbrw(iel) + volume(iel)*(                         &
       rom*xgamma/visct*tinstw(iel)                               &
       -xbeta*rom*xw**2                                           &
       +2.d0*rom/xw*(1.d0-xxf1)/ckwsw2*w1(iel) )

enddo

!===============================================================================
! 8. PRISE EN COMPTE DES TERMES D'ACCUMULATION DE MASSE ET
!         DE LA DEUXIEME PARTIE DES TS UTILISATEURS (PARTIE EXPLICITE)
!         STOCKAGE POUR EXTRAPOLATION EN TEMPS
!      On utilise                       SMBRK,SMBRW
!      En sortie de l'etape on conserve SMBRK,SMBRW,DAM,W2-4

!    Remarque : l'extrapolation telle qu'elle est ecrite n'a pas grand
!               sens si IKECOU=1
!===============================================================================

!     Si on extrapole les T.S.
if(isto2t(iphas).gt.0) then

  do iel = 1, ncel

!       Sauvegarde pour echange
    tuexpk = propce(iel,iptsta)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta) = smbrk(iel)
!       Termes dependant de la variable resolue et theta PROPCE
    smbrk(iel) = iconv(ikiph)*w4(iel)*rtpa(iel,ikiph)             &
         - thets*tuexpk
!       On suppose -DAM > 0 : on implicite
!         le terme utilisateur dependant de la variable resolue
    smbrk(iel) = dam(iel)*rtpa(iel,ikiph) + smbrk(iel)

!       Sauvegarde pour echange
    tuexpw = propce(iel,iptsta+1)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta+1) = smbrw(iel)
!       Termes dependant de la variable resolue et theta PROPCE
    smbrw(iel) = iconv(iomgip)*w4(iel)*rtpa(iel,iomgip)           &
         - thets*tuexpw
!       On suppose -W3 > 0 : on implicite
!         le terme utilisateur dependant de la variable resolue
    smbrw(iel) =  w3(iel)*rtpa(iel,iomgip) + smbrw(iel)

  enddo

!     Si on n'extrapole pas les T.S.
else
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + dam(iel)*rtpa(iel,ikiph)            &
         +iconv(ikiph)*w4(iel)*rtpa(iel,ikiph)
    smbrw(iel) = smbrw(iel) + w3 (iel)*rtpa(iel,iomgip)           &
         +iconv(iomgip)*w4(iel)*rtpa(iel,iomgip)
  enddo
endif

!===============================================================================
! 8.1 PRISE EN COMPTE DES TERMES SOURCES LAGRANGIEN : PARTIE EXPLICITE
!     COUPLAGE RETOUR
!===============================================================================

!     Ordre 2 non pris en compte
if (iilagr.eq.2 .and. ltsdyn.eq.1 .and. iphas.eq.ilphas) then

  do iel = 1,ncel

! Termes sources explicte et implicte sur k

    smbrk(iel)  = smbrk(iel) + tslagr(iel,itske)

! Termes sources explicte sur omega : on reprend la constante CE4 directement
!    du k-eps sans justification ... a creuser si necessaire           !

    smbrw(iel)  = smbrw(iel)                                      &
                + ce4 *tslagr(iel,itske) * propce(iel,ipcroo)     &
                /propce(iel,ipcvto)
  enddo

endif


!===============================================================================
! 9. PRISE EN COMPTE DES TERMES DE CONV/DIFF DANS LE SECOND MEMBRE

!      Tableaux de travail              W7, W8, W1, TINSTK, TINSTW, DRTP
!      Les termes sont stockes dans     W5 et W6, puis ajoutes a SMBRK, SMBRW
!      En sortie de l'etape on conserve W2-6,SMBRK,SMBRW,DAM
!===============================================================================

!     Ceci ne sert a rien si IKECOU n'est pas egal a 1

if (ikecou(iphas).eq.1) then

  do iel = 1, ncel
    w5 (iel) = 0.d0
    w6 (iel) = 0.d0
  enddo

! ---> Traitement de k

  ivar   = ikiph

  ipp    = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef )
  iclvaf = iclrtp(ivar,icoeff)
  chaine = nomvar(ipp)

  if( idiff(ivar).ge. 1 ) then

    do iel = 1, ncel
      xxf1 = xf1(iel)
      sigma = xxf1*ckwsk1 + (1.d0-xxf1)*ckwsk2
      w7(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*propce(iel,ipcvst)/sigma
    enddo
    call viscfa                                                   &
    !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w7     ,                                                       &
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

  iccocg = 1
  inc    = 1
  iconvp = iconv (ivar)
  idiffp = idiff (ivar)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  iwarnp = iwarni(ivar)
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = thetav(ivar)

  call bilsc2                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  ,          &
   w5  ,                                                          &
!        --
   w7     , w8     , w1     , tinstk , tinstw , drtp ,            &
   ra     )

  if (iwarni(ivar).ge.2) then
    isqrt = 1
    call prodsc(ncelet,ncel,isqrt,smbrk,smbrk,rnorm)
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif


! ---> Traitement de omega

  ivar   = iomgip

  ipp    = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef )
  iclvaf = iclrtp(ivar,icoeff)
  chaine = nomvar(ipp)

  if( idiff(ivar).ge. 1 ) then
    do iel = 1, ncel
      xxf1 = xf1(iel)
      sigma = xxf1*ckwsw1 + (1.d0-xxf1)*ckwsw2
      w7(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*propce(iel,ipcvst)/sigma
    enddo
    call viscfa                                                   &
    !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w7     ,                                                       &
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

  iccocg = 1
  inc    = 1
  iconvp = iconv (ivar)
  idiffp = idiff (ivar)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  iwarnp = iwarni(ivar)
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = thetav(ivar)

  call bilsc2                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  ,          &
   w6  ,                                                          &
!        --
   w7     , w8     , w1     , tinstk , tinstw , drtp ,            &
   ra     )

  if (iwarni(ivar).ge.2) then
    isqrt = 1
    call prodsc(ncelet,ncel,isqrt,smbrw,smbrw,rnorm)
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif

  do iel = 1,ncel
    smbrk(iel) = smbrk(iel) + w5(iel)
    smbrw(iel) = smbrw(iel) + w6(iel)
  enddo

endif

!===============================================================================
! 10. AJOUT DES TERMES SOURCES DE MASSE EXPLICITES

!       Les parties implicites eventuelles sont conservees dans W7 et W8
!         et utilisees dans la phase d'implicitation cv/diff

!       Les termes sont stockes dans     SMBRK, SMBRW, W7, W8
!       En sortie de l'etape on conserve W2-8,SMBRK,SMBRW,DAM
!===============================================================================

if (ncesmp.gt.0) then

  do iel = 1, ncel
    w7(iel) = 0.d0
    w8(iel) = 0.d0
  enddo

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

!       On incremente SMBRS par -Gamma RTPA et ROVSDT par Gamma (*theta)
  ivar = ikiph
  call catsma                                                     &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
                                 isto2t(iphas) , thetav(ivar) ,   &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipriph) ,    &
   smbrk  , w7     , tinstk )
  ivar = iomgip
  call catsma                                                     &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
                                 isto2t(iphas) , thetav(ivar) ,   &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipriph) ,    &
   smbrw  , w8     , tinstw )

!       Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if(isto2t(iphas).gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta  ) = propce(iel,iptsta  ) + tinstk(iel)
      propce(iel,iptsta+1) = propce(iel,iptsta+1) + tinstw(iel)
    enddo
!       Sinon on le met directement dans SMBR
  else
    do iel = 1, ncel
      smbrk(iel) = smbrk(iel) + tinstk(iel)
      smbrw(iel) = smbrw(iel) + tinstw(iel)
    enddo
  endif

endif

!     Finalisation des termes sources
if(isto2t(iphas).gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + thetp1 * propce(iel,iptsta)
    smbrw(iel) = smbrw(iel) + thetp1 * propce(iel,iptsta+1)
  enddo
endif

!===============================================================================
! 11. INCREMENTS DES TERMES SOURCES DANS LE SECOND MEMBRE

!       Les termes sont stockes dans     SMBRK, SMBRW
!       En sortie de l'etape on conserve W3-8,SMBRK,SMBRW
!===============================================================================

!     Ordre 2 non pris en compte
if(ikecou(iphas).eq.1) then

  do iel = 1, ncel

    rom = propce(iel,ipcrom)

!   RESOLUTION COUPLEE

    romvsd     = 1.d0/(rom*volume(iel))
    smbrk(iel) = smbrk(iel)*romvsd
    smbrw(iel) = smbrw(iel)*romvsd
    divp23     = d2s3*max(divukw(iel),zero)
    produc     = s2kw(iel)-d2s3*divukw(iel)**2+w2(iel)
    xk         = rtpa(iel,ikiph)
    xw         = rtpa(iel,iomgip)
    xxf1       = xf1(iel)
    xgamma     = xxf1*ckwgm1 + (1.d0-xxf1)*ckwgm2
    xbeta      = xxf1*ckwbt1 + (1.d0-xxf1)*ckwbt2

    a11 = 1.d0/dt(iel)                                            &
         - 1.d0/xw*min(produc,zero)+divp23+cmu*xw
    a12 = cmu*xk
    a21 = 0.d0
    a22 = 1.d0/dt(iel)+xgamma*divp23+2.d0*xbeta*xw

    unsdet = 1.d0/(a11*a22 -a12*a21)

    deltk = ( a22*smbrk(iel) -a12*smbrw(iel) )*unsdet
    deltw = (-a21*smbrk(iel) +a11*smbrw(iel) )*unsdet

!     NOUVEAU TERME SOURCE POUR CODITS

    romvsd = rom*volume(iel)/dt(iel)

    smbrk(iel) = romvsd*deltk
    smbrw(iel) = romvsd*deltw

  enddo

endif

!===============================================================================
! 12. TERMES INSTATIONNAIRES

!     Les termes sont stockes dans     TINSTK, TINSTW
!     En sortie de l'etape on conserve SMBRK, SMBRW,  TINSTK, TINSTW
!===============================================================================

! --- PARTIE EXPLICITE

!     on enleve la convection/diffusion au temps n a SMBRK et SMBRW
!     s'ils ont ete calcules
if (ikecou(iphas).eq.1) then
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) - w5(iel)
    smbrw(iel) = smbrw(iel) - w6(iel)
  enddo
endif

! --- RHO/DT et DIV
!     Extrapolation ou non, meme forme par coherence avec bilsc2

do iel = 1, ncel
  rom = propce(iel,ipcrom)
  romvsd = rom*volume(iel)/dt(iel)
  tinstk(iel) = istat(ikiph)*romvsd                               &
               -iconv(ikiph)*w4(iel)*thetav(ikiph)
  tinstw(iel) = istat(iomgip)*romvsd                              &
               -iconv(iomgip)*w4(iel)*thetav(iomgip)
enddo

! --- Source de masse (le theta est deja inclus par catsma)
if (ncesmp.gt.0) then
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) + w7(iel)
    tinstw(iel) = tinstw(iel) + w8(iel)
  enddo
endif

! --- Termes sources utilisateurs
if(isto2t(iphas).gt.0) then
  thetak = thetav(ikiph)
  thetaw = thetav(iomgip)
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) -dam(iel)*thetak
    tinstw(iel) = tinstw(iel) -w3 (iel)*thetaw
  enddo
else
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) + max(-dam(iel),zero)
    tinstw(iel) = tinstw(iel) + max(-w3 (iel),zero)
  enddo
endif

! --- PRISE EN COMPTE DES TERMES LAGRANGIEN : COUPLAGE RETOUR

!     Ordre 2 non pris en compte
if (iilagr.eq.2 .and. ltsdyn.eq.1 .and. iphas.eq.1) then

  do iel = 1,ncel

! Termes sources implicite sur k

    tinstk(iel) = tinstk(iel) + max(-tslagr(iel,itsli),zero)

! Termes sources implicte sur omega

    tinstw(iel) = tinstw(iel)                                     &
          + max( (-ce4*tslagr(iel,itske)/rtpa(iel,ikiph)) , zero)

  enddo

endif

! Si IKECOU=0, on implicite plus fortement k et omega

if(ikecou(iphas).eq.0)then
  do iel=1,ncel
    xw    = rtpa(iel,iomgip)
    xxf1  = xf1(iel)
    xbeta = xxf1*ckwbt1 + (1.d0-xxf1)*ckwbt2
    rom = propce(iel,ipcrom)
      tinstk(iel) = tinstk(iel) +                                 &
         volume(iel)*cmu*rom*xw
      tinstw(iel) = tinstw(iel) +                                 &
         volume(iel)*xbeta*rom*xw
  enddo
endif


!===============================================================================
! 13. RESOLUTION

!       On utilise                      SMBRK, SMBRW,  TINSTK, TINSTW
!===============================================================================

! ---> Traitement de k

ivar = ikiph
iclvar = iclrtp(ivar,icoef )
iclvaf = iclrtp(ivar,icoeff)

ipp    = ipprtp(ivar)

!     "VITESSE" DE DIFFUSION FACETTE

if( idiff(ivar).ge. 1 ) then

  do iel = 1, ncel
    xxf1 = xf1(iel)
    sigma = xxf1*ckwsk1 + (1.d0-xxf1)*ckwsk2
    w1(iel) = propce(iel,ipcvis)                                  &
                        + idifft(ivar)*propce(iel,ipcvst)/sigma
  enddo
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

!     RESOLUTION POUR K

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
!MO      IPP    =
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)
thetap = thetav(ivar)

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   tinstk , smbrk  , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )


! ---> Traitement de omega

ivar = iomgip
iclvar = iclrtp(ivar,icoef )
iclvaf = iclrtp(ivar,icoeff)

ipp    = ipprtp(ivar)


!     "VITESSE" DE DIFFUSION FACETTE

if( idiff(ivar).ge. 1 ) then
  do iel = 1, ncel
    xxf1 = xf1(iel)
    sigma = xxf1*ckwsw1 + (1.d0-xxf1)*ckwsw2
    w1(iel) = propce(iel,ipcvis)                                  &
                        + idifft(ivar)*propce(iel,ipcvst)/sigma
  enddo
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

! RESOLUTION POUR OMEGA

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
!MO      IPP    =
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)
thetap = thetav(ivar)

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   tinstw , smbrw  , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )


!===============================================================================
! 14. CLIPPING
!===============================================================================

!     Calcul des Min/Max avant clipping, pour affichage
do ii = 1, 2
  if(ii.eq.1) then
    ivar = ikiph
  elseif(ii.eq.2) then
    ivar = iomgip
  endif
  ipp  = ipprtp(ivar)

  vrmin =  grand
  vrmax = -grand
  do iel = 1, ncel
    var = rtp(iel,ivar)
    vrmin = min(vrmin,var)
    vrmax = max(vrmax,var)
  enddo
  if (irangp.ge.0) then
    call parmax (vrmax)
    !==========
    call parmin (vrmin)
    !==========
  endif
  varmna(ipp) = vrmin
  varmxa(ipp) = vrmax

enddo

!     On clippe simplement k et omega par valeur absolue
iclipk = 0
iclipw = 0
do iel = 1, ncel
  xk = rtp(iel,ikiph )
  xw = rtp(iel,iomgip)
  if (abs(xk).le.epz2) then
    iclipk = iclipk + 1
    rtp(iel,ikiph) = max(rtp(iel,ikiph),epz2)
  elseif(xk.le.0.d0) then
    iclipk = iclipk + 1
    rtp(iel,ikiph) = -xk
  endif
  if (abs(xw).le.epz2) then
    iclipw = iclipw + 1
    rtp(iel,iomgip) = max(rtp(iel,iomgip),epz2)
  elseif(xw.le.0.d0) then
    iclipw = iclipw + 1
    rtp(iel,iomgip) = -xw
  endif
enddo

if (irangp.ge.0) then
  call parcpt (iclipk)
  !==========
  call parcpt (iclipw)
  !==========
endif

! ---  Stockage nb de clippings pour listing

iclpmn(ipprtp(ikiph )) = iclipk
iclpmn(ipprtp(iomgip)) = iclipw


!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** PHASE ',I4,' RESOLUTION DU K-OMEGA                     ',/,&
'      ----------------------------------                     ',/)
 1100 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)

#else

 1000 format(/,                                                   &
'   ** PHASE ',I4,' SOLVING K-OMEGA'                           ,/,&
'      ----------------------------'                           ,/)
 1100 format(1X,A8,' : EXPLICIT BALANCE = ',E14.5)

#endif

!----
! FIN
!----

return

end subroutine
