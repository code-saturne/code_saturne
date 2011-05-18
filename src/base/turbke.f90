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

subroutine turbke &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   tslagr , coefa  , coefb  , ckupdc , smacel ,                   &
   viscf  , viscb  , prdv2f ,                                     &
   dam    , xam    ,                                              &
   drtp   , smbrk  , smbre  , tinstk , tinste ,                   &
   divu   , w1     , w2     , w3     , w4     ,                   &
   w5     , w6     , w7     , w8     , w9     ,                   &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS K-EPS 1 PHASE INCOMPRESSIBLE OU
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
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
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
! prdv2f(ncelet    ! tr ! <-- ! tableau de stockage du terme de                !
!                  !    !     ! prod de turbulence pour le v2f                 !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr.(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! tinst.(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! divu(ncelet      ! tr ! --- ! tableau de travail                             !
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
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ncepdp , ncesmp

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
double precision prdv2f(ncelet)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet),   smbrk(ncelet),  smbre(ncelet)
double precision tinstk(ncelet), tinste(ncelet), divu(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra, ifinia
integer          iel   , ifac  , init  , inc   , iccocg, ivar
integer          iivar , iiun
integer          iclip , isqrt
integer          nswrgp, imligp
integer          ipriph, iuiph , iviph , iwiph
integer          ikiph , ieiph , iphiph
integer          icliup, iclivp, icliwp
integer          iclvar, iclvaf
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          ipcrom, ipbrom, ipcvst, ipcvis, iflmas, iflmab
integer          iwarnp, ipp
integer          iptsta
integer          ipcroo, ipbroo, ipcvto, ipcvlo
double precision rnorm , d2s3, divp23
double precision deltk , delte, a11, a12, a22, a21
double precision gravke, epssuk, unsdet, romvsd
double precision prdtur, xk, xeps, xphi, xnu, ttke, ttmin, tt
double precision visct , rom   , ceps1 , ctsqnu
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision thetp1, thetak, thetae, thets, thetap
double precision tuexpk, tuexpe
double precision cmueta, sqrcmu, xs

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

ipriph = ipr
iuiph  = iu
iviph  = iv
iwiph  = iw
ikiph  = ik
ieiph  = iep
iphiph = iphi

icliup = iclrtp(iuiph,icoef)
iclivp = iclrtp(iviph,icoef)
icliwp = iclrtp(iwiph,icoef)

ipcrom = ipproc(irom  )
ipcvst = ipproc(ivisct)
ipcvis = ipproc(iviscl)
iflmas = ipprof(ifluma(iuiph))
iflmab = ipprob(ifluma(iuiph))
ipbrom = ipprob(irom  )

thets  = thetst

ipcroo = ipcrom
ipbroo = ipbrom
ipcvto = ipcvst
ipcvlo = ipcvis
if(isto2t.gt.0) then
  if (iroext.gt.0) then
    ipcroo = ipproc(iroma)
    ipbroo = ipprob(iroma)
  endif
  if(iviext.gt.0) then
    ipcvto = ipproc(ivista)
    ipcvlo = ipproc(ivisla)
  endif
endif

if(isto2t.gt.0) then
  iptsta = ipproc(itstua)
else
  iptsta = 0
endif

if(iwarni(ikiph).ge.1) then
  if (iturb.eq.20) then
    write(nfecra,1000)
  else if (iturb.eq.21) then
    write(nfecra,1001)
  else
    write(nfecra,1002)
  endif
endif

!     Pour la prod lineaire on a besoin de SQRT(Cmu)
sqrcmu = sqrt(cmu)

!===============================================================================
! 2. CALCUL DE SIJ SIJ ET DE DIVU

!      Tableaux de travail              SMBRK,TINSTE,W1,W2,W3,W4,W5,W6
!      SijSij est stocke dans           TINSTK
!      DivU est stocke dans             DIVU
!      En sortie de l'etape on conserve TINSTK, DIVU
!===============================================================================

iccocg = 1
inc = 1

! ON S'APPUIE SUR LES TABLEAUX DE TRAVAIL W4,W5,W6
!     ON UTILISE EGALEMENT SMBRK  ET TINSTE COMME
!                           TABLEAUX DE TRAVAIL TEMPORAIRES
!     W1,W2,W3 SONT UTILISES DANS GRDCEL
!     TINSTK RECOIT LA PRODUCTION
!        TINSTK = {(2 (S11)**2 + 2 (S22)**2 +2 (S33)**2 )
!           +         ((2 S12)**2 + (2 S13)**2 +(2 S23)**2 )}
!               = 2 Sij.Sij


! SMBRK  = DUDX ,W4 = DUDY ,W5 = DUDZ

nswrgp = nswrgr(iuiph)
imligp = imligr(iuiph)
iwarnp = iwarni(ikiph)
epsrgp = epsrgr(iuiph)
climgp = climgr(iuiph)
extrap = extrag(iuiph)

call grdcel                                                       &
!==========
 ( iuiph  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rtpa(1,iuiph)   , coefa(1,icliup) , coefb(1,icliup) ,          &
   smbrk  , w4     , w5     ,                                     &
!        ------   ------   ------
   w1     , w2     , w3     ,                                     &
   ra     )


! TINSTK = (S11)**2
! DIVU = DUDX

do iel = 1, ncel
  tinstk(iel)   = smbrk(iel)**2
  divu  (iel)   = smbrk(iel)
enddo

!               ,W4     = DUDY ,W5     = DUDZ
! TINSTE = DVDX ,SMBRK  = DVDY ,W6     = DVDZ

nswrgp = nswrgr(iviph)
imligp = imligr(iviph)
iwarnp = iwarni(ikiph)
epsrgp = epsrgr(iviph)
climgp = climgr(iviph)
extrap = extrag(iviph)

call grdcel                                                       &
!==========
 ( iviph  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rtpa(1,iviph)   , coefa(1,iclivp) , coefb(1,iclivp) ,          &
   tinste , smbrk  , w6     ,                                     &
!        ------   ------   ------
   w1     , w2     , w3     ,                                     &
   ra     )


! TINSTK = 2 (S11)**2 + 2 (S22)**2
!      + (2 S12)**2
! DIVU = DUDX + DVDY

do iel = 1, ncel
  tinstk (iel)   = 2.d0*(tinstk(iel) + smbrk(iel)**2 )            &
           +  (tinste(iel)+w4(iel))**2
  divu  (iel)   = divu(iel) + smbrk(iel)
enddo

!               ,              ,W5     = DUDZ
!               ,              ,W6     = DVDZ
! W4     = DWDX ,TINSTE = DWDY ,SMBRK  = DWDZ

nswrgp = nswrgr(iwiph)
imligp = imligr(iwiph)
iwarnp = iwarni(ikiph)
epsrgp = epsrgr(iwiph)
climgp = climgr(iwiph)
extrap = extrag(iwiph)

call grdcel                                                       &
!==========
 ( iwiph  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rtpa(1,iwiph)   , coefa(1,icliwp) , coefb(1,icliwp) ,          &
   w4     , tinste , smbrk  ,                                     &
!        ------   ------   ------
   w1     , w2     , w3     ,                                     &
   ra     )

! TINSTK = PRODUCTION = (
!      2 (S11)**2 + 2 (S22)**2 + 2 (S33)**2
!      + (2 S12)**2 + (2 S13)**2 + (2 S23)**2 )
! DIVU = DUDX + DVDY + DWDZ

do iel = 1, ncel
  tinstk (iel)   = tinstk(iel) + 2.d0*smbrk(iel)**2               &
           +  (w4(iel)+w5(iel))**2                                &
           +  (tinste(iel)+w6(iel))**2
  divu  (iel)   = divu(iel) + smbrk(iel)
enddo

! On libere   SMBRK,TINSTE,W1,W2,W3,W4,W5,W6

!===============================================================================
! 3. PRISE EN COMPTE DES TERMES SOURCES UTILISATEURS

!      On passe 2 Sij.Sij = TINSTK et la divergence DIVU
!      Tableaux de travail                        W1, W2, W3, W4, W5, W6
!                                VISCF VISCB XAM DRTP SMBRK SMBRE TINSTE
!      La partie a expliciter est stockee dans    W7, W8
!      La partie a impliciter est stockee dans    DAM, W9
!      En sortie de l'etape on conserve           TINSTK, DIVU,
!                                                 W7 , W8, DAM, W9
!===============================================================================
! viscf viscb xam drtp smbrk smbre tinste
do iel = 1, ncel
  dam(iel) = 0.d0
  w9 (iel) = 0.d0
  w7 (iel) = 0.d0
  w8 (iel) = 0.d0
enddo

call ustske                                                       &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , tinstk , divu   ,          &
   w7     , w8     , dam    , w9     ,                            &
!        ------   ------   ------   ------
   viscf  , viscb  , xam    ,                                     &
   ra     )

! On libere W1, W2, W3, W4, W5, W6
!           VISCF, VISCB, XAM, DRTP, SMBRK, SMBRE, TINSTE
!===============================================================================
! 4. AJOUT DE - 2/3 DIV(U) * DIV(U)

!      En sortie de l'etape on conserve TINSTK, DIVU,
!                                       W7 , W8, DAM, W9
!===============================================================================

! Dans le cas de la production lineaire, seul le terme en divu est
! multiplie par VISCT. Pour les autres modeles, la multiplication par
! VISCT sera faite ulterieurement.
! A ce stade, TINSTK contient S**2
d2s3 = 2.d0/3.d0
if (iturb.eq.21) then
  do iel = 1, ncel
    rom   = propce(iel,ipcroo)
    visct = propce(iel,ipcvto)
    xs = sqrt(tinstk(iel))
    cmueta = cmu*rtpa(iel,ikiph)/rtpa(iel,ieiph)*xs
    cmueta = min(cmueta,sqrcmu)
    tinstk(iel) = rom*cmueta*xs*rtpa(iel,ikiph)                   &
         - d2s3*visct*divu(iel)*divu(iel)
  enddo
else
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) - d2s3*divu(iel)*divu(iel)
  enddo
endif

!===============================================================================
! 5. CALCUL DU TERME DE GRAVITE

!      Les s.m. recoivent production et termes de gravite
!      Tableaux de travail              W1,W2,W3,W4,W5,W6,VISCB
!      Les s.m. sont stockes dans       TINSTK TINSTE
!      En sortie de l'etape on conserve TINSTK, TINSTE,
!                                       DIVU,
!                                       W7 , W8, DAM, W9
!===============================================================================

if (igrake.eq.1 .and. ippmod(iatmos).ge.1) then

    !  Calcul du terme de gravite pour la version atmospherique

    call atprke                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nscal  ,                                                       &
   ipcvto,                                                        &
   ia     ,                                                       &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3    ,                                      &
   w4     , w5     , w6    ,                                      &
   tinstk , tinste ,                                              &
   ra     )

else if (igrake.eq.1) then

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

  call grdcel                                                     &
  !==========
 ( iivar  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   propce(1,ipcroo), propfb(1,ipbroo), viscb  ,                   &
   w4     , w5     , w6     ,                                     &
!        ------   ------   ------
   w1     , w2     , w3     ,                                     &
   ra     )


!      Production et terme de gravite
!        TINSTK=P+G et TINSTE=P+(1-CE3)*G

  if(iscalt.gt.0.and.nscal.ge.iscalt) then
    prdtur = sigmas(iscalt)
  else
    prdtur = 1.d0
  endif

!     En production lineaire, on multiplie tout de suite le terme
!     de gravite par VISCT, car le reste est deja multiplie.
!     Dans les autres cas, la multiplication est faite plus tard.
  if (iturb.eq.21) then
    do iel = 1, ncel
      gravke = -(w4(iel)*gx+w5(iel)*gy+w6(iel)*gz)/               &
           (propce(iel,ipcroo)*prdtur)
      tinste(iel)=tinstk(iel)+propce(iel,ipcvto)*max(gravke,zero)
      tinstk(iel)=tinstk(iel)+propce(iel,ipcvto)*gravke
    enddo
  else
    do iel = 1, ncel
      gravke = -(w4(iel)*gx+w5(iel)*gy+w6(iel)*gz)/               &
           (propce(iel,ipcroo)*prdtur)
      tinste(iel) = tinstk(iel) + max( gravke,zero )
      tinstk(iel) = tinstk(iel) + gravke
    enddo
  endif

else


! --- Production sans termes de gravite
!       TINSTK=TINSTE=P

  do iel = 1, ncel
    tinste(iel) = tinstk(iel)
  enddo

endif

!     En V2F, on stocke TINSTK dans PRDV2F qui sera complete plus loin pour
!     contenir le terme de production complet
if (iturb.eq.50) then
  do iel = 1, ncel
    prdv2f(iel) = tinstk(iel)
  enddo
endif

! On libere W1,W2,W3,W4,W5,W6,VISCB

!===============================================================================
! 6. TERME D'ACCUMULATION DE MASSE -(dRO/dt)*Volume

!      Le terme est stocke dans         W1
!      En sortie de l'etape on conserve W1, TINSTK, TINSTE, DIVU,
!                                       W7 , W8, DAM, W9
!===============================================================================

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w1)

!===============================================================================
! 7. ON FINALISE LE CALCUL DES TERMES SOURCES

!      Les termes sont stockes dans     SMBRK, SMBRE
!      En sortie de l'etape on conserve W1, TINSTK, TINSTE, DIVU,
!                                       SMBRK, SMBRE
!                                       W7 , W8, DAM, W9
!===============================================================================

! SMBRE = CEPS1 EPSILON/K (PROD + G ) - RO0 VOLUME EPSILON EPSILON/K
! SMBRK = PROD + G                  - RO0 VOLUME EPSILON

!     Si on extrapole les termes sources et rho  , il faut ici rho^n
!                                        et visct, il faut ici visct^n

if (iturb.eq.20) then

  do iel = 1, ncel

    visct = propce(iel,ipcvto)
    rom   = propce(iel,ipcroo)

    smbrk(iel) = volume(iel)*(                                    &
         visct*tinstk(iel)                                        &
         -d2s3*rom*rtpa(iel,ikiph)*divu(iel)                      &
         -rom*rtpa(iel,ieiph) )

    smbre(iel) = volume(iel)*rtpa(iel,ieiph)/rtpa(iel,ikiph)*(    &
         ce1*( visct*tinste(iel)                                  &
         -d2s3*rom*rtpa(iel,ikiph)*divu(iel) )                    &
         -ce2*rom*rtpa(iel,ieiph) )

  enddo

else if (iturb.eq.21) then

  do iel = 1, ncel

    rom   = propce(iel,ipcroo)

    smbrk(iel) = volume(iel)*(                                    &
         tinstk(iel)                                              &
         -d2s3*rom*rtpa(iel,ikiph)*divu(iel)                      &
         -rom*rtpa(iel,ieiph) )

    smbre(iel) = volume(iel)*rtpa(iel,ieiph)/rtpa(iel,ikiph)*(    &
         ce1*(tinste(iel)                                         &
         -d2s3*rom*rtpa(iel,ikiph)*divu(iel) )                    &
         -ce2*rom*rtpa(iel,ieiph) )

  enddo

else if (iturb.eq.50) then

  do iel = 1, ncel

    visct = propce(iel,ipcvto)
    rom   = propce(iel,ipcroo)
    xeps = rtpa(iel,ieiph )
    xk   = rtpa(iel,ikiph )
    xphi = rtpa(iel,iphiph)
    xphi = max(xphi,epzero)
    xnu  = propce(iel,ipcvlo)/rom
    ceps1= 1.4d0*(1.d0+cv2fa1*sqrt(1.d0/xphi))
    ttke = xk / xeps
    ttmin = cv2fct*sqrt(xnu/xeps)
    tt = max(ttke,ttmin)

    smbrk(iel) = volume(iel)*(                                    &
         visct*tinstk(iel)                                        &
         -d2s3*rom*rtpa(iel,ikiph)*divu(iel)                      &
         -rom*rtpa(iel,ieiph) )

    smbre(iel) = volume(iel)/tt*(                                 &
         ceps1*( visct*tinste(iel)                                &
         -d2s3*rom*rtpa(iel,ikiph)*divu(iel) )                    &
         -cv2fe2*rom*rtpa(iel,ieiph) )

!     On stocke la partie en Pk dans PRDV2F pour etre reutilise dans RESV2F
    prdv2f(iel) = visct*prdv2f(iel)                               &
         -d2s3*rom*rtpa(iel,ikiph)*divu(iel)

  enddo

endif


!===============================================================================
! 8. PRISE EN COMPTE DES TERMES SOURCES UTILISATEURS
!                        ET ACCUMULATION DE MASSE    : PARTIE EXPLICITE
!      On utilise                       W1,  W7, W8, DAM, W9
!      Les termes sont stockes dans     SMBRK, SMBRE
!      En sortie de l'etape on conserve W1, TINSTK, TINSTE, DIVU,
!                                       SMBRK, SMBRE
!                                       DAM, W9

!    Remarque : l'extrapolation telle qu'elle est ecrite n'a pas grand
!               sens si IKECOU=1
!===============================================================================

!     Si on extrapole les T.S.
if(isto2t.gt.0) then

  do iel = 1, ncel

!       Sauvegarde pour echange
    tuexpk = propce(iel,iptsta)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta) = smbrk(iel) + w7(iel)
!       Termes dependant de la variable resolue et theta PROPCE
    smbrk(iel) = iconv(ikiph)*w1(iel)*rtpa(iel,ikiph)             &
         - thets*tuexpk
!       On suppose -DAM > 0 : on implicite
!         le terme utilisateur dependant de la variable resolue
    smbrk(iel) = dam(iel)*rtpa(iel,ikiph) + smbrk(iel)

!       Sauvegarde pour echange
    tuexpe = propce(iel,iptsta+1)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta+1) = smbre(iel) + w8(iel)
!       Termes dependant de la variable resolue et theta PROPCE
    smbre(iel) = iconv(ieiph)*w1(iel)*rtpa(iel,ieiph)             &
         - thets*tuexpe
!       On suppose -W9 > 0 : on implicite
!         le terme utilisateur dependant de la variable resolue
    smbre(iel) =  w9(iel)*rtpa(iel,ieiph) + smbre(iel)

  enddo

!     Si on n'extrapole pas les T.S.
else
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + dam(iel)*rtpa(iel,ikiph) + w7(iel)  &
         +iconv(ikiph)*w1(iel)*rtpa(iel,ikiph)
    smbre(iel) = smbre(iel) + w9 (iel)*rtpa(iel,ieiph) + w8(iel)  &
         +iconv(ieiph)*w1(iel)*rtpa(iel,ieiph)
  enddo
endif

!===============================================================================
! 8.1 PRISE EN COMPTE DES TERMES SOURCES LAGRANGIEN : PARTIE EXPLICITE
!     COUPLAGE RETOUR
!===============================================================================

!     Ordre 2 non pris en compte
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  do iel = 1,ncel

! Termes sources explicte et implicte sur k

    smbrk(iel)  = smbrk(iel) + tslagr(iel,itske)

! Termes sources explicte sur Eps

    smbre(iel)  = smbre(iel)                                      &
                + ce4 *tslagr(iel,itske) *rtpa(iel,ieiph)         &
                                         /rtpa(iel,ikiph)

  enddo

endif

! ON LIBERE W7, W8

!===============================================================================
! 9. PRISE EN COMPTE DES TERMES DE CONV/DIFF DANS LE SECOND MEMBRE

!      Tableaux de travail              W2, W3, W4, W5, W6, DRTP
!      Les termes sont stockes dans     W7 et W8, puis ajoutes a SMBRK, SMBRE
!      En sortie de l'etape on conserve W1, TINSTK, TINSTE, DIVU,
!                                       SMBRK, SMBRE
!                                       DAM, W7, W8, W9
!===============================================================================

!     Ceci ne sert a rien si IKECOU n'est pas egal a 1

if (ikecou.eq.1) then

  do iel = 1, ncel
    w7 (iel) = 0.d0
    w8 (iel) = 0.d0
  enddo

! ---> Traitement de k

  ivar   = ikiph

  ipp    = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef )
  iclvaf = iclrtp(ivar,icoeff)
  chaine = nomvar(ipp)

  if( idiff(ivar).ge. 1 ) then

    do iel = 1, ncel
      w4(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*propce(iel,ipcvst)/sigmak
    enddo
    call viscfa                                                   &
    !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w4     ,                                                       &
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
   nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  ,          &
   w7  ,                                                          &
!        --
   w2     , w3     , w4     , w5     , w6     , drtp ,            &
   ra     )

  if (iwarni(ivar).ge.2) then
    isqrt = 1
    call prodsc(ncelet,ncel,isqrt,smbrk,smbrk,rnorm)
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif


! ---> Traitement de epsilon

  ivar   = ieiph

  ipp    = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef )
  iclvaf = iclrtp(ivar,icoeff)
  chaine = nomvar(ipp)

  if( idiff(ivar).ge. 1 ) then
    do iel = 1, ncel
      w4(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*propce(iel,ipcvst)/sigmae
    enddo
    call viscfa                                                   &
    !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w4     ,                                                       &
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
   nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  ,          &
   w8  ,                                                          &
!        --
   w2     , w3     , w4     , w5     , w6     , drtp ,            &
   ra     )

  if (iwarni(ivar).ge.2) then
    isqrt = 1
    call prodsc(ncelet,ncel,isqrt,smbre,smbre,rnorm)
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif

  do iel = 1,ncel
    smbrk(iel) = smbrk(iel) + w7(iel)
    smbre(iel) = smbre(iel) + w8(iel)
  enddo

endif

!===============================================================================
! 10. AJOUT DES TERMES SOURCES DE MASSE EXPLICITES

!       Les parties implicites eventuelles sont conservees dans W2 et W3
!         et utilisees dans la phase d'implicitation cv/diff

!       Les termes sont stockes dans     SMBRK, SMBRE, W2, W3
!       En sortie de l'etape on conserve W1, TINSTK, TINSTE, DIVU,
!                                        SMBRK, SMBRE
!                                        DAM, W9, W2, W3
!===============================================================================

if (ncesmp.gt.0) then

  do iel = 1, ncel
    w2(iel) = 0.d0
    w3(iel) = 0.d0
  enddo

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

!       On incremente SMBRS par -Gamma RTPA et ROVSDT par Gamma (*theta)
  ivar = ikiph
  call catsma                                                     &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
                                 isto2t , thetav(ivar) ,   &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipriph) ,    &
   smbrk  , w2     , w4 )
  ivar = ieiph
  call catsma                                                     &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
                                 isto2t , thetav(ivar) ,   &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipriph) ,    &
   smbre  , w3     , w5 )

!       Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if(isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta  ) = propce(iel,iptsta  ) + w4(iel)
      propce(iel,iptsta+1) = propce(iel,iptsta+1) + w5(iel)
    enddo
!       Sinon on le met directement dans SMBR
  else
    do iel = 1, ncel
      smbrk(iel) = smbrk(iel) + w4(iel)
      smbre(iel) = smbre(iel) + w5(iel)
    enddo
  endif

endif

!     ON LIBERE                       W4, W5, W6, TINSTE

!     Finalisation des termes sources
if(isto2t.gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + thetp1 * propce(iel,iptsta)
    smbre(iel) = smbre(iel) + thetp1 * propce(iel,iptsta+1)
  enddo
endif

!===============================================================================
! 11. INCREMENTS DES TERMES SOURCES DANS LE SECOND MEMBRE

!       On utilise                       TINSTK, TINSTE, DIVU
!       Les termes sont stockes dans     SMBRK, SMBRE
!       En sortie de l'etape on conserve W1, SMBRK, SMBRE,
!                                        DAM, W9, W2, W3, W7, W8
!===============================================================================

!     Ordre 2 non pris en compte
if(ikecou.eq.1) then

  if (iturb.eq.20) then

    do iel = 1, ncel

      rom = propce(iel,ipcrom)

!   RESOLUTION COUPLEE

      romvsd=1.d0/(rom*volume(iel))
      smbrk(iel)=smbrk(iel)*romvsd
      smbre(iel)=smbre(iel)*romvsd
      divp23= d2s3*max(divu(iel),zero)

      epssuk = rtpa(iel,ieiph)/rtpa(iel,ikiph)

      a11 = 1.d0/dt(iel)                                          &
           -2.d0*rtpa(iel,ikiph)/rtpa(iel,ieiph)                  &
           *cmu*min(tinstk(iel),zero)+divp23
      a12 = 1.d0
      a21 = -ce1*cmu*tinste(iel)-ce2*epssuk*epssuk
      a22 = 1.d0/dt(iel)+ce1*divp23                               &
           +2.d0*ce2*epssuk

      unsdet = 1.d0/(a11*a22 -a12*a21)

      deltk = ( a22*smbrk(iel) -a12*smbre(iel) )*unsdet
      delte = (-a21*smbrk(iel) +a11*smbre(iel) )*unsdet

!     NOUVEAU TERME SOURCE POUR CODITS

      romvsd = rom*volume(iel)/dt(iel)

      smbrk(iel) = romvsd*deltk
      smbre(iel) = romvsd*delte

    enddo

!     Dans verini on bloque la combinaison ITURB=21/IKECOU=1
  else if (iturb.eq.21) then

    WRITE(NFECRA,*)'IKECOU=1 NON VALIDE EN K-EPS PROD LIN'
    call csexit (1)
!  Section non totalement validee (a priori ca marche, mais pas trop stable) :
!  en fait le v2f est meilleur avec IKECOU=0, on bloque donc la combinaison
!  ITURB=50/IKECOU=1 au niveau de verini. Ces lignes sont donc inaccessibles.
!  On les laisse au cas ou .....
  else if (iturb.eq.50) then

    do iel = 1, ncel

      rom = propce(iel,ipcrom)

!   RESOLUTION COUPLEE

      romvsd=1.d0/(rom*volume(iel))
      smbrk(iel)=smbrk(iel)*romvsd
      smbre(iel)=smbre(iel)*romvsd
      divp23= d2s3*max(divu(iel),zero)

      xeps = rtpa(iel,ieiph )
      xk   = rtpa(iel,ikiph )
      xphi = rtpa(iel,iphiph)
      xphi = max(xphi,epzero)
      xnu  = propce(iel,ipcvis)/propce(iel,ipcrom)
      ctsqnu= cv2fct*sqrt(xnu)
      ceps1= 1.4d0*(1.d0+cv2fa1*sqrt(1.d0/xphi))
      epssuk = xeps/xk
      ttke = xk / xeps
      ttmin = cv2fct*sqrt(xnu/xeps)

      if(ttke.gt.ttmin) then
        a11 = 1.d0/dt(iel)                                        &
             -2.d0*xk/xeps*xphi                                   &
             *cv2fmu*min(tinstk(iel),zero)+divp23
!     Pour A12 on fait comme en k-eps standard pour l'instant,
!     on ne prend pas le terme en P+G ... est-ce judicieux ?
        a12 = 1.d0
        a21 = -ceps1*cv2fmu*xphi*tinste(iel)-cv2fe2*epssuk*epssuk
        a22 = 1.d0/dt(iel)+ceps1*divp23                           &
             +2.d0*cv2fe2*epssuk
      else
        a11 = 1.d0/dt(iel)                                        &
             -cv2fmu*xphi*ctsqnu*min(tinstk(iel),zero)/sqrt(xeps) &
             +divp23
!     Pour A12 on fait comme en k-eps standard pour l'instant,
!     on ne prend pas le terme en P+G ... est-ce judicieux ?
        a12 = 1.d0
!     Le terme en DIVP23 dans A21 n'est pas forcement judicieux
!     (a-t-on besoin du MAX ?)
        a21 = -ceps1*cv2fmu*xphi*tinste(iel)                      &
             +ceps1*sqrt(xeps)/ctsqnu*divp23
        a22 = 1.d0/dt(iel)+1.d0/2.d0*ceps1*divp23*xk              &
             /ctsqnu/sqrt(xeps)                                   &
             +3.d0/2.d0*cv2fe2/ctsqnu*sqrt(xeps)
      endif

      unsdet = 1.d0/(a11*a22 -a12*a21)

      deltk = ( a22*smbrk(iel) -a12*smbre(iel) )*unsdet
      delte = (-a21*smbrk(iel) +a11*smbre(iel) )*unsdet

!     NOUVEAU TERME SOURCE POUR CODITS

      romvsd = rom*volume(iel)/dt(iel)

      smbrk(iel) = romvsd*deltk
      smbre(iel) = romvsd*delte

    enddo

  endif

endif

!     ON LIBERE                       TINSTK, TINSTE, DIVU

!===============================================================================
! 12. TERMES INSTATIONNAIRES

!     On utilise                       W1, W2, W3, W7, W8
!                                      DAM, W9
!     Les termes sont stockes dans     TINSTK, TINSTE
!     En sortie de l'etape on conserve SMBRK, SMBRE,  TINSTK, TINSTE
!===============================================================================

! --- PARTIE EXPLICITE

!     on enleve la convection/diffusion au temps n a SMBRK et SMBRE
!     si on les avait calcules
if (ikecou.eq.1) then
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) - w7(iel)
    smbre(iel) = smbre(iel) - w8(iel)
  enddo
endif

! --- RHO/DT et DIV
!     Extrapolation ou non, meme forme par coherence avec bilsc2

do iel = 1, ncel
  rom = propce(iel,ipcrom)
  romvsd = rom*volume(iel)/dt(iel)
  tinstk(iel) = istat(ikiph)*romvsd                               &
               -iconv(ikiph)*w1(iel)*thetav(ikiph)
  tinste(iel) = istat(ieiph)*romvsd                               &
               -iconv(ieiph)*w1(iel)*thetav(ieiph)
enddo

! --- Source de masse (le theta est deja inclus par catsma)
if (ncesmp.gt.0) then
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) + w2(iel)
    tinste(iel) = tinste(iel) + w3(iel)
  enddo
endif

! --- Termes sources utilisateurs
if(isto2t.gt.0) then
  thetak = thetav(ikiph)
  thetae = thetav(ieiph)
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) -dam(iel)*thetak
    tinste(iel) = tinste(iel) -w9 (iel)*thetae
  enddo
else
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) + max(-dam(iel),zero)
    tinste(iel) = tinste(iel) + max(-w9 (iel),zero)
  enddo
endif

! --- PRISE EN COMPTE DES TERMES LAGRANGIEN : COUPLAGE RETOUR

!     Ordre 2 non pris en compte
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  do iel = 1,ncel

! Termes sources implicite sur k

    tinstk(iel) = tinstk(iel) + max(-tslagr(iel,itsli),zero)

! Termes sources implicte sur Eps

    tinste(iel) = tinste(iel)                                     &
          + max( (-ce4*tslagr(iel,itske)/rtpa(iel,ikiph)) , zero)

  enddo

endif

! Si IKECOU=0, on implicite plus fortement k et eps

if(ikecou.eq.0)then
  if(itytur.eq.2)then
    do iel=1,ncel
      xeps = rtpa(iel,ieiph )
      xk   = rtpa(iel,ikiph )
      rom = propce(iel,ipcrom)
      ttke = xk / xeps
      if(xk.gt.1.d-12) then
        tinstk(iel) = tinstk(iel) +                               &
             rom*volume(iel)/ttke
      endif
      tinste(iel) = tinste(iel) +                                 &
           ce2*rom*volume(iel)/ttke
    enddo
  else if(iturb.eq.50)then
    do iel=1,ncel
      xeps = rtpa(iel,ieiph )
      xk   = rtpa(iel,ikiph )
      rom = propce(iel,ipcrom)
      xnu  = propce(iel,ipcvis)/rom
      ttke = xk / xeps
      ttmin = cv2fct*sqrt(xnu/xeps)
      tt = max(ttke,ttmin)
      if(xk.gt.1.d-12) then
        tinstk(iel) = tinstk(iel) +                               &
             rom*volume(iel)/ttke
      endif
      tinste(iel) = tinste(iel) +                                 &
           cv2fe2*rom*volume(iel)/tt
    enddo

  endif
endif

! ON LIBERE W1, W2, W3, DAM, W9

!===============================================================================
! 13. RESOLUTION

!       On utilise                      SMBRK, SMBRE,  TINSTK, TINSTE
!       Tableaux de travail             W1, W2, W3, W4, W5, W6
!===============================================================================

! ---> Traitement de k

ivar = ikiph
iclvar = iclrtp(ivar,icoef )
iclvaf = iclrtp(ivar,icoeff)

ipp    = ipprtp(ivar)

!     "VITESSE" DE DIFFUSION FACETTE

if( idiff(ivar).ge. 1 ) then

  do iel = 1, ncel
    w1(iel) = propce(iel,ipcvis)                                  &
                        + idifft(ivar)*propce(iel,ipcvst)/sigmak
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
   nvar   , nscal  ,                                              &
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


! ---> Traitement de epsilon

ivar = ieiph
iclvar = iclrtp(ivar,icoef )
iclvaf = iclrtp(ivar,icoeff)

ipp    = ipprtp(ivar)


!     "VITESSE" DE DIFFUSION FACETTE

if( idiff(ivar).ge. 1 ) then
  do iel = 1, ncel
    w1(iel) = propce(iel,ipcvis)                                  &
                        + idifft(ivar)*propce(iel,ipcvst)/sigmae
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

! RESOLUTION POUR EPSILON

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
   nvar   , nscal  ,                                              &
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
   tinste , smbre  , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )


!===============================================================================
! 14. CLIPPING
!===============================================================================

iclip = 1
iwarnp = iwarni(ikiph)
call clipke                                                       &
!==========
 ( ncelet , ncel   , nvar   ,                                     &
   iclip  , iwarnp ,                                              &
   propce , rtp    )


!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION DU K-EPSILON                   ',/,&
'      -----------------------                   ',/)
 1001 format(/,                                                   &
'   ** RESOLUTION DU K-EPSILON A PROD LINEAIRE   ',/,&
'      ---------------------------------------   ',/)
 1002 format(/,                                                   &
'   ** RESOLUTION DU V2F (K ET EPSILON)          ',/,&
'      --------------------------------          ',/)
 1100 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)

#else

 1000 format(/,                                                   &
'   ** SOLVING K-EPSILON'                         ,/,&
'      -----------------'                         ,/)
 1001 format(/,                                                   &
'   ** SOLVING K-EPSILON WITH LINEAR PROD'        ,/,&
'      ----------------------------------'        ,/)
 1002 format(/,                                                   &
'   ** SOLVING V2F (K AND EPSILON)'               ,/,&
'      ---------------------------'               ,/)
 1100 format(1X,A8,' : EXPLICIT BALANCE = ',E14.5)
#endif

!----
! FIN
!----

return

end subroutine
