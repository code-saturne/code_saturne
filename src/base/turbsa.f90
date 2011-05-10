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

subroutine turbsa &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   tslagr , coefa  , coefb  , ckupdc , smacel ,                   &
   itypfb ,                                                       &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbrsa          , tinssa          ,                   &
   divu   , w1     , w2     , w3     , w4     ,                   &
   w5     , w6     , w7     , w8     , w9     ,                   &
   ra     )

!===============================================================================
! Purpose:
! --------

! Solving op the equation of nusa, which is the scalar quantity defined by
! the Spalart-Allmaras model, for 1 phase for 1 time-step.

!-------------------------------------------------------------------------------
! Arguments
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
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr.(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! tinst.(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! divu(ncelet      ! tr ! --- ! tableau de travail                             !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
!We have to know if there is any rough wall
use parall
use pointe

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)
integer          itypfb(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet),   smbrsa(ncelet),tinssa(ncelet)
double precision divu(ncelet)
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
integer          nswrgp, imligp, iphydp
integer          ipriph, iuiph , iviph , iwiph
integer          inuiph
integer          icliup, iclivp, icliwp
integer          iclvar, iclvaf
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          ipcrom, ipbrom, ipcvst, ipcvis, iflmas, iflmab
integer          iwarnp, ipp
integer          iptsta
integer          ipcroo, ipbroo, ipcvto, ipcvlo
integer          maxelt, ils
integer          ipatrg

double precision romvsd
double precision visct , rom
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision thets, thetv, thetp1, thetap
double precision tuexpn
double precision cofbnu
double precision chi  , chi3, taussa, nusa, distbf, fw, fv1, fv2
double precision gsa , rsa , dsigma, cv13
double precision surfn, nu0, dsa0, hssa

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
inuiph = inusa

icliup = iclrtp(iuiph,icoef)
iclivp = iclrtp(iviph,icoef)
icliwp = iclrtp(iwiph,icoef)

ipcrom = ipproc(irom  )
ipcvst = ipproc(ivisct)
ipcvis = ipproc(iviscl)
iflmas = ipprof(ifluma(iuiph))
iflmab = ipprob(ifluma(iuiph))
ipbrom = ipprob(irom  )

! S pour source, V pour variable
!terme source grandeur turbulente
thets  = thetst

ivar   = inuiph
thetv  = thetav(ivar)

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

! extrapolation des TS?
if(isto2t.gt.0) then
  iptsta = ipproc(itstua)
else
  iptsta = 0
endif

! Calculation of some constants
dsigma = 1.d0 / csasig
cv13 = csav1**3

!===============================================================================
! 2. CALCUL DE OmegaIJ OmegaIJ ET DE DIVU ET NORME DE GRAD NUSA

!      Tableaux de travail              SMBRSA,TINSSA,W1,W2,W3,W4,W5,W6
!      SijSij est stocke dans           TINSSA
!      DivU est stocke dans             DIVU
!      Grad nu est stocke dans          SMBRSA
!      En sortie de l'etape on conserve TINSSA, DIVU, SMBRSA
!===============================================================================

iccocg = 1
inc = 1

! ON S'APPUIE SUR LES TABLEAUX DE TRAVAIL W4,W5,W6
!     ON UTILISE EGALEMENT SMBRSA ET TINSSA COMME
!                           TABLEAUX DE TRAVAIL TEMPORAIRES
!     W1,W2,W3 SONT UTILISES DANS GRDCEL
!     TINSSA RECOIT OMEGA**2
!        TINSSA = {          DUDY**2 - 2*DUDY*DVDX + DUDZ**2 -2*DUDZ*DWDX
!                  DVDX**2                         + DVDZ**2 -2*DVDZ*DWDY
!                  DWDX**2 + DWDY**2
!               = 2 Oij.Oij



nswrgp = nswrgr(iuiph)
imligp = imligr(iuiph)
iwarnp = iwarni(inuiph)
epsrgp = epsrgr(iuiph)
climgp = climgr(iuiph)
extrap = extrag(iuiph)
iphydp = 0

! SMBRSA  = DUDX ,W4 = DUDY ,W5 = DUDZ

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   iuiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iuiph)   , coefa(1,icliup) , coefb(1,icliup) ,          &
   smbrsa , w4     , w5     ,                                     &
!  ------   ------   ------
   w1     , w2     , w3     ,                                     &
   ra     )


! TINSSA = DUDY**2 + DUDZ**2
! DIVU   = DUDX

do iel = 1, ncel
  tinssa(iel)   = w4(iel)**2 + w5(iel)**2
  divu  (iel)   = smbrsa(iel)
enddo

!               ,W4     = DUDY ,W5     = DUDZ

! W6     = DVDX ,W7     = DVDY ,W8     = DVDZ

nswrgp = nswrgr(iviph)
imligp = imligr(iviph)
iwarnp = iwarni(inuiph)
epsrgp = epsrgr(iviph)
climgp = climgr(iviph)
extrap = extrag(iviph)
iphydp = 0

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   iviph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iviph)   , coefa(1,iclivp) , coefb(1,iclivp) ,          &
!  dpdx     dpdy     dpdz
   w6     , w7     , w8     ,                                     &
!  ------   ------   ------
   w1     , w2     , w3     ,                                     &
   ra     )

! W6     = DVDX ,W7     = DVDY ,W8     = DVDZ

! TINSSA =           DUDY**2 - 2*DUDY*DVDX + DUDZ**2 -2*DUDZ*DWDX
!          DVDX**2                         + DVDZ**2
! DIVU = DUDX + DVDY

do iel = 1, ncel
  tinssa(iel)   = tinssa(iel) - 2.d0*w4(iel)*w6(iel)              &
           + w6(iel)**2 + w8(iel)**2
  divu  (iel)   = divu(iel) + w7   (iel)
enddo

!               ,              ,W5     = DUDZ
!               ,              ,W8     = DVDZ
! W4     = DWDX ,W6     = DWDY ,W7     = DWDZ

nswrgp = nswrgr(iwiph)
imligp = imligr(iwiph)
iwarnp = iwarni(inuiph)
epsrgp = epsrgr(iwiph)
climgp = climgr(iwiph)
extrap = extrag(iwiph)
iphydp = 0

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   iwiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iwiph)   , coefa(1,icliwp) , coefb(1,icliwp) ,          &
   w4     , w6     , w7     ,                                     &
!  ------   ------   ------
   w1     , w2     , w3     ,                                     &
   ra     )

! TINSSA = 2 OMEGA**2 =
!                            DUDY**2 - 2*DUDY*DVDX + DUDZ**2 -2*DUDZ*DWDX
!                  DVDX**2                         + DVDZ**2 -2*DVDZ*DWDY
!                  DWDX**2 + DWDY**2
! DIVU = DUDX + DVDY + DWDZ

do iel = 1, ncel
  tinssa (iel)   = tinssa(iel) - 2.d0*w8(iel)*W6(iel)             &
           +  w4(iel)**2 +w6(iel)**2
  divu  (iel)   = divu(iel) + w7   (iel)
enddo

! On libere   SMBRSA,W1,W2,W3,W4,W5,W6,W7,W8
! On garde tinssa et divu

! CALCUL DE GRAD nusa
! W4     = DNUDX ,W5 = DNUDY ,W6  = DNUDZ

nswrgp = nswrgr(inuiph)
imligp = imligr(inuiph)
iwarnp = iwarni(inuiph)
epsrgp = epsrgr(inuiph)
climgp = climgr(inuiph)
extrap = extrag(inuiph)
iphydp = 0

iclvar = iclrtp(inuiph,icoef)

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
!  ------
   inuiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w1     , w1     , w1     ,                                     &
!  ------            ------            ------
   rtpa(1,inuiph)  , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   w4     , w5     , w6     ,                                     &
!  ------   ------   ------
   w1     , w2     , w3     ,                                     &
   ra     )
! SMBRSA = GRADnu**2

do iel = 1, ncel
  smbrsa(iel)   = w4(iel)**2 + w5(iel)**2 + w6(iel)**2
enddo

! On libere
! W1,..,W9



!===============================================================================
! 3. PRISE EN COMPTE DES TERMES SOURCES UTILISATEURS
!

!      On passe 2 Omega**2 = TINSSA et la divergence DIVU
!      Tableaux de travail                        W1, W2, W3, W4, W5, W6
!                                VISCF VISCB XAM DRTP SMBRSA W8 W9
!      La partie a expliciter est stockee dans    W7
!      La partie a impliciter est stockee dans    DAM
!      En sortie de l'etape on conserve           TINSSA, DIVU,
!                                                 W7, DAM
!===============================================================================
do iel = 1, ncel
  dam(iel) = 0.d0
  w7 (iel) = 0.d0
enddo

maxelt = max(ncelet, nfac, nfabor)
ils    = idebia
ifinia = ils + maxelt
CALL IASIZE('TURBSA',IFINIA)

call ustssa                                                       &
!==========
 ( ifinia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   maxelt , ia(ils),                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , tinssa , divu   ,          &
   w7     , dam    ,                                              &
!  ------   ------
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   drtp   , smbrsa , w8     , w9     ,                            &
   ra     )

! On libere W1, W2, W3, W4, W5, W6, W8, W9,
!           VISCF, VISCB, XAM, DRTP, SMBRSA,

!===============================================================================
! 4. CALCUL DU TERME DE GRAVITE
!===============================================================================

! Gravity is not taken into account at the moment


!===============================================================================
! 5. TERME D'ACCUMULATION DE MASSE -(dRO/dt)*Volume

!      Le terme est stocke dans         W1
!      En sortie de l'etape on conserve W1, TINSSA, DIVU, SMBRSA, W7, DAM
!===============================================================================

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w1)

!===============================================================================
! 6. ON FINALISE LE CALCUL DES TERMES SOURCES

!      Les termes sont stockes dans     SMBRSA
!      En sortie de l'etape on conserve W1, TINSTSA, DIVU,
!                                       SMBRSA,
!                                       W7 , DAM
!===============================================================================

! Herebelow, we only handle  the case where all the walls have the same roughness
! To extend it, we should be able to link every fluid cell to a boundary face
! (and then give it the appropriate roughness value)

ipatrg = 0
dsa0   = -999.d0

iclvar = iclrtp(inuiph,icoef)
do ifac = 1, nfabor
  if ( itypfb(ifac).eq.iparug ) then
    ipatrg = 1
    cofbnu = coefb(ifac,iclvar)
    ! Roughness of the wall
    dsa0   = distb(ifac) *cofbnu/(1.d0-cofbnu)
    hssa   = exp(8.5d0*xkappa)*dsa0
  endif
  if(ipatrg.ne.0) goto 100
enddo
 100    continue

if(irangp.ge.0) then
  call parcpt(ipatrg)
  if(ipatrg.ne.0) then
    call parsom(dsa0)
    dsa0=dsa0/ipatrg
  endif
endif

!     Si on extrapole les termes sources et rho  , il faut ici rho^n
!                                        et visct, il faut ici visct^n
do iel = 1, ncel

  visct = propce(iel,ipcvto)
  rom   = propce(iel,ipcroo)
  ! Kinematic viscosity
  nu0   = propce(iel,ipcvis)/rom
  distbf= ra(idipar+iel-1)
  ! viscosity of SA
  nusa  = rtpa(iel,inuiph)
  chi   = nusa/nu0
  ! If we have a rough wall
  if(ipatrg.ne.0) then
    distbf = distbf + dsa0
    chi  = chi + 0.5d0* hssa/distbf
  endif
  chi3  = chi**3
  fv1   = chi3/(chi3 + cv13 )
  fv2   = 1.d0 - nusa /(nu0 + nusa*fv1)
  taussa= sqrt(tinssa(iel)*0.5d0)+nusa/(xkappa*distbf)**2*fv2

  ! Computation of fw
  rsa   = min( nusa/(taussa*(xkappa*distbf)**2),10.D0)
  gsa   = rsa + csaw2*(rsa**6-rsa)
  fw    = gsa*( (1.D0+csaw3**6)/(gsa**6+csaw3**6))**(1.D0/6.D0)

  ! SMBRSA = Grad nu . Grad nu
  smbrsa(iel) = volume(iel)*rom*(                                 &
  !  1/SIGMA
  !  -----
     dsigma * csab2*smbrsa(iel)+csab1*taussa*nusa-csaw1*fw*(nusa/distbf)**2)

  ! implicitation of the negative source term of the SA equation.
  ! NB : this term could be negative, and if so, then we explicit it.
  tinssa(iel) = (max(csaw1*fw*nusa/distbf**2-csab1*taussa,0.d0)         &
                      )*rom*volume(iel)

enddo

!===============================================================================
! 7. PRISE EN COMPTE DES TERMES SOURCES UTILISATEURS
!                        ET ACCUMULATION DE MASSE    : PARTIE EXPLICITE
!      On utilise                       W1,  W7, DAM
!      Le terme est stocke dans         SMBRSA
!      En sortie de l'etape on conserve W1, TINSSA, DIVU,
!                                       SMBRSA
!                                       W7, DAM
!===============================================================================

!     Si on extrapole les T.S.
if(isto2t.gt.0) then

  do iel = 1, ncel

!       Sauvegarde de Ts^(n-1) (Term Utilisateur EXPlicite Nusa)
     tuexpn =propce(iel,iptsta)

!       Pour la suite et le pas de temps suivant
!       On stoque les TS explicites du temps n (TS model + TS utilisateur)
    propce(iel,iptsta) = smbrsa(iel) + w7(iel)

!       Termes dependant de la variable resolue et theta PROPCE

!                               Div(rhoU)*nusa^n
!                               -------   ---------------
    smbrsa(iel) = iconv(inuiph)*w1(iel)  *rtpa(iel,inuiph)            &
!        -Thetas*PROPCE^(n-1)
!          ----- ------
         - thets*tuexpn

!       On suppose -DAM > 0 : on implicite
!         le terme utilisateur dependant de la variable resolue

!                 Ts_imp  * nusa^n
!                 --------  ---------------
    smbrsa(iel) = dam(iel)*rtpa(iel,inuiph) + smbrsa(iel)

  enddo

!     Si on n'extrapole pas les T.S. : W7 --> TS explicite
else
  do iel = 1, ncel
    smbrsa(iel) = smbrsa(iel) + dam(iel)*rtpa(iel,inuiph) + w7(iel)  &
         +iconv(inuiph)*w1(iel)*rtpa(iel,inuiph)
  enddo
endif

!===============================================================================
! 8 PRISE EN COMPTE DES TERMES SOURCES LAGRANGIEN : PARTIE EXPLICITE
!     COUPLAGE RETOUR
!===============================================================================

! Not accounted for at the moment.

!===============================================================================
! 9. AJOUT DES TERMES SOURCES DE MASSE EXPLICITES

!       Les parties implicites eventuelles sont conservees dans W2 et W3
!         et utilisees dans la phase d'implicitation cv/diff

!       Les termes sont stockes dans     SMBRSA, W2, W3
!       En sortie de l'etape on conserve W1, TINSSA, DIVU,
!                                        SMBRSA,
!                                        DAM, W9, W2, W3
!===============================================================================

if (ncesmp.gt.0) then

  do iel = 1, ncel
    w2(iel) = 0.d0
    w3(iel) = 0.d0
  enddo

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

!       On incremente SMBRSA par -Gamma RTPA et ROVSDT par Gamma (*theta)
  ivar = inuiph

  call catsma                                                     &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
                                 isto2t , thetv        ,   &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipriph) ,    &
   smbrsa , w2     , w4 )

!       Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if(isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta ) = propce(iel,iptsta ) + w4(iel)
    enddo
!       Sinon on le met directement dans SMBRSA
  else
    do iel = 1, ncel
      smbrsa(iel) = smbrsa(iel) + w4(iel)
    enddo
  endif

endif

!     ON LIBERE                       W4

!     Finalisation des termes sources
if(isto2t.gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
!                               (1+thetas)* PROPCE^n
!                               ------      ------------------
    smbrsa(iel) = smbrsa(iel) + thetp1    * propce(iel,iptsta)
  enddo
endif

!===============================================================================
! 10. TERMES INSTATIONNAIRES

!     On utilise                       W1, W2, W3, W7, W8
!                                      DAM, W9
!     Les termes sont stockes dans     TINSSA
!     En sortie de l'etape on conserve SMBRSA, TINSSA
!===============================================================================

! --- PARTIE EXPLICITE

! --- RHO/DT et DIV
!     Extrapolation ou non, meme forme par coherence avec bilsc2

do iel = 1, ncel
  rom = propce(iel,ipcrom)
  romvsd = rom*volume(iel)/dt(iel)

! TINSSA already contains the negativ implicited source term
  tinssa(iel) = tinssa(iel)                                        &
               +istat(inuiph)*romvsd                               &
               -iconv(inuiph)*w1(iel)*thetv
enddo

! --- Source de masse (le theta est deja inclus par catsma)
if (ncesmp.gt.0) then
  do iel = 1, ncel
    tinssa(iel) = tinssa(iel) + w2(iel)
  enddo
endif

!----------------------------------
! --- Termes sources utilisateurs?
!... Implicitation des TS?
if(isto2t.gt.0) then
  do iel = 1, ncel
    tinssa(iel) = tinssa(iel) -dam(iel)*thetv
  enddo
else
  do iel = 1, ncel
    tinssa(iel) = tinssa(iel) + max(-dam(iel),zero)
  enddo
endif

!===============================================================================
! 11. RESOLUTION

!       On utilise                      SMBRSA, TINSSA,
!       Tableaux de travail             W1, W2, W3, W4, W5, W6
!===============================================================================

! ---> Traitement de nusa

ivar = inuiph
iclvar = iclrtp(ivar,icoef )
iclvaf = iclrtp(ivar,icoeff)

ipp    = ipprtp(ivar)

!    "VITESSE" DE DIFFUSION FACETTE

if( idiff(ivar).ge. 1 ) then

  do iel = 1, ncel
    rom = propce(iel,ipcrom)
    ! diffusibility: 1/SIGMA*(mu_laminaire+ rho*nusa)
    ! nusa  = rtpa(iel,inuiph)
    w1(iel) = dsigma *( propce(iel,ipcvis)                        &
                        + idifft(ivar)*rtpa(iel,inuiph)*rom )
  enddo

  call viscfa                                                     &
  !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  ,                                              &
   ra     )

  ! Be carefull with the walls:
  !  If we have a smooth wall then nusa is zero at the wall
  !  If we have a rough wall then nusa_wall*(1- IprF/d0)=Vipr

  do ifac = 1, nfabor

    iel   = ifabor(ifac)
    surfn = surfbn(ifac)

    ! Smooth wall
    if(    itypfb(ifac).eq.iparoi) then
      viscb(ifac) = dsigma * propce(iel,ipcvis)*surfn/distb(ifac)

    ! Rough wall
    elseif(itypfb(ifac).eq.iparug) then

      rom = propce(iel,ipcrom)
      ! dsa0 is recomputed in case of many different roughness
      cofbnu = coefb(ifac,iclvar)
      ! Roughness of the wall
      dsa0   = distb(ifac) *cofbnu/(1.d0-cofbnu)
      hssa   = exp(8.5d0*xkappa)*dsa0
      ! For rough walls: nusa_F*(IprF/d0+1) = nusa_Ipr
      viscb(ifac) = dsigma * ( propce(iel,ipcvis)                    &
                   + idifft(ivar)*rtpa(iel,inuiph)*rom               &
                   * dsa0/(distb(ifac)+dsa0)            )*surfn/distb(ifac)

    endif

  enddo

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

!     RESOLUTION

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
!  ------   ------
   tinssa , smbrsa , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )


!===============================================================================
! 12. CLIPPING
!===============================================================================

iclip = 0

iwarnp = iwarni(inuiph)
call clipsa                                                       &
!==========
 ( ncelet , ncel   , nvar   , nphas  ,                            &
   iclip  , iwarnp ,                                              &
   propce , rtp    )


!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** PHASE ',I4,' RESOLUTION DE SPALART-ALLMARAS            ',/,&
'      ------------------------------------                   ',/)
 1100 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)

#else

 1000 format(/,                                                   &
'   ** PHASE ',I4,' SOLVING SPALART-A'                         ,/,&
'      ------------------------------'                         ,/)
 1100 format(1X,A8,' : EXPLICIT BALANCE = ',E14.5)
#endif

!----
! FIN
!----

return

end subroutine
