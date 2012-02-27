!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine turbsa &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   tslagr , coefa  , coefb  , ckupdc , smacel ,                   &
   itypfb )

!===============================================================================
! Purpose:
! --------

! Solving op the equation of nusa, which is the scalar quantity defined by
! the Spalart-Allmaras model for 1 time-step.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
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
use mesh
!We have to know if there is any rough wall
use parall
use pointe, only: dispar, coefau, coefbu

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          itypfb(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)

! Local variables

character*80     chaine
integer          iel   , ifac  , init  , inc   , iccocg, ivar
integer          iivar , iiun
integer          iclip , isqrt
integer          nswrgp, imligp
integer          icliup
integer          iclvar, iclvaf
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          ipcrom, ipbrom, ipcvst, ipcvis, iflmas, iflmab
integer          iwarnp, ipp
integer          iptsta
integer          ipcroo, ipbroo, ipcvto, ipcvlo
integer          ipatrg

logical          ilved

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

double precision rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:) :: smbrsa, tinssa, divu
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:,:,:) :: gradv
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4
double precision, allocatable, dimension(:) :: w7

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(dam(ncelet))
allocate(smbrsa(ncelet))
allocate(tinssa(ncelet), divu(ncelet))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet))
allocate(w7(ncelet))

icliup = iclrtp(iu,icoef)

ipcrom = ipproc(irom  )
ipcvst = ipproc(ivisct)
ipcvis = ipproc(iviscl)
iflmas = ipprof(ifluma(iu))
iflmab = ipprob(ifluma(iu))
ipbrom = ipprob(irom  )

! S pour source, V pour variable
!terme source grandeur turbulente
thets  = thetst

ivar   = inusa
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
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(ncelet,3,3))

iccocg = 1
inc = 1

nswrgp = nswrgr(iu)
imligp = imligr(iu)
iwarnp = iwarni(inusa)
epsrgp = epsrgr(iu)
climgp = climgr(iu)
extrap = extrag(iu)

if (ivelco.eq.1) then

  ilved = .false.

  call grdvec &
  !==========
( iu     , imrgra , inc    , iccocg , nswrgp , imligp ,          &
  iwarnp , nfecra ,                                              &
  epsrgp , climgp , extrap ,                                     &
  ilved ,                                                        &
  rtpa(1,iu) ,  coefau , coefbu,                                 &
  gradv  )

else

  call grdvni &
  !==========
( iu  , imrgra , inc    , iccocg , nswrgp , imligp ,             &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
  rtpa(1,iu)   , coefa(1,icliup) , coefb(1,icliup) ,             &
  gradv  )

endif


! TINSSA =  OMEGA**2 = DUDY**2 + DVDX**2 + DUDZ**2 + DWDX**2 + DVDZ**2 + DWDY**2
!                     - 2*DUDY*DVDX - 2*DUDZ*DWDX - 2*DVDZ*DWDY
!
!        = 2 Oij.Oij
! DIVU = DUDX + DVDY + DWDZ

do iel = 1, ncel
  tinssa(iel) = (gradv(iel,2,1) - gradv(iel,1,2))**2   &
              + (gradv(iel,3,1) - gradv(iel,1,3))**2   &
              + (gradv(iel,3,2) - gradv(iel,2,3))**2
  divu(iel) = gradv(iel,1,1) + gradv(iel,2,2) + gradv(iel,3,3)
enddo

! Free memory
deallocate(gradv)

! Allocate a temporary array for the gradient calculation
allocate(grad(ncelet,3))

! CALCUL DE GRAD nusa

nswrgp = nswrgr(inusa)
imligp = imligr(inusa)
iwarnp = iwarni(inusa)
epsrgp = epsrgr(inusa)
climgp = climgr(inusa)
extrap = extrag(inusa)

iclvar = iclrtp(inusa,icoef)

call grdcel &
!==========
 ( inusa , imrgra , inc    , iccocg , nswrgp , imligp ,           &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,inusa)  , coefa(1,iclvar) , coefb(1,iclvar) ,           &
   grad   )

! SMBRSA = GRADnu**2

do iel = 1, ncel
  smbrsa(iel) = grad(iel,1)**2 + grad(iel,2)**2 + grad(iel,3)**2
enddo

! Free memory
deallocate(grad)

!===============================================================================
! 3. PRISE EN COMPTE DES TERMES SOURCES UTILISATEURS
!

!      On passe 2 Omega**2 = TINSSA et la divergence DIVU
!      Tableaux de travail                        W1, W2, W3, W4, W5, W6
!                                VISCF VISCB SMBRSA W8 W9
!      La partie a expliciter est stockee dans    W7
!      La partie a impliciter est stockee dans    DAM
!      En sortie de l'etape on conserve           TINSSA, DIVU,
!                                                 W7, DAM
!===============================================================================
do iel = 1, ncel
  dam(iel) = 0.d0
  w7 (iel) = 0.d0
enddo

call ustssa                                                       &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , tinssa , divu   ,          &
   w7     , dam    )

! On libere W1, W2, W3, W4, W5, W6, W8, W9,
!           VISCF, VISCB, SMBRSA,

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

iclvar = iclrtp(inusa,icoef)
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
  distbf= dispar(iel)
  ! viscosity of SA
  nusa  = rtpa(iel,inusa)
  chi   = nusa/nu0
  ! If we have a rough wall
  if(ipatrg.ne.0) then
    distbf = distbf + dsa0
    chi  = chi + 0.5d0* hssa/distbf
  endif
  chi3  = chi**3
  fv1   = chi3/(chi3 + cv13 )
  fv2   = 1.d0 - nusa /(nu0 + nusa*fv1)
  taussa= max(sqrt(tinssa(iel))+nusa/(xkappa*distbf)**2*fv2, epzero)

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
    smbrsa(iel) = iconv(inusa)*w1(iel)  *rtpa(iel,inusa)            &
!        -Thetas*PROPCE^(n-1)
!          ----- ------
         - thets*tuexpn

!       On suppose -DAM > 0 : on implicite
!         le terme utilisateur dependant de la variable resolue

!                 Ts_imp  * nusa^n
!                 --------  ---------------
    smbrsa(iel) = dam(iel)*rtpa(iel,inusa) + smbrsa(iel)

  enddo

!     Si on n'extrapole pas les T.S. : W7 --> TS explicite
else
  do iel = 1, ncel
    smbrsa(iel) = smbrsa(iel) + dam(iel)*rtpa(iel,inusa) + w7(iel)  &
         +iconv(inusa)*w1(iel)*rtpa(iel,inusa)
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
  ivar = inusa

  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
                                 isto2t , thetv        ,          &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipr) ,       &
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
               +istat(inusa)*romvsd                                &
               -iconv(inusa)*w1(iel)*thetv
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

ivar = inusa
iclvar = iclrtp(ivar,icoef )
iclvaf = iclrtp(ivar,icoeff)

ipp    = ipprtp(ivar)

!    "VITESSE" DE DIFFUSION FACETTE

if( idiff(ivar).ge. 1 ) then

  do iel = 1, ncel
    rom = propce(iel,ipcrom)
    ! diffusibility: 1/SIGMA*(mu_laminaire+ rho*nusa)
    ! nusa  = rtpa(iel,inusa)
    w1(iel) = dsigma *( propce(iel,ipcvis)                        &
                        + idifft(ivar)*rtpa(iel,inusa)*rom )
  enddo

  call viscfa                                                     &
  !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

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
                   + idifft(ivar)*rtpa(iel,inusa)*rom                &
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
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
!  ------   ------
   tinssa , smbrsa , rtp(1,ivar)     ,                            &
   rvoid  )


!===============================================================================
! 12. CLIPPING
!===============================================================================

iclip = 0

iwarnp = iwarni(inusa)
call clipsa                                                       &
!==========
 ( ncelet , ncel   , nvar   ,                                     &
   iclip  , iwarnp ,                                              &
   propce , rtp    )


! Free memory
deallocate(viscf, viscb)
deallocate(dam)
deallocate(smbrsa)
deallocate(tinssa, divu)
deallocate(w1, w2, w3)
deallocate(w4)
deallocate(w7)

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
