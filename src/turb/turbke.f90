!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

!===============================================================================
! Function:
! ---------

!> \file turbke.f90
!>
!> \brief Solving the \f$ k - \epsilon \f$ for incompressible flows or slightly
!> compressible flows for one time step.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        index of the ncepdp cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        mass source type for the variables (cf. ustsma)
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rtp           calculated variables at cell centers
!>                               (at the current time step)
!> \param[in]     rtpa          calculated variables at cell centers
!>                               (at the previous time step)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfa        physical properties at interior face centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[in]     tslagr        coupling term of the lagangian module
!> \param[in]     coefa, coefb  boundary conditions
!>
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        values of the variables associated to the
!>                               mass source
!>                               (for ivar=ipr, smacel is the mass flux)
!> \param[in]     prdv2f        production term stored for the v2f
!_______________________________________________________________________________

subroutine turbke &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   tslagr , coefa  , coefb  , ckupdc , smacel ,                   &
   prdv2f )

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
use lagran
use ppincl
use pointe, only: coefau, coefbu
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision prdv2f(ncelet)

! Local variables

character*80     chaine
integer          iel   , ifac  , init  , inc   , iccocg, ivar
integer          iivar , iiun
integer          iclip , isqrt
integer          nswrgp, imligp
integer          icliup, iclivp, icliwp
integer          iclvar, iclvaf
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          ipcrom, ipbrom, ipcvst, ipcvis, iflmas, iflmab
integer          iwarnp, ipp
integer          iptsta
integer          ipcroo, ipbroo, ipcvto, ipcvlo
integer          iphydp
integer          imucpp

double precision rnorm , d2s3, divp23
double precision deltk , delte, a11, a12, a22, a21
double precision gravke, epssuk, unsdet, romvsd
double precision prdtur, xk, xeps, xphi, xnu, xnut, ttke, ttmin, tt
double precision visct , rom   , ceps1 , ctsqnu
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision thetp1, thetak, thetae, thets, thetap
double precision tuexpk, tuexpe
double precision cmueta, sqrcmu, xs
double precision hint

logical          ilved

double precision rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: usimpk
double precision, allocatable, dimension(:) :: smbrk, smbre, rovsdt
double precision, allocatable, dimension(:) :: tinstk, tinste, divu
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5
double precision, allocatable, dimension(:) :: w7, w8, usimpe
double precision, allocatable, dimension(:) :: w10, w11, w12
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:,:,:), allocatable :: gradv
double precision, allocatable, dimension(:) :: dpvar

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbrk(ncelet), smbre(ncelet), rovsdt(ncelet))
allocate(tinstk(ncelet), tinste(ncelet), divu(ncelet))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet))
allocate(usimpk(ncelet))
allocate(w7(ncelet), w8(ncelet), usimpe(ncelet))
allocate(dpvar(ncelet))

if (iturb.eq.51) then
  allocate(w10(ncelet),w11(ncelet))
endif

icliup = iclrtp(iu,icoef)
iclivp = iclrtp(iv,icoef)
icliwp = iclrtp(iw,icoef)

ipcrom = ipproc(irom  )
ipcvst = ipproc(ivisct)
ipcvis = ipproc(iviscl)
iflmas = ipprof(ifluma(iu))
iflmab = ipprob(ifluma(iu))
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

if(iwarni(ik).ge.1) then
  if (iturb.eq.20) then
    write(nfecra,1000)
  else if (iturb.eq.21) then
    write(nfecra,1001)
  else
    write(nfecra,1002)
  endif
endif

! For the model with linear production, sqrt(Cmu) is required
sqrcmu = sqrt(cmu)

!===============================================================================
! 2. Compute the scalar strain rate SijSij and the trace of the velocity
!    gradient

!      SijSij est stocke dans           tinstk
!      DivU est stocke dans             divu
!      En sortie de l'etape on conserve tinstk, divu
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(ncelet,3,3))

iccocg = 1
inc = 1

nswrgp = nswrgr(iu)
imligp = imligr(iu)
iwarnp = iwarni(ik)
epsrgp = epsrgr(iu)
climgp = climgr(iu)
extrap = extrag(iu)

if (ivelco.eq.1) then

  ilved = .false.

  call grdvec &
  !==========
( iu     , imrgra , inc    , nswrgp , imligp ,                   &
  iwarnp , nfecra ,                                              &
  epsrgp , climgp , extrap ,                                     &
  ilved  ,                                                       &
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

! tinstk = Stain rate = ( 2 (S11)**2 + 2 (S22)**2 + 2 (S33)**2
!                       + (2 S12)**2 + (2 S13)**2 + (2 S23)**2 )
!        = 2 Sij.Sij
! divu = dudx + dvdy + dwdz

do iel = 1, ncel

  tinstk(iel) = 2.d0*( gradv(iel,1,1)**2 + gradv(iel,2,2) **2 &
                     + gradv(iel,3,3)**2 )                    &
              + (gradv(iel,2,1) + gradv(iel,1,2))**2          &
              + (gradv(iel,3,1) + gradv(iel,1,3))**2          &
              + (gradv(iel,3,2) + gradv(iel,2,3))**2

  divu(iel) = gradv(iel,1,1) + gradv(iel,2,2) + gradv(iel,3,3)

enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3. Take user source terms into account

!    The scalar strain rate (tinstk) and the trace of the velocity gradient
!     (divu) are available.
!
!    La partie a expliciter est stockee dans    w7, w8
!    La partie a impliciter est stockee dans    usimpk, usimpe
!    En sortie de l'etape on conserve           tinstk, divu,
!                                               w7 , w8, usimpk, usimpe
!===============================================================================

do iel = 1, ncel
  usimpk(iel) = 0.d0
  usimpe(iel) = 0.d0
  w7(iel) = 0.d0
  w8(iel) = 0.d0
enddo

call ustske                                                       &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , tinstk , divu   ,          &
   w7     , w8     , usimpk , usimpe )

!===============================================================================
! 4. Add - 2/3 div(u) * div(u)

!      En sortie de l'etape on conserve tinstk, divu,
!                                       w7 , w8, usimpk, usimpe
!===============================================================================

! Dans le cas de la production lineaire, seul le terme en divu est
! multiplie par visct. Pour les autres modeles, la multiplication par
! visct sera faite ulterieurement.
! A ce stade, tinstk contient S**2
d2s3 = 2.d0/3.d0
if (iturb.eq.21) then
  do iel = 1, ncel
    rom   = propce(iel,ipcroo)
    visct = propce(iel,ipcvto)
    xs = sqrt(tinstk(iel))
    cmueta = cmu*rtpa(iel,ik)/rtpa(iel,iep)*xs
    cmueta = min(cmueta,sqrcmu)
    tinstk(iel) = rom*cmueta*xs*rtpa(iel,ik)                   &
         - d2s3*visct*divu(iel)*divu(iel)
  enddo
else
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) - d2s3*divu(iel)*divu(iel)
  enddo
endif

!===============================================================================
! 5. Compute the buoyant term

!      Les s.m. recoivent production et termes de gravite
!      Tableaux de travail              viscb
!      Les s.m. sont stockes dans       tinstk, tinste
!      En sortie de l'etape on conserve tinstk, tinste,
!                                       divu,
!                                       w7 , w8, usimpk, usimpe
!===============================================================================

if (igrake.eq.1 .and. ippmod(iatmos).ge.1) then

    !  Calcul du terme de gravite pour la version atmospherique

    call atprke &
    !==========
 ( nscal  ,                                                       &
   ipcvto,                                                        &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   tinstk , tinste )

else if (igrake.eq.1) then

  ! Allocate a temporary for the gradient calculation
  allocate(grad(ncelet,3))

! --- Terme de gravite G = BETA*G*GRAD(SCA)/PRDTUR/RHO
!     Ici on calcule   G =-G*GRAD(RHO)/PRDTUR/RHO

  iccocg = 1
  inc = 1

!     Le choix ci dessous a l'avantage d'etre simple

  nswrgp = nswrgr(ik)
  epsrgp = epsrgr(ik)
  imligp = imligr(ik)
  iwarnp = iwarni(ik)
  climgp = climgr(ik)
  extrap = extrag(ik)

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
   propce(1,ipcroo), propfb(1,ipbroo), viscb  ,                   &
   grad   )


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
      gravke = -(grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz) / &
                (propce(iel,ipcroo)*prdtur)
      tinste(iel) = tinstk(iel) + propce(iel,ipcvto)*max(gravke,zero)
      tinstk(iel) = tinstk(iel) + propce(iel,ipcvto)*gravke
    enddo
  else
    do iel = 1, ncel
      gravke = -(grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz) / &
                (propce(iel,ipcroo)*prdtur)
      tinste(iel) = tinstk(iel) + max( gravke,zero )
      tinstk(iel) = tinstk(iel) + gravke
    enddo
  endif

  ! Free memory
  deallocate(grad)

else


! --- Production sans termes de gravite
!       tinstk=tinste=P

  do iel = 1, ncel
    tinste(iel) = tinstk(iel)
  enddo

endif

! En v2f, on stocke tinstk dans prdv2f qui sera complete plus loin pour
! contenir le terme de production complet
if (itytur.eq.5) then
  do iel = 1, ncel
    prdv2f(iel) = tinstk(iel)
  enddo
endif

! On libere viscb

!===============================================================================
! 7. pre Seulement pour le modele bl-v2/k, calcul de e et ceps2*

!      Les termes sont stockes dans     w10, w11
!      Tableaux de travail              w2, w3, w4, w5, drtp,smbrk,smbre
!                                       viscf, viscb
!      En sortie de l'etape on conserve w10, w11
!===============================================================================

if (iturb.eq.51) then

  ! Calcul du terme CEPS2*: Il est stocke dans w10

  do iel=1,ncel
    visct = propce(iel,ipcvto)
    rom   = propce(iel,ipcroo)
    w3(iel) = visct/rom/sigmak
  enddo

  call viscfa &
  !==========
( imvisf ,        &
  w3     ,        &
  viscf  , viscb  )

  ivar = ik
  iclvar = iclrtp(ivar,icoef)
  iclvaf = iclrtp(ivar,icoeff)

  ! Translate coefa into cofaf and coefb into cofbf
  do ifac = 1, nfabor

    iel = ifabor(ifac)

    hint = w3(iel)/distb(ifac)

    ! Translate coefa into cofaf and coefb into cofbf
    coefa(ifac, iclvaf) = -hint*coefa(ifac,iclvar)
    coefb(ifac, iclvaf) = hint*(1.d0-coefb(ifac,iclvar))

  enddo

  iccocg = 1
  inc = 1
  init = 1

  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  iphydp = 0

  call itrgrp &
  !==========
( nvar   , nscal  ,                                              &
  init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydp , &
  iwarnp , nfecra ,                                              &
  epsrgp , climgp , extrap ,                                     &
  w2     , w2     , w2     ,                                     &
  rtpa(1,ivar)    ,                                              &
  coefa(1,iclvar) , coefb(1,iclvar) ,                            &
  coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
  viscf  , viscb  ,                                              &
  w3     , w3     , w3     ,                                     &
  w10    )

  do iel=1,ncel
    w10(iel) = -w10(iel)/volume(iel)/rtpa(iel,iep)
    w10(iel) = tanh(abs(w10(iel))**1.5d0)
    w10(iel) = cpale2*(1.d0-(cpale2-cpale4)/cpale2*w10(iel)*rtpa(iel,ial)**3)
  enddo

!  Calcul du terme 2*Ceps3*(1-alpha)^3*nu*nut/eps*d2Ui/dxkdxj*d2Ui/dxkdxj:
!   (i.e. E term / k)           : Il est stocke dans W11

  ! Allocate a work array
  allocate(w12(ncelet))

  call tsepls &
  !==========
( dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
  coefa  , coefb  ,                                              &
  w12    )

  do iel=1,ncel

    rom   = propce(iel,ipcroo)
    xnu   = propce(iel,ipcvlo)/rom
    xnut  = propce(iel,ipcvto)/rom
    xeps = rtpa(iel,iep )
    xk   = rtpa(iel,ik )
    xphi = rtpa(iel,iphi)

    ttke = xk/xeps
    ttmin = cpalct*sqrt(xnu/xeps)
    tt = sqrt(ttke**2 + ttmin**2)

    w11(iel) = 2.d0*xnu*xnut*w12(iel)*cpale3/xeps                  &
                *(1.d0-rtpa(iel,ial))**3

  enddo

  ! Free memory
  deallocate(w12)

endif

!===============================================================================
! 7. On finalise le calcul des termes sources

!      Les termes sont stockes dans     smbrk, smbre
!      En sortie de l'etape on conserve tinstk, tinste, divu,
!                                       smbrk, smbre
!                                       w7 , w8, usimpk, usimpe
!===============================================================================

! smbre = ceps1 epsilon/k (prod + g ) - rho0 volume epsilon epsilon/k
! smbrk =                  prod + g   - rho0 volume epsilon

! Si on extrapole les termes sources et rho  , il faut ici rho^n
!                                    et visct, il faut ici visct^n

if (iturb.eq.20) then

  do iel = 1, ncel

    visct = propce(iel,ipcvto)
    rom   = propce(iel,ipcroo)

    smbrk(iel) = volume(iel)*(                                    &
         visct*tinstk(iel)                                        &
         -d2s3*rom*rtpa(iel,ik)*divu(iel)                         &
         -rom*rtpa(iel,iep) )

    smbre(iel) = volume(iel)*rtpa(iel,iep)/rtpa(iel,ik)*(         &
         ce1*( visct*tinste(iel)                                  &
         -d2s3*rom*rtpa(iel,ik)*divu(iel) )                       &
         -ce2*rom*rtpa(iel,iep) )

  enddo

else if (iturb.eq.21) then

  do iel = 1, ncel

    rom   = propce(iel,ipcroo)

    smbrk(iel) = volume(iel)*(                                    &
         tinstk(iel)                                              &
         -d2s3*rom*rtpa(iel,ik)*divu(iel)                         &
         -rom*rtpa(iel,iep) )

    smbre(iel) = volume(iel)*rtpa(iel,iep)/rtpa(iel,ik)*(         &
         ce1*(tinste(iel)                                         &
         -d2s3*rom*rtpa(iel,ik)*divu(iel) )                       &
         -ce2*rom*rtpa(iel,iep) )

  enddo

else if (iturb.eq.50) then

  do iel = 1, ncel

    visct = propce(iel,ipcvto)
    rom   = propce(iel,ipcroo)
    xeps = rtpa(iel,iep )
    xk   = rtpa(iel,ik )
    xphi = rtpa(iel,iphi)
    xphi = max(xphi,epzero)
    xnu  = propce(iel,ipcvlo)/rom
    ceps1= 1.4d0*(1.d0+cv2fa1*sqrt(1.d0/xphi))
    ttke = xk / xeps
    ttmin = cv2fct*sqrt(xnu/xeps)
    tt = max(ttke,ttmin)

    smbrk(iel) = volume(iel)*(                                    &
         visct*tinstk(iel)                                        &
         -d2s3*rom*rtpa(iel,ik)*divu(iel)                         &
         -rom*rtpa(iel,iep) )

    smbre(iel) = volume(iel)/tt*(                                 &
         ceps1*( visct*tinste(iel)                                &
         -d2s3*rom*rtpa(iel,ik)*divu(iel) )                       &
         -cv2fe2*rom*rtpa(iel,iep) )

!     On stocke la partie en Pk dans PRDV2F pour etre reutilise dans RESV2F
    prdv2f(iel) = visct*prdv2f(iel)                               &
         -d2s3*rom*rtpa(iel,ik)*divu(iel)

  enddo

else if (iturb.eq.51) then

  do iel=1,ncel

    visct = propce(iel,ipcvto)
    rom   = propce(iel,ipcroo)
    xeps = rtpa(iel,iep )
    xk   = rtpa(iel,ik )
    xphi = rtpa(iel,iphi)
    xnu  = propce(iel,ipcvlo)/rom
    ttke = xk / xeps
    ttmin = cpalct*sqrt(xnu/xeps)
    tt = sqrt(ttke**2.d0+ttmin**2.d0)

    smbrk(iel) = volume(iel)*(                                    &
         visct*tinstk(iel)                                        &
         -d2s3*rom*rtpa(iel,ik)*divu(iel)                      &
         -rom*rtpa(iel,iep)                                     &
         -rom*w11(iel)*xk )

    smbre(iel) = volume(iel)/tt*(                                 &
         cpale1*( visct*tinste(iel)                               &
         -d2s3*rom*rtpa(iel,ik)*divu(iel) )                    &
         -w10(iel)*rom*rtpa(iel,iep) )

!     On stocke la partie en Pk dans PRDV2F pour etre reutilise dans RESV2F
    prdv2f(iel) = visct*prdv2f(iel)                               &
         -d2s3*rom*rtpa(iel,ik)*divu(iel)

  enddo

endif


!===============================================================================
! 8. Prise en compte des termes sources utilisateurs
!                        et accumulation de masse    : partie explicite
!      On utilise                       w7, w8, usimpk, usimpe
!      Les termes sont stockes dans     smbrk, smbre
!      En sortie de l'etape on conserve tinstk, tinste, divu,
!                                       smbrk, smbre
!                                       usimpk, usimpe

!    Remarque : l'extrapolation telle qu'elle est ecrite n'a pas grand
!               sens si IKECOU=1
!===============================================================================

! Si on extrapole les T.S.
if(isto2t.gt.0) then

  do iel = 1, ncel

    ! Sauvegarde pour echange
    tuexpk = propce(iel,iptsta)
    ! Pour la suite et le pas de temps suivant
    propce(iel,iptsta) = smbrk(iel) + w7(iel)
    ! Termes dependant de la variable resolue et theta PROPCE
    smbrk(iel) = - thets*tuexpk
    ! On suppose -usimpk > 0 : on implicite
    !  le terme utilisateur dependant de la variable resolue
    smbrk(iel) = usimpk(iel)*rtpa(iel,ik) + smbrk(iel)

    ! Sauvegarde pour echange
    tuexpe = propce(iel,iptsta+1)
    ! Pour la suite et le pas de temps suivant
    propce(iel,iptsta+1) = smbre(iel) + w8(iel)
    ! Termes dependant de la variable resolue et theta PROPCE
    smbre(iel) = - thets*tuexpe
    ! On suppose -usimpe > 0 : on implicite
    !  le terme utilisateur dependant de la variable resolue
    smbre(iel) =  usimpe(iel)*rtpa(iel,iep) + smbre(iel)

  enddo

! Si on n'extrapole pas les T.S.
else
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + usimpk(iel)*rtpa(iel,ik) + w7(iel)
    smbre(iel) = smbre(iel) + usimpe (iel)*rtpa(iel,iep) + w8(iel)
  enddo
endif

!===============================================================================
! 8.1 Prise en compte des termes sources lagrangien : partie explicite
!     couplage retour
!===============================================================================

! Ordre 2 non pris en compte
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  do iel = 1,ncel

    ! Termes sources explicte et implicte sur k
    smbrk(iel)  = smbrk(iel) + tslagr(iel,itske)

    ! Termes sources explicte sur Eps
    smbre(iel)  = smbre(iel)                                    &
                + ce4 *tslagr(iel,itske) *rtpa(iel,iep)         &
                                         /rtpa(iel,ik)

  enddo

endif

!===============================================================================
! 9. Prise en compte des termes de conv/diff dans le second membre

!      Tableaux de travail              w4, w5
!      Les termes sont stockes dans     w7 et w8, puis ajoutes a smbrk, smbre
!      En sortie de l'etape on conserve tinstk, tinste, divu,
!                                       smbrk, smbre
!                                       usimpk, w7, w8, usimpe
!===============================================================================

! Ceci ne sert a rien si ikecou n'est pas egal a 1

if (ikecou.eq.1) then

  do iel = 1, ncel
    w7 (iel) = 0.d0
    w8 (iel) = 0.d0
  enddo

! ---> Traitement de k

  ivar   = ik

  ipp    = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef )
  iclvaf = iclrtp(ivar,icoeff)
  chaine = nomvar(ipp)

  if( idiff(ivar).ge. 1 ) then

    do iel = 1, ncel
      if(iturb.eq.51) then
        w4(iel) = propce(iel,ipcvis)/2.d0                         &
             + idifft(ivar)*propce(iel,ipcvst)/sigmak
      else
        w4(iel) = propce(iel,ipcvis)                              &
             + idifft(ivar)*propce(iel,ipcvst)/sigmak
      endif
    enddo
    call viscfa                                                   &
    !==========
 ( imvisf ,                                                       &
   w4     ,                                                       &
   viscf  , viscb  )

    if (iturb.eq.51) then
      ! Translate coefa into cofaf and coefb into cofbf
      do ifac = 1, nfabor

        iel = ifabor(ifac)

        hint = w4(iel)/distb(ifac)

        ! Translate coefa into cofaf and coefb into cofbf
        coefa(ifac, iclvaf) = -hint*coefa(ifac,iclvar)
        coefb(ifac, iclvaf) = hint*(1.d0-coefb(ifac,iclvar))

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
  imucpp = 0
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = thetav(ivar)

  call bilsc2                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp , imucpp ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  , rvoid  , &
   w7     )

  if (iwarni(ivar).ge.2) then
    isqrt = 1
    call prodsc(ncel,isqrt,smbrk,smbrk,rnorm)
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif


! ---> Traitement de epsilon

  ivar   = iep

  ipp    = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef )
  iclvaf = iclrtp(ivar,icoeff)
  chaine = nomvar(ipp)

  if( idiff(ivar).ge. 1 ) then
    do iel = 1, ncel
      if(iturb.eq.51) then
        w4(iel) = propce(iel,ipcvis)/2.0                          &
             + idifft(ivar)*propce(iel,ipcvst)/cpalse
      else
        w4(iel) = propce(iel,ipcvis)                              &
             + idifft(ivar)*propce(iel,ipcvst)/sigmae
      endif
    enddo

    call viscfa &
 ( imvisf ,                                                       &
   w4     ,                                                       &
   viscf  , viscb  )

    if (iturb.eq.51) then
      ! Translate coefa into cofaf and coefb into cofbf
      do ifac = 1, nfabor

        iel = ifabor(ifac)

        hint = w4(iel)/distb(ifac)

        ! Translate coefa into cofaf and coefb into cofbf
        coefa(ifac, iclvaf) = -hint*coefa(ifac,iclvar)
        coefb(ifac, iclvaf) = hint*(1.d0-coefb(ifac,iclvar))

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
  imucpp = 0
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = thetav(ivar)

  call bilsc2                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp , imucpp ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  , rvoid  , &
   w8     )

  if (iwarni(ivar).ge.2) then
    isqrt = 1
    call prodsc(ncel,isqrt,smbre,smbre,rnorm)
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif

  do iel = 1,ncel
    smbrk(iel) = smbrk(iel) + w7(iel)
    smbre(iel) = smbre(iel) + w8(iel)
  enddo

endif

!===============================================================================
! 10. Ajout des termes sources de masse explicites

!       Les parties implicites eventuelles sont conservees dans w2 et w3
!         et utilisees dans la phase d'implicitation cv/diff

!       Les termes sont stockes dans     smbrk, smbre, w2, w3
!       En sortie de l'etape on conserve tinstk, tinste, divu,
!                                        smbrk, smbre
!                                        usimpk, usimpe, w2, w3
!===============================================================================

if (ncesmp.gt.0) then

  do iel = 1, ncel
    w2(iel) = 0.d0
    w3(iel) = 0.d0
  enddo

  ! Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

  ! On incremente smbrs par -Gamma rtpa et rovsdt par Gamma (*theta)
  ivar = ik

  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
   isto2t , thetav(ivar)    ,                                     &
   icetsm , itypsm(1,ivar)  ,                                     &
   volume , rtpa(1,ivar)    , smacel(1,ivar) , smacel(1,ipr) ,    &
   smbrk  , w2     , w4 )

  ivar = iep

  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
   isto2t , thetav(ivar)    ,                                     &
   icetsm , itypsm(1,ivar)  ,                                     &
   volume , rtpa(1,ivar)    , smacel(1,ivar) , smacel(1,ipr) ,    &
   smbre  , w3     , w5 )

  ! Si on extrapole les TS on met Gamma Pinj dans propce
  if(isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta  ) = propce(iel,iptsta  ) + w4(iel)
      propce(iel,iptsta+1) = propce(iel,iptsta+1) + w5(iel)
    enddo
  ! Sinon on le met directement dans smbr
  else
    do iel = 1, ncel
      smbrk(iel) = smbrk(iel) + w4(iel)
      smbre(iel) = smbre(iel) + w5(iel)
    enddo
  endif

endif

! On libere                       w4, w5

! Finalisation des termes sources
if(isto2t.gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + thetp1 * propce(iel,iptsta)
    smbre(iel) = smbre(iel) + thetp1 * propce(iel,iptsta+1)
  enddo
endif

!===============================================================================
! 11. Increments des termes sources dans le second membre

!       On utilise                       tinstk, tinste, divu
!       Les termes sont stockes dans     smbrk, smbre
!       En sortie de l'etape on conserve smbrk, smbre,
!                                        usimpk, usimpe, w2, w3, w7, w8
!===============================================================================

! Ordre 2 non pris en compte
if(ikecou.eq.1) then

  if (iturb.eq.20) then

    do iel = 1, ncel

      rom = propce(iel,ipcrom)

      ! Coupled solving
      romvsd=1.d0/(rom*volume(iel))
      smbrk(iel)=smbrk(iel)*romvsd
      smbre(iel)=smbre(iel)*romvsd
      divp23= d2s3*max(divu(iel),zero)

      epssuk = rtpa(iel,iep)/rtpa(iel,ik)

      a11 = 1.d0/dt(iel)                                          &
           -2.d0*rtpa(iel,ik)/rtpa(iel,iep)                       &
           *cmu*min(tinstk(iel),zero)+divp23
      a12 = 1.d0
      a21 = -ce1*cmu*tinste(iel)-ce2*epssuk*epssuk
      a22 = 1.d0/dt(iel)+ce1*divp23                               &
           +2.d0*ce2*epssuk

      unsdet = 1.d0/(a11*a22 -a12*a21)

      deltk = ( a22*smbrk(iel) -a12*smbre(iel) )*unsdet
      delte = (-a21*smbrk(iel) +a11*smbre(iel) )*unsdet

      ! New source term for the iterative process
      romvsd = rom*volume(iel)/dt(iel)

      smbrk(iel) = romvsd*deltk
      smbre(iel) = romvsd*delte

    enddo

  ! Dans verini on bloque la combinaison iturb=21/ikecou=1
  else if (iturb.eq.21) then

    write(nfecra,*)'ikecou=1 non valide en k-eps prod lin'
    call csexit (1)
  !  Section non totalement validee (a priori ca marche, mais pas trop stable) :
  !  en fait le v2f est meilleur avec ikecou=0, on bloque donc la combinaison
  !  iturb=50/ikecou=1 au niveau de verini. Ces lignes sont donc inaccessibles.
  !  On les laisse au cas ou .....
  else if (iturb.eq.50) then

    do iel = 1, ncel

      rom = propce(iel,ipcrom)

      ! Coupled solving
      romvsd=1.d0/(rom*volume(iel))
      smbrk(iel)=smbrk(iel)*romvsd
      smbre(iel)=smbre(iel)*romvsd
      divp23= d2s3*max(divu(iel),zero)

      xeps = rtpa(iel,iep )
      xk   = rtpa(iel,ik )
      xphi = rtpa(iel,iphi)
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
        ! Pour a12 on fait comme en k-eps standard pour l'instant,
        ! on ne prend pas le terme en P+G ... est-ce judicieux ?
        a12 = 1.d0
        a21 = -ceps1*cv2fmu*xphi*tinste(iel)-cv2fe2*epssuk*epssuk
        a22 = 1.d0/dt(iel)+ceps1*divp23                           &
             +2.d0*cv2fe2*epssuk
      else
        a11 = 1.d0/dt(iel)                                        &
             -cv2fmu*xphi*ctsqnu*min(tinstk(iel),zero)/sqrt(xeps) &
             +divp23
        ! Pour a12 on fait comme en k-eps standard pour l'instant,
        ! on ne prend pas le terme en P+G ... est-ce judicieux ?
        a12 = 1.d0
        ! Le terme en divp23 dans a21 n'est pas forcement judicieux
        ! (a-t-on besoin du MAX ?)
        a21 = -ceps1*cv2fmu*xphi*tinste(iel)                      &
             +ceps1*sqrt(xeps)/ctsqnu*divp23
        a22 = 1.d0/dt(iel)+1.d0/2.d0*ceps1*divp23*xk              &
             /ctsqnu/sqrt(xeps)                                   &
             +3.d0/2.d0*cv2fe2/ctsqnu*sqrt(xeps)
      endif

      unsdet = 1.d0/(a11*a22 -a12*a21)

      deltk = ( a22*smbrk(iel) -a12*smbre(iel) )*unsdet
      delte = (-a21*smbrk(iel) +a11*smbre(iel) )*unsdet

      ! New source term for the iterative process
      romvsd = rom*volume(iel)/dt(iel)

      smbrk(iel) = romvsd*deltk
      smbre(iel) = romvsd*delte

    enddo

  ! Dans verini on bloque la combinaison iturb=51/ikecou=1
  else if (iturb.eq.51) then

    write(nfecra,*)'ikecou=1 non valide en bl-v2/k'
    call csexit (1)

  endif

endif

! On libere                       tinstk, tinste, divu

!===============================================================================
! 12. Termes instationnaires

!     On utilise                       w2, w3, w7, w8
!                                      usimpk, usimpe
!     Les termes sont stockes dans     tinstk, tinste
!     En sortie de l'etape on conserve smbrk, smbre,  tinstk, tinste
!===============================================================================

! --- Explicit part

!     on enleve la convection/diffusion au temps n a smbrk et smbre
!     si on les avait calcules
if (ikecou.eq.1) then
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) - w7(iel)
    smbre(iel) = smbre(iel) - w8(iel)
  enddo
endif

! --- RHO/DT

do iel = 1, ncel
  rom = propce(iel,ipcrom)
  romvsd = rom*volume(iel)/dt(iel)
  tinstk(iel) = istat(ik)*romvsd
  tinste(iel) = istat(iep)*romvsd
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
  thetak = thetav(ik)
  thetae = thetav(iep)
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) -usimpk(iel)*thetak
    tinste(iel) = tinste(iel) -usimpe (iel)*thetae
  enddo
else
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) + max(-usimpk(iel),zero)
    tinste(iel) = tinste(iel) + max(-usimpe (iel),zero)
  enddo
endif

! --- prise en compte des termes lagrangien : couplage retour

! Ordre 2 non pris en compte
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  do iel = 1,ncel

! Termes sources implicite sur k

    tinstk(iel) = tinstk(iel) + max(-tslagr(iel,itsli),zero)

! Termes sources implicte sur Eps

    tinste(iel) = tinste(iel)                                     &
          + max( (-ce4*tslagr(iel,itske)/rtpa(iel,ik)) , zero)

  enddo

endif

! Si ikecou=0, on implicite plus fortement k et eps

if(ikecou.eq.0)then
  if(itytur.eq.2)then
    do iel=1,ncel
      xeps = rtpa(iel,iep )
      xk   = rtpa(iel,ik )
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
      xeps = rtpa(iel,iep )
      xk   = rtpa(iel,ik )
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
  else if(iturb.eq.51)then
    do iel=1,ncel
      xeps = rtpa(iel,iep )
      xk   = rtpa(iel,ik )
      rom = propce(iel,ipcrom)
      xnu  = propce(iel,ipcvis)/rom
      ttke = xk / xeps
      ttmin = cpalct*sqrt(xnu/xeps)
      tt = sqrt(ttke**2.d0+ttmin**2.d0)
      if(xk.gt.1.d-12) then
        tinstk(iel) = tinstk(iel) +                               &
             rom*volume(iel)/ttke
      endif
      tinstk(iel) = tinstk(iel) +                                 &
             rom*w11(iel)*volume(iel)
      tinste(iel) = tinste(iel) +                                 &
           w10(iel)*rom*volume(iel)/tt
    enddo

  endif
endif

! On libere w2, w3, usimpk, usimpe

!===============================================================================
! 13. Solving

!       On utilise                      smbrk, smbre,  tinstk, tinste
!       Tableaux de travail             w1
!===============================================================================

! ---> turbulent kinetic (k) energy treatment

ivar = ik
iclvar = iclrtp(ivar,icoef )
iclvaf = iclrtp(ivar,icoeff)

ipp    = ipprtp(ivar)

! "vitesse" de diffusion facette

if( idiff(ivar).ge. 1 ) then

  do iel = 1, ncel
    if(iturb.eq.51) then
      w1(iel) = propce(iel,ipcvis)/2.d0                           &
                          + idifft(ivar)*propce(iel,ipcvst)/sigmak
    else
      w1(iel) = propce(iel,ipcvis)                                &
                          + idifft(ivar)*propce(iel,ipcvst)/sigmak
    endif
  enddo
  call viscfa                                                     &
  !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

! Solving k

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
imucpp = 0
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

call codits &
!==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp ,                            &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   tinstk , smbrk  , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , rvoid  )

! ---> Turbulent dissipation (epsilon) treatment

ivar = iep
iclvar = iclrtp(ivar,icoef )
iclvaf = iclrtp(ivar,icoeff)

ipp    = ipprtp(ivar)


! "vitesse" de diffusion facette

if( idiff(ivar).ge. 1 ) then
  do iel = 1, ncel
    if(iturb.eq.51) then
      w1(iel) = propce(iel,ipcvis)/2.d0                           &
                          + idifft(ivar)*propce(iel,ipcvst)/cpalse
    else
      w1(iel) = propce(iel,ipcvis)                                &
                          + idifft(ivar)*propce(iel,ipcvst)/sigmae
    endif
  enddo

  call viscfa &
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

! Solving epsilon

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
imucpp = 0
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

call codits &
!==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp ,                            &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   tinste , smbre  , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , rvoid  )

!===============================================================================
! 14. Clipping
!===============================================================================

iclip = 1
iwarnp = iwarni(ik)
call clipke                                                       &
!==========
 ( ncelet , ncel   , nvar   ,                                     &
   iclip  , iwarnp ,                                              &
   propce , rtp    )

! Free memory
deallocate(viscf, viscb)
deallocate(usimpk)
deallocate(smbrk, smbre, rovsdt)
deallocate(tinstk, tinste, divu)
deallocate(w1, w2, w3)
deallocate(w4, w5)
deallocate(w7, w8, usimpe)
deallocate(dpvar)

if (allocated(w10)) deallocate(w10, w11)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                      &
'   ** Resolution du k-epsilon                   ',/,&
'      -----------------------                   ',/)
 1001 format(/,                                      &
'   ** Resolution du k-epsilon a prod lineaire   ',/,&
'      ---------------------------------------   ',/)
 1002 format(/,                                      &
'   ** Resolution du v2f (k et epsilon)          ',/,&
'      --------------------------------          ',/)
 1100 format(1X,A8,' : Bilan explicite = ',E14.5)

#else

 1000 format(/,                                      &
'   ** Solving k-epsilon'                         ,/,&
'      -----------------'                         ,/)
 1001 format(/,                                      &
'   ** solving k-epsilon with linear prod'        ,/,&
'      ----------------------------------'        ,/)
 1002 format(/,                                      &
'   ** Solving v2f (k and epsilon)'               ,/,&
'      ---------------------------'               ,/)
 1100 format(1X,A8,' : Explicit balance = ',E14.5)
#endif

!----
! End
!----

return

end subroutine
