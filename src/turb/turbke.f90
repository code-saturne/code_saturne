!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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
!> \param[in]     tslagr        coupling term of the lagangian module
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        values of the variables associated to the
!>                               mass source
!>                               (for ivar=ipr, smacel is the mass flux)
!> \param[in]     prdv2f        production term stored for the v2f
!_______________________________________________________________________________

subroutine turbke &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce ,                            &
   tslagr , ckupdc , smacel , prdv2f )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use cstnum
use cstphy
use optcal
use lagran
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision tslagr(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision prdv2f(ncelet)

! Local variables

character*80     chaine
integer          iel   , ifac  , init  , inc   , iccocg, ivar
integer          iivar , iiun
integer          iclip , isqrt
integer          nswrgp, imligp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          ipcvst, ipcvis, iflmas, iflmab
integer          iwarnp, ipp
integer          iptsta
integer          ipcvto, ipcvlo
integer          iphydp
integer          imucpp, idftnp, iswdyp

integer          icvflb
integer          ivoid(1)

double precision rnorm , d2s3, d1s3, divp23
double precision deltk , delte, a11, a12, a22, a21
double precision gravke, epssuk, unsdet, romvsd
double precision prdtur, xk, xeps, xphi, xnu, xnut, ttke, ttmin, tt
double precision visct , rho   , ceps1
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
double precision, allocatable, dimension(:) :: prdtke, prdeps
double precision, allocatable, dimension(:) :: strain
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5
double precision, allocatable, dimension(:) :: w7, w8, usimpe
double precision, allocatable, dimension(:) :: w10, w11, w12
double precision, allocatable, dimension(:) :: ce2rc
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:,:,:), allocatable :: gradv
double precision, allocatable, dimension(:) :: dpvar
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, crom, bromo, cromo
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbrk(ncelet), smbre(ncelet), rovsdt(ncelet))
allocate(tinstk(ncelet), tinste(ncelet), divu(ncelet), strain(ncelet))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet))
allocate(usimpk(ncelet))
allocate(w7(ncelet), w8(ncelet), usimpe(ncelet))
allocate(dpvar(ncelet))

if (iturb.eq.20) then
  allocate(prdtke(ncelet), prdeps(ncelet))
endif

if (iturb.eq.51) then
  allocate(w10(ncelet),w11(ncelet))
endif

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

call field_get_val_s(icrom, crom)
ipcvst = ipproc(ivisct)
ipcvis = ipproc(iviscl)

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_val_s(ibrom, brom)

thets  = thetst

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)
call field_get_val_s(icrom, cromo)
call field_get_val_s(ibrom, bromo)
ipcvto = ipcvst
ipcvlo = ipcvis
if(isto2t.gt.0) then
  if (iroext.gt.0) then
    call field_get_val_prev_s(icrom, cromo)
    call field_get_val_prev_s(ibrom, bromo)
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

d2s3 = 2.d0/3.d0
d1s3 = 1.d0/3.d0

!===============================================================================
! 2. Compute the scalar strain rate SijSij and the trace of the velocity
!    gradient

!      (Sij^D) (Sij^D)  is stored in    strain (deviatoric strain tensor rate)
!      tr(Grad u)       is stored in    divu
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(3, 3, ncelet))

iccocg = 1
inc = 1

nswrgp = nswrgr(iu)
imligp = imligr(iu)
iwarnp = iwarni(ik)
epsrgp = epsrgr(iu)
climgp = climgr(iu)
extrap = extrag(iu)

ilved = .false.

! WARNING: gradv(xyz, uvw, iel)
call grdvec &
!==========
( iu     , imrgra , inc    , nswrgp , imligp ,                   &
  iwarnp , epsrgp , climgp ,                                     &
  ilved  ,                                                       &
  rtpa(1,iu) ,  coefau , coefbu,                                 &
  gradv  )

! strain = Stain rate of the deviatoric part of the strain tensor
!        = 2 (Sij^D).(Sij^D)
! divu   = trace of the velocity gradient
!        = dudx + dvdy + dwdz

do iel = 1, ncel

  strain(iel) = 2.d0                                                           &
    *( ( d2s3*gradv(1, 1, iel)-d1s3*gradv(2, 2, iel)-d1s3*gradv(3, 3, iel))**2 &
     + (-d1s3*gradv(1, 1, iel)+d2s3*gradv(2, 2, iel)-d1s3*gradv(3, 3, iel))**2 &
     + (-d1s3*gradv(1, 1, iel)-d1s3*gradv(2, 2, iel)+d2s3*gradv(3, 3, iel))**2 &
     )                                                                         &
    + (gradv(2, 1, iel) + gradv(1, 2, iel))**2                                 &
    + (gradv(3, 1, iel) + gradv(1, 3, iel))**2                                 &
    + (gradv(3, 2, iel) + gradv(2, 3, iel))**2

  divu(iel) = gradv(1, 1, iel) + gradv(2, 2, iel) + gradv(3, 3, iel)

enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3. Instationnary terms (stored in tinstk and tinste)
!===============================================================================

do iel = 1, ncel
  rho = crom(iel)
  romvsd = rho*volume(iel)/dt(iel)
  tinstk(iel) = istat(ik)*romvsd
  tinste(iel) = istat(iep)*romvsd
enddo

!===============================================================================
! 4. Compute the first part of the production term: muT (S^D)**2

!      En sortie de l'etape on conserve strain, divu,
!===============================================================================

! For the Linear Production k-epsilon model,
! the production term is assumed to be asymptotically in S and
! not in mu_TxS**2
if (iturb.eq.21) then
  do iel = 1, ncel
    rho   = cromo(iel)
    visct = propce(iel,ipcvto)
    xs = sqrt(strain(iel))
    cmueta = min(cmu*rtpa(iel,ik)/rtpa(iel,iep)*xs, sqrcmu)
    smbrk(iel) = rho*cmueta*xs*rtpa(iel,ik)
    smbre(iel) = smbrk(iel)
  enddo
else
  do iel = 1, ncel
    visct = propce(iel,ipcvto)
    smbrk(iel) = visct*strain(iel)
    smbre(iel) = smbrk(iel)
  enddo
endif

!=============================================================================
! 5. Take into account rotation/curvature correction, if necessary
!=============================================================================

! Cazalbou correction: the Ceps2 coefficient of destruction term of epsislon
! is modified by rotation and curvature

! Allocate an array for the modified Ceps2 coefficient
allocate(ce2rc(ncel))

if (irccor.eq.1) then

  ! Compute the modified Ceps2 coefficient (w1 array not used)

  call rotcor(dt, rtpa, w1, ce2rc)
  !==========

else

  if (itytur.eq.2) then
    do iel = 1, ncel
      ce2rc(iel) = ce2
    enddo
  elseif(iturb.eq.50) then
    do iel = 1, ncel
      ce2rc(iel) = cv2fe2
    enddo
  elseif(iturb.eq.51) then
    do iel = 1, ncel
      ce2rc(iel) = ccaze2
    enddo
  endif

endif
! ce2rc array is used all along the subroutine. It is deallocated at the end.

!===============================================================================
! 6. Compute the buoyancy term

!      Les s.m. recoivent production et termes de gravite
!      Tableaux de travail              viscb
!      Les s.m. sont stockes dans       smbrk, smbre
!      En sortie de l'etape on conserve smbrk, smbre,
!                                       divu,
!===============================================================================

! Buoyant term for the Atmospheric module
! (function of the potential temperature)
if (igrake.eq.1 .and. ippmod(iatmos).ge.1) then

  call atprke &
  !==========
 ( nscal  ,                                                       &
   rtpa   , propce ,                                              &
   tinstk ,                                                       &
   smbrk  , smbre  )

! --- Buoyancy term     G = Beta*g.Grad(scalar)/prdtur/rho
!     Here is computed  G =-g.grad(rho)/prdtur/rho
else if (igrake.eq.1) then

  ! Allocate a temporary for the gradient calculation
  allocate(grad(ncelet,3))

  iccocg = 1
  inc = 1
  nswrgp = nswrgr(ik)
  epsrgp = epsrgr(ik)
  imligp = imligr(ik)
  iwarnp = iwarni(ik)
  climgp = climgr(ik)
  extrap = extrag(ik)

  ! Dirichlet boundary condition on the gradient of rho
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

  iivar = 0

  call grdcel &
  !==========
 ( iivar  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   cromo  , bromo  , viscb  ,                                     &
   grad   )

  ! Production term due to buoyancy
  !   smbrk = P+G
  !   smbre = P+(1-ce3)*G
  if(iscalt.gt.0.and.nscal.ge.iscalt) then
    prdtur = sigmas(iscalt)
  else
    prdtur = 1.d0
  endif

  ! smbr* store mu_TxS**2
  do iel = 1, ncel
    rho   = cromo(iel)
    visct = propce(iel,ipcvto)
    xeps = rtpa(iel, iep)
    xk   = rtpa(iel, ik)
    ttke = xk / xeps

    gravke = -(grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz) &
           / (rho*prdtur)

    ! Implicit Buoyant terms when negativ
    tinstk(iel) = tinstk(iel) + max(-rho*volume(iel)*cmu*ttke*gravke, 0.d0)

    ! Explicit Buoyant terms
    smbre(iel) = smbre(iel) + visct*max(gravke, zero)
    smbrk(iel) = smbrk(iel) + visct*gravke
  enddo

  ! Free memory
  deallocate(grad)

endif

! En v2f, on stocke la production  dans prdv2f qui sera complete plus loin pour
! contenir le terme de production complet
if (itytur.eq.5) then
  do iel = 1, ncel
    prdv2f(iel) = smbrk(iel)
  enddo
endif

!===============================================================================
! 7. pre Seulement pour le modele bl-v2/k, calcul de e et ceps2*

!      Les termes sont stockes dans     w10, w11
!      Tableaux de travail              w2, w3
!                                       viscf, viscb
!      En sortie de l'etape on conserve w10, w11
!===============================================================================

if (iturb.eq.51) then

  ! Calcul du terme CEPS2*: Il est stocke dans w10

  do iel=1,ncel
    visct = propce(iel,ipcvto)
    rho   = cromo(iel)
    w3(iel) = visct/rho/sigmak
  enddo

  call viscfa &
  !==========
( imvisf ,        &
  w3     ,        &
  viscf  , viscb  )

  ivar = ik
  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)
  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  ! Translate coefa into cofaf and coefb into cofbf
  do ifac = 1, nfabor

    iel = ifabor(ifac)

    hint = w3(iel)/distb(ifac)

    ! Translate coefa into cofaf and coefb into cofbf
    cofafp(ifac) = -hint*coefap(ifac)
    cofbfp(ifac) = hint*(1.d0-coefbp(ifac))

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
( init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydp , &
  iwarnp ,                                                       &
  epsrgp , climgp , extrap ,                                     &
  rvoid  ,                                                       &
  rtpa(1,ivar)    ,                                              &
  coefap , coefbp , cofafp , cofbfp ,                            &
  viscf  , viscb  ,                                              &
  w3     , w3     , w3     ,                                     &
  w10    )

  do iel = 1, ncel
    w10(iel) = -w10(iel)/volume(iel)/rtpa(iel,iep)
    w10(iel) = tanh(abs(w10(iel))**1.5d0)
    w10(iel) = cpale2*(1.d0-(cpale2-cpale4)/cpale2*w10(iel)*rtpa(iel,ial)**3)
  enddo

  ! Calcul du terme 2*Ceps3*(1-alpha)^3*nu*nut/eps*d2Ui/dxkdxj*d2Ui/dxkdxj:
  !  (i.e. E term / k)           : Il est stocke dans w11

  ! Allocate a work array
  allocate(w12(ncelet))

  call tsepls(rtpa, w12)
  !==========

  do iel = 1, ncel

    rho   = cromo(iel)
    xnu   = propce(iel,ipcvlo)/rho
    xnut  = propce(iel,ipcvto)/rho
    xeps = rtpa(iel,iep )
    xk   = rtpa(iel,ik )
    xphi = rtpa(iel,iphi)

    ttke = xk/xeps
    ttmin = cpalct*sqrt(xnu/xeps)
    tt = sqrt(ttke**2 + ttmin**2)

    w11(iel) = 2.d0*xnu*xnut*w12(iel)*cpale3/xeps                  &
                *(1.d0-rtpa(iel,ial))**3

  enddo

  ! Take into account the Cazalbou rotation/curvature correction if necessary
  if (irccor.eq.1) then
     do iel =1, ncel
       w10(iel) = w10(iel)*ce2rc(iel)/ccaze2
       w11(iel) = w11(iel)*ce2rc(iel)/ccaze2
     enddo
  endif

  ! Free memory
  deallocate(w12)

endif

!===============================================================================
! 8. Finalization of explicit and implicit source terms

!===============================================================================

! smbre = ceps1 epsilon/k (prod + g ) - rho0 volume epsilon epsilon/k
! smbrk =                  prod + g   - rho0 volume epsilon

! Si on extrapole les termes sources et rho  , il faut ici rho^n
!                                    et visct, il faut ici visct^n

if (itytur.eq.2) then

  ! Stores the production terms for the k-epsilon coupling option
  if (iturb.eq.20) then
    do iel = 1, ncel
      prdtke(iel) = smbrk(iel)
      prdeps(iel) = smbre(iel)
    enddo
  endif

  do iel = 1, ncel

    rho   = cromo(iel)

    smbrk(iel) = volume(iel)*                                     &
               ( smbrk(iel) - rho*rtpa(iel,iep)                   &
               - d2s3*rho*rtpa(iel,ik)*divu(iel)                  &
               )

    smbre(iel) = volume(iel)*                                              &
               ( rtpa(iel,iep)/rtpa(iel,ik)*( ce1*smbre(iel)               &
                                            - ce2rc(iel)*rho*rtpa(iel,iep) &
                                            )                              &
               - d2s3*rho*ce1*rtpa(iel,iep)*divu(iel)                      &
               )

  enddo

  ! If the solving of k-epsilon is uncoupled, negative source terms are implicited
  if (ikecou.eq.0) then
    do iel = 1, ncel
      xeps = rtpa(iel,iep )
      xk   = rtpa(iel,ik )
      rho = crom(iel)
      ttke = xk / xeps
      tinstk(iel) = tinstk(iel) + rho*volume(iel)/ttke            &
                  + max(d2s3*rho*volume(iel)*divu(iel), 0.d0)
      tinste(iel) = tinste(iel) + ce2rc(iel)*rho*volume(iel)/ttke &
                  + max(d2s3*ce1*rho*volume(iel)*divu(iel), 0.d0)
    enddo
  endif

else if (iturb.eq.50) then

  do iel = 1, ncel

    visct = propce(iel,ipcvto)
    rho  = cromo(iel)
    xnu  = propce(iel,ipcvlo)/rho
    xeps = rtpa(iel,iep)
    xk   = rtpa(iel,ik)
    xphi = rtpa(iel,iphi)
    xphi = max(xphi, epzero)
    ceps1= 1.4d0*(1.d0+cv2fa1*sqrt(1.d0/xphi))
    ttke = xk / xeps
    ttmin = cv2fct*sqrt(xnu/xeps)
    tt = max(ttke, ttmin)

    ! Explicit part
    smbrk(iel) = volume(iel)*                                     &
               ( smbrk(iel) - rho*rtpa(iel,iep)                   &
               - d2s3*rho*rtpa(iel,ik)*divu(iel)                  &
               )

    smbre(iel) = volume(iel)*                                       &
               ( 1.d0/tt*(ceps1*smbre(iel) - ce2rc(iel)*rho*xeps)   &
               - d2s3*rho*ceps1*xk*divu(iel)                        &
               )

    ! On stocke la partie en Pk dans PRDV2F pour etre reutilise dans RESV2F
    prdv2f(iel) = prdv2f(iel)                               &
                - d2s3*rho*rtpa(iel,ik)*divu(iel)!FIXME this term should be removed

    ! Implicit part
    if (xk.gt.1.d-12) then !FIXME make it dimensionless
      tinstk(iel) = tinstk(iel) + rho*volume(iel)/ttke
    endif
    tinstk(iel) = tinstk(iel) + max(d2s3*rho*volume(iel)*divu(iel), 0.d0)
    tinste(iel) = tinste(iel) + ce2rc(iel)*rho*volume(iel)/tt                &
                + max(d2s3*ceps1*ttke/tt*rho*volume(iel)*divu(iel), 0.d0)

  enddo

else if (iturb.eq.51) then

  do iel=1,ncel

    visct = propce(iel,ipcvto)
    rho   = cromo(iel)
    xnu  = propce(iel,ipcvlo)/rho
    xeps = rtpa(iel,iep )
    xk   = rtpa(iel,ik )
    xphi = rtpa(iel,iphi)
    ttke = xk / xeps
    ttmin = cpalct*sqrt(xnu/xeps)
    tt = sqrt(ttke**2.d0+ttmin**2.d0)

    ! Explicit part
    smbrk(iel) = volume(iel)*                                     &
               ( smbrk(iel)                                       &
               - rho*xeps                                         &
               - rho*w11(iel)*xk                                  &
               - d2s3*rho*xk*divu(iel)                            &
               )

    smbre(iel) = volume(iel)*                                               &
               ( 1.d0/tt*(cpale1*smbre(iel) - w10(iel)*rho*xeps)            &
               - d2s3*rho*cpale1*xk/tt*divu(iel)                            &
               )

    ! On stocke la partie en Pk dans PRDV2F pour etre reutilise dans RESV2F
    prdv2f(iel) = prdv2f(iel)                               &
                - d2s3*rho*rtpa(iel,ik)*divu(iel)!FIXME this term should be removed

    ! Implicit part
    if (xk.gt.1.d-12) then !FIXME make it dimensionless
      tinstk(iel) = tinstk(iel) + rho*volume(iel)/ttke
    endif
    tinstk(iel) = tinstk(iel) + max(d2s3*rho*volume(iel)*divu(iel), 0.d0)
    tinstk(iel) = tinstk(iel) + w11(iel)*rho*volume(iel)
    tinste(iel) = tinste(iel) + w10(iel)*rho*volume(iel)/tt                  &
                + max(d2s3*cpale1*ttke/tt*rho*volume(iel)*divu(iel), 0.d0)

  enddo

endif

!===============================================================================
! 9. Take user source terms into account

!    The scalar strain rate (strain) and the trace of the velocity gradient
!     (divu) are available.
!
!    La partie a expliciter est stockee dans    w7, w8
!    La partie a impliciter est stockee dans    usimpk, usimpe
!    En sortie de l'etape on conserve           strain, divu,
!===============================================================================

do iel = 1, ncel
  usimpk(iel) = 0.d0
  usimpe(iel) = 0.d0
  w7(iel) = 0.d0
  w8(iel) = 0.d0
enddo

call ustske &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smacel , strain , divu   ,                            &
   w7     , w8     , usimpk , usimpe )

! If source terms are extrapolated over time
if (isto2t.gt.0) then

  thetak = thetav(ik)
  thetae = thetav(iep)

  do iel = 1, ncel

    ! Recover the value at time (n-1)
    tuexpk = propce(iel,iptsta)
    tuexpe = propce(iel,iptsta+1)
    ! Save the values for the next time-step
    propce(iel,iptsta) = smbrk(iel) + w7(iel)
    propce(iel,iptsta+1) = smbre(iel) + w8(iel)

    ! Explicit Part
    smbrk(iel) = - thets*tuexpk
    smbre(iel) = - thets*tuexpe
    ! It is assumed that (-usimpk > 0) and though this term is implicit
    smbrk(iel) = usimpk(iel)*rtpa(iel,ik) + smbrk(iel)
    smbre(iel) = usimpe(iel)*rtpa(iel,iep) + smbre(iel)

    ! Implicit part
    tinstk(iel) = tinstk(iel) - usimpk(iel)*thetak
    tinste(iel) = tinste(iel) - usimpe(iel)*thetae

  enddo

! If no extrapolation over time
else
  do iel = 1, ncel
    ! Explicit part
    smbrk(iel) = smbrk(iel) + usimpk(iel)*rtpa(iel,ik) + w7(iel)
    smbre(iel) = smbre(iel) + usimpe(iel)*rtpa(iel,iep) + w8(iel)

    ! Implicit part
    tinstk(iel) = tinstk(iel) + max(-usimpk(iel),zero)
    tinste(iel) = tinste(iel) + max(-usimpe(iel),zero)
  enddo
endif

!===============================================================================
! 10. Prise en compte des termes sources lagrangien
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

    ! Termes sources implicite sur k
    tinstk(iel) = tinstk(iel) + max(-tslagr(iel,itsli),zero)

    ! Termes sources implicte sur Eps
    tinste(iel) = tinste(iel) + max((-ce4*tslagr(iel,itske)/rtpa(iel,ik)), zero)

  enddo

endif

!===============================================================================
! 11. Mass source terms (Implicit and explicit parts)

!       En sortie de l'etape on conserve divu,
!                                        smbrk, smbre
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

  ! Implicit part (theta is already taken into account in catsma)
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) + w2(iel)
    tinste(iel) = tinste(iel) + w3(iel)
  enddo

endif

!===============================================================================
! 12.1 Prise en compte des termes de conv/diff dans le second membre pour le
!      couplage renforcé k-epsilon (ikecou == 1)

!      Tableaux de travail              w4, w5
!      Les termes sont stockes dans     w7 et w8, puis ajoutes a smbrk, smbre
!      En sortie de l'etape on conserve divu,
!                                       smbrk, smbre
!                                       w7, w8
!===============================================================================

if (ikecou.eq.1) then

  do iel = 1, ncel
    w7 (iel) = 0.d0
    w8 (iel) = 0.d0
  enddo

  ! ---> Traitement de k
  ivar   = ik
  call field_get_label(ivarfl(ivar), chaine)

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)
  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  if (idiff(ivar).ge. 1) then

    do iel = 1, ncel
      w4(iel) = propce(iel,ipcvis)                              &
              + idifft(ivar)*propce(iel,ipcvst)/sigmak
    enddo

    call viscfa &
    !==========
 ( imvisf ,                                                       &
   w4     ,                                                       &
   viscf  , viscb  )

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
  idftnp = 1 ! no tensorial diffusivity
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = thetav(ivar)
  ! all boundary convective flux with upwind
  icvflb = 0

  call bilsca &
  !==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   iwarnp , imucpp , idftnp ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , rvoid  , rvoid  ,                            &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   w7     )

  if (iwarni(ivar).ge.2) then
    isqrt = 1
    call prodsc(ncel,isqrt,smbrk,smbrk,rnorm)
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif


  ! ---> Traitement de epsilon
  ivar   = iep
  call field_get_label(ivarfl(ivar), chaine)

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)
  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  if (idiff(ivar).ge. 1) then
    do iel = 1, ncel
      w4(iel) = propce(iel,ipcvis)                              &
              + idifft(ivar)*propce(iel,ipcvst)/sigmae
    enddo

    call viscfa &
    !==========
 ( imvisf ,                                                       &
   w4     ,                                                       &
   viscf  , viscb  )

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
  idftnp = 1 ! no tensorial diffusivity
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = thetav(ivar)
  ! all boundary convective flux with upwind
  icvflb = 0

  call bilsca &
  !==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   iwarnp , imucpp , idftnp ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , rvoid  , rvoid  ,                            &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
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
! 12.2 k-Epsilon coupling (ikecou == 1)

!===============================================================================

! Ordre 2 non pris en compte
if (ikecou.eq.1) then

  if (iturb.eq.20) then

    do iel = 1, ncel

      rho = crom(iel)
      visct = propce(iel,ipcvto)

      ! Coupled solving
      romvsd = 1.d0/(rho*volume(iel))
      smbrk(iel)=smbrk(iel)*romvsd
      smbre(iel)=smbre(iel)*romvsd
      divp23= d2s3*max(divu(iel),zero)

      epssuk = rtpa(iel,iep)/rtpa(iel,ik)

      a11 = 1.d0/dt(iel)                                          &
           -2.d0*rtpa(iel,ik)/rtpa(iel,iep)                       &
           *cmu*min(prdtke(iel)/visct,zero)+divp23
      a12 = 1.d0
      a21 = -ce1*cmu*prdeps(iel)/visct-ce2rc(iel)*epssuk*epssuk
      a22 = 1.d0/dt(iel)+ce1*divp23                               &
           +2.d0*ce2rc(iel)*epssuk

      unsdet = 1.d0/(a11*a22 -a12*a21)

      deltk = ( a22*smbrk(iel) -a12*smbre(iel) )*unsdet
      delte = (-a21*smbrk(iel) +a11*smbre(iel) )*unsdet

      ! New source term for the iterative process
      romvsd = rho*volume(iel)/dt(iel)

      smbrk(iel) = romvsd*deltk
      smbre(iel) = romvsd*delte

    enddo

    ! on enleve la convection/diffusion au temps n a smbrk et smbre
    ! si on les avait calcules
    do iel = 1, ncel
      smbrk(iel) = smbrk(iel) - w7(iel)
      smbre(iel) = smbre(iel) - w8(iel)
    enddo

  ! Dans verini on bloque la combinaison iturb!=20/ikecou=1
  else

    write(nfecra,*)'ikecou=1 non valide avec ce modele de turbulence'
    call csexit (1)

  endif

endif

!===============================================================================
! 13. Finalization of the Right Hand Side when activating 2nd time order
!===============================================================================

if (isto2t.gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + thetp1 * propce(iel,iptsta)
    smbre(iel) = smbre(iel) + thetp1 * propce(iel,iptsta+1)
  enddo
endif

!===============================================================================
! 14. Solving

!       On utilise                      smbrk, smbre,  tinstk, tinste
!       Tableaux de travail             w1
!===============================================================================

! ---> turbulent kinetic (k) energy treatment
ivar = ik
ipp    = ipprtp(ivar)

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

! Face viscosity
if (idiff(ivar).ge.1) then

  do iel = 1, ncel
    if (iturb.eq.51) then
      w1(iel) = propce(iel,ipcvis)/2.d0                           &
              + idifft(ivar)*propce(iel,ipcvst)/sigmak
    else
      w1(iel) = propce(iel,ipcvis)                                &
              + idifft(ivar)*propce(iel,ipcvst)/sigmak
    endif
  enddo

  call viscfa &
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
idftnp = idften(ivar)
iswdyp = iswdyn(ivar)
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
! all boundary convective flux with upwind
icvflb = 0

call codits &
!==========
 ( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   tinstk , smbrk  , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , rvoid  )

! ---> Turbulent dissipation (epsilon) treatment
ivar = iep
ipp    = ipprtp(ivar)

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

! Face viscosity
if (idiff(ivar).ge.1) then
  do iel = 1, ncel
    if (iturb.eq.51) then
      w1(iel) = propce(iel,ipcvis)/2.d0                           &
              + idifft(ivar)*propce(iel,ipcvst)/cpalse
    else
      w1(iel) = propce(iel,ipcvis)                                &
              + idifft(ivar)*propce(iel,ipcvst)/sigmae
    endif
  enddo

  call viscfa &
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
idftnp = idften(ivar)
iswdyp = iswdyn(ivar)
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
! all boundary convective flux with upwind
icvflb = 0

call codits &
!==========
 ( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   tinste , smbre  , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , rvoid  )

!===============================================================================
! 15. Clipping
!===============================================================================

iclip = 1
iwarnp = iwarni(ik)
call clipke &
!==========
 ( ncelet , ncel   , nvar   ,                                     &
   iclip  , iwarnp ,                                              &
   propce , rtp    )

! Free memory
deallocate(viscf, viscb)
deallocate(usimpk)
deallocate(smbrk, smbre, rovsdt)
deallocate(tinstk, tinste, divu, strain)
deallocate(w1, w2, w3)
deallocate(w4, w5)
deallocate(w7, w8, usimpe)
deallocate(dpvar)
deallocate(ce2rc)

if (allocated(w10)) deallocate(w10, w11)
if (allocated(prdtke)) deallocate(prdtke, prdeps)

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
