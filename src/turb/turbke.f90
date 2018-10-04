!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
!>
!> Please refer to the
!> <a href="../../theory.pdf#keps"><b>\f$k - \epsilon \f$ model</b></a>
!> section of the theory guide for more informations, as well as the
!> <a href="../../theory.pdf#turbke"><b>turbke</b></a> section.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        index of the ncepdp cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        mass source type for the variables
!>                              (cf. cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
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
   dt     ,                                                       &
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
use field_operator
use cs_c_bindings
use atincl, only: kopint

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision tslagr(ncelet,*)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision prdv2f(ncelet)

! Local variables

character(len=80) :: chaine
integer          iel   , ifac  , init  , inc   , iccocg, ivar
integer          f_id0 , iiun
integer          iclip
integer          nswrgp, imligp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          iflmas, iflmab
integer          iwarnp
integer          istprv, f_oi_id
integer          iphydp, iprev
integer          imucpp, idftnp, iswdyp
integer          key_t_ext_id
integer          iroext
integer          iviext

integer          icvflb, imasac
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
double precision turb_schmidt
double precision xnoral, xnal(3)
double precision utaurf,ut2,ypa,ya,xunorm, limiter, nu0,alpha
double precision normp

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
double precision, dimension(:), pointer :: cvar_k, cvara_k
double precision, dimension(:), pointer :: cvar_ep, cvara_ep
double precision, dimension(:), pointer :: cvara_al, cvara_phi
double precision, dimension(:), pointer :: cpro_pcvto, cpro_pcvlo
double precision, dimension(:), pointer :: viscl, cvisct
double precision, dimension(:), pointer :: c_st_k_p, c_st_eps_p
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_al

type(var_cal_opt) :: vcopt_k, vcopt_e

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

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

! Map field arrays
call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

call field_get_val_s(ivisct, cvisct)
call field_get_val_s(iviscl, viscl)

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_val_s(ivarfl(ik), cvar_k)
call field_get_val_prev_s(ivarfl(ik), cvara_k)
call field_get_val_s(ivarfl(iep), cvar_ep)
call field_get_val_prev_s(ivarfl(iep), cvara_ep)
if (iturb.eq.50.or.iturb.eq.51) call field_get_val_prev_s(ivarfl(iphi), cvara_phi)
if (iturb.eq.51) call field_get_val_prev_s(ivarfl(ial), cvara_al)

thets  = thetst

call field_get_key_int(ivarfl(ik), kstprv, istprv)
if (istprv.ge.0) then
  call field_get_val_s(istprv, c_st_k_p)
  call field_get_key_int(ivarfl(iep), kstprv, istprv)
  if (istprv.ge.0) then
    call field_get_val_s(istprv, c_st_eps_p)
  endif
  if (istprv.ge.0) istprv = 1
endif

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)
call field_get_val_s(icrom, cromo)
call field_get_val_s(ibrom, bromo)
call field_get_val_s(ivisct, cpro_pcvto)
call field_get_val_s(iviscl, cpro_pcvlo)
if (istprv.ge.0) then
  call field_get_key_int(icrom, key_t_ext_id, iroext)
  if (iroext.gt.0) then
    call field_get_val_prev_s(icrom, cromo)
    call field_get_val_prev_s(ibrom, bromo)
  endif
  call field_get_key_int(iviscl, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_get_val_prev_s(iviscl, cpro_pcvlo)
  endif
  call field_get_key_int(ivisct, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_get_val_prev_s(ivisct, cpro_pcvto)
  endif
endif

call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt_k)
call field_get_key_struct_var_cal_opt(ivarfl(iep), vcopt_e)

if (vcopt_k%iwarni.ge.1) then
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
! 1.1 Advanced reinit
!===============================================================================

! Automatic reinitialization at the end of the first iteration:
! wall distance y^+ is computed with -C log(1-alpha), where C=CL*Ceta*L*kappa,
! then y so we have an idea of the wall distance in complex geometries.
! Then U is initialized with a Reichard layer,
! Epsilon by 1/(kappa y), clipped next to the wall at its value for y^+=15.
! k is given by a blending between eps/(2 nu)*y^2 and utau/sqrt(Cmu)
! The blending function is chosen so that the asymptotic behavior
! and give the correct peak of k (not the same blending than for the EBRSM
! because k profile is not the same for k-omega).

!TODO FIXME: Are the BC uncompatible?
if (ntcabs.eq.1.and.reinit_turb.eq.1.and.iturb.eq.51) then

  allocate(grad(3,ncelet))

  ! Compute the gradient of Alpha
  iprev  = 0
  inc    = 1
  iccocg = 1

  call field_gradient_scalar(ivarfl(ial), iprev, imrgra, inc,     &
                             iccocg,                              &
                             grad)

  call field_get_val_s(ivarfl(iep), cvar_ep)
  call field_get_val_s(ivarfl(ial), cvar_al)
  call field_get_val_v(ivarfl(iu), vel)

  utaurf = 0.05d0*uref
  nu0 = viscl0 / ro0

  do iel = 1, ncel
    ! Compute the velocity magnitude
    xunorm = vel(1,iel)**2 + vel(2,iel)**2 + vel(3,iel)**2
    xunorm = sqrt(xunorm)

    ! y+ is bounded by 400, because in the Reichard profile,
    ! it corresponds to saturation (u>uref)
    cvar_al(iel) = max(min(cvar_al(iel),(1.d0-exp(-400.d0/50.d0))),0.d0)

    call field_current_to_previous(ivarfl(ial))

    ! Compute the magnitude of the alpha gradient
    xnoral = ( grad(1,iel)*grad(1,iel)          &
           +   grad(2,iel)*grad(2,iel)          &
           +   grad(3,iel)*grad(3,iel) )
    xnoral = sqrt(xnoral)
   ! Compute the unitary vector of Alpha
    if (xnoral.le.epzero/cell_f_vol(iel)**(1.d0/3.d0)) then
      xnal(1) = 1.d0/sqrt(3.d0)
      xnal(2) = 1.d0/sqrt(3.d0)
      xnal(3) = 1.d0/sqrt(3.d0)
    else
      xnal(1) = grad(1,iel)/xnoral
      xnal(2) = grad(2,iel)/xnoral
      xnal(3) = grad(3,iel)/xnoral
    endif

    alpha = cvar_al(iel)

    ! Compute YA, therefore alpha is given by 1-exp(-YA/(50 nu/utau))
    ! NB: y^+ = 50 give the best compromise
    ya = -dlog(1.d0-cvar_al(iel))*50.d0*nu0/utaurf
    ypa = ya/(nu0/utaurf)

    ! Velocity magnitude is imposed (limited only), the direction is
    ! conserved
    if (xunorm.le.1.d-12*uref) then
      limiter = 1.d0
    else
      limiter = min( utaurf/xunorm*(2.5d0*dlog(1.d0+0.4d0*ypa)  &
                    +7.8d0*(1.d0-dexp(-ypa/11.d0)               &
                    -(ypa/11.d0)*dexp(-0.33d0*ypa))),           &
                    1.d0)
    endif

    vel(1,iel) = limiter*vel(1,iel)
    vel(2,iel) = limiter*vel(2,iel)
    vel(3,iel) = limiter*vel(3,iel)

    ut2 = 0.05d0*uref

    cvar_ep(iel) = utaurf**3*min(1.d0/(xkappa*15.d0*nu0/utaurf), &
                                 1.d0/(xkappa*ya))
    cvar_k(iel) =  cvar_ep(iel)/2.d0/nu0*ya**2                 &
                   * exp(-ypa/25.d0)**2                        &
                 + ut2**2/sqrt(cmu)*(1.d0-exp(-ypa/25.d0))**2
  enddo

  call field_current_to_previous(ivarfl(iu))
  call field_current_to_previous(ivarfl(ik))
  call field_current_to_previous(ivarfl(iep)) !TODO phi ?

  deallocate(grad)
endif


!===============================================================================
! 2. Compute the scalar strain rate SijSij and the trace of the velocity
!    gradient

!      (Sij^D) (Sij^D)  is stored in    strain (deviatoric strain tensor rate)
!      tr(Grad u)       is stored in    divu
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(3, 3, ncelet))

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

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
! 3. Unsteady terms (stored in tinstk and tinste)
!===============================================================================

do iel = 1, ncel
  rho = crom(iel)
  romvsd = rho*cell_f_vol(iel)/dt(iel)
  tinstk(iel) = vcopt_k%istat*romvsd
  tinste(iel) = vcopt_e%istat*romvsd
enddo

!===============================================================================
! 4. Compute the first part of the production term: muT (S^D)**2

!      Going out of the step we keep strain, divu,
!===============================================================================

! For the Linear Production k-epsilon model,
! the production term is assumed to be asymptotically in S and
! not in mu_TxS**2
if (iturb.eq.21) then
  do iel = 1, ncel
    rho   = cromo(iel)
    visct = cpro_pcvto(iel)
    xs = sqrt(strain(iel))
    cmueta = min(cmu*cvara_k(iel)/cvara_ep(iel)*xs, sqrcmu)
    smbrk(iel) = rho*cmueta*xs*cvara_k(iel)
    smbre(iel) = smbrk(iel)
  enddo
else
  do iel = 1, ncel
    visct = cpro_pcvto(iel)
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

  call rotcor(dt, w1, ce2rc)

else

  if (itytur.eq.2) then
    do iel = 1, ncel
      ce2rc(iel) = ce2
    enddo
  elseif (iturb.eq.50) then
    do iel = 1, ncel
      ce2rc(iel) = cv2fe2
    enddo
  elseif (iturb.eq.51) then
    do iel = 1, ncel
      ce2rc(iel) = ccaze2
    enddo
  endif

endif
! ce2rc array is used all along the subroutine. It is deallocated at the end.

!===============================================================================
! 6. Compute the buoyancy term

!      The mass sources receive production and gravity terms
!      Work arrays                      viscb
!      The mass sources are stored in   smbrk, smbre
!      Going out of the step we keep    smbrk, smbre,
!                                       divu,
!===============================================================================

! Buoyant term for the Atmospheric module
! (function of the potential temperature)
if (igrake.eq.1 .and. ippmod(iatmos).ge.1) then

  call atprke(nscal, tinstk, smbrk, smbre)

! --- Buoyancy term     G = Beta*g.Grad(scalar)/prdtur/rho
!     Here is computed  G =-g.grad(rho)/prdtur/rho
else if (igrake.eq.1) then

  ! Allocate a temporary for the gradient calculation
  allocate(grad(3,ncelet))

  iccocg = 1
  inc = 1
  nswrgp = vcopt_k%nswrgr
  epsrgp = vcopt_k%epsrgr
  imligp = vcopt_k%imligr
  iwarnp = vcopt_k%iwarni
  climgp = vcopt_k%climgr
  extrap = vcopt_k%extrag

  ! Dirichlet boundary condition on the gradient of rho
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

  f_id0 = -1

  call gradient_s                                                 &
 ( f_id0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , epsrgp , climgp , extrap ,                            &
   cromo  , bromo  , viscb  ,                                     &
   grad   )

  ! Production term due to buoyancy
  !   smbrk = P+G
  !   smbre = P+(1-ce3)*G
  if (iscalt.gt.0.and.nscal.ge.iscalt) then
    call field_get_key_double(ivarfl(isca(iscalt)), ksigmas, turb_schmidt)
    prdtur = turb_schmidt
  else
    prdtur = 1.d0
  endif

  ! smbr* store mu_TxS**2
  do iel = 1, ncel
    rho   = cromo(iel)
    visct = cpro_pcvto(iel)
    xeps = cvara_ep(iel)
    xk   = cvara_k(iel)
    ttke = xk / xeps

    gravke = -(grad(1,iel)*gx + grad(2,iel)*gy + grad(3,iel)*gz) &
           / (rho*prdtur)

    ! Implicit Buoyant terms when negativ
    tinstk(iel) = tinstk(iel) + max(-rho*cell_f_vol(iel)*cmu*ttke*gravke, 0.d0)

    ! Explicit Buoyant terms
    smbre(iel) = smbre(iel) + visct*max(gravke, zero)
    smbrk(iel) = smbrk(iel) + visct*gravke
  enddo

  ! Free memory
  deallocate(grad)

endif

! In v2f, we store the production in prdv2f which will be complete further
! for containing the complete production term
if (itytur.eq.5) then
  do iel = 1, ncel
    prdv2f(iel) = smbrk(iel)
  enddo
endif

!===============================================================================
! 7. pre Only for the bl-v2/k model, calculation of E and Ceps2*

!      The terms are stored in          w10, w11
!      Work arrays                      w2, w3
!                                       viscf, viscb
!      Going out of the step we keep w10, w11
!===============================================================================

if (iturb.eq.51) then

  ! Calculation of Ceps2*: it is stored in w10

  do iel=1,ncel
    visct = cpro_pcvto(iel)
    rho   = cromo(iel)
    w3(iel) = visct/rho/sigmak
  enddo

  call viscfa(imvisf, w3, viscf, viscb)

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

  nswrgp = vcopt_k%nswrgr
  imligp = vcopt_k%imligr
  iwarnp = vcopt_k%iwarni
  epsrgp = vcopt_k%epsrgr
  climgp = vcopt_k%climgr
  extrap = vcopt_k%extrag
  iphydp = 0

  call itrgrp &
( ivarfl(ivar), init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydp , &
  iwarnp ,                                                                     &
  epsrgp , climgp , extrap ,                                                   &
  rvoid  ,                                                                     &
  cvara_k         ,                                                            &
  coefap , coefbp , cofafp , cofbfp ,                                          &
  viscf  , viscb  ,                                                            &
  w3     ,                                                                     &
  w10    )

  do iel = 1, ncel
    w10(iel) = -w10(iel)/volume(iel)/cvara_ep(iel)
    w10(iel) = tanh(abs(w10(iel))**1.5d0)
    w10(iel) = cpale2*(1.d0-(cpale2-cpale4)/cpale2*w10(iel)*cvara_al(iel)**3)
  enddo

  ! Calculation of 2*Ceps3*(1-alpha)^3*nu*nut/eps*d2Ui/dxkdxj*d2Ui/dxkdxj:
  !  (i.e. E term / k)           : it is stored in w11

  ! Allocate a work array
  allocate(w12(ncelet))

  call tsepls(w12)

  do iel = 1, ncel

    rho   = cromo(iel)
    xnu   = cpro_pcvlo(iel)/rho
    xnut  = cpro_pcvto(iel)/rho
    xeps = cvara_ep(iel)
    xk   = cvara_k(iel)
    xphi = cvara_phi(iel)

    ttke = xk/xeps
    ttmin = cpalct*sqrt(xnu/xeps)
    tt = sqrt(ttke**2 + ttmin**2)

    w11(iel) = 2.d0*xnu*xnut*w12(iel)*cpale3/xeps                  &
                *(1.d0-cvara_al(iel))**3

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

! If we extrapolate the source terms and rho, we need here rho^n
!                                    and visct, we need here visct^n

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

    smbrk(iel) = cell_f_vol(iel)*                                     &
               ( smbrk(iel) - rho*cvara_ep(iel)                       &
               - d2s3*rho*cvara_k(iel)*divu(iel)                      &
               )

    smbre(iel) = cell_f_vol(iel)*                                              &
               ( cvara_ep(iel)/cvara_k(iel)*( ce1*smbre(iel)                       &
                                            - ce2rc(iel)*rho*cvara_ep(iel)     &
                                            )                              &
               - d2s3*rho*ce1*cvara_ep(iel)*divu(iel)                          &
               )

  enddo

  ! If the solving of k-epsilon is uncoupled, negative source terms are implicited
  if (ikecou.eq.0) then
    do iel = 1, ncel
      xeps = cvara_ep(iel)
      xk   = cvara_k(iel)
      rho = crom(iel)
      ttke = xk / xeps
      tinstk(iel) = tinstk(iel) + rho*cell_f_vol(iel)/ttke            &
                  + max(d2s3*rho*cell_f_vol(iel)*divu(iel), 0.d0)
      tinste(iel) = tinste(iel) + ce2rc(iel)*rho*cell_f_vol(iel)/ttke &
                  + max(d2s3*ce1*rho*cell_f_vol(iel)*divu(iel), 0.d0)
    enddo
  endif

else if (iturb.eq.50) then

  do iel = 1, ncel

    visct = cpro_pcvto(iel)
    rho  = cromo(iel)
    xnu  = cpro_pcvlo(iel)/rho
    xeps = cvara_ep(iel)
    xk   = cvara_k(iel)
    xphi = cvara_phi(iel)
    xphi = max(xphi, epzero)
    ceps1= 1.4d0*(1.d0+cv2fa1*sqrt(1.d0/xphi))
    ttke = xk / xeps
    ttmin = cv2fct*sqrt(xnu/xeps)
    tt = max(ttke, ttmin)

    ! Explicit part
    smbrk(iel) = cell_f_vol(iel)*                                     &
               ( smbrk(iel) - rho*cvara_ep(iel)                       &
               - d2s3*rho*cvara_k(iel)*divu(iel)                      &
               )

    smbre(iel) = cell_f_vol(iel)*                                       &
               ( 1.d0/tt*(ceps1*smbre(iel) - ce2rc(iel)*rho*xeps)   &
               - d2s3*rho*ceps1*xk*divu(iel)                        &
               )

    ! We store the part with Pk in prdv2f which will be reused in resv2f
    prdv2f(iel) = prdv2f(iel)                               &
                - d2s3*rho*cvara_k(iel)*divu(iel)!FIXME this term should be removed

    ! Implicit part
    tinstk(iel) = tinstk(iel) + rho*cell_f_vol(iel)/max(ttke, epzero * ttmin)

    tinstk(iel) = tinstk(iel) + max(d2s3*rho*cell_f_vol(iel)*divu(iel), 0.d0)
    tinste(iel) = tinste(iel) + ce2rc(iel)*rho*cell_f_vol(iel)/tt                &
                + max(d2s3*ceps1*ttke/tt*rho*cell_f_vol(iel)*divu(iel), 0.d0)

  enddo

else if (iturb.eq.51) then

  do iel=1,ncel

    visct = cpro_pcvto(iel)
    rho   = cromo(iel)
    xnu  = cpro_pcvlo(iel)/rho
    xeps = cvara_ep(iel)
    xk   = cvara_k(iel)
    xphi = cvara_phi(iel)
    ttke = xk / xeps
    ttmin = cpalct*sqrt(xnu/xeps)
    tt = sqrt(ttke**2.d0+ttmin**2.d0)

    ! Explicit part
    smbrk(iel) = cell_f_vol(iel)*                                     &
               ( smbrk(iel)                                       &
               - rho*xeps                                         &
               - rho*w11(iel)*xk                                  &
               - d2s3*rho*xk*divu(iel)                            &
               )

    smbre(iel) = cell_f_vol(iel)*                                               &
               ( 1.d0/tt*(cpale1*smbre(iel) - w10(iel)*rho*xeps)            &
               - d2s3*rho*cpale1*xk/tt*divu(iel)                            &
               )

    ! We store the part with Pk in prdv2f which will be reused in resv2f
    prdv2f(iel) = prdv2f(iel)                               &
                - d2s3*rho*cvara_k(iel)*divu(iel)!FIXME this term should be removed

    ! Implicit part
    tinstk(iel) = tinstk(iel) + rho*cell_f_vol(iel)/max(ttke, epzero * ttmin)

    tinstk(iel) = tinstk(iel) + max(d2s3*rho*cell_f_vol(iel)*divu(iel), 0.d0)
    ! Note that w11 is positiv
    tinstk(iel) = tinstk(iel) + w11(iel)*rho*cell_f_vol(iel)
    ! Note that w10 is positiv
    tinste(iel) = tinste(iel) + w10(iel)*rho*cell_f_vol(iel)/tt                  &
                + max(d2s3*cpale1*ttke/tt*rho*cell_f_vol(iel)*divu(iel), 0.d0)

  enddo

endif

!===============================================================================
! 9. Take user source terms into account

!    The scalar strain rate (strain) and the trace of the velocity gradient
!     (divu) are available.
!
!    The part to be explicit is stored in       w7, w8
!    The part to be implicit is stored in       usimpk, usimpe
!    Going out of the step we keep              strain, divu,
!===============================================================================

do iel = 1, ncel
  usimpk(iel) = 0.d0
  usimpe(iel) = 0.d0
  w7(iel) = 0.d0
  w8(iel) = 0.d0
enddo

call cs_user_turbulence_source_terms &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivarfl(ik)      ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   w7     , usimpk )

call cs_user_turbulence_source_terms &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivarfl(iep)     ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   w8     , usimpe )

if (ippmod(iatmos).ge.0) then
  ! Nudging towards optimal interpolation for k
  call field_get_key_int(ivarfl(ik), kopint, f_oi_id)
  if (f_oi_id.ge.0) then
    call cs_at_data_assim_source_term(ivarfl(ik), w7, usimpk)
  endif

  ! Nudging towards optimal interpolation for epsilon
  call field_get_key_int(ivarfl(iep), kopint, f_oi_id)
  if (f_oi_id.ge.0) then
    call cs_at_data_assim_source_term(ivarfl(ik), w8, usimpe)
  endif
endif

! If source terms are extrapolated over time
if (istprv.ge.0) then

  thetak = vcopt_k%thetav
  thetae = vcopt_e%thetav

  do iel = 1, ncel

    ! Recover the value at time (n-1)
    tuexpk = c_st_k_p(iel)
    tuexpe = c_st_eps_p(iel)
    ! Save the values for the next time-step
    c_st_k_p(iel) = smbrk(iel) + w7(iel)
    c_st_eps_p(iel) = smbre(iel) + w8(iel)

    ! Explicit Part
    smbrk(iel) = - thets*tuexpk
    smbre(iel) = - thets*tuexpe
    ! It is assumed that (-usimpk > 0) and though this term is implicit
    smbrk(iel) = usimpk(iel)*cvara_k(iel) + smbrk(iel)
    smbre(iel) = usimpe(iel)*cvara_ep(iel) + smbre(iel)

    ! Implicit part
    tinstk(iel) = tinstk(iel) - usimpk(iel)*thetak
    tinste(iel) = tinste(iel) - usimpe(iel)*thetae

  enddo

! If no extrapolation over time
else
  do iel = 1, ncel
    ! Explicit part
    smbrk(iel) = smbrk(iel) + usimpk(iel)*cvara_k(iel) + w7(iel)
    smbre(iel) = smbre(iel) + usimpe(iel)*cvara_ep(iel) + w8(iel)

    ! Implicit part
    tinstk(iel) = tinstk(iel) + max(-usimpk(iel),zero)
    tinste(iel) = tinste(iel) + max(-usimpe(iel),zero)
  enddo
endif

!===============================================================================
! 10. Taking into account the lagrangian source terms
!     output coupling
!===============================================================================

! Second order is not taken into account
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  do iel = 1,ncel

    ! Explicit and implicit source terms on k
    smbrk(iel)  = smbrk(iel) + tslagr(iel,itske)

    ! Explicit source terms on Eps
    smbre(iel)  = smbre(iel)                                    &
                + ce4 *tslagr(iel,itske) *cvara_ep(iel)             &
                                         /cvara_k(iel)

    ! Implicit source terms on k
    tinstk(iel) = tinstk(iel) + max(-tslagr(iel,itsli),zero)

    ! Implicit souce terms on Eps
    tinste(iel) = tinste(iel) + max((-ce4*tslagr(iel,itske)/cvara_k(iel)), zero)

  enddo

endif

!===============================================================================
! 11. Mass source terms (Implicit and explicit parts)

!       Going out of the step we keep divu,
!                                     smbrk, smbre
!===============================================================================

if (ncesmp.gt.0) then

  do iel = 1, ncel
    w2(iel) = 0.d0
    w3(iel) = 0.d0
  enddo

  ! Integer equal to 1 (for navsto: nb of sur-iter)
  iiun = 1

  ! We incremente smbrs with -Gamma.var_prev and rovsdt with Gamma
  ivar = ik

  call catsma &
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
   isto2t ,                                                       &
   icetsm , itypsm(:,ivar)  ,                                     &
   volume , cvara_k         , smacel(:,ivar) , smacel(1,ipr) ,    &
   smbrk  , w2     , w4 )

  ivar = iep

  call catsma &
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
   isto2t ,                                                       &
   icetsm , itypsm(:,ivar)  ,                                     &
   volume , cvara_ep        , smacel(:,ivar) , smacel(1,ipr) ,    &
   smbre  , w3     , w5 )

  ! If we extrapolate the source terms we put Gamma Pinj in c_st
  if (istprv.ge.0) then
    do iel = 1, ncel
      c_st_k_p(iel)   = c_st_k_p(iel)   + w4(iel)
      c_st_eps_p(iel) = c_st_eps_p(iel) + w5(iel)
    enddo
  ! Otherwise we put it directly in smbr
  else
    do iel = 1, ncel
      smbrk(iel) = smbrk(iel) + w4(iel)
      smbre(iel) = smbre(iel) + w5(iel)
    enddo
  endif

  ! Implicit part
  do iel = 1, ncel
    tinstk(iel) = tinstk(iel) + w2(iel)
    tinste(iel) = tinste(iel) + w3(iel)
  enddo

endif

!===============================================================================
! 12.1 Taking into account the terms of conv/diff in the second member for the
!      strengthened coupling k-epsilon (ikecou == 1)

!      Work table                       w4, w5
!      The terms are stored in          w7 et w8, then added to smbrk, smbre
!      Going out of the step we keep    divu,
!                                       smbrk, smbre
!                                       w7, w8
!===============================================================================

if (ikecou.eq.1) then

  do iel = 1, ncel
    w7 (iel) = 0.d0
    w8 (iel) = 0.d0
  enddo

  ! ---> Treatment of k
  ivar   = ik
  call field_get_label(ivarfl(ivar), chaine)

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)
  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  if (vcopt_k%idiff.ge. 1) then

    do iel = 1, ncel
      w4(iel) = viscl(iel) + vcopt_k%idifft*cvisct(iel)/sigmak
    enddo

    call viscfa(imvisf, w4, viscf, viscb)

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
  iconvp = vcopt_k%iconv
  imasac = 1
  idiffp = vcopt_k%idiff
  nswrgp = vcopt_k%nswrgr
  imligp = vcopt_k%imligr
  ircflp = vcopt_k%ircflu
  ischcp = vcopt_k%ischcv
  isstpp = vcopt_k%isstpc
  iwarnp = vcopt_k%iwarni
  imucpp = 0
  idftnp = ISOTROPIC_DIFFUSION
  blencp = vcopt_k%blencv
  epsrgp = vcopt_k%epsrgr
  climgp = vcopt_k%climgr
  extrap = vcopt_k%extrag
  relaxp = vcopt_k%relaxv
  thetap = vcopt_k%thetav
  ! all boundary convective flux with upwind
  icvflb = 0

  call bilsca &
 ( idtvar , ivarfl(ivar)    , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                            &
   iwarnp , imucpp , idftnp , imasac ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,                   &
   cvara_k         , cvara_k         ,                                     &
   coefap , coefbp , cofafp , cofbfp ,                                     &
   imasfl , bmasfl ,                                                       &
   viscf  , viscb  , rvoid  , rvoid  ,                                     &
   rvoid  , rvoid  ,                                                       &
   icvflb , ivoid  ,                                                       &
   w7     )

  if (vcopt_k%iwarni.ge.2) then
    rnorm = sqrt(cs_gdot(ncel,smbrk,smbrk))
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif


  ! ---> Treatment of epsilon
  ivar   = iep
  call field_get_label(ivarfl(ivar), chaine)

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)
  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  if (vcopt_e%idiff.ge. 1) then
    do iel = 1, ncel
      w4(iel) = viscl(iel)                                        &
              + vcopt_e%idifft*cvisct(iel)/sigmae
    enddo

    call viscfa(imvisf, w4, viscf, viscb)

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
  iconvp = vcopt_e%iconv
  imasac = 1
  idiffp = vcopt_e%idiff
  nswrgp = vcopt_e%nswrgr
  imligp = vcopt_e%imligr
  ircflp = vcopt_e%ircflu
  ischcp = vcopt_e%ischcv
  isstpp = vcopt_e%isstpc
  iwarnp = vcopt_e%iwarni
  imucpp = 0
  idftnp = ISOTROPIC_DIFFUSION
  blencp = vcopt_e%blencv
  epsrgp = vcopt_e%epsrgr
  climgp = vcopt_e%climgr
  extrap = vcopt_e%extrag
  relaxp = vcopt_e%relaxv
  thetap = vcopt_e%thetav
  ! all boundary convective flux with upwind
  icvflb = 0

  call bilsca &
 ( idtvar , ivarfl(ivar)    , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                            &
   iwarnp , imucpp , idftnp , imasac ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,                   &
   cvara_ep        , cvara_ep        ,                                     &
   coefap , coefbp , cofafp , cofbfp ,                                     &
   imasfl , bmasfl ,                                                       &
   viscf  , viscb  , rvoid  , rvoid  ,                                     &
   rvoid  , rvoid  ,                                                       &
   icvflb , ivoid  ,                                                       &
   w8     )

  if (vcopt_e%iwarni.ge.2) then
    rnorm = sqrt(cs_gdot(ncel,smbre,smbre))
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

! Second order is not taken into account
if (ikecou.eq.1) then

  if (iturb.eq.20) then

    do iel = 1, ncel

      rho = crom(iel)
      visct = cpro_pcvto(iel)

      ! Coupled solving
      romvsd = 1.d0/(rho*volume(iel))
      smbrk(iel)=smbrk(iel)*romvsd
      smbre(iel)=smbre(iel)*romvsd
      divp23= d2s3*max(divu(iel),zero)

      epssuk = cvara_ep(iel)/cvara_k(iel)

      a11 = 1.d0/dt(iel)                                          &
           -2.d0*cvara_k(iel)/cvara_ep(iel)                               &
           *cmu*min(prdtke(iel)/visct,zero)+divp23
      a12 = 1.d0
      a21 = -ce1*cmu*prdeps(iel)/visct-ce2rc(iel)*epssuk*epssuk
      a22 = 1.d0/dt(iel)+ce1*divp23                               &
           +2.d0*ce2rc(iel)*epssuk

      unsdet = 1.d0/(a11*a22 -a12*a21)

      deltk = ( a22*smbrk(iel) -a12*smbre(iel) )*unsdet
      delte = (-a21*smbrk(iel) +a11*smbre(iel) )*unsdet

      ! New source term for the iterative process
      romvsd = rho*cell_f_vol(iel)/dt(iel)

      smbrk(iel) = romvsd*deltk
      smbre(iel) = romvsd*delte

    enddo

    ! we remove the convection/diffusion at time n from smbrk and smbre
    ! if they were calculated
    do iel = 1, ncel
      smbrk(iel) = smbrk(iel) - w7(iel)
      smbre(iel) = smbre(iel) - w8(iel)
    enddo

  ! In verini we block the mix iturb!=20/ikecou=1
  else

    write(nfecra,*)'ikecou=1 is not validated with this turbulent model'
    call csexit (1)

  endif

endif

!===============================================================================
! 13. Finalization of the Right Hand Side when activating 2nd time order
!===============================================================================

if (istprv.ge.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + thetp1 * c_st_k_p(iel)
    smbre(iel) = smbre(iel) + thetp1 * c_st_eps_p(iel)
  enddo
endif

!===============================================================================
! 14. Solving

!       We use                          smbrk, smbre,  tinstk, tinste
!       Work table                      w1
!===============================================================================

! ---> turbulent kinetic (k) energy treatment
ivar = ik

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

! Face viscosity
if (vcopt_k%idiff.ge.1) then

  do iel = 1, ncel
    if (iturb.eq.51) then
      w1(iel) = viscl(iel)/2.d0                                   &
              + vcopt_k%idifft*cvisct(iel)/sigmak
    else
      w1(iel) = viscl(iel)                                        &
              + vcopt_k%idifft*cvisct(iel)/sigmak
    endif
  enddo

  call viscfa(imvisf, w1, viscf, viscb)

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

! Solving k
iconvp = vcopt_k%iconv
idiffp = vcopt_k%idiff
ndircp = vcopt_k%ndircl
nswrsp = vcopt_k%nswrsm
nswrgp = vcopt_k%nswrgr
imligp = vcopt_k%imligr
ircflp = vcopt_k%ircflu
ischcp = vcopt_k%ischcv
isstpp = vcopt_k%isstpc
iescap = 0
imucpp = 0
idftnp = vcopt_k%idften
iswdyp = vcopt_k%iswdyn
iwarnp = vcopt_k%iwarni
blencp = vcopt_k%blencv
epsilp = vcopt_k%epsilo
epsrsp = vcopt_k%epsrsm
epsrgp = vcopt_k%epsrgr
climgp = vcopt_k%climgr
extrap = vcopt_k%extrag
relaxp = vcopt_k%relaxv
thetap = vcopt_k%thetav
! all boundary convective flux with upwind
icvflb = 0
normp = -1.d0
init   = 1

call codits &
 ( idtvar , init   , ivarfl(ivar)    , iconvp , idiffp , ndircp , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp , normp  ,                                              &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   cvara_k         , cvara_k         ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   tinstk , smbrk  , cvar_k          , dpvar  ,                   &
   rvoid  , rvoid  )

! ---> Turbulent dissipation (epsilon) treatment
ivar = iep

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

! Face viscosity
if (vcopt_e%idiff.ge.1) then
  do iel = 1, ncel
    if (iturb.eq.51) then
      w1(iel) = viscl(iel)/2.d0                                   &
              + vcopt_e%idifft*cvisct(iel)/cpalse
    else
      w1(iel) = viscl(iel)                                        &
              + vcopt_e%idifft*cvisct(iel)/sigmae
    endif
  enddo

  call viscfa(imvisf, w1, viscf, viscb)

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

! Solving epsilon
iconvp = vcopt_e%iconv
idiffp = vcopt_e%idiff
ndircp = vcopt_e%ndircl
nswrsp = vcopt_e%nswrsm
nswrgp = vcopt_e%nswrgr
imligp = vcopt_e%imligr
ircflp = vcopt_e%ircflu
ischcp = vcopt_e%ischcv
isstpp = vcopt_e%isstpc
iescap = 0
imucpp = 0
idftnp = vcopt_e%idften
iswdyp = vcopt_e%iswdyn
iwarnp = vcopt_e%iwarni
blencp = vcopt_e%blencv
epsilp = vcopt_e%epsilo
epsrsp = vcopt_e%epsrsm
epsrgp = vcopt_e%epsrgr
climgp = vcopt_e%climgr
extrap = vcopt_e%extrag
relaxp = vcopt_e%relaxv
thetap = vcopt_e%thetav
! all boundary convective flux with upwind
icvflb = 0
normp = -1.d0

call codits &
 ( idtvar , init   , ivarfl(ivar)    , iconvp , idiffp , ndircp , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp , normp  ,                                              &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   cvara_ep        , cvara_ep        ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   tinste , smbre  , cvar_ep         , dpvar  ,                   &
   rvoid  , rvoid  )

!===============================================================================
! 15. Clipping
!===============================================================================

iclip = 1
call clipke(ncelet, ncel, iclip)

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
'   ** Resolution du k-epsilon a prod lineaire',/,&
'      ---------------------------------------   ',/)
 1002 format(/,                                      &
'   ** Resolution du v2f (k et epsilon)               ',/,&
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
