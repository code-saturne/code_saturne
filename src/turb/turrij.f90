!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file turrij.f90
!>
!> \brief Solving the \f$ R_{ij} - \epsilon \f$ for incompressible flows or
!>  slightly compressible flows for one time step.
!>
!> Please refer to the
!> <a href="../../theory.pdf#rijeps"><b>\f$ R_{ij} - \epsilon \f$ model</b></a>
!> section of the theory guide for more informations, as well as the
!> <a href="../../theory.pdf#turrij"><b>turrij</b></a> section.
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
!> \param[in]     dt            time step (per cell)
!> \param[in]     tslagr        coupling term of the lagangian module
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        values of the variables associated to the
!>                               mass source
!>                               (for ivar=ipr, smacel is the mass flux)
!_______________________________________________________________________________

subroutine turrij &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   tslagr ,                                                       &
   ckupdc , smacel )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use lagran
use numvar
use entsor
use cstphy
use cstnum
use optcal
use ppincl
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp)

integer, dimension(ncesmp,nvar), target :: itypsm

double precision dt(ncelet)
double precision ckupdc(6,ncepdp)

double precision, dimension(ncesmp,nvar), target ::  smacel
double precision, dimension(ncelet,ntersl), target :: tslagr
double precision, dimension(:), pointer :: bromo, cromo

! Local variables

integer          ifac  , iel   , ivar  , isou, jsou
integer          iflmas, iflmab
integer          inc   , iccocg
integer          iwarnp, iclip
integer          imrgrp, nswrgp, imligp
integer          f_id0 , f_id, st_prv_id
integer          iprev
integer          key_t_ext_id
integer          iroext
integer          f_id_phij
integer          icvflb
integer          ivoid(1)
double precision epsrgp, climgp
double precision rhothe
double precision utaurf, ut2, ypa, ya, tke, xunorm, limiter, nu0, alpha
double precision xnoral, xnal(3)
double precision k, p, thets, thetv, tuexpr, thetp1
double precision d1s3, d2s3

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbr, rovsdt
double precision, allocatable, dimension(:,:,:), target :: gradv
double precision, allocatable, dimension(:,:), target :: produc
double precision, allocatable, dimension(:,:) :: gradro
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:,:) :: smbrts, gatinj
double precision, allocatable, dimension(:,:,:) :: rovsdtts
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb

double precision, pointer, dimension(:) :: tslagi
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:), pointer :: cvara_scalt
double precision, dimension(:), pointer :: cvar_ep, cvar_al
double precision, dimension(:,:), pointer :: cvara_rij, cvar_rij, vel
double precision, dimension(:,:), pointer :: c_st_prv, lagr_st_rij
double precision, dimension(:,:), pointer :: cpro_produc
double precision, dimension(:,:,:), pointer :: cpro_gradv
double precision, dimension(:,:), pointer :: cpro_press_correl
double precision, dimension(:), pointer :: cpro_beta

double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:,:), pointer :: coefa_rij, cofaf_rij
double precision, dimension(:,:,:), pointer :: coefb_rij, cofbf_rij

type(var_cal_opt) :: vcopt
type(var_cal_opt), target :: vcopt_loc
type(var_cal_opt), pointer :: p_k_value
type(c_ptr) :: c_k_value

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbr(ncelet), rovsdt(ncelet))

call field_get_id_try("velocity_gradient", f_id)
if (f_id.ge.0) then
  call field_get_val_t(f_id, cpro_gradv)
else
  allocate(gradv(3, 3, ncelet))
  cpro_gradv => gradv
endif

allocate(smbrts(6,ncelet))
allocate(rovsdtts(6,6,ncelet))

! Allocate other arrays, depending on user options
call field_get_id_try("rij_production", f_id)
if (f_id.ge.0) then
  call field_get_val_v(f_id, cpro_produc)
else
  allocate(produc(6,ncelet))
  cpro_produc => produc
endif

call field_get_id_try("rij_pressure_strain_correlation", f_id_phij)
if (f_id_phij.ge.0) then
  call field_get_val_v(f_id_phij, cpro_press_correl)
endif
call field_get_val_s(ivarfl(iep), cvar_ep)

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

call field_get_val_v(ivarfl(irij), cvar_rij)
call field_get_val_prev_v(ivarfl(irij), cvara_rij)
call field_get_key_struct_var_cal_opt(ivarfl(irij), vcopt)

thets  = thetst
thetv  = vcopt%thetav

if (vcopt%iwarni.ge.1) then
  if (iturb.eq.30) then
    write(nfecra,1000)
  elseif (iturb.eq.31) then
    write(nfecra,1001)
  else
    write(nfecra,1002)
  endif
endif

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

call field_get_key_int(ivarfl(irij), kstprv, st_prv_id)
if (st_prv_id .ge. 0) then
  call field_get_val_v(st_prv_id, c_st_prv)
else
  c_st_prv=> null()
endif

! Mass flux
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

!===============================================================================
! 1.1 Advanced init for EBRSM
!===============================================================================

! Automatic reinitialization at the end of the first iteration:
! wall distance y^+ is computed with -C log(1-alpha), where C=CL*Ceta*L*kappa,
! then y so we have an idea of the wall distance in complex geometries.
! Then U is initialized with a Reichard layer,
! Epsilon by 1/(kappa y), clipped next to the wall at its value for y^+=15
! k is given by a blending between eps/(2 nu)*y^2 and utau/sqrt(Cmu).
! The blending function is chosen so that the asymptotic behavior
! and give the correct peak of k.

!TODO FIXME Are the BC uncompatible?
if (ntcabs.eq.1.and.reinit_turb.eq.1.and.iturb.eq.32) then

  allocate(grad(3,ncelet))

  ! Compute the gradient of Alpha
  iprev  = 0
  inc    = 1
  iccocg = 1

  call field_gradient_scalar(ivarfl(ial), iprev, 0, inc, iccocg, grad)

  call field_get_val_v(ivarfl(irij), cvar_rij)
  utaurf=0.05d0*uref

  call field_get_val_s(ivarfl(iep), cvar_ep)
  call field_get_val_s(ivarfl(ial), cvar_al)
  call field_get_val_v(ivarfl(iu), vel)

  nu0 = viscl0/ro0

  do iel = 1, ncel
    ! Compute the velocity magnitude
    xunorm = vel(1,iel)**2 + vel(2,iel)**2 + vel(3,iel)**2
    xunorm = sqrt(xunorm)

    ! y+ is bounded by 400, because in the Reichard profile,
    ! it corresponds to saturation (u>uref)
    cvar_al(iel) = max(min(cvar_al(iel),(1.d0-exp(-400.d0/50.d0))),0.d0)
    ! Compute the magnitude of the Alpha gradient
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
    ! Velocity magnitude is imposed (limitted only), the direction is
    ! conserved
    if (xunorm.le.1.d-12*uref) then
      limiter = 1.d0
    else
      limiter = min( utaurf/xunorm*(2.5d0*dlog(1.d0+0.4d0*ypa) &
                    +7.8d0*(1.d0-dexp(-ypa/11.d0)              &
                    -(ypa/11.d0)*dexp(-0.33d0*ypa))),          &
                    1.d0)
    endif

    vel(1,iel) = limiter*vel(1,iel)
    vel(2,iel) = limiter*vel(2,iel)
    vel(3,iel) = limiter*vel(3,iel)

    ut2 = 0.05d0*uref
    cvar_ep(iel) = utaurf**3*min(1.d0/(xkappa*15.d0*nu0/utaurf), &
                                 1.d0/(xkappa*ya))
    tke =  cvar_ep(iel)/2.d0/nu0*ya**2                  &
           * exp(-ypa/25.d0)**2                         &
          + ut2**2/0.3d0*(1.d0-exp(-ypa/25.d0))**2

    cvar_rij(1,iel) = alpha**3       *2.d0/3.d0        *tke &
                    + (1.d0-alpha**3)*(1.d0-xnal(1)**2)*tke
    cvar_rij(2,iel) = alpha**3       *2.d0/3.d0        *tke &
                    + (1.d0-alpha**3)*(1.d0-xnal(2)**2)*tke
    cvar_rij(3,iel) = alpha**3       *2.d0/3.d0        *tke &
                    + (1.d0-alpha**3)*(1.d0-xnal(3)**2)*tke
    cvar_rij(4,iel) = -(1.d0-alpha**3)*(xnal(1)*xnal(2))*tke
    cvar_rij(5,iel) = -(1.d0-alpha**3)*(xnal(2)*xnal(3))*tke
    cvar_rij(6,iel) = -(1.d0-alpha**3)*(xnal(1)*xnal(3))*tke
  enddo

  call field_current_to_previous(ivarfl(iu))
  call field_current_to_previous(ivarfl(irij))

  deallocate(grad)
endif

!===============================================================================
! 2.1 Compute the velocity gradient
!===============================================================================

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, 0, inc, cpro_gradv)

!===============================================================================
! 2.2 Compute the production term for Rij
!===============================================================================

do iel = 1, ncel

   ! Pij = - (Rik dUk/dXj + dUk/dXi Rkj)
   ! Pij is stored as (P11, P22, P33, P12, P23, P13)
   cpro_produc(1,iel) = &
                  - 2.0d0*(cvara_rij(1,iel)*cpro_gradv(1, 1, iel) +           &
                           cvara_rij(4,iel)*cpro_gradv(2, 1, iel) +           &
                           cvara_rij(6,iel)*cpro_gradv(3, 1, iel) )

   cpro_produc(4,iel) = &
                  - (cvara_rij(4,iel)*cpro_gradv(1, 1, iel) +                 &
                     cvara_rij(2,iel)*cpro_gradv(2, 1, iel) +                 &
                     cvara_rij(5,iel)*cpro_gradv(3, 1, iel) )                 &
                  - (cvara_rij(1,iel)*cpro_gradv(1, 2, iel) +                 &
                     cvara_rij(4,iel)*cpro_gradv(2, 2, iel) +                 &
                     cvara_rij(6,iel)*cpro_gradv(3, 2, iel) )

   cpro_produc(6,iel) = &
                  - (cvara_rij(6,iel)*cpro_gradv(1, 1, iel) +                 &
                     cvara_rij(5,iel)*cpro_gradv(2, 1, iel) +                 &
                     cvara_rij(3,iel)*cpro_gradv(3, 1, iel) )                 &
                  - (cvara_rij(1,iel)*cpro_gradv(1, 3, iel) +                 &
                     cvara_rij(4,iel)*cpro_gradv(2, 3, iel) +                 &
                     cvara_rij(6,iel)*cpro_gradv(3, 3, iel) )

   cpro_produc(2,iel) = &
                  - 2.0d0*(cvara_rij(4,iel)*cpro_gradv(1, 2, iel) +           &
                           cvara_rij(2,iel)*cpro_gradv(2, 2, iel) +           &
                           cvara_rij(5,iel)*cpro_gradv(3, 2, iel) )

   cpro_produc(5,iel) = &
                  - (cvara_rij(6,iel)*cpro_gradv(1, 2, iel) +                 &
                     cvara_rij(5,iel)*cpro_gradv(2, 2, iel) +                 &
                     cvara_rij(3,iel)*cpro_gradv(3, 2, iel) )                 &
                  - (cvara_rij(4,iel)*cpro_gradv(1, 3, iel) +                 &
                     cvara_rij(2,iel)*cpro_gradv(2, 3, iel) +                 &
                     cvara_rij(5,iel)*cpro_gradv(3, 3, iel) )

   cpro_produc(3,iel) = &
                  - 2.0d0*(cvara_rij(6,iel)*cpro_gradv(1, 3, iel) +           &
                           cvara_rij(5,iel)*cpro_gradv(2, 3, iel) +           &
                           cvara_rij(3,iel)*cpro_gradv(3, 3, iel) )

enddo

!===============================================================================
! 2.3 Compute the pressure correlation  term for Rij
!===============================================================================
! Phi,ij = Phi1,ij+Phi2,ij
! Phi,ij = -C1 k/eps (Rij-2/3k dij) - C2 (Pij-2/3P dij)
! Phi,ij is stored as (Phi11, Phi22, Phi33, Phi12, Phi23, Phi13)

! TODO : coherent with the model
if (f_id_phij.ge.0) then
  d1s3 = 1.0d0/3.0d0
  d2s3 = 2.0d0/3.0d0
  do iel = 1, ncel
    k=0.5*(cvara_rij(1,iel)+cvara_rij(2,iel)+cvara_rij(3,iel))
    p=0.5*(cpro_produc(1,iel)+cpro_produc(2,iel)+cpro_produc(3,iel))
    do isou=1,3
      cpro_press_correl(isou, iel) = -crij1*cvar_ep(iel)/k*(cvara_rij(isou,iel)-d2s3*k)  &
                                     -crij2*(cpro_produc(isou,iel)-d2s3*p)
    enddo
    do isou=4,6
      cpro_press_correl(isou, iel) = -crij1*cvar_ep(iel)/k*(cvara_rij(isou,iel))  &
                                     -crij2*(cpro_produc(isou,iel))
    enddo
  enddo
endif

!===============================================================================
! 3. Compute the density gradient for buoyant terms
!===============================================================================

! Note that the buoyant term is normally expressed in temr of
! (u'T') or (u'rho') here modelled with a GGDH:
!   (u'rho') = C * k/eps * R_ij Grad_j(rho)

! Buoyant term for the Atmospheric module
! (function of the potential temperature)
if (igrari.eq.1 .and. ippmod(iatmos).ge.1) then
  ! Allocate a temporary array for the gradient calculation
  ! Warning, grad(theta) here
  allocate(gradro(3,ncelet))

  call field_get_val_s(icrom, cromo)

  call field_get_val_prev_s(ivarfl(isca(iscalt)), cvara_scalt)

  inc = 1
  iccocg = 1

  call field_gradient_scalar(ivarfl(isca(iscalt)), 1, 0, inc, iccocg, gradro)

  ! gradro stores: - rho grad(theta)/theta
  ! grad(rho) and grad(theta) have opposite signs
  do iel = 1, ncel
    rhothe = cromo(iel)/cvara_scalt(iel)
    gradro(1, iel) = -rhothe*gradro(1, iel)
    gradro(2, iel) = -rhothe*gradro(2, iel)
    gradro(3, iel) = -rhothe*gradro(3, iel)
  enddo

else if (igrari.eq.1) then
  ! Allocate a temporary array for the gradient calculation
  allocate(gradro(3,ncelet))

  ! Boussinesq approximation, only for the thermal scalar for the moment
  if (idilat.eq.0)  then

    iccocg = 1
    inc = 1

    ! Use the current value...
    call field_gradient_scalar(ivarfl(isca(iscalt)), 0, 0, inc, iccocg, gradro)

    !FIXME make it dependant on the scalar and use is_buoyant field
    call field_get_val_s(ibeta, cpro_beta)

    ! gradro stores: - beta grad(T)
    ! grad(rho) and grad(T) have opposite signs
    do iel = 1, ncel
      gradro(1, iel) = - ro0 * cpro_beta(iel) * gradro(1, iel)
      gradro(2, iel) = - ro0 * cpro_beta(iel) * gradro(2, iel)
      gradro(3, iel) = - ro0 * cpro_beta(iel) * gradro(3, iel)
    enddo

  else

    ! Boundary conditions: Dirichlet romb
    ! We use viscb to store the relative coefficient of rom
    ! We impose in Dirichlet (coefa) the value romb

    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

    ! The choice below has the advantage to be simple

    imrgrp = vcopt%imrgra
    nswrgp = vcopt%nswrgr
    imligp = vcopt%imligr
    iwarnp = vcopt%iwarni
    epsrgp = vcopt%epsrgr
    climgp = vcopt%climgr

    f_id0 = -1
    iccocg = 1

    call field_get_key_int(icrom, key_t_ext_id, iroext)
    ! If we extrapolate the source terms and rho, we use cpdt rho^n
    if(isto2t.gt.0.and.iroext.gt.0) then
      call field_get_val_prev_s(icrom, cromo)
      call field_get_val_prev_s(ibrom, bromo)
    else
      call field_get_val_s(icrom, cromo)
      call field_get_val_s(ibrom, bromo)
    endif

    call gradient_s(f_id0, imrgrp, inc, iccocg, nswrgp, imligp,      &
                    iwarnp, epsrgp, climgp, cromo, bromo, viscb,     &
                    gradro)

  endif
endif

!===============================================================================
! 4.1 Prepare to solve Rij
!     We solve the equation in a routine similar to covofi.f90
!===============================================================================

!===============================================================================
! 4.1.1 Source terms for Rij
!===============================================================================

do iel = 1, ncel
  do isou = 1 ,6
    smbrts(isou,iel) = 0.d0
  enddo
enddo
do iel = 1, ncel
  do isou = 1, 6
    do jsou = 1, 6
      rovsdtts(isou,jsou,iel) = 0.d0
    enddo
  enddo
enddo

! User source terms
!------------------

call cs_user_turbulence_source_terms2(nvar, nscal, ncepdp, ncesmp,   &
                                      ivarfl(irij),                  &
                                      icepdc, icetsm, itypsm,        &
                                      ckupdc, smacel,                &
                                      smbrts, rovsdtts)

! C version
call user_source_terms(ivarfl(irij), smbrts, rovsdtts)

! If we extrapolate the source terms
if (st_prv_id.ge.0) then
  !     S as Source, V as Variable
  do iel = 1, ncel
    do isou = 1, 6
      ! Save for exchange
      tuexpr = c_st_prv(isou,iel)
      ! For continuation and the next time step
      c_st_prv(isou,iel) = smbrts(isou,iel)

      smbrts(isou,iel) = - thets*tuexpr
      ! Right hand side of the previous time step
      ! We suppose -rovsdt > 0: we implicit
      !    the user source term (the rest)
      do jsou = 1, 6
        smbrts(isou,iel) =   smbrts(isou,iel) &
                           + rovsdtts(jsou,isou,iel)*cvara_rij(jsou,iel)
        ! Diagonal
        rovsdtts(jsou,isou,iel) = - thetv*rovsdtts(jsou,isou,iel)
      enddo
    enddo
  enddo
else
  do iel = 1, ncel
    do isou = 1, 6
      do jsou = 1, 6
        smbrts(isou,iel) =   smbrts(isou,iel) &
                           + rovsdtts(jsou,isou,iel)*cvara_rij(jsou,iel)
      enddo
      rovsdtts(isou,isou,iel) = max(-rovsdtts(isou,isou,iel), 0.d0)
    enddo
  enddo
endif

! Lagrangian source terms
!------------------------

! 2nd order is not taken into account
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  call field_get_val_v_by_name('rij_st_lagr', lagr_st_rij)
  tslagi => tslagr(1:ncelet,itsli)

  do iel = 1,ncel
    do isou = 1, 6
      smbrts(isou, iel) = smbrts(isou, iel) + lagr_st_rij(isou,iel)
      rovsdtts(isou,isou,iel) = rovsdtts(isou,isou, iel) + max(-tslagi(iel),zero)
    enddo
  enddo

endif

! Mass source term
!-----------------

if (ncesmp.gt.0) then

  allocate(gatinj(6,ncelet))

  ! We increment smbr with -Gamma.var_prev and rovsdr with Gamma
  call catsmt(ncesmp, 1, icetsm, itypsm(:,ivar),                            &
              cell_f_vol, cvara_rij, smacel(:,ivar), smacel(:,ipr),  &
              smbrts, rovsdtts, gatinj)

  do isou = 1, 6

    ! If we extrapolate the source terms we put Gamma Pinj in the previous st
    if (st_prv_id.ge.0) then
      do iel = 1, ncel
        c_st_prv(isou,iel) = c_st_prv(isou,iel) + gatinj(isou,iel)
      enddo
    ! Otherwise we put it directly in the RHS
    else
      do iel = 1, ncel
        smbrts(isou, iel) = smbrts(isou, iel) + gatinj(isou,iel)
      enddo
    endif

  enddo

  deallocate(gatinj)

endif

!===============================================================================
! 4.1.2 Unsteady term
!===============================================================================

! ---> Added in the matrix diagonal

if (vcopt%istat .eq. 1) then
  do iel = 1, ncel
    do isou = 1, 6
      rovsdtts(isou,isou,iel) =   rovsdtts(isou,isou,iel)                     &
                                + (crom(iel)/dt(iel))*cell_f_vol(iel)
    enddo
  enddo
endif

ivar = irij

!===============================================================================
! 4.1.3 Rij-epsilon model-specific terms
!===============================================================================

allocate(viscce(6,ncelet))
allocate(weighf(2,nfac))
allocate(weighb(nfabor))

! Rij-epsilon standard (LRR)
if (iturb.eq.30) then

  if (irijco.eq.1) then
    call resrij2(ivar,                                       &
                 cpro_gradv, cpro_produc, gradro,                 &
                 viscf, viscb, viscce,                       &
                 smbrts, rovsdtts,                           &
                 weighf, weighb)
  else
    call resrij(ivar,                                        &
                cpro_produc, gradro,                         &
                viscf, viscb, viscce,                        &
                smbrts, rovsdtts,                            &
                weighf, weighb)
  endif

! Rij-epsilon SSG or EBRSM
elseif (iturb.eq.31.or.iturb.eq.32) then

  call resssg2(ivar,                                         &
               cpro_gradv, cpro_produc, gradro,                   &
               viscf, viscb, viscce,                         &
               smbrts, rovsdtts,                             &
               weighf, weighb)

endif

!===============================================================================
! 4.2 Solve Rij
!===============================================================================

if (st_prv_id.ge.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    do isou = 1, 6
      smbrts(isou,iel) = smbrts(isou,iel) + thetp1*c_st_prv(isou,iel)
    enddo
  enddo
endif

! all boundary convective flux with upwind
icvflb = 0

call field_get_coefa_v(ivarfl(ivar), coefa_rij)
call field_get_coefb_v(ivarfl(ivar), coefb_rij)
call field_get_coefaf_v(ivarfl(ivar), cofaf_rij)
call field_get_coefbf_v(ivarfl(ivar), cofbf_rij)

vcopt_loc = vcopt

vcopt_loc%istat  = -1
vcopt_loc%idifft = -1
vcopt_loc%iwgrec = 0 ! Warning, may be overwritten if a field
vcopt_loc%thetav = thetv
vcopt_loc%blend_st = 0 ! Warning, may be overwritten if a field

p_k_value => vcopt_loc
c_k_value = equation_param_from_vcopt(c_loc(p_k_value))

call cs_equation_iterative_solve_tensor(idtvar, ivarfl(ivar), c_null_char,     &
                                        c_k_value, cvara_rij, cvara_rij,       &
                                        coefa_rij, coefb_rij,                  &
                                        cofaf_rij, cofbf_rij,                  &
                                        imasfl, bmasfl, viscf,                 &
                                        viscb, viscf, viscb, viscce,           &
                                        weighf, weighb, icvflb, ivoid,         &
                                        rovsdtts, smbrts, cvar_rij)

!===============================================================================
! 5. Solve Epsilon
!===============================================================================

call reseps(nvar, ncesmp, icetsm, itypsm,                    &
            dt, cpro_gradv, cpro_produc, gradro,                  &
            smacel, viscf, viscb,                            &
            smbr, rovsdt)

!===============================================================================
! 6. Clipping
!===============================================================================

if (iturb.eq.32) then
  iclip = 1
else
  iclip = 2
endif

if (irijco.eq.1) then
  call clprij2(ncel)
else
  call clprij(ncel, iclip)
endif

! Free memory
deallocate(viscf, viscb)
deallocate(smbr, rovsdt)
if (allocated(gradro)) deallocate(gradro)
if (allocated(produc)) deallocate(produc)
if (allocated(gradv)) deallocate(gradv)

deallocate(smbrts)
deallocate(rovsdtts)

!--------
! Formats
!--------

 1000 format(/,                                      &
'   ** Solving Rij-EPSILON LRR'                   ,/,&
'      -----------------------'                   ,/)
 1001 format(/,                                      &
'   ** Solving Rij-EPSILON SSG'                   ,/,&
'      -----------------------'                   ,/)
 1002 format(/,                                      &
'   ** Solving Rij-EPSILON EBRSM'                 ,/,&
'      -------------------------'                 ,/)

!----
! End
!----

return

end subroutine
