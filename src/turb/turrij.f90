!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
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
!> \param[in]     itypsm        mass source type for the variables (cf. cs_user_mass_source_terms)
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
use pointe, only: rvoid1, rvoid2
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
double precision ckupdc(ncepdp,6)

double precision, dimension(ncesmp,nvar), target ::  smacel
double precision, dimension(ncelet,ntersl), target :: tslagr
double precision, dimension(:), pointer :: bromo, cromo

! Local variables

integer          ifac  , iel   , ivar  , isou  , ii
integer          inc   , iccocg
integer          iwarnp, iclip
integer          nswrgp, imligp
integer          f_id0
integer          iitsla
integer          iprev
double precision epsrgp, climgp, extrap
double precision rhothe
double precision utaurf,ut2,ypa,ya,tke,xunorm, limiter, nu0,alpha
double precision xnoral, xnal(3)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbr, rovsdt
double precision, allocatable, dimension(:,:,:) :: gradv
double precision, allocatable, dimension(:,:) :: produc
double precision, allocatable, dimension(:,:) :: gradro
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:,:) :: smbrts
double precision, allocatable, dimension(:,:,:) ::rovsdtts

double precision, pointer, dimension(:) :: tslage, tslagi
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: cvara_r12, cvara_r13, cvara_r23
double precision, dimension(:), pointer :: cvara_scalt
double precision, dimension(:), pointer :: cvar_ep, cvar_al
double precision, dimension(:,:), pointer :: cvara_rij, cvar_rij, vel
double precision, dimension(:,:), pointer :: tslage2

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

tslagi  => rvoid1
tslage  => rvoid1
tslage2 => rvoid2

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbr(ncelet), rovsdt(ncelet))
allocate(gradv(3, 3, ncelet))
if (irijco.eq.1) then
  allocate(smbrts(6,ncelet))
  allocate(rovsdtts(6,6,ncelet))
endif

! Allocate other arrays, depending on user options
if (iturb.eq.30) then
  allocate(produc(6,ncelet))
endif

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

if (iturb.eq.30) then
  if (irijco.eq.1) then
    call field_get_val_prev_v(ivarfl(irij), cvara_rij)
  else
    call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
    call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
    call field_get_val_prev_s(ivarfl(ir33), cvara_r33)
    call field_get_val_prev_s(ivarfl(ir12), cvara_r12)
    call field_get_val_prev_s(ivarfl(ir13), cvara_r13)
    call field_get_val_prev_s(ivarfl(ir23), cvara_r23)
  endif
endif

call field_get_key_struct_var_cal_opt(ivarfl(iep), vcopt)

if(vcopt%iwarni.ge.1) then
  if (iturb.eq.30) then
    write(nfecra,1000)
  elseif (iturb.eq.31) then
    write(nfecra,1001)
  else
    write(nfecra,1002)
  endif
endif

!===============================================================================
! 1.1 Advanced init for EBRSM
!===============================================================================
! Automatic reinitialization at the end of the first iteration:
! wall distance y^+ is computed with -C log(1-alpha), where C=CL*Ceta*L*kappa, then y
! so we have an idea of the wall distance in complexe geometries.
! Then U is initialized with a Reichard lay
! Epsilon by 1/(kappa y), clipped next to the wall at its value for y^+=15
! k is given by a blending between eps/(2 nu)*y^2 and utau/sqrt(Cmu)
! The blending function is chosen so that the asymptotic behavior
! and give the correct pic of k

!TODO FIXME Are the BC uncompatible?
if (ntcabs.eq.1.and.reinit_turb.eq.1.and.iturb.eq.32) then

  allocate(grad(3,ncelet))

  ! Compute the gradient of Alpha
  iprev  = 0
  inc    = 1
  iccocg = 1

  call field_gradient_scalar(ivarfl(ial), iprev, imrgra, inc,     &
                             iccocg,                              &
                             grad)

  if (irijco.eq.1) then
    call field_get_val_v(ivarfl(irij), cvar_rij)
    call field_get_val_prev_v(ivarfl(irij), cvara_rij)
  else
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
    call field_get_val_s(ivarfl(ir12), cvar_r12)
    call field_get_val_s(ivarfl(ir13), cvar_r13)
    call field_get_val_s(ivarfl(ir23), cvar_r23)
    call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
    call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
    call field_get_val_prev_s(ivarfl(ir33), cvara_r33)
    call field_get_val_prev_s(ivarfl(ir12), cvara_r12)
    call field_get_val_prev_s(ivarfl(ir13), cvara_r13)
    call field_get_val_prev_s(ivarfl(ir23), cvara_r23)
  endif

  utaurf=0.05d0*uref

  call field_get_val_s(ivarfl(iep), cvar_ep)
  call field_get_val_s(ivarfl(ial), cvar_al)
  call field_get_val_v(ivarfl(iu), vel)

  nu0 = viscl0/ro0

  do iel = 1, ncel
    ! Compute the velocity magnitude
    xunorm = vel(1,iel)**2 + vel(2,iel)**2 + vel(3,iel)**2
    xunorm = sqrt(xunorm)

    ! y+ is bounded by 400, because in the Reichard profile, it corresponds to saturation (u>uref)
    cvar_al(iel) = max(min(cvar_al(iel),(1.d0-exp(-400.d0/50.d0))) &
    ,0.d0)
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
      limiter = min(utaurf/xunorm*(2.5d0*dlog(1.d0+0.4d0*ypa)            &
      +7.8d0*(1.d0-dexp(-ypa/11.d0)          &
      -(ypa/11.d0)*dexp(-0.33d0*ypa))),      &
      1.d0)
    endif

    vel(1,iel) = limiter*vel(1,iel)
    vel(2,iel) = limiter*vel(2,iel)
    vel(3,iel) = limiter*vel(3,iel)

    ut2 = 0.05d0*uref
    cvar_ep(iel) = utaurf**3*min(1.d0/(xkappa*15.d0*nu0/utaurf), &
    1.d0/(xkappa*ya))
    tke = cvar_ep(iel)/2.d0/nu0*ya**2             &
    * exp(-ypa/25.d0)**2                         &
    + ut2**2/0.3d0*(1.d0-exp(-ypa/25.d0))**2

    if (irijco.eq.0) then
      cvar_r11(iel) = alpha**3       *2.d0/3.d0        *tke &
                    + (1.d0-alpha**3)*(1.d0-xnal(1)**2)*tke
      cvar_r22(iel) = alpha**3       *2.d0/3.d0        *tke &
                    + (1.d0-alpha**3)*(1.d0-xnal(2)**2)*tke
      cvar_r33(iel) = alpha**3       *2.d0/3.d0        *tke &
                    + (1.d0-alpha**3)*(1.d0-xnal(3)**2)*tke
      cvar_r12(iel) = -(1.d0-alpha**3)*(xnal(1)*xnal(2))*tke
      cvar_r23(iel) = -(1.d0-alpha**3)*(xnal(2)*xnal(3))*tke
      cvar_r13(iel) = -(1.d0-alpha**3)*(xnal(1)*xnal(3))*tke
    else
      cvar_rij(1,iel) = alpha**3       *2.d0/3.d0        *tke &
                    + (1.d0-alpha**3)*(1.d0-xnal(1)**2)*tke
      cvar_rij(2,iel) = alpha**3       *2.d0/3.d0        *tke &
                    + (1.d0-alpha**3)*(1.d0-xnal(2)**2)*tke
      cvar_rij(3,iel) = alpha**3       *2.d0/3.d0        *tke &
                    + (1.d0-alpha**3)*(1.d0-xnal(3)**2)*tke
      cvar_rij(4,iel) = -(1.d0-alpha**3)*(xnal(1)*xnal(2))*tke
      cvar_rij(5,iel) = -(1.d0-alpha**3)*(xnal(2)*xnal(3))*tke
      cvar_rij(6,iel) = -(1.d0-alpha**3)*(xnal(1)*xnal(3))*tke
    end if
  enddo

  call field_current_to_previous(ivarfl(iu))
  if (irijco.eq.0) then
    call field_current_to_previous(ivarfl(ir11))
    call field_current_to_previous(ivarfl(ir22))
    call field_current_to_previous(ivarfl(ir33))
    call field_current_to_previous(ivarfl(ir12))
    call field_current_to_previous(ivarfl(ir23))
    call field_current_to_previous(ivarfl(ir13))
  else
    call field_current_to_previous(ivarfl(irij))
  endif

  deallocate(grad)
end if


!===============================================================================
! 2.1 Compute the velocity gradient
!===============================================================================

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

!===============================================================================
! 2.2 Compute the production term for Rij LRR (iturb =30)
!===============================================================================

if (iturb.eq.30) then
    do ii = 1 , 6
      do iel = 1, ncel
        produc(ii,iel) = 0.0d0
      enddo
    enddo
  if (irijco.eq.1) then
    do iel = 1 , ncel

      ! grad u

      produc(1,iel) = produc(1,iel)                                  &
                    - 2.0d0*(cvara_rij(1,iel)*gradv(1, 1, iel) +           &
                             cvara_rij(4,iel)*gradv(2, 1, iel) +           &
                             cvara_rij(6,iel)*gradv(3, 1, iel) )

      produc(4,iel) = produc(4,iel)                                  &
                    - (cvara_rij(4,iel)*gradv(1, 1, iel) +                 &
                       cvara_rij(2,iel)*gradv(2, 1, iel) +                 &
                       cvara_rij(5,iel)*gradv(3, 1, iel) )

      produc(6,iel) = produc(6,iel)                                  &
                    - (cvara_rij(6,iel)*gradv(1, 1, iel) +                 &
                       cvara_rij(5,iel)*gradv(2, 1, iel) +                 &
                       cvara_rij(3,iel)*gradv(3, 1, iel) )

      ! grad v

      produc(2,iel) = produc(2,iel)                                  &
                    - 2.0d0*(cvara_rij(4,iel)*gradv(1, 2, iel) +           &
                             cvara_rij(2,iel)*gradv(2, 2, iel) +           &
                             cvara_rij(5,iel)*gradv(3, 2, iel) )

      produc(4,iel) = produc(4,iel)                                  &
                    - (cvara_rij(1,iel)*gradv(1, 2, iel) +                 &
                       cvara_rij(4,iel)*gradv(2, 2, iel) +                 &
                       cvara_rij(6,iel)*gradv(3, 2, iel) )

      produc(5,iel) = produc(5,iel)                                  &
                    - (cvara_rij(6,iel)*gradv(1, 2, iel) +                 &
                       cvara_rij(5,iel)*gradv(2, 2, iel) +                 &
                       cvara_rij(3,iel)*gradv(3, 2, iel) )

      ! grad w

      produc(3,iel) = produc(3,iel)                                  &
                    - 2.0d0*(cvara_rij(6,iel)*gradv(1, 3, iel) +           &
                             cvara_rij(5,iel)*gradv(2, 3, iel) +           &
                             cvara_rij(3,iel)*gradv(3, 3, iel) )

      produc(6,iel) = produc(6,iel)                                  &
                    - (cvara_rij(1,iel)*gradv(1, 3, iel) +                 &
                       cvara_rij(4,iel)*gradv(2, 3, iel) +                 &
                       cvara_rij(6,iel)*gradv(3, 3, iel) )

      produc(5,iel) = produc(5,iel)                                  &
                    - (cvara_rij(4,iel)*gradv(1, 3, iel) +                 &
                       cvara_rij(2,iel)*gradv(2, 3, iel) +                 &
                       cvara_rij(5,iel)*gradv(3, 3, iel) )

    enddo
  else
    do iel = 1 , ncel

      ! grad u

      produc(1,iel) = produc(1,iel)                                  &
                    - 2.0d0*(cvara_r11(iel)*gradv(1, 1, iel) +           &
                             cvara_r12(iel)*gradv(2, 1, iel) +           &
                             cvara_r13(iel)*gradv(3, 1, iel) )

      produc(4,iel) = produc(4,iel)                                  &
                    - (cvara_r12(iel)*gradv(1, 1, iel) +                 &
                       cvara_r22(iel)*gradv(2, 1, iel) +                 &
                       cvara_r23(iel)*gradv(3, 1, iel) )

      produc(6,iel) = produc(6,iel)                                  &
                    - (cvara_r13(iel)*gradv(1, 1, iel) +                 &
                       cvara_r23(iel)*gradv(2, 1, iel) +                 &
                       cvara_r33(iel)*gradv(3, 1, iel) )

      ! grad v

      produc(2,iel) = produc(2,iel)                                  &
                    - 2.0d0*(cvara_r12(iel)*gradv(1, 2, iel) +           &
                             cvara_r22(iel)*gradv(2, 2, iel) +           &
                             cvara_r23(iel)*gradv(3, 2, iel) )

      produc(4,iel) = produc(4,iel)                                  &
                    - (cvara_r11(iel)*gradv(1, 2, iel) +                 &
                       cvara_r12(iel)*gradv(2, 2, iel) +                 &
                       cvara_r13(iel)*gradv(3, 2, iel) )

      produc(5,iel) = produc(5,iel)                                  &
                    - (cvara_r13(iel)*gradv(1, 2, iel) +                 &
                       cvara_r23(iel)*gradv(2, 2, iel) +                 &
                       cvara_r33(iel)*gradv(3, 2, iel) )

      ! grad w

      produc(3,iel) = produc(3,iel)                                  &
                    - 2.0d0*(cvara_r13(iel)*gradv(1, 3, iel) +           &
                             cvara_r23(iel)*gradv(2, 3, iel) +           &
                             cvara_r33(iel)*gradv(3, 3, iel) )

      produc(6,iel) = produc(6,iel)                                  &
                    - (cvara_r11(iel)*gradv(1, 3, iel) +                 &
                       cvara_r12(iel)*gradv(2, 3, iel) +                 &
                       cvara_r13(iel)*gradv(3, 3, iel) )

      produc(5,iel) = produc(5,iel)                                  &
                    - (cvara_r12(iel)*gradv(1, 3, iel) +                 &
                       cvara_r22(iel)*gradv(2, 3, iel) +                 &
                       cvara_r23(iel)*gradv(3, 3, iel) )

    enddo
  endif
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

  iprev = 1
  inc = 1
  iccocg = 1

  call field_gradient_scalar(ivarfl(isca(iscalt)), 1, imrgra, inc, &
                             iccocg,                               &
                             gradro)

  ! gradro stores: - rho grad(theta)/theta
  ! grad(ro) and grad(teta) have opposite sign
  do iel = 1, ncel
    rhothe = cromo(iel)/cvara_scalt(iel)
    gradro(1, iel) = -rhothe*gradro(1, iel)
    gradro(2, iel) = -rhothe*gradro(2, iel)
    gradro(3, iel) = -rhothe*gradro(3, iel)
  enddo

else if (igrari.eq.1) then
  ! Allocate a temporary array for the gradient calculation
  allocate(gradro(3,ncelet))

! Boundary conditions: Dirichlet romb
!   We use viscb to store the relative coefficient of rom
!   We impose in Dirichlet (coefa) the value romb

  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

! The choice below has the advantage to be simple
  call field_get_key_struct_var_cal_opt(ivarfl(ir11), vcopt)

  nswrgp = vcopt%nswrgr
  imligp = vcopt%imligr
  iwarnp = vcopt%iwarni
  epsrgp = vcopt%epsrgr
  climgp = vcopt%climgr
  extrap = vcopt%extrag

  f_id0 = -1
  iccocg = 1

  ! If we extrapolate the source terms and rho, we use cpdt rho^n
  if(isto2t.gt.0.and.iroext.gt.0) then
    call field_get_val_prev_s(icrom, cromo)
    call field_get_val_prev_s(ibrom, bromo)
  else
    call field_get_val_s(icrom, cromo)
    call field_get_val_s(ibrom, bromo)
  endif

  call gradient_s                                                 &
 ( f_id0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , epsrgp , climgp , extrap ,                            &
   cromo  , bromo  , viscb           ,                            &
   gradro )

endif

!===============================================================================
! 4.  Loop on the variables Rij (6 variables)
!     The order is R11 R22 R33 R12 R13 R23 (The place of those variables
!      is IR11.    ..
!     We solve the equation in a routine similar to covofi.f90
!===============================================================================

if (irijco.eq.1) then
  ivar = irij

  if (iilagr.eq.2) then
    tslage2 => tslagr(1:ncelet,itsr11:itsr13)
    tslagi  => tslagr(1:ncelet,itsli)
  endif


  ! Rij-epsilon standard (LRR)
  if (iturb.eq.30) then !TODO
    call resrij2 &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   produc , gradro ,                                              &
   ckupdc , smacel ,                                              &
   viscf  , viscb  ,                                              &
   tslage , tslagi ,                                              &
   smbrts   , rovsdtts )

  ! Rij-epsilon SSG or EBRSM
  elseif (iturb.eq.31.or.iturb.eq.32) then

    call resssg2 &
  ( nvar    , nscal  , ncepdp , ncesmp ,                            &
    ivar    ,                                                       &
    icepdc  , icetsm , itypsm ,                                     &
    dt      ,                                                       &
    gradv   , gradro ,                                              &
    ckupdc  , smacel ,                                              &
    viscf   , viscb  ,                                              &
    tslage2 , tslagi ,                                              &
    smbrts  , rovsdtts )
  endif
else
  do isou = 1, 6
    if    (isou.eq.1) then
      ivar   = ir11
    elseif(isou.eq.2) then
      ivar   = ir22
    elseif(isou.eq.3) then
      ivar   = ir33
    elseif(isou.eq.4) then
      ivar   = ir12
    elseif(isou.eq.5) then
      ivar   = ir23
    elseif(isou.eq.6) then
      ivar   = ir13
    endif

    if (iilagr.eq.2) then
      iitsla = itsr11 + (isou-1)
      tslage => tslagr(1:ncelet,iitsla)
      tslagi => tslagr(1:ncelet,itsli)
    endif

    ! Rij-epsilon standard (LRR)
    if (iturb.eq.30) then
      call resrij &
   ( nvar   , nscal  , ncepdp , ncesmp ,                            &
     ivar   , isou   ,                                              &
     icepdc , icetsm , itypsm ,                                     &
     dt     ,                                                       &
     produc , gradro ,                                              &
     ckupdc , smacel ,                                              &
     viscf  , viscb  ,                                              &
     tslage , tslagi ,                                              &
     smbr   , rovsdt )

    ! Rij-epsilon SSG or EBRSM
    elseif (iturb.eq.31.or.iturb.eq.32) then
        call resssg &
      ( nvar   , nscal  , ncepdp , ncesmp ,                            &
        ivar   , isou   ,                                              &
        icepdc , icetsm , itypsm ,                                     &
        dt     ,                                                       &
        gradv  , gradro ,                                              &
        ckupdc , smacel ,                                              &
        viscf  , viscb  ,                                              &
        tslage , tslagi ,                                              &
        smbr   , rovsdt )
    endif

  enddo
endif

!===============================================================================
! 5. Solve Epsilon
!===============================================================================

call reseps &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   gradv  , produc , gradro ,                                     &
   ckupdc , smacel ,                                              &
   viscf  , viscb  ,                                              &
   tslagr ,                                                       &
   smbr   , rovsdt )

!===============================================================================
! 6. Clipping
!===============================================================================

if (iturb.eq.32) then
  iclip = 1
else
  iclip = 2
endif

if (irijco.eq.1) then
  call clprij2(ncelet, ncel, iclip)
else
  call clprij(ncelet, ncel, iclip)
endif

! Free memory
deallocate(viscf, viscb)
deallocate(smbr, rovsdt)
if (allocated(gradro)) deallocate(gradro)
if (allocated(produc)) deallocate(produc)
deallocate(gradv)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                      &
'   ** Resolution du Rij-EPSILON LRR'                   ,/,&
'      -----------------------------'             ,/)
 1001 format(/,                                      &
'   ** Resolution du Rij-EPSILON SSG'                   ,/,&
'      -----------------------------'             ,/)
 1002 format(/,                                      &
'   ** Resolution du Rij-EPSILON EBRSM'                 ,/,&
'      -------------------------------'           ,/)

#else

 1000 format(/,                                      &
'   ** Solving Rij-EPSILON LRR'                   ,/,&
'      -----------------------'                   ,/)
 1001 format(/,                                      &
'   ** Solving Rij-EPSILON SSG'                   ,/,&
'      -----------------------'                   ,/)
 1002 format(/,                                      &
'   ** Solving Rij-EPSILON EBRSM'                 ,/,&
'      -------------------------'                 ,/)

#endif

!----
! End
!----

return

end subroutine
