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
! Function :
! --------

!> \file driflu.f90
!>
!> \brief Compute the modified convective flux for scalars with a drift.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iflid         index of the current drift scalar field
!> \param[in]     dt            time step (per cell)
!> \param[in,out] imasfl        scalar mass flux at interior face centers
!> \param[in,out] bmasfl        scalar mass flux at boundary face centers
!> \param[in,out] divflu        divergence of drift flux
!______________________________________________________________________________

subroutine driflu &
( iflid  ,                                                       &
  dt     ,                                                       &
  imasfl , bmasfl ,                                              &
  divflu  )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: nvar, ndimfb
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use ppthch
use pointe
use field
use mesh
use parall
use period
use cpincl
use cs_coal_incl
use ppincl
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iflid

double precision dt(ncelet)
double precision imasfl(nfac), bmasfl(nfabor)
double precision divflu(ncelet)

! Local variables

integer          ivar
integer          ifac  , iel
integer          init  , inc   , iccocg
integer          ifcvsl, iflmas, iflmab
integer          imrgrp, nswrgp, imligp, iwarnp
integer          iconvp, idiffp, imvisp
integer          isou  , jsou
integer          f_id  , f_id0
integer          iflmb0, iphydp, ivisep, itypfl
integer          keysca, iscal, keydri, iscdri, icvflb
integer          keyccl
integer          icla, jcla, jvar
integer          ivoid(1)
integer          ielup, id_x1, id_vdp_i, imasac, id_pro, id_drift

double precision epsrgp, climgp, extrap
double precision visls_0
double precision thetap
double precision rhovdt
double precision omegaa
double precision rho

double precision rvoid(1)

character(len=80) :: fname

double precision, dimension(:), allocatable :: w1, viscce, rtrace
double precision, dimension(:), allocatable :: coefap, coefbp
double precision, dimension(:), allocatable :: cofafp, cofbfp
double precision, dimension(:,:), allocatable :: coefa1
double precision, dimension(:,:,:), allocatable :: coefb1
double precision, dimension(:,:), allocatable :: dudt
double precision, dimension(:), allocatable :: viscf, viscb
double precision, dimension(:), allocatable :: flumas, flumab

double precision, dimension(:), pointer :: cpro_taup
double precision, dimension(:), pointer :: cpro_taufpt
double precision, dimension(:,:), pointer :: cpro_drift
double precision, dimension(:,:), pointer :: coefav, cofafv
double precision, dimension(:,:,:), pointer :: coefbv, cofbfv
double precision, dimension(:), pointer :: imasfl_mix, bmasfl_mix
double precision, dimension(:,:), pointer :: vdp_i
double precision, dimension(:), pointer :: x2
double precision, dimension(:), pointer :: imasfl_gas, bmasfl_gas
double precision, dimension(:), pointer :: cpro_x1, bpro_x1
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:,:), pointer :: vel, vela
double precision, dimension(:), pointer :: cvar_k
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:), pointer :: visct, cpro_viscls
double precision, dimension(:), pointer :: cvara_var

type(var_cal_opt) :: vcopt, vcopt_u
type(var_cal_opt), target   :: vcopt_loc
type(var_cal_opt), pointer  :: p_k_value
type(c_ptr)                 :: c_k_value

!===============================================================================

!===============================================================================
! 0. Key words for fields
!===============================================================================

! Key id for scalar id
call field_get_key_id("scalar_id", keysca)
call field_get_key_int(iflid, keysca, iscal)

ivar = isca(iscal)
f_id0 = -1

call field_get_val_prev_s(ivarfl(ivar), cvara_var)

! Eventual scalar class:
!   0: scalar belonging to no class
!  -1: deduced continuous (gas) class
!  >0: particle class
call field_get_key_id("scalar_class", keyccl)
call field_get_key_int(ivarfl(ivar), keyccl, icla)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)
call field_get_key_int(iflid, keydri, iscdri)

! Massic fraction of gas
call field_get_id_try("x_c", id_x1)
if (id_x1.ne.-1) then
  call field_get_val_s(id_x1, cpro_x1)

  ! Mass fraction of the gas at the boundary
  call field_get_val_s_by_name("b_x_c", bpro_x1)

  ! Get the mass flux of the continuous phase (gas)
  ! that is the negative scalar class
  do jvar = 1, nvar
    call field_get_key_int(ivarfl(jvar), keyccl, jcla)
    if (jcla.eq.-1) then
      call field_get_key_int(ivarfl(jvar), kimasf, iflmas)
      call field_get_key_int(ivarfl(jvar), kbmasf, iflmab)
    endif
  enddo
  call field_get_val_s(iflmas, imasfl_gas)
  ! Pointer to the Boundary mass flux
  call field_get_val_s(iflmab, bmasfl_gas)

endif

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)

! Pointers to the mass fluxes of the mix (based on mix velocity)
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl_mix)
call field_get_val_s(iflmab, bmasfl_mix)

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)
call field_get_val_prev_v(ivarfl(iu), vela)

!===============================================================================
! 1. Initialization
!===============================================================================

! --- Physical properties
call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)
call field_get_val_s(ivisct, visct)

if (itytur.eq.3) then
  call field_get_val_v(ivarfl(irij), cvar_rij)
elseif (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then
  call field_get_val_s(ivarfl(ik), cvar_k)
endif

! --- Brownian diffusivity
call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_viscls)
endif

! Name of the scalar
call field_get_name(ivarfl(ivar), fname)

! Vector containing all the additional convective terms
allocate(dudt(3, ncelet))
allocate(w1(ncelet), viscce(ncelet))
allocate(coefap(ndimfb), coefbp(ndimfb))
allocate(cofafp(ndimfb), cofbfp(ndimfb))
allocate(coefa1(3, ndimfb), coefb1(3, 3, ndimfb))
allocate(viscf(nfac), viscb(nfabor))
allocate(flumas(nfac), flumab(nfabor))

if (btest(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)) then
  ! Index of the corresponding relaxation time (cpro_taup)
  call field_get_id_try('drift_tau_'//trim(fname), id_pro)
  if (id_pro.ne.-1) call field_get_val_s(id_pro, cpro_taup)

  ! Index of the corresponding relaxation time (cpro_taup)
  call field_get_id_try('drift_vel_'//trim(fname), id_drift)
  if (id_drift.ne.-1) call field_get_val_v(id_drift, cpro_drift)

  ! Index of the corresponding interaction time particle--eddies (cpro_taufpt)
  if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS)) then
    call field_get_id('drift_turb_tau_'//trim(fname), f_id)
    call field_get_val_s(f_id, cpro_taufpt)
  endif

  ! Initialization of the convection flux for the current particle class
  do ifac = 1, nfac
    viscf(ifac) = 0.d0
    flumas(ifac) = 0.d0
  enddo

  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
    flumab(ifac) = 0.d0
  enddo

  ! Initialization of the gas "class" convective flux by the
  ! first particle "class":
  !  it is initialized by the mass flux of the bulk
  if (icla.eq.1.and.id_x1.ne.-1) then
    do ifac = 1, nfac
      imasfl_gas(ifac) = imasfl_mix(ifac)
    enddo

    do ifac = 1, nfabor
      bmasfl_gas(ifac) = bmasfl_mix(ifac)
    enddo
  endif

!===============================================================================
! 2. initialize the additional convective flux with the gravity term
!===============================================================================

  ! Test if a deviation velocity of particles class exists

  if (icla.ge.1) then

    write(fname,'(a,i2.2)')'vd_p_' ,icla
    call field_get_id_try(fname, id_vdp_i)

    if (id_vdp_i.ne.-1) then
      call field_get_val_v(id_vdp_i, vdp_i)

      do iel = 1, ncel

        rho = crom(iel)
        cpro_drift(1, iel) = rho*vdp_i(1, iel)
        cpro_drift(2, iel) = rho*vdp_i(2, iel)
        cpro_drift(3, iel) = rho*vdp_i(3, iel)

      enddo
    endif

  else if (icla.ge.0.and.id_pro.ne.-1.and.id_drift.ne.-1) then

    do iel = 1, ncel

      rho = crom(iel)
      cpro_drift(1, iel) = rho*cpro_taup(iel)*gx
      cpro_drift(2, iel) = rho*cpro_taup(iel)*gy
      cpro_drift(3, iel) = rho*cpro_taup(iel)*gz

    enddo

  endif

!===============================================================================
! 3. Computation of the turbophoresis and the thermophoresis terms
!===============================================================================

  ! Initialized to 0
  do iel = 1, ncel
    viscce(iel) = 0.d0
  enddo

  if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS).and.iturb.ne.0) then

    ! The diagonal part is easy to implicit (Grad (K) . n = (K_j - K_i)/IJ)
    ! Compute the K=1/3*trace(R) coefficient (diffusion of Zaichik)

    if (itytur.eq.3) then

      allocate(rtrace(ncel))
      do iel = 1, ncel
        rtrace(iel) = cvar_rij(1,iel) + cvar_rij(2,iel) + cvar_rij(3,iel)
      enddo
      do iel = 1, ncel
        ! Correction by Omega
        omegaa = cpro_taup(iel)/cpro_taufpt(iel)
        ! FIXME: use idifft or not?
        viscce(iel) = 1.d0/3.d0 * cpro_taup(iel)/(1.d0+omegaa)*rtrace(iel)
      enddo
      deallocate(rtrace)

    elseif (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then

      do iel = 1, ncel
        ! Correction by Omega
        omegaa = cpro_taup(iel)/cpro_taufpt(iel)
        viscce(iel) = 2.d0/3.d0*cpro_taup(iel)/(1.d0+omegaa)*cvar_k(iel)
      enddo

    endif

  endif


  if (btest(iscdri, DRIFT_SCALAR_THERMOPHORESIS)) then

    ! cpro_viscls(iel): contains the Brownian motion
    !-----------------------------------------------

    if (ifcvsl.ge.0) then

      do iel = 1, ncel
        viscce(iel) = viscce(iel) + cpro_viscls(iel)/crom(iel)
      enddo

    else

      call field_get_key_double(ivarfl(ivar), kvisl0, visls_0)

      do iel = 1, ncel
        viscce(iel) = viscce(iel) + visls_0/crom(iel)
      enddo

    endif

  endif

  if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS).or.    &
      btest(iscdri, DRIFT_SCALAR_THERMOPHORESIS)) then

    iphydp = 0
    inc    = 1
    iccocg = 1
    imrgrp = vcopt%imrgra
    nswrgp = vcopt%nswrgr
    imligp = vcopt%imligr
    iwarnp = vcopt%iwarni
    epsrgp = vcopt%epsrgr
    climgp = vcopt%climgr
    imvisp = vcopt%imvisf
    extrap = 0

    ! Face diffusivity of rho to compute rho*(Grad K . n)_face
    do iel = 1, ncel
      w1(iel) = crom(iel)
    enddo

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(w1)
    endif

    call viscfa &
    ( imvisp ,            &
      w1     ,            &
      viscf  , viscb  )

    ! Homogeneous Neumann BC
    do ifac = 1, nfabor
      ! BCs for gradients
      coefap(ifac) = 0.d0
      coefbp(ifac) = 1.d0
      ! BCs for fluxes
      cofafp(ifac) = 0.d0
      cofbfp(ifac) = 0.d0
    enddo

    init   = 0

    ! The computed convective flux has the dimension of rho*velocity
    call itrmas &
   ( f_id0  , init , inc , imrgrp , iccocg , nswrgp , imligp , iphydp ,      &
     0      , iwarnp ,                                                       &
     epsrgp , climgp , extrap ,                                              &
     rvoid  ,                                                                &
     viscce ,                                                                &
     coefap , coefbp ,                                                       &
     cofafp , cofbfp ,                                                       &
     viscf  , viscb  ,                                                       &
     w1     ,                                                                &
     flumas , flumab )

  ! TODO add extradiagonal part
  endif

!===============================================================================
! 4. Centrifugal force (particular derivative Du/Dt)
!===============================================================================

  if (btest(iscdri, DRIFT_SCALAR_CENTRIFUGALFORCE)) then

    do iel = 1, ncel

      rhovdt = crom(iel)*volume(iel)/dt(iel)

      dudt(1,iel) = -rhovdt*(vel(1,iel)-vela(1,iel))
      dudt(2,iel) = -rhovdt*(vel(2,iel)-vela(2,iel))
      dudt(3,iel) = -rhovdt*(vel(3,iel)-vela(3,iel))
    enddo

    iconvp = 1
    idiffp = 0
    inc    = 1
    ivisep = 0
    icvflb = 0

    ! Reset viscf and viscb
    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo

    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

    ! Get Boundary conditions of the velocity
    call field_get_coefa_v (ivarfl(iu), coefav)
    call field_get_coefb_v (ivarfl(iu), coefbv)
    call field_get_coefaf_v(ivarfl(iu), cofafv)
    call field_get_coefbf_v(ivarfl(iu), cofbfv)

    ! The added convective scalar mass flux is:
    !      (thetap*Y_\face-imasac*Y_\celli)*mf.
    ! When building the implicit part of the rhs, one
    ! has to impose 1 on mass accumulation.
    imasac = 1

    vcopt_loc = vcopt_u

    vcopt_loc%iconv  = iconvp
    vcopt_loc%istat  = -1
    vcopt_loc%idiff  = idiffp
    vcopt_loc%idifft = -1
    vcopt_loc%iswdyn = -1
    vcopt_loc%nswrsm = -1
    vcopt_loc%iwgrec = 0
    vcopt_loc%blend_st = 0 ! Warning, may be overwritten if a field
    vcopt_loc%epsilo = -1
    vcopt_loc%epsrsm = -1

    p_k_value => vcopt_loc
    c_k_value = equation_param_from_vcopt(c_loc(p_k_value))

    ! Warning: cs_balance_vector adds "-( grad(u) . rho u)"
    call cs_balance_vector &
   ( idtvar , ivarfl(iu)      , imasac , inc    , ivisep ,                   &
     c_k_value                , vel    , vel    , coefav , coefbv , cofafv , &
     cofbfv , imasfl_mix      , bmasfl_mix      , viscf  , viscb  , rvoid  , &
     rvoid  , rvoid  , rvoid  , rvoid  ,                                     &
     icvflb , ivoid  , dudt   )

    do iel = 1, ncel
      cpro_drift(1, iel) = cpro_drift(1, iel)                       &
                         + cpro_taup(iel)*dudt(1, iel)/volume(iel)
      cpro_drift(2, iel) = cpro_drift(2, iel)                       &
                         + cpro_taup(iel)*dudt(2, iel)/volume(iel)
      cpro_drift(3, iel) = cpro_drift(3, iel)                       &
                         + cpro_taup(iel)*dudt(3, iel)/volume(iel)
    enddo

  endif

!===============================================================================
! 5. Electrophoresis term
!===============================================================================

  if (btest(iscdri, DRIFT_SCALAR_ELECTROPHORESIS)) then

    !TODO
    call csexit(1)

  endif

!===============================================================================
! 6. Finalization of the mass flux of the current class
!===============================================================================

  ! For all scalar with a drift excpted the gas phase which is deduced
  ! And for those whom mass flux is imposed elsewhere
  if (icla.ge.0.and.(.not.btest(iscdri, DRIFT_SCALAR_IMPOSED_MASS_FLUX))) then

    ! Homogeneous Neumann at the boundary
    if (btest(iscdri, DRIFT_SCALAR_ZERO_BNDY_FLUX)) then
      do ifac = 1, nfabor
        do isou = 1, 3
          coefa1(isou, ifac) = 0.d0
          do jsou = 1, 3
            coefb1(isou, jsou, ifac) = 0.d0
          enddo
        enddo
      enddo
    else if (btest(iscdri, DRIFT_SCALAR_ZERO_BNDY_FLUX_AT_WALLS)) then
      do ifac = 1, nfabor
        do isou = 1, 3
          coefa1(isou, ifac) = 0.d0
          do jsou = 1, 3
            coefb1(isou, jsou, ifac) = 0.d0
          enddo
          if (itypfb(ifac).ne.iparoi .and. itypfb(ifac).ne.iparug) then
            coefb1(isou, isou, ifac) = 1.d0
          endif
        enddo
      enddo
    else
      do ifac = 1, nfabor
        do isou = 1, 3
          coefa1(isou, ifac) = 0.d0
          do jsou = 1, 3
            coefb1(isou, jsou, ifac) = 0.d0
          enddo
          coefb1(isou, isou, ifac) = 1.d0
        enddo
      enddo
    endif

    init   = 0
    inc    = 1
    iflmb0 = 0
    itypfl = 0 ! drift has already been multiplied by rho
    imrgrp = vcopt%imrgra
    nswrgp = vcopt%nswrgr
    imligp = vcopt%imligr
    iwarnp = vcopt%iwarni
    epsrgp = vcopt%epsrgr
    climgp = vcopt%climgr
    extrap = 0

    call inimav &
     ( f_id0  , itypfl ,                                              &
       iflmb0 , init   , inc    , imrgrp , nswrgp , imligp ,          &
       iwarnp ,                                                       &
       epsrgp , climgp ,                                              &
       crom   , brom   ,                                              &
       cpro_drift  ,                                                  &
       coefa1 , coefb1 ,                                              &
       flumas , flumab )

    ! Update the convective flux, exception for the Gas "class"
    do ifac = 1, nfac
      imasfl(ifac) = imasfl_mix(ifac) + flumas(ifac)
    enddo

    do ifac = 1, nfabor
      bmasfl(ifac) = bmasfl_mix(ifac) + flumab(ifac)
    enddo


!===============================================================================
! 7. Deduce the convective flux of the gas "class" by removing the flux
!     of the current particle "class":
!  (rho x1 V1)_ij = (rho Vs)_ij - sum_classes (rho x2 V2)_ij
!===============================================================================

    if (icla.ge.1) then

      write(fname,'(a,i2.2)')'x_p_' ,icla
      call field_get_id_try(fname, f_id)

      if (f_id.ne.-1) then
        call field_get_val_s(f_id, x2)

        do ifac = 1, nfac

          ! Upwind value of x2 at the face, consistent with the
          !  other transport equations
          ielup = ifacel(2,ifac)
          if (imasfl(ifac).ge.0.d0) ielup = ifacel(1,ifac)

          imasfl_gas(ifac) = imasfl_gas(ifac) - x2(ielup)*imasfl(ifac)
        enddo

        do ifac = 1, nfabor
          ! TODO Upwind value of x2 at the face, consistent with the
          !  other transport equations
          !if (bmasfl(ifac).ge.0.d0)
          ielup = ifabor(ifac)
          bmasfl_gas(ifac) = bmasfl_gas(ifac) - x2(ielup)*bmasfl(ifac)
        enddo
      endif
    endif

  ! Finalize the convective flux of the gas "class" by scaling by x1
  !  (rho x1 V1)_ij = (rho Vs)_ij - sum_classes (rho x2 V2)_ij
  ! Warning, x1 at the face must be computed so that it is consistent
  ! with an upwind scheme on (rho V1)
  else if (icla.eq.-1.and.id_x1.ne.-1) then

    do ifac = 1, nfac

      ! Upwind value of x2 at the face, consistent with the
      !  other transport equations
      ielup = ifacel(2,ifac)
      if (imasfl_gas(ifac).ge.0.d0) ielup = ifacel(1,ifac)

      imasfl_gas(ifac) = imasfl_gas(ifac)/cpro_x1(ielup)
    enddo

    do ifac = 1, nfabor
      ! Upwind value of x1 at the face, consistent with the
      !  other transport equations
      ielup = ifabor(ifac)
      if (bmasfl_gas(ifac).lt.0.d0) then
        bmasfl_gas(ifac) = bmasfl_gas(ifac)/bpro_x1(ifac)
      else
        bmasfl_gas(ifac) = bmasfl_gas(ifac)/cpro_x1(ielup)
      endif
    enddo
  endif

endif

!===============================================================================
! 8. Mass aggregation term of the additional part "div(rho(u_p-u_f))"
!===============================================================================

init = 1
iconvp = vcopt%iconv
thetap = vcopt%thetav!FIXME not multiplied by theta?

! recompute the difference between mixture and the class
if (btest(iscdri, DRIFT_SCALAR_IMPOSED_MASS_FLUX)) then
  do ifac = 1, nfac
    flumas(ifac) = - imasfl_mix(ifac)
  enddo

  do ifac = 1, nfabor
    flumab(ifac) = - bmasfl_mix(ifac)
  enddo
else
  do ifac = 1, nfac
    flumas(ifac) = imasfl(ifac) - imasfl_mix(ifac)
  enddo

  do ifac = 1, nfabor
    flumab(ifac) = bmasfl(ifac) - bmasfl_mix(ifac)
  enddo
endif

call divmas(init, flumas, flumab, divflu)

! Free memory
deallocate(viscce)
deallocate(dudt)
deallocate(w1)
deallocate(viscf, viscb)
deallocate(flumas, flumab)
deallocate(coefap, coefbp)
deallocate(cofafp, cofbfp)
deallocate(coefa1, coefb1)

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
