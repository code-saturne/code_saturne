!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
!> \param[in,out] rovsdt        unsteady term and mass aggregation term
!> \param[in,out] smbrs         right hand side for the scalar iscal
!______________________________________________________________________________

subroutine driflu &
( iflid  ,                                                       &
  dt     ,                                                       &
  imasfl , bmasfl ,                                              &
  rovsdt , smbrs  )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb, nvar
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
double precision imasfl(nfac), bmasfl(ndimfb)
double precision rovsdt(ncelet), smbrs(ncelet)

! Local variables

integer          ivar
integer          ifac  , iel
integer          init  , inc   , iccocg
integer          ifcvsl, iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp
integer          ircflp, ischcp, isstpp
integer          isou  , jsou
integer          f_id  , f_id0
integer          iflmb0, idftnp, iphydp, ivisep, itypfl
integer          keysca, iscal, keydri, iscdri, icvflb
integer          keyccl
integer          icla, jcla, jvar
integer          ivoid(1)
integer          ielup, id_x1, id_vdp_i, imasac

double precision epsrgp, climgp, extrap, blencp
double precision thetap
double precision rhovdt
double precision omegaa
double precision rho
double precision relaxp

double precision rvoid(1)

character(len=80) :: fname

double precision, dimension(:), allocatable :: w1, viscce
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
double precision, dimension(:), pointer :: x1, b_x1
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:,:), pointer :: vel, vela
double precision, dimension(:), pointer :: cvar_k
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: visct, cpro_viscls
double precision, dimension(:), pointer :: cvara_var

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
  call field_get_val_s(id_x1, x1)

  ! Mass fraction of the gas at the boundary
  call field_get_val_s_by_name("b_x_c", b_x1)

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
call field_get_val_s(iprpfl(ivisct), visct)

if (itytur.eq.3) then
  call field_get_val_s(ivarfl(ir11), cvar_r11)
  call field_get_val_s(ivarfl(ir22), cvar_r22)
  call field_get_val_s(ivarfl(ir33), cvar_r33)
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
  call field_get_id('drift_tau_'//trim(fname), f_id)
  call field_get_val_s(f_id, cpro_taup)

  ! Index of the corresponding relaxation time (cpro_taup)
  call field_get_id('drift_vel_'//trim(fname), f_id)
  call field_get_val_v(f_id, cpro_drift)

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
  if (icla.eq.1) then
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
  write(fname,'(a,i2.2)')'vd_p_' ,icla
  call field_get_id_try(fname, id_vdp_i)

  if (icla.ge.1.and.id_vdp_i.ne.-1) then

    call field_get_val_v(id_vdp_i, vdp_i)

    do iel = 1, ncel

      rho = crom(iel)
      cpro_drift(1, iel) = rho*vdp_i(1, iel)
      cpro_drift(2, iel) = rho*vdp_i(2, iel)
      cpro_drift(3, iel) = rho*vdp_i(3, iel)

    enddo

  else if (icla.ge.0) then

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

    ! Compute the K=1/3*trace(K) coefficient (diffusion of Zaichik)

    if (itytur.eq.3) then

      do iel = 1, ncel

        ! Correction by Omega
        omegaa = cpro_taup(iel)/cpro_taufpt(iel)
        ! FIXME: use idifft or not?
        viscce(iel) = 1.d0/3.d0                                        &
                    * cpro_taup(iel)/(1.d0+omegaa)*( cvar_r11(iel)     &
                                                   + cvar_r22(iel)     &
                                                   + cvar_r33(iel) )
      enddo

    elseif (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then

      do iel = 1, ncel

        ! Correction by Omega
        omegaa = cpro_taup(iel)/cpro_taufpt(iel)

        viscce(iel) = 2.d0/3.d0*cpro_taup(iel)/(1.d0+omegaa)*cvar_k(iel)
      enddo

    else

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

      do iel = 1, ncel
        viscce(iel) = viscce(iel) + visls0(iscal)/crom(iel)
      enddo

    endif

  endif

  if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS).or.    &
      btest(iscdri, DRIFT_SCALAR_THERMOPHORESIS)) then

    iphydp = 0
    inc    = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    ! Face diffusivity of rho to compute rho*(Grad K . n)_face
    do iel = 1, ncel
      w1(iel) = crom(iel)
    enddo

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(w1)
    endif

    call viscfa &
    ( imvisf ,            &
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
   ( f_id0  , init , inc , imrgra , iccocg , nswrgp , imligp , iphydp ,      &
     iwarnp ,                                                                &
     epsrgp , climgp , extrap ,                                              &
     rvoid  ,                                                                &
     viscce ,                                                                &
     coefap , coefbp ,                                                       &
     cofafp , cofbfp ,                                                       &
     viscf  , viscb  ,                                                       &
     w1     , w1     , w1     ,                                              &
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
    nswrgp = nswrgr(iu)
    imligp = imligr(iu)
    ircflp = ircflu(iu)
    ischcp = ischcv(iu)
    isstpp = isstpc(iu)
    inc    = 1
    ivisep = 0
    iwarnp = iwarni(iu)
    idftnp = idften(iu)
    blencp = blencv(iu)
    epsrgp = epsrgr(iu)
    climgp = climgr(iu)
    thetap = thetav(iu)
    relaxp = relaxv(iu)
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

    ! Warning: bilsc adds "-( grad(u) . rho u)"
    call bilscv &
   ( idtvar , iu     , iconvp , idiffp , nswrgp , imligp , ircflp , &
     ischcp , isstpp , inc    , imrgra , ivisep ,                   &
     iwarnp , idftnp , imasac ,                                     &
     blencp , epsrgp , climgp , relaxp , thetap ,                   &
     vel    , vel    ,                                              &
     coefav , coefbv , cofafv , cofbfv ,                            &
     imasfl_mix , bmasfl_mix ,                                      &
     viscf  , viscb  , rvoid  , rvoid  ,                            &
     icvflb , ivoid  ,                                              &
     dudt   )

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
  if (icla.ge.0) then

    ! Zero additional flux at the boundary
    do ifac = 1, nfabor

      do isou = 1, 3
        coefa1(isou, ifac) = 0.d0
        do jsou = 1, 3
          coefb1(isou, jsou, ifac) = 0.d0
        enddo
      enddo

    enddo

    init   = 0
    inc    = 1
    iflmb0 = 0
    itypfl = 0 ! drift has already been multiplied by rho
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    call inimav &
     ( f_id0  , itypfl ,                                              &
       iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
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
      call field_get_val_s_by_name(fname, x2)

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

  ! Finalize the convective flux of the gas "class" by scaling by x1
  !  (rho x1 V1)_ij = (rho Vs)_ij - sum_classes (rho x2 V2)_ij
  ! Warning, x1 at the face must be computed so that it is consistent
  ! with an upwind scheme on (rho V1)
  else if (icla.eq.-1) then

    do ifac = 1, nfac

      ! Upwind value of x2 at the face, consistent with the
      !  other transport equations
      ielup = ifacel(2,ifac)
      if (imasfl_gas(ifac).ge.0.d0) ielup = ifacel(1,ifac)

      imasfl_gas(ifac) = imasfl_gas(ifac)/x1(ielup)
    enddo

    do ifac = 1, nfabor
      ! Upwind value of x1 at the face, consistent with the
      !  other transport equations
      ielup = ifabor(ifac)
      if (bmasfl_gas(ifac).lt.0.d0) then
        bmasfl_gas(ifac) = bmasfl_gas(ifac)/b_x1(ifac)
      else
        bmasfl_gas(ifac) = bmasfl_gas(ifac)/x1(ielup)
      endif
    enddo
  endif

endif

!===============================================================================
! 8. Mass aggregation term of the additional part "div(rho(u_p-u_f))"
!===============================================================================

init = 1
iconvp = iconv(ivar)
thetap = thetav(ivar)

! recompute the difference between mixture and the class
do ifac = 1, nfac
  flumas(ifac) = imasfl(ifac) - imasfl_mix(ifac)
enddo

do ifac = 1, nfabor
  flumab(ifac) = bmasfl(ifac) - bmasfl_mix(ifac)
enddo

call divmas(init, flumas, flumab, w1)

! NB: if the porosity module is swiched on, the the porosity is already
! taken into account in w1

! --> mass aggregation term
do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel) + iconvp*thetap*w1(iel)
  smbrs(iel) = smbrs(iel) - iconvp*w1(iel)*cvara_var(iel)
enddo

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
