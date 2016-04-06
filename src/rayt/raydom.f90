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

!> \file raydom.f90
!> \brief Main subroutine for solving of the radiative transfer equation.
!>
!> Two types of method are available:
!> - Discretes Ordinates Methods (DOM)
!> - P-1 approximation (only recommended for the CP)
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     itypfb        boundary face types
!> \param[in]     izfrad        zone number of boundary faces
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine raydom &
 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   izfrad ,                                                       &
   dt )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use cs_fuel_incl
use ppincl
use cpincl
use radiat
use ihmpre
use dimens, only: ndimfb
use mesh
use field
use cs_c_bindings
use pointe, only:pmapper_double_r1

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          itypfb(ndimfb)
integer          izfrad(ndimfb)

double precision dt(ncelet)

! Local variables

integer          iappel
integer          ifac   , iel    , iok    , izone  , isou   , jsou
integer          inc    , iccocg , iwarnp , imligp , nswrgp
integer          icla   , ipcla  , f_id0
integer          idverl
integer          iflux(nozrdm)
integer          ngg, i
double precision epsrgp, climgp, extrap
double precision aa, bb, ckmin, unspi, xlimit, flunmn
double precision flux(nozrdm)
double precision vv, sf, xlc, xkmin, pp, xptk

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbrs, rovsdt
double precision, allocatable, dimension(:) :: dcp
double precision, allocatable, dimension(:) :: ckmel
double precision, allocatable, dimension(:,:,:) :: grad
double precision, allocatable, dimension(:,:) :: tempk
double precision, allocatable, dimension(:,:) :: q, coefaq
double precision, allocatable, dimension(:,:,:) :: coefbq
double precision, allocatable, dimension(:) :: coefap, coefbp
double precision, allocatable, dimension(:) :: cofafp, cofbfp
double precision, allocatable, dimension(:) :: flurds, flurdb
double precision, allocatable, dimension(:,:) :: kgi,agi,agbi,iabparh2,iempexh2,iempimh2
double precision, allocatable, dimension(:)   :: iqxpar,iqypar,iqzpar
double precision, allocatable, dimension(:)   :: iabgaz,iabpar,iemgex,iempex,ilutot
double precision, allocatable, dimension(:)   :: iqpato,iemgim,iempim,tparo

double precision, dimension(:), pointer :: bqinci, b_temp
double precision, dimension(:), pointer :: bxlam, bepa, beps, bfnet
double precision, dimension(:), pointer :: cvara_scalt
double precision, dimension(:), pointer :: cvara_yfol
double precision, dimension(:), pointer :: cpro_cp
double precision, dimension(:,:), pointer :: bqinsp
double precision, dimension(:), pointer :: cpro_qx, cpro_qy, cpro_qz
double precision, dimension(:), pointer :: cpro_cak1, cpro_tsri1, cpro_tsre1
double precision, dimension(:), pointer :: cpro_abso1, cpro_emi1, cpro_temp2
double precision, dimension(:), pointer :: cpro_cak, cpro_tsri, cpro_tsre
double precision, dimension(:), pointer :: cpro_abso, cpro_emi
double precision, dimension(:), pointer :: cpro_x2, cpro_lumin

integer    ipadom
data       ipadom /0/
save       ipadom

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================
!---> Number of passes
ipadom = ipadom + 1

! Allocate temporary arrays for the radiative equations resolution
allocate(viscf(nfac), viscb(ndimfb))
allocate(smbrs(ncelet), rovsdt(ncelet))

! Allocate specific arrays for the radiative transfer module
allocate(tempk(ncelet,nrphas))
allocate(coefap(ndimfb), coefbp(ndimfb))
allocate(cofafp(ndimfb), cofbfp(ndimfb))
allocate(flurds(nfac), flurdb(ndimfb))

! Allocate work arrays
allocate(ckmel(ncelet)) ! Absorption coeffcient of the bulk phase
allocate(dcp(ncelet))   ! Specific heat capacity of the bulk phase
allocate(tparo(nfabor))

! Map field arrays
call field_get_val_s(itempb,b_temp)     ! Boundary temperature
call field_get_val_s(iqinci, bqinci)    ! Irradiating flux density
call field_get_val_s(ixlam, bxlam)      ! Heat conduction coeffcient of walls
call field_get_val_s(iepa, bepa)        ! Thickness of walls
call field_get_val_s(ieps, beps)        ! Emmisivity of walls
call field_get_val_s(ifnet, bfnet)      ! Radiosity at walls

if (icp.gt.0) then
  call field_get_val_s(iprpfl(icp), cpro_cp)
endif
! ADF model parameters
if (imoadf.ge.1.or.imfsck.eq.1) then
  call field_get_val_v(iqinsp,bqinsp)   ! Irradiating spectral flux density
endif

allocate(kgi(ncelet,nwsgg),agi(ncelet,nwsgg))    ! Radiation coeffcient kgi and
! the corresponding weight agi of the i-th grey gas
allocate(iqxpar(ncelet),iqypar(ncelet),iqzpar(ncelet)) ! Flux density components
allocate(iabgaz(ncelet),iabpar(ncelet)) ! Radiation absorbed by
! the gasphase and the solid phase (all particles classes)
allocate(iemgex(ncelet),iempex(ncelet)) ! Emmitted radtion of the gasphase and
! the solid phase (all particle classes)
allocate(iabparh2(ncelet,nclacp),iempexh2(ncelet,nclacp)) ! Absorbed and emmitted
! radiation of a single size class (needed to calculate the source terms of the
! particle enthalpy equation)
allocate(iemgim(ncelet),iempim(ncelet)) ! Implicit source terms of the bulk phase
! enthalpie equation
allocate(iempimh2(ncelet,nclacp))       ! Implicit source term of the particle
! enthalpie equation
allocate(ilutot(ncelet))                ! Total emitted intensity
allocate(iqpato(nfabor))                ! Irradiating  flux density at walls
! Careful: Should not be mixed up with bqinci
allocate(agbi(nfabor,nwsgg))            ! Weight of the i-th grey gas at walls

! Wall temperature

if (itpscl.eq.2) then
  xptk = tkelvi
else
  xptk =0
endif

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then
    tparo(ifac) = b_temp(ifac) + xptk
  else
    tparo(ifac) = 0.d0
  endif
enddo

! FSCK model parameters
if (ipadom.eq.1) then
  ! Weight of the i-the gaussian quadrature
  allocate(wq(nwsgg))
endif

!===============================================================================
! 1. Initializations
!===============================================================================

if (ipadom.gt.1 .and. mod(ntcabs,nfreqr).ne.0) return

write(nfecra,1000)

!---> Constants initialization
unspi = 1.d0/pi

call field_get_val_s(iprpfl(icak(1)),cpro_cak1)
call field_get_val_s(iprpfl(itsri(1)),cpro_tsri1)
call field_get_val_s(iprpfl(itsre(1)),cpro_tsre1)
call field_get_val_s(iprpfl(iabso(1)),cpro_abso1)
call field_get_val_s(iprpfl(iemi(1)),cpro_emi1)
call field_get_val_s(iprpfl(iqx),cpro_qx)
call field_get_val_s(iprpfl(iqy),cpro_qy)
call field_get_val_s(iprpfl(iqz),cpro_qz)
call field_get_val_s(iprpfl(ilumin),cpro_lumin)

!--> Working arrays
do iel = 1, ncel
  cpro_cak1(iel)  = 0.d0  ! Radiation coefficient k of the gas phase
  cpro_tsri1(iel) = 0.d0  ! TS implicit due to emission
  cpro_tsre1(iel) = 0.d0  ! TS explicit due to emission and absorption
  cpro_abso1(iel) = 0.d0  ! Absortion: Sum,i((kg,i+kp) * Integral(Ii)dOmega)
  cpro_emi1(iel)  = 0.d0  ! Emission:  Sum,i((kg,i+kp) * stephn * T^4 *agi)
!
  cpro_qx(iel) = 0.d0 ! X-component of the radiative flux vector
  cpro_qy(iel) = 0.d0 ! Y-compnent of the radiative flux vector
  cpro_qz(iel) = 0.d0 ! Z-Component of the radiative flux vector

  iabgaz(iel) = 0.d0 ! Absorption of the gas phase: kg,i * Integral(Ii) * dOmega
  iabpar(iel) = 0.d0 ! Absortion of particles: kp * Integral(Ii) * dOmega
  iemgex(iel) = 0.d0 ! Emission of the gas phase: kg,i * stephn * T^4 *agi
  iempex(iel) = 0.d0 ! Emission of particles: kp * stephn * T^4 *agi
  iemgim(iel) = 0.d0 ! Gas phase related implicit source term in the bulk phase enthalpy eqn.
  iempim(iel) = 0.d0 ! Particle related implicit source term in the bulk phase enthalpy eqn.
  ilutot(iel) = 0.d0 ! Total emitted intensity
  ckmel(iel)  = 0.d0 ! Radiation coeffcient of the bulk phase

  if (icp.gt.0) then
    dcp(iel) = 1.d0/cpro_cp(iel)
  else
    dcp(iel) = 1.d0/cp0
  endif

  do i=1,nwsgg
    kgi(iel,i)= 0.d0
    agi(iel,i)= 1.d0 ! In case of grey gas radiation properties (kgi!=f(lambda))
                     ! agi must be set to 1.
  enddo
enddo

do ifac = 1, nfabor
  iqpato(ifac) = 0.d0
  do i = 1, nwsgg
    agbi(ifac,i) = 1.d0 ! In case of grey gas radiation properties (kgi!=f(lambda))
                        ! agbi must be set to 1.
  enddo
enddo

do iel = 1,ncel
  do icla = 1, nclacp
    iabparh2(iel,icla) = 0.d0
    iempexh2(iel,icla) = 0.d0
    iempimh2(iel,icla) = 0.d0
  enddo
enddo

if (ipadom.eq.1) then
  do i = 1,nwsgg
    wq(i) = 1.d0      ! Must be set to 1 in case of using the standard as well as
                      ! the ADF radiation models
  enddo
endif
!=============================================================================
! 2. Temperature storing (in Kelvin) in tempk(iel, irphas)
!=============================================================================

!---> Temperature transport
if (itherm.eq.1) then

  call field_get_val_prev_s(ivarfl(isca(iscalt)), cvara_scalt)

  ! iscalt is in Celsius
  if (itpscl.eq.2) then
    do iel = 1, ncel
      tempk(iel,1) = cvara_scalt(iel) + tkelvi
    enddo
  else
    do iel = 1, ncel
      tempk(iel,1) = cvara_scalt(iel)
    enddo
  endif

!---> Enthalpy transport (flurdb is a temporary array)
else if (itherm.eq.2) then

  call field_get_val_prev_s(ivarfl(isca(iscalt)), cvara_scalt)

  call c_h_to_t(cvara_scalt, tempk(:,1))

  if (ippmod(iccoal).ge.0) then

    ! Particules' temperature
    do icla = 1, nclacp
      call field_get_val_s(iprpfl(itemp2(icla)),cpro_temp2)
      ipcla = 1+icla
      do iel = 1, ncel
        tempk(iel,ipcla) = cpro_temp2(iel)
      enddo
    enddo

  ! Fuel
  else if (ippmod(icfuel).ge.0) then

    do icla = 1, nclafu
      call field_get_val_s(iprpfl(itemp2(icla)),cpro_temp2)
      ipcla = 1+icla
      do iel = 1, ncel
        tempk(iel,ipcla) = cpro_temp2(iel)
      enddo
    enddo

  endif

else
  write(nfecra,3500) itherm
  call csexit (1)
endif

!=============================================================================
! 3. Absorption coefficient
!=============================================================================

!--> Initialization to a non-admissible value for testing after usray3
do iel = 1, ncel
  cpro_cak1(iel) = -grand
enddo

!--> Absorption coefficient for different modules

! Warning: for the approximation P-1, the absorption coefficient is required
!          for boundary conditions

if (ippmod(iphpar).ge.2) then

  call ppcabs(tempk, kgi, agi, agbi)

else

  !---> Reading of User datas

  !   - Interface Code_Saturne
  !     ======================

  if (iihmpr.eq.1) then

    call uiray3(cpro_cak1, ncel, imodak) !FIXME for ADF

    if (iirayo.eq.2 .and. ippmod(iphpar).le.1 .and. ipadom.le.3) then
      sf = 0.d0
      vv = 0.d0

      ! Compute the characteristic length of the computational domain
      do ifac = 1, nfabor
        sf = sf + sqrt(surfbo(1,ifac)**2 +                      &
                       surfbo(2,ifac)**2 +                      &
                       surfbo(3,ifac)**2 )
      enddo
      if (irangp.ge.0) then
        call parsom(sf)
      endif

      do iel = 1, ncel
        vv = vv + volume(iel)
      enddo
      if (irangp.ge.0) then
        call parsom(vv)
      endif

      xlc = 3.6d0 * vv / sf

      !  Clipping on ck
      xkmin = 1.d0 / xlc

      iok = 0
      do iel = 1, ncel
        if (cpro_cak1(iel).lt.xkmin) then
          iok = iok +1
        endif
      enddo

      ! Warning if the optical thickness is too big
      pp = xnp1mx/100.0d0
      if (dble(iok).gt.pp*dble(ncel)) then
        write(nfecra,6000) xkmin, dble(iok)/dble(ncel)*100.d0,  &
                           xnp1mx
      endif
    endif

  endif

  ! Only necessary when grey gas radiation properties are applied.
  ! In case of the ADF model this test doesnt make sence.

  if (imoadf.eq.0.and.imfsck.eq.0) then
    call usray3 &
  ( nvar   , nscal  , iappel ,                                     &
    itypfb ,                                                       &
    izfrad ,                                                       &
    dt     ,                                                       &
    cpro_cak1)
  endif

endif

!--> General checking

!--> Test if the radiation coeffcient has been assigned
if (iirayo.ge.1) then
  if (imoadf.eq.0.and.imfsck.eq.0) then
    ckmin = cpro_cak1(1)
    do iel = 1, ncel
      ckmin = min(ckmin,cpro_cak1(iel))
    enddo

    if (irangp.ge.0) then
      call parmin(ckmin)
    endif

    if (ckmin.lt.0.d0) then
      if (iirayo.eq.2) then
        write(nfecra,2020)
      else if (iirayo.eq.1) then
        write(nfecra,2010)
      endif
      call csexit (1)
    endif
  else
    ckmin = 0.d0
    do iel = 1, ncel
      do ngg = 1, nwsgg
        ckmin = min(ckmin, kgi(iel,ngg))
      enddo
    enddo

    if (irangp.ge.0) then
      call parmin(ckmin)
    endif

    if (ckmin.lt.0.d0) then
      if (iirayo.eq.2) then
        write(nfecra,2020)
      else if (iirayo.eq.1) then
        write(nfecra,2010)
      endif
      call csexit (1)
    endif
  endif
endif

!---> Check of a transparent case
idverl = idiver
!=============================================================================
! 4. Solving the ETR
!=============================================================================
! Loop over all grey gases. Remember: In case of the basic radiation models of
! Code_Saturne nwsgg=1

do ngg = 1, nwsgg

  if (imoadf.ge.1.or.imfsck.eq.1) then
    do iel = 1, ncel
      cpro_cak1(iel) = kgi(iel, ngg) ! TODO merge the two arrays
    enddo

  else
    aa = 0.d0
    do iel = 1, ncel
      aa = max(aa, cpro_cak1(iel))
    enddo
    if (irangp.ge.0) then
      call parmax(aa)
    endif
    if (aa.le.epzero) then
      write(nfecra,1100)
      idverl = -1
    endif
  endif

!===============================================================================
! 4.1 Radiative P-1 model
!===============================================================================

  if (iirayo.eq.2) then

    !--> Gas phase: Explicit source term in the transport eqn. of theta4

    do iel = 1, ncel
      smbrs(iel) = 3.d0*cpro_cak1(iel)*(tempk(iel,1)**4)          &
                 * agi(iel, ngg)*volume(iel)
    enddo

    !--> Solid phase:

    ! Coal particles: Explicit source term in the transport eqn. of theta4
    if (ippmod(iccoal).ge.0) then
      do icla = 1, nclacp
        ipcla = 1+icla
        call field_get_val_s(iprpfl(ix2(icla)),cpro_x2)
        call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
        do iel = 1,ncel
          smbrs(iel) = smbrs(iel)                               &
                     + (3.d0*cpro_x2(iel)                       &
                       *  cpro_cak(iel)                         &
                       * (tempk(iel,ipcla)**4)*agi(iel,ngg)     &
                       *  volume(iel))
        enddo
      enddo

    ! Fuel droplets: Explicit source term in the transport eqn. of theta4
    else if (ippmod(icfuel).ge.0) then
      do icla = 1, nclafu
        ipcla = 1+icla
        call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
        call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
        do iel = 1,ncel
          smbrs(iel) =  smbrs(iel)                              &
                     + (3.d0*cvara_yfol(iel)                    &
                       *  cpro_cak(iel)                         &
                       * (tempk(iel,ipcla)**4)*agi(iel,ngg)     &
                       *  volume(iel) )
        enddo
      enddo
    endif

    !--> Gas phase: Implicit source term in the transport eqn. of theta4

    do iel = 1, ncel
      rovsdt(iel) =  3.d0*cpro_cak1(iel)*volume(iel)
    enddo

    !--> Solid phase:

    ! Coal particles: Implicit source term in the transport eqn. of theta4
    if (ippmod(iccoal).ge.0) then
      do icla = 1, nclacp
        ipcla = 1+icla
        call field_get_val_s(iprpfl(ix2(icla)),cpro_x2)
        call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
        do iel = 1,ncel
          rovsdt(iel) = rovsdt(iel)                                      &
                      + (3.d0*cpro_x2(iel)                               &
                        * cpro_cak(iel) * volume(iel) )
        enddo
      enddo

    ! Fuel droplets: Implicit source term in the transport eqn. of theta4
    else if (ippmod(icfuel).ge.0) then
      do icla = 1, nclafu
        ipcla = 1+icla
        call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
        call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
        do iel = 1,ncel
          rovsdt(iel) = rovsdt(iel)                                      &
                      + (3.d0*cvara_yfol(iel)                            &
                        * cpro_cak(iel) * volume(iel) )
        enddo
      enddo

    endif

    ! Radiation coeffcient of the bulk phase
    ! Gas phase:
    do iel = 1, ncel
      ckmel(iel) = cpro_cak1(iel)
    enddo

    if (ippmod(iccoal).ge.0) then
      ! Solid phase:
      ! Coal particles
      do icla = 1, nclacp
        ipcla = 1+icla
        call field_get_val_s(iprpfl(ix2(icla)),cpro_x2)
        call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
        do iel = 1, ncel
          ckmel(iel) = ckmel(iel)                                  &
                     + ( cpro_x2(iel)                              &
                       * cpro_cak(iel) )
        enddo
      enddo
    ! Fuel droplets
    else if (ippmod(icfuel).ge.0) then
      do icla = 1, nclafu
        ipcla = 1+icla
        call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
        call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
        do iel = 1, ncel
          ckmel(iel) = ckmel(iel)                                  &
                     + ( cvara_yfol(iel)                           &
                       * cpro_cak(iel) )
        enddo
      enddo
    endif

    ! Test if ckmel is gt zero
    do iel = 1, ncel
      if (ckmel(iel).le.0.d0) then
        write(nfecra,7000)
        call csexit (1)
      endif
    enddo

    ! Update Boundary condition coefficients

    call raycll &
    ( itypfb ,                                                       &
      izfrad ,                                                       &
      coefap , coefbp ,                                              &
      cofafp , cofbfp ,                                              &
      tparo  , bqinci , beps   ,                                     &
      ckmel, agbi, ngg )

    ! Solving

    call raypun &
    ( itypfb ,                                                       &
      coefap , coefbp ,                                              &
      cofafp , cofbfp ,                                              &
      flurds , flurdb ,                                              &
      viscf  , viscb  ,                                              &
      smbrs  , rovsdt ,                                              &
      cpro_abso1 , cpro_emi1 ,                                       &
      cpro_tsre1 ,                                                   &
      iqxpar , iqypar , iqzpar ,                                     &
      bqinci , beps   , tparo  ,                                     &
      ckmel  , agbi   , ngg    )

  !===============================================================================
  ! 4.2 Solving of the radiative transfert equation (DOM)
  !===============================================================================

  else if (iirayo.eq.1) then

    !--> Gas phase: Explicit source term of the ETR
    do iel = 1, ncel
      smbrs(iel) = stephn*cpro_cak1(iel)              &
                 *(tempk(iel,1)**4)*agi(iel,ngg)*volume(iel)*unspi
    enddo

    !--> Solid phase:
    ! Coal particles: Explicit source term of the ETR
    if (ippmod(iccoal).ge.0) then
      do icla = 1, nclacp
        ipcla = 1+icla
        call field_get_val_s(iprpfl(ix2(icla)),cpro_x2)
        call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
        do iel = 1, ncel
          smbrs(iel) = smbrs(iel)                                 &
                     + cpro_x2(iel)                               &
                       * agi(iel,ngg)*stephn                      &
                       *cpro_cak(iel)                             &
                       *(tempk(iel,ipcla)**4)                     &
                       * volume(iel)*unspi
        enddo
      enddo
    ! Fuel droplets: Explicit source term of the ETR
    elseif (ippmod(icfuel).ge.0) then
      do icla = 1,nclafu
        ipcla = 1+icla
        call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
        call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
        do iel = 1,ncel
          smbrs(iel) = smbrs(iel)                         &
                     + cvara_yfol(iel)                    &
                       *agi(iel,ngg)*stephn               &
                       *cpro_cak(iel)                     &
                       *(tempk(iel,ipcla)**4)             &
                       *volume(iel)* unspi
        enddo
      enddo
    endif

    !--> Gas phase: Implicit source term of the ETR
    do iel = 1, ncel
      rovsdt(iel) = cpro_cak1(iel) * volume(iel)
    enddo

    !--> Solid phase
    ! Coal particles: Implicit source term of the ETR
    if (ippmod(iccoal).ge.0) then
      do icla = 1, nclacp
        ipcla = 1+icla
        call field_get_val_s(iprpfl(ix2(icla)),cpro_x2)
        call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel)                                      &
                      + cpro_x2(iel)                                     &
                        * cpro_cak(iel) * volume(iel)
        enddo
      enddo
    ! Fuel droplets: Implicit source term of the ETR
    elseif (ippmod(icfuel).ge.0) then
      do icla = 1, nclafu
        ipcla = 1+icla
        call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
        call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel)                                      &
                      + cvara_yfol(iel)                                  &
                        * cpro_cak(iel) * volume(iel)
        enddo
      enddo
    endif

    ! Update Boundary condiction coefficients

    call raycll &
    ( itypfb ,                                                       &
      izfrad ,                                                       &
      coefap , coefbp ,                                              &
      cofafp , cofbfp ,                                              &
      tparo  , bqinci , beps   ,                                     &
      ckmel, agbi, ngg )

    ! Solving

    call raysol &
    ( coefap , coefbp ,                                              &
      cofafp , cofbfp ,                                              &
      flurds , flurdb ,                                              &
      viscf  , viscb  ,                                              &
      smbrs  , rovsdt ,                                              &
      cpro_tsre1 ,                                                 &
      iqxpar , iqypar , iqzpar  ,                                    &
      bqinci , bfnet  , ngg)

  endif

  ! Summing up the quantities of each grey gas
  do iel = 1, ncel

    ! Absorption
    iabgaz(iel) = iabgaz(iel)+(cpro_cak1(iel) *          &
                               cpro_tsre1(iel) *         &
                               wq(ngg))
  enddo

  if (ippmod(iccoal).ge.0) then
    do icla = 1, nclacp
      ipcla = 1+icla
      call field_get_val_s(iprpfl(ix2(icla)),cpro_x2)
      call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
      do iel = 1, ncel
        iabpar(iel) = iabpar(iel) + (cpro_x2(iel)                   &
                                  * cpro_cak(iel)                   &
                                  * cpro_tsre1(iel)                  &
                                  * wq(ngg))
        iabparh2(iel,icla)        = iabparh2(iel,icla)              &
                                  + (cpro_cak(iel)                  &
                                  * cpro_tsre1(iel)                  &
                                  * wq(ngg))
      enddo
    enddo
  else if (ippmod(icfuel).ge.0) then
    do icla = 1, nclafu
      ipcla = 1+icla
      call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
      call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
      do iel = 1, ncel
        iabpar(iel) = iabpar(iel) + (cvara_yfol(iel)                &
                                  * cpro_cak(iel)                   &
                                  * cpro_tsre1(iel)                  &
                                  * wq(ngg))
        iabparh2(iel,icla)        = iabparh2(iel,icla)              &
                                  + (cpro_cak(iel)                  &
                                  * cpro_tsre1(iel)                  &
                                  * wq(ngg))
      enddo
    enddo
  endif

  ! Emission
  do iel = 1, ncel
    iemgex(iel)=iemgex(iel)-(cpro_cak1(iel)                          &
                           *agi(iel,ngg)*4.d0*stephn                &
                           *(tempk(iel,1)**4)*wq(ngg))
    iemgim(iel)=iemgim(iel)-(16.d0*dcp(iel)*cpro_cak1(iel)           &
                           * agi(iel,ngg)*stephn*(tempk(iel,1)**3)*wq(ngg))

  enddo
  if (ippmod(iccoal).ge.0) then
    do icla = 1, nclacp
      ipcla = 1+icla
      call field_get_val_s(iprpfl(ix2(icla)),cpro_x2)
      call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
      do iel = 1, ncel
        iempex(iel) = iempex(iel) -(4.0d0*cpro_x2(iel)                       &
                                  *stephn * cpro_cak(iel)                    &
                                  *(tempk(iel,ipcla)**4)*agi(iel,ngg)*wq(ngg))
        iempexh2(iel,icla) = iempexh2(iel,icla)                              &
                             -(4.0d0*stephn * cpro_cak(iel)                  &
                             *(tempk(iel,ipcla)**4)*agi(iel,ngg)*wq(ngg))
        iempim(iel) = iempim(iel) -(16.d0*cpro_cak(iel)                      &
                                  *cpro_x2(iel)                              &
                                  *stephn*(tempk(iel,ipcla)**3)*agi(iel,ngg) &
                                  /cp2ch(ichcor(icla))*wq(ngg))
        iempimh2(iel,icla) = iempimh2(iel,icla)                              &
                             -(16.d0*cpro_cak(iel)                           &
                             *stephn*(tempk(iel,ipcla)**3)*agi(iel,ngg)      &
                             /cp2ch(ichcor(icla))*wq(ngg))

      enddo
    enddo
  else if (ippmod(icfuel).ge.0) then
    do icla = 1, nclafu
      ipcla = 1+icla
      call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
      call field_get_val_s(iprpfl(icak(ipcla)),cpro_cak)
      do iel = 1, ncel
        iempex(iel) = iempex(iel)  -(4.0d0*cvara_yfol(iel)                    &
                                   *stephn * cpro_cak(iel)                    &
                                   *(tempk(iel,ipcla)**4)*agi(iel,ngg)*wq(ngg))
        iempexh2(iel,icla) = iempexh2(iel,icla)                               &
                             -(4.0d0*stephn * cpro_cak(iel)                   &
                             *(tempk(iel,ipcla)**4)*agi(iel,ngg)*wq(ngg))
        iempim(iel) = iempim(iel) -(16.d0*cpro_cak(iel)                       &
                                  *cvara_yfol(iel)*stephn                     &
                                  *(tempk(iel,ipcla)**3)*agi(iel,ngg)/cp2fol  &
                                  *wq(ngg))
        iempimh2(iel,icla) = iempimh2(iel,icla)                               &
                             -(16.d0*cpro_cak(iel)                            &
                             *stephn*(tempk(iel,ipcla)**3)*agi(iel,ngg)/cp2fol&
                             *wq(ngg))
      enddo
    enddo
  endif

  do iel = 1, ncel
    ! Emitted intensity
    ilutot(iel)=ilutot(iel)+(cpro_tsre1(iel)*wq(ngg))

    ! Flux vector components
    cpro_qx(iel) = cpro_qx(iel)+(iqxpar(iel)*wq(ngg))
    cpro_qy(iel) = cpro_qy(iel)+(iqypar(iel)*wq(ngg))
    cpro_qz(iel) = cpro_qz(iel)+(iqzpar(iel)*wq(ngg))
  enddo

  ! If the ADF model is activated we have to sum up the spectral flux densities
  if (imoadf.ge.1) then
    do ifac =1, nfabor
      iqpato(ifac)=iqpato(ifac)+(bqinsp(ngg,ifac)*wq(ngg))
    enddo
  endif
enddo

!The total radiative flux is copied in bqinci
!a) for post-processing reasons and
!b) in order to calculate bfnet
if (imoadf.ge.1) then
  do ifac=1,nfabor
    bqinci(ifac) = iqpato(ifac)
  enddo
endif

!===============================================================================
! 5.  Storing of the total emitted intensity
!===============================================================================
!                             /    ->  ->
!                        SA= /  L( X , S ). DOMEGA
!                           /4.PI

do iel=1,ncel
  cpro_lumin(iel) = ilutot(iel)
enddo

!===============================================================================
! 6. Net radiative flux at walls: computation and integration
!===============================================================================

!--> Initialization to a non-admissible value for testing after usray5
do ifac = 1,nfabor
  bfnet(ifac) = -grand
enddo

!---> Reading of User datas
!CAREFUL: The user has access to the radiation coeffcient (array cpro_cak1)
!in usray5. However, only when the standard radiation models of code_saturne are
!applied, this table contains the true value given by the user. Thus, the usage
!of the radiation coeffcient in usray5 must be done carefully. In its present
!version usray5 does NOT use the radiation coeffcient, and thus, usray5 can still
!be called here, even if the ADF model is activated.

call usray5 &
( nvar   , nscal  ,                                              &
  itypfb ,                                                       &
  izfrad ,                                                       &
  dt     ,                                                       &
  coefap , coefbp ,                                              &
  cofafp , cofbfp ,                                              &
  tparo  , bqinci ,                                              &
  bfnet  , bxlam  , bepa   , beps   ,                            &
  cpro_cak1 )

!---> Check flunet
iok = 0
xlimit = -grand*0.1d0
flunmn = grand

do ifac = 1, nfabor
  if (bfnet(ifac).le.xlimit) then
    iok = iok + 1
    flunmn = min(flunmn,bfnet(ifac))
    if (irangp.ge.0) then
      call parmin(flunmn)
    endif
    write(nfecra,4000)ifac,izfrad(ifac),itypfb(ifac)
  endif
enddo

if (iok.ne.0) then
  write(nfecra,4100) flunmn
  call csexit (1)
endif

!--> Integration du flux net sur les differentes zones de frontieres
!     IFLUX sert en parallele pour reperer les zones existantes

do izone = 1, nozrdm
  flux(izone) = 0.d0
  iflux(izone) = 0
enddo
do ifac = 1,nfabor
  izone = izfrad(ifac)
  flux(izone) = flux(izone) + bfnet(ifac)*surfbn(ifac)
  iflux(izone) = 1
enddo
if(irangp.ge.0) then
  call parrsm(nozarm,flux )
  call parimx(nozarm,iflux)
endif


write(nfecra,5000)
write(nfecra,5010)
do izone = 1, nozarm
  if(iflux(izone).eq.1) then
    write(nfecra,5020) izone,flux(izone)
  endif
enddo
write(nfecra,5000)


!--> Integration de la densite de flux net aux frontieres

aa = 0.d0
do ifac = 1,nfabor
  aa =  aa + bfnet(ifac) * surfbn(ifac)
enddo
if (irangp.ge.0) then
  call parsom(aa)
endif
write(nfecra,5030) aa

!===============================================================================
! 7. Implicit and explicit radiative source terms
!===============================================================================


!===============================================================================
! 7.1 Semi-analitical radiative source termes
!===============================================================================

if (idverl.ge.0) then

  do iel = 1, ncel
    ! Absorption of the gas is copied into abso1
    cpro_abso1(iel) = iabgaz(iel)
    ! Emission of the gas phase is copied into emi1
    cpro_emi1(iel) = iemgex(iel)
  enddo
  if (ippmod(iccoal).ge.0.or.ippmod(icfuel).ge.0) then
    do iel = 1, ncel
      ! Absoprtion of particles is added to abso1
      cpro_abso1(iel) = cpro_abso1(iel) + iabpar(iel)
      ! Emission of particles is added to emi1
      cpro_emi1(iel) = cpro_emi1(iel) + iempex(iel)
    enddo
  endif
  do iel = 1, ncel
    ! Emission + Absorption of gas and particles --> TSexplicit
    cpro_tsre1(iel) = cpro_abso1(iel) + cpro_emi1(iel)
  enddo

  do iel = 1, ncel
    ! TSimplicit of the gas phase
    cpro_tsri1(iel) = iemgim(iel)
  enddo
  if (ippmod(iccoal).ge.0.or.ippmod(icfuel).ge.0) then
    do iel = 1, ncel
      ! TSimplicit of the solid phase is added to stri1
      cpro_tsri1(iel) = cpro_tsri1(iel) + iempim(iel)
    enddo
  endif

  ! In order to determine the source terms of the particle enthalpy tranport eqn.,
  ! we have to copie the approriate determined aboce into the corressponding tables
  if (ippmod(iccoal).ge.0.or.ippmod(icfuel).ge.0) then
    do icla = 1,nclacp
      ipcla = 1+icla
      call field_get_val_s(iprpfl(itsri(ipcla)),cpro_tsri)
      call field_get_val_s(iprpfl(itsre(ipcla)),cpro_tsre)
      call field_get_val_s(iprpfl(iabso(ipcla)),cpro_abso)
      call field_get_val_s(iprpfl(iemi(ipcla)),cpro_emi)
      do iel = 1, ncel
        cpro_abso(iel) = iabparh2(iel,icla)
        cpro_emi(iel)  = iempexh2(iel,icla)
        cpro_tsre(iel) = iabparh2(iel,icla)+iempexh2(iel,icla)
        cpro_tsri(iel) = iempimh2(iel,icla)
      enddo
    enddo
  endif
else
  do iel = 1, ncel
    cpro_abso1(iel)  = 0.d0
    cpro_emi1(iel)  = 0.d0
    cpro_tsre1(iel) = 0.d0
    cpro_tsri1(iel) = 0.d0
  enddo
endif

!===============================================================================
! 6.2 Explicit conservative radiative source terms
!===============================================================================

! coefap and coefbp are NOW Boundary conditions on the divergence

if (idverl.eq.1 .or. idverl.eq.2) then

  ! Allocate  temporary arrays for gradient computation

  allocate(q(3,ncelet))

  do iel = 1, ncel
    q(1,iel) = cpro_qx(iel)
    q(2,iel) = cpro_qy(iel)
    q(3,iel) = cpro_qz(iel)
  enddo

  allocate(coefaq(3,nfabor))
  allocate(coefbq(3,3,nfabor))

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    coefaq(1,ifac) = bfnet(ifac)*surfbo(1,ifac) / surfbn(ifac)
    coefaq(2,ifac) = bfnet(ifac)*surfbo(2,ifac) / surfbn(ifac)
    coefaq(3,ifac) = bfnet(ifac)*surfbo(3,ifac) / surfbn(ifac)
  enddo

  do ifac = 1, nfabor
    do isou = 1, 3
      do jsou = 1, 3
        coefbq(isou,jsou,ifac) = zero
      enddo
    enddo
  enddo

  allocate(grad(3,3,ncelet))

  ! Donnees pour le calcul de la divergence
  inc     = 1
  iccocg  = 1
  imligp  = -1
  iwarnp  = iimlum
  epsrgp  = 1.d-8
  climgp  = 1.5d0
  extrap  = 0.d0
  nswrgp  = 100

  f_id0 = -1

  call cgdvec                                                     &
  ( f_id0  , imrgra , inc    , nswrgp , iwarnp , imligp ,         &
    epsrgp , climgp , coefaq , coefbq , q      , grad   )

  do iel = 1,ncel
    cpro_tsre1(iel) = - grad(1,1,iel)                              &
                      - grad(2,2,iel)                              &
                      - grad(3,3,iel)
  enddo

  ! Free memory
  deallocate(grad)
  deallocate(coefbq)
  deallocate(coefaq)
  deallocate(q)

! Fin du calcul de la divergence
endif

!===============================================================================
! 7.3 Explicit radiative semi-analytical corrected source term
!===============================================================================

if (idverl.eq.2) then

  !---> Comparison of the semi-analytical and conservative source terms
  aa = 0.d0
  do iel = 1, ncel
    aa = aa + cpro_tsre1(iel) * volume(iel)
  enddo

  bb = 0.d0
  do iel = 1,ncel
    bb = bb + (cpro_abso1(iel)+cpro_emi1(iel))*volume(iel)
  enddo

  if(irangp.ge.0) then
    call parsom(aa)
    call parsom(bb)
  endif

  aa = aa/bb

  !---> Correction of the semi-analytical source term by the conservative source
  ! term
  do iel = 1,ncel
    cpro_tsre1(iel) = ( cpro_abso1(iel)       &
                      + cpro_emi1(iel))       &
                       *aa
  enddo

endif

!===============================================================================
! 7.4 Finalization of explicit source terms
!===============================================================================

if (idverl.ge.0) then

  !--> Integration volumique du terme source explicite
  !    Le resultat de cette integration DOIT etre le meme que l'integration
  !    surfacique de la densite de flux net radiatif faite plus haut
  !    si  IDVERL = 1 ou 2

  aa = 0.d0
  do iel = 1, ncel
    aa = aa + cpro_tsre1(iel) * volume(iel)
  enddo

  if (irangp.ge.0) then
    call parsom(aa)
  endif

  write(nfecra,5040) aa
  write(nfecra,5050)
  write(nfecra,5000)
!--> Correction du terme source explicite dans raysca pour permettre un
!    post-processing correct du terme source explicite
!    lorsque la variable transportee est la temperature
!    (pour les calculs en combustion la variable transportee est toujours
!    l'enthalpie)
else
  write(nfecra,5000)
endif

! Free memory
deallocate(viscf, viscb)
deallocate(smbrs, rovsdt)
deallocate(tempk)
deallocate(coefap, coefbp)
deallocate(cofafp, cofbfp)
deallocate(flurds, flurdb)
deallocate(ckmel, dcp, tparo)
deallocate(kgi,agi,agbi)
deallocate(iqxpar,iqypar,iqzpar)
deallocate(iabgaz,iabpar,iemgex,iempex,ilutot)
deallocate(iemgim,iempim)
deallocate(iabparh2,iempexh2,iempimh2)
if (ntcabs.eq.ntmabs) deallocate(wq)
!--------
! Formats
!--------

 1000 FORMAT (/, 3X,'** INFORMATIONS SUR LE TERME SOURCE RADIATIF',/,   &
           3X,'   -----------------------------------------' )
 1100 FORMAT (/, 3X,'   Calcul effectue en rayonnement transparent'  ,/)

 2010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE RAYONNEMENT EST ACTIVE AVEC LE MODELE DOM.           ',/,&
'@      LA VALEUR MINIMALE DU COEFFICIENT D ABSORPTION A EST  ',/,&
'@      EGALE A ', E14.5                                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE RAYONNEMENT EST ACTIVE AVEC LE MODELE P-1.           ',/,&
'@      LE COEFFICIENT D''ABSORBTION DOIT ETRE STRICTEMENT    ',/,&
'@      SUPERIEUR A ZERO.                                     ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3500 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE RAYONNEMENT EST ACTIVE.                              ',/,&
'@                                                            ',/,&
'@    Le modele ITHERM devrait valoir 1 (temperature) ou      ',/,&
'@      2 (enthalpie). On a ITHERM = ',i10                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4000 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT (FLUNET    NON RENSEIGNE)       ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )
 4100 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    LE FLUNET    N''EST PAS RENSEIGNEE POUR CERTAINES       ',/,&
'@        FACES DE BORD                                       ',/,&
'@                                                            ',/,&
'@        Valeur minimale ',E14.5                              ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usray5.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 5000 format('-----------------------------------------------------',   &
          '--------------')

 5010 format('Zone         Flux net radiatif (Watt) (normale',          &
          ' unitaire sortante)')

 5020 format(i6,13x,e11.4)

 5030 format('Flux net radiatif sur toutes les frontieres  Fnet = ',    &
           e11.4,' Watt')

 5040 format('Integrale volumique du terme source radiatif Srad = ',    &
           e11.4,' Watt')

 5050 format('(Si IDIVER = 1 ou 2 alors on doit avoir Srad = -Fnet)')

 6000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT APPROXIMATION P-1  (RAYDOM)     ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    LA LONGUEUR OPTIQUE DU MILIEU SEMI-TRANSPARENT          ',/,&
'@      DOIT AU MOINS ETRE DE L''ORDRE DE L''UNITE POUR ETRE  ',/,&
'@      DANS LE DOMAINE D''APPLICATION DE L''APPROXIMATION P-1',/,&
'@    CELA NE SEMBLE PAS ETRE LE CAS ICI.                     ',/,&
'@                                                            ',/,&
'@    LE COEFFICIENT D''ABSORPTION MINIMUM POUR ASSURER CETTE ',/,&
'@      LONGUEUR OPTIQUE EST XKMIN = ',e11.4                   ,/,&
'@    CETTE VALEUR N''EST PAS ATTEINTE POUR ', e11.4,'%       ',/,&
'@      DES CELLULES DU MAILLAGE.                             ',/,&
'@    LE POURCENTAGE DE CELLULES DU MAILLAGE POUR LESQUELLES  ',/,&
'@      ON ADMET QUE CETTE CONDITION SOIT VIOLEE EST IMPOSE   ',/,&
'@      PAR DEFAUT OU DANS USINI1 A XNP1MX = ', e11.4,'%      ',/,&
'@                                                            ',/,&
'@    Verifier les valeurs du coefficient d''absorption CK    ',/,&
'@      dans l''interface ou le modifier dans USRAY3.         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : P-1 Radiation model (Subroutine RAYDOM)     ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    The local radiation coeffcient of the bulk phase ckmel  ',/,&
'@    takes the value 0 somewhere. This often occurs during   ',/,&
'@    the very first iterations of the simulation.            ',/,&
'@    Thus, make sure the coal and/or the char mass fraction  ',/,&
'@    have been initialzed to values different from zero.     ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

end subroutine
