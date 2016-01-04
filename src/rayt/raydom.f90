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

subroutine raydom &
 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   izfrad ,                                                       &
   dt     , propce )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  Enveloppe principale du module de resolution de l'equation
!  des transferts radiatifs

!  Deux methodes sont disponibles :

!    1) La methode : "Discretes Ordinates Methods" (DOM)
!    2) L'approximation P-1 (recommande uniquement pour le CP)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfrad(nfabor    ! te ! <-- ! numero de zone des faces de bord               !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          itypfb(ndimfb)
integer          izfrad(ndimfb)

double precision dt(ncelet)
double precision propce(ncelet,*)

! Local variables

integer          iappel
integer          ifac   , iel    , iok    , izone  , isou   , jsou
integer          inc    , iccocg , iwarnp , imligp , nswrgp
integer          mode   , icla   , ipcla  , f_id0
integer          idverl
integer          iflux(nozrdm)
integer          ngg, i
double precision epsrgp, climgp, extrap
double precision aa, bb, ckmin, ckmax, unspi, xlimit, flunmn
double precision flux(nozrdm)
double precision vv, sf, xlc, xkmin, pp

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
double precision, allocatable, dimension(:)   :: iqpato,iemgim,iempim

double precision, dimension(:), pointer :: tparo, bqinci
double precision, dimension(:), pointer :: bxlam, bepa, beps, bfnet
double precision, dimension(:), pointer :: cvara_scalt
double precision, dimension(:), pointer :: cvara_yfol
double precision, dimension(:), pointer :: cpro_cp
double precision, dimension(:,:), pointer :: bqinsp

integer    ipadom
data       ipadom /0/
save       ipadom

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

! Allocate temporary arrays for the radiative equations resolution
allocate(viscf(nfac), viscb(ndimfb))
allocate(smbrs(ncelet), rovsdt(ncelet))

! Allocate specific arrays for the radiative transfert module
allocate(tempk(ncelet,nrphas))
allocate(coefap(ndimfb), coefbp(ndimfb))
allocate(cofafp(ndimfb), cofbfp(ndimfb))
allocate(flurds(nfac), flurdb(ndimfb))

! Allocate work arrays
allocate(ckmel(ncelet)) ! Absorption coeffcient of the bulk phase
allocate(dcp(ncelet))   ! Specific heat capacity of the bulk phase

! Map field arrays
call field_get_val_s(itparo, tparo)     ! Wall temperature
call field_get_val_s(iqinci, bqinci)    ! Irradiating flux density
call field_get_val_s(ixlam, bxlam)      ! Heat conduction coeffcient of walls
call field_get_val_s(iepa, bepa)        ! Thickness of walls
call field_get_val_s(ieps, beps)        ! Emmisivity of walls
call field_get_val_s(ifnet, bfnet)      ! Radiosity at walls

if (icp.gt.0) then
  call field_get_val_s(iprpfl(icp), cpro_cp)
endif
! ADF model parameters
if (imoadf.ge.1) then
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

!===============================================================================
! 1. Initializations
!===============================================================================

!---> Number of passes
ipadom = ipadom + 1
if (ipadom.gt.1 .and. mod(ntcabs,nfreqr).ne.0) return

write(nfecra,1000)

!---> Constants initialization
unspi = 1.d0/pi

!--> Working arrays
do iel = 1, ncel
  propce(iel,ipproc(icak(1)))  = 0.d0  ! Radiation coefficient k of the gas phase
  propce(iel,ipproc(itsri(1))) = 0.d0  ! TS implicit due to emission
  propce(iel,ipproc(itsre(1))) = 0.d0  ! TS explicit due to emission and absorption
  propce(iel,ipproc(iabso(1)))  = 0.d0  ! Absortion: Sum,i((kg,i+kp) * Integral(Ii)dOmega)
  propce(iel,ipproc(iemi(1)))  = 0.d0  ! Emission:  Sum,i((kg,i+kp) * stephn * T^4 *agi)
!
  propce(iel,ipproc(iqx)) = 0.d0 ! X-component of the radiative flux vector
  propce(iel,ipproc(iqy)) = 0.d0 ! Y-compnent of the radiative flux vector
  propce(iel,ipproc(iqz)) = 0.d0 ! Z-Component of the radiative flux vector

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

  mode = 1

  if (ippmod(iphpar).le.1) then

    call usray4 &
    ( nvar   , nscal  ,                                            &
      mode   ,                                                     &
      itypfb ,                                                     &
      dt     ,                                                     &
      tparo  , flurdb , tempk(1,1)  )

  else

    call ppray4 &
  ( mode   ,                                                       &
    itypfb ,                                                       &
    propce ,                                                       &
    tparo  , flurdb , tempk(1,1)  )

  endif

  if (ippmod(iccoal).ge.0) then

    ! Particules' temperature
    do icla = 1, nclacp
      ipcla = 1+icla
      do iel = 1, ncel
        tempk(iel,ipcla) = propce(iel,ipproc(itemp2(icla)))
      enddo
    enddo

  ! Fuel
  else if (ippmod(icfuel).ge.0) then

    do icla = 1, nclafu
      ipcla = 1+icla
      do iel = 1, ncel
        tempk(iel,ipcla) = propce(iel,ipproc(itemp2(icla)))
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
  propce(iel,ipproc(icak(1))) = -grand
enddo

!--> Absorption coefficient for different modules

! Warning: for the approximation P-1, the absorption coefficient is required
!          for boundary conditions

if (ippmod(iphpar).ge.2) then

  call ppcabs(propce, tempk, kgi, agi, agbi)

else

  !---> Reading of User datas

  !   - Interface Code_Saturne
  !     ======================

  if (iihmpr.eq.1) then

    call uiray3(propce(1,ipproc(icak(1))), ncel, imodak) !FIXME for ADF

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
        if (propce(iel,ipproc(icak(1))).lt.xkmin) then
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

  if (imoadf.eq.0) then
    call usray3 &
  ( nvar   , nscal  , iappel ,                                     &
    itypfb ,                                                       &
    izfrad ,                                                       &
    dt     ,                                                       &
    propce(1,ipproc(icak(1))))
  endif

endif

!--> General checking

!--> Test if the radiation coeffcient has been assigned
if (iirayo.ge.1) then
  if (imoadf.eq.0) then
    ckmin = propce(1,ipproc(icak(1)))
    do iel = 1, ncel
      ckmin = min(ckmin,propce(iel,ipproc(icak(1))))
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

!=============================================================================
! 4. Solving the ETR
!=============================================================================
! Loop over all grey gases. Remember: In case of the basic radiation models of
! Code_Saturne nwsgg=1

do ngg = 1, nwsgg
  !---> Check of a transparent case
  idverl = idiver

  if (imoadf.ge.1) then
    do iel = 1, ncel
      propce(iel, ipproc(icak(1))) = kgi(iel, ngg) ! TODO merge the two arrays
    enddo

  else
    aa = 0.d0
    do iel = 1, ncel
      aa = max(aa, propce(iel,ipproc(icak(1))))
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
      smbrs(iel) = 3.d0*propce(iel,ipproc(icak(1)))*(tempk(iel,1)**4)          &
                 * agi(iel, ngg)*volume(iel)
    enddo

    !--> Solid phase:

    ! Coal particles: Explicit source term in the transport eqn. of theta4
    if (ippmod(iccoal).ge.0) then
      do icla = 1, nclacp
        ipcla = 1+icla
        do iel = 1,ncel
          smbrs(iel) = smbrs(iel)                               &
                     + (3.d0*propce(iel,ipproc(ix2(icla)))      &
                       *  propce(iel,ipproc(icak(ipcla)))       &
                       * (tempk(iel,ipcla)**4)*agi(iel,ngg)     &
                       *  volume(iel))
        enddo
      enddo

    ! Fuel droplets: Explicit source term in the transport eqn. of theta4
    else if (ippmod(icfuel).ge.0) then
      do icla = 1, nclafu
        ipcla = 1+icla
        call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
        do iel = 1,ncel
          smbrs(iel) =  smbrs(iel)                              &
                     + (3.d0*cvara_yfol(iel)                    &
                       *  propce(iel,ipproc(icak(ipcla)))       &
                       * (tempk(iel,ipcla)**4)*agi(iel,ngg)     &
                       *  volume(iel) )
        enddo
      enddo
    endif

    !--> Gas phase: Implicit source term in the transport eqn. of theta4

    do iel = 1, ncel
      rovsdt(iel) =  3.d0*propce(iel,ipproc(icak(1)))*volume(iel)
    enddo

    !--> Solid phase:

    ! Coal particles: Implicit source term in the transport eqn. of theta4
    if (ippmod(iccoal).ge.0) then
      do icla = 1, nclacp
        ipcla = 1+icla
        do iel = 1,ncel
          rovsdt(iel) = rovsdt(iel)                                      &
                      + (3.d0*propce(iel,ipproc(ix2(icla)))              &
                        * propce(iel,ipproc(icak(ipcla))) * volume(iel) )
        enddo
      enddo

    ! Fuel droplets: Implicit source term in the transport eqn. of theta4
    else if (ippmod(icfuel).ge.0) then
      do icla = 1, nclafu
        ipcla = 1+icla
        call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
        do iel = 1,ncel
          rovsdt(iel) = rovsdt(iel)                                      &
                      + (3.d0*cvara_yfol(iel)                            &
                        * propce(iel,ipproc(icak(ipcla))) * volume(iel) )
        enddo
      enddo

    endif

    ! Radiation coeffcient of the bulk phase
    ! Gas phase:
    do iel = 1, ncel
      ckmel(iel) = propce(iel,ipproc(icak(1)))
    enddo

    if (ippmod(iccoal).ge.0) then
      ! Solid phase:
      ! Coal particles
      do icla = 1, nclacp
        ipcla = 1+icla
        do iel = 1, ncel
          ckmel(iel) = ckmel(iel)                                  &
                     + ( propce(iel,ipproc(ix2(icla)))             &
                       * propce(iel,ipproc(icak(ipcla))) )
        enddo
      enddo
    ! Fuel droplets
    else if (ippmod(icfuel).ge.0) then
      do icla = 1, nclafu
        ipcla = 1+icla
        call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
        do iel = 1, ncel
          ckmel(iel) = ckmel(iel)                                  &
                     + ( cvara_yfol(iel)                           &
                       * propce(iel,ipproc(icak(ipcla))) )
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

    ! Update Boundary condiction coefficients

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
      propce(1,ipproc(iabso(1))),propce(1,ipproc(iemi(1))),          &
      propce(1,ipproc(itsre(1))) ,                                   &
      iqxpar , iqypar , iqzpar ,                                     &
      bqinci , beps   , tparo  ,                                     &
      ckmel  , agbi   , ngg    )

  !===============================================================================
  ! 4.2 Solving of the radiative transfert equation (DOM)
  !===============================================================================

  else if (iirayo.eq.1) then

    !--> Gas phase: Explicit source term of the ETR
    do iel = 1, ncel
      smbrs(iel) = stephn*propce(iel,ipproc(icak(1)))              &
                 *(tempk(iel,1)**4)*agi(iel,ngg)*volume(iel)*unspi
    enddo

    !--> Solid phase:
    ! Coal particles: Explicit source term of the ETR
    if (ippmod(iccoal).ge.0) then
      do icla = 1, nclacp
        ipcla = 1+icla
        do iel = 1, ncel
          smbrs(iel) = smbrs(iel)                                 &
                     + propce(iel,ipproc(ix2(icla)))              &
                       * agi(iel,ngg)*stephn                      &
                       *propce(iel,ipproc(icak(ipcla)))           &
                       *(tempk(iel,ipcla)**4)                     &
                       * volume(iel)*unspi
        enddo
      enddo
    ! Fuel droplets: Explicit source term of the ETR
    elseif (ippmod(icfuel).ge.0) then
      do icla = 1,nclafu
        ipcla = 1+icla
        call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
        do iel = 1,ncel
          smbrs(iel) = smbrs(iel)                         &
                     + cvara_yfol(iel)                    &
                       *agi(iel,ngg)*stephn               &
                       *propce(iel,ipproc(icak(ipcla)))   &
                       *(tempk(iel,ipcla)**4)             &
                       *volume(iel)* unspi
        enddo
      enddo
    endif

    !--> Gas phase: Implicit source term of the ETR
    do iel = 1, ncel
      rovsdt(iel) = propce(iel,ipproc(icak(1))) * volume(iel)
    enddo

    !--> Solid phase
    ! Coal particles: Implicit source term of the ETR
    if (ippmod(iccoal).ge.0) then
      do icla = 1, nclacp
        ipcla = 1+icla
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel)                                      &
                      + propce(iel,ipproc(ix2(icla)))                    &
                        * propce(iel,ipproc(icak(ipcla))) * volume(iel)
        enddo
      enddo
    ! Fuel droplets: Implicit source term of the ETR
    elseif (ippmod(icfuel).ge.0) then
      do icla = 1, nclafu
        ipcla = 1+icla
        call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel)                                      &
                      + cvara_yfol(iel)                                  &
                        * propce(iel,ipproc(icak(ipcla))) * volume(iel)
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
      propce(1,ipproc(itsre(1))),                                    &
      iqxpar , iqypar , iqzpar  ,                                    &
      bqinci , bfnet  , ngg)

  endif

  ! Summing up the quantities of each grey gas
  do iel = 1, ncel

    ! Absorption
    iabgaz(iel) = iabgaz(iel)+(propce(iel,ipproc(icak(1))) *          &
                               propce(iel,ipproc(itsre(1))))
  enddo

  if (ippmod(iccoal).ge.0) then
    do icla = 1, nclacp
      ipcla = 1+icla
      do iel = 1, ncel
        iabpar(iel) = iabpar(iel) + (propce(iel,ipproc(ix2(icla)))  &
                                  * propce(iel,ipproc(icak(ipcla))) &
                                  * propce(iel,ipproc(itsre(1))))
        iabparh2(iel,icla)        = iabparh2(iel,icla)              &
                                  + propce(iel,ipproc(icak(ipcla))) &
                                  * propce(iel,ipproc(itsre(1)))
      enddo
    enddo
  else if (ippmod(icfuel).ge.0) then
    do icla = 1, nclafu
      ipcla = 1+icla
      call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
      do iel = 1, ncel
        iabpar(iel) = iabpar(iel) + (cvara_yfol(iel)                &
                                  * propce(iel,ipproc(icak(ipcla))) &
                                  * propce(iel,ipproc(itsre(1))))
        iabparh2(iel,icla)        = iabparh2(iel,icla)              &
                                  + propce(iel,ipproc(icak(ipcla))) &
                                  * propce(iel,ipproc(itsre(1)))
      enddo
    enddo
  endif

  ! Emission
  do iel = 1, ncel
    iemgex(iel)=iemgex(iel)-(propce(iel,ipproc(icak(1)))            &
                           *agi(iel,ngg)*4.d0*stephn                &
                           *(tempk(iel,1)**4))
    iemgim(iel)=iemgim(iel)-(16.d0*dcp(iel)*propce(iel,ipproc(icak(1)))      &
                           * agi(iel,ngg)*stephn*(tempk(iel,1)**3))

  enddo
  if (ippmod(iccoal).ge.0) then
    do icla = 1, nclacp
      ipcla = 1+icla
      do iel = 1, ncel
        iempex(iel) = iempex(iel) -(4.0d0*propce(iel,ipproc(ix2(icla)))      &
                                  *stephn * propce(iel,ipproc(icak(ipcla)))  &
                                  *(tempk(iel,ipcla)**4)*agi(iel,ngg))
        iempexh2(iel,icla) = iempexh2(iel,icla)                              &
                             -(4.0d0*stephn * propce(iel,ipproc(icak(ipcla)))&
                             *(tempk(iel,ipcla)**4)*agi(iel,ngg))
        iempim(iel) = iempim(iel) -16.d0*propce(iel,ipproc(icak(ipcla)))     &
                                  *propce(iel,ipproc(ix2(icla)))             &
                                  *stephn*(tempk(iel,ipcla)**3)*agi(iel,ngg) &
                                  /cp2ch(ichcor(icla))
        iempimh2(iel,icla) = iempimh2(iel,icla)                              &
                             -16.d0*propce(iel,ipproc(icak(ipcla)))          &
                             *stephn*(tempk(iel,ipcla)**3)*agi(iel,ngg)      &
                             /cp2ch(ichcor(icla))

      enddo
    enddo
  else if (ippmod(icfuel).ge.0) then
    do icla = 1, nclafu
      ipcla = 1+icla
      call field_get_val_prev_s(ivarfl(isca(iyfol(icla))), cvara_yfol)
      do iel = 1, ncel
        iempex(iel) = iempex(iel)   -(4.0d0*cvara_yfol(iel)                   &
                                    *stephn * propce(iel,ipproc(icak(ipcla))) &
                                    *(tempk(iel,ipcla)**4)*agi(iel,ngg))
        iempexh2(iel,icla) = iempexh2(iel,icla)                               &
                             -(4.0d0*stephn * propce(iel,ipproc(icak(ipcla))) &
                             *(tempk(iel,ipcla)**4)*agi(iel,ngg))
        iempim(iel) = iempim(iel) -16.d0*propce(iel,ipproc(icak(ipcla)))      &
                                  *cvara_yfol(iel)*stephn                     &
                                  *(tempk(iel,ipcla)**3) *agi(iel,ngg)/cp2fol
        iempimh2(iel,icla) = iempimh2(iel,icla)                               &
                             -16.d0*propce(iel,ipproc(icak(ipcla)))           &
                             *stephn*(tempk(iel,ipcla)**3)*agi(iel,ngg)/cp2fol
      enddo
    enddo
  endif

  do iel = 1, ncel
    ! Emitted intensity
    ilutot(iel)=ilutot(iel)+propce(iel,ipproc(itsre(1)))

    ! Flux vector components
    propce(iel,ipproc(iqx))=propce(iel,ipproc(iqx))+iqxpar(iel)
    propce(iel,ipproc(iqy))=propce(iel,ipproc(iqy))+iqypar(iel)
    propce(iel,ipproc(iqz))=propce(iel,ipproc(iqz))+iqzpar(iel)
  enddo

  ! If the ADF model is activated we have to sum up the spectral flux densities
  if (imoadf.ge.1) then
    do ifac =1, nfabor
      iqpato(ifac)=iqpato(ifac)+bqinsp(ngg,ifac)
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
  propce(iel,ipproc(ilumin))  = ilutot(iel)
enddo

!===============================================================================
! 6. Net radiative flux at walls: computation and integration
!===============================================================================

!--> Initialization to a non-admissible value for testing after usray5
do ifac = 1,nfabor
  bfnet(ifac) = -grand
enddo

!---> Reading of User datas
!CAREFUL: The user has acces to the radiation coeffcient propce(1,ipproc(icak(1)))
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
  propce(1,ipproc(icak(1)))  )

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
    ! Absoprtion of the gas is copied into iabso(1)
    propce(iel,ipproc(iabso(1))) = iabgaz(iel)
    ! Emission of the gas phase is copied into iemi(1)
    propce(iel,ipproc(iemi(1))) = iemgex(iel)
  enddo
  if (ippmod(iccoal).ge.0.or.ippmod(icfuel).ge.0) then
    do iel = 1, ncel
      ! Absoprtion of particles is added to iabso(1)
      propce(iel,ipproc(iabso(1))) = propce(iel,ipproc(iabso(1))) + iabpar(iel)
      ! Emission of particles is added to iemi(1)
      propce(iel,ipproc(iemi(1))) = propce(iel,ipproc(iemi(1))) + iempex(iel)
    enddo
  endif
  do iel = 1, ncel
    ! Emission + Absorption of gas and particles --> TSexplicit
    propce(iel,ipproc(itsre(1))) = propce(iel,ipproc(iabso(1)))                &
                                 + propce(iel,ipproc(iemi(1)))
  enddo

  do iel = 1, ncel
    ! TSimplicit of the gas phase
    propce(iel,ipproc(itsri(1))) = iemgim(iel)
  enddo
  if (ippmod(iccoal).ge.0.or.ippmod(icfuel).ge.0) then
    do iel = 1, ncel
      ! TSimplicit of the solid phase is added to istri(1)
      propce(iel,ipproc(itsri(1))) = propce(iel,ipproc(itsri(1))) + iempim(iel)
    enddo
  endif

  ! In order to determine the source terms of the particle enthalpy tranport eqn.,
  ! we have to copie the approriate determined aboce into the corressponding tables
  if (ippmod(iccoal).ge.0.or.ippmod(icfuel).ge.0) then
    do icla = 1,nclacp
      ipcla = 1+icla
        do iel = 1, ncel
          propce(iel,ipproc(iabso(ipcla))) = iabparh2(iel,icla)
          propce(iel,ipproc(iemi(ipcla))) = iempexh2(iel,icla)
          propce(iel,ipproc(itsre(ipcla)))= iabparh2(iel,icla)+iempexh2(iel,icla)
          propce(iel,ipproc(itsri(ipcla)))= iempimh2(iel,icla)
        enddo
    enddo
  endif
else
  do iel = 1, ncel
    propce(iel,ipproc(iabso(1)))  = 0.d0
    propce(iel,ipproc(iemi(1)))  = 0.d0
    propce(iel,ipproc(itsre(1))) = 0.d0
    propce(iel,ipproc(itsri(1))) = 0.d0
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
    q(1,iel) = propce(iel,ipproc(iqx))
    q(2,iel) = propce(iel,ipproc(iqy))
    q(3,iel) = propce(iel,ipproc(iqz))
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
    propce(iel,ipproc(itsre(1))) = - grad(1,1,iel)                &
                                   - grad(2,2,iel)                &
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
! 7.3 Explicite radiative semi-analytical corrected source term
!===============================================================================

if (idverl.eq.2) then

  !---> Comparison of the semi-analytical and conservative source terms
  aa = 0.d0
  do iel = 1, ncel
    aa = aa + propce(iel,ipproc(itsre(1))) * volume(iel)
  enddo

  bb = 0.d0
  do iel = 1,ncel
    bb = bb                                                                    &
       + (propce(iel,ipproc(iabso(1)))+propce(iel,ipproc(iemi(1))))*volume(iel)
  enddo

  if(irangp.ge.0) then
    call parsom(aa)
    call parsom(bb)
  endif

  aa = aa/bb

  !---> Correction of the semi-analytical source term by the conservative source
  ! term
  do iel = 1,ncel
    propce(iel,ipproc(itsre(1))) = ( propce(iel,ipproc(iabso(1)))              &
                                   + propce(iel,ipproc(iemi(1))))             &
                                 * aa
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
    aa = aa + propce(iel,ipproc(itsre(1))) * volume(iel)
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
deallocate(ckmel)
deallocate(tempk)
deallocate(coefap, coefbp)
deallocate(cofafp, cofbfp)
deallocate(flurds, flurdb)
deallocate(kgi,agi,agbi)
deallocate(iqxpar,iqypar,iqzpar)
deallocate(iabgaz,iabpar,iemgex,iempex,ilutot)
deallocate(iemgim,iempim)
deallocate(iabparh2,iempexh2,iempimh2)

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

