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
!> \file attssc.f90
!>
!> \brief Additional right-hand side source terms for scalar equations
!> taking into account dry and humid atmospheric variables.
!> If 1D atmospheric radiative module is used (iatra1 = 1) additional source
!> terms for the thermal scalar equation to take into account the radiative
!> forcing.
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   iscal           scalar number
!> \param[in]   crvexp          explicit part of the right hand side
!_______________________________________________________________________________

subroutine attssc ( iscal, crvexp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use atincl
use field
use radiat
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision crvexp(ncelet)

! Local variables

character(len=80) :: chaine
integer          ivar,  iel

double precision pp, dum

double precision, dimension(:), allocatable :: ray3Di, ray3Dst
double precision, dimension(:,:), allocatable, save :: grad1, grad2
double precision, dimension(:), allocatable, save :: r3

double precision, save :: qliqmax,r3max
logical, save :: r3_is_defined = .false.
integer, save :: treated_scalars = 0

double precision, dimension(:), allocatable :: pphy
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cvar_ntdrp
double precision, dimension(:), pointer :: cvar_pottemp
double precision, dimension(:), pointer :: cpro_tempc
double precision, dimension(:), pointer :: cpro_liqwt
double precision, dimension(:), pointer :: cpro_rad_cool
double precision, dimension(:), pointer :: cpro_met_p

!===============================================================================
! 1. Initialization
!===============================================================================

! variable number from scalar number
ivar = isca(iscal)

! variable name
call field_get_name(ivarfl(ivar), chaine)

! density
call field_get_val_s(icrom, crom)

if (imeteo.ge.2) then
  call field_get_val_s_by_name('meteo_pressure', cpro_met_p)
endif

!===============================================================================
! 2. Taking into account radiative forcing for the 1d radiative module
!    (if the 3D module is not activated)
!===============================================================================

if (ippmod(iatmos).ge.1.and.iatra1.ge.1.and.iirayo.eq.0) then

  call field_get_val_s(ivarfl(isca(iscalt)), cvar_pottemp)
  call field_get_val_s(itempc, cpro_tempc)

  !   2.1 Source terms in the equation of the liquid potential temperature

  if (ivar.eq.isca(iscalt)) then

    allocate(ray3Di(ncel))
    allocate(ray3Dst(ncel))

    ! Call the 1D radiative model
    ! Compute the divergence of the ir and solar radiative fluxes:
    call atr1vf

    ! Cressman interpolation of the 1D radiative fluxes on the 3D mesh:
    ! Infra red
    call mscrss(idrayi, 1, ray3Di)

    ! Sun
    call mscrss(idrayst, 1, ray3Dst)

    ! Store radiative fluxes for droplet nucleation model
    ! for humid atmosphere
    if (ippmod(iatmos).eq.2.and.modsedi.eq.1.and.modnuc.gt.0) then
      call field_get_val_s_by_name('radiative_cooling', cpro_rad_cool)
      do iel = 1, ncel
        cpro_rad_cool(iel) = (ray3Dst(iel)-ray3Di(iel) )
      enddo
    endif

    ! Explicit source term for the thermal scalar equation:

    do iel = 1, ncel
      crvexp(iel) = crvexp(iel) +                                   &
        cp0*cell_f_vol(iel)*crom(iel)*(-ray3Di(iel) + ray3Dst(iel)) &
        ! Conversion Temperature -> Potential Temperature
        * cvar_pottemp(iel) / (cpro_tempc(iel) + tkelvi)

    enddo

    deallocate(ray3Di)
    deallocate(ray3Dst)

  endif

endif

!===============================================================================
! 3. Take into source terms fort thetal, qw and nc due to sedimentation of drops
!===============================================================================
! FIXME gravity is assumed to follow z-axis direction

if (ippmod(iatmos).eq.2.and.modsedi.eq.1) then ! for humid atmo. physics only

  call field_get_val_s(iliqwt, cpro_liqwt)
  call field_get_val_s(itempc, cpro_tempc)

  ! Test minimum liquid water to carry out drop sedimentation
  qliqmax = 0.d0
  do iel = 1, ncel
    qliqmax = max(cpro_liqwt(iel),qliqmax)
  enddo
  if (irangp.ge.0) call parmax(qliqmax)

  if (qliqmax.gt.1d-8) then

    if (.not.r3_is_defined)then

      call field_get_val_s(ivarfl(isca(intdrp)), cvar_ntdrp)

      allocate(r3(ncelet))
      ! Compute the mean value: (<r^3>)**1/3
      call define_r3
      r3_is_defined = .true.

      allocate(grad1(3,ncelet), grad2(3,ncelet))

      call grad_sed(grad1, grad2)

    endif ! r3_not_defined

    ivar = isca(iscal)
    if (ivar.eq.isca(iscalt)) then

      do iel = 1, ncel
        if (imeteo.eq.0) then
          call atmstd(xyzcen(3,iel),pp,dum,dum)
        else if (imeteo.eq.1) then
          call intprf &
               ( nbmett, nbmetm,                                        &
                 ztmet , tmmet , phmet , xyzcen(3,iel) , ttcabs, pp )
        else
          pp = cpro_met_p(iel)
        endif

        crvexp(iel) = crvexp(iel) -clatev*(ps/pp)**(rair/cp0)           &
                    *cell_f_vol(iel)*grad1(3,iel)
      enddo
      treated_scalars = treated_scalars + 1

    elseif (ivar.eq.isca(iymw)) then

      do iel = 1, ncel
        crvexp(iel) = crvexp(iel) - cell_f_vol(iel)*grad1(3,iel)
      enddo

      treated_scalars = treated_scalars + 1

    elseif (ivar.eq.isca(intdrp)) then

      do iel = 1, ncel
        crvexp(iel) = crvexp(iel) + cell_f_vol(iel)*grad2(3,iel)
      enddo

      treated_scalars = treated_scalars + 1

    endif

    treated_scalars = mod(treated_scalars, 3)

    if (treated_scalars.eq.0) then ! keeping same gradients for 3 atm. var.
      deallocate(r3)
      r3_is_defined = .false.
      deallocate(grad1)
      deallocate(grad2)
    endif
  endif ! qliqmax.gt.1.d-8
endif ! for humid atmosphere physics only

!--------
! Formats
!--------

return

!===============================================================================

contains

  !-----------------------------------------------------------------------------

  !> \brief Compute the mean volumic radius of the droplets

  subroutine define_r3

    !===========================================================================

    !===========================================================================
    ! Module files
    !===========================================================================

    use cstnum, only: pi

    !===========================================================================

    implicit none


    ! Local variables

    double precision rho, qliq, nc

    double precision rho_water
    parameter (rho_water=1.d+3) ! FIXME should be defined somewhere else
    double precision conversion
    parameter (conversion=1d+6)! passing from 1/cm**3 to 1/m**3
    !===========================================================================

    r3max = 0.d0
    do iel = 1, ncel
      rho = crom(iel)
      qliq = cpro_liqwt(iel)
      nc = cvar_ntdrp(iel)
      if (qliq.ge.1d-8) then!FIXME useless
        nc = max(nc,1.d0)
        r3(iel) = (0.75d0/pi*(rho*qliq)/(rho_water*nc*conversion))**(1.d0/3.d0)
      else
        r3(iel) = 0.d0
      endif
      r3max = max(r3(iel),r3max)
    enddo

    if (irangp.ge.0) call parmax (r3max)
  end subroutine define_r3

  !-----------------------------------------------------------------------------

  !> \brief Computation of the gradient of the two following quantities
  !> 1) rho*qliq*V(r3)*exp(5*sc^2)
  !> 2) nc*V(r3)*exp(-sc^2)
  !> FIXME describe quantities

  !> \param[out]    grad1
  !> \param[out]    grad2

  subroutine grad_sed(grad1, grad2)

    !===========================================================================

    !===========================================================================
    ! Module files
    !===========================================================================

    use cs_c_bindings
    use pointe, only: itypfb

    !===========================================================================

    implicit none

    ! Arguments

    double precision grad1(3,ncelet), grad2(3,ncelet)

    ! Local variables

    double precision climgp, epsrgp, depo

    integer    iccocg, ii, iifld, imligp, inc, iwarnp, imrgrp, nswrgp, ifac, iel
    double precision, dimension(:), allocatable :: local_coefa, local_coefb
    double precision, dimension(:), allocatable :: local_field, sed_vel
    double precision, dimension(:), allocatable :: pres

    double precision, dimension(:), pointer :: rugd
    double precision, dimension(:), pointer :: bcfnns, ustar, rugt
    double precision, dimension(:), pointer :: cpro_tempc

    type(var_cal_opt) :: vcopt

    !===========================================================================

    ! no droplets case

    if (r3max.lt.1.d-10) then
      do iel = 1, ncel
        do ii = 1, 3
          grad1(ii,iel) = 0.d0
          grad2(ii,iel) = 0.d0
        enddo
      enddo

      return
    endif

    ! compute sedimentation velocity

    allocate(sed_vel(ncel))
    do iel = 1, ncel
      sed_vel(iel) = 1.19d+8 * (r3(iel))**2  ! taup g, with taup = cuning * d^2 * rhop / (18 * mu) ...
    enddo

    ! take into account deposition if enabled

    if (moddep.gt.0) then
      allocate(pres(ncel))
      call field_get_val_s(itempc, cpro_tempc)

      do iel = 1, ncel
        if (imeteo.eq.0) then
          call atmstd(xyzcen(3,iel),pres(iel),dum,dum)
        else if (imeteo.eq.1) then
          call intprf                                                      &
              ( nbmett, nbmetm,                                            &
                ztmet , tmmet , phmet , xyzcen(3,iel) , ttcabs, pres(iel) )
        else
          pres(iel) = cpro_met_p(iel)
        endif
      enddo

      call field_get_val_s_by_name('non_neutral_scalar_correction', bcfnns)
      call field_get_val_s_by_name('ustar', ustar)
      call field_get_val_s_by_name('boundary_roughness', rugd)
      call field_get_val_s_by_name('boundary_thermal_roughness', rugt)

      do ifac = 1, nfabor
        if (itypfb(ifac).eq.iparug.or.itypfb(ifac).eq.iparoi) then
          iel  = ifabor(ifac)
          if (r3(iel).gt.0.d0) then
            if (rugd(ifac).gt.0.d0) then
              call deposition_vel(cpro_tempc(iel), crom(iel), pres(iel),         &
                                  bcfnns(ifac), ustar(ifac), rugt(ifac),   &
                                  r3(iel), sed_vel(iel), depo)

              sed_vel(iel) = sed_vel(iel) + depo
            endif
          endif
        endif ! itypfb.eq.iparug
      enddo

      deallocate(pres)
    endif ! moddep.gt.0

    ! options for gradient calculation
    iccocg = 1
    inc = 1
    iifld = -1
    call field_get_key_struct_var_cal_opt(ivarfl(isca(iymw)), vcopt)
    imrgrp = vcopt%imrgra
    nswrgp = vcopt%nswrgr
    epsrgp = vcopt%epsrgr
    imligp = vcopt%imligr
    iwarnp = vcopt%iwarni
    climgp = vcopt%climgr

    ! homogeneous Neumann BCs for gradient computation
    allocate(local_coefa(nfabor))
    do ifac = 1, nfabor
      local_coefa(ifac) = 0.d0
    enddo
    allocate(local_coefb(nfabor))
    do ifac = 1, nfabor
      local_coefb(ifac) = 1.d0
    enddo

    allocate(local_field(ncelet))

    ! Computation of the gradient of rho*qliq*V(r3)*exp(5*sc^2)
    ! it corresponds to div(qliq* exp() rho vel_d)
    do iel = 1, ncel
      local_field(iel) = crom(iel)        & ! mass density of the air kg/m3
                        *cpro_liqwt(iel)  & ! total liquid water content kg/kg
                        *sed_vel(iel)     & ! deposition velocity m/s
                        *exp(5.d0*sigc**2)  ! coefficient coming from log-norm
                                            ! law of the droplet spectrum
    enddo

    call gradient_s                                                 &
   ( iifld  , imrgrp , inc    , iccocg , nswrgp ,imligp,            &
     iwarnp , epsrgp , climgp ,                                     &
     local_field     , local_coefa , local_coefb ,                  &
     grad1   )

    ! Computation of the gradient of Nc*V(r3)*exp(-sc^2)

    do iel = 1, ncel
      local_field(iel) = crom(iel)        & ! mass density of the air kg/m3
                        *cvar_ntdrp(iel)  & ! number of droplets 1/cm**3
                        *sed_vel(iel)     & ! deposition velocity m/s
                        *exp(-sigc**2)      ! coefficient coming from log-normal
                                            ! law of the droplet spectrum
    enddo

    call gradient_s                                                 &
   ( iifld  , imrgrp , inc    , iccocg , nswrgp ,imligp,            &
     iwarnp , epsrgp , climgp ,                                     &
     local_field     , local_coefa , local_coefb ,                  &
     grad2   )

    deallocate(sed_vel)

    deallocate(local_coefa)
    deallocate(local_coefb)
    deallocate(local_field)

  end subroutine grad_sed

  !-----------------------------------------------------------------------------

  !> \brief Compute deposition velocity

  !> \param[in]       tempc       in Celsius
  !> \param[in]       rom
  !> \param[in]       pres
  !> \param[in]       cfnns       non neutral correction coefficient for scalars
  !> \param[in]       ustar
  !> \param[in]       rugt
  !> \param[in]       rcloudvolmoy
  !> \param[in]       wg
  !> \param[in]       depo

  subroutine deposition_vel(tempc, rom , pres,          &
                            cfnns, ustar, rugt,        &
                            rcloudvolmoy, wg, depo)

    !===========================================================================

    use cstnum

    implicit none

    !===========================================================================

    ! Arguments

    double precision tempc, rom, pres
    double precision cfnns, ustar, rugt
    double precision rcloudvolmoy, depo

    ! Local variables

    double precision temp
    double precision ckarm, eps0, cbolz, gamma, alpha, arecep
    double precision dp
    double precision muair, nuair, dbrow, cebro, lpm, ccunning
    double precision wg, ather
    double precision raero, st, ceimp, ceint, rsurf,rhoeau
    double precision dzmin

    !===========================================================================

    ! deposition is computed only for first level
    ckarm  = 0.4d0
    eps0   = 3.d0
    cbolz  = 1.38d-23
    gamma  = 0.56d0
    alpha  = 1.5d0
    arecep = 0.01d0
    rhoeau = 1000.d0
    dzmin  = 4.d0

    dp = dzmin/2.d0
    temp = tempc + tkelvi

    muair = 1.83d-5*(416.16d0/(temp+120.d0))                    &
           *((temp/296.16d0)**1.5d0)

    nuair = muair/rom

    lpm = (2.d0*muair/pres)*((0.125d0*pi*rair*temp)**(0.5))
    ccunning = 1.d0 + (lpm/rcloudvolmoy)*(1.257d0 + 0.4d0       &
              *exp(-1.1d0*rcloudvolmoy/lpm))

    dbrow = cbolz*temp*ccunning/                                &
            (6.d0*pi*muair*rcloudvolmoy)

    cebro = nuair**((-1.d0)*gamma)/dbrow

    ather = ckarm/log((dp+rugt)/rugt)
    raero = 1.d0 / (ather * ustar * cfnns)

    st = wg*ustar/(9.81d0*arecep)
    ceimp = (st/(st+alpha))**(2.)
    ceint = 2.d0*((rcloudvolmoy/arecep)**(2.))

    ! Fog or cloud droplet deposition
    if (ustar.gt.0.d0) then
      rsurf = 1.d0 / (eps0*ustar*(ceimp+ceint+cebro)*exp(-sqrt(st)))
      depo = 1.d0 / (raero + rsurf)
    else
      depo = 0.d0
    endif

  end subroutine deposition_vel

  !-----------------------------------------------------------------------------

end subroutine attssc
