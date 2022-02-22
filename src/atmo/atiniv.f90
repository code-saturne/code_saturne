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
!> \file atiniv.f90
!> \brief Initialisation of calculation variables for the atmospheric module,
!> it is the counterpart of usiniv.f90.
!>
!> Initialise for example the meteorological field for each cell of
!> the domain by interpolation of the data from the meteo file
!> First stage, before GUI.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!-------------------------------------------------------------------------------

subroutine atiniv0

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use ppincl
use atincl
use mesh
use atchem
use sshaerosol
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer          imode, iel
integer          k,ii, isc
integer          fid_axz
double precision d2s3
double precision zent,xuent,xvent, xwent, xkent,xeent,tpent,qvent,ncent
double precision xcent
double precision r_nt
double precision vel_dir(3), shear_dir(3)

type(var_cal_opt) :: vcopt_p, vcopt_u

double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi
double precision, dimension(:), pointer :: cvar_fb, cvar_omg, cvar_nusa
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:), pointer :: cvar_despgi, cvar_sc
double precision, dimension(:), pointer :: cvar_scalt, cvar_totwt, cvar_ntdrp
double precision, dimension(:,:), pointer :: cpro_met_vel
double precision, dimension(:), pointer :: cpro_met_potemp
double precision, dimension(:), pointer :: cpro_met_qv, cpro_met_nc
double precision, dimension(:), pointer :: cpro_met_k, cpro_met_eps
double precision, dimension(:), pointer :: cpro_met_axz

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

call field_get_id_try('meteo_shear_anisotropy', fid_axz)
if (fid_axz.ne.-1) then
  call field_get_val_s(fid_axz, cpro_met_axz)
endif

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

d2s3 = 2.d0/3.d0

if (itytur.eq.2) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (itytur.eq.3) then
  call field_get_val_v(ivarfl(irij), cvar_rij)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (iturb.eq.50) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
  call field_get_val_s(ivarfl(iphi), cvar_phi)
  call field_get_val_s(ivarfl(ifb), cvar_fb)
elseif (iturb.eq.60) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iomg), cvar_omg)
elseif (iturb.eq.70) then
  call field_get_val_s(ivarfl(inusa), cvar_nusa)
endif

if (imeteo.ge.2) then
  call field_get_val_s_by_name('meteo_pot_temperature', cpro_met_potemp)
  call field_get_val_v_by_name('meteo_velocity', cpro_met_vel)
  call field_get_val_s_by_name('meteo_tke', cpro_met_k)
  call field_get_val_s_by_name('meteo_eps', cpro_met_eps)
  if (ippmod(iatmos).eq.2) then
    call field_get_val_s_by_name('meteo_humidity', cpro_met_qv)
    call field_get_val_s_by_name('meteo_drop_nb', cpro_met_nc)
  endif
endif

!===============================================================================
! 2. READING THE METEO PROFILE FILE (IF IMETEO = 1 DEFAULT OPTION):
!===============================================================================

! Meteo file
if (imeteo.eq.1) then
  imode = 1
  call atlecm(imode)

! Recomputed from cs_glog_atmo_option values
else if (imeteo.eq.2) then
  call cs_atmo_compute_meteo_profiles()
endif

if (iatra1.gt.0) then

  imode = 1
  call usatdv(imode)

endif

! Atmospheric gaseous chemistry
if (ifilechemistry.ge.1) then

  ! Second reading of chemical profiles file
  imode = 1
  call atlecc(imode)

  ! Volume initilization with profiles for species present
  ! in the chemical profiles file
  if (isuite.eq.0 .or. (isuite.ne.0.and.init_at_chem.eq.1)) then

    do k = 1, nespgi
      call field_get_val_s(ivarfl(isca(isca_chem(idespgi(k)))), cvar_despgi)

      do iel = 1, ncel
        zent = xyzcen(3,iel)
        call intprf                                                         &
        (nbchmz, nbchim,                                                    &
         zproc, tchem, espnum(1+(k-1)*nbchim*nbchmz), zent  , ttcabs, xcent )
        ! The first nespg user scalars are supposed to be chemical species
        cvar_despgi(iel) = xcent
      enddo

    enddo
  endif

endif

! Atmospheric gaseous chemistry
if (ichemistry.ge.1) then

  ! Computation of the conversion factor matrix used for
  ! the reaction rates jaccobian matrix
  do ii = 1, nespg
    do k = 1, nespg
      conv_factor_jac((chempoint(k)-1)*nespg+chempoint(ii)) = dmmk(ii)/dmmk(k)
    enddo
  enddo

endif

! Atmospheric aerosol chemistry
if (iaerosol.ne.CS_ATMO_AEROSOL_OFF) then

  ! Reading intial concentrations and numbers from file
  !   or from the aerosol library
  call atleca()

  ! Initialization
  if (isuite.eq.0 .or. init_at_chem.eq.1) then

    ! Writing
    if (vcopt_u%iwarni.ge.1.or.vcopt_p%iwarni.ge.1) then
      write(nfecra,2001)
    endif

    do ii = 1, nlayer_aer*n_aer + n_aer
      isc = isca_chem(nespg + ii)
      call field_get_val_s(ivarfl(isca(isc)), cvar_sc)
      do iel = 1, ncel
        cvar_sc(iel) = dlconc0(ii)
      enddo
    enddo
  endif

  ! Do not free memory allocated in atleca, dlconc0 is used in attycl

endif

! Check simulation times used by atmo
! radiative transfer or chemistry models
if (     (iatra1.eq.1.or.ifilechemistry.ge.1)              &
    .and.(syear.eq.-1.or.squant.eq.-1.or.shour.eq.-1 &
          .or.smin.eq.-1.or.ssec.le.-1.d0)) then
  if (iatra1.eq.1) write(nfecra,1000)
  if (ifilechemistry.ge.1) write(nfecra,1001)
  call csexit (1)
endif

! Check radiative module latitude / longitude
if (iatra1.eq.1 .and. (xlat.ge.rinfin*0.5 .or. xlon.ge.rinfin*0.5)) then
  write(nfecra,1002)
  call csexit (1)
endif

! Check latitude / longitude from meteo file
if (imeteo.eq.1) then
  if (maxval(xmet).ge.rinfin*0.5 .or. maxval(ymet).ge.rinfin*0.5) then
    write(nfecra,1003)
    call csexit (1)
  endif
endif

! Check latitude / longitude from chemistry file
if (ifilechemistry.ge.1) then
  if (maxval(xchem).ge.rinfin*0.5 .or. maxval(ychem).ge.rinfin*0.5) then
    write(nfecra,1004)
    call csexit (1)
  endif
endif


!===============================================================================
! 3. Dry atmosphere: default initialization of potential temperature
!===============================================================================

! Only if the simulation is not a restart from another one
if (isuite.eq.0) then

  if (initmeteo.eq.1) then

    if (ippmod(iatmos).eq.1) then
      call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
    else if (ippmod(iatmos).eq.2) then
      call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
      call field_get_val_s(ivarfl(isca(iymw)), cvar_totwt)
      call field_get_val_s(ivarfl(isca(intdrp)), cvar_ntdrp)
    endif

    if (imeteo.eq.0) then

      if (ippmod(iatmos).eq.1) then
        ! The thermal scalar is potential temperature
        do iel = 1, ncel
          cvar_scalt(iel) = t0
        enddo
      endif

      if (ippmod(iatmos).eq.2) then
        ! The thermal scalar is liquid potential temperature
        do iel = 1, ncel
          cvar_scalt(iel) = t0
          cvar_totwt(iel) = 0.d0
          cvar_ntdrp(iel) = 0.d0
        enddo
      endif

    ! Only if meteo file is present:
    else

      ! Writing
      if (vcopt_u%iwarni.ge.1.or.vcopt_p%iwarni.ge.1) then
        write(nfecra,2000)
      endif

      do iel = 1, ncel

        zent = xyzcen(3,iel)

        ! Meteo file
        if (imeteo.eq.1) then
          call intprf &
            (nbmetd, nbmetm,                                               &
            zdmet, tmmet, umet , zent  , ttcabs, xuent )

          call intprf &
            (nbmetd, nbmetm,                                               &
            zdmet, tmmet, vmet , zent  , ttcabs, xvent )

          call intprf &
            (nbmetd, nbmetm,                                               &
            zdmet, tmmet, ekmet, zent  , ttcabs, xkent )

          call intprf &
            (nbmetd, nbmetm,                                               &
            zdmet, tmmet, epmet, zent  , ttcabs, xeent )

          xwent = 0.d0
        else
          xuent = cpro_met_vel(1, iel)
          xvent = cpro_met_vel(2, iel)
          xwent = cpro_met_vel(3, iel)
          xkent = cpro_met_k(iel)
          xeent = cpro_met_eps(iel)
        endif

        vel(1,iel) = xuent
        vel(2,iel) = xvent
        vel(3,iel) = xwent

        ! Velocity direction normalized
        vel_dir(1) = xuent
        vel_dir(2) = xvent
        vel_dir(3) = xwent
        call vector_normalize(vel_dir, vel_dir)
        shear_dir(1) = 0.d0
        shear_dir(2) = 0.d0
        if (fid_axz.eq.-1) then
          shear_dir(3) = -sqrt(cmu) ! Rxz/k
        else
          shear_dir(3) = cpro_met_axz(iel) ! Rxz/k
        endif

        ! ITYTUR est un indicateur qui vaut ITURB/10
        if    (itytur.eq.2) then

          cvar_k(iel)  = xkent
          cvar_ep(iel) = xeent

        elseif (itytur.eq.3) then

          r_nt = - sqrt(cmu) * xkent
          cvar_rij(1,iel) = d2s3*xkent
          cvar_rij(2,iel) = d2s3*xkent
          cvar_rij(3,iel) = d2s3*xkent
          ! Rxy
          cvar_rij(4,iel) = xkent * &
             (vel_dir(1)*shear_dir(2)+vel_dir(2)*shear_dir(1))
          ! Ryz
          cvar_rij(5,iel) = xkent * &
             (vel_dir(2)*shear_dir(3)+vel_dir(3)*shear_dir(2))
          ! Rxz
          cvar_rij(6,iel) = xkent * &
             (vel_dir(1)*shear_dir(3)+vel_dir(3)*shear_dir(1))
          cvar_ep(iel)  = xeent

        elseif (iturb.eq.50) then

          cvar_k(iel)   = xkent
          cvar_ep(iel)  = xeent
          cvar_phi(iel) = d2s3
          cvar_fb(iel)  = 0.d0

        elseif (iturb.eq.60) then

          cvar_k(iel)   = xkent
          cvar_omg(iel) = xeent/cmu/xkent

        elseif (iturb.eq.70) then

          cvar_nusa(iel) = cmu*xkent**2/xeent

        endif

        if (ippmod(iatmos).eq.1) then
          if (imeteo.eq.1) then
            ! The thermal scalar is potential temperature
            call intprf &
              (nbmett, nbmetm,                                               &
              ztmet, tmmet, tpmet, zent  , ttcabs, tpent )
          else
            tpent = cpro_met_potemp(iel)
          endif

          cvar_scalt(iel) = tpent
        endif

        if (ippmod(iatmos).eq.2) then

          if (imeteo.eq.1) then
            ! The thermal scalar is liquid potential temperature
            call intprf(nbmett, nbmetm,   &
                        ztmet, tmmet, tpmet, zent, ttcabs, tpent)

            call intprf(nbmett, nbmetm,   &
                        ztmet, tmmet, qvmet, zent, ttcabs, qvent)

            call intprf(nbmett, nbmetm,  &
                        ztmet, tmmet, ncmet, zent, ttcabs, ncent)
          else
            tpent = cpro_met_potemp(iel)
            qvent = cpro_met_qv(iel)
            ncent = cpro_met_nc(iel)
          endif
          cvar_scalt(iel) = tpent
          cvar_totwt(iel) = qvent
          cvar_ntdrp(iel) = ncent
        endif

      enddo

    endif
  endif

endif

!--------
! Formats
!--------

 2000 format(/,                                                   &
'   ** INIT DYNAMIC VARIABLES FROM METEO FILE'                ,/,&
'      --------------------------------------'                ,/)

 2001 format(/,                                                   &
'   ** INIT ATMO CHEMISTRY VARIABLE FROM FILE'                ,/,&
'      --------------------------------------'                ,/)

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                ATMOSPHERIC  MODULE                         ',/,&
'@                RADITIVE MODEL (IATRA1)                     ',/,&
'@                                                            ',/,&
'@    The simulation time is wrong                            ',/,&
'@    Check variables syear, squant, shour, smin, ssec        ',/,&
'@                                                            ',/,&
'@    By decreasing priority, these variables can be defined  ',/,&
'@    in cs_user_parameters or the meteo file                 ',/,&
'@    or the chemistry file                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@    The simulation time is wrong                            ',/,&
'@    Check variables syear, squant, shour, smin, ssec        ',/,&
'@                                                            ',/,&
'@    By decreasing priority, these variables can be defined  ',/,&
'@    in cs_user_parameters or the meteo file                 ',/,&
'@    or the chemistry file                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                ATMOSPHERIC  MODULE                         ',/,&
'@                RADITIVE MODEL (IATRA1)                     ',/,&
'@                                                            ',/,&
'@    Wrong xlat and xlon coordinates.                        ',/,&
'@                                                            ',/,&
'@    See cs_user_parameters.f90                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1003 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC MODULE                                    ',/,&
'@                                                            ',/,&
'@    Wrong coordinates xmet, ymet for the meteo profile.     ',/,&
'@                                                            ',/,&
'@    See cs_user_parameters.f90                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1004 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@    Wrong xchem, ychem coordinates for the concentration    ',/,&
'@    profiles (chemistry model).                             ',/,&
'@                                                            ',/,&
'@    See cs_user_parameters.f90                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return

end subroutine atiniv0

!-------------------------------------------------------------------------------
!> \file atiniv.f90
!> \brief Initialisation of calculation variables for the atmospheric module,
!> it is the counterpart of usiniv.f90.
!>
!> Initialise for example the meteorological field for each cell of
!> the domain by interpolation of the data from the meteo file
!> Second stage.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!-------------------------------------------------------------------------------

subroutine atiniv1

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use ppincl
use atincl
use mesh
use atchem
use sshaerosol
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

double precision d2s3

double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi
double precision, dimension(:), pointer :: cvar_fb, cvar_omg, cvar_nusa
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:,:), pointer :: cpro_met_vel
double precision, dimension(:), pointer :: cpro_met_potemp
double precision, dimension(:), pointer :: cpro_met_qv, cpro_met_nc
double precision, dimension(:), pointer :: cpro_met_k, cpro_met_eps

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

d2s3 = 2.d0/3.d0

if (itytur.eq.2) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (itytur.eq.3) then
  call field_get_val_v(ivarfl(irij), cvar_rij)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (iturb.eq.50) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
  call field_get_val_s(ivarfl(iphi), cvar_phi)
  call field_get_val_s(ivarfl(ifb), cvar_fb)
elseif (iturb.eq.60) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iomg), cvar_omg)
elseif (iturb.eq.70) then
  call field_get_val_s(ivarfl(inusa), cvar_nusa)
endif

if (imeteo.ge.2) then
  call field_get_val_s_by_name('meteo_pot_temperature', cpro_met_potemp)
  call field_get_val_v_by_name('meteo_velocity', cpro_met_vel)
  call field_get_val_s_by_name('meteo_tke', cpro_met_k)
  call field_get_val_s_by_name('meteo_eps', cpro_met_eps)
  if (ippmod(iatmos).eq.2) then
    call field_get_val_s_by_name('meteo_humidity', cpro_met_qv)
    call field_get_val_s_by_name('meteo_drop_nb', cpro_met_nc)
  endif
endif

!===============================================================================
! 2. Checks of the reading the meteo profile file (if imeteo = 1 default option)
!===============================================================================

! Check simulation times used by atmo
! radiative transfer or chemistry models
if (     (iatra1.eq.1.or.ifilechemistry.ge.1)              &
    .and.(syear.eq.-1.or.squant.eq.-1.or.shour.eq.-1 &
          .or.smin.eq.-1.or.ssec.le.-1.d0)) then
  if (iatra1.eq.1) write(nfecra,1000)
  if (ifilechemistry.ge.1) write(nfecra,1001)
  call csexit (1)
endif

! Check radiative module latitude / longitude
if (iatra1.eq.1 .and. (xlat.ge.rinfin*0.5 .or. xlon.ge.rinfin*0.5)) then
  write(nfecra,1002)
  call csexit (1)
endif

! Check latitude / longitude from meteo file
if (imeteo.eq.1) then
  if (maxval(xmet).ge.rinfin*0.5 .or. maxval(ymet).ge.rinfin*0.5) then
    write(nfecra,1003)
    call csexit (1)
  endif
endif

! Check latitude / longitude from chemistry file
if (ifilechemistry.ge.1) then
  if (maxval(xchem).ge.rinfin*0.5 .or. maxval(ychem).ge.rinfin*0.5) then
    write(nfecra,1004)
    call csexit (1)
  endif
endif

!--------
! Formats
!--------

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                ATMOSPHERIC  MODULE                         ',/,&
'@                RADITIVE MODEL (IATRA1)                     ',/,&
'@                                                            ',/,&
'@    The simulation time is wrong                            ',/,&
'@    Check variables syear, squant, shour, smin, ssec        ',/,&
'@                                                            ',/,&
'@    By decreasing priority, these variables can be defined  ',/,&
'@    in cs_user_parameters or the meteo file                 ',/,&
'@    or the chemistry file                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@    The simulation time is wrong                            ',/,&
'@    Check variables syear, squant, shour, smin, ssec        ',/,&
'@                                                            ',/,&
'@    By decreasing priority, these variables can be defined  ',/,&
'@    in cs_user_parameters or the meteo file                 ',/,&
'@    or the chemistry file                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                ATMOSPHERIC  MODULE                         ',/,&
'@                RADITIVE MODEL (IATRA1)                     ',/,&
'@                                                            ',/,&
'@    Wrong xlat and xlon coordinates.                        ',/,&
'@                                                            ',/,&
'@    See cs_user_parameters.f90                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1003 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC MODULE                                    ',/,&
'@                                                            ',/,&
'@    Wrong coordinates xmet, ymet for the meteo profile.     ',/,&
'@                                                            ',/,&
'@    See cs_user_parameters.f90                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1004 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@    Wrong xchem, ychem coordinates for the concentration    ',/,&
'@    profiles (chemistry model).                             ',/,&
'@                                                            ',/,&
'@    See cs_user_parameters.f90                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return

end subroutine atiniv1
