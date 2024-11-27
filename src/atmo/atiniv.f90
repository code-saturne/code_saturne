!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

procedure() :: atlecm, atlecc, atleca

! Local variables

integer          imode, iel
integer          k,ii, isc
double precision zent
double precision xcent

type(var_cal_opt) :: vcopt_p, vcopt_u

double precision, dimension(:), pointer :: cvar_despgi, cvar_sc

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

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

! Atmospheric gaseous chemistry
if (ichemistry.ge.1) then

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
         zproc, tchem, espnum, zent  , ttcabs, xcent )
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
if (     (iatra1.eq.1.or.ichemistry.ge.1)              &
    .and.(syear.eq.-1.or.squant.eq.-1.or.shour.eq.-1 &
          .or.smin.eq.-1.or.ssec.le.-1.d0)) then
  if (iatra1.eq.1) write(nfecra,1000)
  if (ichemistry.ge.1) write(nfecra,1001)
  call csexit (1)
endif

! Check radiative module latitude / longitude
if (iatra1.eq.1 .and. (xlat.ge.rinfin*0.5 .or. xlon.ge.rinfin*0.5)) then
  write(nfecra,1002)
  call csexit (1)
endif

! Check latitude / longitude from meteo file
if (imeteo.eq.1) then
  if (maxval(xyp_met(1,:)).ge.rinfin*0.5 .or. maxval(xyp_met(2,:)).ge.rinfin*0.5) then
    write(nfecra,1003)
    call csexit (1)
  endif
endif

! Check latitude / longitude from chemistry file
if (ichemistry.ge.1) then
  if (maxval(xchem).ge.rinfin*0.5 .or. maxval(ychem).ge.rinfin*0.5) then
    write(nfecra,1004)
    call csexit (1)
  endif
endif

!--------
! Formats
!--------

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
'@    Check your data and parameters (GUI and user functions) ',/,&
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
'@    Check your data and parameters (GUI and user functions) ',/,&
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
'@    Check your data and parameters (GUI and user functions) ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return

end subroutine atiniv0

