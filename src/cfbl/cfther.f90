!-------------------------------------------------------------------------------

!VERS

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
! File description:
! ----------------

!> \file cfther.f90
!>
!> Define thermodynamic laws (especially for the compressible flow scheme).
!>
!> Only the perfect gas law is available for now. The molar mass has to be
!> provided either in the GUI or in cs_user_parameters.f90.

!===============================================================================


!-------------------------------------------------------------------------------
!> \brief This subroutine is a driver allowing to call the appropriate
!> thermodynamicalfunctions depending on the quantities provided by the user.
!> Hence it is only used during the initialization step and at the boundaries
!> of type supersonic inlet. It is described in the following how to select the
!> quantity to be returned.
!>
!> When calling the user subroutine, the integer 'iccfth' specifies which
!> calculation has to be performed (and which quantity has to be returned).
!> The values for 'iccfth' for each case are provided below.
!>
!>   The variables are referred to using a different index i:
!>
!>     - pressure: 2
!>     - density: 3
!>     - temperature: 5
!>     - internal energy: 7
!>     - entropy: 13
!>
!>   iccfth is as follows, depending on which quantity needs to be computed:
!>     - variables at cell centers from variable i and variable j (i<j):
!>           iccfth = i*j*10000
!>     - variables at boundary faces from variable i and variable j (i<j):
!>           iccfth = i*j*10000+900
!>
!> Detailed values of iccfth and corresponding computations:
!>
!>   Values at the cell centers:
!>
!>     - temperature and energy from pressure and density: iccfth =  60000
!>     - density and energy from pressure and temperature: iccfth =  100000
!>     - density and temperature from pressure and energy: iccfth =  140000
!>     - pressure and energy from density and temperature: iccfth =  150000
!>     - pressure and temperature from density and energy: iccfth =  210000
!>
!>   Values at the faces for boundary conditions:
!>     - temperature and energy from pressure and density: iccfth = 60900
!>     - density and energy from pressure and temperature: iccfth = 100900
!>     - density and temperature from pressure and energy: iccfth = 140900
!>     - pressure and energy from density and temperature: iccfth = 150900
!>     - pressure and temperature from density and energy: iccfth = 210900
!>
!> \param[in]     nvar          total number of variables
!> \param[in]     iccfth        id of computation
!> \param[in]     imodif        modification indicator
!>                              - 0 variables and properties not modified
!>                              - > 0 variables and properties modified for a
!>                                bulk computation / face number for a B.C.
!> \param[out]    output1       computed therm. property 1 at cell centers
!> \param[out]    output2       computed therm. property 2 at cell centers
!> \param[in,out] bval          variable values at boundary faces
!-------------------------------------------------------------------------------

subroutine cfther(nvar, iccfth, imodif, output1, output2, bval)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use parall
use pointe
use entsor
use ppppar
use ppthch
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar, iccfth, imodif

double precision output1(*), output2(*), bval(nfabor,nvar)

! Local variables

integer          ifac0, l_size, iel, ifac, itk, ien

double precision, dimension(:), pointer :: crom, brom
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_pr, cvar_tk, cvar_en

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

!===============================================================================
! 0. Initialization.
!===============================================================================

! Rank of the variables in their associated arrays
if (iccfth.ge.0) then
  call field_get_val_s(icrom, crom)
  call field_get_val_s(ibrom, brom)
  itk = isca(itempk)
  ien = isca(ienerg)
  call field_get_val_s(ivarfl(itk), cvar_tk)
  call field_get_val_s(ivarfl(ien), cvar_en)
endif

call field_get_val_s(ivarfl(ipr), cvar_pr)

! For calculation of values at the cell faces,
! ifac0 is the number of the current face
ifac0 = imodif

! Calculation of temperature and energy from pressure and density
if (iccfth.eq.60000) then

  call cf_check_density(crom, ncel)

  call cf_thermo_te_from_dp(cvar_pr, crom, output1, output2, vel, ncel)

  ! Transfer to cell fields
  if (imodif.gt.0) then
    do iel = 1, ncel
      cvar_tk(iel) = output1(iel)
      cvar_en(iel) = output2(iel)
    enddo
  endif


! Calculation of density and energy from pressure and temperature:
elseif (iccfth.eq.100000) then

  call cf_check_temperature(cvar_tk, ncel)

  call cf_thermo_de_from_pt(cvar_pr, cvar_tk, output1, output2, vel, ncel)

  ! Transfer to the cell fields
  if (imodif.gt.0) then
    do iel = 1, ncel
      crom(iel) = output1(iel)
      cvar_en(iel) = output2(iel)
    enddo
  endif

! Calculation of density and temperature from pressure and energy
elseif (iccfth.eq.140000) then

  call cf_thermo_dt_from_pe(cvar_pr, cvar_en, output1, output2, vel, ncel)

  ! Transfer to cell fields
  if (imodif.gt.0) then
    do iel = 1, ncel
      crom(iel) = output1(iel)
      cvar_tk(iel) = output2(iel)
    enddo
  endif


! Calculation of pressure and energy from density and temperature
elseif (iccfth.eq.150000) then

  call cf_thermo_pe_from_dt(crom, cvar_tk, cvar_pr, cvar_en, vel, ncel )

  ! Transfer to cell fields
  if (imodif.gt.0) then
    do iel = 1, ncel
      cvar_pr(iel) = output1(iel)
      cvar_en(iel) = output2(iel)
    enddo
  endif

! Calculation of pressure and temperature from density and energy
elseif (iccfth.eq.210000) then

  call cf_thermo_pt_from_de(crom, cvar_en, output1, output2,    &
                            vel, ncel)

  ! Transfer to cell fields
  if (imodif.gt.0) then
    do iel = 1, ncel
      cvar_pr(iel) = output1(iel)
      cvar_tk(iel) = output2(iel)
    enddo
  endif

! Calculation of temperature and energy from pressure and density
! (it is postulated that the pressure and density values are strictly positive)
elseif (iccfth.eq.60900) then

  ifac = ifac0
  l_size = 1

  call cf_thermo_te_from_dp_ni(bval(ifac,ipr), brom(ifac:ifac), bval(ifac,itk), &
                               bval(ifac,ien), bval(ifac,iu), bval(ifac,iv),    &
                               bval(ifac,iw), l_size)


! Calculation of density and energy from pressure and temperature
elseif (iccfth.eq.100900) then

  ifac = ifac0
  l_size = 1

  call cf_thermo_de_from_pt_ni(bval(ifac,ipr), bval(ifac, itk), brom(ifac:ifac), &
                               bval(ifac,ien), bval(ifac,iu), bval(ifac,iv),     &
                               bval(ifac,iw), l_size)


! Calculation of density and temperature from pressure and total energy
elseif (iccfth.eq.140900) then

  ifac = ifac0
  l_size = 1

  call cf_thermo_dt_from_pe_ni(bval(ifac,ipr), bval(ifac,ien), brom(ifac:ifac), &
                               bval(ifac,itk), bval(ifac,iu), bval(ifac,iv),    &
                               bval(ifac,iw), l_size)

! Calculation of pressure and energy from density and temperature
elseif (iccfth.eq.150900) then

  ifac = ifac0
  l_size = 1

  call cf_thermo_pe_from_dt_ni(brom(ifac:ifac), bval(ifac,itk), bval(ifac,ipr),  &
                               bval(ifac,ien), bval(ifac,iu), bval(ifac,iv),     &
                               bval(ifac,iw), l_size)

! Calculation of pressure and temperature from density and energy
elseif (iccfth.eq.210900) then

  ifac = ifac0
  l_size = 1

  call cf_thermo_pt_from_de_ni(brom(ifac:ifac), bval(ifac,ien), bval(ifac,ipr),  &
                               bval(ifac,itk), bval(ifac,iu), bval(ifac,iv),     &
                               bval(ifac,iw), l_size)

endif

return
end subroutine cfther

!===============================================================================

!-------------------------------------------------------------------------------
! Each calculation has to be explicitly implemented in the appropriate
! subroutine below (already done for perfect gas).
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> \brief Set variability of isobaric specific heat and isochoric specific heat
!> according to the chosen thermodynamic law.
!-------------------------------------------------------------------------------

subroutine cf_set_thermo_options

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use ppincl

!===============================================================================

implicit none

! Arguments

!===============================================================================

if (ieos.eq.1) then

   ! --- Calculation options: constant Cp and Cv (perfect gas)
   ! The value for the isobaric specific heat Cp0 must be provided in
   ! the user subroutine 'usipsu'. The value for the isochoric
   ! specific heat Cv0 is calculated in a subsequent section (from Cp0)

    icp = 0
    icv = 0

endif

end subroutine cf_set_thermo_options

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Retrieve molar mass.

!> \param[out]   xmasml    molar mass
!-------------------------------------------------------------------------------

subroutine cf_get_molar_mass( xmasml )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use ppincl
use cstphy

!===============================================================================

implicit none

! Arguments

double precision xmasml

!===============================================================================

if (ieos.eq.1) then

   ! Molar mass of the gas (kg/mol)
   xmasml = xmasmr

endif

end subroutine cf_get_molar_mass

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute \f$\gamma\f$.

!> \param[out]   gamma    ratio of specific heat
!-------------------------------------------------------------------------------

subroutine cf_thermo_gamma(gamma)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use cstphy
use ppincl
use ppthch
use entsor

!===============================================================================

implicit none

! Arguments

double precision gamma

! Local variables

double precision gamagp, xmasml

!===============================================================================

call cf_get_molar_mass(xmasml)

! Gamagp is supposed to be superior or equal to 1.
! It is computed at each call, even if this may seem costly,
! to be coherent with the "constant gamma" case for which this
! constant is not saved. A ''save'' instruction and a test would
! be sufficient to avoid computing gamagp at each call if necessary.

if (ieos.eq.1) then

  gamagp = 1.d0 + rr/(xmasml*cp0-rr)

  if (gamagp.lt.1.d0) then
    write(nfecra,1010) gamagp
    call csexit (1)
  endif

endif

gamma = gamagp

!--------
! Formats
!--------

 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in thermodynamic computations      ',/,  &
'@       (cfther.f90), for perfect gas with constant gamma.',/,   &
'@',/,                                                            &
'@     Gamma = ',e12.4   ,/,                                      &
'@     Gamma must be a real number greater or equal to 1.',/,     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

return
end subroutine cf_thermo_gamma

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Initialize density, total energy and isochoric specific heat
!> according to the chosen thermodynamic law using the default parameters.

!> \param[in]     ncel    number of cells
!> \param[in]     ncelet  total number of cells on the local rank
!-------------------------------------------------------------------------------

subroutine cf_thermo_default_init(ncel, ncelet)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstphy
use numvar
use optcal
use ppincl
use ppthch
use field

!===============================================================================

implicit none

! Arguments

integer ncel, ncelet

! Local variables

integer ien, iel

double precision xmasml

double precision, dimension(:), pointer :: crom, cvar_en

!===============================================================================

! Default initializations
! t0 is positive (this assumption has been checked in verini)

ien = isca(ienerg)
call field_get_val_s(ivarfl(ien), cvar_en)
call field_get_val_s(icrom, crom)

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  cv0 = cp0 - rr/xmasml

  if (isuite.eq.0) then
    do iel = 1, ncel
      crom(iel) = p0 * xmasml/(rr*t0)
      cvar_en(iel) = cv0 * t0
    enddo
  endif

endif

end subroutine cf_thermo_default_init

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Check positivity of the pressure.

!> \param[in]     pres    array of pressure values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_check_pressure(pres, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use cstnum
use numvar
use optcal
use parall
use entsor

!===============================================================================

implicit none

! Arguments

integer l_size

double precision pres(*)

! Local variables

integer iel, ierr

!===============================================================================

! If the pressure is lower or equal to zero: clipping, write and stop.
! Indeed, if this is the case, the thermodynamic computations will most
! probably fail. This call is done at the end of the density calculation (after
! a classical clipping and before parallel communications).

ierr = 0
do iel = 1, l_size
  if (pres(iel).le.0.d0) then
    pres(iel) = epzero
    ierr = ierr + 1
  endif
enddo

if (irangp.ge.0) then
  call parcpt (ierr)
endif

if (ierr.gt.0) then
  ntmabs = ntcabs
  write(nfecra,8000)ierr, epzero
endif

!--------
! Formats
!--------

8000 format (                                                    &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in thermodynamic computations      ',/,  &
'@       (cfther.f90), for perfect gas with constant gamma.',/,   &
'@',/,                                                            &
'@     Negative values of the density were encountered ',/,       &
'@     in ',i10   ,' cells.',/,                                   &
'@     The density was clipped at ',e12.4  ,/                     &
'@     The run was stopped.',/,                                   &
'@',/,                                                            &
'@     If it is desired to continue the run in spite of this ',/, &
'@     behavior, it is possible to force a standard clipping ',/, &
'@     by setting a minimum value for the density variable in',/, &
'@     the GUI or in the user subroutine ''usipsu'' (set the ',/, &
'@     scamin value associated to the variable ',/,               &
'@     isca(irho).',/,                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

end subroutine cf_check_pressure

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Check strict positivity of the internal energy.

!> \param[in]     ener    array of total energy values
!> \param[in]     l_size  l_size of the array
!> \param[in]     vel     array of velocity values
!-------------------------------------------------------------------------------

subroutine cf_check_internal_energy(ener, l_size, vel)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use cstnum
use numvar
use optcal
use parall
use ppincl
use entsor

!===============================================================================

implicit none

! Arguments

integer l_size

double precision ener(*), vel(3,*)

! Local variables

integer iel, ierr

double precision enint

!===============================================================================

! If the internal energy <= zero: clipping, write and stop
!   Indeed, if this is the case, the thermodynamic computations will
!   most probably fail.

ierr = 0
do iel = 1, l_size
  enint = ener(iel) - 0.5d0*(  vel(1,iel)**2  &
                             + vel(2,iel)**2  &
                             + vel(3,iel)**2)
  if (enint.le.0.d0) then
    ener(iel) = epzero + 0.5d0*(  vel(1,iel)**2  &
                                + vel(2,iel)**2  &
                                + vel(3,iel)**2)
    ierr = ierr + 1
  endif
enddo

if (irangp.ge.0) then
  call parcpt (ierr)
endif

if (ierr.gt.0) then
  ntmabs = ntcabs
  write(nfecra,8100) ierr, epzero
endif

!--------
! Formats
!--------

8100 format (                                                    &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in thermodynamic computations      ',/,  &
'@       (cfther.f90), for perfect gas with constant gamma.',/,   &
'@',/,                                                            &
'@     Negative values of the internal energy were encountered',/,&
'@     in ',i10   ,' cells.',/,                                   &
'@     The internal energy  was clipped at ',e12.4  ,/            &
'@     The run was stopped.',/,                                   &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

end subroutine cf_check_internal_energy

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Check strict positivity of the density given by the user.

!> \param[in]     dens    array of density values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_check_density(dens, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use entsor

!===============================================================================

implicit none

! Arguments

integer l_size

double precision dens(*)

! Local variables

integer iel, ierr

!===============================================================================

! Verification of the values of the density
! Stop if a negative value is detected (since the density has been
! provided by the user, one potential cause is a wrong user
! initialization)

ierr = 0
do iel = 1, l_size
  if (dens(iel).le.0.d0) then
    write(nfecra,3010) dens(iel),iel
    ierr = ierr + 1
  endif
enddo

if (irangp.ge.0) then
  call parcpt (ierr)
endif

if (ierr.gt.0) then
  call csexit (1)
endif

!--------
! Formats
!--------

3010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in thermodynamic computations      ',/,  &
'@       (cfther.f90), for perfect gas with constant gamma.',/,   &
'@',/,                                                            &
'@     The computation of temperature failed.',/,                 &
'@',/,                                                            &
'@     Density = ',e12.4   ,' in cell ',i10  ,/,                  &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

end subroutine cf_check_density

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Check strict positivity of the temperature (Celsius) given by the user.

!> \param[in]     temp    array of temperature values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_check_temperature(temp, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use entsor

!===============================================================================

implicit none

! Arguments

integer l_size

double precision temp(*)

! Local variables

integer iel, ierr

!===============================================================================

! Verification of the values of the temperature
! Stop if a negative value is detected (since the temperature has been
! provided by the user, one potential cause is a wrong user
! initialization)

ierr = 0
do iel = 1, l_size
  if (temp(iel).le.0.d0) then
    write(nfecra,2010) temp(iel), iel
    ierr = ierr + 1
  endif
enddo

if (irangp.ge.0) then
  call parcpt(ierr)
endif

if (ierr.gt.0) then
  call csexit (1)
endif

!--------
! Formats
!--------

2010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in thermodynamic computations      ',/,  &
'@       (cfther.f90), for perfect gas with constant gamma.',/,   &
'@',/,                                                            &
'@     The computation of density failed.',/,                     &
'@',/,                                                            &
'@     Temperature = ',e12.4   ,' in cell ',i10  ,/,              &
'@     Temperature must be strictly positive.',/,                 &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

end subroutine cf_check_temperature

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute temperature and total energy from density and pressure.

!> \param[in]     pres    array of pressure values
!> \param[in]     dens    array of density values
!> \param[out]    temp    array of temperature values
!> \param[out]    ener    array of total energy values
!> \param[in]     velx    array of velocity x-component values
!> \param[in]     velz    array of velocity y-component values
!> \param[in]     vely    array of velocity z-component values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_te_from_dp_ni(pres, dens, temp, ener, velx, vely, velz, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use dimens
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision pres(*), dens(*), temp(*), ener(*)
double precision velx(*), vely(*), velz(*)

! local variables

integer ii

double precision xmasml

!===============================================================================

! calculation of temperature and energy from pressure and density

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  do ii = 1, l_size
    ! temperature
    temp(ii) = xmasml * pres(ii) / (rr*dens(ii))
    ! total energy
    ener(ii) =  cv0*temp(ii)                                             &
              + 0.5d0*( velx(ii)**2 + vely(ii)**2 + velz(ii)**2 )
  enddo

endif

end subroutine cf_thermo_te_from_dp_ni

!-------------------------------------------------------------------------------
!> \brief Compute temperature and total energy from density and pressure.

!> \param[in]     pres    array of pressure values
!> \param[in]     dens    array of density values
!> \param[out]    temp    array of temperature values
!> \param[out]    ener    array of total energy values
!> \param[in]     vel     array of velocity component values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_te_from_dp(pres, dens, temp, ener, vel, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use dimens
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision pres(*), dens(*), temp(*), ener(*)
double precision vel(3,*)

! local variables

integer ii

double precision xmasml

!===============================================================================

! calculation of temperature and energy from pressure and density

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  do ii = 1, l_size
    ! temperature
    temp(ii) = xmasml * pres(ii) / (rr*dens(ii))
    ! total energy
    ener(ii) =  cv0*temp(ii)                                             &
              + 0.5d0*( vel(1,ii)**2 + vel(2,ii)**2 + vel(3,ii)**2 )
  enddo

endif

end subroutine cf_thermo_te_from_dp

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute density and total energy from pressure and temperature.

!> \param[in]     pres    array of pressure values
!> \param[in]     temp    array of temperature values
!> \param[out]    dens    array of density values
!> \param[out]    ener    array of total energy values
!> \param[in]     velx    array of velocity x-component values
!> \param[in]     velz    array of velocity y-component values
!> \param[in]     vely    array of velocity z-component values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_de_from_pt_ni(pres, temp, dens, ener, velx, vely, velz, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision pres(*), temp(*), dens(*), ener(*)
double precision velx(*), vely(*), velz(*)

! Local variables

integer iel

double precision xmasml

!===============================================================================

! Calculation of density and energy from pressure and temperature

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  do iel = 1, l_size
    ! Temperature
    dens(iel) = xmasml * pres(iel) / (rr*temp(iel))
    ! Total energy
    ener(iel) =  cv0*temp(iel)                                             &
               + 0.5d0*( velx(iel)**2 + vely(iel)**2 + velz(iel)**2 )
  enddo

endif

end subroutine cf_thermo_de_from_pt_ni

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute density and total energy from pressure and temperature
!>        (interleaved version).

!> \param[in]     pres    array of pressure values
!> \param[in]     temp    array of temperature values
!> \param[out]    dens    array of density values
!> \param[out]    ener    array of total energy values
!> \param[in]     vel     array of velocity component values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_de_from_pt(pres, temp, dens, ener, vel, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision pres(*), temp(*), dens(*), ener(*)
double precision vel(3,*)

! Local variables

integer iel

double precision xmasml

!===============================================================================

! Calculation of density and energy from pressure and temperature

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  do iel = 1, l_size
    ! Temperature
    dens(iel) = xmasml * pres(iel) / (rr*temp(iel))
    ! Total energy
    ener(iel) =  cv0*temp(iel)                                             &
               + 0.5d0*( vel(1,iel)**2 + vel(2,iel)**2 + vel(3,iel)**2 )
  enddo

endif

end subroutine cf_thermo_de_from_pt

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute density and temperature from pressure and total energy.

!> \param[in]     pres    array of pressure values
!> \param[in]     ener    array of total energy values
!> \param[out]    dens    array of density values
!> \param[out]    temp    array of temperature values
!> \param[in]     velx    array of velocity x-component values
!> \param[in]     velz    array of velocity y-component values
!> \param[in]     vely    array of velocity z-component values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_dt_from_pe_ni(pres, ener, dens, temp, velx, vely, velz, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision pres(*), ener(*), dens(*), temp(*)
double precision velx(*), vely(*), velz(*)

! Local variables

integer iel

double precision enint, gamagp, xmasml

!===============================================================================

! Calculation of density and temperature from pressure and energy

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  call cf_thermo_gamma(gamagp)

  do iel = 1, l_size
    ! Internal energy (to avoid the need to divide by the temperature
    ! to compute density)
    enint =  ener(iel)                                      &
           - 0.5d0*( velx(iel)**2                           &
           + vely(iel)**2                                   &
           + velz(iel)**2 )
    ! Density
    dens(iel) = pres(iel) / ( (gamagp-1.d0) * enint )
    ! Temperature
    temp(iel) = xmasml * (gamagp-1.d0) * enint / rr
  enddo

endif

end subroutine cf_thermo_dt_from_pe_ni

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute density and temperature from pressure and total energy
!>        (interleaved version).

!> \param[in]     pres    array of pressure values
!> \param[in]     ener    array of total energy values
!> \param[out]    dens    array of density values
!> \param[out]    temp    array of temperature values
!> \param[in]     vel     array of velocity component values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_dt_from_pe(pres, ener, dens, temp, vel, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision pres(*), ener(*), dens(*), temp(*)
double precision vel(3,*)

! Local variables

integer iel

double precision enint, gamagp, xmasml

!===============================================================================

! Calculation of density and temperature from pressure and energy

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  call cf_thermo_gamma(gamagp)

  do iel = 1, l_size
    ! Internal energy (to avoid the need to divide by the temperature
    ! to compute density)
    enint =  ener(iel)                                      &
           - 0.5d0*( vel(1,iel)**2                           &
           + vel(2,iel)**2                                   &
           + vel(3,iel)**2 )
    ! Density
    dens(iel) = pres(iel) / ( (gamagp-1.d0) * enint )
    ! Temperature
    temp(iel) = xmasml * (gamagp-1.d0) * enint / rr
  enddo

endif

end subroutine cf_thermo_dt_from_pe

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute pressure and total energy from density and temperature.

!> \param[in]     dens    array of density values
!> \param[in]     temp    array of temperature values
!> \param[out]    pres    array of pressure values
!> \param[out]    ener    array of total energy values
!> \param[in]     velx    array of velocity x-component values
!> \param[in]     velz    array of velocity y-component values
!> \param[in]     vely    array of velocity z-component values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_pe_from_dt_ni(dens, temp, pres, ener, velx, vely, velz, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision dens(*), temp(*), pres(*), ener(*)
double precision velx(*), vely(*), velz(*)

! Local variables

integer iel

double precision xmasml

!===============================================================================

! Calculation of pressure and energy from density and temperature

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  do iel = 1, l_size
    ! Pressure
    pres(iel) = dens(iel)*temp(iel)*rr/xmasml
    ! Total energy
    ener(iel) =  cv0*temp(iel)                                          &
               + 0.5d0*( velx(iel)**2 + vely(iel)**2 + velz(iel)**2 )
  enddo

endif

end subroutine cf_thermo_pe_from_dt_ni

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute pressure and total energy from density and temperature
!>        (interleaved version).

!> \param[in]     dens    array of density values
!> \param[in]     temp    array of temperature values
!> \param[out]    pres    array of pressure values
!> \param[out]    ener    array of total energy values
!> \param[in]     vel     array of velocity component values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_pe_from_dt(dens, temp, pres, ener, vel, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision dens(*), temp(*), pres(*), ener(*)
double precision vel(3,*)

! Local variables

integer iel

double precision xmasml

!===============================================================================

! Calculation of pressure and energy from density and temperature

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  do iel = 1, l_size
    ! Pressure
    pres(iel) = dens(iel)*temp(iel)*rr/xmasml
    ! Total energy
    ener(iel) =  cv0*temp(iel)                                          &
               + 0.5d0*( vel(1,iel)**2 + vel(2,iel)**2 + vel(3,iel)**2 )
  enddo

endif

end subroutine cf_thermo_pe_from_dt

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute pressure and temperature from density and total energy.

!> \param[in]     dens    array of density values
!> \param[in]     temp    array of temperature values
!> \param[out]    pres    array of pressure values
!> \param[out]    ener    array of total energy values
!> \param[in]     velx    array of velocity x-component values
!> \param[in]     velz    array of velocity y-component values
!> \param[in]     vely    array of velocity z-component values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_pt_from_de_ni(dens, ener, pres, temp, velx, vely, velz, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision dens(*), ener(*), pres(*), temp(*)
double precision velx(*), vely(*), velz(*)

! Local variables

integer iel

double precision enint, gamagp, xmasml

!===============================================================================

! Calculation of pressure and temperature from density and energy

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  call cf_thermo_gamma(gamagp)

  do iel = 1, l_size
    ! Internal energy (to avoid the need to divide by the temperature
    ! to compute density)
    enint =  ener(iel)                                        &
           - 0.5d0*( velx(iel)**2                             &
           + vely(iel)**2                             &
           + velz(iel)**2 )

    ! Pressure
    pres(iel) = (gamagp-1.d0) * dens(iel) * enint
    ! Temperature
    temp(iel) = xmasml * (gamagp-1.d0) * enint / rr
  enddo

endif

end subroutine cf_thermo_pt_from_de_ni

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute pressure and temperature from density and total energy.

!> \param[in]     dens    array of density values
!> \param[in]     temp    array of temperature values
!> \param[out]    pres    array of pressure values
!> \param[out]    ener    array of total energy values
!> \param[in]     vel     array of velocity component values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_pt_from_de(dens, ener, pres, temp, vel, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision dens(*), ener(*), pres(*), temp(*)
double precision vel(3,*)

! Local variables

integer iel

double precision enint, gamagp, xmasml

!===============================================================================

! Calculation of pressure and temperature from density and energy

call cf_get_molar_mass(xmasml)

if (ieos.eq.1) then

  call cf_thermo_gamma(gamagp)

  do iel = 1, l_size
    ! Internal energy (to avoid the need to divide by the temperature
    ! to compute density)
    enint =  ener(iel)                                        &
           - 0.5d0*( vel(1,iel)**2                             &
           + vel(2,iel)**2                             &
           + vel(3,iel)**2 )

    ! Pressure
    pres(iel) = (gamagp-1.d0) * dens(iel) * enint
    ! Temperature
    temp(iel) = xmasml * (gamagp-1.d0) * enint / rr
  enddo

endif

end subroutine cf_thermo_pt_from_de

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute square of sound velocity:
!>
!> \f[c^2  = \left(\frac{\partial p}{\partial \rho}\right)_s\f]
!>
!> for perfect gas, the explicit formula is:
!>
!> \f[c^2  = \gamma \frac{p}{\rho}\f]

!> \param[in]    pres    array of pressure values
!> \param[in]    dens    array of density values
!> \param[out]   c2      array of the values of the square of sound velocity
!> \param[in]    l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_c_square(pres, dens, c2, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision pres(*), dens(*), c2(*)

! Local variables

integer iel

double precision gamagp

!===============================================================================

if (ieos.eq.1) then

  call cf_thermo_gamma(gamagp)

  do iel = 1, l_size
    c2(iel) = gamagp * pres(iel) / dens(iel)
  enddo

endif

end subroutine cf_thermo_c_square

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute \f$\beta\f$:
!>
!> \f[ \beta = \left(\frac{\partial p}{\partial s}\right)_\rho \f]
!>
!> for a perfect gas, the explicit formula is:
!>
!> \f[ \beta = \rho^\gamma \f]

!> \param[in]    dens    array of density values
!> \param[out]   beta    array of beta values
!> \param[in]    l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_beta(dens, beta, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision dens(*), beta(*)

! Local variables

integer iel

double precision gamagp


if (ieos.eq.1) then

  call cf_thermo_gamma(gamagp)

  do iel = 1, l_size
    beta(iel) = dens(iel)**gamagp
  enddo

endif

end subroutine cf_thermo_beta

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute the isochoric specific heat:
!>
!> \f[C_v = \left(\frac{\partial e}{\partial T}\right)_\rho\f]

!> \param[in]     cp      array of isobaric specific heat values
!> \param[in]     cv      array of isochoric specific heat values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_cv(cp, cv, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision cp(*), cv(*)

! Local variables

! Constant quantity if perfect gas chosen,
! nothing to be done.

end subroutine cf_thermo_cv

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute entropy from pressure and density:
!>
!> \f[s = \frac{p}{\rho^\gamma}\f]

!> \param[in]     dens    array of density values
!> \param[in]     pres    array of pressure values
!> \param[out]    entr    array of total energy values
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_s_from_dp(dens, pres, entr, l_size)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use parall
use ppincl
use ppthch

!===============================================================================

implicit none

! Arguments

integer l_size

double precision dens(*), pres(*), entr(*)

! Local variables

integer iel

double precision gamagp

!===============================================================================

if (ieos.eq.1) then

  call cf_check_density(dens, l_size)

  call cf_thermo_gamma(gamagp)

  do iel = 1, l_size
    entr(iel) = pres(iel) / (dens(iel)**gamagp)
  enddo

endif

end subroutine cf_thermo_s_from_dp

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute wall boundary conditions.

!> \param[in,out] wbfb    output work array
!> \param[in]     ifac    boundary face indice
!-------------------------------------------------------------------------------

subroutine cf_thermo_wall_bc(wbfb, ifac)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use cstnum
use paramx
use mesh
use numvar
use parall
use ppincl
use ppthch
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer ifac

double precision wbfb(nfabor)

! Local variables

integer iel

double precision gamagp, xmach

double precision, dimension(:), pointer :: crom
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_pr

!===============================================================================

if (ieos.eq.1) then

  ! Calculation of the boundary conditions on the wall face ifac

  call cf_thermo_gamma(gamagp)

  ! Map field arrays
  call field_get_val_v(ivarfl(iu), vel)
  call field_get_val_s(icrom, crom)

  call field_get_val_s(ivarfl(ipr), cvar_pr)

  iel = ifabor(ifac)

  ! Calculation of the Mach number at the boundary face, using the
  ! cell center velocity projected on the vector normal to the boundary
  xmach = ( vel(1,iel)*surfbo(1,ifac)                          &
          + vel(2,iel)*surfbo(2,ifac)                          &
          + vel(3,iel)*surfbo(3,ifac) ) / surfbn(ifac)         &
         / sqrt( gamagp * cvar_pr(iel) / crom(iel) )

  ! Pressure

  ! A Neumann boundary condition is used. This does not allow to use
  ! the Rusanov scheme, but some stabilization effect is expected.
  ! A test based on the value of coefb at the previous time step
  ! is implemented to avoid oscillating between a rarefaction
  ! situation and a shock configuration from one time step to the
  ! next.

  ! Rarefaction !FIXME with the new cofaf cofbf
  if (xmach.lt.0.d0.and.wbfb(ifac).le.1.d0) then

     if (xmach.gt.2.d0/(1.d0-gamagp)) then
       wbfb(ifac) = (1.d0 + (gamagp-1.d0)/2.d0 * xmach)    &
                    **(2.d0*gamagp/(gamagp-1.d0))
     else
        ! In case the rarefaction is too strong, a zero Dirichlet value
        ! is used for pressure (the value of wbfb is used here as an
        ! indicator)
        wbfb(ifac) = rinfin
     endif

  ! Shock
  elseif (xmach.gt.0.d0.and.wbfb(ifac).ge.1.d0) then

    wbfb(ifac) = 1.d0 + gamagp*xmach                                    &
                       *((gamagp+1.d0)/4.d0*xmach                       &
                         + sqrt(1.d0 + (gamagp+1.d0)**2/16.d0*xmach**2) )

  ! Oscillation between rarefaction and shock or zero Mach number
  else
    wbfb(ifac) = 1.d0
  endif

endif

end subroutine cf_thermo_wall_bc

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute subsonic outlet boundary conditions.

!> \param[in,out] bval    variable values at boundary faces
!> \param[in]     ifac    boundary face indice
!-------------------------------------------------------------------------------

subroutine cf_thermo_subsonic_outlet_bc(bval, ifac)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use mesh
use paramx
use numvar
use parall
use ppincl
use ppthch
use field

!===============================================================================

implicit none

! Arguments

integer ifac

double precision bval(nfabor,*)

! Local variables

integer iel, ien

double precision gamagp
double precision roi, ro1, pri, ei, uni, un1, uns
double precision ci, c1, mi, a, b, sigma1, pinf

double precision, dimension(:), pointer :: crom, brom
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_pr, cvar_en

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

!===============================================================================

! Calculation of the boundary conditions on the subsonic outlet face ifac

ien = isca(ienerg)

if (ieos.eq.1) then

  call cf_thermo_gamma(gamagp)

  call field_get_val_s(icrom, crom)
  call field_get_val_s(ibrom, brom)

  call field_get_val_s(ivarfl(ipr), cvar_pr)
  call field_get_val_s(ivarfl(ien), cvar_en)

  iel = ifabor(ifac)

  pinf = bval(ifac,ipr)
  pri  = cvar_pr(iel)
  roi  = crom(iel)

  ci   = sqrt(gamagp * pri / roi)
  uni  = ( vel(1,iel) * surfbo(1,ifac)                             &
         + vel(2,iel) * surfbo(2,ifac)                             &
         + vel(3,iel) * surfbo(3,ifac) ) / surfbn(ifac)
  ei   = cvar_en(iel) - 0.5d0 * uni**2

  ! Rarefaction case
  if (pinf.le.pri) then

    ! Computation of the velocity in state 1 using Riemann invariants
    ! of the 1-rarefaction
    a = 2 * ci / (gamagp - 1.d0) * (1.d0 -  (pinf / pri) ** ( (gamagp - 1.d0) &
                                                             /(2.d0 * gamagp)))
    un1 = uni + a

    ! Computation of the density in state 1 using Rieman invariants
    ! of the 1-rarefaction
    ro1 = roi * (pinf/pri)**(1.d0 / gamagp)

    ! Subsonic inlet - state 2 should be imposed but too few information
    ! is available to compute it
    ! for want of anything better, state 1 is imposed
    if (un1.lt.0.d0) then

      ! Density
      brom(ifac) = ro1
      ! Velocity
      bval(ifac,iu) = vel(1,iel) + a * surfbo(1,ifac) / surfbn(ifac)
      bval(ifac,iv) = vel(2,iel) + a * surfbo(2,ifac) / surfbn(ifac)
      bval(ifac,iw) = vel(3,iel) + a * surfbo(3,ifac) / surfbn(ifac)
      ! Total energy
      bval(ifac,ien) = pinf / ((gamagp - 1.d0) * ro1)                       &
                     + 0.5d0 * (bval(ifac,iu)**2                            &
                     + bval(ifac,iv)**2                                     &
                     + bval(ifac,iw)**2)

      ! Outlet
    else

      ! Computation of the sound speed in state 1
      c1 = sqrt(gamagp * pinf / ro1)

      ! Subsonic outlet - state 1 is imposed
      if ((un1-c1).lt.0.d0) then

        ! Density
        brom(ifac) = ro1
        ! Velocity
        bval(ifac,iu) = vel(1,iel) + a * surfbo(1,ifac) / surfbn(ifac)
        bval(ifac,iv) = vel(2,iel) + a * surfbo(2,ifac) / surfbn(ifac)
        bval(ifac,iw) = vel(3,iel) + a * surfbo(3,ifac) / surfbn(ifac)
        ! Total energy
        bval(ifac,ien) = pinf / ((gamagp - 1.d0) * ro1)                     &
                       + 0.5d0 * (bval(ifac,iu)**2                          &
                       + bval(ifac,iv)**2                         &
                       + bval(ifac,iw)**2)

        ! Sonic outlet
      else if ((uni-ci).lt.0.d0) then

        ! Mach number in the domain
        mi = uni / ci

        b = (gamagp - 1.d0) / (gamagp + 1.d0) * (mi + 2.d0 / (gamagp - 1))

        ! Sonic state pressure
        bval(ifac,ipr) = pri * b ** (2.d0 * gamagp / (gamagp - 1.d0))
        ! Sonic state density
        brom(ifac) = roi * b ** (2.d0 / (gamagp - 1.d0))
        ! Sonic state velocity
        uns = b * ci
        bval(ifac,iu) = uns * surfbo(1,ifac) / surfbn(ifac)
        bval(ifac,iv) = uns * surfbo(2,ifac) / surfbn(ifac)
        bval(ifac,iw) = uns * surfbo(3,ifac) / surfbn(ifac)
        ! Sonic state energy
        bval(ifac,isca(ienerg)) = bval(ifac,ipr)/((gamagp - 1.d0) * brom(ifac))&
                                + 0.5d0 * uns**2

        ! Supersonic outlet
      else

        ! pb = pri
        bval(ifac,ipr) = pri
        ! ub = uni
        bval(ifac,iu) = vel(1,iel)
        bval(ifac,iv) = vel(2,iel)
        bval(ifac,iw) = vel(3,iel)
        ! rob = roi
        brom(ifac) = roi
        ! eb = ei
        bval(ifac,isca(ienerg)) = cvar_en(iel)

      endif

    endif

    ! Shock case
  else

    ! Computation of the density in state 1 with Rankine-Hugoniot relations
    ro1 = roi * ((gamagp - 1.d0) * pri  + (gamagp + 1.d0) * pinf)   &
              / ((gamagp - 1.d0) * pinf + (gamagp + 1.d0) * pri )

    ! Computation of the velocity in state 1 with Rankine-Hugoniot relations
    ! un1 = un2
    a = sqrt( (pinf - pri) * (1.d0/roi - 1.d0/ro1) )
    un1 = uni - a

    ! Subsonic inlet - state 2 should be imposed but too few information
    ! is available to compute it
    ! for want of anything better, state 1 is imposed
    if (un1.le.0d0) then

      ! Density
      brom(ifac) = ro1
      ! Velocity
      bval(ifac,iu) = vel(1,iel) - a * surfbo(1,ifac) / surfbn(ifac)
      bval(ifac,iv) = vel(2,iel) - a * surfbo(2,ifac) / surfbn(ifac)
      bval(ifac,iw) = vel(3,iel) - a * surfbo(3,ifac) / surfbn(ifac)
      ! Total energy
      bval(ifac,ien) = pinf / ((gamagp-1.d0) * brom(ifac))           &
                     + 0.5d0 * (bval(ifac,iu)**2                     &
                     + bval(ifac,iv)**2                              &
                     + bval(ifac,iw)**2)

    ! Outlet
    else

      ! Computation of the shock velocity
      sigma1 = (roi * uni - ro1 * un1) / (roi - ro1)

      ! Subsonic outlet - state 1 is imposed
      if (sigma1.le.0.d0) then

        ! Density
        brom(ifac) = ro1
        ! Velocity
        bval(ifac,iu) = vel(1,iel) - a * surfbo(1,ifac) / surfbn(ifac)
        bval(ifac,iv) = vel(2,iel) - a * surfbo(2,ifac) / surfbn(ifac)
        bval(ifac,iw) = vel(3,iel) - a * surfbo(3,ifac) / surfbn(ifac)
        ! Total energy
        bval(ifac,ien) =  pinf / ((gamagp-1.d0) * brom(ifac))  &
                        + 0.5d0 * (bval(ifac,iu)**2            &
                        + bval(ifac,iv)**2                     &
                        + bval(ifac,iw)**2)

        ! Supersonic outlet
      else

        ! pb = pri
        bval(ifac,ipr) = pri
        ! unb = uni
        bval(ifac,iu) = vel(1,iel)
        bval(ifac,iv) = vel(2,iel)
        bval(ifac,iw) = vel(3,iel)
        ! rob = roi
        brom(ifac) = roi
        ! eb = ei
        bval(ifac,isca(ienerg)) = cvar_en(iel)

      endif ! test on shock speed sign

    endif ! test on state 1 velocity sign

  endif ! test on pinf-pri sign

endif

end subroutine cf_thermo_subsonic_outlet_bc

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute inlet boundary condition with total pressure and total
!> enthalpy imposed.

!> \param[in,out] bval    variable values at boundary faces
!> \param[in]     ifac    boundary face number
!-------------------------------------------------------------------------------

subroutine cf_thermo_ph_inlet_bc(bval, ifac)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use mesh
use numvar
use optcal
use parall
use ppincl
use ppthch
use field
use cstnum
use entsor

!===============================================================================

implicit none

! Arguments

integer ifac

double precision bval(nfabor,*)

! Local variables

integer iel, niter, nitermax

double precision gamagp, bMach, eps, pstat, old_pstat, ptot, res, rhotot
double precision roi, ro1, pri, ei, uni, un1, y, uns, bc, cosalp, norm
double precision ci, c1, mi, a, sigma1, utxi, utyi, utzi
double precision dir(3)

double precision, dimension(:), pointer :: crom, brom
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_pr, cvar_en

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

!===============================================================================

if (ieos.eq.1) then

  ! Calculation of the boundary conditions on the inlet face ifac

  call cf_thermo_gamma(gamagp)

  call field_get_val_s(icrom, crom)
  call field_get_val_s(ibrom, brom)

  call field_get_val_s(ivarfl(ipr), cvar_pr)
  call field_get_val_s(ivarfl(isca(ienerg)), cvar_en)

  iel  = ifabor(ifac)

  niter = 0

  roi  = crom(iel)
  pri  = cvar_pr(iel)

  ! Normalize the direction vector given by the user
  norm = sqrt(bval(ifac,iu)**2 + bval(ifac,iv)**2 + bval(ifac,iw)**2)
  if (norm.lt.epzero) then
    write(nfecra,9010)ifac
    call csexit (1)
  endif

  dir(1) = bval(ifac,iu) / norm
  dir(2) = bval(ifac,iv) / norm
  dir(3) = bval(ifac,iw) / norm

  ! Angle between the imposed direction and the inlet normal
  cosalp = ( dir(1)*surfbo(1,ifac) + dir(2)*surfbo(2,ifac)  &
           + dir(3)*surfbo(3,ifac) ) /surfbn(ifac)

  ! If direction vector is outward, warn the user
  if (cosalp.gt.epzero) then
    write(nfecra,9020)ifac
  endif

  ! Computation of the sound speed inside the domain
  ci = sqrt(gamagp * pri / roi)

  uni = ( vel(1,iel) * surfbo(1,ifac)                             &
        + vel(2,iel) * surfbo(2,ifac)                             &
        + vel(3,iel) * surfbo(3,ifac) ) / surfbn(ifac)

  bMach = uni / ci

  utxi = vel(1,iel) - uni * surfbo(1,ifac) * surfbn(ifac)
  utyi = vel(2,iel) - uni * surfbo(2,ifac) * surfbn(ifac)
  utzi = vel(3,iel) - uni * surfbo(3,ifac) * surfbn(ifac)

  ei   = cvar_en(iel) - 0.5d0 * (  vel(1,iel)**2           &
                                 + vel(2,iel)**2           &
                                 + vel(3,iel)**2 )

  ptot = bval(ifac,ipr)
  rhotot = gamagp / (gamagp - 1.d0) * ptot / bval(ifac,isca(ienerg))
  old_pstat = ptot

  nitermax = 100
  eps = epsrsm(ipr)
  res = 1.d0

  do while (niter.le.100.and.res.gt.eps)

    pstat =  ptot*(1.d0+(gamagp - 1.d0)*0.5d0*bMach**2)**(gamagp/(1.d0-gamagp))
    y = pri / pstat

    ! 1-shock
    if (y.lt.1.d0) then

      ! Computation of the density in state 1 with Rankine-Hugoniot relations
      ro1 = roi * ((gamagp - 1.d0) * pri   + (gamagp + 1.d0) * pstat)  &
                / ((gamagp - 1.d0) * pstat + (gamagp + 1.d0) * pri)

      ! Computation of the velocity in state 1 with Rankine-Hugoniot relations
      ! un1 = un2
      un1 = uni - sqrt( (pstat - pri) * (1.d0/roi - 1.d0/ro1) )

      ! Subsonic inlet
      if (un1.le.0.d0) then

        ! unb = u2
        bval(ifac,iu) = un1 / cosalp * dir(1)
        bval(ifac,iv) = un1 / cosalp * dir(2)
        bval(ifac,iw) = un1 / cosalp * dir(3)
        ! rob = ro2
        brom(ifac) = (pstat / ptot)**(1.d0/gamagp) * rhotot
        ! eb = e2
        bval(ifac,isca(ienerg)) = pstat / ((gamagp - 1.d0) * brom(ifac))       &
             + 0.5d0 * (bval(ifac,iu)**2 + bval(ifac,iv)**2 + bval(ifac,iw)**2)
        ! Outlet
      else
        ! Computation of the shock velocity
        sigma1 = (roi * uni - ro1 * un1) / (roi - ro1)

        ! subsonic outlet
        if (sigma1.le.0.d0) then

          ! unb = u1
          bval(ifac,iu) = utxi + un1 * surfbo(1,ifac) / surfbn(ifac)
          bval(ifac,iv) = utyi + un1 * surfbo(2,ifac) / surfbn(ifac)
          bval(ifac,iw) = utzi + un1 * surfbo(3,ifac) / surfbn(ifac)
          ! rob = ro1
          brom(ifac) = ro1
          ! eb = e1
          bval(ifac,isca(ienerg)) = ei                                         &
                           - 0.5d0 * (pstat + pri) * (1.d0 / ro1 - 1.d0 / roi) &
                           + 0.5d0 * (un1**2 + utxi**2 + utyi**2 + utzi**2)

        ! supersonic outlet
        else

          ! pb = pri
          pstat = pri
          ! unb = uni
          bval(ifac,iu) = vel(1,iel)
          bval(ifac,iv) = vel(2,iel)
          bval(ifac,iw) = vel(3,iel)
          ! rob = roi
          brom(ifac) = roi
          ! eb = ei
          bval(ifac,isca(ienerg)) = cvar_en(iel)

        endif

      endif

      ! 1-rarefaction
    else

      ! Computation of the velocity in state 1 using Riemann invariants
      ! of the 1-rarefaction
      un1 = uni +  2 * ci / (gamagp - 1.d0)                                    &
                 * (1.d0 - (pstat / pri) ** ((gamagp - 1.d0) / (2.d0 * gamagp)))

      ! Computation of the density in state 1 using Riemann invariants
      ! of the 1-rarefaction
      ro1 = (pstat / pri) ** (1.d0 / gamagp) * roi

      ! Subsonic inlet
      if (un1.le.0.d0) then

        ! unb = u2
        bval(ifac,iu) = un1 / cosalp * dir(1)
        bval(ifac,iv) = un1 / cosalp * dir(2)
        bval(ifac,iw) = un1 / cosalp * dir(3)
        ! rob = ro2
        brom(ifac) = (pstat / ptot)**(1.d0/gamagp) * rhotot
        ! eb = e2
        bval(ifac,isca(ienerg)) = pstat / ((gamagp - 1.d0) * brom(ifac))      &
                                + 0.5d0 * (  bval(ifac,iu)**2                 &
                                           + bval(ifac,iv)**2                 &
                                           + bval(ifac,iw)**2 )
        ! Outlet
      else

        ! Computation of the sound speed in state 1
        c1 = sqrt(gamagp * pstat / ro1)

        ! Subsonic outlet
        if ((un1 - c1).lt.0.d0) then

          ! unb = u1
          bval(ifac,iu) = utxi + un1 * surfbo(1,ifac) / surfbn(ifac)
          bval(ifac,iv) = utyi + un1 * surfbo(2,ifac) / surfbn(ifac)
          bval(ifac,iw) = utzi + un1 * surfbo(3,ifac) / surfbn(ifac)
          ! rob = ro1
          brom(ifac) = ro1
          ! eb = e1
          bval(ifac,isca(ienerg)) = pstat / (ro1 * (gamagp - 1.d0))            &
                                  + 0.5d0 * (un1**2 + utxi**2 + utyi**2 + utzi**2)

          ! Supersonic outlet
        else if ((uni - ci).ge.0.d0) then

          ! pb = pri
          pstat = pri
          ! ub = uni
          bval(ifac,iu) = vel(1,iel)
          bval(ifac,iv) = vel(2,iel)
          bval(ifac,iw) = vel(3,iel)
          ! rob = roi
          brom(ifac) = roi
          ! eb = ei
          bval(ifac,isca(ienerg)) = cvar_en(iel)

          ! Outlet in sonic state
        else

          ! Mach number in the domain
          mi = uni / ci

          a = (gamagp - 1.d0) / (gamagp + 1.d0) * (mi + 2.d0 / (gamagp - 1))

          ! Sonic state pressure
          pstat = pri * a ** (2.d0 * gamagp / (gamagp - 1.d0))
          ! Sonic state density
          brom(ifac) = roi * a ** (2.d0 / (gamagp - 1.d0))
          ! Sonic state velocity
          uns = a * ci
          bval(ifac,iu) = uns * surfbo(1,ifac) / surfbn(ifac)
          bval(ifac,iv) = uns * surfbo(2,ifac) / surfbn(ifac)
          bval(ifac,iw) = uns * surfbo(3,ifac) / surfbn(ifac)
          ! Sonic state energy
          bval(ifac,isca(ienerg)) =  pstat / ((gamagp - 1.d0) * brom(ifac)) &
                                   + 0.5d0 * uns**2

        endif

      endif

    endif

    bc = sqrt(gamagp * pstat / brom(ifac))
    bMach = ( bval(ifac,iu) * surfbo(1,ifac)                             &
            + bval(ifac,iv) * surfbo(2,ifac)                             &
            + bval(ifac,iw) * surfbo(3,ifac) ) / surfbn(ifac) / bc

    bval(ifac,ipr) = pstat

    ! Pressure residual
    res = abs((pstat - old_pstat) / ptot)

    ! Prepare next iteration
    old_pstat = pstat
    niter = niter + 1

  enddo

  ! Warn the user if fixed point algorithm did not converge
  if (niter.eq.101) then
    write(nfecra,3000)ifac,res
  endif

endif

!--------
! Formats
!--------

3000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    in thermodynamics computations',/,              &
'@    =======',/,                                                 &
'@     Error encountered in thermodynamic computations      ',/,  &
'@       (cfther.f90), for perfect gas with constant gamma.',/,   &
'@',/,                                                            &
'@     Fixed point algorithm did not converge when',/,            &
'@     computing the subsonic inlet boundary condition',/,        &
'@     with total pressure and enthalpy imposed.',/,              &
'@',/,                                                            &
'@     At boundary face ',i10   , /,                              &
'@     boundary Mach number residual = ', e12.4 ,/,               &
'@     maximum number of iterations (100) was reached.' ,/,       &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

9010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    in thermodynamics computations',/,              &
'@    =======',/,                                                 &
'@     Error encountered in thermodynamic computations      ',/,  &
'@       (cfther.f90), for perfect gas with constant gamma.',/,   &
'@',/,                                                            &
'@     The computation of the subsonic inlet boundary',/,         &
'@     condition with imposed total pressure and',/,              &
'@     total enthalpy failed.',/,                                 &
'@',/,                                                            &
'@     At boundary face ',i10   , /,                              &
'@     The direction vector given by the user can''t be null.',/, &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

9020 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    in thermodynamics computations',/,              &
'@    =======',/,                                                 &
'@     Error encountered in thermodynamic computations      ',/,  &
'@       (cfther.f90), for perfect gas with constant gamma.',/,   &
'@',/,                                                            &
'@       in the computation of the subsonic inlet with',/,        &
'@       imposed total pressure and total enthalpy.',/,           &
'@',/,                                                            &
'@     At boundary face ',i10   , /,                              &
'@     The direction vector points outward the fluid domain.',/,  &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

end subroutine cf_thermo_ph_inlet_bc

!===============================================================================

!-------------------------------------------------------------------------------
!> \brief Compute epsilon sup:
!>
!> \f[\epsilon_{\textrm{sup}} = e - C_v T\f]
!>
!> for perfect gas: \f[\epsilon_{\textrm{sup}} = 0\f]

!> \param[out]    output  output work array
!> \param[in]     l_size  l_size of the array
!-------------------------------------------------------------------------------

subroutine cf_thermo_eps_sup(output, l_size)

!===============================================================================
! Module files
!===============================================================================

use ppincl

!===============================================================================

implicit none

! Arguments

integer l_size

double precision output(*)

! Local variables

integer ii

!===============================================================================

if (ieos.eq.1) then

  ! It is zero for a perfect gas

  do ii = 1, l_size
    output(ii) = 0.d0
  enddo

endif

end subroutine cf_thermo_eps_sup

!===============================================================================
