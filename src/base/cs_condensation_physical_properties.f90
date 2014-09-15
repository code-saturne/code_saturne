!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

!> \file cs_condensation_physical_properties.f90
!>
!> \brief This subroutine fills physical properties which are variable in time
!>        for the condensation modelling.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________


subroutine cs_condensation_physical_properties

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use pointe
use albase
use parall
use period
use ihmpre
use ppppar
use ppthch
use ppincl
use mesh
use field
use cavitation
use cs_c_bindings
use dimens, only: nscal
use cs_f_interfaces

!===============================================================================

implicit none

! Arguments

! Local variables

character*80     chaine

integer          iel   , iscal, iesp, jesp, ierror
integer          f_id  , ii

double precision xsum_mu, xsum_lambda, phi_mu, phi_lambda, tk, x_k
double precision mu_i, mu_j, lambda_i, lambda_j

type(severe_acc_species_prop), pointer :: s_j, s_i
type(severe_acc_species_prop), target :: s_h2o_g
type(severe_acc_species_prop), dimension(:), allocatable, target :: s_k

double precision, allocatable, dimension(:) :: lambd_m, mix_mol_mas

double precision, dimension(:), pointer :: cpro_rho
double precision, dimension(:), pointer :: cpro_viscl, cpro_cp
double precision, dimension(:), pointer :: cpro_venth, cpro_vyk
double precision, dimension(:), pointer :: cvar_enth , cvar_yk
double precision, dimension(:), pointer :: cvar_yi, cvar_yj
double precision, dimension(:), pointer :: y_h2o_g


!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

ierror = 0
if (irovar.le.0) ierror = 1
if (ivivar.le.0) ierror = 1
if (icp.le.0)    ierror = 1

if (ierror.gt.0) then
  call csexit(1)
endif


call field_get_val_s(ivarfl(isca(iscalt)), cvar_enth)

! --- Density value
call field_get_val_s(icrom, cpro_rho)
! --- Molecular dynamic viscosity value
call field_get_val_s(iprpfl(iviscl), cpro_viscl)
! --- Specific heat
if (icp.gt.0) then
  call field_get_val_s(iprpfl(icp), cpro_cp)

! --- Stop if Cp is not variable
else
  write(nfecra,1000) icp
  call csexit (1)
endif

! --- Lambda/CP
call field_get_val_s(iprpfl(ivisls(iscalt)), cpro_venth)

! Condensable gas H20 (vapor)
call field_get_val_s_by_name("y_h2o_g", y_h2o_g)
call field_get_id("y_h2o_g", f_id)
call field_get_key_struct_severe_acc_species_prop(f_id, s_h2o_g)

allocate(lambd_m(ncelet), mix_mol_mas(ncelet))

allocate(s_k(nscasp+1))
!===============================================================================
!2. Define the physical properties for the mixinf flow with:
!   - the density variable and function of the temperature and species scalars
!   - the specific heat and dynamic viscosity function of the species scalars
!   - the conductivity and diffusivity coefficients of the scalars are defined
!     as below.
!===============================================================================

! Compute the mass fraction of H2O vapor
! Deduced from the non-condensable species (O2, N2)
!--------------------------------------------------

! Initialization
do iel = 1, ncel
  y_h2o_g(iel) = 1.d0

  ! Mixing Specific heat function
  cpro_cp(iel) = 0.d0

  ! Mixing molar mass
  mix_mol_mas(iel) = 0.d0

  ! Mixing molecular diffusivity
  cpro_viscl(iel) = 0.d0

  lambd_m(iel) = 0.d0
enddo

do iesp = 1, nscasp
  ! mass fraction array of the different
  ! species (O2, N2)
  !-------------------------------------------------
  call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)

  call field_get_key_struct_severe_acc_species_prop( &
       ivarfl(isca(iscasp(iesp))), s_k(iesp))

  do iel = 1, ncel
    y_h2o_g(iel) = y_h2o_g(iel)-cvar_yk(iel)
    mix_mol_mas(iel) = mix_mol_mas(iel) + cvar_yk(iel)/s_k(iesp)%mol_mas
  enddo
enddo

! Clipping
do iel = 1, ncel
  y_h2o_g(iel) = max(y_h2o_g(iel), 0.d0)
enddo

!Finalize the computation of the Mixing molar mass
s_k(nscasp+1) = s_h2o_g
do iel = 1, ncel
  mix_mol_mas(iel) = mix_mol_mas(iel) + y_h2o_g(iel)/s_h2o_g%mol_mas
  mix_mol_mas(iel) = 1.d0/mix_mol_mas(iel)
enddo

! Mixing Specific heat function
!------------------------------
do iesp = 1, nscasp
  call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)

  do iel = 1, ncel
    cpro_cp(iel) = cpro_cp(iel) + cvar_yk(iel)*s_k(iesp)%cp
  enddo
enddo

! Finalization
do iel = 1, ncel
  cpro_cp(iel) = cpro_cp(iel) + y_h2o_g(iel)*s_h2o_g%cp
enddo

!==================================================
! Mixing density function of the constant pressure
! and the species and temperature scalars as below:
!             ----------------------
!        rho = p0/( R T (1/sum[Yi/Mi]) )
!             ----------------------
! with:
!         i ={1, .. ,N} species
!         Yi : mass fraction of each species
!         Mi : molar fraction [kg/mole]
!         R  = prefect gas constant [J/mole/K]
!         p0 = atmos. pressure (Pa)
!==================================================

do iel = 1, ncel
  ! Evaluate the temperature thanks to the enthalpy
  tk = cvar_enth(iel)/ cpro_cp(iel)
  cpro_rho(iel) = p0*mix_mol_mas(iel)/(rr*tk)

enddo

!==================================================
! Bulk viscosity computation
! and
! Bulk enthalpy Conductivity
! computation
!==================================================

! Loop over ALL the species
do iesp = 1, nscasp+1

  s_i => s_k(iesp)

  ! H2O_g
  if (iesp.eq.nscasp+1) then
    cvar_yi => y_h2o_g
  ! Non condensable species
  else
    call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yi)
  endif

  do iel = 1, ncel

    tk = cvar_enth(iel) / cpro_cp(iel)

    ! The viscosity law for each species is defined
    ! as below:
    ! 1/. for O2 and N2 species:
    !      mu = mu_a T + mu_b, with T (°K)
    ! 2/. for H20_v species:
    !      mu = mu_a (T - tkelvi), with dT (°C)
    !              ------------------------
    ! The conductivity expression for each species is
    ! defined as:
    ! 1/. for O2 and N2 species :
    !      lambda = lambda_a . T + lambda_b, with T (°K)
    ! 2/. for H20g species:
    !      lambda = lambda_a . (T-tkelvi) + lambda_b, with dT (°C)
    if (iesp.eq.nscasp+1) then
      mu_i     = s_i%mu_a    *(tk-tkelvi) + s_i%mu_b
      lambda_i = s_i%lambda_a*(tk-tkelvi) + s_i%lambda_b
    else
      mu_i     = s_i%mu_a     * tk + s_i%mu_b
      lambda_i = s_i%lambda_a * tk + s_i%lambda_b
    endif

    xsum_mu = 0.d0
    xsum_lambda = 0.d0

    ! Loop over ALL the species
    do jesp = 1, nscasp+1

      s_j => s_k(jesp)

      ! H2O_g
      if (jesp.eq.nscasp+1) then
        cvar_yj => y_h2o_g
        mu_j     = s_j%mu_a    *(tk-tkelvi) + s_j%mu_b
        lambda_j = s_j%lambda_a*(tk-tkelvi) + s_j%lambda_b

      ! Non condensable species
      else
        call field_get_val_s(ivarfl(isca(iscasp(jesp))), cvar_yj)
        mu_j     = s_j%mu_a     * tk + s_j%mu_b
        lambda_j = s_j%lambda_a * tk + s_j%lambda_b
      endif

      phi_mu = (1.d0/sqrt(8.d0))                                 &
          *(1.d0 +  s_i%mol_mas / s_j%mol_mas)**(-0.5d0)      &
          *(1.d0 + (mu_i        / mu_j       )**(+0.5d0)      &
                 * (s_j%mol_mas / s_i%mol_mas)**(+0.25d0))**2


      phi_lambda = (1.d0/sqrt(8.d0))                                 &
          *(1.d0 +  s_i%mol_mas / s_j%mol_mas)**(-0.5d0)      &
          *(1.d0 + (lambda_i    / lambda_j   )**(+0.5d0)      &
                 * (s_j%mol_mas / s_i%mol_mas)**(+0.25d0))**2

      x_k = cvar_yj(iel)*mix_mol_mas(iel)/s_j%mol_mas
      xsum_mu = xsum_mu + x_k * phi_mu
      xsum_lambda = xsum_lambda + x_k * phi_lambda

    enddo

    ! Mixing viscosity defined as function of the scalars
    !----------------------------------------------------
    x_k = cvar_yi(iel)*mix_mol_mas(iel)/s_i%mol_mas
    cpro_viscl(iel) = cpro_viscl(iel) + x_k * mu_i / xsum_mu
    lambd_m(iel) = lambd_m(iel) + x_k * lambda_i / xsum_lambda

  enddo
enddo

! Same diffusivity for all the scalars except the enthalpy
do ii = 1, nscapp

  iscal = iscapp(ii)

  if (iscal.ne.iscalt) then

    call field_get_val_s(iprpfl(ivisls(iscal)), cpro_vyk)

    do iel = 1, ncel
      cpro_vyk(iel)= cpro_viscl(iel)
    enddo

  endif

enddo

! --- Lambda/Cp of the thermal scalar
do iel = 1, ncel
  cpro_venth(iel) = lambd_m(iel)/cpro_cp(iel)
enddo

deallocate(lambd_m)
deallocate(s_k)

!===============================================================================
! 3. Checking of the user values
!===============================================================================

!--------
! Formats
!--------

 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@      usipph specifies that the specific heat is uniform',/,    &
'@        icp = ',i10   ,' while',/,                              &
'@      usphyv prescribes a variable specific heat.',/,           &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipph or usphyv.',/,                                &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

!----
! End
!----

return
end subroutine
