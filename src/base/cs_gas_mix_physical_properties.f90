!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file cs_gas_mix_physical_properties.f90
!>
!> \brief This subroutine fills physical properties which are variable in time
!>        for the gas mixtures modelling with or without steam inside the fluid
!>        domain. In presence of steam, this one is deduced from the
!>        noncondensable gases transported as scalars
!>        (by means of the mass fraction of each species).
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________


subroutine cs_gas_mix_physical_properties

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
use cs_f_interfaces

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iel   , iscal, ifcvsl, iesp, jesp, ierror, f_id

character(len=80) :: name_i, name_j, name_d

double precision xsum_mu, xsum_lambda, phi_mu, phi_lambda, x_k
double precision mu_i, mu_j, lambda_i, lambda_j

type(gas_mix_species_prop), pointer :: s_j, s_i
type(gas_mix_species_prop), target :: s_d
type(gas_mix_species_prop), dimension(:), allocatable, target :: s_k

double precision, allocatable, dimension(:), target :: lam_loc
double precision, allocatable, dimension(:), target :: tk_loc

double precision, dimension(:), pointer :: cpro_rho
double precision, dimension(:), pointer :: cpro_viscl, cpro_cp
double precision, dimension(:), pointer :: cpro_venth, cpro_vyk
double precision, dimension(:), pointer :: cvar_enth , cvar_yk, tempk
double precision, dimension(:), pointer :: cvar_yi, cvar_yj
double precision, dimension(:), pointer :: y_d, ya_d
double precision, dimension(:), pointer :: mix_mol_mas
double precision, dimension(:), pointer :: lambda

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

ierror = 0
if (irovar.le.0) ierror = 1
if (ivivar.le.0) ierror = 1
if (icp.lt.0)    ierror = 1

if (ierror.gt.0) then
  call csexit(1)
endif

! In compressible, the density is updated after the pressure step (cfmspr)
if (ippmod(icompf).lt.0) then
  call field_get_val_s(ivarfl(isca(iscalt)), cvar_enth)
  ! Density value
  call field_get_val_s(icrom, cpro_rho)
  allocate(tk_loc(ncel))
  tempk => tk_loc
else
  call field_get_val_s(ivarfl(isca(itempk)), tempk)
endif

! Molecular dynamic viscosity value
call field_get_val_s(iviscl, cpro_viscl)
! Specific heat value
if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)

! Stop if Cp is not variable
else
  write(nfecra,1000) icp
  call csexit (1)
endif

! Lambda/Cp value
call field_get_key_int(ivarfl(isca(iscalt)), kivisl, ifcvsl)
call field_get_val_s(ifcvsl, cpro_venth)

call field_get_val_s(igmxml, mix_mol_mas)

! Deduce mass fraction (y_d) which is
! y_h2o_g in presence of steam or
! y_he/y_h2 with noncondensable gases
if (ippmod(igmix).eq.0) then
  name_d = "y_he"
elseif (ippmod(igmix).eq.1) then
  name_d = "y_h2"
elseif (ippmod(igmix).ge.2.and.ippmod(igmix).lt.5) then
  name_d = "y_h2o_g"
else ! ippmod(igmix).eq.5
  name_d = "y_o2"
endif
call field_get_val_s_by_name(name_d, y_d)
call field_get_val_prev_s_by_name(name_d, ya_d)
call field_get_id(name_d, f_id)
call field_get_key_struct_gas_mix_species_prop(f_id, s_d)

if (ippmod(icompf).lt.0) then
  allocate(lam_loc(ncelet))
  lambda => lam_loc
else
  call field_get_key_int(ivarfl(isca(itempk)), kivisl, ifcvsl)
  call field_get_val_s(ifcvsl, lambda)
endif

allocate(s_k(nscasp+1))

!===============================================================================
!2. Define the physical properties for the gas mixture with:
!   - the density (rho_m) and specific heat (cp_m) of the gas mixture function
!     temperature and species scalar (yk),
!   - the dynamic viscosity (mu_m) and conductivity coefficient (lbd_m) of
!     the gas mixture function ot the enthalpy and species scalars,
!   - the diffusivity coefficients of the scalars (Dk, D_enh).
!===============================================================================

!Storage the previous value of the deduced mass fraction ya_d
do iel=1, ncelet
  ya_d(iel) = y_d(iel)
enddo

!-----------------------------------------
! Compute the mass fraction (y_d) deduced
! from the mass fraction (yk) transported
!-----------------------------------------

! Initialization
do iel = 1, ncel
  y_d(iel) = 1.d0

  ! Mixture specific heat
  cpro_cp(iel) = 0.d0

  ! Mixture molar mass
  mix_mol_mas(iel) = 0.d0

  ! Mixture molecular diffusivity
  cpro_viscl(iel) = 0.d0

  ! Thermal conductivity
  lambda(iel) = 0.d0
enddo

do iesp = 1, nscasp
  ! Mass fraction array of the different species
  call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)

  call field_get_key_struct_gas_mix_species_prop( &
       ivarfl(isca(iscasp(iesp))), s_k(iesp))

  do iel = 1, ncel
    y_d(iel) = y_d(iel)-cvar_yk(iel)
    mix_mol_mas(iel) = mix_mol_mas(iel) + cvar_yk(iel)/s_k(iesp)%mol_mas
  enddo
enddo

! Clipping
do iel = 1, ncel
  y_d(iel) = max(y_d(iel), 0.d0)
enddo

!Finalize the computation of the Mixture molar mass
s_k(nscasp+1) = s_d
do iel = 1, ncel
  mix_mol_mas(iel) = mix_mol_mas(iel) + y_d(iel)/s_d%mol_mas
  mix_mol_mas(iel) = 1.d0/mix_mol_mas(iel)
enddo

!==============================================================
! Mixture specific heat function of species specific heat (cpk)
! and mass fraction of each gas species (yk), as below:
!             -----------------------------
! - noncondensable gases and the mass fraction deduced:
!             cp_m(iel) = Sum( yk.cpk)_k[0, nscasp]
!                       + y_d.cp_d
!             -----------------------------
! remark: the mass fraction deduced depending of the
!         modelling chosen by the user.
!         with:
!             - igmix = 0 or 1, a noncondensable gas
!             - igmix > 2     , a condensable gas (steam)
!==============================================================
do iesp = 1, nscasp
  call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)

  do iel = 1, ncel
    cpro_cp(iel) = cpro_cp(iel) + cvar_yk(iel)*s_k(iesp)%cp
  enddo
enddo

! Finalization
do iel = 1, ncel
  cpro_cp(iel) = cpro_cp(iel) + y_d(iel)*s_d%cp
enddo

!===========================================================
! gas mixture density function of the temperature, pressure
! and the species scalars with taking into account the
! dilatable effects, as below:
!             ----------------------
! - with inlet/outlet conditions:
!   [idilat=2]: rho= p0/(R. temp(1/sum[y_i/M_i]))
! - with only wall conditions:
!   [idilat=3]: rho= pther/(R. temp(1/sum[y_i/M_i]))
!             ----------------------
!         i ={1, .. ,N} : species scalar number
!         y_i           : mass fraction of each species
!         M_i           : molar fraction [kg/mole]
!         R             : ideal gas constant [J/mole/K]
!         p0            : atmos. pressure (Pa)
!         pther         : pressure (Pa) integrated on the
!                         fluid domain
!===========================================================

if (ippmod(icompf).lt.0) then
  do iel = 1, ncel
    ! Evaluate the temperature thanks to the enthalpy
    tempk(iel) = cvar_enth(iel)/ cpro_cp(iel)
    if (idilat.eq.3) then
      cpro_rho(iel) = pther*mix_mol_mas(iel)/(cs_physical_constants_r*tempk(iel))
    else
      cpro_rho(iel) = p0*mix_mol_mas(iel)/(cs_physical_constants_r*tempk(iel))
    endif
  enddo
endif

!==================================================
! Dynamic viscosity and conductivity coefficient
! the physical properties associated to the gas
! mixture with or without condensable gas.
!==================================================

! Loop over ALL the species
do iesp = 1, nscasp+1

  s_i => s_k(iesp)

  ! Mass fraction deduced
  ! (as steam or noncondensable gas)
  if (iesp.eq.nscasp+1) then
    cvar_yi => y_d
    name_i = name_d
  ! Noncondensable species
  else
    call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yi)
    call field_get_name (ivarfl(isca(iscasp(iesp))), name_i)
  endif

  do iel = 1, ncel

    ! Viscosity and conductivity laws
    ! for each mass fraction species

    if (ivsuth.eq.0) then
        ! With a linear law
    call cs_local_physical_properties &
    !================================
   ( mu_i, lambda_i, tempk(iel), tkelvi, s_i, name_i)
    else
      ! Or : with a Sutherland law
      call cs_local_physical_properties_suth                        &
      !================================
     ( mu_i, lambda_i, tempk(iel), s_i, name_i)
    endif

    xsum_mu = 0.d0
    xsum_lambda = 0.d0

    ! Loop over ALL the species
    do jesp = 1, nscasp+1

      s_j => s_k(jesp)

      ! Mass fraction deduced
      ! (as steam or noncondensable gas)
      if (jesp.eq.nscasp+1) then
        cvar_yj => y_d
        name_j = name_d
      ! Noncondensable species
      else
        call field_get_val_s(ivarfl(isca(iscasp(jesp))), cvar_yj)
        call field_get_name (ivarfl(isca(iscasp(jesp))), name_j)
      endif

      if (ivsuth.eq.0) then
          ! With a linear law
      call cs_local_physical_properties &
      !================================
     ( mu_j, lambda_j, tempk(iel), tkelvi, s_j, name_j)
      else
        ! Or : with a Sutherland law
        call cs_local_physical_properties_suth                        &
        !================================
       ( mu_j, lambda_j, tempk(iel), s_j, name_j)
      endif

      phi_mu = (1.d0/sqrt(8.d0))                              &
           *(1.d0 +  s_i%mol_mas / s_j%mol_mas)**(-0.5d0)     &
           *(1.d0 + (mu_i        / mu_j       )**(+0.5d0)     &
           * (s_j%mol_mas / s_i%mol_mas)**(+0.25d0))**2


      phi_lambda = (1.d0/sqrt(8.d0))                          &
           *(1.d0 +  s_i%mol_mas / s_j%mol_mas)**(-0.5d0)     &
           *(1.d0 + (lambda_i    / lambda_j   )**(+0.5d0)     &
           * (s_j%mol_mas / s_i%mol_mas)**(+0.25d0))**2

      x_k = cvar_yj(iel)*mix_mol_mas(iel)/s_j%mol_mas
      xsum_mu = xsum_mu + x_k * phi_mu

      xsum_lambda = xsum_lambda + x_k * phi_lambda

    enddo

    ! Mixture viscosity defined as function of the scalars
    !-----------------------------------------------------
    x_k = cvar_yi(iel)*mix_mol_mas(iel)/s_i%mol_mas
    cpro_viscl(iel) = cpro_viscl(iel) + x_k * mu_i / xsum_mu

    lambda(iel) = lambda(iel) + x_k * lambda_i / xsum_lambda

  enddo
enddo

!=====================================================
! Dynamic viscosity and conductivity coefficient
! the physical properties filled for the gas mixture
!=====================================================
! Same diffusivity for all the scalars except the enthalpy
do iesp = 1, nscasp

  iscal = iscasp(iesp)

  call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
  call field_get_val_s(ifcvsl, cpro_vyk)

  do iel = 1, ncel
    cpro_vyk(iel)= cpro_viscl(iel)
  enddo

enddo

if(ippmod(icompf).lt.0) then
  ! --- Lambda/Cp of the thermal scalar
  do iel = 1, ncel
    cpro_venth(iel) = lambda(iel)/cpro_cp(iel)
  enddo

  ! deallocate local arrays if not compressible
  deallocate(tk_loc)
  deallocate(lam_loc)
  tempk => null()
  lambda => null()
endif

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
'@      usipsu specifies that the specific heat is uniform',/,    &
'@        icp = ',i10   ,' while',/,                              &
'@      usphyv prescribes a variable specific heat.',/,           &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or usphyv.',/,                                &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

!----
! End
!----

return
end subroutine cs_gas_mix_physical_properties

!===============================================================================
! Purpose:
! -------

!> \file cs_gas_mix_physical_properties.f90
!>
!> \brief This user subroutine is used to compute the dynamic viscosity and
!> conductivity coefficient associated to each gas species.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    mu            dynamic viscosity associated to the gas species
!> \param[out]    lambda        conductivity coefficient of the gas species
!> \param[in]     tk            temperature variable in kelvin
!> \param[in]     tkelvin       reference temperature value
!> \param[in]     spro          constants used for the physcial laws
!> \param[in]     name          name of the field associated to the gas species
!_______________________________________________________________________________

subroutine cs_local_physical_properties(mu, lambda, tk, tkelvin, spro, name)
!===============================================================================

use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

double precision mu, lambda
double precision tk, tkelvin

character(len=80) :: name

type(gas_mix_species_prop) spro

!===============================================================================
! The viscosity law for each species is defined
! as below:
!              ----------------------------------
! 1/. for Steam species:
!      mu = mu_a.(tk - tkelvi), with a law in (°C) unit
! 2/. for Helium species:
!      mu = mu_a .(tk/tkelvi)**c + mu_b, with nodim. law
! 3/. for Hydrogen species:
!      mu = mu_a .(tk-tkelvi) + mu_b, with t (°C)
! 4/. for Oxygen and Nitrogen species:
!      mu = mu_a .(tk) + mu_b, with t in (°K)
!
! The conductivity expression for each species is
! defined as:
!              ----------------------------------
! 1/. for Steam species:
!      lambda = lambda_a .(tk-tkelvi) + lambda_b, with a law in (°C) unit
! 2/. for Helium species:
!     lambda = lambda_a .(tk/tkelvi)**c + lambda_b, with nodim. law
! 3/. for Hydrogen species:
!      lambda = lambda_a .tk + lambda_b, with tk (°K)
! 4/. for Oxygen and Nitrogen species :
!      lambda = lambda_a .(tk) + lambda_b, with t (°K)
!===============================================================================

if (name.eq.'y_h2o_g') then
  mu     = spro%mu_a    *(tk-tkelvin) + spro%mu_b
  lambda = spro%lambda_a*(tk-tkelvin) + spro%lambda_b
elseif(name.eq.'y_he') then
  mu     = spro%mu_a     * (tk/tkelvin)**0.7d0
  lambda = spro%lambda_a * (tk/tkelvin)**0.7d0
elseif(name.eq.'y_h2') then
  mu     = spro%mu_a     * (tk-tkelvin) + spro%mu_b
  lambda = spro%lambda_a * tk + spro%lambda_b
elseif (name.eq.'y_o2'.or.name.eq.'y_n2') then
  mu     = spro%mu_a     * tk + spro%mu_b
  lambda = spro%lambda_a * tk + spro%lambda_b
else
  call csexit(1)
endif

!----
! End
!----
return

end subroutine cs_local_physical_properties

subroutine cs_local_physical_properties_suth(mu, lambda, tk,spro,name)
!===============================================================================
use field
use cs_c_bindings
use cstphy
use ppthch
!===============================================================================

implicit none

! Arguments

double precision mu, lambda
double precision tk

character(len=80) :: name
type(gas_mix_species_prop)  spro

! Local variables

double precision muref, lamref
double precision trefmu, treflam, smu, slam

!===============================================================================
! Sutherland law for viscosity and thermal conductivity
! The viscosity law for each specie is defined
! as below:
!              ----------------------------------
! mu = muref*(T/Tref)**(3/2)*(Tref+S1)/(T+S1)

! The conductivity expression for each specie is
! defined as:
!              ----------------------------------
!  lambda = lambdaref*(T/Tref)**(3/2)*(Tref+S2)/(T+S2)
!             ------------------------------------
! S1 and S2 are respectively Sutherland temperature for conductivity and
! Sutherland temperature for viscosity of the considered specie
! Tref is a reference temperature, equal to 273K for a perfect gas.
! For steam (H20), Tref has not the same value in the two formulae.
! Available species : O2, N2, H2, H20 and  He
! The values for the parameters come from F.M. White's book "Viscous Fluid Flow"
!================================================================================

if (name.ne.'y_h2o_g' .and. name.ne.'y_he' .and. name.ne.'y_o2'         &
     .and. name.ne.'y_n2' .and. name.ne.'y_h2') then
  call csexit(1)
endif

muref = spro%muref
lamref = spro%lamref
trefmu = spro%trefmu
treflam = spro%treflam
smu = spro%smu
slam = spro%slam

mu =  muref * (tk / trefmu)**1.5d0                              &
     * ((trefmu+smu) / (tk+smu))
lambda = lamref  * (tk / treflam)**1.5d0                        &
     * ((treflam+slam) / (tk+slam))

!----
! End
!----
return

end subroutine cs_local_physical_properties_suth
