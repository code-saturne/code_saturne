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
! Function:
! ---------

!> \file uchida_model.f90
!>
!> \brief The COPAIN correlations used to approximate the condensation source
!> term and the thermal exchange coefficient to impose at the wall where
!> condensation occurs.
!>
!>
!> This subroutine is used to compute at each cell the
!> \f$ \Lambda _{\mbox{cond}} \f$ and \f$ h_{\mbox{f,bf}} \f$.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ifbpcd        index of faces with condensation source terms
!> \param[in]     izzftcd       faces zone with condensation source terms imposed
!>                              (at previous and current time steps)
!> \param[out]    gam_s         value associated to each variable in the
!>                              condensation source terms (Lambda_cond)
!> \param[out]    hpcond        value associated to the fluid exchange coeff.
!>                              evaluated with empiric law at the BC face with
!>                              condensation
!_______________________________________________________________________________

subroutine condensation_uchida_model &
 ( nvar   , nfbpcd , ifbpcd , izzftcd ,  &
   gam_s  , hpcond , regime)

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
use field
use mesh
use cs_c_bindings
use cs_f_interfaces
use cs_nz_condensation, only: thermal_condensation_flux, iztag1d, flthr, dflthr

use condensation_module

!===============================================================================

implicit none

! Arguments

integer          nvar, nfbpcd, ifbpcd(nfbpcd), izzftcd(nfbpcd)
integer          regime

double precision gam_s(nfbpcd,nvar)
double precision hpcond(nfbpcd)

! Local variables

integer          ii, iz, iel, ifac
integer          ivar, f_id, ifcvsl

double precision flux, h_uchida, h_cond
double precision sink_term, gamma_cond, pressure
double precision lambda, theta, tinf, psat, t_wall
double precision mol_mas_int
double precision Prdtl, xnu, drho, gravity
double precision rho_wall, rho_av
double precision x_vapint,x_ncond_int, x_inc
double precision h1min, h1max
double precision h2min, h2max
double precision flmin, flmax

type(gas_mix_species_prop) s_h2o_g
type(var_cal_opt) :: vcopt

double precision, allocatable, dimension(:) :: mix_mol_mas, mol_mas_ncond
double precision, allocatable, dimension(:) :: x_h2o_g, diff_m
double precision, dimension(:), pointer :: cpro_rho, cpro_viscl, cpro_cp, cpro_venth
double precision, dimension(:), pointer :: cvar_enth
double precision, dimension(:), pointer :: y_h2o_g

if (regime > 0) then
  write(nfecra,*) "***********************************************************************"
  write(nfecra,*) "*** Stop : icondb_regime > 0 non authorized with (icondb_model = 2) ***"
  write(nfecra,*) "***********************************************************************"
  call csexit(1)
endif

!===============================================================================
! Allocate a temporary array for cells selection
allocate(mix_mol_mas(ncelet), diff_m(ncelet), mol_mas_ncond(ncelet) )
allocate(x_h2o_g(ncelet))

!===============================================================================
! 0 - Initialization
!===============================================================================

h1min = 1.d+20 ; h1max = -1.d+20
h2min = 1.d+20 ; h2max = -1.d+20
flmin = 1.d+20 ; flmax = -1.d+20

sink_term = 0.d0

call field_get_id("y_h2o_g", f_id)
call field_get_key_struct_gas_mix_species_prop(f_id, s_h2o_g)

!===============================================================================
! 2 - Pointers to phys. properties for mixing flow at cell center.
!===============================================================================

! Coefficient exchange of the enthalpy scalar
ivar = isca(iscalt)

call field_get_val_s(ivarfl(ivar), cvar_enth)

! Condensable gas H20 (vapor)
call field_get_val_s_by_name("y_h2o_g", y_h2o_g)

!-- Pointer to density value
call field_get_val_s(icrom, cpro_rho)
! --- Molecular dynamic viscosity value
call field_get_val_s(iviscl, cpro_viscl)

!-- Specific heat
if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
else
  write(nfecra,1000) icp
  call csexit (1)
endif

!-- (Lambda/Cp) of the thermal scalar
call field_get_key_int (ivarfl(isca(iscalt)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_venth)
else
  write(nfecra,1010) iscalt
  call csexit (1)
endif

! General properties
gravity = sqrt(gx**2+gy**2+gz**2)
if (idilat.eq.3) then
  pressure = pther
else
  pressure = p0
endif

!===============================================================================
! 3 - Compute mass fraction of H2O, etc.
!===============================================================================

! mix_mol_mas = molecular weight of the mixture
! mol_mas_ncond = molecular weight of the non condensable gas
! x_ncond = mole fraction of non condensable gas
! x_h2o_g = mole fraction of steam
! diff_m = "binary" diffusivity of steam into the non condensable gas
call compute_mix_properties(ncel, mix_mol_mas, mol_mas_ncond, x_h2o_g, diff_m)

!===============================================================================
! 4 - Loop on the (ii) boundary faces with condensation
!     -------------------------------------------------
!     -> the sink source term gam_s(ii) associated to the condensation
!        is computed at each face associated to the cold wall
!        (loss of vapor mass flow rate due to condensation)
!     -> the exchange coefficient hpcond(ii) associated to the phase change
!        is computed at each face and imposed at the boundary condition
!        for the enthalpy scalar.
!        --------------------------------------------------------------
!          the correlations used are those from the copain modelling
!===============================================================================

do ii = 1, nfbpcd

  ifac= ifbpcd(ii)
  iel = ifabor(ifac)

  iz  = izzftcd(ii)


  xnu = cpro_viscl(iel)/cpro_rho(iel)

  call get_wall_temperature(iz, ii, t_wall)
  call get_temperature(cvar_enth(iel), cpro_cp(iel), tinf)

  call compute_psat(t_wall, psat)
  x_vapint = psat/pressure

  call compute_prandtl(cpro_viscl(iel), cpro_venth(iel), Prdtl)
  lambda = cpro_venth(iel)*cpro_cp(iel)

  x_ncond_int = 1.d0 - x_vapint
  mol_mas_int = x_vapint*s_h2o_g%mol_mas + x_ncond_int*mol_mas_ncond(iel)
  rho_wall = pressure*mol_mas_int / (cs_physical_constants_r*t_wall)
  rho_av = 0.5d0 * (rho_wall + cpro_rho(iel))
  drho = dabs(rho_wall - cpro_rho(iel)) / rho_av
  theta = 1.d0

  ! Condensation case
  if (x_h2o_g(iel).gt.x_vapint) then
    x_inc = 1.d0-x_h2o_g(iel)

! Total heat flux (latent + sensible) from 2*UCHIDA correlation
    h_uchida = 380.0d0 * ( (1.0d0-y_h2o_g(iel)) / y_h2o_g(iel) )**(-0.7d0)

! Sensible component (Mc Adams)
    theta = 1.d0 +0.625d0*(x_ncond_int-x_inc)/x_ncond_int
    hpcond(ii)= 0.13d0*theta*lambda                     &
           *( gravity*drho*Prdtl/(xnu**2) )**(1.d0/3.d0) &
           /cpro_cp(iel)

! Latent component deduced from (total - sensible)
    h_cond = h_uchida - hpcond(ii) * cpro_cp(iel)

    ! Computation of sink source term gam_s(ii)
    !         per unit of area (kg/s/m2)
    sink_term = h_cond * (tinf - t_wall) / lcond
    gam_s(ii,ipr) = gam_s(ii, ipr) - sink_term

    ! Heat flux
    flux = h_uchida * (tinf - t_wall)

  else ! no condensation

    theta = 1.d0
    hpcond(ii) = 0.13d0*theta*lambda                     &
           *( gravity*drho*Prdtl/(xnu**2) )**(1.d0/3.d0) &
           /cpro_cp(iel)
    flux = hpcond(ii)*cpro_cp(iel)*(tinf-t_wall)
    sink_term = 0.d0
    h_cond = 0.0d0
  endif

  ! Store it for postprocessing
  thermal_condensation_flux(ii) = flux

  h1min = min(h1min, h_cond)
  h1max = max(h1max, h_cond)
  h2min = min(h2min, hpcond(ii)*cpro_cp(iel))
  h2max = max(h2max, hpcond(ii)*cpro_cp(iel))
  flmin = min(flmin,flux)
  flmax = max(flmax,flux)

  !For 1D thermal conduction model :
  !store the flux and its derivative
  if(iztag1d(iz).eq.1) then
    flthr(ii) = flux
   dflthr(ii) = 0.d0
  endif

enddo !  End of the Loop on the (ii) boundary faces with condensation

! Deallocate the temporary array
deallocate(mix_mol_mas, mol_mas_ncond)
deallocate(x_h2o_g, diff_m)

!===============================================================================
! 5 - Print in the log every ntlist of the min/max values
!         for exchange coefficient and flux computed
!===============================================================================

gamma_cond = 0.d0
do ii = 1, nfbpcd
  ifac = ifbpcd(ii)
  gamma_cond = gamma_cond + gam_s(ii, ipr)
enddo

if (irangp.ge.0) then
  call parmin(h1min)
  call parmin(h2min)
  call parmax(h1max)
  call parmax(h2max)
  call parmin(flmin )
  call parmax(flmax )
  call parsom(gamma_cond)
endif

if (mod(ntcabs,ntlist).eq.0) then
  write(nfecra,*) ' Minmax values of latent htc: ', &
                    ntcabs, h1min, h1max
  write(nfecra,*) ' Minmax values of sensible htc: ', &
                    ntcabs, h2min, h2max
  write(nfecra,*) ' Minmax values of thermal flux        : ', &
                    ntcabs, flmin, flmax
  write(nfecra,1061) gamma_cond
endif

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)
 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@      usipsu indique que la chaleur specifique est uniforme ',/,&
'@        ICP = ',I10   ,' alors que                          ',/,&
'@      copain model impose une chaleur specifique variable.  ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipsu ou copain model.                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@    La diffusivite du scalaire ',i10,' est uniforme alors   ',/,&
'@      que copain model impose une diffusivite variable.     ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipsu ou copain model.                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
#else
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
'@      copain model prescribes a variable specific heat.',/,     &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu',/,                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@    For scalar', i10,/,                                         &
'@      usipsu specifies that the diffusivity is uniform',/,      &
'@        ivislc(',i10   ,') = ',i10   ,' while',/,               &
'@      copain model prescribes a variable diffusivity.',/,       &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or copain model',/,                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
#endif
 1061 format(/,                                                   &
' ** Condensation source terms (Gamma) added:',E14.5           ,/,&
'    ----------------------------------------',/)

!----
! End
!----

return
end subroutine condensation_uchida_model
