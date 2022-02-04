!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
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

!> \file copain_model.f90
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

subroutine condensation_copain_model &
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
use cs_nz_condensation, only: thermal_condensation_flux, izcophc, izcophg, &
                              iztag1d, zxrefcond, zprojcond, flthr, dflthr

use condensation_module


!===============================================================================

implicit none

! Arguments

integer          nvar, nfbpcd, ifbpcd(nfbpcd), izzftcd(nfbpcd)
integer          regime

double precision gam_s(nfbpcd,nvar)
double precision hpcond(nfbpcd)

! Local variables

integer          ii, iz, iel, ifac, iesp
integer          ivar, f_id, ifcvsl
integer          ustar_id, conv_regime

double precision flux, pressure
double precision sink_term, gamma_cond
double precision lambda, theta, tinf, psat
double precision Sh_z, Nu_z, Gr_z, schdt, Prdtl, Re_z
double precision xnu
double precision x_inc,y_ncond
double precision x_vapint,x_ncond_int,y_ncond_int
double precision distbf
double precision drho
double precision mix_mol_mas_int
double precision gravity
double precision dplus, yplus, sigmat, ypth
double precision hcond,hcdt,hcdcop,hw_cop,hflui,hpflui,hw_enth
double precision h1min, h1max, h2min, h2max
double precision h3min, h3max, h4min, h4max
double precision flmin, flmax
double precision t_wall
double precision uk, rough_t
double precision lcar, u_ref

type(gas_mix_species_prop) s_h2o_g
type(var_cal_opt) :: vcopt

double precision, allocatable, dimension(:) :: mix_mol_mas, mol_mas_ncond
double precision, allocatable, dimension(:) :: x_h2o_g, diff_m
double precision, dimension(:), pointer :: cpro_rho, cpro_viscl, cpro_cp, cpro_venth
double precision, dimension(:), pointer :: cvar_enth, cvar_yk
double precision, dimension(:), pointer :: y_h2o_g
double precision, dimension(:), pointer :: yplbr
double precision, dimension(:), pointer :: bpro_ustar
double precision, dimension(:,:), pointer :: cvar_vel

!===============================================================================
! Allocate a temporary array for cells selection
allocate(mix_mol_mas(ncelet), mol_mas_ncond(ncelet))
allocate(x_h2o_g(ncelet), diff_m(ncelet))

!===============================================================================
! 0 - Initialization
!===============================================================================

h1min = 1.d+20 ; h1max = -1.d+20
h2min = 1.d+20 ; h2max = -1.d+20
h3min = 1.d+20 ; h3max = -1.d+20
h4min = 1.d+20 ; h4max = -1.d+20
flmin = 1.d+20 ; flmax = -1.d+20

sink_term = 0.d0

yplbr => null()

if (iyplbr.ge.0) call field_get_val_s(iyplbr, yplbr)

call field_get_id("y_h2o_g", f_id)
call field_get_key_struct_gas_mix_species_prop(f_id, s_h2o_g)

call field_get_id_try('ustar', ustar_id)
if (ustar_id.ge.0) then
  call field_get_val_s(ustar_id, bpro_ustar)
endif

! Convection regime (1: natural, 2: forced, 3: mixed)
conv_regime = regime + 1

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
!-- Stop if Cp is not variable
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

call field_get_val_v(ivarfl(iu), cvar_vel)       ! Velocity

if (idilat == 3) then
  pressure = pther
else
  pressure = p0
endif
gravity = sqrt(gx**2+gy**2+gz**2)

!===============================================================================
! 3 - Compute properties of the mix (molecular weight, mole fractions, diffusivity)
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

lcar = 1.0d0
do ii = 1, nfbpcd

  ifac= ifbpcd(ii)
  iel = ifabor(ifac)

  iz  = izzftcd(ii)
  distbf = distb(ifac) ! wall distance

  Re_z = 0.0d0
  if (conv_regime > 1) then
    call normalize_vector(3, zprojcond(1:3, iz))
    call compute_characteristic_length(cdgfbo(1:3, ifac), zxrefcond(1:3, iz), &
                                       zprojcond(1:3, iz), lcar)
    call compute_tangential_velocity(cvar_vel(1:3, iel), surfbo(1:3, ifac), 1.0d0/surfbn(ifac), u_ref)
    Re_z = cpro_rho(iel) * u_ref * lcar / cpro_viscl(iel)
  endif

  call get_wall_temperature(iz, ii, t_wall)

  !-- kinematic viscosity --------------------------
  xnu = cpro_viscl(iel)/cpro_rho(iel)

  ! Thermal roughness
  rough_t = 0.d0
  ! Turbulent friction velocity
  if (ustar_id.ge.0) then
    uk = bpro_ustar(ifac)
  else
    uk = 0.d0
  endif
  dplus = 0.d0
  yplus = yplbr(ifac)
  sigmat = 0.9d0

  ! Bulk temperature is taken at first cell center
  call get_temperature(cvar_enth(iel), cpro_cp(iel), tinf)
  call compute_drho(tinf, t_wall, drho)
  call compute_grashof(gravity, drho, lcar, xnu, Gr_z)

  call compute_psat(t_wall, psat)
  x_vapint = psat / pressure

  sink_term = 0.d0

  ! if (Xv > Xi,v) we have condensation  -----------
  if (x_h2o_g(iel).gt.x_vapint) then

    call compute_schmidt(xnu, diff_m(iel), schdt)
    ! Condensation exchange coefficient based on turbulent wall law
    call hturbp(iwalfs,xnu,schdt,sigmat,rough_t,uk,yplus,dplus,hflui,ypth)
    hcdt =  diff_m(iel)/distbf*hflui

    h3max = max(h3max,hcdt)
    h3min = min(h3min,hcdt)

    ! Condensation exchange coefficient based on COPAIN correlation
    x_inc = 1.d0-x_h2o_g(iel)
    x_ncond_int = 1.d0 - x_vapint
    mix_mol_mas_int = x_vapint*s_h2o_g%mol_mas + x_ncond_int*mol_mas_ncond(iel)
    y_ncond =  1.d0 - y_h2o_g(iel)
    y_ncond_int = x_ncond_int*mol_mas_ncond(iel)/mix_mol_mas_int

    ! Corrective factor to account for suction and film effects
    theta = 1.d0 +0.625d0*(x_ncond_int-x_inc)/x_ncond_int
    call compute_exchange_adimensional(theta, Re_z, Gr_z, schdt, conv_regime, Sh_z)
    call hcd_copain(diff_m(iel), Sh_z, lcar, hcdcop)

    h4max = max(h4max,hcdcop)
    h4min = min(h4min,hcdcop)

    ! Choose between turbulent wall law and COPAIN correlation
    if (izcophc(iz).eq.1) then
      hcond = hcdt
    else if (izcophc(iz).eq.2) then
      hcond = hcdcop
    else if (izcophc(iz).eq.3) then
      hcond = max(hcdt,hcdcop)
    endif

    ! Final value of condensation exchange coefficient
    hcond = cpro_rho(iel) * hcond * (y_ncond_int - y_ncond) / y_ncond_int * &
            lcond / (tinf - t_wall)

    sink_term = hcond * (tinf - t_wall) / lcond
    gam_s(ii,ipr) = gam_s(ii, ipr) - sink_term

  else ! no condensation

    sink_term = 0.d0
    theta = 1.0d0

  endif

  !-- Convective exchange coefficient from COPAIN correlation
  lambda = cpro_venth(iel)*cpro_cp(iel)
  call compute_prandtl(cpro_viscl(iel), cpro_venth(iel), Prdtl)
  call compute_exchange_adimensional(theta, Re_z, Gr_z, Prdtl, conv_regime, Nu_z)
  hw_cop = lambda * Nu_z / (cpro_cp(iel) * lcar)

  h2max = max(h2max,hw_cop)
  h2min = min(h2min,hw_cop)

  ! Convective exchange coefficient from turbulent wall law
  call hturbp(iwalfs,xnu,Prdtl,sigmat,rough_t,uk,yplus,dplus,hpflui,ypth)
  hw_enth = cpro_venth(iel)/distbf*hpflui

  h1max = max(h1max,hw_enth)
  h1min = min(h1min,hw_enth)

  if (izcophg(iz).eq.1) then
    hpcond(ii) = hw_enth
  else if (izcophg(iz).eq.2) then
    hpcond(ii) = hw_cop
  else if (izcophg(iz).eq.3) then
    hpcond(ii) = max(hw_cop,hw_enth)
  endif

  !Total thermal flux = condensation flux + convective flux
  flux = hpcond(ii)*cpro_cp(iel)*(tinf-t_wall)  &
       + sink_term*lcond

  ! Store it for postprocessing
  thermal_condensation_flux(ii) = flux

  flmin = min(flmin,flux)
  flmax = max(flmax,flux)

  !===================================================
  !== With the 1D thermal conduction model is used  ==
  !==       (iagt1d:=1), we stored the flux         ==
  !==         and its derivative.                   ==
  !===================================================
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
  call parmin(h2min)
  call parmin(h3min)
  call parmin(h4min)
  call parmax(h2max)
  call parmax(h3max)
  call parmax(h4max)
  call parmin(flmin )
  call parmax(flmax )
  call parsom(gamma_cond)
endif

if (mod(ntcabs,ntlist).eq.0) then
  write(nfecra,*) ' Minmax values of h_gas turb-enth hp  : ', &
                    ntcabs, h1min, h1max
  write(nfecra,*) ' Minmax values of copain hcop hp      : ', &
                    ntcabs, h2min, h2max
  write(nfecra,*) ' Minmax values of h_gas hturb-mix sink: ', &
                    ntcabs, h3min, h3max
  write(nfecra,*) ' Minmax values of copain hcop sink    : ', &
                    ntcabs, h4min, h4max
  write(nfecra,*) ' Minmax values of thermal flux        : ', &
                    ntcabs, flmin, flmax
endif

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

if (vcopt%iwarni.ge.1) then
  write(nfecra,1061) gamma_cond
endif

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
 1061 format(/,                                                   &
' ** Condensation source terms (Gamma) added:',E14.5           ,/,&
'    ----------------------------------------',/)

!----
! End
!----

return
end subroutine condensation_copain_model

subroutine hcd_copain(mass_diffusivity, sherwood, length, hcd)
  implicit none
  double precision, intent(in) :: mass_diffusivity, sherwood, length
  double precision, intent(out) :: hcd
  hcd = mass_diffusivity * sherwood / length
end subroutine hcd_copain


subroutine compute_drho(t_inf, t_wall, drho)
  implicit none
  double precision, intent(in) :: t_inf, t_wall
  double precision, intent(out) :: drho
  drho = abs((t_inf-t_wall)/(t_inf) )
end subroutine compute_drho
