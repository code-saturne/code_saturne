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
! ANY WAp0 / (RRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
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

subroutine condensation_dehbi_model &
 ( nvar   , nfbpcd , ifbpcd , izzftcd ,  &
   p_gam_s  , hpcond , regime) &
  bind(C, name="condensation_dehbi_model")

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
use cs_nz_condensation, only: flthr, dflthr, thermal_condensation_flux, iztag1d, ztpar, zxrefcond, zprojcond
use cs_nz_tagmr, only: ztpar0, ztmur

use condensation_module


!===============================================================================

implicit none

! Arguments

integer(c_int), value :: nvar, nfbpcd, regime
integer(c_int), dimension(*) :: ifbpcd, izzftcd

real(kind=c_double), dimension(*) :: hpcond
type(c_ptr) :: p_gam_s
double precision, dimension(:,:), pointer :: gam_s

! Local variables

integer          ii, iz, iel, ifac, iesp, idims
integer          ivar, f_id, ifcvsl, idir
integer          conv_regime

double precision gravity, pressure
double precision flux, h_dehbi
double precision t_wall
double precision sink_term, gamma_cond
double precision lambda, theta, tinf, psat
double precision Sh_z, Nu_z, Gr_z, schdt, Prdtl
double precision Sh_z_NC, Sh_z_FC, Re_z
double precision Nu_z_NC, Nu_z_FC
double precision lcar, xnu, mu, u_ref, u_square, u_norm
double precision x_inc,y_ncond, x_k
double precision x_vapint,x_ncond_int,y_ncond_int, B
double precision ratio_tkpr,drho, rho_av
double precision mol_mas_int, rho_wall
double precision hcond,hcdcop
double precision h1min, h1max, h2min, h2max
double precision h3min, h3max, h4min, h4max
double precision flmin, flmax

type(gas_mix_species_prop) s_h2o_g, s_k
type(var_cal_opt) :: vcopt

double precision, allocatable, dimension(:) :: mix_mol_mas, mol_mas_ncond
double precision, allocatable, dimension(:) :: x_ncond, x_h2o_g, diff_m
double precision, dimension(:), pointer :: cpro_rho, cpro_viscl, cpro_cp, cpro_venth
double precision, dimension(:), pointer :: cvar_enth, cvar_yk
double precision, dimension(:), pointer :: y_h2o_g
double precision, dimension(:), pointer :: yplbr
double precision, dimension(:,:), pointer :: cvar_vel

!===============================================================================
! C pointer to table bindings
call c_f_pointer(p_gam_s, gam_s, [nfbpcd,nvar])

!===============================================================================
! Allocate a temporary array for cells selection
allocate(mix_mol_mas(ncelet), mol_mas_ncond(ncelet), x_ncond(ncelet))
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

!===============================================================================
! 1 - Constant physical values of the Copain modelling
!     with (C_k) constants of the saturated pressure and
!     (lcar) the characteristic length scales
!===============================================================================

! Convection regime (1: natural, 2: forced, 3:mixed)
conv_regime = regime + 1

!===============================================================================
! 2 - Pointers to phys. properties for mixing flow at cell center.
!===============================================================================

ivar = isca(iscalt)

call field_get_val_s(ivarfl(ivar), cvar_enth)
call field_get_val_s_by_name("y_h2o_g", y_h2o_g)
call field_get_val_s(icrom, cpro_rho)
call field_get_val_s(iviscl, cpro_viscl)
call field_get_val_v(ivarfl(iu), cvar_vel)
if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
else
  write(nfecra,1000) icp
  call csexit (1)
endif
call field_get_key_int (ivarfl(isca(iscalt)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_venth) !lambda/cp
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
! 3 - Compute mixture properties
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
!          the correlations used are those from Dehbi 2015
!===============================================================================

lcar = 1.0d0
do ii = 1, nfbpcd

  ifac= ifbpcd(ii)
  iel = ifabor(ifac)

  iz  = izzftcd(ii)

  Re_z = 0.0d0
  ! Compute reference distance (vertical distance since start of condenser plate)
  if (conv_regime > 1) then
    call normalize_vector(3, zprojcond(1:3, iz))
    call compute_characteristic_length(cdgfbo(1:3, ifac), zxrefcond(1:3, iz), &
                                       zprojcond(1:3, iz), lcar)
    call compute_tangential_velocity(cvar_vel(1:3, iel), surfbo(1:3, ifac), 1.0d0/surfbn(ifac), u_ref)
    Re_z = cpro_rho(iel) * u_ref * lcar / cpro_viscl(iel)
  endif

  call get_wall_temperature(iz, ii, t_wall)

  xnu = cpro_viscl(iel)/cpro_rho(iel)
  mu  = cpro_viscl(iel)

  call compute_prandtl(cpro_viscl(iel), cpro_venth(iel), Prdtl)

  theta = 1.d0
  call get_temperature(cvar_enth(iel), cpro_cp(iel), tinf)

  call compute_psat(t_wall, psat)
  x_vapint = psat/pressure

  ! if (Xv > Xi,v) we have condensation
  if (x_h2o_g(iel).gt.x_vapint) then

    lambda = cpro_venth(iel)*cpro_cp(iel)

    x_ncond_int = 1.d0 - x_vapint
    x_inc = 1.d0-x_h2o_g(iel)
    mol_mas_int = x_vapint*s_h2o_g%mol_mas + x_ncond_int*mol_mas_ncond(iel)
    y_ncond =  1.d0 - y_h2o_g(iel)
    y_ncond_int = x_ncond_int*mol_mas_ncond(iel)/mol_mas_int

    rho_wall = pressure*mol_mas_int / (cs_physical_constants_r*t_wall)
    rho_av = 0.5d0 * (rho_wall + cpro_rho(iel))
    drho = dabs(rho_wall - cpro_rho(iel)) / rho_av

    B = (y_ncond - y_ncond_int) / (y_ncond_int)
    theta = 1.33d0 * LOG(1.0d0 + B) / B ! theta for global flux (from dehbi)
    xnu = mu / rho_av
    call compute_schmidt(xnu, diff_m(iel), schdt)
    call compute_grashof(gravity, drho, lcar, xnu, Gr_z)

    ! McAdams/schlichting correlation for Sherwood number
    call compute_exchange_adimensional(theta, Re_z, Gr_z, schdt, conv_regime, Sh_z)
    call compute_exchange_adimensional(theta, Re_z, Gr_z, Prdtl, conv_regime, Nu_z)

    ! Sensible heat flux, computed with COPAIN expression
    hpcond(ii) = Nu_z * lambda / (lcar * cpro_cp(iel))

    hcdcop = diff_m(iel)*Sh_z/lcar
    h_dehbi = hcdcop * rho_av * (y_ncond_int -y_ncond) / y_ncond_int &
              * lcond / (tinf - t_wall)

    h3max = max(h3max,h_dehbi)
    h3min = min(h3min,h_dehbi)

    hcond = h_dehbi - hpcond(ii) * cpro_cp(iel)

    !=================================================
    !== Computation of sink source term gam_s(ii)
    !==         per unit of area (kg/s/m2)
    !==     ----------------------------------------
    !==       (if Xv > Xi,v we have condensation )
    !=================================================

    sink_term = hcond * (tinf - t_wall) / lcond
    gam_s(ii,ipr) = gam_s(ii, ipr) - sink_term

    flux = h_dehbi * (tinf - t_wall)

  else

    !================================================
    !==   The sink source term (kg.s-1.m-2) is equal
    !==   to zero :
    !==           gam_s(ii) = 0.
    !==   (if Xv < Xi,v we do not have condensation )
    !================================================



    sink_term = 0.d0
    hcond = 0.0d0
    !-- Grasholf number based on temperature if no condensation ---
    drho = abs((tinf-t_wall)/tinf )
    theta = 1.0d0
    call compute_grashof(gravity, drho, lcar, xnu, Gr_z)
    call compute_exchange_adimensional(theta, Re_z, Gr_z, Prdtl, conv_regime, Nu_z)
    hpcond(ii) = Nu_z * lambda / (lcar * cpro_cp(iel))
    flux = hpcond(ii)*cpro_cp(iel)*(tinf-t_wall)
  endif

  h2max = max(h2max,hpcond(ii)*cpro_cp(iel))
  h2min = min(h2min,hpcond(ii)*cpro_cp(iel))

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
deallocate(mix_mol_mas, x_ncond, mol_mas_ncond)
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
  call parmax(h2max)
  call parmax(h3max)
  call parmin(flmin )
  call parmax(flmax )
  call parsom(gamma_cond)
endif

if (mod(ntcabs,ntlist).eq.0) then
  write(nfecra,*) ' Minmax values of sensible htc: ', &
                    ntcabs, h2min, h2max
  write(nfecra,*) ' Minmax values of dehbi htc: ', &
                    ntcabs, h3min, h3max
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
end subroutine condensation_dehbi_model
