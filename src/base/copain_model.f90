!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
   gam_s  , hpcond )

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
use pointe, only: thermal_condensation_flux, flthr, dflthr
use parall
use period
use field
use mesh
use cs_c_bindings
use cs_f_interfaces
use cs_nz_condensation, only: izcophc, izcophg, iztag1d, ztpar
use cs_nz_tagmr, only: ztpar0, ztmur


!===============================================================================

implicit none

! Arguments

integer          nvar, nfbpcd, ifbpcd(nfbpcd), izzftcd(nfbpcd)

double precision gam_s(nfbpcd,nvar)
double precision hpcond(nfbpcd)

! Local variables

integer          ii, iz, iel, ifac, iesp
integer          ivar, f_id, ifcvsl

double precision flux
double precision sink_term, gamma_cond
double precision lambda, theta, tinf, psat
double precision Sh_z, Nu_z, Gr_z, schdt, Prdtl
double precision lcar, l_cell, xnu
double precision x_inc,y_ncond, x_k
double precision x_vapint,x_ncond_int,y_ncond_int
double precision distbf
double precision ratio_tkpr,drho
double precision xkloc
double precision xmab,xvab,patm,a1
double precision g_n
double precision dplus, yplus, sigmat, ypth
double precision hcond,hcdt,hcdcop,hw_cop,hflui,hpflui,hw_enth
double precision h1min, h1max, h2min, h2max
double precision h3min, h3max, h4min, h4max
double precision flmin, flmax
double precision C_k1,C_k2,C_k3,C_k4,C_k5,C_k6
double precision C_k7,C_k8,C_k9,t_wall
double precision pr_c,t_c,dtheta
double precision lcond

type(gas_mix_species_prop) s_h2o_g, s_k
type(var_cal_opt) :: vcopt

double precision, allocatable, dimension(:) :: mix_mol_mas, mol_mas_ncond
double precision, allocatable, dimension(:) :: x_ncond, x_h2o_g, diff_m
double precision, dimension(:), pointer :: cpro_rho, cpro_viscl, cpro_cp, cpro_venth
double precision, dimension(:), pointer :: cvar_enth, cvar_yk
double precision, dimension(:), pointer :: y_h2o_g
double precision, dimension(:), pointer :: yplbr

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
!     (lcar, l_cell) the characteristic length scales
!===============================================================================

pr_c = 221.2d+5 ; t_c  = 647.3d0
patm = 101320.d0

C_k1 =  -7.691234564d0 ; C_k2 = -26.08023696d0
C_k3 =-168.1706546d0   ; C_k4 =  64.23285504d0
C_k5 =-118.9646225d0   ; C_k6 =   4.16711732d0
C_k7 =  20.9750676d0   ; C_k8 =  -1.d+9
C_k9 =   6.d0

lcar = 1.d0 ;l_cell = 1.d0

!- Define the vaporization Latent heat (lambda_cond)
!------------------------------------------------------
lcond = 2278.0d+3

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

!===============================================================================
! 3 - Compute mass fraction of H2O, etc.
!===============================================================================

! Compute the mass fraction of H2O vapor
! Deduced from the non-condensable species (O2, N2)
!--------------------------------------------------

! Initialization
do iel = 1, ncel
  y_h2o_g(iel) = 1.d0
  ! Mixture molar mass
  mix_mol_mas(iel) = 0.d0
  ! Molar mass location of noncondensable gas
  mol_mas_ncond(iel) = 0.d0
  ! Molar mass fraction of noncondensable gas
  x_ncond(iel) = 0.d0
  ! Mixture molecular diffusivity
  ! between the vapor and the noncondensable gases
  diff_m(iel) = 0.d0
enddo

do iesp = 1, nscasp
  ! mass fraction array of the different
  ! species (O2, N2)
  !-------------------------------------------------
  call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)

  call field_get_key_struct_gas_mix_species_prop( &
       ivarfl(isca(iscasp(iesp))), s_k)

  do iel = 1, ncel
    y_h2o_g(iel) = y_h2o_g(iel)-cvar_yk(iel)
    mix_mol_mas(iel) = mix_mol_mas(iel) + cvar_yk(iel)/s_k%mol_mas
  enddo
enddo

!Finalize the computation of the Mixture molar mass
do iel = 1, ncel
  mix_mol_mas(iel) = mix_mol_mas(iel) + y_h2o_g(iel)/s_h2o_g%mol_mas
  mix_mol_mas(iel) = 1.d0/mix_mol_mas(iel)
enddo

! Warning: only on non-condensable gases
do iesp = 1, nscasp
  call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)
  call field_get_key_struct_gas_mix_species_prop( &
       ivarfl(isca(iscasp(iesp))), s_k)
  do iel = 1, ncel
    x_k = cvar_yk(iel)*mix_mol_mas(iel)/s_k%mol_mas
    mol_mas_ncond(iel) = mol_mas_ncond(iel) + x_k*s_k%mol_mas
    x_ncond(iel) = x_ncond(iel) + x_k
  enddo
enddo

! Molar fraction of h20 vapor
do iel = 1, ncel
  x_h2o_g(iel) = y_h2o_g(iel)*mix_mol_mas(iel)/s_h2o_g%mol_mas
enddo

! Finalize the molar mass of non-condensable species
do iel = 1, ncel
  mol_mas_ncond(iel) = mol_mas_ncond(iel)/x_ncond(iel)
enddo

!-- Compute the mixing molecular diffusivity
!-- between the vapor and the nocondensable gases
!-------------------------------------------------

! Warning: only on non-condensable gases
do iesp = 1, nscasp

  call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)
  call field_get_key_struct_gas_mix_species_prop( &
       ivarfl(isca(iscasp(iesp))), s_k)

  do iel = 1, ncel

    x_k = cvar_yk(iel)*mix_mol_mas(iel)/s_k%mol_mas

    if (idilat.eq.3) then
      ratio_tkpr = ((cvar_enth(iel)/cpro_cp(iel))**1.75d0)/pther
    else
      ratio_tkpr = ((cvar_enth(iel)/cpro_cp(iel))**1.75d0)/p0
    endif

    xmab = sqrt(2.d0/( 1.d0/(s_h2o_g%mol_mas*1000.d0) &
                      +1.d0/(    s_k%mol_mas*1000.d0) ) )
    xvab = ( s_h2o_g%vol_dif**(1.d0/3.d0) &
            +    s_k%vol_dif**(1.d0/3.d0) )**2.d0
    a1   = 1.43d-7/(xmab*xvab)*patm

    diff_m(iel) = diff_m(iel) + x_k/(a1*ratio_tkpr)
  enddo
enddo

! Finalization
do iel = 1, ncel
  diff_m(iel) = x_ncond(iel) / diff_m(iel)
enddo

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
  !-- Geometric quantities -------------------------
  distbf = distb(ifac)

  !-- If the 1D thermal conduction model is activated,
  !-- the wall temperature is in unit (Celsius °C)
  !---------------------------------------------------
  if(iztag1d(iz).eq.1) then
    if(isuite.eq.0.and.ntcabs.eq.1) then
      t_wall = ztpar0(iz)
    else
      t_wall = ztmur(ii,1)
    endif
  else
    t_wall = ztpar(iz)
  endif

  !-- kinematic viscosity --------------------------
  xnu = cpro_viscl(iel)/cpro_rho(iel)

  !-- Compute the t_in temperature far away
  !-- from the wall. t_inf := t_atmos, is the
  !-- temperature associate to the cell.
  !--      t_inf given in unit (°C)
  !-------------------------------------------------

  tinf = cvar_enth(iel)/cpro_cp(iel) - tkelvi

  !-- the saturated pressure associated
  !-- dtheta = dT/T_c
  !-------------------------------------------------
  dtheta = (t_wall+tkelvi)/t_c
  psat  = pr_c*exp((1.d0/dtheta)       &
        *( C_k1*(1.d0-dtheta)          &
          +C_k2*(1.d0-dtheta)**2       &
          +C_k3*(1.d0-dtheta)**3       &
          +C_k4*(1.d0-dtheta)**4       &
          +C_k5*(1.d0-dtheta)**5)      &
        / (1.d0+C_k6*(1.d0-dtheta)     &
               +C_k7*(1.d0-dtheta)**2) &
      -(1.d0-dtheta)/(C_k8*(1.d0-dtheta)**2+C_k9))
  !-- Prandtl number -------------------------------
  Prdtl = cpro_viscl(iel)/ cpro_venth(iel)
  !-- Grasholf number ------------------------------
  drho = abs((tinf-t_wall)/(tinf+tkelvi) )
  Gr_z = sqrt(gx**2+gy**2+gz**2)*drho*lcar**3/(xnu**2)
  !-- Nusselt number of the natural convection -----
  theta = 1.d0
  Nu_z = theta*(0.13d0*(Gr_z*Prdtl)**(1.d0/3.d0))

  !-- Molar fraction of the noncondensable gases (Xinc,inf)
  !-------------------------------------------------
  x_inc = 1.d0-x_h2o_g(iel)

  !-- Molar fraction of the interface vapor
  !-------------------------------------------------

  if (idilat.eq.3) then
    x_vapint = psat/pther
  else
    x_vapint = psat/p0
  endif

  ! if (Xv > Xi,v) we have condensation  -----------
  if (x_h2o_g(iel).gt.x_vapint) then
    ! noncondensable Molar fraction at the interface
    x_ncond_int = 1.d0 - x_vapint

    ! Mixture molar mass at the interface -----------
    xkloc = x_vapint*s_h2o_g%mol_mas + x_ncond_int*mol_mas_ncond(iel)

    ! Noncondensable mass fraction
    y_ncond =  1.d0 - y_h2o_g(iel)

    ! Mass fraction of noncondensable gas at interface
    y_ncond_int = x_ncond_int*mol_mas_ncond(iel)/xkloc

    ! Computation of the corrective factor theta
    ! using for Nusselt and Sherwood numbers
    !-----------------------------------------------
    theta = 1.d0 +0.625d0*(x_ncond_int-x_inc)/x_ncond_int

    ! Nusselt number of natural convection
    !-----------------------------------------------
    Nu_z = theta*(0.13d0*(Gr_z*Prdtl)**(1.d0/3.d0))

    ! Sink source term (kg/s) -- Initialization
    !-----------------------------------------------
    sink_term = 0.d0

    !-- Schmidt number
    schdt = xnu/diff_m(iel)

    !-- The McAdams correlation gives the following
    !-- expression for the Sherwood number
    !-------------------------------------------------
    Sh_z  = theta*(0.13d0*(Gr_z*schdt)**(1.d0/3.d0))

    !-- [1]. Computation of the (hcdt) turbulent
    !-- exchange coefficient of condensation based on
    !-- the schmidt number and the molecular
    !-- diffusivity of the mixing gas
    !-------------------------------------------------
    dplus = 0.d0
    yplus = yplbr(ifac)
    sigmat = 0.9d0
    call hturbp(iwalfs,schdt,sigmat,yplus,dplus,hflui,ypth)

    hcdt =  diff_m(iel)/distbf*hflui

    h3max = max(h3max,hcdt)
    h3min = min(h3min,hcdt)

    !-- [2]. Computation of the (hcdcop) exchange
    !-- coefficient of condensation based on
    !--           the COPAIN correlation
    !-------------------------------------------------
    hcdcop = diff_m(iel)*Sh_z/(l_cell)

    h4max = max(h4max,hcdcop)
    h4min = min(h4min,hcdcop)

    !-- Computation of (hcond) the exchange
    !-- coefficient of condensation function of the
    !-- the user parameter choice (icophc).
    !-------------------------------------------------
    if (izcophc(iz).eq.1) then

      !-- [1]. Choose the (hcdt) turbulent exch. coeff
      !-----------------------------------------------
      hcond = hcdt
    else if (izcophc(iz).eq.2) then

      !-- [2]. Choose the (hcdcop) COPAIN  exch. coeff
      !-----------------------------------------------
      hcond = hcdcop
    else if (izcophc(iz).eq.3) then

      !-- [3]. Choose the maximum value between the
      !--    compute turbulent coefficient and the
      !--    theorical one given by Copain correlation
      !-----------------------------------------------
      hcond = max(hcdt,hcdcop)
    endif

    !=================================================
    !== Computation of sink source term gam_s(ii)
    !==         per unit of area (kg/s/m2)
    !==     ----------------------------------------
    !==       (if Xv > Xi,v we have condensation )
    !=================================================

    sink_term = cpro_rho(iel)*hcond   &
              *((y_ncond_int-y_ncond)/y_ncond_int)
    gam_s(ii,ipr) = gam_s(ii, ipr) - sink_term

  else

    !================================================
    !==   The sink source term (kg.s-1.m-2) is equal
    !==   to zero :
    !==           gam_s(ii) = 0.
    !==   (if Xv < Xi,v we do not have condensation )
    !================================================

    sink_term = 0.d0
  endif

  !===================================================
  !== Computation of the thermal exchange
  !== coefficient (kg.m-1.s-1) to the wall surface
  !===================================================

  lambda = cpro_venth(iel)*cpro_cp(iel)

  !-- Nusselt number
  Nu_z = theta*(0.13d0*(Gr_z*Prdtl)**(1.d0/3.d0))
  !-- the gravity term value
  g_n  = sqrt(gx**2+gy**2+gz**2)


  !--  Compute the thermal exchange coefficient
  !--      with the Copain correlation
  !---------------------------------------------------
  hw_cop = 0.13d0*theta*lambda                     &
         *( g_n*drho*Prdtl/(xnu**2) )**(1.d0/3.d0) &
         /cpro_cp(iel)

  h2max = max(h2max,hw_cop)
  h2min = min(h2min,hw_cop)


  !-- Computation of (hw_cop) the thermal exchange
  !-- coefficient to the wall surface, function of
  !-- the user parameter choice (icophg).
  !-------------------------------------------------
  dplus = 0.d0
  yplus = yplbr(ifac)
  sigmat = 0.9d0
  Prdtl = cpro_viscl(iel)/cpro_venth(iel)

  call hturbp(iwalfs,Prdtl,sigmat,yplus,dplus,hpflui,ypth)
  hw_enth = cpro_venth(iel)/distbf*hpflui

  h1max = max(h1max,hw_enth)
  h1min = min(h1min,hw_enth)
  if (izcophg(iz).eq.1) then

    !-- [1]. Choose the (hw_enth) exch. coeff
    hpcond(ii) = hw_enth
  else if (izcophg(iz).eq.2) then

    !-- [2]. Choose the (hp_cop) exch. coeff
    !-----------------------------------------------
    hpcond(ii) = hw_cop
  else if (izcophg(iz).eq.3) then

    !-- [3]. Choose the maximal value of
    !--         (hp_cop) exch. coeff
    !-----------------------------------------------
    hpcond(ii) = max(hw_cop,hw_enth)
  endif

  !===================================================
  !== Computation of the thermal flux to the face   ==
  !===================================================

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
end subroutine condensation_copain_model
