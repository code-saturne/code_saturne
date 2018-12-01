!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file metal_structures_copain_model.f90
!>
!> \brief The COPAIN modelling to estimate the heat and mass transfer
!> associated to the steam condensation phenomena at each cell corresponding to
!> the metal structures volume identified by geometric criteria.
!>
!> This subroutine is used to compute at each cell the
!> \f$ \Gamma^{v}_{\mbox{cond}} \f$ and \f$ \phi^{v}_{\mbox{cond}} \f$.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ncmast       number of cells with metal structures condensation
!> \param[in]     ltmast       index of cells with metal structures condensation
!> \param[in]     tmet         temperature imposed at the cold wall
!>                             as constant or variable in time
!>                             with a 0-D metal thermal model
!> \param[out]    gam_ms       value associated to each variable in the
!>                             condensation source terms \f$ \Gamma^v_{cond}\f$
!> \param[out]    flux_ms      value associated to the heat transfer estimated
!>                             for the metal structures \f$\phi^v_{cond}\f$.
!_______________________________________________________________________________

subroutine metal_structures_copain_model &
 ( ncmast ,  ltmast ,              &
   tmet   ,                        &
   gam_ms, flux_ms )

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
use cs_tagms, only:s_metal, t_metal, tmet0

!===============================================================================

implicit none

! Arguments

integer          ncmast ,  ltmast(ncelet)

double precision tmet
double precision gam_ms(ncelet)
double precision flux_ms(ncelet)

! Local variables

integer          icmet,  iel, iesp
integer          ivar, f_id, ifcvsl

double precision flux
double precision sink_term, gamma_cond
double precision lambda, theta, tinf, psat
double precision Sh_z, Nu_z, Gr_z, schdt, Prdtl
double precision lcar, l_cell, xnu
double precision x_inc,y_ncond, x_k
double precision x_vapint,x_ncond_int,y_ncond_int
double precision ratio_tkpr,drho
double precision xkloc
double precision xmab,xvab,patm,a1
double precision flmin, flmax
double precision C_k1,C_k2,C_k3,C_k4,C_k5,C_k6
double precision C_k7,C_k8,C_k9,t_wall
double precision pr_c,t_c,dtheta
double precision lcond
double precision surfbm, volm

type(gas_mix_species_prop) s_h2o_g, s_k
type(var_cal_opt) :: vcopt

double precision, allocatable, dimension(:) :: mix_mol_mas, mol_mas_ncond
double precision, allocatable, dimension(:) :: x_ncond, x_h2o_g, diff_m
double precision, dimension(:), pointer :: cpro_rho, cpro_viscl, cpro_cp, cpro_venth
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: cvar_enth, cvar_yk
double precision, dimension(:), pointer :: y_h2o_g

!===============================================================================
! Allocate a temporary array for cells selection
allocate(mix_mol_mas(ncelet), mol_mas_ncond(ncelet), x_ncond(ncelet))
allocate(x_h2o_g(ncelet), diff_m(ncelet))

!===============================================================================
! 0 - Initialization
!===============================================================================
flmin = 1.d+20 ; flmax = -1.d+20

sink_term = 0.d0

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

drho  = 0.d0 ; Prdtl = 1.d0
schdt = 1.d0 ; xnu   = 1.d0

! Define the vaporization Latent heat (Lcond)
! -------------------------------------------
lcond = 2278.0d+3

! Define the global metal structures volume
! -----------------------------------------
do icmet = 1, ncmast
  iel =  ltmast(icmet)
  volm = volm + volume(iel)
enddo

if (irangp.ge.0) then
  call parsom(volm)
endif

!===============================================================================
! 2 - Pointers to phys. properties for gas mixture at cell center.
!===============================================================================

! Coefficient exchange of the enthalpy scalar
ivar = isca(iscalt)
call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

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
! Deduced from the noncondensable species (O2, N2)
!-------------------------------------------------

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

! Warning: only on noncondensable gases
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

! Finalize the molar mass of noncondensable species
do iel = 1, ncel
  mol_mas_ncond(iel) = mol_mas_ncond(iel)/x_ncond(iel)
enddo

!-- Compute the mixing molecular diffusivity
!-- between the vapor and the nocondensable gases
!-------------------------------------------------

! Warning: only on noncondensable gases
do iesp = 1, nscasp

  call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)
  call field_get_key_struct_gas_mix_species_prop( &
       ivarfl(isca(iscasp(iesp))), s_k)

  do iel = 1, ncel

    x_k = cvar_yk(iel)*mix_mol_mas(iel)/s_k%mol_mas

    ratio_tkpr = ((cvar_enth(iel)/cpro_cp(iel))**1.75d0)/pther

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
! 4 - Loop on the (ii) cells with metal structures condensation
!     ---------------------------------------------------------
!     -> the sink source term gam_ms(ii) associated to the condensation
!        is computed at each cell associated to the metal structure volume
!        (loss of vapor mass flow rate due to condensation)
!     -> the heat transfer flux_ms(ii) associated to the phase change
!        is computed at each cell and adding in the explicit enthalpy
!        source term as a volume heat transfer flux.
!        ------------------------------------------------------------
!         the correlations used are similar to the copain modelling
!===============================================================================

do icmet = 1, ncmast

  iel =  ltmast(icmet)
  !-- Geometric quantities ---------------------------
  surfbm = s_metal*volume(iel)/volm

  !-- If the 0-D thermal conduction model is activated,
  !-- the wall temperature is in unit (Celsius °C)
  !---------------------------------------------------
  if(itagms.eq.1) then
    if(isuite.eq.0.and.ntcabs.eq.1) then
      t_wall = tmet0
    else
      t_wall = t_metal(iel,1)
    endif
  else
    t_wall = tmet
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

  x_vapint = psat/pther

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

    !=================================================
    !== Computation of sink source term gam_s(ii)
    !==         per unit of area (kg/s/m2)
    !==     ----------------------------------------
    !==       (if Xv > Xi,v we have condensation )
    !=================================================

    sink_term = cpro_rho(iel)*diff_m(iel)*Sh_z     &
              *((y_ncond_int-y_ncond)/             &
                 y_ncond_int )/(l_cell)
    gam_ms(iel) = gam_ms(iel) - sink_term

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

  !===================================================
  !== Computation of the thermal flux to the face   ==
  !===================================================
  flux =lambda*Nu_z*(tinf-t_wall)*surfbm

  flmin = min(flmin,flux)
  flmax = max(flmax,flux)
  !===================================================
  !==  With the 0-D thermal model is used           ==
  !==  we stored the flux to pass to the model      ==
  !===================================================
  if(itagms.eq.1) then
    flux_ms(iel) = flux
  endif
  !===================================================

enddo !  End of the Loop on the (ii) cells with metal structures condensation

! Deallocate the temporary array
deallocate(mix_mol_mas, x_ncond, mol_mas_ncond)
deallocate(x_h2o_g, diff_m)

!===============================================================================
! 5 - Print in the log every ntlist of the min/max values
!         for exchange coefficient and flux computed
!===============================================================================

gamma_cond = 0.d0
do icmet = 1, ncmast
  iel =  ltmast(icmet)
  gamma_cond = gamma_cond + gam_ms(iel)
enddo

if (irangp.ge.0) then
  call parmin(flmin )
  call parmax(flmax )
  call parsom(gamma_cond)
endif

if (mod(ntcabs,ntlist).eq.0) then
  write(nfecra,*) ' Minmax values of metal structures thermal flux : ', &
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
end subroutine metal_structures_copain_model
