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

subroutine iniusi

!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE APPELANT LES ROUTINES UTILISATEUR POUR L'ENTREE DES
!   PARAMETRES DE CALCUL : ON PASSE ICI POUR TOUT CALCUL

! CETTE ROUTINE PERMET DE CACHER A L'UTILISATEUR LES APPELS
!   A VARPOS ET AU LECTEUR XML DE L'IHM

! LE DECOUPAGE DE L'ANCIEN USINI1 PERMET EGALEMENT DE MIEUX
!   CONTROLER LES ZONES OU SONT INITIALISES LES VARIABLES (PAR
!   LE BIAIS DE PARAMETRES PASSES EN ARGUMENT)


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use atincl
use atsoil, only:nbrsol, tab_sol
use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use albase
use parall
use period
use pointe, only:ibm_porosity_mode
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use radiat
use cs_coal_incl
use cs_c_bindings
use cs_cf_bindings
use field
use cdomod

!===============================================================================

implicit none

! Arguments

! Local variables

integer          nmodpp
integer          nscmax
integer          l_size, f_id
integer          error, n_elts
double precision l_cp(1), l_xmasm(1), l_cv(1)
integer, dimension(:), pointer :: elt_ids

!===============================================================================

procedure() :: varpos, usppmo, uialin, cscpva, usipph, cfnmtd, fldvar, csivis
procedure() :: atini1, solcat, csidtv, csiphy, fldprp, cstime, usipsu
procedure() :: indsui

interface

  subroutine cs_function_default_define()  &
       bind(C, name='cs_function_default_define')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_function_default_define

  subroutine cs_gui_ale_diffusion_type()  &
       bind(C, name='cs_gui_ale_diffusion_type')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_ale_diffusion_type

  subroutine cs_gui_checkpoint_parameters()  &
       bind(C, name='cs_gui_checkpoint_parameters')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_checkpoint_parameters

  subroutine cs_combustion_initialize()  &
       bind(C, name='cs_combustion_initialize')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_combustion_initialize

  subroutine cs_gui_combustion_ref_values()  &
       bind(C, name='cs_gui_combustion_ref_values')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_combustion_ref_values

  subroutine cs_gui_mobile_mesh_structures_add()  &
       bind(C, name='cs_gui_mobile_mesh_structures_add')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_mobile_mesh_structures_add

  subroutine cs_gui_physical_constants()  &
       bind(C, name='cs_gui_physical_constants')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_physical_constants

  subroutine cs_gui_physical_properties()  &
       bind(C, name='cs_gui_physical_properties')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_physical_properties

  subroutine cs_gui_define_fans()  &
       bind(C, name='cs_gui_define_fans')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_define_fans

  subroutine cs_gui_equation_parameters()  &
       bind(C, name='cs_gui_equation_parameters')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_equation_parameters

  subroutine cs_gui_error_estimator()  &
       bind(C, name='cs_gui_error_estimator')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_error_estimator

  subroutine cs_gui_numerical_options()  &
       bind(C, name='cs_gui_numerical_options')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_numerical_options

  subroutine cs_gui_output_boundary()  &
       bind(C, name='cs_gui_output_boundary')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_output_boundary

  subroutine cs_gui_porous_model()  &
       bind(C, name='cs_gui_porous_model')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_porous_model

  subroutine cs_gui_radiative_transfer_parameters()  &
       bind(C, name='cs_gui_radiative_transfer_parameters')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_radiative_transfer_parameters

  subroutine cs_gui_thermal_model()  &
       bind(C, name='cs_gui_thermal_model')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_thermal_model

  subroutine cs_gui_scalar_model_settings()  &
       bind(C, name='cs_gui_scalar_model_settings')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_scalar_model_settings

  subroutine cs_mobile_structures_setup()  &
       bind(C, name='cs_mobile_structures_setup')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_mobile_structures_setup

  subroutine cs_velocity_pressure_set_solid()  &
       bind(C, name='cs_velocity_pressure_set_solid')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_velocity_pressure_set_solid

  subroutine cs_pressure_correction_model_activate()  &
       bind(C, name='cs_pressure_correction_model_activate')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_pressure_correction_model_activate

  subroutine cs_runaway_check_define_field_max(f_id, value)  &
       bind(C, name='cs_runaway_check_define_field_max')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: f_id
    real(c_double), value :: value
  end subroutine cs_runaway_check_define_field_max

  subroutine cs_user_radiative_transfer_parameters()  &
       bind(C, name='cs_user_radiative_transfer_parameters')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_user_radiative_transfer_parameters

  ! Interface to C function to initialize CDO model structures

  subroutine cs_f_domain_setup_init_model_context()  &
       bind(C, name='cs_f_domain_setup_init_model_context')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_f_domain_setup_init_model_context

  ! Interface to C function building properties

  subroutine cs_create_added_properties() &
    bind(C, name='cs_create_added_properties')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_create_added_properties

end interface

!===============================================================================

elt_ids => null()

! Check for restart and read matching time steps and notebook values

call parameters_read_restart_info

!===============================================================================
! 1. Initialize model settings
!===============================================================================

! Flow model selection through GUI

call cs_gui_physical_model_select

call cfnmtd(ficfpp, len(ficfpp))

! Flow model selection through user Fortran subroutine

call usppmo(1)
call cs_wall_condensation_set_onoff_state(icondb, icondv)

! Other models selection through GUI

! ALE parameters
call uialin (nalinf, nalimx, epalim)

! thermal model
call cs_gui_thermal_model

! turbulence model choice
call cs_gui_turb_model

! constant or variable specific heat, volume viscosity, ...
call csfpva

! Other models selection through user Fortran subroutine

call usipph(1, iturb, itherm, iale)

call csidtv()

call csiphy()

! Gravity and Coriolis
! Presence or not of gravity may be needed to determine whether some fields
! are created, so this is called before cs_user_model (to provide a
! user call site allowing to modify GUI-defined value programatically
! before property fields are created).
call cs_gui_physical_constants

! Activate radiative transfer model

! This module must be activated early so as to reserve the associated
! variables in some physical models.

call cs_gui_radiative_transfer_parameters
call cs_user_radiative_transfer_parameters ! deprecated, kept for compatibility

! Flow and other models selection through user C function
call cs_user_model

! Initialize some model structures if needed

call cs_combustion_initialize
call cs_gui_combustion_ref_values

! Set type and order of the turbulence model
call cs_set_type_order_turbulence_model()

! If CDO has been activated, initialize the context structures for
! models which have been activated
if (icdo.ge.1) then
   call cs_f_domain_setup_init_model_context
endif

! Activate CDO for ALE
if (iale.eq.2) then
  call cs_ale_activate
endif

if (iale.ge.1) then
  call cs_gui_mobile_mesh_structures_add
endif

! Read thermochemical data for specific physics
call pplecd

! Other model parameters, including user-defined scalars

call cs_gui_user_variables
call cs_gui_user_arrays
call cs_gui_calculator_functions

! Solid zones

call cs_velocity_pressure_set_solid

!===============================================================================
! 2. Initialize parameters for specific physics
!===============================================================================

call cs_rad_transfer_options

! Define fields for variables, check and build iscapp
! and computes the number of user scalars (nscaus)
if (icdo.lt.2) then
  call fldvar(nmodpp)

  ! Activate the pressure correction model only if CDO mode is not stand-alone
  call cs_pressure_correction_model_activate
endif

if (iale.ge.1) then
  call cs_gui_ale_diffusion_type
endif
call csivis

nscmax = nscamx

!===============================================================================
! Specific physics modules
!===============================================================================
! Note: part of what is inside ppini1 could be moved here
! so that usipsu / cs_user_parameters can be used by the user to modify default
! settings

! Atmospheric flows
if (ippmod(iatmos).ge.0) then
  call atini1
endif

! Compressible
call cscfgp
call field_get_id_try('velocity', f_id)
if (f_id .ge. 0) then
  if (ippmod(icompf).ge.0) then
    call cs_runaway_check_define_field_max(f_id, 1.0d5)
  else
    call cs_runaway_check_define_field_max(f_id, 1.0d4)
  endif
endif

! Atmospheric module
if (ippmod(iatmos).ge.0) then
  ! Some advanced init/allocation for the soil model
  if (iatsoil.ge.0) then

    ! Get the number of sol only (warning, number of element of the zone not yet
    ! computed)
    call atmo_get_soil_zone(n_elts, nbrsol, elt_ids)

    ! Allocation of table of values
    allocate(tab_sol(nbrsol),stat = error)

    if (error /= 0) then
      write(nfecra,*) "Allocation error of atmodsol::tab_sol"
      call csexit(1)
    endif

    ! First pass, default soil categories parameters function of the number
    ! of soil. Can be modified by the user
    call solcat(1)

  endif
endif

!===============================================================================
! 3. INITIALISATION DE PARAMETRES "GLOBAUX"
!===============================================================================

! --- Parametres globaux

!     Pas de temps
!     Couplage vitesse/pression
!     Prise en compte de la pression hydrostatique
!     Estimateurs (pas encore dans l'IHM)


!   - Interface code_saturne
!     ======================

! Postprocessing

call cs_gui_output_boundary

! Define main properties (pointers, checks, ipp) if not in CDO mode only
if (icdo.lt.2) then
  call fldprp
endif

!===============================================================================
! 4. INITIALISATION DE PARAMETRES UTILISATEUR SUPPLEMENTAIRES
!===============================================================================

! --- Format des fichiers aval (entsor.h)
! --- Options du calcul (optcal.h)
! --- Constantes physiques (cstphy.h)

!   - Interface code_saturne
!     ======================

! Restart, read auxiliary file, frozen velocity field

call cs_gui_checkpoint_parameters()

! Time step (only ntmabs, dtref)
call cstime()

! Local numerical options

call cs_gui_equation_parameters

! If CDO mode only, no pressure is defined at this stage
if (icdo.lt.2) then
  call cs_gui_numerical_options
endif

! Physical properties
call cs_gui_physical_properties

! Turbulence reference values (uref, almax)
call cs_gui_turb_ref_values

! Set turbulence constants according to model choices.
! This can be overwritten by the user in cs_user_parameters()
call cs_f_turb_complete_constants(-1)

! Scamin, scamax, turbulent flux model, diffusivities
call cs_gui_scalar_model_settings()

! Porosity model
call cs_gui_porous_model()

! Init fan
call cs_gui_define_fans()

! Init error estimator
call cs_gui_error_estimator()

! Initialize base evaluation functions
call cs_function_default_define()

!   - User functions
!     ==============

call usipsu(nmodpp)
call user_parameters

! If time step is local or variable, pass information to C layer, as it
! may be needed for some field (or moment) definitions.
if (idtvar.ne.0) then
  call time_step_define_variable(1)
endif
if (idtvar.eq.2.or.idtvar.eq.-1) then
  call time_step_define_local(1)
endif

call indsui(isuite)

! Default values of physical properties for the compressible model
if (ippmod(icompf).ge.0) then
  ! ieos has been set above in uippmo with the GUI or in usppmo without the GUI.
  ! The variability of the thermal conductivity
  ! (diffusivity_id for itempk) and the volume viscosity (iviscv) has
  ! been set in fldprp.

  ! Compute cv0 according to chosen EOS.
  l_size = 1
  l_cp(1) = cp0 ! dummy argument in stiffened gas
  l_xmasm(1) = xmasmr ! dummy argument in stiffened gas
  call cs_cf_thermo_cv(l_cp, l_xmasm, l_cv, l_size)
  cv0 = l_cv(1)
endif

if (ibm_porosity_mode.gt.0) then
  iporos = 3
endif

! Choose the porous model
call cs_porous_model_set_model(iporos)

! --- Varpos
! If CDO mode only, skip this stage
if (icdo.lt.2) then
  call varpos
endif
! --- Internal coupling
call cs_gui_internal_coupling
call cs_user_internal_coupling

call cs_internal_coupling_setup

! Mobile structures
! (after call to cs_gui_mobile_mesh_structures_add and possible
! call by user to cs_mobile_structures_add_n_structures)
if (iale.ge.1) then
  call cs_mobile_structures_setup
endif

!----
! Formats
!----

return
end subroutine
