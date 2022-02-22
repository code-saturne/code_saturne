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

!> \file ssh-aerosol.f90
!> \brief Module for aerosol chemistry in the atmospheric module

module sshaerosol

  !=============================================================================

  use, intrinsic :: iso_c_binding
  use ppppar, only: nozppm

  !=============================================================================
  !> \defgroup at_aerosol_chemistry Aerosol chemistry parameters for the
  !>                                atmospheric module

  !> \addtogroup at_aerosol_chemistry
  !> \{

  !> Enum and flag to activate or not the aerosol model
  enum, bind(C)
    enumerator :: CS_ATMO_AEROSOL_OFF = 0
    enumerator :: CS_ATMO_AEROSOL_SSH = 1
  end enum
  integer(kind=kind(CS_ATMO_AEROSOL_OFF)), pointer, save :: iaerosol

  !> Flag to desactivate gaseous chemistry
  logical(c_bool), pointer, save :: nogaseouschemistry

  !> Flag to initialize gas species with the aerosol library
  logical(c_bool), pointer, save :: init_gas_with_lib

  !> Flag to initialize aerosols with the aerosol library
  logical(c_bool), pointer, save :: init_aero_with_lib

  ! Number of aerosol layers
  integer(c_int), pointer, save :: nlayer_aer

  !> Number of aerosols
  integer(c_int), pointer, save :: n_aer

  !> Initial gaseous and particulate concentrations
  !> and aerosol number read in file
  double precision, save, allocatable, dimension(:) :: dlconc0

  !> read zone boundary conditions from file
  integer, save :: iprofa(nozppm)

  !> \}

  !=============================================================================

contains

  !-------------------------------------------------------------------------------

  !> \brief Get the aerosols concentrations and numbers from SSH-aerosol

  !> \param[out]    array         array with the aerosols
  !>                              concentrations (microg / m^3)
  !>                              and numbers (molecules / m^3)

  subroutine sshaerosol_get_aero(array)

    use, intrinsic :: iso_c_binding
    use cs_c_bindings

    implicit none

    ! Arguments

    double precision, intent(out) :: array(n_aer*(1+nlayer_aer))

    ! Local variables

    real(kind=c_double) :: c_array(n_aer*(1+nlayer_aer))

    call cs_atmo_aerosol_get_aero(c_array)

    array(1:n_aer*(1+nlayer_aer)) = c_array(1:n_aer*(1+nlayer_aer))

  end subroutine sshaerosol_get_aero

  !-----------------------------------------------------------------------------

  !> \brief Get the gas species concentrations from SSH-aerosol

  !> \param[out]    array         array with the gas species concentrations
  !>                              (microg / m^3)

  subroutine sshaerosol_get_gas(array)

    use, intrinsic :: iso_c_binding
    use cs_c_bindings
    use atchem, only : nespg

    implicit none

    ! Arguments

    double precision, intent(out) :: array(nespg)

    ! Local variables

    real(kind=c_double) :: c_array(nespg)

    call cs_atmo_aerosol_get_gas(c_array)

    array(1:nespg) = c_array(1:nespg)

  end subroutine sshaerosol_get_gas

  !=============================================================================

  !> \brief Map pointers to arrays

  subroutine init_aerosol_pointers

    use, intrinsic :: iso_c_binding
    use cs_c_bindings
    use atchem, only : nespg, isca_chem, dmmk, chempoint

    implicit none

    ! Local variables

    type(c_ptr) :: c_species_to_scalar_id, c_molar_mass, c_chempoint

    call cs_f_atmo_chem_arrays_get_pointers(c_species_to_scalar_id, &
                                            c_molar_mass, &
                                            c_chempoint)

    call c_f_pointer(c_species_to_scalar_id, isca_chem,  &
                     [nespg+(nlayer_aer+1)*n_aer])
    call c_f_pointer(c_molar_mass, dmmk, [nespg])
    call c_f_pointer(c_chempoint, chempoint, [nespg])

  end subroutine init_aerosol_pointers

  !=============================================================================

end module sshaerosol
