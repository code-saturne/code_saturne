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

! Module for cooling towers

module ctincl

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use cstphy

  implicit none

  !=============================================================================

  !> \defgroup ctincl Module for cooling towers constants

  !> \addtogroup ctincl
  !> \{

  !=============================================================================

  !> Initial absolute humidity in the cooling tower
  real(c_double), pointer, save :: humidity0

  !> Cp of dry air
  real(c_double), pointer, save :: cp_a

  !> Cp of water vapour
  real(c_double), pointer, save :: cp_v

  !> Cp of liquid water
  real(c_double), pointer, save :: cp_l

  !> Enthalpy of vapourisation of water
  real(c_double), pointer, save :: hv0

  !> Density of liquid water
  real(c_double), pointer, save :: rho_l

  !> Conductivity of liquid water
  real(c_double), pointer, save :: lambda_l

  !> Conductivity of humid air
  real(c_double), pointer, save :: lambda_h

  !> Droplet diameter for rain zones and liquid water initialization
  real(c_double), pointer, save :: droplet_diam

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS


    ! Interface to C function retrieving pointers to members of the
    ! global fluid properties structure
    subroutine cs_ctwr_glob_properties_get_pointer( &
        humidity0, &
        cp_a, cp_v, cp_l, hv0, rho_l, &
        lambda_h, lambda_l, &
        droplet_diam) &
        bind(C, name='cs_ctwr_glob_properties_get_pointer')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: humidity0
      type(c_ptr), intent(out) :: cp_a, cp_v, cp_l, hv0, rho_l
      type(c_ptr), intent(out) :: lambda_h, lambda_l
      type(c_ptr), intent(out) :: droplet_diam
    end subroutine cs_ctwr_glob_properties_get_pointer

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

contains

  subroutine ctwr_properties_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_humidity0
    type(c_ptr) :: c_cp_a, c_cp_v, c_cp_l, c_hv0, c_rho_l
    type(c_ptr) :: c_lambda_h, c_lambda_l
    type(c_ptr) :: c_droplet_diam

    call cs_ctwr_glob_properties_get_pointer( &
      c_humidity0, &
      c_cp_a, c_cp_v, c_cp_l, c_hv0, c_rho_l, &
      c_lambda_h, c_lambda_l, &
      c_droplet_diam)

    call c_f_pointer(c_humidity0   , humidity0   )
    call c_f_pointer(c_cp_a        , cp_a        )
    call c_f_pointer(c_cp_v        , cp_v        )
    call c_f_pointer(c_cp_l        , cp_l        )
    call c_f_pointer(c_hv0         , hv0         )
    call c_f_pointer(c_rho_l       , rho_l       )
    call c_f_pointer(c_lambda_h    , lambda_h    )
    call c_f_pointer(c_lambda_l    , lambda_l    )
    call c_f_pointer(c_droplet_diam, droplet_diam)

  end subroutine ctwr_properties_init

  !=============================================================================

end module ctincl
