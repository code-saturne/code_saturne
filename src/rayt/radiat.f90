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

!> \file radiat.f90
!> Module for Radiation

module radiat

  !===========================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !> \defgroup radiat Module for Radiative transfer

  !> \addtogroup radiat
  !> \{

  !===========================================================================

  !> Activation of the radiative transfer module:
  !>  - 0: not activated
  !>  - 1: DOM
  !>  - 2: P1
  integer(c_int), pointer, save :: iirayo

  !>Spectral radiation models
  !>Number of ETRs to solve.
  integer(c_int), pointer, save :: nwsgg

  !> Time step interval for radiative properties updating from the library
  integer(c_int), pointer, save :: nt_rad_prp

  !> Atmospheric radiation model:
  !> - Direct Solar the first bit (H2O band)
  !> - Direct Solar O3 band for the second bit
  !> - diFfuse Solar for the third bit (SIR H2O band)
  !> - diFfuse Solar O3 for the fourth bit (SUV O3 band)
  !> - InfraRed for the fifth bitPeriod of the radiation module.
  integer(c_int), pointer, save :: rad_atmo_model

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    ! Interface to C function finalizing quadrature
    subroutine cs_rad_transfer_finalize() &
      bind(C, name='cs_rad_transfer_finalize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_rad_transfer_finalize

    !---------------------------------------------------------------------------

    ! Interface to C function to retrieve pointers
    subroutine cs_rad_transfer_get_pointers(p_iirayo,     p_nwsgg,             &
                                            p_nt_rad_prp, p_rad_atmo_model)    &
      bind(C, name='cs_rad_transfer_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none

      type(c_ptr), intent(out) :: p_iirayo, p_nwsgg
      type(c_ptr), intent(out) :: p_nt_rad_prp, p_rad_atmo_model

    end subroutine cs_rad_transfer_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function for restart

    subroutine cs_rad_transfer_read() &
      bind(C, name='cs_rad_transfer_read')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_rad_transfer_read

    !---------------------------------------------------------------------------

    ! Interface to C function for checkpoint

    subroutine cs_rad_transfer_write() &
      bind(C, name='cs_rad_transfer_write')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_rad_transfer_write

    !---------------------------------------------------------------------------

    ! Interface to C function defining options

    subroutine cs_rad_transfer_options()                               &
      bind(C, name='cs_rad_transfer_options')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_rad_transfer_options

    !---------------------------------------------------------------------------

    ! Interface to C function handling resolution

    subroutine cs_rad_transfer_solve(bc_type)  &
      bind(C, name='cs_rad_transfer_solve')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), dimension(*) :: bc_type
    end subroutine cs_rad_transfer_solve

    !---------------------------------------------------------------------------

    ! Interface to C function handling source terms

    subroutine cs_rad_transfer_source_terms(smbrs, rovsdt)   &
      bind(C, name='cs_rad_transfer_source_terms')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*) :: smbrs, rovsdt
    end subroutine cs_rad_transfer_source_terms

    !---------------------------------------------------------------------------

    ! Interface to C function

    function cs_rad_time_is_active() result(is_active)  &
      bind(C, name='cs_rad_time_is_active')
      use, intrinsic :: iso_c_binding
      implicit none
      logical(kind=c_bool) :: is_active
    end function cs_rad_time_is_active

    !---------------------------------------------------------------------------

  end interface

contains

  !=============================================================================

  ! get C pointers

  subroutine radiat_init

    use ppppar
    use ppincl
    use optcal
    use ppcpfu
    use numvar

    type(c_ptr) :: p_iirayo, p_nwsgg, p_nt_rad_prp, p_rad_atmo_model

    call cs_rad_transfer_get_pointers(p_iirayo,     p_nwsgg,                   &
                                      p_nt_rad_prp, p_rad_atmo_model)

    call c_f_pointer(p_iirayo, iirayo)
    call c_f_pointer(p_nwsgg, nwsgg)
    call c_f_pointer(p_nt_rad_prp, nt_rad_prp)
    call c_f_pointer(p_rad_atmo_model, rad_atmo_model)

  end subroutine radiat_init

  !=============================================================================

  ! Free related arrays

  subroutine radiat_finalize

    call cs_rad_transfer_finalize

  end subroutine radiat_finalize

  !=============================================================================

end module radiat
