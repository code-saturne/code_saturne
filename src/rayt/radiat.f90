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

  !> Period of the radiation module.
  !> The radiation module is called every \ref nfreqr time steps (more precisely,
  !> every time \ref optcal::ntcabs "ntcabs" is a multiple of \ref nfreqr).
  !> Also, in order to have proper initialization of the variables, whatever
  !> the value of \ref nfreqr, the radiation module is called at
  !> the first time step of a calculation (restart or not).
  !> Useful if and only if the radiation module is activated}
  integer(c_int), pointer, save ::           nfreqr

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
    subroutine cs_rad_transfer_get_pointers(p_iirayo, p_nfreqr)  &
      bind(C, name='cs_rad_transfer_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none

      type(c_ptr), intent(out) :: p_iirayo, p_nfreqr

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

    ! Interface to C function handling source terms

    subroutine cs_rad_transfer_bcs(nvar, bc_type, icodcl, dt, rcodcl)  &
      bind(C, name='cs_rad_transfer_bcs')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: nvar
      integer(kind=c_int), dimension(*) :: bc_type, icodcl
      real(kind=c_double), dimension(*) :: dt, rcodcl
    end subroutine cs_rad_transfer_bcs

    !---------------------------------------------------------------------------

    ! Interface to C function handling resolution

    subroutine cs_rad_transfer_solve(verbosity, bc_type,     &
                                     cp2fol, cp2ch, ichcor)  &
      bind(C, name='cs_rad_transfer_solve')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: verbosity
      integer(kind=c_int), dimension(*) :: bc_type, ichcor
      real(kind=c_double), value :: cp2fol
      real(kind=c_double), dimension(*) :: cp2ch
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

    type(c_ptr) :: p_iirayo, p_nfreqr

    call cs_rad_transfer_get_pointers(p_iirayo, p_nfreqr)

    call c_f_pointer(p_iirayo, iirayo)
    call c_f_pointer(p_nfreqr, nfreqr)

  end subroutine radiat_init

  !=============================================================================

  ! Free related arrays

  subroutine radiat_finalize

    call cs_rad_transfer_finalize

  end subroutine radiat_finalize

  !=============================================================================

end module radiat
