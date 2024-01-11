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

!> \file cplsat.f90
!> \brief Module for code/code coupling

module cplsat

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup cplsat Module for code/code coupling
  !> code / code - management of key parameters

  !> \addtogroup cplsat
  !> \{

  !> number of couplings code_saturne / code_saturne
  integer, save :: nbrcpl = 0
  !> indicator coupling face / face only
  !integer, save :: ifaccp = 0
  integer(c_int), pointer, save :: ifaccp
  !> maximum permissible number of coupling
  integer   nbcpmx
  parameter(nbcpmx=10)
  !> turbulence model of the remote instance
  integer, save :: iturcp(nbcpmx)
  !> indicator to update location of the coupling
  integer, save :: imajcp(nbcpmx)
  !> indicator of calulation in relative reference frame
  integer, save :: icormx(nbcpmx)
  !> number of variables to send/receive
  integer, save :: nvarcp(nbcpmx)
  !> size of exchange tables
  integer, save :: nvarto(nbcpmx)
  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    subroutine cs_f_sat_coupling_get_pointers(ifaccp)     &
      bind(C, name='cs_f_sat_coupling_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ifaccp
    end subroutine cs_f_sat_coupling_get_pointers

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran code_saturne coupling API.
  !> This maps Fortran pointers to global C structure members.

  subroutine cplsat_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ifaccp

    call cs_f_sat_coupling_get_pointers(c_ifaccp)

    call c_f_pointer(c_ifaccp, ifaccp)

  end subroutine cplsat_init

  !=============================================================================

end module cplsat
