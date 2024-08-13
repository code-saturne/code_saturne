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

!> \file pointe.f90
!> \brief Module for pointer variables

module pointe

  !=============================================================================

  use, intrinsic :: iso_c_binding
  use paramx

  implicit none

  !=============================================================================
  !> \defgroup pointer_variables Module for pointer variables

  !> \addtogroup pointer_variables
  !> \{

  !> \defgroup fortran_pointer_containers Containers to Fortran array pointers.
  !> An array of one of these derived types can be used to manage a set of
  !> pointers (Fortran not allow arrays of pointers directly, so this
  !> technique is a classical workaround.

  ! Note also that Fortran bounds remapping could in theory be used to
  ! handle different pointer shapes with a single type.

  !> \addtogroup fortran_pointer_containers
  !> \{

  !> container for rank 1 double precision array pointer.
  type pmapper_double_r1
    double precision, dimension(:),  pointer :: p !< rank 1 array pointer
  end type pmapper_double_r1

  !> \}

  !=============================================================================

  !> \defgroup coupled_case Specific arrays for the coupled case

  !> \addtogroup coupled_case
  !> \{

  !> \anchor itypfb
  !> boundary condition type at the boundary face \c ifac
  !> (see \ref cs_user_boundary_conditions)
  integer, dimension(:), pointer, save :: itypfb

  !> to identify boundary zones associated with boundary faces
  !> (specific physics models)
  integer, dimension(:), pointer :: izfppp

  !> \}

  !> \}

contains

  !=============================================================================

  subroutine boundary_conditions_init() &
    bind(C, name='cs_f_boundary_conditions_init')

    use, intrinsic :: iso_c_binding
    use mesh
    use cs_c_bindings

    implicit none

    ! Local variables

    type(c_ptr) :: c_itypfb, c_izfppp

    call cs_f_boundary_conditions_create

    call cs_f_boundary_conditions_get_pointers(c_itypfb, c_izfppp)

    call c_f_pointer(c_itypfb, itypfb, [nfabor])
    call c_f_pointer(c_izfppp, izfppp, [nfabor])

  end subroutine boundary_conditions_init

  !=============================================================================

end module pointe
