!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file albase.f90
!> Module for Arbitrary Lagrangian Eulerian method (ALE)

module albase

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  !> \defgroup albase Module for Arbitrary Lagrangian Eulerian method (ALE)

  !> \addtogroup albase
  !> \{

  !> Activates (> 0) or not (= 0), activate the ALE module:
  !>  - 0: not activated
  !>  - 1: legacy solver
  !>  - 2: CDO solver
  integer(c_int), pointer, save :: iale
  !> the number of sub-iterations of initialization of the fluid
  integer(c_int), save :: nalinf = 0
  !> maximum number of implicitation iterations of the structure displacement
  integer(c_int), save :: nalimx = 1
  !> relative precision of implicitation of the structure displacement
  real(c_double), save :: epalim = 1.d-5
  !> iteration (yes=1, no=0) to initialize ALE
  integer(c_int), save :: italin = -999

  bind(C, name='cs_glob_mobile_structures_i_max') :: nalimx
  bind(C, name='cs_glob_mobile_structures_i_eps') :: epalim
  bind(C, name='cs_glob_ale_n_ini_f') :: nalinf
  bind(C, name='cs_glob_ale_need_init') :: italin

  !> indicator of imposed displacement
  integer(c_int), pointer :: impale(:)
  !> defines the mesh velocity from the color of the boundary faces,
  !> or more generally from their properties (colors, groups, ...),
  !> from the boundary conditions defined in cs user boundary conditions,
  !> or even from their coordinates.
  integer(c_int), pointer :: ialtyb(:)
  !> Pointer to field over vertices: mesh displacement
  integer, save :: fdiale
  !> \}

contains

  !=============================================================================

  subroutine map_ale

    use, intrinsic :: iso_c_binding

    !---------------------------------------------------------------------------

    interface

      subroutine cs_f_ale_get_pointers(iale) &
        bind(C, name='cs_f_ale_get_pointers')
        use, intrinsic :: iso_c_binding
        implicit none
        type(c_ptr), intent(out) :: iale
      end subroutine cs_f_ale_get_pointers

    end interface

    ! Local variables

    type(c_ptr) :: c_iale

    call cs_f_ale_get_pointers(c_iale)

    call c_f_pointer(c_iale, iale)

  end subroutine map_ale

  !=============================================================================

  !> \brief Map Fortran ALE boundary conditions info.

  subroutine ale_models_bc_maps

    use, intrinsic :: iso_c_binding
    use mesh, only:nfabor, nnod

    implicit none

    interface

      subroutine cs_f_ale_bc_get_pointers(p_impale, p_ale_bc_type) &
        bind(C, name='cs_f_ale_bc_get_pointers')
        use, intrinsic :: iso_c_binding
        implicit none
        type(c_ptr), intent(out) :: p_impale, p_ale_bc_type
      end subroutine cs_f_ale_bc_get_pointers

    end interface

    ! Arguments

    ! Local variables
    type(c_ptr) :: p_impale, p_ale_bc_type

    call cs_f_ale_bc_get_pointers(p_impale, p_ale_bc_type)

    call c_f_pointer(p_impale, impale, [nnod])
    call c_f_pointer(p_ale_bc_type, ialtyb, [nfabor])

  end subroutine ale_models_bc_maps

  !=============================================================================

end module albase
