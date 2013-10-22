!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

! Module for turbomachinery capabilities

module turbomachinery

  !=============================================================================

  ! Type of turbomachinery computation :
  !   none (0), frozen rotor (1), transient (2)

  integer, save :: iturbo

  ! Rotation axis

  double precision, save :: rotax(3)

  ! Flag on cells related to the rotor

  logical, dimension(:), pointer :: irotce

  !=============================================================================

  interface

    ! Interface to C function mapping some data for turbomachinery

    subroutine map_turbomachinery_module (iturbo2, rotax2, irotce2) &
      bind(C, name='cs_f_map_turbomachinery_module')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(out) :: iturbo2
      real(c_double), dimension(1:3), intent(out) :: rotax2
      type(c_ptr), intent(out) :: irotce2
    end subroutine map_turbomachinery_module

    !---------------------------------------------------------------------------

    ! Interface to C function updating mesh for turbomachinery

    subroutine turbomachinery_update_mesh(t_cur_mob, t_elapsed)            &
      bind(C, name='cs_turbomachinery_update_mesh')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), value :: t_cur_mob
      real(kind=c_double), intent(out) :: t_elapsed
    end subroutine turbomachinery_update_mesh

    !---------------------------------------------------------------------------

    ! Interface to C function resetting face fields for turbomachinery

    subroutine turbomachinery_reinit_i_face_fields() &
      bind(C, name='cs_turbomachinery_reinit_i_face_fields')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine turbomachinery_reinit_i_face_fields

    !---------------------------------------------------------------------------

    ! Interface to C function resizing owned cell fields for turbomachinery

    subroutine turbomachinery_resize_cell_fields() &
      bind(C, name='cs_turbomachinery_resize_cell_fields')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine turbomachinery_resize_cell_fields

  end interface

contains

  !=============================================================================

  ! Map turbomachinery module components to global c turbomachinery structure

  subroutine turbomachinery_init_mapping

    use, intrinsic :: iso_c_binding
    use mesh

    implicit none

    ! Local variables

    type(c_ptr) :: c_p

    ! Synchronize the list of rotating cells

    call map_turbomachinery_module(iturbo, rotax, c_p)

    call c_f_pointer(c_p, irotce, [ncelet])

    return

  end subroutine turbomachinery_init_mapping

  !=============================================================================

  ! Sync turbomachinery module components to global c turbomachinery structure

  subroutine turbomachinery_update

    use, intrinsic :: iso_c_binding
    use mesh

    implicit none

    ! Local variables

    type(c_ptr) :: c_p
    integer(c_int) :: iturbo2
    real(c_double), dimension(1:3) :: rotax2

    call map_turbomachinery_module(iturbo2, rotax2, c_p)

    call turbomachinery_resize_cell_fields

    call c_f_pointer(c_p, irotce, [ncelet])

    return

  end subroutine turbomachinery_update

  !=============================================================================

end module turbomachinery
