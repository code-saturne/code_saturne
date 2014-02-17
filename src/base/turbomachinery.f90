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

  integer, dimension(:), pointer :: irotce

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

  ! Initialization of turbomachinery module variables

  subroutine turbomachinery_init

    use, intrinsic :: iso_c_binding
    use mesh
    use cstphy
    use cplsat, only: imobil

    implicit none

    ! Local variables

    integer iel
    type(c_ptr) :: c_p

    ! Map turbomachinery module components to global c turbomachinery structure

    call map_turbomachinery_module(iturbo, rotax, c_p)

    call c_f_pointer(c_p, irotce, [ncelet])

    ! In case of relative frame computation or turbomachinery computation
    ! with code/code coupling, the data management is merged with the one
    ! of turbomachinery module

    if (iturbo.eq.0) then
      if (icorio.eq.1.or.imobil.eq.1) then

        rotax(1) = omegax
        rotax(2) = omegay
        rotax(3) = omegaz

        allocate(irotce(ncelet))
        do iel = 1, ncelet
          irotce(iel) = 1
        enddo
      endif
    endif

    return

  end subroutine turbomachinery_init

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

  ! Finalization of turbomachinery module variables

  subroutine turbomachinery_finalize

    use cstphy, only: icorio
    use cplsat, only: imobil

    if (iturbo.eq.0) then
      if (icorio.eq.1.or.imobil.eq.1)  deallocate(irotce)
    endif

    return

  end subroutine turbomachinery_finalize

  !=============================================================================

end module turbomachinery
