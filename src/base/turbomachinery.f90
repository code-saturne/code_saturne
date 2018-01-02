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

!> \file turbomachinery.f90
!> \brief Module for turbomachinery computations

module turbomachinery

  !=============================================================================
  !> \defgroup at_turbomachinery Module for turbomachinery computations

  !> \addtogroup at_turbomachinery
  !> \{

  !> Type of turbomachinery computation :
  !>   none (0), frozen rotor (1), transient (2)

  integer, save :: iturbo

  !> Rotor identifier list (1:ncel)

  integer, dimension(:), pointer :: irotce

  ! Elapsed time for logging

  double precision, save :: rs_ell(2)

  ! Arrays associated to wall BC update

  double precision, dimension(:), allocatable :: coftur, hfltur

  !> \}

  !=============================================================================

  interface

    ! Interface to C function mapping some data for turbomachinery

    subroutine map_turbomachinery_module (iturbo2, irotce2) &
      bind(C, name='cs_f_map_turbomachinery_module')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(out) :: iturbo2
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

    call map_turbomachinery_module(iturbo, c_p)

    call c_f_pointer(c_p, irotce, [ncelet])

    ! In case of relative frame computation or turbomachinery computation
    ! with code/code coupling, the data management is merged with the one
    ! of turbomachinery module

    if (iturbo.eq.0) then
      if (icorio.eq.1.or.imobil.eq.1) then
        allocate(irotce(ncelet))
        do iel = 1, ncelet
          irotce(iel) = 1
        enddo
      endif
    endif

    ! Initialization of elapsed times

    rs_ell(1) = 0.d0
    rs_ell(2) = 0.d0

    ! Allocate arrays for wall velocity BC update

    if (iturbo.eq.2)  allocate(coftur(nfabor), hfltur(nfabor))

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

    call map_turbomachinery_module(iturbo2, c_p)

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

    if (iturbo.eq.2)  deallocate(coftur, hfltur)

    return

  end subroutine turbomachinery_finalize

  !=============================================================================

end module turbomachinery
