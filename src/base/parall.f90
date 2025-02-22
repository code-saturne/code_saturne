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

!> \file parall.f90
!> \brief Module for basic MPI and OpenMP parallelism-related values

module parall

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup parall Module for basic MPI and OpenMP parallelism-related values

  !> \addtogroup parall
  !> \{

  !> process rank
  !> - -1 in sequential mode
  !> - r (0 < r < n_processes) in distributed parallel run
  integer, save ::  irangp = -1

  !> number of processes (=1 if sequental)
  integer, save ::  nrangp = 1

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \brief Compute the global maximum of a real number in case of parellism.

    !> \param[in, out]   max  local max in, global max out

    subroutine parmax(max)  &
      bind(C, name='cs_f_parall_max_r')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), intent(inout) :: max
    end subroutine parmax

    !---------------------------------------------------------------------------

    !> \brief Compute the global minimum of a real number in case of parellism.

    !> \param[in, out]   min  local min in, global min out

    subroutine parmin(min)  &
      bind(C, name='cs_f_parall_min_r')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), intent(inout) :: min
    end subroutine parmin

    !---------------------------------------------------------------------------

    !> \brief Compute the global sum of an integer in case of parellism.

    !> Note that for counters, on very large meshes, if the sum exceeds
    !> 2**31, the result will be false on most machines. To avoid this,
    !> using the C API (with counters as cs_gnum_t) is preferred.

    !> \param[in, out]   sum  local sum in, global sum out

    subroutine parcpt(count)  &
      bind(C, name='cs_f_parall_sum_i')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(inout) :: count
    end subroutine parcpt

    !---------------------------------------------------------------------------

    !> \brief Compute the global sum of a real number in case of parellism.

    !> \param[in, out]   sum  local sum in, global sum out

    subroutine parsom(sum)  &
      bind(C, name='cs_f_parall_sum_r')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), intent(inout) :: sum
    end subroutine parsom

    !---------------------------------------------------------------------------

    !> \brief Given an (id, rank, value) tuple, return the local id, rank,
    !>        and value corresponding to the global minimum value.

    !> \param[in, out]   elt_id   element id for which the value is the smallest
    !>                            (local in, global out)
    !> \param[in, out]   rank_id  rank id for which the value is the smallest
    !>                            (local in, global out)
    !> \param[in]        val      associated local minimum value

    subroutine parfpt(elt_id, rank_id, val)  &
      bind(C, name='cs_parall_min_id_rank_r')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(inout) :: elt_id, rank_id
      real(c_double), value :: val
    end subroutine parfpt

    !---------------------------------------------------------------------------

  end interface

end module parall


