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

!> \file parall.f90
!> \brief Module for basic MPI and OpenMP parallelism-related values

module parall

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup parall Module for basic MPI and OpenMP parallelism-related values

  !> \addtogroup parall
  !> \{

  !> thr_n_min : minimum number of elements for loops on threads
  integer   thr_n_min
  parameter(thr_n_min = 128)

  !> process rank
  !> - -1 in sequential mode
  !> - r (0 < r < n_processes) in distributed parallel run
  integer, save ::  irangp

  !> number of processes (=1 if sequental)
  integer, save ::  nrangp

  ! Global dimensions (i.e. independent of parallel partitioning)

  !> global number of cells
  integer(kind=8), save :: ncelgb
  !> global number of interior faces
  integer(kind=8), save :: nfacgb
  !> global number of boundary faces
  integer(kind=8), save :: nfbrgb
  !> global number of vertices
  integer(kind=8), save :: nsomgb

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \brief Compute the global maximum of an integer in case of parellism.

    !> \param[in, out]   max  local max in, global max out

    subroutine parcmx(count)  &
      bind(C, name='cs_f_parall_max_i')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(inout) :: count
    end subroutine parcmx

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

    !> \brief Compute the global minimum of an integer in case of parellism.

    !> \param[in, out]   min  local min in, global min out

    subroutine parcmn(count)  &
      bind(C, name='cs_f_parall_min_i')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(inout) :: count
    end subroutine parcmn

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

    !> \brief Compute the global maxima of an array of integers
    !> in case of parellism.

    !> \param[in]        n_elts  size of array
    !> \param[in, out]   max  local max in, global max out

    subroutine parimx(n_elts, array)  &
      bind(C, name='cs_f_parall_max_n_i')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n_elts
      integer(c_int), dimension(*), intent(inout) :: array
    end subroutine parimx

    !---------------------------------------------------------------------------

    !> \brief Compute the global maxima of an array of real numbers
    !> in case of parellism.

    !> \param[in]        n_elts  size of array
    !> \param[in, out]   max     local max in, global max out

    subroutine parrmx(n_elts, array)  &
      bind(C, name='cs_f_parall_max_n_r')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n_elts
      real(c_double), dimension(*), intent(inout) :: array
    end subroutine parrmx

    !---------------------------------------------------------------------------

    !> \brief Compute the global minima of an array of integers
    !> in case of parellism.

    !> \param[in]        n_elts  size of array
    !> \param[in, out]   min  local min in, global min out

    subroutine parimn(n_elts, array)  &
      bind(C, name='cs_f_parall_min_n_i')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n_elts
      integer(c_int), dimension(*), intent(inout) :: array
    end subroutine parimn

    !---------------------------------------------------------------------------

    !> \brief Compute the global minima of an array of real numbers
    !> in case of parellism.

    !> \param[in]        n_elts  size of array
    !> \param[in, out]   min     local min in, global min out

    subroutine parrmn(n_elts, array)  &
      bind(C, name='cs_f_parall_min_n_r')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n_elts
      real(c_double), dimension(*), intent(inout) :: array
    end subroutine parrmn

    !---------------------------------------------------------------------------

    !> \brief Compute the global sums of an array of integers
    !> in case of parellism.

    !> Note that for counters, on very large meshes, if a sum exceeds
    !> 2**31, the resuly will be false on most machines. To avoid this,
    !> using the C API (with counters as cs_gnum_t) is preferred.

    !> \param[in]        n_elts  size of array
    !> \param[in, out]   sum  local sum in, global sum out

    subroutine parism(n_elts, array)  &
      bind(C, name='cs_f_parall_sum_n_i')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n_elts
      integer(c_int), dimension(*), intent(inout) :: array
    end subroutine parism

    !---------------------------------------------------------------------------

    !> \brief Compute the global sums of an array of real numbers
    !> in case of parellism.

    !> \param[in]        n_elts  size of array
    !> \param[in, out]   sum  local sum in, global sum out

    subroutine parrsm(n_elts, array)  &
      bind(C, name='cs_f_parall_sum_n_r')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n_elts
      real(c_double), dimension(*), intent(inout) :: array
    end subroutine parrsm

    !---------------------------------------------------------------------------

    !> \brief Broadcast an integer in case of parellism.

    !> \param[in]        root_rank  rank of the sending process
    !> \param[in, out]   val        value to broadcast
    !>                              (input on root_rank, output on others)

    subroutine parall_bcast_i(root_rank, val)  &
      bind(C, name='cs_f_parall_bcast_i')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: root_rank
      integer(c_int), intent(inout) :: val
    end subroutine parall_bcast_i

    !---------------------------------------------------------------------------

    !> \brief Broadcast a real number in case of parellism.

    !> \param[in]        root_rank  rank of the sending process
    !> \param[in, out]   val        value to broadcast
    !>                              (input on root_rank, output on others)

    subroutine parall_bcast_r(root_rank, val)  &
      bind(C, name='cs_f_parall_bcast_r')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: root_rank
      real(c_double), intent(inout) :: val
    end subroutine parall_bcast_r

    !---------------------------------------------------------------------------

    !> \brief Broadcast an array of integers in case of parellism.

    !> \param[in]        root_rank  rank of the sending process
    !> \param[in]        n_elts     size of array
    !> \param[in, out]   array      array to broadcast
    !>                              (input on root_rank, output on others)

    subroutine parbci(root_rank, n_elts, array)  &
      bind(C, name='cs_f_parall_bcast_n_i')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: root_rank, n_elts
      integer(c_int), dimension(*), intent(inout) :: array
    end subroutine parbci

    !---------------------------------------------------------------------------

    !> \brief Broadcast an array of real numbers in case of parellism.

    !> \param[in]        root_rank  rank of the sending process
    !> \param[in]        n_elts     size of array
    !> \param[in, out]   array      array to broadcast
    !>                              (input on root_rank, output on others)

    subroutine parbcr(root_rank, n_elts, array)  &
      bind(C, name='cs_f_parall_bcast_n_r')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: root_rank, n_elts
      real(c_double), dimension(*), intent(inout) :: array
    end subroutine parbcr

    !---------------------------------------------------------------------------

    !> \brief Maximum value of a real and the value of related array on all
    !> default communicator processes.

    !> \param[in]       n             size of the related array
    !> \param[in, out]  max           local max in, global max out
    !> \param[in, out]  max_loc_vals  array values at location of local max in,
    !>                                and at location of global max out

    subroutine parmxl(n, max, max_loc_vals)  &
      bind(C, name='cs_parall_max_loc_vals')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n
      real(c_double), intent(inout) :: max
      real(c_double), dimension(*), intent(inout) :: max_loc_vals
    end subroutine parmxl

    !---------------------------------------------------------------------------

    !> \brief Minimum value of a real and the value of related array on all
    !> default communicator processes.

    !> \param[in]       n             size of the related array
    !> \param[in, out]  min           local min in, global min out
    !> \param[in, out]  min_loc_vals  array values at location of local min in,
    !>                                and at location of global min out

    subroutine parmnl(n, min, min_loc_vals)  &
      bind(C, name='cs_parall_min_loc_vals')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n
      real(c_double), intent(inout) :: min
      real(c_double), dimension(*), intent(inout) :: min_loc_vals
    end subroutine parmnl

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

    !> \brief Build a global array from each local array in each domain.

    !> Local arrays are appened in order of owning MPI rank.
    !> The size of each local array may be different.

    !> Use of this function may be quite practical, but should be limited
    !> to user functions, as it may limit scalability (especially as regards
    !> memory usage).

    !> \param[in]   n_elts    size of the local array
    !> \param[in]   n_g_elts  size of the global array
    !> \param[in]   array     local array (size: n_elts)
    !> \param[out]  g_array   global array  (size: n_g_elts)

    subroutine cs_parall_allgather_r(n_elts, n_g_elts, array, g_array)  &
      bind(C, name='cs_parall_allgather_r')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n_elts, n_g_elts
      real(c_double), dimension(*), intent(in) :: array
      real(c_double), dimension(*), intent(inout) :: g_array
    end subroutine cs_parall_allgather_r

    !---------------------------------------------------------------------------

    !> \brief Set a barrier on all default communicator processes.

    !> The this function will exit only once all processes have called
    !> it. This is not recommended for production code, but may be useful
    !> when debugging.

    subroutine parbar()  &
      bind(C, name='cs_f_parall_barrier')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine parbar

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Build a global array from each local array in each domain.

  !> Local arrays are appened in order of owning MPI rank.
  !> The size of each local array may be different.

  !> Use of this function may be quite practical, but should be limited
  !> to user functions, as it may limit scalability (especially as regards
  !> memory usage).

  !> \param[in]   n_elts    size of the local array
  !> \param[in]   n_g_elts  size of the global array
  !> \param[in]   array     local array (size: n_elts)
  !> \param[out]  g_array   global array  (size: n_g_elts)

  subroutine paragv(n_elts, n_g_elts, array, g_array)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: n_elts, n_g_elts
    real(c_double), dimension(:), intent(in) :: array
    real(c_double), dimension(:), intent(inout) :: g_array
    call cs_parall_allgather_r(n_elts, n_g_elts, array(:), g_array(:))
  end subroutine paragv

  !=============================================================================

end module parall


