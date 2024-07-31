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

!> \file entsor.f90
!> \brief Module for input/output

module entsor

  !=============================================================================
  use, intrinsic :: iso_c_binding
  use paramx

  implicit none

  !=============================================================================

  !> \defgroup entsor Module for input/output

  !> \addtogroup entsor
  !> \{

  !> standard output
  integer, save :: nfecra

  !> field key for output label (\ref label field keyword).
  integer, save :: keylbl = -1

  !> field key for logging (\ref log field keyword).
  integer, save :: keylog = -1

  !> <a name="keyvis"></a>
  !> field key for postprocessing output (\ref post_vis field keyword).
  integer, save :: keyvis = -1

  !> \}

  !> \defgroup userfile Additional user files

  !> \addtogroup userfile
  !> \{

  !> name of the thermochemical data file for combustion.
  !>
  !> Useful in case of gas combustion.
  character(len=64), save :: ficrad

  !> logical unit of the thermochemical data file.
  !> Useful in case of gas combustion;
  integer, save :: impfpp = 25

  !> \}

  !=============================================================================

contains

  !=============================================================================

  !> \brief Flush Fortran log

  subroutine flush_nfecra() bind(C, name='cs_f_flush_logs')
    flush(nfecra)
  end subroutine flush_nfecra

  !=============================================================================

end module entsor
