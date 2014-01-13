!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

!> \file mltgrd.f90
!> \brief Module for multigrid parameters

module mltgrd

  !=============================================================================

  use paramx

  implicit none

  !> \defgroup mltgrd Module for multigrid parameters

  !> \addtogroup mltgrd
  !> \{

  !=============================================================================

  !> for the multi-grid method, maximum number of cells on the coarsest grid.
  !> useful if and only if imgr(ivar) = 1 for at least one variable ivar
  integer, save :: ncegrm

  !> when using the multi-grid method, maximum number of grid levels.
  !> useful if and only if imgr(ivar) = 1 for at least one variable ivar
  integer, save :: ngrmax


  !> mean number of cells under which merging should take place
  integer, save :: mltmmn

  !> global number of cells under which merging should take place
  integer, save :: mltmgl

  !> number of ranks over which merging takes place (stride)
  integer, save :: mltmst

  !> number of active ranks under which no merging is done
  integer, save :: mltmmr

  !> multigrid coarsening type (face traversal order)
  integer, save :: mlttyp

  !> maximum aggregation per grid level
  integer, save :: nagmx0(nvarmx)

  !> if > 0, activates post-processing output of aggregation,
  !> by projecting the coarse cell numbers (modulo ncpmgr(ivar))
  !> on the finest mesh
  integer, save :: ncpmgr(nvarmx)

  !=============================================================================
  !> \}

end module mltgrd
