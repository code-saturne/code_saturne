!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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
!> Module for multigrid parameters

module mltgrd

  !=============================================================================

  use paramx

  implicit none

  !=============================================================================

  ! Multigrid
  ! -----------
  !   ncegrm : maximum number of cells on coarsest mesh
  !   ngrmax : maximum number of grids

  !   mltmmn : mean number of cells under which merging should take place
  !   mltmgl : global number of cells under which merging should take place
  !   mltmmr : number of active ranks under which no merging is done
  !   mltmst : number of ranks over which merging takes place (stride)
  !   mlttyp : multigrid coarsening type (face traversal order)
  !
  !   nagmx0 : maximum aggregation per grid level
  !   ncpmgr : if > 0, activates post-processing output of aggregation,
  !            by projecting the coarse cell numbers (modulo ncpmgr(ivar))
  !            on the finest mesh

  integer, save :: ncegrm, ngrmax
  integer, save :: mltmmn, mltmgl, mltmst, mltmmr, mlttyp
  integer, save :: nagmx0(nvarmx), ncpmgr(nvarmx)

  !=============================================================================

end module mltgrd
