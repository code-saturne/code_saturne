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

!> \file dimens.f90
!> Module for dimensions

module dimens

  !=============================================================================

  ! Mesh and field data

  !=============================================================================

  integer, save :: nvar, nscal, nvisls

  integer, save :: ncofab

  integer, save :: nproce, nprofb

  ! Fake dimension for arrays propfb, coefa and coefb
  ! where nfabor = 0 (to avoid issues with array bounds when
  ! multidimensional arrays have size nfabor in one dimension)

  integer, save :: ndimfb

  !=============================================================================

end module dimens
