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

!> \file mesh.f90
!> Module for mesh-related arrays

module mesh

  !=============================================================================

  implicit none

  !=============================================================================

  ! Mesh Fortran structure, pointers to the C structure

  ! ndim : spatial dimension

  integer :: ndim
  parameter(ndim=3)

  integer, save :: ncelet  !< number of extended (real + ghost) cells
  integer, save :: ncel    !< number of cells
  integer, save :: nfac    !< number of interior faces
  integer, save :: nfabor  !< number of boundary faces
  integer, save :: nnod    !< number of vertices

  integer, save :: ncelbr  !< number of cells with faces on boundary

  integer, save :: lndfac  !< size of nodfac indexed array
  integer, save :: lndfbr  !< size of nodfbr indexed array

  integer, save :: nfml    !< number of families (group classes)

  !> interior faces -> cells connectivity
  integer, dimension(:,:), pointer :: ifacel

  !> boundary  faces -> cells connectivity
  integer, dimension(:), pointer :: ifabor

  !> interior face -> vertex index
  integer, dimension(:), pointer :: ipnfac

  !> interior face -> vertex connectivity
  integer, dimension(:), pointer :: nodfac

  !> boundary face -> vertex index
  integer, dimension(:), pointer :: ipnfbr

  !> boundary face -> vertex connectivity
  integer, dimension(:), pointer :: nodfbr

  integer, dimension(:), pointer :: ifmfbr    !< boundary face family numbers
  integer, dimension(:), pointer :: ifmcel    !< cell family numbers

  !> list of cells adjacent to boundary faces
  integer, dimension(:), pointer :: icelbr

  !> symmetry marker (0 to cancel mass flux)
  integer, dimension(:), pointer :: isympa

  !> cell centers
  double precision, dimension(:,:), pointer :: xyzcen

  !> interior face surface vectors
  double precision, dimension(:,:), pointer :: surfac

  !> boundary face surface vectors
  double precision, dimension(:,:), pointer :: surfbo

  !> interior face centers of gravity
  double precision, dimension(:,:), pointer :: cdgfac

  !> boundary face centers of gravity
  double precision, dimension(:,:), pointer :: cdgfbo

  !> vertex coordinates
  double precision, dimension(:,:), pointer :: xyznod

  !> cell volumes
  double precision, dimension(:), pointer :: volume

  !> interior face surfaces
  double precision, dimension(:), pointer :: surfan

  !> boundary face surfaces
  double precision, dimension(:), pointer :: surfbn

  !> distance IJ.Nij
  double precision, dimension(:), pointer :: dist

  !> distance IF.N for boundary faces
  double precision, dimension(:), pointer :: distb

  !> weighting (Aij=pond Ai+(1-pond)Aj)
  double precision, dimension(:), pointer :: pond

  !> vector I'J' for interior faces
  double precision, dimension(:,:), pointer :: dijpf

  !> vector II' for interior faces
  double precision, dimension(:,:), pointer :: diipb

  !> vector OF for interior faces
  double precision, dimension(:,:), pointer :: dofij

  !=============================================================================

end module mesh
