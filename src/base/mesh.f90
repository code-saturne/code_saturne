!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

! Module for mesh-related arrays

module mesh

  !=============================================================================

  ! Mesh Fortran structure, pointers to the C structure

  ! ndim : spatial dimension

  integer :: ndim
  parameter(ndim=3)

  ! ncelet : number of extended (real + ghost) cells
  ! ncel   : number of cells
  ! nfac   : number of interior faces
  ! nfabor : number of boundary faces
  ! nnod   : number of vertices

  integer, save :: ncelet, ncel, nfac, nfabor, nnod

  ! ncelbr : number of cells with faces on boundary

  integer, save :: ncelbr

  ! lndfac : size of nodfac indexed array
  ! lndfbr : size of nodfbr indexed array

  integer, save :: lndfac, lndfbr

  ! nfml   : number of families (group classes)
  ! nprfml : number of properties per family (group class)

  integer, save :: nprfml, nfml

  ! ifacel : interior faces -> cells connectivity
  ! ifabor : boundary faces -> cells connectivity

  integer, dimension(:,:), pointer :: ifacel
  integer, dimension(:), pointer :: ifabor

  ! ipnfac : position du premier noeud de chaque face interne dans nodfac
  ! nodfac : connectivite faces internes/noeuds
  ! ipnfbr : position du premier noeud de chaque face de bord dans nodfbr
  ! nodfbr : connectivite faces de bord/noeuds

  integer, dimension(:), pointer :: ipnfac
  integer, dimension(:), pointer :: nodfac
  integer, dimension(:), pointer :: ipnfbr
  integer, dimension(:), pointer :: nodfbr

  ! ifmfbr : boundary face family numbers
  ! ifmcel : cell family numbers
  ! iprfml : property numbers per family

  integer, dimension(:), pointer :: ifmfbr
  integer, dimension(:), pointer :: ifmcel
  integer, dimension(:,:), pointer :: iprfml

  ! icelbr : list of cells adjacent to boundary faces

  integer, dimension(:), pointer :: icelbr

  ! xyzcen : cell centers
  ! surfac : interior faces surface vectors
  ! surfbo : boundary faces surface vectors
  ! cdgfac : interior faces centers of gravity
  ! cdgfbo : boundary faces centers of gravity
  ! xyznod : vertex coordinates (optional)
  ! volume : cell volumes
  ! surfan : interior face surfaces
  ! surfbn : boundary face surfaces
  ! dist   : distance IJ.Nij
  ! distb  : likewise for border faces
  ! pond   : weighting (Aij=pond Ai+(1-pond)Aj)
  ! dijpf  : vector I'J'
  ! diipb  : likewise for border faces
  ! dofij  : vector OF at interior faces

  double precision, dimension(:,:), pointer :: xyzcen
  double precision, dimension(:,:), pointer :: surfac
  double precision, dimension(:,:), pointer :: surfbo
  double precision, dimension(:,:), pointer :: cdgfac
  double precision, dimension(:,:), pointer :: cdgfbo
  double precision, dimension(:,:), pointer :: xyznod
  double precision, dimension(:), pointer :: volume
  double precision, dimension(:), pointer :: surfan
  double precision, dimension(:), pointer :: surfbn
  double precision, dimension(:), pointer :: dist
  double precision, dimension(:), pointer :: distb
  double precision, dimension(:), pointer :: pond
  double precision, dimension(:,:), pointer :: dijpf
  double precision, dimension(:,:), pointer :: diipb
  double precision, dimension(:,:), pointer :: dofij

  !=============================================================================

end module mesh
