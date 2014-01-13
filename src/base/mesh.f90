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

!> \file mesh.f90
!> \brief Module for mesh-related arrays

module mesh

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup mesh Mesh Fortran structure, pointers to the C structure

  !> \addtogroup mesh
  !> \{

  !< spatial dimension (3)
  integer :: ndim
  parameter(ndim=3)

  !< number of extended (real + ghost of the 'halo') cells. See \ref note_1
  integer, save :: ncelet

  !< number of real cells in the mesh
  integer, save :: ncel

  !< number of internal faces  (see \ref note_2)
  integer, save :: nfac

  !< number of boundary faces (see \ref note_2)
  integer, save :: nfabor

  !< number of vertices in the mesh
  integer, save :: nnod

  !< number of cells with at least one boundary
  integer, save :: ncelbr

  !< size of the array \c nodfac of internal faces - nodes connectivity
  !< (see \ref note_3)
  integer, save :: lndfac

  !< size of the array \c nodfbr of boundary faces - nodes connectivity
  !< (see \ref note_3)
  integer, save :: lndfbr

  !< Number of referenced families of entities (boundary faces, elements, ...)
  integer, save :: nfml

  !< Index-numbers of the two (only) neighbouring cells for each internal face
  integer, dimension(:,:), pointer :: ifacel

  !> index-number of the (unique) neighbouring cell for each boundary face
  integer, dimension(:), pointer :: ifabor

  !> position of the first node of the each internal face in the array nodfac
  !> (see \ref note_3)
  integer, dimension(:), pointer :: ipnfac

  !> index-numbers of the nodes of each internal face
  !> (see \ref note_3)
  integer, dimension(:), pointer :: nodfac

  !> position of the first node of the each boundary face in the array nodfbr
  !> (see \ref note_3)
  integer, dimension(:), pointer :: ipnfbr

  !> index-numbers of the nodes of each boundary face
  !> (see \ref note_3)
  integer, dimension(:), pointer :: nodfbr

  !> family number of the boundary faces. See \ref note_1
  integer, dimension(:), pointer :: ifmfbr

  !> family number of the elements. See \ref note_1
  integer, dimension(:), pointer :: ifmcel

  !> list of cells having at least one boundary face
  integer, dimension(:), pointer :: icelbr

  !> integer to mark out the "symmetry" (itypfb=isymet) boundary faces
  !> where the mass flow has to be canceled when the ALE module is switched
  !> off (these faces are impermeable).
  !> For instance, if the face ifac is symmetry face,
  !> isympa(ifac)=0, otherwise isympa(ifac)=1.
  integer, dimension(:), pointer :: isympa

  !> coordinate of the cell centers
  double precision, dimension(:,:), pointer :: xyzcen

  !> surface vector of the internal faces. Its norm is the surface of the face
  !> and it is oriented from \c ifacel(1,.) to \c ifacel(2,.)
  double precision, dimension(:,:), pointer :: surfac


  !> surface vector of the boundary faces. Its norm is the surface of the face
  !> and it is oriented outwards
  double precision, dimension(:,:), pointer :: surfbo

  !> coordinates of the centres of the internal faces
  double precision, dimension(:,:), pointer :: cdgfac

  !> coordinates of the centres of the boundary faces
  double precision, dimension(:,:), pointer :: cdgfbo

  !> coordinates of the mesh vertices
  double precision, dimension(:,:), pointer :: xyznod

  !> volume of each cell
  double precision, dimension(:), pointer :: volume

  !> norm of the surface vector of the internal faces
  double precision, dimension(:), pointer :: surfan

  !> norm of the surface of the boundary faces
  double precision, dimension(:), pointer :: surfbn

  !> distance IJ.Nij
  !> for every internal face, dot product of the vectors
  !> \f$ \vect{IJ}\f$ and \f$\vect{n}\f$.  I and J are respectively
  !> the centres of the first and the second neighbouring cell.
  !> The vector \f$\vect{n}\f$ is the unit vector normal to the face
  !> and oriented from the first to the second cell
  double precision, dimension(:), pointer :: dist

  !> distance IF.N for boundary faces
  !> For every boundary face, dot product between the vectors
  !> \f$\vect{IF}\f$ and \f$\vect{n}\f$.
  !> I is the center of the neighbouring cell. F is the face center.
  !> The vector \f$\vect{n}\f$ is the unit vector normal to the face and
  !> oriented to the exterior of the domain
  double precision, dimension(:), pointer :: distb

  !> weighting (Aij=pond Ai+(1-pond)Aj)
  !> for every internal face,
  !> \f$\displaystyle\frac{\vect{FJ}.\vect{n}}{\vect{IJ}.\vect{n}}\f$.
  !> With regard to the mesh quality, its ideal value is 0.5
  double precision, dimension(:), pointer :: pond

  !> vector I'J' for interior faces
  !> for every internal face, the three components of the vector
  !> \f$\vect{I'J'}\f$, where I' and J' are
  !> respectively the orthogonal projections of the neighbouring cell
  !> centres I and J on a straight line orthogonal to the face and passing
  !> through its center
  double precision, dimension(:,:), pointer :: dijpf

  !> vector II' for interior faces
  !> for every boundary face, the three components of the vector
  !> \f$\vect{II'}\f$. I' is the orthogonal projection of I,
  !> center of the neighbouring cell, on the
  !> straight line perpendicular to the face and passign through its center
  double precision, dimension(:,:), pointer :: diipb

  !> vector OF for interior faces
  !> for every internal face, the three components of the vector
  !> \f$\vect{OF}\f$. O is the intersection
  !> point between the face and the straight line joining the centres
  !> of the two neighbouring cells. F is the face center
  double precision, dimension(:,:), pointer :: dofij

  !=============================================================================

  !> \}

end module mesh
