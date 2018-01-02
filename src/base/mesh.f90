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

!> \file mesh.f90
!> \brief Module for mesh-related arrays

module mesh

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup mesh Mesh Fortran structure, pointers to the C structure

  !> \addtogroup mesh
  !> \{

  !> \anchor ndim
  !> spatial dimension (3)
  integer :: ndim
  parameter(ndim=3)

  !> \anchor ncelet
  !> number of extended (real + ghost of the 'halo') cells. See \ref note_1
  integer, save :: ncelet

  !> \anchor ncel
  !> number of real cells in the mesh
  integer, save :: ncel

  !> \anchor nfac
  !> number of internal faces  (see \ref note_2)
  integer, save :: nfac

  !> \anchor nfabor
  !> number of boundary faces (see \ref note_2)
  integer, save :: nfabor

  !> \anchor nnod
  !> number of vertices in the mesh
  integer, save :: nnod

  !> \anchor lndfac
  !> size of the array \c nodfac of internal faces - nodes connectivity
  !> (see \ref note_3)
  integer, save :: lndfac

  !> \anchor lndfbr
  !> size of the array \c nodfbr of boundary faces - nodes connectivity
  !> (see \ref note_3)
  integer, save :: lndfbr

  !> \anchor nfml
  !> Number of referenced families of entities (boundary faces, elements, ...)
  integer, save :: nfml

  ! pointer to C array used by ifacel (0 to n-1 numbering)
  integer, dimension(:,:), pointer :: ifacel_0

  ! pointer to C array used by ifabor (0 to n-1 numbering)
  integer, dimension(:), pointer :: ifabor_0

  ! pointer to C array used by ipnfac (0 to n-1 numbering)
  integer, dimension(:), pointer :: ipnfac_0

  ! pointer to C array used by nodfac (0 to n-1 numbering)
  integer, dimension(:), pointer :: nodfac_0

  ! pointer to C array used by ipnfbr (0 to n-1 numbering)
  integer, dimension(:), pointer :: ipnfbr_0

  ! pointer to C array used by nodfbr (0 to n-1 numbering)
  integer, dimension(:), pointer :: nodfbr_0

  !> \anchor ifmfbr
  !> family number of the boundary faces. See \ref note_1
  integer, dimension(:), pointer :: ifmfbr

  !> \anchor ifmcel
  !> family number of the elements. See \ref note_1
  integer, dimension(:), pointer :: ifmcel

  !> \anchor isympa
  !> integer to mark out the "symmetry" (itypfb=isymet) boundary faces
  !> where the mass flow has to be canceled when the ALE module is switched
  !> off (these faces are impermeable).
  !> For instance, if the face ifac is symmetry face,
  !> isympa(ifac)=0, otherwise isympa(ifac)=1.
  integer, dimension(:), pointer :: isympa

  !> \anchor xyzcen
  !> coordinate of the cell centers
  double precision, dimension(:,:), pointer :: xyzcen

  !> \anchor surfac
  !> surface vector of the internal faces. Its norm is the surface of the face
  !> and it is oriented from \c ifacel(1,.) to \c ifacel(2,.)
  double precision, dimension(:,:), pointer :: surfac

  !> \anchor surfbo
  !> surface vector of the boundary faces. Its norm is the surface of the face
  !> and it is oriented outwards
  double precision, dimension(:,:), pointer :: surfbo

  !> \anchor suffac
  !> fluid surface vector of the internal faces. Its norm is the surface of the face
  !> and it is oriented from \c ifacel(1,.) to \c ifacel(2,.)
  double precision, dimension(:,:), pointer :: suffac

  !> \anchor suffbo
  !> surface vector of the boundary faces. Its norm is the surface of the face
  !> and it is oriented outwards
  double precision, dimension(:,:), pointer :: suffbo

  ! isolid_0
  ! integer to mark out the "solid" cells (where the fluid volume is 0).
  ! Only available when iporos > 0.
  ! When iporos = 0, this array has a unique value (isolid_0(1:1)=0).
  integer, dimension(:), pointer :: isolid_0

  !> \anchor cdgfac
  !> coordinates of the centers of the internal faces
  double precision, dimension(:,:), pointer :: cdgfac

  !> \anchor cdgfbo
  !> coordinates of the centers of the boundary faces
  double precision, dimension(:,:), pointer :: cdgfbo

  !> \anchor xyznod
  !> coordinates of the mesh vertices
  double precision, dimension(:,:), pointer :: xyznod

  !> \anchor volume
  !> volume of each cell
  double precision, dimension(:), pointer :: volume

  !> \anchor cell_f_vol
  !> fluid volume of each cell
  double precision, dimension(:), pointer :: cell_f_vol

  !> \anchor surfan
  !> norm of the surface vector of the internal faces
  double precision, dimension(:), pointer :: surfan

  !> \anchor surfbn
  !> norm of the surface of the boundary faces
  double precision, dimension(:), pointer :: surfbn

  !> \anchor suffan
  !> norm of the fluid surface vector of the internal faces
  double precision, dimension(:), pointer :: suffan

  !> \anchor suffbn
  !> norm of the fluid surface of the boundary faces
  double precision, dimension(:), pointer :: suffbn

  !> \anchor dist
  !> for every internal face, dot product of the vectors
  !> \f$ \vect{IJ}\f$ and \f$\vect{n}\f$.  I and J are respectively
  !> the centers of the first and the second neighboring cell.
  !> The vector \f$\vect{n}\f$ is the unit vector normal to the face
  !> and oriented from the first to the second cell
  double precision, dimension(:), pointer :: dist

  !> \anchor distb
  !> For every boundary face, dot product between the vectors
  !> \f$\vect{IF}\f$ and \f$\vect{n}\f$.
  !> I is the center of the neighboring cell. F is the face center.
  !> The vector \f$\vect{n}\f$ is the unit vector normal to the face and
  !> oriented to the exterior of the domain
  double precision, dimension(:), pointer :: distb

  !> \anchor pond
  !> weighting (Aij=pond Ai+(1-pond)Aj)
  !> for every internal face,
  !> \f$\displaystyle\frac{\vect{FJ}.\vect{n}}{\vect{IJ}.\vect{n}}\f$.
  !> With regard to the mesh quality, its ideal value is 0.5
  double precision, dimension(:), pointer :: pond

  !> \anchor dijpf
  !> vector I'J' for interior faces
  !> for every internal face, the three components of the vector
  !> \f$\vect{I'J'}\f$, where I' and J' are
  !> respectively the orthogonal projections of the neighboring cell
  !> centers I and J on a straight line orthogonal to the face and passing
  !> through its center
  double precision, dimension(:,:), pointer :: dijpf

  !> \anchor diipb
  !> vector II' for interior faces
  !> for every boundary face, the three components of the vector
  !> \f$\vect{II'}\f$. I' is the orthogonal projection of I,
  !> center of the neighboring cell, on the
  !> straight line perpendicular to the face and passign through its center
  double precision, dimension(:,:), pointer :: diipb

  !> \anchor dofij
  !> vector OF for interior faces
  !> for every internal face, the three components of the vector
  !> \f$\vect{OF}\f$. O is the intersection
  !> point between the face and the straight line joining the centers
  !> of the two neighboring cells. F is the face center
  double precision, dimension(:,:), pointer :: dofij

  !=============================================================================

contains

  !=============================================================================

  !> \anchor ifacel
  !> Index-numbers of the two (only) neighboring cells for each internal face

  elemental pure function ifacel(iside, ifac) result(icel)

    implicit none

    ! Parameters

    integer, intent(in) :: iside, ifac
    integer             :: icel

    ! Function body

    icel = ifacel_0(iside, ifac) + 1

  end function ifacel

  !=============================================================================

  !> \anchor ifabor
  !> index-number of the (unique) neighboring cell for each boundary face

  elemental pure function ifabor(ifac) result(icel)

    implicit none

    ! Parameters

    integer, intent(in) :: ifac
    integer             :: icel

    ! Function body

    icel = ifabor_0(ifac) + 1

  end function ifabor

  !=============================================================================

  !> \anchor ipnfac
  !> position of the first node of the each internal face in the array
  !> returned by \ref nodfac (see \ref note_3)

  elemental pure function ipnfac(ifac) result(ipn)

    implicit none

    ! Parameters

    integer, intent(in) :: ifac
    integer             :: ipn

    ! Function body

    ipn = ipnfac_0(ifac) + 1

  end function ipnfac

  !=============================================================================

  !> \anchor nodfac
  !> indexed-numbers of the nodes of each internal face
  !> (see \ref note_3)

  elemental pure function nodfac(ipn) result(inod)

    implicit none

    ! Parameters

    integer, intent(in) :: ipn
    integer             :: inod

    ! Function body

    inod = nodfac_0(ipn) + 1

  end function nodfac

  !=============================================================================

  !> \anchor ipnfbr
  !> position of the first node of the each boundary face in the array returned
  !> by \ref nodfbr (see \ref note_3)

  elemental pure function ipnfbr(ifac) result(ipn)

    implicit none

    ! Parameters

    integer, intent(in) :: ifac
    integer             :: ipn

    ! Function body

    ipn = ipnfbr_0(ifac) + 1

  end function ipnfbr

  !=============================================================================

  !> \anchor nodfbr
  !> indexed-numbers of the nodes of each boundary face
  !> (see \ref note_3)

  elemental pure function nodfbr(ipn) result(inod)

    implicit none

    ! Parameters

    integer, intent(in) :: ipn
    integer             :: inod

    ! Function body

    inod = nodfbr_0(ipn) + 1

  end function nodfbr

  !=============================================================================

  !> \anchor isolid
  !> integer to mark out the "solid" cells (where the fluid volume is 0).
  !> Only available when \ref optcal::iporos "iposros" > 0.
  !> When \ref iporos = 0, this array has a unique value (isolid_0(1:1)=0).

  elemental pure function isolid(iporos, iel) result(isol)

    implicit none

    ! Parameters

    integer, intent(in) :: iporos, iel
    integer             :: isol

    ! Function body

    isol = isolid_0(min(iporos, 1)*(iel - 1) + 1)

  end function isolid

  !> \}

end module mesh
