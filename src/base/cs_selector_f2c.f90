!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!===============================================================================
! Function:
! ---------

!> \file cs_selector_f2c.f90

!> \brief Build the list of interior faces matching a criteria string.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     fstr          criteria string
!> \param[out]    facnb         number of selected faces
!> \param[out]    faces         selected faces
!_______________________________________________________________________________

subroutine getfac &
 ( fstr , facnb, faces)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

! Arguments

character*(*) fstr
integer       faces(*), facnb

! Local variables

character(len=len_trim(fstr)+1, kind=c_char) :: c_crit

!===============================================================================

interface

  subroutine cs_selector_get_i_face_num_list(criteria,                  &
                                             n_faces, i_face_num_list)  &
    bind(C, name='cs_selector_get_i_face_num_list')
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=1), dimension(*), intent(in) :: criteria
    integer(c_int), intent(out) :: n_faces
    integer(c_int), dimension(*), intent(out) :: i_face_num_list
  end subroutine cs_selector_get_i_face_num_list

end interface

c_crit = trim(fstr)//c_null_char

call cs_selector_get_i_face_num_list(c_crit, facnb, faces)

return

end subroutine getfac

!===============================================================================

!===============================================================================
! Function:
! ---------

!> \brief Build the list of boundary faces matching a criteria string.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     fstr          criteria string
!> \param[out]    facnb         number of selected faces
!> \param[out]    faces         selected faces
!_______________________________________________________________________________

subroutine getfbr &
 ( fstr , facnb, faces)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

! Arguments

character*(*)    fstr
integer      faces(*), facnb

! Local variables

character(len=len_trim(fstr)+1, kind=c_char) :: c_crit

!===============================================================================

interface

  subroutine cs_selector_get_b_face_num_list(criteria,                  &
                                             n_faces, b_face_num_list)  &
    bind(C, name='cs_selector_get_b_face_num_list')
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=1), dimension(*), intent(in) :: criteria
    integer(c_int), intent(out) :: n_faces
    integer(c_int), dimension(*), intent(out) :: b_face_num_list
  end subroutine cs_selector_get_b_face_num_list

end interface

c_crit = trim(fstr)//c_null_char

call cs_selector_get_b_face_num_list(c_crit, facnb, faces)

return

end subroutine getfbr

!===============================================================================

!===============================================================================
! Function:
! ---------

!> \brief Build the list of cells matching a criteria string.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     fstr          criteria string
!> \param[out]    cellnb        number of selected cells
!> \param[out]    cells         selected cells
!_______________________________________________________________________________


subroutine getcel &
 ( fstr , cellnb, cells)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

! Arguments

character*(*) fstr
integer       cells(*), cellnb

! Local variables

character(len=len_trim(fstr)+1, kind=c_char) :: c_crit

!===============================================================================

interface

  subroutine cs_selector_get_cell_num_list(criteria, n_faces, cell_num_list)  &
    bind(C, name='cs_selector_get_cell_num_list')
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=1), dimension(*), intent(in) :: criteria
    integer(c_int), intent(out) :: n_faces
    integer(c_int), dimension(*), intent(out) :: cell_num_list
  end subroutine cs_selector_get_cell_num_list

end interface

c_crit = trim(fstr)//c_null_char

call cs_selector_get_cell_num_list(c_crit, cellnb, cells)

return

end subroutine getcel

!===============================================================================

!===============================================================================
! Function:
! ---------

!> \brief Build the lists of faces at the boundary of cells matching a
!> criteria string.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     fstr          criteria string
!> \param[out]    ifaces        selected interior faces
!> \param[out]    bfaces        selected boundary faces
!> \param[out]    ifacnb        number fo selected interior faces
!> \param[out]    bfacnb        number fo selected boundary faces
!_______________________________________________________________________________

subroutine getceb &
 ( fstr , ifacnb, bfacnb, ifaces, bfaces )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

! Arguments

character*(*) fstr
integer       ifaces(*), bfaces(*), ifacnb, bfacnb

! Local variables

integer :: ii
character(len=len_trim(fstr)+1, kind=c_char) :: c_crit

!===============================================================================

interface

  subroutine cs_selector_get_cells_boundary(criteria, n_i_faces, n_b_faces,   &
                                            i_face_list, b_face_list)         &
    bind(C, name='cs_selector_get_cells_boundary')
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=1), dimension(*), intent(in) :: criteria
    integer(c_int), intent(out) :: n_i_faces, n_b_faces
    integer(c_int), dimension(*), intent(out) :: i_face_list, b_face_list
  end subroutine cs_selector_get_cells_boundary

end interface

c_crit = trim(fstr)//c_null_char

call cs_selector_get_cells_boundary(c_crit, ifacnb, bfacnb, ifaces, bfaces)

do ii = 1, ifacnb
  ifaces(ii) = ifaces(ii) + 1
enddo
do ii = 1, bfacnb
  bfaces(ii) = bfaces(ii) + 1
enddo

return

end subroutine getceb

!===============================================================================

!===============================================================================
! Function:
! ---------

!> \brief Build the list of mesh element families matching a criteria string.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     fstr          criteria string
!> \param[out]    families      selected families
!> \param[out]    famnb         number of selected families
!_______________________________________________________________________________

subroutine getfam &
 ( fstr , famnb, families)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

! Arguments

character*(*) fstr
integer       families(*), famnb

! Local variables

character(len=len_trim(fstr)+1, kind=c_char) :: c_crit

!===============================================================================

interface

  subroutine cs_selector_get_family_list(criteria, n_families, family_list)   &
    bind(C, name='cs_selector_get_family_list')
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=1), dimension(*), intent(in) :: criteria
    integer(c_int), intent(out) :: n_families
    integer(c_int), dimension(*), intent(out) :: family_list
  end subroutine cs_selector_get_family_list

end interface

c_crit = trim(fstr)//c_null_char

call cs_selector_get_family_list(c_crit, famnb, families)

return

end subroutine getfam

!===============================================================================
! Function:
! ---------

!> \brief Build the list of interior faces belonging to a given periodicity.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     perio_num     periodicity number
!> \param[out]    n_faces       number of faces
!> \param[out]    face_list     faces list
!_______________________________________________________________________________

subroutine getfpe &
 ( perio_num, n_faces, face_list )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

! Arguments

integer       face_list(*), perio_num, n_faces

! Local variables

integer       ii

!===============================================================================

interface

  subroutine cs_selector_get_perio_face_list(perio_num, n_faces, face_list)   &
    bind(C, name='cs_selector_get_perio_face_list')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: perio_num
    integer(c_int), intent(out) :: n_faces
    integer(c_int), dimension(*), intent(out) :: face_list
  end subroutine cs_selector_get_perio_face_list

end interface

call cs_selector_get_perio_face_list(perio_num, n_faces, face_list)

do ii = 1, n_faces
  face_list(ii) = face_list(ii) + 1
enddo

return

end subroutine getfpe

!===============================================================================


