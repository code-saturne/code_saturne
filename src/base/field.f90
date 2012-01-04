!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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


! Module for field-related operations

module field

  !=============================================================================

  ! Temporary pointer values used for mapping

  double precision, dimension(:),   pointer :: field_tmp_scal => null()
  double precision, dimension(:,:), pointer :: field_tmp_vect => null()

contains

  !=============================================================================

  ! Define a field.

  subroutine flddef (name, iexten, itycat, ityloc, idim, ilved, iprev, ifield)

    implicit none

    ! Arguments

    character(len=*), intent(in) :: name    ! Field name
    integer, intent(in)          :: iexten  ! 1: intensive; 2: extensive
    integer, intent(in)          :: itycat  ! Field category (may be added)
                                            !   4: variable
                                            !   8: property
                                            !  16: postprocess
                                            !  32: accumulator
                                            !  64: user
    integer, intent(in)          :: ityloc  ! Location type
                                            !   0: none
                                            !   1: cells
                                            !   2: interior faces
                                            !   3: interior faces
                                            !   4: vertices
    integer, intent(in)          :: idim    ! Field dimension
    integer, intent(in)          :: ilved   ! 0: not intereaved; 1: interleaved
    integer, intent(in)          :: iprev   ! 0: no previous values, 1: previous

    integer, intent(out)         :: ifield  ! Id of defined field

    ! Local variables

    integer :: lname

    lname = len(name)

    call fldde1(name, lname, iexten, itycat, ityloc, idim, ilved, iprev, ifield)
    !==========

    return

  end subroutine flddef

  !=============================================================================

  ! Return an id associated with a given field name if present.

  ! If the field has not been defined previously, -1 is returned.

  subroutine fldfid (name, ifield)

    implicit none

    ! Arguments

    character(len=*), intent(in) :: name    ! Field name

    integer, intent(out)         :: ifield  ! Id of field

    ! Local variables

    integer :: lname

    lname = len(name)

    call fldfi1(name, lname, ifield)
    !==========

    return

  end subroutine fldfid

  !=============================================================================

  ! Return an id associated with a given key name if present.

  ! If the key has not been defined previously, -1 is returned.

  subroutine fldkid (name, ikey)

    implicit none

    ! Arguments

    character(len=*), intent(in) :: name    ! Key name

    integer, intent(out)         :: ikey    ! Id of key

    ! Local variables

    integer :: lname

    lname = len(name)

    call fldki1(name, lname, ikey)
    !==========

    return

  end subroutine fldkid

  !=============================================================================

  ! Assign a character string for a given key to a field.

  ! If the key id is not valid, or the value type or field category is not
  ! compatible, a fatal error is provoked.

  subroutine fldsks (ifield, ikey, str)

    implicit none

    ! Arguments

    integer, intent(in)          :: ifield   ! Field id
    integer, intent(in)          :: ikey     ! Key id
    character(len=*), intent(in) :: str      ! Associated string

    ! Local variables

    integer :: lstr

    lstr = len(str)

    call fldsk1(ifield, ikey, str, lstr)
    !==========

    return

  end subroutine fldsks

  !=============================================================================

  ! Return a character string for a given key associated with a field.

  ! If the key id is not valid, or the value type or field category is not
  ! compatible, a fatal error is provoked.

  subroutine fldgks (ifield, ikey, str)

    implicit none

    ! Arguments

    integer, intent(in)           :: ifield   ! Field id
    integer, intent(in)           :: ikey     ! Key id
    character(len=*), intent(out) :: str      ! Associated string

    ! Local variables

    integer :: lstr

    lstr = len(str)

    call fldgk1(ifield, ikey, str, lstr)
    !==========

    return

  end subroutine fldgks

  !=============================================================================

  ! Remove character X, x, U, u, or 1 from a Fortran character string if the
  ! compared strings are identical except for the last character, respectively
  ! Y, y, V,v, or 2 and Z, z, W, w, or 3.

  subroutine fldsnv (name1, name2, name3)

    implicit none

    ! Arguments

    character(len=*), intent(inout) :: name1  ! Name of base character string
    character(len=*), intent(in)    :: name2  ! Name of second character string
    character(len=*), intent(in)    :: name3  ! Name of third character string

    ! Local variables

    integer :: ii, jj
    integer :: lnmvar, lnmva2, lnmva3

    lnmvar = len(name1)
    lnmva2 = len(name2)
    lnmva3 = len(name3)

    if ((lnmvar .eq. lnmva2) .and. (lnmvar .eq. lnmva3)) then

      do ii = lnmvar, 1, -1
        if (    name1(ii:ii)  .ne. ' '         &
            .or. name2(ii:ii) .ne. ' '         &
            .or. name3(ii:ii) .ne. ' ') exit
      enddo

      if (ii .gt. 1) then

        jj = ii

        ! Handle the case where the next-to-last character changes, such
        ! as with VelocityX1, VelocityX2, ... in case of a calculation
        ! with multiple phases.

        if (      (ii .gt. 2)                                       &
            .and. (name1(ii:ii) .eq. name2(ii:ii))                  &
            .and. (name1(ii:ii) .eq. name3(ii:ii))) then
          ii = jj-1
        endif

        ! Remove the character related to the spatial axis

        if (      name1(ii:ii) .eq. 'X'                              &
            .and. name2(ii:ii) .eq. 'Y'                              &
            .and. name3(ii:ii) .eq. 'Z') then
          name1(ii:ii) = ' '
        else if (      name1(ii:ii) .eq. 'x'                         &
                 .and. name2(ii:ii) .eq. 'y'                         &
                 .and. name3(ii:ii) .eq. 'z') then
          name1(ii:ii) = ' '
        else if (      name1(ii:ii) .eq. 'U'                         &
                 .and. name2(ii:ii) .eq. 'V'                         &
                 .and. name3(ii:ii) .eq. 'W') then
          name1(ii:ii) = ' '
        else if (      name1(ii:ii) .eq. 'u'                         &
                 .and. name2(ii:ii) .eq. 'v'                         &
                 .and. name3(ii:ii) .eq. 'w') then
          name1(ii:ii) = ' '
        else if (      name1(ii:ii) .eq. '1'                         &
                 .and. name2(ii:ii) .eq. '2'                         &
                 .and. name3(ii:ii) .eq. '3') then
          name1(ii:ii) = ' '
        endif

        ! If the next-to last character was removed, shift the last one.

        if (ii .eq. jj+1) then
          name1(ii:ii) = name1(jj:jj)
          name1(jj:jj) = ' '
        endif

      endif

    endif

    return

  end subroutine fldsnv

  !=============================================================================

  ! Return a pointer to scalar field's values

  ! If the field id is not valid, a fatal error is provoked.

  ! Note that this function is not thread-safe.

  subroutine fldpts (ifield, iprev, val)

    implicit none

    ! Arguments

    integer, intent(in)                     :: ifield   ! Field id
    integer, intent(in)                     :: iprev    ! If 1, previous values
    double precision, dimension(:), pointer :: val      ! Associated pointer

    call fldps1(ifield, iprev)

    val => field_tmp_scal
    field_tmp_scal => null()

    return

  end subroutine fldpts

  !=============================================================================

  ! Return a pointer to vector field's values

  ! If the field id is not valid, a fatal error is provoked.

  ! Note that this function is not thread-safe.

  subroutine fldptv (ifield, iprev, val)

    implicit none

    ! Arguments

    integer, intent(in)                       :: ifield   ! Field id
    integer, intent(in)                       :: iprev    ! If 1, prev. values
    double precision, dimension(:,:), pointer :: val      ! Associated pointer

    call fldpv1(ifield, iprev)

    val => field_tmp_vect
    field_tmp_vect => null()

    return

  end subroutine fldptv

  !=============================================================================

end module field

!===============================================================================

! Subroutines defined outside of module so that their names are not mangled,
! as they must be callable from C code.

!===============================================================================

! Set global temporary scalar field pointer to null.

subroutine fldps2

  use field

  implicit none

  field_tmp_scal => null()

end subroutine fldps2

!===============================================================================

! Set global temporary scalar field pointer to a given array

subroutine fldps3(nval, val)

  use field

  implicit none

  ! Arguments

  integer, intent(in)                    :: nval
  double precision, dimension(*), target :: val

  ! Local variables

  field_tmp_scal => val(1:nval)

end subroutine fldps3

!===============================================================================

! Set global temporary vector field pointer to null.

subroutine fldpv2

  use field

  implicit none

  field_tmp_vect => null()

end subroutine fldpv2

!===============================================================================

! Set global temporary vector field pointer to a given array

subroutine fldpv3(nval1, nval2, val)

  use field

  implicit none

  ! Arguments

  integer, intent(in)                           :: nval1, nval2
  double precision, dimension(nval1, *), target :: val

  ! Local variables

  field_tmp_vect => val(1:nval1, 1:nval2)

end subroutine fldpv3
