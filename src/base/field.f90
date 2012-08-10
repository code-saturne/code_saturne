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

  integer FIELD_INTENSIVE, FIELD_EXTENSIVE
  integer FIELD_VARIABLE, FIELD_PROPERTY
  integer FIELD_POSTPROCESS, FIELD_ACCUMULATOR, FIELD_USER

  parameter (FIELD_INTENSIVE=1)
  parameter (FIELD_EXTENSIVE=2)
  parameter (FIELD_VARIABLE=4)
  parameter (FIELD_PROPERTY=8)
  parameter (FIELD_POSTPROCESS=16)
  parameter (FIELD_ACCUMULATOR=32)
  parameter (FIELD_USER=64)

  !=============================================================================

  interface

    ! Interface to C function creating a field descriptor

    function cs_field_create(name, type_flag, location_id, dim, interleaved, &
                             has_previous) result(f) &
      bind(C, name='cs_field_create')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name
      integer(c_int), value                                    :: type_flag
      integer(c_int), value                                    :: location_id
      integer(c_int), value                                    :: dim
      logical(c_bool), value                                   :: interleaved
      logical(c_bool), value                                   :: has_previous
      type(c_ptr)                                              :: f
    end function cs_field_create

    !---------------------------------------------------------------------------

    ! Interface to C function obtaining a field's id by its name

    function cs_f_field_id_by_name(name) result(id) &
      bind(C, name='cs_f_field_id_by_name')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name
      integer(c_int)                                           :: id
    end function cs_f_field_id_by_name

    !---------------------------------------------------------------------------

    ! Interface to C function obtaining field's pointer by its id

    function cs_field_by_id(id) result(f) &
      bind(C, name='cs_field_by_id')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: id
      type(c_ptr)           :: f
    end function cs_field_by_id

    !---------------------------------------------------------------------------

    ! Interface to C function allocating boundary condition coefficients

    subroutine cs_field_allocate_bc_coeffs(f, have_flux_bc)  &
      bind(C, name='cs_field_allocate_bc_coeffs')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value           :: f
      logical(c_bool), value       :: have_flux_bc
    end subroutine cs_field_allocate_bc_coeffs

    !---------------------------------------------------------------------------

    ! Interface to C function returning field's value pointer and dimensions.

    ! If the field id is not valid, a fatal error is provoked.

    subroutine cs_f_field_var_ptr_by_id(id, p_type, p_rank, f_dim, c_p)  &
      bind(C, name='cs_f_field_var_ptr_by_id')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value        :: id
      integer(c_int), value        :: p_type
      integer(c_int), value        :: p_rank
      integer(c_int), dimension(2) :: f_dim
      type(c_ptr), intent(out)     :: c_p
    end subroutine cs_f_field_var_ptr_by_id

    !---------------------------------------------------------------------------

    ! Interface to C function returning field's boundary condition
    ! coefficient valuesvalue pointer and dimensions.

    ! If the field id is not valid, a fatal error is provoked.

    subroutine cs_f_field_bc_coeffs_ptr_by_id(id, p_type, p_rank, f_dim, c_p)  &
      bind(C, name='cs_f_field_bc_coeffs_ptr_by_id')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value        :: id
      integer(c_int), value        :: p_type
      integer(c_int), value        :: p_rank
      integer(c_int), dimension(3) :: f_dim
      type(c_ptr), intent(out)     :: c_p
    end subroutine cs_f_field_bc_coeffs_ptr_by_id

  end interface

  !=============================================================================

contains

  !=============================================================================

  ! Define a field.

  subroutine field_create(name, type_flag, location_id, dim,   &
                          interleaved, has_previous,           &
                          id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name         ! Field name
    integer, intent(in)          :: type_flag    ! Field category (may be added)
    integer, intent(in)          :: location_id  ! Location type
                                                 !   0: none
                                                 !   1: cells
                                                 !   2: interior faces
                                                 !   3: interior faces
                                                 !   4: vertices
    integer, intent(in)          :: dim          ! Field dimension
    logical, intent(in)          :: interleaved  ! true if values interleaved
                                                 ! ignored if < 2 components
    logical, intent(in)          :: has_previous ! true if values at previous
                                                 ! time step are maintained

    integer, intent(out)         :: id           ! Id of defined field

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: c_type_flag
    integer(c_int) :: c_location_id
    integer(c_int) :: c_dim
    logical(c_bool) :: c_interleaved
    logical(c_bool) :: c_has_previous
    type(c_ptr)     :: f

    c_name = trim(name)//c_null_char
    c_type_flag = type_flag
    c_location_id = location_id
    c_dim = dim
    c_interleaved = interleaved
    c_has_previous = has_previous

    f = cs_field_create(c_name, c_type_flag, c_location_id, c_dim, &
                        c_interleaved, c_has_previous)
    id = cs_f_field_id_by_name(c_name)

    return

  end subroutine field_create

  !=============================================================================

  ! Return an id associated with a given field name if present.

  ! If the field has not been defined previously, -1 is returned.

  subroutine field_id_get(name, id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name         ! Field name
    integer, intent(out)         :: id           ! Id of defined field

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name

    c_name = trim(name)//c_null_char

    id = cs_f_field_id_by_name(c_name)

    return

  end subroutine field_id_get

  !=============================================================================

  ! Allocate boundary condition coefficient arrays if applicable.

  subroutine field_allocate_bc_coeffs(id, have_flux_bc)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in) :: id           ! Id of defined field
    logical, intent(in) :: have_flux_bc

    ! Local variables

    integer(c_int) :: c_id
    logical(c_bool) :: c_have_flux_bc
    type(c_ptr)     :: f

    f = cs_field_by_id(id)
    c_have_flux_bc = have_flux_bc
    call cs_field_allocate_bc_coeffs(f, c_have_flux_bc)

    return

  end subroutine field_allocate_bc_coeffs

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

  !> \brief Return pointer to the values array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to scalar field values

  subroutine field_val_get_s (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                  :: field_id
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(2) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 1
    p_rank = 1

    call cs_f_field_var_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_val_get_s

  !=============================================================================

  !> \brief Return pointer to the values array of a given vector field

  !> \param[in]     field_id  id of given field (which must be vectorial)
  !> \param[out]    pointer to vector field values

  subroutine field_val_get_v (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                    :: field_id
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(2) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 1
    p_rank = 2

    call cs_f_field_var_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_val_get_v

  !=============================================================================

  !> \brief Return pointer to the previous values array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to previous scalr field values

  subroutine field_val_prev_get_s (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                  :: field_id
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 2
    p_rank = 1

    call cs_f_field_var_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_val_prev_get_s

  !=============================================================================

  !> \brief Return pointer to the previous values array of a given vector field

  !> \param[in]     field_id  id of given field (which must be vectorial)
  !> \param[out]    pointer to previous vector field values

  subroutine field_val_prev_get_v (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                    :: field_id
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 2
    p_rank = 2

    call cs_f_field_var_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_val_prev_get_v

  !=============================================================================

  !> \brief Return pointer to the coefa array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to scalar field BC coefa values

  subroutine field_coefa_get_s (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                  :: field_id
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 1
    p_rank = 1

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_coefa_get_s

  !=============================================================================

  !> \brief Return pointer to the coefa array of a given vector field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to vector field BC coefa values

  subroutine field_coefa_get_v (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                    :: field_id
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 1
    p_rank = 2

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_coefa_get_v

  !=============================================================================

  !> \brief Return pointer to the coefb array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to scalar field BC coefa values

  subroutine field_coefb_get_s (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                  :: field_id
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 2
    p_rank = 1

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_coefb_get_s

  !=============================================================================

  !> \brief Return pointer to the coefb array of a given uncoupled vector field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to vector field BC coefa values

  subroutine field_coefb_get_uv (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                    :: field_id
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 2
    p_rank = 2

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_coefb_get_uv

  !=============================================================================

  !> \brief Return pointer to the coefb array of a given coupled vector field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to vector field BC coefa values

  subroutine field_coefb_get_v (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                      :: field_id
    double precision, dimension(:,:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 2
    p_rank = 3

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2), f_dim(3)])

  end subroutine field_coefb_get_v

  !=============================================================================

  !> \brief Return pointer to the coefa array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to scalar field BC coefa values

  subroutine field_coefaf_get_s (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                  :: field_id
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 3
    p_rank = 1

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_coefaf_get_s

  !=============================================================================

  !> \brief Return pointer to the coefa array of a given vector field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to vector field BC coefa values

  subroutine field_coefaf_get_v (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                    :: field_id
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 3
    p_rank = 2

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_coefaf_get_v

  !=============================================================================

  !> \brief Return pointer to the coefb array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to scalar field BC coefa values

  subroutine field_coefbf_get_s (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                  :: field_id
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 4
    p_rank = 1

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_coefbf_get_s

  !=============================================================================

  !> \brief Return pointer to the coefb array of a given uncoupled vector field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to vector field BC coefa values

  subroutine field_coefbf_get_uv (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                    :: field_id
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 4
    p_rank = 2

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_coefbf_get_uv

  !=============================================================================

  !> \brief Return pointer to the coefb array of a given coupled vector field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    pointer to vector field BC coefa values

  subroutine field_coefbf_get_v (field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                      :: field_id
    double precision, dimension(:,:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 4
    p_rank = 3

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2), f_dim(3)])

  end subroutine field_coefbf_get_v

  !=============================================================================

end module field
