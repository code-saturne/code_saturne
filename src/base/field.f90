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


!> \file field.f90
!> Module for field-related operations

module field

  !=============================================================================

  implicit none

  !=============================================================================

  integer :: FIELD_INTENSIVE, FIELD_EXTENSIVE, FIELD_STEADY
  integer :: FIELD_VARIABLE, FIELD_PROPERTY
  integer :: FIELD_POSTPROCESS, FIELD_ACCUMULATOR, FIELD_USER

  integer :: FIELD_OK, FIELD_INVALID_KEY_NAME, FIELD_INVALID_KEY_ID,   &
             FIELD_INVALID_CATEGORY, FIELD_INVALID_TYPE

  parameter (FIELD_INTENSIVE=1)
  parameter (FIELD_EXTENSIVE=2)
  parameter (FIELD_STEADY=4)
  parameter (FIELD_VARIABLE=8)
  parameter (FIELD_PROPERTY=16)
  parameter (FIELD_POSTPROCESS=32)
  parameter (FIELD_ACCUMULATOR=64)
  parameter (FIELD_USER=128)

  parameter (FIELD_OK=0)
  parameter (FIELD_INVALID_KEY_NAME=1)
  parameter (FIELD_INVALID_KEY_ID=2)
  parameter (FIELD_INVALID_CATEGORY=3)
  parameter (FIELD_INVALID_TYPE=4)

  !=============================================================================

  interface

    ! Interface to C function allocating field values

    !> \brief Allocate arrays for all defined fields based on their location.

    !> Location sized must thus be known.

    !> Fields that do not own their data should all have been mapped at this
    !> stage, and are checked.

    subroutine field_allocate_or_map_all()  &
      bind(C, name='cs_field_allocate_or_map_all')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine field_allocate_or_map_all

    !---------------------------------------------------------------------------

    ! Interface to C function assigning integer value to a key

    !> \brief Assign a floating point value for a given key to a field.

    !> If the key id is not valid, or the value type or field category is not
    !> compatible, a fatal error is provoked.

    !> \param[in]   f_id     field id
    !> \param[in]   k_id     id of associated key
    !> \param[in]   k_value  value associated with key

    subroutine field_set_key_int(f_id, k_id, k_value)  &
      bind(C, name='cs_f_field_set_key_int')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id, k_id, k_value
    end subroutine field_set_key_int

    !---------------------------------------------------------------------------

    ! Interface to C function assigning integer bit values to a key

    !> \brief Set integer bits matching a mask to 1 for a given key for a field.

    !> If the key id is not valid, or the value type or field category is not
    !> compatible, a fatal error is provoked.

    !> \param[in]   f_id  field id
    !> \param[in]   k_id  id of associated key
    !> \param[in]   mask  associated mask

    subroutine field_set_key_int_bits(f_id, k_id, k_value)  &
      bind(C, name='cs_f_field_set_key_int_bits')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id, k_id, k_value
    end subroutine field_set_key_int_bits

    !---------------------------------------------------------------------------

    ! Interface to C function assigning integer bit values to a key

    !> \brief Set integer bits matching a mask to 0 for a given key for a field.

    !> If the key id is not valid, or the value type or field category is not
    !> compatible, a fatal error is provoked.

    !> \param[in]   f_id  field id
    !> \param[in]   k_id  id of associated key
    !> \param[in]   mask  associated mask

    subroutine field_clear_key_int_bits(f_id, k_id, k_value)  &
      bind(C, name='cs_f_field_clear_key_int_bits')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id, k_id, k_value
    end subroutine field_clear_key_int_bits

    !---------------------------------------------------------------------------

    ! Interface to C function assigning floating-point value to a key

    !> \brief Assign a floating point value for a given key to a field.

    !> If the key id is not valid, or the value type or field category is not
    !> compatible, a fatal error is provoked.

    !> \param[in]   f_id     field id
    !> \param[in]   k_id     id of associated key
    !> \param[in]   k_value  value associated with key

    subroutine field_set_key_double(f_id, k_id, k_value)  &
      bind(C, name='cs_f_field_set_key_double')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id, k_id
      real(c_double), value :: k_value
    end subroutine field_set_key_double

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function creating a field descriptor

    function cs_field_create(name, type_flag, location_id, dim, &
                             has_previous) result(f) &
      bind(C, name='cs_field_create')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name
      integer(c_int), value                                    :: type_flag
      integer(c_int), value                                    :: location_id
      integer(c_int), value                                    :: dim
      logical(c_bool), value                                   :: has_previous
      type(c_ptr)                                              :: f
    end function cs_field_create

    !---------------------------------------------------------------------------

    ! Interface to C function returning or creating a field descriptor

    function cs_field_find_or_create(name, type_flag, location_id, &
                                     dim) result(f) &
      bind(C, name='cs_field_find_or_create')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name
      integer(c_int), value                                    :: type_flag
      integer(c_int), value                                    :: location_id
      integer(c_int), value                                    :: dim
      type(c_ptr)                                              :: f
    end function cs_field_find_or_create

    !---------------------------------------------------------------------------

    ! Interface to C function obtaining the number of fields

    function cs_f_field_n_fields() result(id) &
      bind(C, name='cs_f_field_n_fields')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int)                                           :: id
    end function cs_f_field_n_fields

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

    ! Interface to C function obtaining a field's location

    function cs_f_field_location(f) result(f_loc) &
      bind(C, name='cs_f_field_location')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: f
      integer(c_int)        :: f_loc
    end function cs_f_field_location

    !---------------------------------------------------------------------------

    ! Interface to C function obtaining a field's id by its name

    function cs_f_field_id_by_name_try(name) result(id) &
      bind(C, name='cs_f_field_id_by_name_try')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name
      integer(c_int)                                           :: id
    end function cs_f_field_id_by_name_try

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

    ! Interface to C function returning a given field name pointer and length.

    subroutine cs_f_field_get_name(f_id, f_name_max, f_name, f_name_len)  &
      bind(C, name='cs_f_field_get_name')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value       :: f_id
      integer(c_int), value       :: f_name_max
      type(c_ptr), intent(out)    :: f_name
      integer(c_int), intent(out) :: f_name_len
    end subroutine cs_f_field_get_name

    !---------------------------------------------------------------------------

    ! Interface to C function returning a given field's dimension info

    subroutine cs_f_field_get_dimension(f_id, f_dim)  &
      bind(C, name='cs_f_field_get_dimension')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      integer(c_int), dimension(1), intent(out) :: f_dim
    end subroutine cs_f_field_get_dimension

    !---------------------------------------------------------------------------

    ! Interface to C function returning a given field's ownership info

    subroutine cs_f_field_get_ownership(f_id, f_is_owner)  &
      bind(C, name='cs_f_field_get_ownership')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      logical(c_bool), intent(out) :: f_is_owner
    end subroutine cs_f_field_get_ownership

    !---------------------------------------------------------------------------

    ! Interface to C function returning a given field's type info

    subroutine cs_f_field_get_type(f_id, f_type)  &
      bind(C, name='cs_f_field_get_type')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      integer(c_int), intent(out) :: f_type
    end subroutine cs_f_field_get_type

    !---------------------------------------------------------------------------

    ! Interface to C function indicating if a field maintains a previous time

    function cs_f_field_have_previous(f_id) result(have_previous)  &
      bind(C, name='cs_f_field_have_previous')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      integer(c_int) :: have_previous
    end function cs_f_field_have_previous

    !---------------------------------------------------------------------------

    ! Interface to C function changing a fields handling of previous values

    subroutine cs_f_field_set_n_previous(f_id, n_previous)  &
      bind(C, name='cs_f_field_set_n_previous')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id, n_previous
    end subroutine cs_f_field_set_n_previous

    !---------------------------------------------------------------------------

    ! Interface to C function allocating field values

    subroutine cs_field_allocate_values(f)  &
      bind(C, name='cs_field_allocate_values')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: f
    end subroutine cs_field_allocate_values

    !---------------------------------------------------------------------------

    ! Interface to C function mapping field values

    subroutine cs_field_map_values(f, var, var_prev)  &
      bind(C, name='cs_field_map_values')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: f
      real(kind=c_double), dimension(*) :: var, var_prev
    end subroutine cs_field_map_values

    !---------------------------------------------------------------------------

    ! Interface to C function allocating boundary condition coefficients

    subroutine cs_field_allocate_bc_coeffs(f, have_flux_bc, have_mom_bc,  &
                                           have_conv_bc)                  &
      bind(C, name='cs_field_allocate_bc_coeffs')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value           :: f
      logical(c_bool), value       :: have_flux_bc
      logical(c_bool), value       :: have_mom_bc
      logical(c_bool), value       :: have_conv_bc
    end subroutine cs_field_allocate_bc_coeffs

    !---------------------------------------------------------------------------

    ! Interface to C function initializing boundary condition coefficients

    subroutine cs_field_init_bc_coeffs(f)                      &
      bind(C, name='cs_field_init_bc_coeffs')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value           :: f
    end subroutine cs_field_init_bc_coeffs

    !---------------------------------------------------------------------------

    ! Interface to C function copying current to previous values

    subroutine cs_field_current_to_previous(f)  &
      bind(C, name='cs_field_current_to_previous')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: f
    end subroutine cs_field_current_to_previous

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

    ! Interface to C function returning field's value pointer and dimensions.

    ! If the field id is not valid, a fatal error is provoked.

    subroutine cs_f_field_var_ptr_by_id_try(id, p_type, p_rank, f_dim, c_p)  &
      bind(C, name='cs_f_field_var_ptr_by_id_try')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value        :: id
      integer(c_int), value        :: p_type
      integer(c_int), value        :: p_rank
      integer(c_int), dimension(2) :: f_dim
      type(c_ptr), intent(out)     :: c_p
    end subroutine cs_f_field_var_ptr_by_id_try

    !---------------------------------------------------------------------------

    ! Interface to C function returning field's boundary condition
    ! coefficient values pointer and dimensions.

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

    !---------------------------------------------------------------------------

    ! Interface to C function obtaining a field key id by its name

    function cs_f_field_key_id(name) result(id) &
      bind(C, name='cs_field_key_id')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name
      integer(c_int)                                           :: id
    end function cs_f_field_key_id

    !---------------------------------------------------------------------------

    ! Interface to C function obtaining a field key id by its name

    function cs_f_field_key_id_try(name) result(id) &
      bind(C, name='cs_field_key_id_try')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name
      integer(c_int)                                           :: id
    end function cs_f_field_key_id_try

    !---------------------------------------------------------------------------

    ! Interface to C function querying if key value was defined

    function cs_field_is_key_set(f, k_id) result(is_set) &
      bind(C, name='cs_field_is_key_set')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: f
      integer(c_int), value :: k_id
      logical(c_bool)       :: is_set
    end function cs_field_is_key_set

    !---------------------------------------------------------------------------

    ! Interface to C function querying if key value was locked

    function cs_field_is_key_locked(f, k_id) result(is_locked) &
      bind(C, name='cs_field_is_key_locked')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: f
      integer(c_int), value :: k_id
      logical(c_bool)       :: is_locked
    end function cs_field_is_key_locked

    !---------------------------------------------------------------------------

    ! Interface to C function locking a key value

    subroutine cs_field_lock_key(f, k_id) &
      bind(C, name='cs_field_lock_key')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: f
      integer(c_int), value :: k_id
    end subroutine cs_field_lock_key

    !---------------------------------------------------------------------------

    ! Interface to C function assigning a character string for a given key
    ! to a field.

    subroutine cs_f_field_set_key_str(f_id, c_id, str) &
      bind(C, name='cs_f_field_set_key_str')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value                                    :: f_id, c_id
      character(kind=c_char, len=1), dimension(*), intent(in)  :: str
    end subroutine cs_f_field_set_key_str

    !---------------------------------------------------------------------------

    ! Interface to C function returning an integer for a given key associated
    ! with a field

    function cs_field_get_key_int(f, k_id) result(k_value) &
      bind(C, name='cs_field_get_key_int')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: f
      integer(c_int), value :: k_id
      integer(c_int)        :: k_value
    end function cs_field_get_key_int

    !---------------------------------------------------------------------------

    ! Interface to C function returning an floating-point valuer for a given
    ! key associated with a field

    function cs_field_get_key_double(f, k_id) result(k_value) &
      bind(C, name='cs_field_get_key_double')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: f
      integer(c_int), value :: k_id
      real(c_double)        :: k_value
    end function cs_field_get_key_double

    !---------------------------------------------------------------------------

    ! Interface to C function returning a string for a given key associated
    ! with a field.

    subroutine cs_f_field_get_key_str(f_id, k_id, str_max, str, str_len) &
      bind(C, name='cs_f_field_get_key_str')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value       :: f_id, k_id, str_max
      type(c_ptr), intent(out)    :: str
      integer(c_int), intent(out) :: str_len
    end subroutine cs_f_field_get_key_str

    !---------------------------------------------------------------------------

    ! Interface to C function copying a structure associated with a field.

    subroutine cs_f_field_set_key_struct(f_id, k_id, k_value) &
      bind(C, name='cs_f_field_set_key_struct')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: f_id, k_id
      type(c_ptr), value                :: k_value
    end subroutine cs_f_field_set_key_struct

    !---------------------------------------------------------------------------

    ! Interface to C function copying a structure associated with a field.

    subroutine cs_f_field_get_key_struct(f_id, k_id, k_value) &
      bind(C, name='cs_f_field_get_key_struct')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: f_id, k_id
      type(c_ptr), value                :: k_value
    end subroutine cs_f_field_get_key_struct

    !---------------------------------------------------------------------------

    ! Interface to C function returning a label associated with a field.

    subroutine cs_f_field_get_label(f_id, str_max, str, str_len) &
      bind(C, name='cs_f_field_get_label')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value       :: f_id, str_max
      type(c_ptr), intent(out)    :: str
      integer(c_int), intent(out) :: str_len
    end subroutine cs_f_field_get_label

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief  Define a field.

  !> \param[in]  name           field name
  !> \param[in]  type_flag      field categories (may be added)
  !> \param[in]  location_id    field location type:
  !>                              0: none
  !>                              1: cells
  !>                              2: interior faces
  !>                              3: interior faces
  !>                              4: vertices
  !> \param[in]  dim            field dimension
  !> \param[in]  has_previous   .true. if values at previous
  !>                            time step are maintained
  !> \param[out] id             id of defined field

  subroutine field_create(name, type_flag, location_id, dim, has_previous, &
                          id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name
    integer, intent(in)          :: type_flag
    integer, intent(in)          :: location_id
    integer, intent(in)          :: dim
    logical, intent(in)          :: has_previous
    integer, intent(out)         :: id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: c_type_flag
    integer(c_int) :: c_location_id
    integer(c_int) :: c_dim
    logical(c_bool) :: c_has_previous
    type(c_ptr)     :: f

    c_name = trim(name)//c_null_char
    c_type_flag = type_flag
    c_location_id = location_id
    c_dim = dim
    c_has_previous = has_previous

    f = cs_field_create(c_name, c_type_flag, c_location_id, c_dim, &
                        c_has_previous)
    id = cs_f_field_id_by_name(c_name)

    return

  end subroutine field_create

  !=============================================================================

  !> \brief  Return the id of a field matching a given name and attributes,
  !>         creating it if necessary.

  !> If a field with the same name but different attributes is present,
  !> this is considered an error.

  !> The default number of time values associated with a field created through
  !> this function is 1. To modify it, use \ref cs_field_set_n_time_vals.

  !> \param[in]  name           field name
  !> \param[in]  type_flag      field categories (may be added)
  !> \param[in]  location_id    field location type:
  !>                              0: none
  !>                              1: cells
  !>                              2: interior faces
  !>                              3: interior faces
  !>                              4: vertices
  !> \param[in]  dim            field dimension
  !> \param[out] id             id of defined field

  subroutine field_find_or_create(name, type_flag, location_id, dim, id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name
    integer, intent(in)          :: type_flag
    integer, intent(in)          :: location_id
    integer, intent(in)          :: dim
    integer, intent(out)         :: id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: c_type_flag
    integer(c_int) :: c_location_id
    integer(c_int) :: c_dim
    type(c_ptr)     :: f

    c_name = trim(name)//c_null_char
    c_type_flag = type_flag
    c_location_id = location_id
    c_dim = dim

    f = cs_field_find_or_create(c_name, c_type_flag, c_location_id, c_dim)
    id = cs_f_field_id_by_name(c_name)

    return

  end subroutine field_find_or_create

  !=============================================================================

  !> \brief  Return the number of defined fields.

  !> \param[out] nfld           number of field

  subroutine field_get_n_fields(nfld)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(out)         :: nfld

    nfld = cs_f_field_n_fields()

    return

  end subroutine field_get_n_fields

  !=============================================================================

  !> \brief  Return an id associated with a given field name.

  !> \param[in]  name           field name
  !> \param[out] id             id of field

  subroutine field_get_id(name, id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name
    integer, intent(out)         :: id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name

    c_name = trim(name)//c_null_char

    id = cs_f_field_id_by_name(c_name)

    return

  end subroutine field_get_id

  !=============================================================================

  !> \brief  Return  the location of a given field.

  !> \param[in]  f_id           field id
  !> \param[out] f_loc          location of the field

  subroutine field_get_location(f_id, f_loc)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)         :: f_id
    integer, intent(out)        :: f_loc

    ! Local variables

    integer(c_int) :: cf_id
    type(c_ptr)    :: f

    cf_id = f_id
    f = cs_field_by_id(cf_id)

    f_loc = cs_f_field_location(f)

    return

  end subroutine field_get_location

  !=============================================================================

  !> \brief  Return an id associated with a given field name if present.

  !> If the field has not been defined previously, -1 is returned.

  !> \param[in]  name           field name
  !> \param[out] id             id of field

  subroutine field_get_id_try(name, id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name
    integer, intent(out)         :: id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name

    c_name = trim(name)//c_null_char

    id = cs_f_field_id_by_name_try(c_name)

    return

  end subroutine field_get_id_try

  !=============================================================================

  !> \brief Return a given field's name.

  !> \param[in]   f_id  field id
  !> \param[out]  name  field's name

  subroutine field_get_name(f_id, name)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)           :: f_id
    character(len=*), intent(out) :: name

    ! Local variables

    integer :: i
    integer(c_int) :: c_f_id, name_max, c_name_len
    type(c_ptr) :: c_name_p
    character(kind=c_char, len=1), dimension(:), pointer :: c_name

    c_f_id = f_id
    name_max = len(name)

    call cs_f_field_get_name(f_id, name_max, c_name_p, c_name_len)
    call c_f_pointer(c_name_p, c_name, [c_name_len])

    do i = 1, c_name_len
      name(i:i) = c_name(i)
    enddo
    do i = c_name_len + 1, name_max
      name(i:i) = ' '
    enddo

    return

  end subroutine field_get_name

  !=============================================================================

  !> \brief Return a given field's dimension.

  !> \param[in]   f_id   field id
  !> \param[out]  f_dim  number of field components (dimension)

  subroutine field_get_dim(f_id, f_dim)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)  :: f_id
    integer, intent(out) :: f_dim

    ! Local variables

    integer(c_int) :: c_f_id
    integer(c_int), dimension(1) :: c_dim

    c_f_id = f_id

    call cs_f_field_get_dimension(c_f_id, c_dim)

    f_dim = c_dim(1)

    return

  end subroutine field_get_dim

  !=============================================================================

  !> \brief Return the field ownership flag.

  !> \param[in]   f_id         field id
  !> \param[out]  f_is_owner  true if field is owner, false otherwise

  subroutine field_get_ownership(f_id, f_is_owner)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)  :: f_id
    logical, intent(out) :: f_is_owner

    ! Local variables

    integer(c_int) :: c_f_id
    logical(c_bool) :: c_is_owner

    c_f_id = f_id

    call cs_f_field_get_ownership(c_f_id, c_is_owner)

    f_is_owner = c_is_owner

    return

  end subroutine field_get_ownership

  !=============================================================================

  !> \brief Return a given field's type.

  !> \param[in]   f_id         field id
  !> \param[out]  f_type       field type flag

  subroutine field_get_type(f_id, f_type)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)  :: f_id
    integer, intent(out) :: f_type

    ! Local variables

    integer(c_int) :: c_f_id
    integer(c_int) :: c_type

    c_f_id = f_id

    call cs_f_field_get_type(c_f_id, c_type)

    f_type = c_type

    return

  end subroutine field_get_type

  !=============================================================================

  !> \brief Indicate if a field maintains values at the previous time step

  !> \param[in]   f_id           field id
  !> \param[out]  have_previous  true if previous values are maintained

  subroutine field_have_previous(f_id, have_previous)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)  :: f_id
    logical, intent(out) :: have_previous

    ! Local variables

    integer(c_int) :: c_f_id, c_have_prev

    c_f_id = f_id

    c_have_prev = cs_f_field_have_previous(c_f_id)

    if (c_have_prev .eq. 0) then
      have_previous = .false.
    else
      have_previous = .true.
    endif

    return

  end subroutine field_have_previous

  !=============================================================================

  !> \brief Interface to C function locking a key value

  !> \param[in]   f_id        field id
  !> \param[in]   k_id        key id

  subroutine field_lock_key(f_id, k_id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in) :: f_id           ! Id of defined field
    integer, intent(in) :: k_id

    ! Local variables

    integer(c_int) :: cf_id, ck_id
    type(c_ptr)    :: f

    cf_id = f_id
    ck_id = k_id
    f = cs_field_by_id(cf_id)

    call cs_field_lock_key(f, ck_id)

  end subroutine field_lock_key

  !=============================================================================

  !> \brief Modify a field's handling of values at the previous time step

  !> \param[in]   f_id        field id
  !> \param[out]  n_previous  number of previous values (0 or 1)

  subroutine field_set_n_previous(f_id, n_previous)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)  :: f_id, n_previous

    ! Local variables

    integer(c_int) :: c_f_id, c_n_previous

    c_f_id = f_id
    c_n_previous = n_previous

    call cs_f_field_set_n_previous(c_f_id, c_n_previous)

    return

  end subroutine field_set_n_previous

  !=============================================================================

  !> \brief Allocate field's value arrays.

  !> \param[in]  id  field id

  subroutine field_allocate_values(id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in) :: id           ! Id of defined field

    ! Local variables

    integer(c_int) :: c_id
    type(c_ptr)    :: f

    c_id = id

    f = cs_field_by_id(c_id)
    call cs_field_allocate_values(f)

    return

  end subroutine field_allocate_values

  !=============================================================================

  !> \brief  Map existing value arrays to field descriptor.

  !> \param[in]  id       field id
  !> \param[in]  val      pointer to array of values
  !> \param[in]  val_pre  pointer to array of previous values (ignored if
  !>                      field was defined with have_previous = .false.)

  subroutine field_map_values(id, val, val_pre)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                        :: id
    double precision, intent(in), dimension(*) :: val, val_pre

    ! Local variables

    integer(c_int) :: c_id
    type(c_ptr)    :: f

    c_id = id

    f = cs_field_by_id(c_id)
    call cs_field_map_values(f, val, val_pre)

    return

  end subroutine field_map_values

  !=============================================================================

  !> \brief Allocate boundary condition coefficient arrays if applicable.

  !> \param[in]  id            field id
  !> \param[in]  have_flux_bc  if .true., flux BC coefficients
  !>                           (coefaf and coefbf) are added
  !> \param[in]  have_mom_bc   if .true., BC coefficients used in divergence
  !>                           term (coefad and coefbd) are added
  !> \param[in]  have_conv_bc   if .true., BC coefficients used in convection
  !>                           term (coefac and coefbc) are added

  subroutine field_allocate_bc_coeffs(id, have_flux_bc, have_mom_bc, have_conv_bc)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in) :: id
    logical, intent(in) :: have_flux_bc
    logical, intent(in) :: have_mom_bc
    logical, intent(in) :: have_conv_bc

    ! Local variables

    integer(c_int) :: c_id
    logical(c_bool) :: c_have_flux_bc
    logical(c_bool) :: c_have_mom_bc
    logical(c_bool) :: c_have_conv_bc
    type(c_ptr)     :: f

    c_id = id

    f = cs_field_by_id(c_id)
    c_have_flux_bc = have_flux_bc
    c_have_mom_bc = have_mom_bc
    c_have_conv_bc = have_conv_bc
    call cs_field_allocate_bc_coeffs(f, c_have_flux_bc, c_have_mom_bc, c_have_conv_bc)

    return

  end subroutine field_allocate_bc_coeffs

  !=============================================================================

  !> \brief Initialize boundary condition coefficient arrays if applicable.

  !> \param[in]  id            field id

  subroutine field_init_bc_coeffs(id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in) :: id

    ! Local variables

    integer(c_int) :: c_id
    type(c_ptr)    :: f

    c_id = id

    f = cs_field_by_id(c_id)
    call cs_field_init_bc_coeffs(f)

    return

  end subroutine field_init_bc_coeffs

  !=============================================================================

  !> \brief  Copy current values to previous values

  !> \param[in]  id  field id

  subroutine field_current_to_previous(id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in) :: id

    ! Local variables

    integer(c_int) :: c_id
    type(c_ptr)    :: f

    c_id = id

    f = cs_field_by_id(c_id)
    call cs_field_current_to_previous(f)

    return

  end subroutine field_current_to_previous

  !=============================================================================

  !> \brief  Query if a given key has been set for a field.

  !> If the key id is not valid, or the field category is not
  !> compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[in]   k_id     id of associated key
  !> \param[out]  is_set   is .true. if the field is set

  subroutine field_is_key_set(f_id, k_id, is_set)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)   :: f_id, k_id
    logical, intent(out)  :: is_set

    ! Local variables

    integer(c_int) :: c_f_id, c_k_id
    logical(c_bool) :: c_is_set
    type(c_ptr) :: f

    is_set = .false.

    c_f_id = f_id
    c_k_id = k_id
    f = cs_field_by_id(c_f_id)
    c_is_set = cs_field_is_key_set(f, k_id)
    if (c_is_set.eqv..true.) is_set = .true.

  end subroutine field_is_key_set

  !=============================================================================

  !> \brief  Return an id associated with a given key name if present.

  !> If the key has not been defined previously, -1 is returned.

  !> \param[in]   name  key name
  !> \param[out]  id    associated key id

  subroutine field_get_key_id(name, id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name
    integer, intent(out)         :: id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int)                               :: c_id

    c_name = trim(name)//c_null_char

    c_id = cs_f_field_key_id_try(c_name)
    id = c_id

    return

  end subroutine field_get_key_id

  !=============================================================================

  !> \brief Return an integer value for a given key associated with a field.

  !> If the key id is not valid, or the value type or field category is not
  !> compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[in]   k_id     id of associated key
  !> \param[out]  k_value  integer value associated with key id for this field

  subroutine field_get_key_int(f_id, k_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)   :: f_id, k_id
    integer, intent(out)  :: k_value

    ! Local variables

    integer(c_int) :: c_f_id, c_k_id, c_k_value
    type(c_ptr) :: f

    c_f_id = f_id
    c_k_id = k_id
    f = cs_field_by_id(c_f_id)
    c_k_value = cs_field_get_key_int(f, c_k_id)
    k_value = c_k_value

    return

  end subroutine field_get_key_int

  !=============================================================================

  !> \brief Return an integer value for a given key associated with a field.

  !> If the key id is not valid, or the value type or field category is not
  !> compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[in]   k_name   key name
  !> \param[out]  k_value  integer value associated with key id for this field

  subroutine field_get_key_int_by_name(f_id, k_name, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)   :: f_id
    character(len=*), intent(in) :: k_name
    integer, intent(out)  :: k_value

    ! Local variables

    integer(c_int) :: c_f_id, c_k_id, c_k_value
    character(len=len_trim(k_name)+1, kind=c_char) :: c_k_name
    type(c_ptr) :: f

    c_k_name = trim(k_name)//c_null_char

    c_k_id = cs_f_field_key_id_try(c_k_name)
    c_f_id = f_id
    f = cs_field_by_id(c_f_id)
    c_k_value = cs_field_get_key_int(f, c_k_id)
    k_value = c_k_value

    return

  end subroutine field_get_key_int_by_name

  !=============================================================================

  !> \brief Return a floating-point value for a given key associated with a
  !> field.

  !> If the key id is not valid, or the value type or field category is not
  !> compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[in]   k_id     id of associated key
  !> \param[out]  k_value  integer value associated with key id for this field

  subroutine field_get_key_double(f_id, k_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)            :: f_id, k_id
    double precision, intent(out)  :: k_value

    ! Local variables

    integer(c_int) :: c_f_id, c_k_id
    real(c_double) :: c_k_value
    type(c_ptr) :: f

    c_f_id = f_id
    c_k_id = k_id
    f = cs_field_by_id(c_f_id)
    c_k_value = cs_field_get_key_double(f, k_id)
    k_value = c_k_value

    return

  end subroutine field_get_key_double

  !=============================================================================

  !> \brief Assign a character string for a given key to a field.

  !> If the key id is not valid, or the value type or field category is not
  !> compatible, a fatal error is provoked.

  !> \param[in]   f_id  field id
  !> \param[in]   k_id  id of associated key
  !> \param[in]   str   string associated with key

  subroutine field_set_key_str(f_id, k_id, str)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)          :: f_id, k_id
    character(len=*), intent(in) :: str

    ! Local variables

    integer(c_int) :: c_f_id, c_k_id
    character(len=len_trim(str)+1, kind=c_char) :: c_str

    c_f_id = f_id
    c_k_id = k_id
    c_str = trim(str)//c_null_char

    call cs_f_field_set_key_str(c_f_id, c_k_id, c_str)

    return

  end subroutine field_set_key_str

  !=============================================================================

  !> \brief Return a character string for a given key associated with a field.

  !> If the key id is not valid, or the value type or field category is not
  !> compatible, a fatal error is provoked.

  !> \param[in]   f_id  field id
  !> \param[in]   k_id  id of associated key
  !> \param[out]  str   string associated with key

  subroutine field_get_key_str(f_id, k_id, str)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)           :: f_id, k_id
    character(len=*), intent(out) :: str

    ! Local variables

    integer :: i
    integer(c_int) :: c_f_id, c_k_id, c_str_max, c_str_len
    type(c_ptr) :: c_str_p
    character(kind=c_char, len=1), dimension(:), pointer :: c_str

    c_f_id = f_id
    c_k_id = k_id
    c_str_max = len(str)

    call cs_f_field_get_key_str(c_f_id, c_k_id, c_str_max, c_str_p, c_str_len)
    call c_f_pointer(c_str_p, c_str, [c_str_len])

    do i = 1, c_str_len
      str(i:i) = c_str(i)
    enddo
    do i = c_str_len + 1, c_str_max
      str(i:i) = ' '
    enddo

    return

  end subroutine field_get_key_str

  !=============================================================================

  ! Remove character X, x, U, u, or 1 from a Fortran character string if the
  ! compared strings are identical except for the last character, respectively
  ! Y, y, V,v, or 2 and Z, z, W, w, or 3.

  subroutine fldsnv(name1, name2, name3)

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

  !> \brief Return a label associated with a field.

  !> If the "label" key has been set for this field, its associated string
  !> is returned. Otherwise, the field's name is returned.

  !> \param[in]   f_id  field id
  !> \param[out]  str   string associated with key

  subroutine field_get_label(f_id, str)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)           :: f_id
    character(len=*), intent(out) :: str

    ! Local variables

    integer :: i
    integer(c_int) :: c_f_id, c_str_max, c_str_len
    type(c_ptr) :: c_str_p
    character(kind=c_char, len=1), dimension(:), pointer :: c_str

    c_f_id = f_id
    c_str_max = len(str)

    call cs_f_field_get_label(c_f_id, c_str_max, c_str_p, c_str_len)
    call c_f_pointer(c_str_p, c_str, [c_str_len])

    do i = 1, c_str_len
      str(i:i) = c_str(i)
    enddo
    do i = c_str_len + 1, c_str_max
      str(i:i) = ' '
    enddo

    return

  end subroutine field_get_label

  !=============================================================================

  !> \brief Return pointer to the values array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to scalar field values

  subroutine field_get_val_s(field_id, p)

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

  end subroutine field_get_val_s

  !=============================================================================

  !> \brief Return pointer to the values array of a given scalar field

  !> \param[in]     name      name of given field (which must be scalar)
  !> \param[out]    p         pointer to scalar field values

  subroutine field_get_val_s_by_name(name, p)

    use, intrinsic :: iso_c_binding
    implicit none

    character(len=*), intent(in)                         :: name
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(2) :: f_dim
    type(c_ptr) :: c_p

    c_name = trim(name)//c_null_char

    f_id = cs_f_field_id_by_name(c_name)
    p_type = 1
    p_rank = 1

    call cs_f_field_var_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_get_val_s_by_name

  !=============================================================================

  !> \brief Return pointer to the array's previous values of a given scalar
  !> field

  !> \param[in]     name      name of given field (which must be scalar)
  !> \param[out]    p         pointer to scalar field values at the previous
  !>                          iteration

  subroutine field_get_val_prev_s_by_name(name, p)

    use, intrinsic :: iso_c_binding
    implicit none

    character(len=*), intent(in)                         :: name
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(2) :: f_dim
    type(c_ptr) :: c_p

    c_name = trim(name)//c_null_char

    f_id = cs_f_field_id_by_name(c_name)

    p_type = 2
    p_rank = 1

    call cs_f_field_var_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_get_val_prev_s_by_name

  !=============================================================================


  !> \brief Return pointer to the values array of a given vector field

  !> \param[in]     field_id  id of given field (which must be vectorial)
  !> \param[out]    p         pointer to vector field values

  subroutine field_get_val_v(field_id, p)

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

  end subroutine field_get_val_v

!=============================================================================

  !> \brief Return pointer to the values array of a given vector field

  !> \param[in]     name      name of given field (which must be vectorial)
  !> \param[out]    p         pointer to scalar field values

  subroutine field_get_val_v_by_name(name, p)

    use, intrinsic :: iso_c_binding
    implicit none

    character(len=*), intent(in)                           :: name
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(2) :: f_dim
    type(c_ptr) :: c_p

    c_name = trim(name)//c_null_char

    f_id = cs_f_field_id_by_name(c_name)
    p_type = 1
    p_rank = 2

    call cs_f_field_var_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_get_val_v_by_name

  !=============================================================================

  !> \brief Return pointer to the array's previous values of a given vector
  !>        field

  !> \param[in]     name      name of given field (which must be vectorial)
  !> \param[out]    p         pointer to vector field values at the previous
  !>                          iteration

  subroutine field_get_val_prev_v_by_name(name, p)

    use, intrinsic :: iso_c_binding
    implicit none

    character(len=*), intent(in)                           :: name
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(2) :: f_dim
    type(c_ptr) :: c_p

    c_name = trim(name)//c_null_char

    f_id = cs_f_field_id_by_name(c_name)
    p_type = 2
    p_rank = 2

    call cs_f_field_var_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_get_val_prev_v_by_name

  !=============================================================================

  !> \brief Return pointer to the previous values array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to previous scalar field values

  subroutine field_get_val_prev_s(field_id, p)

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

  end subroutine field_get_val_prev_s

  !=============================================================================

  !> \brief Return pointer to the previous values array of a given scalar field
  !> if it exists, to the current value otherwise

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to previous scalar field values

  subroutine field_get_val_prev_s_try(field_id, p)

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

    call cs_f_field_var_ptr_by_id_try(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_get_val_prev_s_try

  !=============================================================================

  !> \brief Return pointer to the previous values array of a given vector field

  !> \param[in]     field_id  id of given field (which must be vectorial)
  !> \param[out]    p         pointer to previous vector field values

  subroutine field_get_val_prev_v(field_id, p)

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

  end subroutine field_get_val_prev_v

  !=============================================================================

  !> \brief Return pointer to the coefa array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to scalar field BC coefa values

  subroutine field_get_coefa_s(field_id, p)

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

  end subroutine field_get_coefa_s

  !=============================================================================

  !> \brief Return pointer to the coefa array of a given vector field

  !> \param[in]     field_id  id of given field (which must be a vector)
  !> \param[out]    p         pointer to vector field BC coefa values

  subroutine field_get_coefa_v(field_id, p)

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

  end subroutine field_get_coefa_v

  !=============================================================================

  !> \brief Return pointer to the coefad array of a given scalar field
  !>        (used in the divergence operator such as div(Rij))

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to vector field BC coefa values

  subroutine field_get_coefad_s(field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                  :: field_id
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 5
    p_rank = 1

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_get_coefad_s

  !=============================================================================

  !> \brief Return pointer to the coefa array of a given vector field
  !>        (used in the divergence operator such as div(u'T'))

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to vector field BC coefa values

  subroutine field_get_coefad_v(field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                    :: field_id
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 5
    p_rank = 2

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_get_coefad_v

  !=============================================================================

  !> \brief Return pointer to the coefac array of a given vector field

  !> \param[in]     field_id  id of given field (which must be a vector)
  !> \param[out]    p         pointer to vector field BC coefa values

  subroutine field_get_coefac_v(field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                    :: field_id
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 7
    p_rank = 2

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_get_coefac_v

  !=============================================================================

  !> \brief Return pointer to the coefb array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to scalar field BC coefa values

  subroutine field_get_coefb_s(field_id, p)

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

  end subroutine field_get_coefb_s

  !=============================================================================

  !> \brief Return pointer to the coefbc array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to scalar field BC coefbc values

  subroutine field_get_coefbc_s(field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                  :: field_id
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 8
    p_rank = 1

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_get_coefbc_s

  !=============================================================================

  !> \brief Return pointer to the coefb array of a given uncoupled vector field

  !> \param[in]     field_id  id of given field (which must be a vector)
  !> \param[out]    p         pointer to vector field BC coefa values

  subroutine field_get_coefb_uv(field_id, p)

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

  end subroutine field_get_coefb_uv

  !=============================================================================

  !> \brief Return pointer to the coefbc array of a given uncoupled vector field

  !> \param[in]     field_id  id of given field (which must be a vector)
  !> \param[out]    p         pointer to vector field BC coefbc values

  subroutine field_get_coefbc_uv(field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                    :: field_id
    double precision, dimension(:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 8
    p_rank = 2

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2)])

  end subroutine field_get_coefbc_uv

  !=============================================================================

  !> \brief Return pointer to the coefb array of a given coupled vector field

  !> \param[in]     field_id  id of given field (which must be a vector)
  !> \param[out]    p         pointer to vector field BC coefa values

  subroutine field_get_coefb_v(field_id, p)

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

  end subroutine field_get_coefb_v

  !=============================================================================

  !> \brief Return pointer to the coefbc array of a given coupled vector field

  !> \param[in]     field_id  id of given field (which must be a vector)
  !> \param[out]    p         pointer to vector field BC coefa values

  subroutine field_get_coefbc_v(field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                      :: field_id
    double precision, dimension(:,:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 8
    p_rank = 3

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2), f_dim(3)])

  end subroutine field_get_coefbc_v

  !=============================================================================

  !> \brief Return pointer to the coefaf array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to scalar field BC coefa values

  subroutine field_get_coefaf_s(field_id, p)

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

  end subroutine field_get_coefaf_s

  !=============================================================================

  !> \brief Return pointer to the coefac array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to scalar field BC coefa values

  subroutine field_get_coefac_s(field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                  :: field_id
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 7
    p_rank = 1

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_get_coefac_s

  !=============================================================================

  !> \brief Return pointer to the coefaf array of a given vector field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to vector field BC coefa values

  subroutine field_get_coefaf_v(field_id, p)

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

  end subroutine field_get_coefaf_v

  !=============================================================================

  !> \brief Return pointer to the coefbf array of a given scalar field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to scalar field BC coefb values

  subroutine field_get_coefbf_s(field_id, p)

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

  end subroutine field_get_coefbf_s

  !=============================================================================

  !> \brief Return pointer to the coefbf array of a given uncoupled vector field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to vector field BC coefb values

  subroutine field_get_coefbf_uv(field_id, p)

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

  end subroutine field_get_coefbf_uv

  !=============================================================================

  !> \brief Return pointer to the coefbf array of a given coupled vector field

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to vector field BC coefb values

  subroutine field_get_coefbf_v(field_id, p)

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

  end subroutine field_get_coefbf_v

  !=============================================================================

  !> \brief Return pointer to the coefbd array of a given scalar field
  !>        (used in the divergence operator such as div(Rij))

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to vector field BC coefb values

  subroutine field_get_coefbd_s(field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                  :: field_id
    double precision, dimension(:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 6
    p_rank = 1

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1)])

  end subroutine field_get_coefbd_s

  !=============================================================================

  !> \brief Return pointer to the coefbd array of a given coupled vector field
  !>        (used in the divergence operator such as div(u'T'))

  !> \param[in]     field_id  id of given field (which must be scalar)
  !> \param[out]    p         pointer to vector field BC coefa values

  subroutine field_get_coefbd_v(field_id, p)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                      :: field_id
    double precision, dimension(:,:,:), pointer, intent(out) :: p

    ! Local variables

    integer(c_int) :: f_id, p_type, p_rank
    integer(c_int), dimension(3) :: f_dim
    type(c_ptr) :: c_p

    f_id = field_id
    p_type = 6
    p_rank = 3

    call cs_f_field_bc_coeffs_ptr_by_id(f_id, p_type, p_rank, f_dim, c_p)
    call c_f_pointer(c_p, p, [f_dim(1), f_dim(2), f_dim(3)])

  end subroutine field_get_coefbd_v

  !=============================================================================

end module field
