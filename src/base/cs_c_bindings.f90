!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!> \file cs_c_bindings.f90
!> Definition of C function and subroutine bindings.

module cs_c_bindings

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use field

  implicit none

  !=============================================================================

  integer :: MESH_LOCATION_NONE, MESH_LOCATION_CELLS
  integer :: MESH_LOCATION_INTERIOR_FACES, MESH_LOCATION_BOUNDARY_FACES
  integer :: MESH_LOCATION_VERTICES, MESH_LOCATION_PARTICLES
  integer :: MESH_LOCATION_OTHER

  integer :: POST_ON_LOCATION, POST_BOUNDARY_NR, POST_MONITOR

  integer :: RESTART_VAL_TYPE_INT_T, RESTART_VAL_TYPE_REAL_T

  integer :: RESTART_DISABLED, RESTART_MAIN, RESTART_AUXILIARY
  integer :: RESTART_RAD_TRANSFER, RESTART_LAGR, RESTART_LAGR_STAT
  integer :: RESTART_1D_WALL_THERMAL, RESTART_LES_INFLOW

  integer :: VOLUME_ZONE_INITIALIZATION, VOLUME_ZONE_POROSITY
  integer :: VOLUME_ZONE_HEAD_LOSS
  integer :: VOLUME_ZONE_SOURCE_TERM, VOLUME_ZONE_MASS_SOURCE_TERM

  parameter (MESH_LOCATION_NONE=0)
  parameter (MESH_LOCATION_CELLS=1)
  parameter (MESH_LOCATION_INTERIOR_FACES=2)
  parameter (MESH_LOCATION_BOUNDARY_FACES=3)
  parameter (MESH_LOCATION_VERTICES=4)
  parameter (MESH_LOCATION_PARTICLES=5)
  parameter (MESH_LOCATION_OTHER=6)

  parameter (POST_ON_LOCATION=1)
  parameter (POST_BOUNDARY_NR=2)
  parameter (POST_MONITOR=4)

  parameter (RESTART_VAL_TYPE_INT_T=1)
  parameter (RESTART_VAL_TYPE_REAL_T=3)

  parameter (RESTART_DISABLED=-1)
  parameter (RESTART_MAIN=0)
  parameter (RESTART_AUXILIARY=1)
  parameter (RESTART_RAD_TRANSFER=2)
  parameter (RESTART_LAGR=3)
  parameter (RESTART_LAGR_STAT=4)
  parameter (RESTART_1D_WALL_THERMAL=5)
  parameter (RESTART_LES_INFLOW=6)

  parameter (VOLUME_ZONE_INITIALIZATION=1)
  parameter (VOLUME_ZONE_POROSITY=2)
  parameter (VOLUME_ZONE_HEAD_LOSS=4)
  parameter (VOLUME_ZONE_SOURCE_TERM=8)
  parameter (VOLUME_ZONE_MASS_SOURCE_TERM=16)

  !-----------------------------------------------------------------------------

  type, bind(c)  :: var_cal_opt
    integer(c_int) :: iwarni
    integer(c_int) :: iconv
    integer(c_int) :: istat
    integer(c_int) :: idircl
    integer(c_int) :: ndircl
    integer(c_int) :: idiff
    integer(c_int) :: idifft
    integer(c_int) :: idften
    integer(c_int) :: iswdyn
    integer(c_int) :: ischcv
    integer(c_int) :: ibdtso
    integer(c_int) :: isstpc
    integer(c_int) :: nswrgr
    integer(c_int) :: nswrsm
    integer(c_int) :: imvisf
    integer(c_int) :: imrgra
    integer(c_int) :: imligr
    integer(c_int) :: ircflu
    integer(c_int) :: iwgrec
    integer(c_int) :: icoupl
    real(c_double) :: thetav
    real(c_double) :: blencv
    real(c_double) :: blend_st
    real(c_double) :: epsilo
    real(c_double) :: epsrsm
    real(c_double) :: epsrgr
    real(c_double) :: climgr
    real(c_double) :: relaxv
  end type var_cal_opt

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    ! Interface to C function copying a var_cal_opt structure associated
    ! with a field.

    subroutine cs_f_field_set_key_struct_var_cal_opt(f_id, k_value) &
      bind(C, name='cs_f_field_set_key_struct_var_cal_opt')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: f_id
      type(c_ptr), value                :: k_value
    end subroutine cs_f_field_set_key_struct_var_cal_opt

    !---------------------------------------------------------------------------

    ! Interface to C function setting a var_cal_opt structure associated
    ! with a field.

    subroutine cs_f_field_get_key_struct_var_cal_opt(f_id, k_value) &
      bind(C, name='cs_f_field_get_key_struct_var_cal_opt')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: f_id
      type(c_ptr), value                :: k_value
    end subroutine cs_f_field_get_key_struct_var_cal_opt

    !---------------------------------------------------------------------------

    ! Interface to C function returninng a pointer to a cs_equation_param_t
    ! structure based on a given var_cal_opt structure.

    function equation_param_from_vcopt(k_value) result(eqp) &
      bind(C, name='cs_f_equation_param_from_var_cal_opt')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: k_value
      type(c_ptr)        :: eqp
    end function equation_param_from_vcopt

    !---------------------------------------------------------------------------

    ! Interface to C exit routine function.

    subroutine csexit(status) &
      bind(C, name='cs_exit')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: status
    end subroutine csexit

    !---------------------------------------------------------------------------

    ! Interface to C function activating default log.

    function cs_log_default_is_active() result(is_active) &
      bind(C, name='cs_log_default_is_active')
      use, intrinsic :: iso_c_binding
      implicit none
      logical(kind=c_bool) :: is_active
    end function cs_log_default_is_active

    !---------------------------------------------------------------------------

    ! Scalar clipping

    subroutine cs_f_scalar_clipping(f_id)  &
      bind(C, name='cs_f_scalar_clipping')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: f_id
    end subroutine cs_f_scalar_clipping

    !---------------------------------------------------------------------------

    ! Temporal and z-axis interpolation for meteorological profiles

    function cs_intprf(nprofz, nproft, profz, proft,              &
                       profv, xz, t) result (var)                 &
         bind(C, name='cs_intprf')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: nprofz, nproft
      real(kind=c_double), dimension(nprofz), intent(in) :: profz
      real(kind=c_double), dimension(nproft), intent(in) :: proft
      real(kind=c_double), dimension(nprofz,nproft), intent(in) :: profv
      real(kind=c_double), intent(in), value :: xz, t
      real(kind=c_double) :: var
    end function cs_intprf

    !---------------------------------------------------------------------------

    ! Z-axis interpolation for meteorological profiles

    subroutine cs_intprz(nprofz, profz, profv, xz, z_lv, var)  &
      bind(C, name='cs_intprz')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: nprofz
      real(kind=c_double), dimension(nprofz), intent(in) :: profz, profv
      real(kind=c_double), intent(in), value :: xz
      integer(c_int), dimension(2), intent(out) :: z_lv
      real(kind=c_double), intent(out) :: var
    end subroutine cs_intprz

    !---------------------------------------------------------------------------

    !> \brief Compute filters for dynamic models.

    !> \param[in]   dim            stride of array to filter
    !> \param[in]   val            array of values to filter
    !> \param[out]  f_val          array of filtered values

    subroutine les_filter(stride, val, f_val)  &
      bind(C, name='cs_les_filter')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: stride
      real(kind=c_double), dimension(*) :: val
      real(kind=c_double), dimension(*), intent(out) :: f_val
    end subroutine les_filter

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function creating the bc type and face zone arrays

    subroutine cs_f_boundary_conditions_create() &
      bind(C, name='cs_boundary_conditions_create')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_f_boundary_conditions_create

    !---------------------------------------------------------------------------

    ! Interface to C function to get the bc type array pointer

    subroutine cs_f_boundary_conditions_get_pointers(itypfb) &
      bind(C, name='cs_f_boundary_conditions_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: itypfb
    end subroutine cs_f_boundary_conditions_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function which returns if the restart present
    function cs_restart_present() result(flag) &
      bind(C, name='cs_restart_present')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: flag
    end function cs_restart_present

    !---------------------------------------------------------------------------

    ! Interface to C function creating a variable field

    function cs_variable_field_create(name, label,                   &
                                      location_id, dim) result(id)   &
      bind(C, name='cs_variable_field_create')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name, label
      integer(c_int), value                                    :: location_id
      integer(c_int), value                                    :: dim
      integer(c_int)                                           :: id
    end function cs_variable_field_create

    !---------------------------------------------------------------------------

    !> \brief Get the gas concentrations from aerosol code

    subroutine cs_atmo_aerosol_get_gas(array)   &
      bind(C, name='cs_atmo_aerosol_get_gas')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(out) :: array
    end subroutine cs_atmo_aerosol_get_gas

    !---------------------------------------------------------------------------

    !> \brief Return pointers to atmo chemistry arrays

    subroutine cs_f_atmo_chem_arrays_get_pointers(isca_chem, dmmk,  &
                                                  chempoint)        &
      bind(C, name='cs_f_atmo_chem_arrays_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: isca_chem, dmmk, chempoint
    end subroutine cs_f_atmo_chem_arrays_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function for atmo

    subroutine raysze(xlat, xlong, jour, heurtu, imer, albe, za, muzero, &
                      omega, fo) &
      bind(C, name='cs_atmo_compute_solar_angles')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), value :: xlat, xlong, jour, heurtu
      integer(kind=c_int), value :: imer
      real(kind=c_double), intent(inout) :: albe, za, muzero, omega, fo
    end subroutine raysze

    !---------------------------------------------------------------------------

    !> \brief Initialize C chemistry structure from Fortran

    subroutine cs_f_atmo_chem_initialize_species_to_fid(species_fid) &
      bind(C, name='cs_f_atmo_chem_initialize_species_to_fid')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), dimension(*), intent(in) :: species_fid
    end subroutine cs_f_atmo_chem_initialize_species_to_fid

    !---------------------------------------------------------------------------

    ! Computes the explicit chemical source term for atmospheric chemistry
    ! in case of a semi-coupled resolution

     subroutine chem_source_terms(iscal, st_exp, st_imp) &
       bind(C, name='cs_atmo_chem_source_terms')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: iscal
       real(kind=c_double), dimension(*), intent(inout) :: st_exp, st_imp
     end subroutine chem_source_terms

    !---------------------------------------------------------------------------

    ! Interface to C function for data assimilation (atmospheric module)

    subroutine cs_at_data_assim_initialize()                        &
      bind(C, name='cs_at_data_assim_initialize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_at_data_assim_initialize

    !---------------------------------------------------------------------------

    ! Interface to C function for data assimilation (atmospheric module)

    function cs_at_opt_interp_is_p1_proj_needed() result (ineeded)   &
      bind(C, name='cs_at_opt_interp_is_p1_proj_needed')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: ineeded
    end function cs_at_opt_interp_is_p1_proj_needed

    !---------------------------------------------------------------------------

    ! Interface to C function for data assimilation (atmospheric module).

    subroutine cs_at_data_assim_source_term(f_id, exp_st, imp_st)   &
      bind(C, name='cs_at_data_assim_source_term')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      real(kind=c_double), dimension(*), intent(inout) :: exp_st
      real(kind=c_double), dimension(*), intent(inout) :: imp_st
    end subroutine cs_at_data_assim_source_term

    !---------------------------------------------------------------------------

    ! Interface to C function computing standard atmospheric profile

    subroutine atmstd(z_ref, p_ref, t_ref, z, p, t, r) &
      bind(C, name='cs_atmo_profile_std')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), intent(in), value :: z_ref, p_ref, t_ref, z
      real(kind=c_double), intent(out) :: p, t, r
    end subroutine atmstd

    !---------------------------------------------------------------------------

    ! Interface to C function computing etheta and eq variable
    ! knowing the saturation.

    subroutine atprke(tinstk, smbrk, smbre)  &
      bind(C, name='cs_atprke')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(inout) :: tinstk, smbrk, smbre
    end subroutine atprke

    !---------------------------------------------------------------------------
    ! Interface to C function to compute the number of aerosols

    subroutine cs_atmo_aerosol_ssh_set_t_p_h(t, p, h) &
       bind(C, name='cs_atmo_aerosol_ssh_set_t_p_h')
       use, intrinsic :: iso_c_binding
       implicit none
       real(kind=c_double), intent(inout) :: t, p, h
    end subroutine cs_atmo_aerosol_ssh_set_t_p_h

    !---------------------------------------------------------------------------

    ! Interface to C function updating scalar array ghost values.

    subroutine synsca(var)  &
      bind(C, name='cs_mesh_sync_var_scal')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), dimension(*) :: var
    end subroutine synsca

    !---------------------------------------------------------------------------

    subroutine atlecc (imode) &
      bind(C, name="cs_f_read_chemistry_profile")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: imode
    end subroutine atlecc

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Temporal and z-axis interpolation for meteorological profiles

  !> An optimized linear interpolation is used.

  subroutine intprf(nprofz, nproft, profz, proft,              &
                    profv, xz, temps, var)
    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), intent(in), value :: nprofz, nproft
    real(kind=c_double), dimension(nprofz), intent(in) :: profz
    real(kind=c_double), dimension(nproft), intent(in) :: proft
    real(kind=c_double), dimension(nprofz, nproft), intent(in) :: profv
    real(kind=c_double), intent(in), value :: xz, temps
    real(kind=c_double), intent(out) :: var

    var = cs_intprf(nprofz, nproft, profz, proft, profv, xz, temps)

  end subroutine intprf

  !=============================================================================

  !> \brief z-axis interpolation for meteorological profiles

  !> An optimized linear interpolation is used.

  subroutine intprz(nprofz, profz, profv, xz, iz1, iz2, var)
    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), intent(in), value :: nprofz
    real(kind=c_double), dimension(nprofz), intent(in) :: profz, profv
    real(kind=c_double), intent(in), value :: xz
    integer(c_int), intent(out) :: iz1, iz2
    real(kind=c_double), intent(out) :: var

    integer(c_int), dimension(2) :: z_lv

    call cs_intprz(nprofz, profz, profv, xz, z_lv, var)
    iz1 = z_lv(1) + 1
    iz2 = z_lv(2) + 1

  end subroutine intprz

  !=============================================================================

  !> \brief Assign a var_cal_opt for a cs_var_cal_opt_t key to a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[in]   k_value  structure associated with key

  subroutine field_set_key_struct_var_cal_opt(f_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                   :: f_id
    type(var_cal_opt), intent(in), target :: k_value

    ! Local variables

    integer(c_int)                 :: c_f_id
    type(var_cal_opt),pointer      :: p_k_value
    type(c_ptr)                    :: c_k_value
    character(len=11+1, kind=c_char) :: c_name

    integer(c_int), save           :: c_k_id = -1

    if (c_k_id .eq. -1) then
      c_name = "var_cal_opt"//c_null_char
      c_k_id = cs_f_field_key_id(c_name)
    endif

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_set_key_struct_var_cal_opt(c_f_id, c_k_value)

    return

  end subroutine field_set_key_struct_var_cal_opt

  !=============================================================================

  !> \brief Return a pointer to the var_cal_opt structure for cs_var_cal_opt key
  !> associated with a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[out]  k_value  integer value associated with key id for this field

  subroutine field_get_key_struct_var_cal_opt(f_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                      :: f_id
    type(var_cal_opt), intent(out), target :: k_value

    ! Local variables

    integer(c_int)                 :: c_f_id
    type(var_cal_opt),pointer      :: p_k_value
    type(c_ptr)                    :: c_k_value

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_get_key_struct_var_cal_opt(c_f_id, c_k_value)

    return

  end subroutine field_get_key_struct_var_cal_opt

  !=============================================================================

  !> \brief  Add field defining a general solved variable, with default options.

  !> \param[in]  name           field name
  !> \param[in]  label          field default label, or empty
  !> \param[in]  location_id    field location type:
  !>                              0: none
  !>                              1: cells
  !>                              2: interior faces
  !>                              3: interior faces
  !>                              4: vertices
  !> \param[in]  dim            field dimension
  !> \param[out] id             id of defined field

  subroutine variable_field_create(name, label, location_id, dim, id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name, label
    integer, intent(in)          :: location_id, dim
    integer, intent(out)         :: id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    character(len=len_trim(label)+1, kind=c_char) :: c_label
    integer(c_int) :: c_location_id, c_dim, c_id

    c_name = trim(name)//c_null_char
    c_label = trim(label)//c_null_char
    c_location_id = location_id
    c_dim = dim

    c_id = cs_variable_field_create(c_name, c_label, c_location_id, c_dim)

    id = c_id

    return

  end subroutine variable_field_create

  !=============================================================================

  !> brief Clipping scalar field.
  ! \param[in]   iscal

  subroutine clpsca(iscal)

    use, intrinsic :: iso_c_binding
    use numvar
    implicit none

    ! Arguments
    integer, intent(in) :: iscal

    ! Local variables

    integer(c_int) :: f_id

    f_id = ivarfl(isca(iscal))

    call cs_f_scalar_clipping(f_id)

  end subroutine clpsca

  !=============================================================================

end module cs_c_bindings
