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

!> \file cstphy.f90
!> \brief Module for physical constants

module cstphy

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  !> \defgroup cstphy Module for physical constants

  !> \addtogroup cstphy
  !> \{

  !> Temperature in Kelvin correponding to 0 degrees Celsius (= +273,15)
  double precision :: tkelvi
  parameter(tkelvi = 273.15d0)

  !> Stephan constant for the radiative module \f$\sigma\f$
  !> in \f$W.m^{-2}.K^{-4}\f$
  double precision :: stephn
  parameter(stephn = 5.6703d-8)

  !> Gravity
  real(c_double), pointer, save :: gy, gz

  !> Perfect gas constant for air (mixture)
  real(c_double), pointer, save :: rair

  !> ratio gas constant h2o/ dry air
  real(c_double), pointer, save :: rvsra

  !> reference pressure for the total pressure.\n
  real(c_double), pointer, save :: p0

  !> reference temperature.
  !>
  !> Useful for the specific physics gas or coal combustion (initialization
  !> of the density), for the electricity modules to initialize the domain
  !> temperature and for the compressible module (initializations).
  !> It must be given in Kelvin.
  real(c_double), pointer, save :: t0

  !> reference specific heat.
  !>
  !> Useful with a thermal model.
  !> Unless the user specifies the specific heat in the user function
  !> \ref cs_user_physical_properties with the compressible module or
  !>  coal combustion, \ref cp0 is also needed even when there is no user scalar.
  !> \note None of the scalars from the specific physics is a temperature.
  !> \note When using the Graphical Interface, \ref cp0 is also used to
  !> calculate the diffusivity of the thermal scalars,
  !> based on their conductivity; it is therefore needed, unless the
  !> diffusivity is also specified in \ref cs_user_physical_properties.
  real(c_double), pointer, save :: cp0

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical constants structure

    subroutine cs_f_physical_constants_get_pointers(gz)     &
      bind(C, name='cs_f_physical_constants_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: gz
    end subroutine cs_f_physical_constants_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global fluid properties structure

    subroutine cs_f_fluid_properties_get_pointers(p0,      &
                                                  t0,      &
                                                  cp0,     &
                                                  rair,    &
                                                  rvsra)   &
      bind(C, name='cs_f_fluid_properties_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p0, t0, cp0
      type(c_ptr), intent(out) :: rair, rvsra
    end subroutine cs_f_fluid_properties_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran physical constants API.
  !> This maps Fortran pointers to global C structure members.

  subroutine physical_constants_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_gz

    call cs_f_physical_constants_get_pointers(c_gz)

    call c_f_pointer(c_gz, gz)

  end subroutine physical_constants_init

  !> \brief Initialize Fortran fluid properties API.
  !> This maps Fortran pointers to global C structure members.

  subroutine fluid_properties_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_p0, c_t0, c_cp0
    type(c_ptr) :: c_rair, c_rvsra

    call cs_f_fluid_properties_get_pointers(c_p0, c_t0, c_cp0,              &
                                            c_rair, c_rvsra)

    call c_f_pointer(c_p0, p0)
    call c_f_pointer(c_t0, t0)
    call c_f_pointer(c_cp0, cp0)
    call c_f_pointer(c_rair, rair)
    call c_f_pointer(c_rvsra, rvsra)

  end subroutine fluid_properties_init

  !=============================================================================

end module cstphy
