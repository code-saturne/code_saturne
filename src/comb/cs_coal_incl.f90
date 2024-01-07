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

!> \file cs_coal_incl.f90
!> Module for coal combustion

module cs_coal_incl

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use ppppar

  implicit none

  ! Combustion du coke par H2O

  integer(c_int), pointer, save :: ihth2o, ighh2o(:), ipci(:)

  ! Modele de NOx

  real(c_double), pointer, save :: qpr(:), fn(:)

  ! Nouveau modele NOx

  real(c_double), pointer, save :: repnck(:), repnle(:), repnlo(:)
  real(c_double), pointer, save :: yhcnle(:), yhcnlo(:), ynh3le(:), ynh3lo(:)
  real(c_double), pointer, save :: yhcnc1(:), ynoch1(:)
  real(c_double), pointer, save :: yhcnc2(:), ynoch2(:)
  real(c_double), pointer, save :: wmchx1, wmchx2

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_coal_incl_get_pointers(                                 &
         p_ihth2o, p_ighh2o, p_ipci, p_qpr, p_fn, p_yhcnle, p_yhcnlo,       &
         p_ynh3le, p_ynh3lo, p_repnle, p_repnlo, p_repnck,                  &
         p_yhcnc1, p_ynoch1, p_yhcnc2, p_ynoch2, p_wmchx1, p_wmchx2)        &
      bind(C, name='cs_f_coal_incl_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) ::                                           &
         p_ihth2o, p_ighh2o, p_ipci, p_qpr, p_fn, p_yhcnle, p_yhcnlo,       &
         p_ynh3le, p_ynh3lo, p_repnle, p_repnlo, p_repnck,                  &
         p_yhcnc1, p_ynoch1, p_yhcnc2, p_ynoch2, p_wmchx1, p_wmchx2
    end subroutine cs_f_coal_incl_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran combustion models properties API.
  !> This maps Fortran pointers to global C variables.

  subroutine cs_f_coal_incl_init() &
    bind(C, name='cs_f_coal_incl_init')

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: p_ihth2o, p_ighh2o, p_ipci, p_qpr, p_fn,  &
                   p_yhcnle, p_yhcnlo, p_ynh3le, p_ynh3lo,   &
                    p_repnle, p_repnlo, p_repnck,            &
                   p_yhcnc1, p_ynoch1, p_yhcnc2, p_ynoch2,   &
                   p_wmchx1, p_wmchx2

    call cs_f_coal_incl_get_pointers(                                       &
         p_ihth2o, p_ighh2o, p_ipci, p_qpr, p_fn, p_yhcnle, p_yhcnlo,       &
         p_ynh3le, p_ynh3lo, p_repnle, p_repnlo, p_repnck,                  &
         p_yhcnc1, p_ynoch1, p_yhcnc2, p_ynoch2, p_wmchx1, p_wmchx2)

    call c_f_pointer(p_ihth2o, ihth2o)
    call c_f_pointer(p_ighh2o, ighh2o, [nclcpm])
    call c_f_pointer(p_ipci, ipci, [ncharm])

    call c_f_pointer(p_qpr, qpr, [ncharm])
    call c_f_pointer(p_fn, fn, [ncharm])
    call c_f_pointer(p_yhcnle, yhcnle, [ncharm])
    call c_f_pointer(p_yhcnlo, yhcnlo, [ncharm])
    call c_f_pointer(p_ynh3le, ynh3le, [ncharm])
    call c_f_pointer(p_ynh3lo, ynh3lo, [ncharm])
    call c_f_pointer(p_repnle, repnle, [ncharm])
    call c_f_pointer(p_repnlo, repnlo, [ncharm])
    call c_f_pointer(p_repnck, repnck, [ncharm])
    call c_f_pointer(p_yhcnc1, yhcnc1, [ncharm])
    call c_f_pointer(p_ynoch1, ynoch1, [ncharm])
    call c_f_pointer(p_yhcnc2, yhcnc2, [ncharm])
    call c_f_pointer(p_ynoch2, ynoch2, [ncharm])
    call c_f_pointer(p_wmchx1, wmchx1)
    call c_f_pointer(p_wmchx2, wmchx2)

  end subroutine cs_f_coal_incl_init

  !=============================================================================

end module cs_coal_incl
