!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file vof.f90
!> Module for Volume-Of-Fluid method

module vof

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  !> \addtogroup vof
  !> \{

  !----------------------------------------------------------------------------
  ! Homogeneous mixture physical properties
  !----------------------------------------------------------------------------

  !> \addtogroup vof_mixture_properties
  !> \{

  !> reference density of fluid 1 (kg/m3).
  !> By convention, liquid phase for cavitation model.
  real(c_double), pointer, save :: rho1

  !> reference density of fluid 2 (kg/m3).
  !> By convention, gas phase for cavitation model.
  real(c_double), pointer, save :: rho2

  !> reference molecular viscosity of fluid 1 (kg/(m s))
  real(c_double), pointer, save :: mu1

  !> reference molecular viscosity of fluid 2 (kg/(m s))
  real(c_double), pointer, save :: mu2

  !> surface tension (N/m)
  real(c_double), pointer, save :: sigmaS

  !> Drift flux factor
  real(c_double), pointer, save :: cdrift

  !> volume fraction gradient factor in drift velocity
  real(c_double), pointer, save :: kdrift

  !> \}

  !=============================================================================

  interface

     !---------------------------------------------------------------------------

     !> \cond DOXYGEN_SHOULD_SKIP_THIS

     !---------------------------------------------------------------------------

     ! Interface to C function retrieving pointers to VOF model indicator
     ! and parameters

     subroutine cs_f_vof_get_pointers(ivofmt, rho1, rho2, mu1, mu2, &
       sigmaS, idrift, cdrift, kdrift) &
       bind(C, name='cs_f_vof_get_pointers')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(out) :: ivofmt, rho1, rho2, mu1, mu2, &
       sigmaS, idrift, cdrift, kdrift
     end subroutine cs_f_vof_get_pointers

     !---------------------------------------------------------------------------

     ! Interface to C function updating mixture physical properties in VOF model

     subroutine vof_compute_linear_rho_mu() &
       bind(C, name='cs_f_vof_compute_linear_rho_mu')
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine vof_compute_linear_rho_mu

     !---------------------------------------------------------------------------

     ! Interface to C function updating mixture physical properties and mass
     ! fluxes in VOF model

     subroutine vof_update_phys_prop() &
       bind(C, name='cs_f_vof_update_phys_prop')
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine vof_update_phys_prop

     !---------------------------------------------------------------------------

     ! Interface to C function computing the surface tension momentum source
     ! term in VOF model

     subroutine vof_surface_tension(stf) &
       bind(C, name='cs_f_vof_surface_tension')
       use, intrinsic :: iso_c_binding
       implicit none
       real(kind=c_double), dimension(*), intent(inout) :: stf
     end subroutine vof_surface_tension

     !---------------------------------------------------------------------------

     ! Interface to C function logging mass budget in VOF model

     subroutine vof_log_mass_budget() &
       bind(C, name='cs_f_vof_log_mass_budget')
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine vof_log_mass_budget

     !---------------------------------------------------------------------------

     ! Interface to C function computing drift flux in VOF model

     subroutine vof_update_drift_flux() &
       bind(C, name='cs_f_vof_update_drift_flux')
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine vof_update_drift_flux

     !---------------------------------------------------------------------------

     ! Interface to C function cs_vof_source_term

     subroutine vof_drift_term(imrgra, nswrgp, imligp, iwarnp, epsrgp,    &
                               climgp, pvar  , pvara , smbrs )            &
       bind(C, name='cs_vof_drift_term')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: imrgra, imligp, iwarnp, nswrgp
       real(kind=c_double), value :: epsrgp, climgp
       real(kind=c_double), dimension(*), intent(in) :: pvar, pvara
       real(kind=c_double), dimension(*), intent(inout) :: smbrs
     end subroutine vof_drift_term

     !---------------------------------------------------------------------------

     !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

     !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran VOF model API.
  !> This maps Fortran pointers to global C structure members and indicator.

  subroutine vof_model_init

    use, intrinsic :: iso_c_binding
    use optcal, only:ivofmt, idrift

    implicit none

    ! Local variables

    type(c_ptr) :: c_ivofmt, c_rho1, c_rho2, c_mu1, c_mu2, &
         c_sigmaS, c_idrift, c_cdrift, c_kdrift

    call cs_f_vof_get_pointers(c_ivofmt, c_rho1, c_rho2, c_mu1, c_mu2,       &
                               c_sigmaS, c_idrift, c_cdrift, c_kdrift)

    call c_f_pointer(c_ivofmt, ivofmt)
    call c_f_pointer(c_rho1, rho1)
    call c_f_pointer(c_rho2, rho2)
    call c_f_pointer(c_mu1, mu1)
    call c_f_pointer(c_mu2, mu2)
    call c_f_pointer(c_sigmaS, sigmaS)
    call c_f_pointer(c_idrift, idrift)
    call c_f_pointer(c_cdrift, cdrift)
    call c_f_pointer(c_kdrift, kdrift)

  end subroutine vof_model_init

end module vof
