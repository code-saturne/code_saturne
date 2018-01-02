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

!> \file cs_cf_bindings.f90
!> Definition of C functions and subroutine bindings for compressible
!> flow module.

module cs_cf_bindings

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    ! Interface to C function allowing to call the appropriate
    ! thermodynamical functions depending on the quantities provided
    ! by the user.

    subroutine cs_cf_thermo(iccfth, face_id, bc_en, bc_pr, bc_tk, bc_vel) &
      bind(C, name='cs_cf_thermo')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: iccfth, face_id
      real(kind=c_double), dimension(*) :: bc_en, bc_pr, bc_tk
      real(kind=c_double), dimension(3, *) :: bc_vel
    end subroutine cs_cf_thermo

    !---------------------------------------------------------------------------

    ! Interface to C function allowing to set the variability of the isobaric
    ! specific heat and the isochoric specific heat according to the
    ! thermodynamic law chosen by the user.

    subroutine cs_cf_set_thermo_options() &
      bind(C, name='cs_cf_set_thermo_options')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_cf_set_thermo_options

    !---------------------------------------------------------------------------

    ! Interface to C function initializing density, total energy and isochoric
    ! specific heat.

    subroutine cs_cf_thermo_default_init() &
      bind(C, name='cs_cf_thermo_default_init')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_cf_thermo_default_init

    !---------------------------------------------------------------------------

    ! Interface to C function checking the positivity of the pressure.

    subroutine cs_cf_check_pressure(pres, l_size) &
      bind(C, name='cs_cf_check_pressure')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: pres
    end subroutine cs_cf_check_pressure

    !---------------------------------------------------------------------------

    ! Interface to C function checking the positivity of internal energy.

    subroutine cs_cf_check_internal_energy(ener, l_size, vel) &
      bind(C, name='cs_cf_check_internal_energy')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: ener
      real(kind=c_double), dimension(3, *) :: vel
    end subroutine cs_cf_check_internal_energy

    !---------------------------------------------------------------------------

    ! Interface to C function checking the positivity of the density.

    subroutine cs_cf_check_density(dens, l_size) &
      bind(C, name='cs_cf_check_density')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: dens
    end subroutine cs_cf_check_density

    !---------------------------------------------------------------------------

    ! Interface to C function checking the positivity of the temperature.

    subroutine cs_cf_check_temperature(temp, l_size) &
      bind(C, name='cs_cf_check_temperature')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: temp
    end subroutine cs_cf_check_temperature

    !---------------------------------------------------------------------------

    ! Interface to C function computing temperature and total energy from
    ! density and pressure.

    subroutine cs_cf_thermo_te_from_dp(cp, cv, pres, dens, temp, ener, vel, &
                                       l_size) &
      bind(C, name='cs_cf_thermo_te_from_dp')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: cp, cv, pres, dens, temp, ener
      real(kind=c_double), dimension(3, *) :: vel
    end subroutine cs_cf_thermo_te_from_dp

    !---------------------------------------------------------------------------

    ! Interface to C function computing density and total energy from pressure
    ! and temperature.

    subroutine cs_cf_thermo_de_from_pt(cp, cv, pres, temp, dens, ener, vel, &
                                       l_size) &
      bind(C, name='cs_cf_thermo_de_from_pt')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: cp, cv, pres, temp, dens, ener
      real(kind=c_double), dimension(3, *) :: vel
    end subroutine cs_cf_thermo_de_from_pt

    !---------------------------------------------------------------------------

    ! Interface to C function computing density and temperature from pressure
    ! and total energy.

    subroutine cs_cf_thermo_dt_from_pe(cp, cv, pres, ener, dens, temp, vel, &
                                       l_size) &
      bind(C, name='cs_cf_thermo_dt_from_pe')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: cp, cv, pres, ener, dens, temp
      real(kind=c_double), dimension(3, *) :: vel
    end subroutine cs_cf_thermo_dt_from_pe

    !---------------------------------------------------------------------------

    ! Interface to C function computing pressure and total energy from
    ! density and temperature.

    subroutine cs_cf_thermo_pe_from_dt(cp, cv, dens, temp, pres, ener, vel, &
                                       l_size) &
      bind(C, name='cs_cf_thermo_pe_from_dt')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: cp, cv, dens, temp, pres, ener
      real(kind=c_double), dimension(3, *) :: vel
    end subroutine cs_cf_thermo_pe_from_dt

    !---------------------------------------------------------------------------

    ! Interface to C function computing pressure and temperature from
    ! from density and total energy.

    subroutine cs_cf_thermo_pt_from_de(cp, cv, dens, ener, pres, temp, vel, &
                                       l_size) &
      bind(C, name='cs_cf_thermo_pt_from_de')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: cp, cv, dens, ener, pres, temp
      real(kind=c_double), dimension(3, *) :: vel
    end subroutine cs_cf_thermo_pt_from_de

    !---------------------------------------------------------------------------

    ! Interface to C function computing the sound velocity square

    subroutine cs_cf_thermo_c_square(cp, cv, pres, dens, c2, l_size) &
      bind(C, name='cs_cf_thermo_c_square')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: cp, cv, pres, dens, c2
    end subroutine cs_cf_thermo_c_square

    !---------------------------------------------------------------------------

    ! Interface to C function computing the thermal expansion coefficient.

    subroutine cs_cf_thermo_beta(cp, cv, dens, beta, l_size) &
      bind(C, name='cs_cf_thermo_beta')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: cp, cv, dens, beta
    end subroutine cs_cf_thermo_beta

    !---------------------------------------------------------------------------

    ! Interface to C function computing the isochoric specific heat.

    subroutine cs_cf_thermo_cv(cp, xmasml, cv, l_size) &
      bind(C, name='cs_cf_thermo_cv')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: cp, xmasml, cv
    end subroutine cs_cf_thermo_cv

    !---------------------------------------------------------------------------

    ! Interface to C function computing the entropy from the density and the
    ! pressure.

    subroutine cs_cf_thermo_s_from_dp(cp, cv, dens, pres, entr, l_size) &
      bind(C, name='cs_cf_thermo_s_from_dp')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: cp, cv, dens, pres, entr
    end subroutine cs_cf_thermo_s_from_dp

    !---------------------------------------------------------------------------

    ! Interface to C function computing the wall boundary condition values.

    subroutine cs_cf_thermo_wall_bc(wbfa, wbfb, face_id) &
      bind(C, name='cs_cf_thermo_wall_bc')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: face_id
      real(kind=c_double), dimension(*) :: wbfa, wbfb
    end subroutine cs_cf_thermo_wall_bc

    !---------------------------------------------------------------------------

    ! Interface to C function computing the subsonic outlet boundary condition
    ! values.

    subroutine cs_cf_thermo_subsonic_outlet_bc(bc_en, bc_pr, bc_vel, face_id) &
      bind(C, name='cs_cf_thermo_subsonic_outlet_bc')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: face_id
      real(kind=c_double), dimension(*) :: bc_en, bc_pr
      real(kind=c_double), dimension(3, *) :: bc_vel
    end subroutine cs_cf_thermo_subsonic_outlet_bc

    !---------------------------------------------------------------------------

    ! Interface to C function computing the values of the inlet boundary
    ! condition with total pressure and total enthalpy imposed.

    subroutine cs_cf_thermo_ph_inlet_bc(bc_en, bc_pr, bc_vel, face_id) &
      bind(C, name='cs_cf_thermo_ph_inlet_bc')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: face_id
      real(kind=c_double), dimension(*) :: bc_en, bc_pr
      real(kind=c_double), dimension(3, *) :: bc_vel
    end subroutine cs_cf_thermo_ph_inlet_bc

    !---------------------------------------------------------------------------

    ! Interface to C function computing epsilon sup.

    subroutine cs_cf_thermo_eps_sup(dens, eps_sup, l_size) &
      bind(C, name='cs_cf_thermo_eps_sup')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: l_size
      real(kind=c_double), dimension(*) :: dens, eps_sup
    end subroutine cs_cf_thermo_eps_sup

  end interface

  !=============================================================================

  end module cs_cf_bindings
