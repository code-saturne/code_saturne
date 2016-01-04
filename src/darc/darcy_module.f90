!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

!> \file darcy_module.f90
!> Module for Darcy calculation options

module darcy_module

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use paramx

  implicit none

  !=============================================================================

  ! Darcy
  integer :: darcy_anisotropic_permeability
  integer :: darcy_unsteady, darcy_anisotropic_diffusion
  integer :: darcy_convergence_criterion
  integer :: darcy_gravity
  double  precision :: darcy_gravity_x
  double  precision :: darcy_gravity_y
  double  precision :: darcy_gravity_z

end module darcy_module
