!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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
!> \brief Module for Darcy calculation options

module darcy_module

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup darcy_module Module for variable numbering

  !> \addtogroup darcy_module
  !> \{

  !-----------------------------------------------------------------------------
  ! Darcy module variables
  !-----------------------------------------------------------------------------

  !> \anchor darcy_anisotropic_permeability
  !> Set permeability to isotropic (0) or anisotropic (1) for all soils
  integer :: darcy_anisotropic_permeability

  !> \anchor darcy_anisotropic_diffusion
  !> Set dispersion to isotropic (0) or anisotropic (1) for all solutes
  integer :: darcy_anisotropic_diffusion

  !> \anchor darcy_unsteady
  !> Set if the transport part is based on a steady (0) or unsteady (1)
  !> flow field
  integer :: darcy_unsteady

  !> \anchor darcy_convergence_criterion
  !> Set convergence criteron of the Newton scheme
  !> - 0: over pressure (recommanded)
  !> - 1: over velocity
  integer :: darcy_convergence_criterion

  !> \anchor darcy_gravity
  !> Set gravity to pass from pressure head \f$ h \f$ to hydraulic
  !> head \f$ H \f$
  integer :: darcy_gravity

  !> \anchor darcy_gravity_x
  !> Set direction \f$ x \f$ to pass from pressure head \f$ h \f$ to hydraulic
  !> head
  !> then \f$ H = h + x \f$
  double precision :: darcy_gravity_x

  !> \anchor darcy_gravity_y
  !> Set direction \f$ y \f$ to pass from pressure head \f$ h \f$ to hydraulic
  !> head
  !> then \f$ H = h + y \f$
  double precision :: darcy_gravity_y

  !> \anchor darcy_gravity_z
  !> Set direction \f$ z \f$ to pass from pressure head \f$ h \f$ to hydraulic
  !> head
  !> then \f$ H = h + z \f$
  double precision :: darcy_gravity_z

  !> \}

end module darcy_module
