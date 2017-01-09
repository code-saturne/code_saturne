!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
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

! Module for cooling towers

module ctincl

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup ctincl Module for cooling towers constants

  !> \addtogroup ctincl
  !> \{


  !=============================================================================


  !> Cp of dry air
  double precision, save :: cp_a

  !> Cp of water vapour
  double precision, save :: cp_v

  !> Cp of liquid water
  double precision, save :: cp_l

  !> Enthalpy of vapourisation of water
  double precision, save :: hv0

  !> Density of liquid water
  double precision, save :: rho_l

  !> Conductivity of liquid water
  double precision, save :: lambda_l

  !> Conductivity of humid air
  double precision, save :: lambda_h

  !> Initial absolute humidity in the cooling tower
  double precision, save :: humidity0

  !> \}

  !=============================================================================

end module ctincl
