!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2025 EDF S.A.
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

!> \file paramx.f90
!> \brief Module for definition of general parameters

module paramx

  !=============================================================================

  implicit none

  !> \defgroup paramx Module for definition of general parameters

  !> \addtogroup paramx
  !> \{

  !=============================================================================

  !> maximum number of scalars solutions of an
  !> advection equation, apart from the variables of the turbulence model
  !> \f$ (k, \varepsilon, R_{ij}, \omega, \varphi, \overline{f}, \alpha, \nu_t\f$)
  !> , that is to say
  !> the temperature and other scalars (passive or not, user-defined or not)

  !> \anchor iparoi
  !> if \ref itypfb=iparoi: smooth solid wall face, impermeable and with friction.
  integer   iparoi

  !> \anchor iparug
  !> if \ref itypfb=iparug: rough solid wall face, impermeable and with friction.
  integer   iparug

  parameter(iparoi=5, iparug=6)

  !=============================================================================

  !> \}

end module paramx
