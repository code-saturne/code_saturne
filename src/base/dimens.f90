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

!> \file dimens.f90
!> \brief Module for dimensions

module dimens

implicit none

  !=============================================================================

  ! Mesh and field data

  !=============================================================================

  !> \defgroup dimens Module for dimensions

  !> \addtogroup dimens
  !> \{

  !> number of solved variables
  integer, save :: nvar = 0

  !> number of solved user scalars
  !> effective number of scalars solutions of an
  !> advection equation, apart from the variables of the turbulence model
  !> (\f$ k \f$, \f$ \varepsilon \f$, \f$ R_{ij} \f$, \f$ \omega \f$,
  !> \f$ \varphi \f$, \f$ \overline{f} \f$,
  !> \f$ \alpha \f$, \f$ \nu_T \f$), that is to say the temperature and other scalars
  !> (passive or not, user-defined or not)
  !> These scalars can be divided into two distinct groups: \ref nscaus
  !> user-defined scalars and \ref nscapp scalars related to a
  !> "specific physics".
  integer, save :: nscal = 0

  !> \}

  !=============================================================================

end module dimens
