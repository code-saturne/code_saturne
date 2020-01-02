!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file cstnum.f90
!> \brief Module for numerical constants

module cstnum

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup cstnum Module for numerical constants

  !> \addtogroup cstnum
  !> \{

  !> epsilon \f$ 10^{-12}\f$
  double precision  epzero
  !> infinite \f$ 10^{+30}\f$
  double precision  rinfin
  !> big value \f$ 10^{+12}\f$
  double precision  grand
  !> zero \f$ 0 \f$
  double precision  zero
  !> \f$ \pi \f$ value with 16 digits
  double precision  pi
  parameter        (epzero=1.d-12, rinfin=1.d+30, grand=1.d+12,  &
                    zero=0.d0, pi=3.141592653589793d0)
  !> \}
  !=============================================================================

end module cstnum
