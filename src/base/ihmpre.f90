!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file ihmpre.f90
!> \brief Module for GUI parameter file flag
!> We could avoid this module by querying a C structure

module ihmpre

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup dimens Module for dimensions

  !> \addtogroup ihmpre
  !> \{


  !> indicator of the use of the GUI
  !> (We could avoid this module by querying a C structure)

  integer, save :: iihmpr

  !> \}
  !=============================================================================

end module ihmpre
