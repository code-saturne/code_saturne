!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!> \file cfpoin.f90
!> Module for fuel combustion

module cfpoin

  !=============================================================================

  implicit none

  !=============================================================================

  ! Tableau Dimension       Description
  ! ifbet  ! nfabor        ! indicateur flux thermique au bord impose
  !                          (il faut annuler des contributions de bord
  !                           de l'eq de E)
  ! ifbrus ! nfabor        ! indicateur flux de bord calcule par rusanov
  !                          (il faut annuler des contributions de bord
  !                           de l'eq de Qdm et de l'eq de E)

  integer, allocatable, dimension(:) :: ifbet , ifbrus

contains

  !=============================================================================

  subroutine init_compf ( nfabor)

    implicit none

    integer nfabor

    allocate(ifbet(nfabor))
    allocate(ifbrus(nfabor))

  end subroutine init_compf

  !=============================================================================

  subroutine finalize_compf

    implicit none

    deallocate(ifbet, ifbrus)

  end subroutine finalize_compf

end module cfpoin
