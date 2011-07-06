!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

! Module for fuel combustion

module cfpoin

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
