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

!> \file alaste.f90
!> Module for ALE with code_aster coupling

module alaste

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !  Methode ALE - mouvement de structures en couplage avec code_aster

  ! nbaste : nombre de structures mobiles

  integer(c_int), save :: nbaste = 0

  bind(C, name='cs_glob_ast_coupling_n_couplings') :: nbaste

  !=============================================================================

end module alaste


