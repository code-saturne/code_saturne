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

!> \file alaste.f90
!> Module for ALE with Code_Aster coupling

module alaste

  !=============================================================================

  implicit none

  ! Nombre de structures max en ALE et Couplage Code_Aster

  integer nastmx
  parameter (nastmx=200)

  !  Methode ALE - mouvement de structures en couplage avec Code_Aster

  ! ntcast : numero d'iteration de couplage avec Code_Aster
  ! nbaste : nombre de structures mobiles
  ! nbfast : nombre de faces couplees
  ! nbnast : nombre de noeuds couples
  ! asddlf : blocage des ddl de force
  ! asddlc : blocage des ddl cinematiques

  integer, save ::  ntcast
  integer, save ::  nbaste, nbfast, nbnast
  integer, save ::  iforas
  integer, save ::  asddlf(3,nastmx), asddlc(3,nastmx)

  !=============================================================================

end module alaste


