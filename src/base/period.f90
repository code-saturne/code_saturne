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

!> \file period.f90
!> Module for periodicity flags

module period

  !=============================================================================

  implicit none

  !=============================================================================

  ! iperio : indique s'il y a au moins une periodicite
  !          valeur admissibles : 0 ou 1 ; valeur par defaut : 0
  ! iperot : indique le nombre de periodicites de rotation
  !          (complete automatiquement)
  !          valeur par defaut : 0

  integer, save :: iperio, iperot, igrper

  !=============================================================================

end module period
