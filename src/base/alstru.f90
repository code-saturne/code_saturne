!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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


!> \file alstru.f90
!> Module for ALE structure movement with internal coupling

module alstru

  !=============================================================================

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup alstru Module for ALE structure movement with internal coupling

  !> \addtogroup alstru
  !> \{

  !> number of structures, automatically computed
  integer, save :: nbstru

  !> \anchor xmstru mass matrix of the structure (kg)
  !> (for \ref xmstru "xmstru"(i,j,k), i and j are the array of mass structure
  !> and k is the index of the structure)
  double precision, save :: xmstru(3,3,nstrmx)

  !> \anchor xcstru damping matric coefficient of the structure (kg/s)
  double precision, save :: xcstru(3,3,nstrmx)

  !> \anchor xkstru spring matrix constant of the structure (kg/s2 = N/m)
  double precision, save :: xkstru(3,3,nstrmx)

  !> \anchor xstr displacement vector of the structure compared to its position
  !> in the initial mesh (m)
  double precision, save :: xstr(3,nstrmx)

  !> \anchor xsta value of \ref xstr at the previous time step (m)
  double precision, save :: xsta(3,nstrmx)

  !> predicted value of \ref xstr (m)
  double precision, save :: xstp(3,nstrmx)

  !> \anchor xstreq equilibrum position of a structure (m)
  double precision, save :: xstreq(3,nstrmx)

  !> \anchor xpstr velocity vector of the structure (m/s)
  double precision, save :: xpstr(3,nstrmx)

  !> \ref xpstr at previous time step (m/s)
  double precision, save :: xpsta(3,nstrmx)

  !> \anchor xppstr acceleration vector of the structure (m/s2)
  double precision, save :: xppstr(3,nstrmx)

  !> \ref xppstr at previous time step (m/s2)
  double precision, save :: xppsta(3,nstrmx)

  !> \anchor forstr force vector acting on the structure (N)
  double precision, save :: forstr(3,nstrmx)

  !> \ref forstr at previous time step (N)
  double precision, save :: forsta(3,nstrmx)

  !> predicted force vector acting on the structure (N)
  double precision, save :: forstp(3,nstrmx)

  !> time step used to solved structure movement
  !> (can be different from th fluid time step)
  double precision, save :: dtstr(nstrmx)

  !> coefficient for the predicted displacement
  double precision, save :: aexxst

  !> coefficient for the predicted displacement
  double precision, save :: bexxst

  !> coefficient for the predicted force
  double precision, save :: cfopre

  !> alpha coefficient for the Newmark hht methode
  double precision, save :: alpnmk

  !> beta coefficient for the Newmark hht methode
  double precision, save :: betnmk

  !> gamma coefficient for the Newmark hht methode
  double precision, save :: gamnmk

  !> \}

  !=============================================================================

end module alstru

