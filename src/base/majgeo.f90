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

subroutine majgeo &
 ( ncel2  , ncele2 , nfabo2 ,                                     &
   ifabo2 ,                                                       &
   xyzce2 , cdgfb2 , srfbn2 )                                     &

 bind(C, name="cs_f_majgeo")

!===============================================================================
! Purpose:
! -------

! Pass mesh information from C to Fortran and compute additional Fortran arrays

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncel2            ! i  ! <-- ! nombre de cellules                             !
! ncele2           ! i  ! <-- ! nombre d'elements halo compris                 !
! nfabo2           ! i  ! <-- ! nombre de faces de bord                        !
! ifabo2           ! ia ! <-- ! boundary face->cells connectivity              !
! xyzce2           ! ra ! <-- ! cell centers                                   !
! cdgfb2           ! ra ! <-- ! boundary face centers                          !
! srfbn2           ! ra ! <-- ! boundary face surfaces                         !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use dimens
use paramx
use entsor
use parall
use cstphy
use mesh

!===============================================================================

use, intrinsic :: iso_c_binding
implicit none

! Arguments

integer(c_int), intent(in) :: ncel2, ncele2, nfabo2

integer(c_int), dimension(nfabo2), target :: ifabo2

real(c_double), dimension(3,ncele2), target :: xyzce2
real(c_double), dimension(3,nfabo2), target :: cdgfb2
real(c_double), dimension(nfabo2), target :: srfbn2

! Local variables

!===============================================================================
! 1. Update number of cells, faces, and vertices
!===============================================================================

ncel = ncel2
ncelet = ncele2

nfabor = nfabo2

!===============================================================================
! 2. Define pointers on mesh structure
!===============================================================================

ifabor_0 => ifabo2(1:nfabor)

xyzcen => xyzce2(1:3,1:ncelet)

!===============================================================================
! 3. Define pointers on mesh quantities
!===============================================================================

cdgfbo => cdgfb2(1:3,1:nfabor)

surfbn => srfbn2(1:nfabor)

!===============================================================================

return
end subroutine

!===============================================================================

subroutine cs_f_get_dimens &
 (nvar2, nscal2)                                     &

 bind(C, name="cs_f_get_dimens")

!===============================================================================
! Purpose:
! -------

! Pass variable information Fortran to C

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar2            ! i  ! <-- ! number of variables                            !
! nscal2           ! i  ! <-- ! number of scalars                              !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use dimens
use paramx

!===============================================================================

use, intrinsic :: iso_c_binding
implicit none

! Arguments

integer(c_int) :: nvar2, nscal2

! Local variables

!===============================================================================
! Update number of variables and scalars
!===============================================================================

nvar2 = nvar
nscal2 = nscal

!===============================================================================

return
end subroutine cs_f_get_dimens
