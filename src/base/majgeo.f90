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
   xyzce2 , surfb2 , suffb2 , cdgfb2 ,                            &
   volum2 , volf2  , srfbn2 , sffbn2 , distb2  )                  &

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
! ncelb2           ! i  ! <-- ! number of boundary cells
! ifabo2           ! ia ! <-- ! boundary face->cells connectivity              !
! xyzce2           ! ra ! <-- ! cell centers                                   !
! surfb2           ! ra ! <-- ! boundary face normals                          !
! suffb2           ! ra ! <-- ! boundary fluid face normals                    !
! cdgfb2           ! ra ! <-- ! boundary face centers                          !
! volum2           ! ra ! <-- ! cell volumes                                   !
! srfbn2           ! ra ! <-- ! boundary face surfaces                         !
! sffbn2           ! ra ! <-- ! boundary fluid face surfaces                   !
! distb2           ! ra ! <-- ! likewise for boundary faces                    !
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
real(c_double), dimension(3,nfabo2), target :: surfb2, cdgfb2
real(c_double), dimension(3,nfabo2), target :: suffb2
real(c_double), dimension(ncele2), target :: volum2
real(c_double), dimension(ncele2), target :: volf2
real(c_double), dimension(nfabo2), target :: srfbn2, sffbn2, distb2

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

surfbo => surfb2(1:3,1:nfabor)
suffbo => suffb2(1:3,1:nfabor)
cdgfbo => cdgfb2(1:3,1:nfabor)

volume => volum2(1:ncelet)
cell_f_vol => volf2(1:ncelet)

surfbn => srfbn2(1:nfabor)
suffbn => sffbn2(1:nfabor)

distb => distb2(1:nfabor)

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
