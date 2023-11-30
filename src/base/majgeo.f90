!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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
 ( ncel2  , ncele2 , nfac2  , nfabo2 ,                            &
   iface2 , ifabo2 , isoli2 ,                                     &
   volmn2 , volmx2 , voltt2 ,                                     &
   xyzce2 , surfa2 , surfb2 , suffa2 , suffb2 ,                   &
   cdgfa2 , cdgfb2 ,                                              &
   volum2 , volf2  , srfan2 , srfbn2 , sffan2 , sffbn2 ,          &
   dist2  , distb2 , pond2  ,                                     &
   dijpf2 , diipb2 , dofij2 )                                     &

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
! nfac2            ! i  ! <-- ! nombre de faces internes                       !
! nfabo2           ! i  ! <-- ! nombre de faces de bord                        !
! ncelb2           ! i  ! <-- ! number of boundary cells
! iface2           ! ia ! <-- ! interior face->cells connectivity              !
! ifabo2           ! ia ! <-- ! boundary face->cells connectivity              !
! ifmfb2           ! ia ! <-- ! boundary face family number                    !
! ifmce2           ! ia ! <-- ! cell family number                             !
! isoli2           ! ia ! <-- ! solid cell flag                                !
! volmn2           ! r  ! <-- ! Minimum control volume                         !
! volmx2           ! r  ! <-- ! Maximum control volume                         !
! voltt2           ! r  ! <-- ! Total   control volume                         !
! xyzce2           ! ra ! <-- ! cell centers                                   !
! surfa2           ! ra ! <-- ! interior face normals                          !
! surfb2           ! ra ! <-- ! boundary face normals                          !
! suffa2           ! ra ! <-- ! interior fluid face normals                    !
! suffb2           ! ra ! <-- ! boundary fluid face normals                    !
! cdgfa2           ! ra ! <-- ! interior face centers                          !
! cdgfb2           ! ra ! <-- ! boundary face centers                          !
! volum2           ! ra ! <-- ! cell volumes                                   !
! srfan2           ! ra ! <-- ! interior face surfaces                         !
! srfbn2           ! ra ! <-- ! boundary face surfaces                         !
! sffan2           ! ra ! <-- ! interior fluid face surfaces                   !
! sffbn2           ! ra ! <-- ! boundary fluid face surfaces                   !
! dist2            ! ra ! <-- ! distance IJ.Nij                                !
! distb2           ! ra ! <-- ! likewise for boundary faces                    !
! pond2            ! ra ! <-- ! weighting (Aij=pond Ai+(1-pond)Aj)             !
! dijpf2           ! ra ! <-- ! vector I'J'                                    !
! diipb2           ! ra ! <-- ! likewise for boundary faces                    !
! dofij2           ! ra ! <-- ! vector OF at interior faces                    !
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

integer(c_int), intent(in) :: ncel2, ncele2, nfac2, nfabo2

integer(c_int), dimension(2,nfac2), target :: iface2
integer(c_int), dimension(nfabo2), target :: ifabo2
integer(c_int), dimension(*), target :: isoli2

real(c_double) :: volmn2, volmx2, voltt2

real(c_double), dimension(3,ncele2), target :: xyzce2
real(c_double), dimension(3,nfac2), target :: surfa2, cdgfa2, dijpf2, dofij2
real(c_double), dimension(3,nfac2), target :: suffa2
real(c_double), dimension(3,nfabo2), target :: surfb2, cdgfb2, diipb2
real(c_double), dimension(3,nfabo2), target :: suffb2
real(c_double), dimension(ncele2), target :: volum2
real(c_double), dimension(ncele2), target :: volf2
real(c_double), dimension(nfac2), target :: srfan2, sffan2, dist2, pond2
real(c_double), dimension(nfabo2), target :: srfbn2, sffbn2, distb2

! Local variables

integer(c_int), pointer :: iporo2
type(c_ptr) :: c_iporos

!===============================================================================

interface

  ! Interface to C function retrieving pointers to mesh quantity options

  subroutine cs_f_porous_model_get_pointers(iporos)  &
    bind(C, name='cs_f_porous_model_get_pointers')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), intent(out) :: iporos
  end subroutine cs_f_porous_model_get_pointers

end interface

!===============================================================================
! 1. Update number of cells, faces, and vertices
!===============================================================================

ncel = ncel2
ncelet = ncele2

nfac = nfac2
nfabor = nfabo2

! Now update ndimfb
if (nfabor.eq.0) then
  ndimfb = 1
else
  ndimfb = nfabor
endif

!===============================================================================
! 2. Define pointers on mesh structure
!===============================================================================

ifacel_0 => iface2(1:2,1:nfac)
ifabor_0 => ifabo2(1:nfabor)

xyzcen => xyzce2(1:3,1:ncelet)

!===============================================================================
! 3. Define pointers on mesh quantities
!===============================================================================

call cs_f_porous_model_get_pointers(c_iporos)
call c_f_pointer(c_iporos, iporo2)

if (iporo2.eq.0) then
  isolid_0 => isoli2(1:1)
else
  isolid_0 => isoli2(1:ncelet)
endif

surfac => surfa2(1:3,1:nfac)
surfbo => surfb2(1:3,1:nfabor)
suffac => suffa2(1:3,1:nfac)
suffbo => suffb2(1:3,1:nfabor)
cdgfac => cdgfa2(1:3,1:nfac)
cdgfbo => cdgfb2(1:3,1:nfabor)

volume => volum2(1:ncelet)
cell_f_vol => volf2(1:ncelet)

surfan => srfan2(1:nfac)
surfbn => srfbn2(1:nfabor)
suffan => sffan2(1:nfac)
suffbn => sffbn2(1:nfabor)

dist => dist2(1:nfac)
distb => distb2(1:nfabor)

pond => pond2(1:nfac)

dijpf => dijpf2(1:3,1:nfac)
diipb => diipb2(1:3,1:nfabor)
dofij => dofij2(1:3,1:nfac)

!===============================================================================
! 4. Define cstphy variables
!===============================================================================

volmin = volmn2
volmax = volmx2
voltot = voltt2

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
