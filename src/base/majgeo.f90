!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
 ( ncel2  , ncele2 , nfac2  , nfabo2 , nsom2  ,                   &
   lndfa2 , lndfb2 , ncelg2 , nfacg2 , nfbrg2 , nsomg2 ,          &
   iface2 , ifabo2 , ifmfb2 , ifmce2 ,                            &
   ipnfa2 , nodfa2 , ipnfb2 , nodfb2 , isymp2 , isoli2 ,          &
   volmn2 , volmx2 , voltt2 ,                                     &
   xyzce2 , surfa2 , surfb2 , suffa2 , suffb2 ,                   &
   cdgfa2 , cdgfb2 , xyzno2 ,                                     &
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
! nsom2            ! i  ! <-- ! nombre de sommets                              !
! lndfa2           ! i  ! <-- ! taille de lndfac                               !
! lndfb2           ! i  ! <-- ! taille de lndfbr                               !
! ncelb2           ! i  ! <-- ! number of boundary cells
! ncelg2           ! i  ! <-- ! nombre global de cellules                      !
! nfacg2           ! i  ! <-- ! nombre global de faces internes                !
! nfbrg2           ! i  ! <-- ! nombre global de faces de bord                 !
! nsomg2           ! i  ! <-- ! nombre global de sommets                       !
! nthdi2           ! i  ! <-- ! nb. max de threads par groupe de faces inter   !
! nthdb2           ! i  ! <-- ! nb. max de threads par groupe de faces de bord !
! ngrpi2           ! i  ! <-- ! nb. groupes de faces interieures               !
! ngrpb2           ! i  ! <-- ! nb. groupes de faces de bord                   !
! idxfi            ! ia ! <-- ! index pour faces internes                      !
! idxfb            ! ia ! <-- ! index pour faces de bord                       !
! iface2           ! ia ! <-- ! interior face->cells connectivity              !
! ifabo2           ! ia ! <-- ! boundary face->cells connectivity              !
! ifmfb2           ! ia ! <-- ! boundary face family number                    !
! ifmce2           ! ia ! <-- ! cell family number                             !
! isymp2           ! ia ! <-- ! boundary face symmetry flag                    !
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
! xyzno2           ! ra ! <-- ! vertex coordinates                             !
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

integer, intent(in) :: ncel2, ncele2, nfac2, nfabo2, nsom2
integer, intent(in) :: lndfa2, lndfb2
integer(kind=8), intent(in) :: ncelg2, nfacg2 , nfbrg2, nsomg2

integer, dimension(2,nfac2), target :: iface2
integer, dimension(ncele2), target :: ifmce2
integer, dimension(nfabo2), target :: ifabo2, ifmfb2
integer, dimension(nfac2+1), target :: ipnfa2
integer, dimension(lndfa2), target :: nodfa2
integer, dimension(nfabo2+1), target :: ipnfb2
integer, dimension(lndfb2), target :: nodfb2
integer, dimension(nfabo2), target :: isymp2
integer, dimension(*), target :: isoli2

double precision :: volmn2, volmx2, voltt2

double precision, dimension(3,ncele2), target :: xyzce2
double precision, dimension(3,nfac2), target :: surfa2, cdgfa2, dijpf2, dofij2
double precision, dimension(3,nfac2), target :: suffa2
double precision, dimension(3,nfabo2), target :: surfb2, cdgfb2, diipb2
double precision, dimension(3,nfabo2), target :: suffb2
double precision, dimension(3,nsom2), target :: xyzno2
double precision, dimension(ncele2), target :: volum2
double precision, dimension(ncele2), target :: volf2
double precision, dimension(nfac2), target :: srfan2, sffan2, dist2, pond2
double precision, dimension(nfabo2), target :: srfbn2, sffbn2, distb2

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

lndfac = lndfa2
lndfbr = lndfb2

! Now update ndimfb
if (nfabor.eq.0) then
  ndimfb = 1
else
  ndimfb = nfabor
endif

nnod = nsom2

!===============================================================================
! 2. Global sizes
!===============================================================================

ncelgb = ncelg2
nfacgb = nfacg2
nfbrgb = nfbrg2
nsomgb = nsomg2

!===============================================================================
! 3. Define pointers on mesh structure
!===============================================================================

ifacel_0 => iface2(1:2,1:nfac)
ifabor_0 => ifabo2(1:nfabor)

ifmfbr => ifmfb2(1:nfabor)
ifmcel => ifmce2(1:ncelet)

ipnfac_0 => ipnfa2(1:nfac+1)
nodfac_0 => nodfa2(1:lndfac)
ipnfbr_0 => ipnfb2(1:nfabor+1)
nodfbr_0 => nodfb2(1:lndfbr)

xyzcen => xyzce2(1:3,1:ncelet)

!===============================================================================
! 4. Define pointers on mesh quantities
!===============================================================================

call cs_f_porous_model_get_pointers(c_iporos)
call c_f_pointer(c_iporos, iporo2)

isympa => isymp2(1:nfabor)
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

xyznod => xyzno2(1:3,1:nnod)

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
! 5. Define cstphy variables
!===============================================================================

volmin = volmn2
volmax = volmx2
voltot = voltt2

!===============================================================================

return
end subroutine
