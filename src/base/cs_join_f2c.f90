!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine defjoi &
!================

 ( numjoi, critjo, fract, plane, iwarnj )

!===============================================================================
! Purpose:
! -------

! Definition of mesh face joining.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! numjoi           ! i  ! --> ! joining operation number                       !
! critjo           ! a  ! <-- ! selection criteria for the border faces to     !
!                  !    !     ! transform                                      !
! fract            ! r  ! <-- ! fraction parameter                             !
! plane            ! r  ! <-- ! face coplanarity parameter                     !
! iwarnj           ! i  ! <-- ! verbosity level                                !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

!===============================================================================

! Arguments

character*(*)    critjo
integer          numjoi, iwarnj
double precision fract, plane

! Local variables

integer       lcritj

!===============================================================================

lcritj = len(critjo)

call defjo1(numjoi, critjo, fract, plane, iwarnj, lcritj)
!==========

return

end subroutine

