!-------------------------------------------------------------------------------

!VERS


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

subroutine ussatc &
!================

( )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Define couplings with Code_Saturne itself.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "entsor.h"
include "parall.h"

!===============================================================================

! Arguments

! Local variables

character*32     namsat
integer          numsat, nbcsat, ii
integer          iwarns

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

numsat = -1
iwarns = 1

nbcsat = 2

! In the case of a coupling between two Code_Saturne instances, the
! 'numsat' and 'namsat' arguments of 'defsat' are ignored.

! In case of multiple couplings, a coupling will be matched with available
! Code_Saturne instances prioritarily based on the 'namsat' (Code_Saturne
! instance name) argument, then on the 'numsat' (Code_Saturne instance
! application number) argument.

! If 'namsat' is empty, matching will be based on 'numsat' only.

! The arguments to defsat are:
!   numsat <-- matching Code_Saturne application id, or -1
!   namsat <-- matching Code_Saturne application name
!   crtcsu <-- cell selection criteria for support
!   crtfsu <-- boundary face selection criteria for support (not functional)
!   crtccp <-- cell selection criteria for coupled cells
!   crtfcp <-- boundary face selection criteria for coupled faces
!   iwarns <-- verbosity level

! Loop on Code_Saturne couplings

do ii = 1, nbcsat

  ! Example: coupling  with instance named 'SATURNE_01'
  !    - coupled faces of color 3 or 4
  !    - all cells available as localization support for instance 'SATURNE_01'

  if (ii .eq. 1) then

    namsat = 'SATURNE_01'

    call defsat(numsat, namsat, 'all[]', ' ', ' ', '3 or 4', iwarns)
    !==========

  ! Example: coupling  with instance named 'SATURNE_02'
  !    - coupled faces of group 'coupled_faces'
  !    - coupled cells (every cell overlapping the distant mesh)
  !    - all cells available as localization support for instance 'SATURNE_02'

  else if (ii .eq. 1) then

    namsat = 'SATURNE_02'

    call defsat(numsat, namsat, 'all[]', ' ', 'all[]', 'coupled_faces', iwarns)
    !==========

  endif

enddo

return
end subroutine
