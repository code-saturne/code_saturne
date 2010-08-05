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

subroutine ussyrc &
!================

 ( )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Define couplings with SYRTHES code.

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

include "paramx.f90"
include "entsor.f90"
include "parall.f90"

!===============================================================================

! Arguments

! Local variables

character*32     namsyr
character        cprjsy
integer          nbcsyr, ii
integer          iwarns

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

iwarns = 1

nbcsyr = 3

! In the case of a single Code_Saturne and single SYRTHES instance, the
! 'namsyr' argument of 'defsyr' is ignored.

! In case of multiple couplings, a coupling will be matched with available
! SYRTHES instances based on the 'namsyr' (SYRTHES instance name) argument.

! The arguments to defsyr are:
!   namsyr <-- matching SYRTHES application name
!   cprjsy <-- ' ' : standard 3D coupling
!              'x', 'y', or 'z': projection axis for coupling with 2D SYRTHES.
!   critsu <-- surface selection criteria
!   critvl <-- volume selection criteria (only for SYRTHES 4)
!   iwarns <-- verbosity level

! Loop on SYRTHES couplings

do ii = 1, nbcsyr

  ! Example: 3D surface coupling at faces of color 3 with instance
  !          named 'SYRTHES_01'

  if (ii .eq. 1) then

    namsyr = 'SYRTHES_01'

    cprjsy = ' '

    call defsyr(namsyr, cprjsy, '3', ' ', iwarns)
    !==========

  ! Example: 2D surface coupling at faces of group 'Wall' with instance
  !          named 'SYRTHES_02'

  else if (ii .eq. 2) then

    namsyr = 'SYRTHES_02'

    cprjsy = 'z'

    call defsyr(namsyr, cprjsy, 'Wall', ' ', iwarns)
    !==========

  ! Example: 3D volume coupling at box with corners (0, 0, 0) and (1, 1, 1)
  !          with instance named 'Solid'

  else if (ii .eq. 3) then

    namsyr = 'Solid'

    cprjsy = ' '

    call defsyr(namsyr, cprjsy,  &
    !==========
                ' ', 'box[0., 0., 0., 1., 1., 1.]', iwarns)


  endif

enddo

return

end subroutine

