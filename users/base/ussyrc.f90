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
! FONCTION :
! ----------

!     DEFINITION DE COUPLAGE(S) AVEC LE CODE SYRTHES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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

! Variables locales

character        cprjsy
integer          numsyr, nbcsyr, ii
integer          iwarns

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

iwarns = 1

numsyr = -1
nbcsyr = 3

do ii = 1, nbcsyr

!       IPRJSY : ' ' : couplage 3D standard
!                'x', 'y', ou 'z' : axe de projection pour couplage avec
!                                   SYRTHES 2D.
!       IWARNS : niveau d'impression associe

!       Exemple : couplage surfacique 3D aux les faces de couleur 3
!                 avec une instance nommee SYRTHES_01

  if (ii .eq. 1) then

    CPRJSY= ' '

    CALL DEFSYR(NUMSYR, 'SYRTHES_01', CPRJSY, '3', ' ', IWARNS)
    !==========

!       Exemple : couplage surfacique 2D aux les faces du groupe 'wall'
!                 avec une instance nommee SYRTHES_02

  else if (ii .eq. 2) then

    CPRJSY= 'z'

    CALL DEFSYR(NUMSYR, 'SYRTHES_02', CPRJSY,                     &
    !==========
                'wall', ' ', IWARNS)

!       Exemple : couplage volumique 3D au niveau de la boite
!                 de coins (0, 0, 0) et (1, 1, 1), avec une instance
!                 nommee SOLID

  else if (ii .eq. 3) then

    CPRJSY= ' '

    CALL DEFSYR(NUMSYR, 'SOLID', CPRJSY,                          &
    !==========
                ' ', 'box[0., 0., 0., 1., 1., 1.]', IWARNS)


  endif

enddo

return

end subroutine

