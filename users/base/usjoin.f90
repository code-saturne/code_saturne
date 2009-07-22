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

subroutine usjoin &
!================

 ( )

!===============================================================================
! FONCTION :
! ----------

!     DEFINITION DE RECOLLEMENT DE MAILLAGES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "entsor.h"
include "parall.h"

!===============================================================================

! Arguments

! Variables locales

integer          nbjoin, ii
integer          iwarnj
double precision fract, plane, rtf, mtf, etf

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

! Parameters (default values)
! ---------------------------

fract = 0.15d0  ! The initial tolerance radius associated to each
                ! vertex is equal to the lenght of the shortest
                ! incident edge, multiplied by this fraction.

plane = 30.0    ! When subdividing faces, 2 faces are considered
                ! coplanar and may be joined if angle between their
                ! unit normals (cosine) does not exceed this parameter.

iwarnj = 1      ! associated verbosity level

! Advanced parameters
! -------------------

etf = 0.50d0    ! Edge equivalence tolerance factor
                ! Used to locally modify the tolerance associated to each
                ! vertex BEFORE adding equivalences between vertices, and
                ! after edge intersections.
                !   = 0 => add no equivalence (may produce very small faces);
                !   < 1 => reduce the number of equivalences between
                !          vertices sharing the same edge (more stable);
                !   = 1 => no change;
                !   > 1 => increase the number of equivalences between
                !          vertices sharing the same edge (more merges).
                !          Not recommmended.

rtf = 0.85d0    ! Reduction tolerance factor during vertices merge
                ! Used when the combination of merges would lead to a
                ! resulting merged vertex from a set of vertices not lying
                ! within the initial tolerance radius of at least one of
                ! its parent vertices.
                ! new tol. = tol * coef. Values between [0.0, 1.0[

mtf = 1.00d0    ! Merge tolerance factor
                ! Used to locally modify the tolerance associated to each
                ! vertex AFTER adding equivalences between vertices.
                !   = 0 => add no equivalence (may produce very small faces);
                !   < 1 => reduce the number of equivalences between
                !          vertices sharing the same edge (more stable);
                !   = 1 => no change;
                !   > 1 => increase the number of equivalences between
                !          vertices sharing the same edge (more merges).
                !          Not recommmended.

! -------------------
! Joinings definition
! -------------------

nbjoin = 1 ! Number of joinings

do ii = 1, nbjoin

  if (ii .eq. 1) then

    call defjoi('98 or 99', fract, plane, rtf, mtf, etf, iwarnj)
    !==========

  endif

enddo

return

end subroutine

