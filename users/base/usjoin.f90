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

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

! Parameters (default values)
! ---------------------------

fract = 0.15d0  ! fraction parameter useful to build initial tolerance
                ! around each selected vertices

plane = 0.80d0  ! coplanarity coefficient for splitting faces

iwarnj = 1      ! associated verbose level

! Advanced parameters
! -------------------

etf = 0.50d0    ! edge equivalence tolerance factor
                ! Coef. used to modify locally the tolerance associated to
                ! each vertex  BEFORE adding equivalences between vertices
                ! after edge intersections.
                !  If coef = 0.0 => add no equivalence
                !  If coef < 1.0 => reduce the number of equivalences between
                !  vertices sharing the same edge. (more stable)
                !  If coef = 1.0 => no change
                !  If coef > 1.0 => increase the number of equivalences between
                !  vertices sharing the same edge. (more merges) Not advised

rtf = 0.85d0    ! reduction tolerance factor during vertices merge
                ! Coef. used when there is a conflict between the resulting
                ! merged vertex from a set of vertices and the tolerance
                ! associated to each vertex of the set.
                ! new tol. = tol * coef. Values between [0.0, 1.0[

mtf = 1.00d0    ! merge tolerance factor
                ! Coef. used to modify locally the tolerance associated
                ! to each vertex BEFORE adding equivalences between vertices
                ! after edge intersections.
                !  If coef = 0.0 => add no equivalence
                !  If coef < 1.0 => reduce the number of equivalences between
                !  vertices sharing the same edge
                !  If coef = 1.0 => no change
                !  If coef > 1.0 => increase the number of equivalences between
                ! vertices sharing the same edge. (more merges) Not advised

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

