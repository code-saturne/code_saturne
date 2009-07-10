!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

 ( critjo, fract, plane, rtf, mtf, etf, iwarni )

!===============================================================================
! FONCTION :
! ----------

!     DEFINITION DE RECOLLEMENT DE MAILLAGES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! critjo           ! a  ! <-- ! critere de selection des faces de              !
!                  !    !     ! bord a recoller                                !
! fract            ! r  ! <-- ! parametre fraction                             !
! plane            ! r  ! <-- ! coefficient de coplaneite                      !
! rtf              ! r  ! <-- ! reduction of tolerance factor                  !
! mtf              ! r  ! <-- ! merge tolerance coefficient                    !
! etf              ! r  ! <-- ! equivalence tolerance coefficient              !
! iwarni           ! e  ! <-- ! niveau d'impression                            !
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

!===============================================================================

! Arguments

character*(*)    critjo
integer          iwarni
double precision fract, plane, rtf, mtf, etf

! Variables locales

integer       lcritj

!===============================================================================

lcritj = len(critjo)

call defjo1(critjo, fract, plane, rtf, mtf, etf, iwarni, lcritj)
!==========

return

end subroutine

