!-------------------------------------------------------------------------------

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

subroutine defjoi &
!================

 ( numjoi, critjo, fract, plane, iwarnj )

!===============================================================================
! FONCTION :
! ----------

!     DEFINITION DE RECOLLEMENT DE MAILLAGES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! numjoi           ! e  ! <-- ! join number                                    !
! critjo           ! a  ! <-- ! selection criteria for the border faces to     !
!                  !    !     ! transform                                      !
! fract            ! r  ! <-- ! fraction parameter                             !
! plane            ! r  ! <-- ! face coplanarity parameter                     !
! iwarnj           ! e  ! <-- ! level of display                               !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

!===============================================================================

! Arguments

character*(*)    critjo
integer          iwarnj, numjoi
double precision fract, plane

! Variables locales

integer       lcritj

!===============================================================================

lcritj = len(critjo)

call defjo1(numjoi, critjo, fract, plane, iwarnj, lcritj)
!==========

return

end subroutine

