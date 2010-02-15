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

subroutine defsyr &
!================

 ( numsyr, nomsyr , cproj, critbo, critvl, iwarni )

!===============================================================================
! FONCTION :
! ----------

!     DEFINITION DE COUPLAGE(S) AVEC LE CODE SYRTHES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! numsyr           ! a  ! <-- ! numero du cas syrthes associe                  !
! nomsyr           ! a  ! <-- ! nom du cas syrthes associe                     !
! cproj            ! c  ! <-- ! direction de projection (' ' pour              !
!                  !    !     ! 3d, 'x', 'y', ou 'z' pour 2d)                  !
! critbo           ! a  ! <-- ! critere de selection des faces de              !
!                  !    !     ! bord couplees (ou vide)                        !
! critvl           ! a  ! <-- ! critere de selection des cellules              !
!                  !    !     ! volumiques couplees (ou vide)                  !
! iwarni           ! e  ! <-- ! niveau d'impression                            !
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

character*(*) nomsyr, critbo, critvl
character     cproj
integer       numsyr, iwarni

! Variables locales

integer       lnomsy, lcritb, lcritv

!===============================================================================

lnomsy = len(nomsyr)
lcritb = len(critbo)
lcritv = len(critvl)

call defsy1(numsyr, nomsyr, cproj,  critbo, critvl, iwarni,       &
!==========
            lnomsy, lcritb, lcritv)


return

end subroutine

