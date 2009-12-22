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

subroutine defsat &
!================

 ( numsat, nomsat, ccesup, cfbsup, ccecpl, cfbcpl, iwarni )

!===============================================================================
! FONCTION :
! ----------

!     DEFINITION DE COUPLAGE(S) AVEC CODE_SATURNE

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! numsat           ! a  ! <-- ! numero du cas saturne associe                  !
! nomsat           ! a  ! <-- ! nom du cas satthes associe                     !
! ccesup           ! a  ! <-- ! critere de selection des cellules              !
!                  !    !     ! support (ou vide)                              !
! cfbsup           ! a  ! <-- ! critere de selection des faces de              !
!                  !    !     ! bord suppport (ou vide)                        !
! ccecpl           ! a  ! <-- ! critere de selection des cellulces             !
!                  !    !     ! couplees (ou vide)                             !
! cfbcpl           ! a  ! <-- ! critere de selection des faces de              !
!                  !    !     ! bord couplees (ou vide)                        !
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

character*(*) nomsat, ccesup, cfbsup, ccecpl, cfbcpl
integer       numsat, iwarni

! Variables locales

integer       lnomsa, lccesu, lcfbsu, lccecp, lcfbcp

!===============================================================================

lnomsa = len(nomsat)
lcfbcp = len(cfbcpl)
lccecp = len(ccecpl)
lcfbsu = len(cfbsup)
lccesu = len(ccesup)

call defsa1(numsat, nomsat, ccesup, cfbsup, ccecpl, cfbcpl, &
!==========
            lnomsa, lccesu, lcfbsu, lccecp, lcfbcp, iwarni)


return
end subroutine
