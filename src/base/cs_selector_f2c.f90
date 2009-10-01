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

subroutine getfbr &
!=================

 ( fstr , facnb, faces)


!===============================================================================
! FONCTION :
! --------

! CONSTRUCTION DE LA LISTE DES FACES DE BORD CORRESPONDANT A
! A L'INTERPRETATION DE LA CHAINE DE CARACTERE FSTR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! fstr             ! a  ! <-- ! chaine a interpreter                           !
! faces            ! e  !  <- ! faces selectionnees                            !
! facnb            ! e  !  <- ! nombre de faces selectionnees                  !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

!===============================================================================

! Arguments

character*(*)    fstr
integer      faces(*), facnb

! VARIABLES LOCALES

integer          lenstr

!===============================================================================

lenstr=len(fstr)
call csgfbr(fstr, lenstr, facnb, faces)

return

end subroutine
subroutine getcel &
!=================

 ( fstr , cellnb, cells)


!===============================================================================
! FONCTION :
! --------

! CONSTRUCTION DE LA LISTE DES CELLULES CORRESPONDANT A
! A L'INTERPRETATION DE LA CHAINE DE CARACTERE FSTR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! fstr             ! a  ! <-- ! chaine a interpreter                           !
! cells            ! e  !  <- ! faces selectionnees                            !
! cellnb           ! e  !  <- ! nombre de faces selectionnees                  !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

!===============================================================================

! Arguments

character*(*) fstr
integer       cells(*), cellnb

! VARIABLES LOCALES

integer          lenstr

!===============================================================================

lenstr=len(fstr)
call csgcel(fstr, lenstr, cellnb, cells)

return

end subroutine
subroutine getfac &
!=================

 ( fstr , facnb, faces)


!===============================================================================
! FONCTION :
! --------

! CONSTRUCTION DE LA LISTE DES FACES INTERNES CORRESPONDANT A
! A L'INTERPRETATION DE LA CHAINE DE CARACTERE FSTR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! fstr             ! a  ! <-- ! chaine a interpreter                           !
! cells            ! e  !  <- ! faces selectionnees                            !
! cellnb           ! e  !  <- ! nombre de faces selectionnees                  !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

!===============================================================================

! Arguments

character*(*) fstr
integer       faces(*), facnb

! VARIABLES LOCALES

integer          lenstr

!===============================================================================

lenstr=len(fstr)
call csgfac(fstr, lenstr, facnb, faces)

return

end subroutine
