!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine csprnt &
!================

 ( chaine , taille )

!===============================================================================
!  FONCTION  :
!  ---------

! IMPRESSION D'UNE CHAINE DE CARACTERES (ISSUE D'UNE FONCTION C)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! chaine           ! c  ! <-- ! chaine a imprimer (tableau)                    !
! taille           ! e  ! <-- ! nombre de caracteres                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor

!===============================================================================

implicit none

character     chaine(*)
integer       taille

character     chloc*16384
integer       ii

!===============================================================================
! 1. IMPRESSION
!===============================================================================

taille = min(taille, 16384 - 1)

do ii = 1, taille
   chloc(ii:ii) = chaine(ii)
enddo

write(nfecra, '(a)', advance='no') chloc(1:taille)

return

end subroutine
