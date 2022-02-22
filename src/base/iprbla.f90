!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

integer function iprbla &
!======================

 ( chaine, lch )

!==============================================================================

!  FONCTION :
!  --------

!    DETERMINER LA POSITION DU PREMIER CARACTERE NON BLANC DANS
!         CHAINE DE LONGUEUR LCH AVEC LA CONVENTION DE ZERO SI
!         LA CHAINE EST BLANCHE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! chaine           ! e  ! <-- ! chaine a verifier                              !
! lch              ! e  ! <-- ! longueur de la chaine                          !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

! Arguments

character chaine*(*)
integer   lch

! Local variables

integer   ii

!===============================================================================

!---------------
! POSITIONNEMENT
!---------------

do 10  ii = 1, lch
   IF ( CHAINE (II:II) .NE. ' ' ) THEN
      iprbla = ii
      goto 20
   endif
   10 continue
iprbla = 0

   20 continue

!----
! FIN
!----

end function
