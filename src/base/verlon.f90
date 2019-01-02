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

subroutine verlon &
!================

 ( chaine, ii1, ii2, lpos )

!==============================================================================
!  FONCTION :
!  --------

!  VERIFICATION DE LA LONGUEUR D'UNE CHAINE DE CARACTERES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! chaine           ! a  ! <-- ! chaine de caracteres a verifier                !
! ii1              ! e  ! --> ! position  premier caractere non blanc          !
! ii2              ! e  ! --> ! position  dernier caractere non blanc          !
! lpos             ! e  ! --> ! longueur effective de la chaine                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

implicit none

! Arguments

character chaine*(*)
integer   ii1,ii2,lpos

! Local variables

integer   n1,iprbla,idrbla

!===============================================================================

ii1  = 0
ii2  = 0
lpos = 0
n1   = len ( chaine )
if ( n1 .le. 0 ) return

ii1  = iprbla ( chaine, n1 )
ii2  = idrbla ( chaine, n1 )
lpos = ii2 - ii1 + 1

end subroutine
