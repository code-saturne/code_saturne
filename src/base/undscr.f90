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

subroutine undscr &
!================

 ( ideb   , ifin   , chaine )

!===============================================================================
!  FONCTION  :
!  ---------

! ON REMPLACE LES CARACTERES GENANTS PAR DES UNDERSCORES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ideb             ! i  ! <-- ! debut de la chaine                             !
! ifin             ! i  ! <-- ! fin   de la chaine                             !
! chaine           ! c  ! <-- ! chaine de caracteres                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

! Arguments

character chaine*(*)
integer   ideb, ifin

! VARIABLES LOCALES

integer   nn

!===============================================================================

!     Si chaine vide, on ne fait rien et en particulier
!       on ne modifie pas CHAINE(0:0)

if(ideb.le.0.and.ifin.le.0) return

do nn = ideb,ifin
  IF(CHAINE(NN:NN).eq.' ') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'(') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.')') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'[') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.']') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'+') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'-') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'@') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'!') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'#') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'*') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'^') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'$') CHAINE(NN:NN) = '_'
  IF(CHAINE(NN:NN).eq.'/') CHAINE(NN:NN) = '_'
enddo

end subroutine
