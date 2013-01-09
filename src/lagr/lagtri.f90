!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine lagtri
!================



!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

! La liste des particules a visualiser est donnee par l'utilisateur
! dans USLAG1. Cette subroutine verifie cette liste :
! 1. elle suprime les trous dans le tableau,
! 2. suprime les doublons,
! 3. range les particules a visualiser par ordre croissant.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!                  !    !     !                                                !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use lagpar
use lagran

!===============================================================================

implicit none

! Local variables

integer i , i1 , i2 , ii , ind , j , tlist(nliste)

!===============================================================================

!... VERIFICATION DES NUMEROS DES PARTICULES DESIREES

ind = 0
do i = 1,nliste
  tlist (i) = -1
  if ( liste (i).gt.0 ) then
    ind = ind + 1
    tlist (ind) = liste (i)
  endif
enddo

!... RANGEMENT DU TABLEAU LISTE SANS TROU

do i = 1,nliste
  liste (i) = tlist (i)
enddo

!... RANGEMENT DU TABLEAU LISTE SANS REPETITION

i1 = 0
do i = 1,ind
  ii = liste (i)
  if ( ii .gt. 0 ) then
    do j = i+1, ind
      if ( ii.eq.liste (j) ) then
        liste (j) = -1
        i1 = i1 + 1
      endif
    enddo
  endif
enddo

if ( i1.gt.0 ) then
  i2 = 0
  do i = 1,ind
    tlist (i) = -1
    if ( liste (i).gt.0 ) then
      i2 = i2 + 1
      tlist (i2) = liste (i)
    endif
  enddo
  ind = i2
  do i = 1,ind
    liste (i) = tlist (i)
  enddo
  do i = ind+1,nliste
    liste (i) = -1
  enddo
endif

!... RANGEMENT DU TABLEAU LISTE PAR ORDRE CROISSANT

i = 1
  50  continue
i2 = i + 1
if (  liste (i2).gt.0  .and.                                      &
      liste (i2).lt.liste (i) ) then
  i1 = liste (i)
  liste (i) = liste (i2)
  liste (i2) = i1
  i = 1
  goto 50
else
  i = i2
  if ( i.lt.ind ) goto 50
endif

if ( ind.gt.nbvis ) then
  do i = ind+1,nliste
    liste (i) = -1
  enddo
endif

!----
! FIN
!----

end subroutine
