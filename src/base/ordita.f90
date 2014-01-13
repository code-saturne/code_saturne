!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine ordita &
!================

 ( nfabor , ifabor, iclass)

!===============================================================================
! FONCTION :
! ---------

! ORDONNE LES NFABOR PAR IFABOR DECROISSANT

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nfabor           ! e  ! <-- ! nb d'el a ordonner                             !
! ifabor(nfabor    ! te ! <-- ! critere d'ordonnancement                       !
! iclass(nfabor    ! te ! --> ! table de renum                                 !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

integer          nfabor
integer          ifabor(nfabor)
integer          iclass(nfabor)

integer          jj, itmp, ii

!===============================================================================

if (nfabor.eq.0) return

do jj = 1, nfabor
  iclass(jj) = jj
enddo
do jj = nfabor/2, 1, -1
  call tdesi1(jj,nfabor,nfabor,ifabor,iclass)
enddo
do ii = nfabor, 3, -1
  itmp       = iclass(1)
  iclass(1)  = iclass(ii)
  iclass(ii) = itmp
  call tdesi1(1,nfabor,ii-1,ifabor,iclass)
enddo
itmp      = iclass(1)
iclass(1) = iclass(2)
iclass(2) = itmp

return
end subroutine
