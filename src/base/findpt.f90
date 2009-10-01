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

subroutine findpt &
!================

 ( ncelet , ncel   , xyzcen ,                                     &
   xx     , yy     , zz     , node   , ndrang)

!===============================================================================

! FONCTION
! --------
!  RECHERCHE DU NOEUD LE PLUS PROCHE DU POINT XX,YY,ZZ

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! xyzcen(ndim      ! tr !  -->! table des coordonnees du                       !
!        ncelet    !    !     !           centre des volumes                   !
! xx,yy,zz         ! tr !  -->! coordonnees du noeud cherche                   !
! node             ! e  ! --> ! noeud cherche (numerotation globale)           !
!                  !    !     !  zero si plantage                              !
! ndrang           ! e  ! --> ! rang du processus associe                      !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "cstnum.h"
include "parall.h"

!===============================================================================

integer          ncelet, ncel, node, ndrang
double precision xyzcen(3,ncelet)
double precision xx, yy, zz

integer          ii
double precision xx1, yy1, zz1, dis2, dis2mn

!===============================================================================
! 1. INITIALISATION
!===============================================================================

node = int((ncel+1)/2)

xx1 = xyzcen(1,node)
yy1 = xyzcen(2,node)
zz1 = xyzcen(3,node)
dis2mn = (xx-xx1)**2+(yy-yy1)**2+(zz-zz1)**2

do ii = 1, ncel
   xx1 = xyzcen(1,ii)
   yy1 = xyzcen(2,ii)
   zz1 = xyzcen(3,ii)
   dis2 = (xx-xx1)**2+(yy-yy1)**2+(zz-zz1)**2
   if (dis2.lt.dis2mn) then
      node = ii
      dis2mn = dis2
   endif
enddo

if (irangp.ge.0) then
   call parfpt (node, ndrang, dis2mn)
   !==========
else
  ndrang = -1
endif

return
end subroutine
