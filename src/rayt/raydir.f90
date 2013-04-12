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

subroutine raydir &
!================

  ( sx,sy,sz,angsol,ndirs )

!===============================================================================
! FONCTION :
! ----------


!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  CALCUL DES COSINUS DIRECTEURS DE LA DIRECTION
!  DE PROPAGATION DU RAYONNEMENT

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!sx,sy,sz          ! r  ! --> ! cosinus directeurs du rayonnement              !
!angsol            ! r  ! --> ! angle solide associe a une direction           !
!  ndirs           ! e  ! ->  ! nombre de directions par 1/8 de                !
!                  !    !     ! sphere                                         !
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

! Arguments

integer          ndirs
double precision sx(ndirs),sy(ndirs),sz(ndirs),angsol(ndirs)

! Local variables

double precision  vec(10), poids(5)

!===============================================================================

! Quadrature T2 : 32 directions

if (ndirs.eq.4) then

  vec(1) = 0.2357022604
  vec(2) = 0.9428090416
  vec(3) = 0.5773502692
  poids(1) = 0.5512855984
  poids(2) = 0.3398369095

  sx(1)   = vec(1)
  sx(2)   = vec(2)
  sx(3)   = vec(3)
  sx(4)   = vec(1)
  sy(1:2) = vec(1)
  sy(3)   = vec(3)
  sy(4)   = vec(2)
  sz(1)   = vec(2)
  sz(2)   = vec(1)
  sz(3)   = vec(3)
  sz(4)   = vec(1)
  angsol(1) = poids(2)
  angsol(2) = poids(2)
  angsol(3) = poids(1)
  angsol(4) = poids(2)

! Quadrature T4 : 128 directions

elseif (ndirs.eq.16) then

  vec(1)  = 0.0990147543
  vec(2)  = 0.4923659639
  vec(3)  = 0.2357022604
  vec(4)  = 0.1230914910
  vec(5)  = 0.8616404369
  vec(6)  = 0.6804138174
  vec(7)  = 0.5773502692
  vec(8)  = 0.2721655270
  vec(9)  = 0.9901475430
  vec(10) = 0.9428090416

  poids(1) =0.0526559083
  poids(2) =0.0995720042
  poids(3) =0.0880369928
  poids(4) =0.1320249278
  poids(5) =0.1552108150

  sx(1)   = vec(1)
  sx(2)   = vec(2)
  sx(3)   = vec(3)
  sx(4)   = vec(4)
  sx(5)   = vec(5)
  sx(6)   = vec(6)
  sx(7)   = vec(7)
  sx(8)   = vec(8)
  sx(9)   = vec(4)
  sx(10)  = vec(9)
  sx(11)  = vec(10)
  sx(12)  = vec(5)
  sx(13)  = vec(6)
  sx(14)  = vec(2)
  sx(15)  = vec(3)
  sx(16)  = vec(1)
  sy(1)   = vec(1)
  sy(2)   = vec(4)
  sy(3)   = vec(3)
  sy(4)   = vec(2)
  sy(5)   = vec(4)
  sy(6)   = vec(8)
  sy(7)   = vec(7)
  sy(8)   = vec(6)
  sy(9)   = vec(5)
  sy(10)  = vec(1)
  sy(11)  = vec(3)
  sy(12)  = vec(2)
  sy(13)  = vec(6)
  sy(14)  = vec(5)
  sy(15)  = vec(10)
  sy(16)  = vec(9)
  sz(1)   = vec(9)
  sz(2)   = vec(5)
  sz(3)   = vec(10)
  sz(4)   = vec(5)
  sz(5)   = vec(2)
  sz(6)   = vec(6)
  sz(7)   = vec(7)
  sz(8)   = vec(6)
  sz(9)   = vec(2)
  sz(10)  = vec(1)
  sz(11)  = vec(3)
  sz(12)  = vec(4)
  sz(13)  = vec(8)
  sz(14)  = vec(4)
  sz(15)  = vec(3)
  sz(16)  = vec(1)

  angsol(1) = poids(1)
  angsol(2) = poids(2)
  angsol(3) = poids(3)
  angsol(4) = poids(2)
  angsol(5) = poids(2)
  angsol(6) = poids(4)
  angsol(7) = poids(5)
  angsol(8) = poids(4)
  angsol(9) = poids(2)
  angsol(10) = poids(1)
  angsol(11) = poids(3)
  angsol(12) = poids(2)
  angsol(13) = poids(4)
  angsol(14) = poids(2)
  angsol(15) = poids(3)
  angsol(16) = poids(1)

endif

return

end subroutine
