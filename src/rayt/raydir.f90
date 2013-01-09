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

  ( sx,sy,sz,ndirs )

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
double precision sx(ndirs),sy(ndirs),sz(ndirs)

! Local variables

double precision  teta(6,3)
double precision  phi1,teta1
integer           ii,kk,nray,jray,iia,iib

!===============================================================================

data teta                                                                 &
 / .785398d0,  .785398d0,  .785398d0, 1.256637d0,  .300162d0,  .448799d0, &
   .785398d0, 1.256637d0, 1.427996d0,  .999598d0, 1.374447d0, 1.121997d0, &
   .785398d0,  .314159d0,  .142800d0,  .205855d0,  .571199d0,  .785398d0 /

!===============================================================================

!    0 - INITIALISATION
!        --------------
!  FORMULE : NDIRS = 3*NRAY -2

if (ndirs.eq.4)  nray = 2
if (ndirs.eq.16) nray = 6

!    1 - CALCUL  DES COSINUS DIRECTEURS
!        ------------------------------

sz(1 ) = cos(atan(tan(teta(1,3))/cos(teta(1,1))))
sy(1 ) = sz(1 )
sx(1 ) = sz(1 )

jray = 1

do ii = 0,2

  iia = 3 + ii
  iib = 1 + ii
  if (iia.gt.3) iia = iia-3

  do kk = 2,nray
    jray = jray+1
    teta1 = teta(kk,iib)
    phi1 = atan(tan(teta(kk,iia))/cos(teta(kk,iib)))
    sx(jray) = sin(phi1)*cos(teta1)
    sy(jray) = sin(phi1)*sin(teta1)
    sz(jray) = cos(phi1)
  enddo

enddo

return

end subroutine
