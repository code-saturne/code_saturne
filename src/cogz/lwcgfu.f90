!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine lwcgfu &
!================

 ( gfunc , f    , fm    , yfp2m  , fp2m )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DES VALEURS DE LA FONCTION G
! SUIVANT LES PARAMETRES F, FM, FP2M. YP2M

! LE RESULTAT EST :
! ---------------
!    CALCUL DE LA VALEUR DE G AU POINT F


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!                  !    !     !                                                !
! gfunc            ! r  ! --> ! valeur de g au point f                         !
! f                ! r  ! <-- ! valeur de la fraction de melange               !
! fm               ! r  ! <-- ! moyenne de la fraction de melange              !
! fp2m             ! r  ! <-- ! variance de la fraction de melange             !
! yp2m             ! r  ! <-- ! variance de la fraction massique               !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================

! Arguments

double precision gfunc , f, fm, yfp2m, fp2m


! Local variables

double precision epsi


!===============================================================================

!===============================================================================
! 1.  CALCULS PRELIMINAIRES
!===============================================================================

! --> Initialisation

gfunc = 0.d0
epsi = 1.d-09

!===============================================================================
! 2.  CALCUL DES VALEURS DE LA FONCTION G
!===============================================================================

if (fp2m .le. epsi) then
  gfunc = 1.d0
else
  gfunc = (f-fm) * sqrt( 1.d0 + yfp2m/fp2m )
endif

!----
! FIN
!----

end subroutine
