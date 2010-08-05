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

subroutine csinit &
!================

 ( argifo , irgpar , nrgpar , nthpar )

!===============================================================================
!  FONCTION  :
!  ---------

! INIT DU LISTING ET DE PARAMETRES ASSOCIES AU PREPROCESSEUR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! argifo           ! e  ! <-- ! valeur de ifoenv ;                             !
!                  !    !     ! format de communication avec                   !
!                  !    !     ! le preprocesseur                               !
!                  !    !     !   0 : pas de communications                    !
!                  !    !     !   1 : communication par fichiers               !
! irgpar           ! e  ! <-- ! rang si parallele ; -1 si sequentiel           !
! nrgpar           ! e  ! <-- ! nombre de processus ; 1 si sequentiel          !
! nthpar           ! e  ! <-- ! nombre de threads                              !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "optcal.f90"
include "entsor.f90"
include "parall.f90"

!===============================================================================

integer          argifo, irgpar, nrgpar, nthpar

!===============================================================================

!===============================================================================
! Initialisation du common IPARAL
!===============================================================================

irangp = irgpar
nrangp = nrgpar

nthrdi = 1
nthrdb = 1
ngrpi = 1
ngrpb = 1

!===============================================================================
! Initialisation des paramètres de lecture des données Préprocesseur
!===============================================================================

ifoenv = argifo

return
end subroutine
