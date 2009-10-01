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

subroutine memcs1 &
!================

 ( idbia0 , idbra0 ,                                              &
   ncesup , nfbsup , ncecpl , nfbcpl , ncencp , nfbncp ,          &
   ilcesu , ilfbsu , ilcecp , ilfbcp , ilcenc , ilfbnc ,          &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE POUR LE COUPLAGE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!                  !    !     !                                                !
!__________________.____._______________.________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================


!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ncesup , nfbsup
integer          ncecpl , nfbcpl , ncencp , nfbncp
integer          ilcesu , ilfbsu
integer          ilcecp , ilfbcp , ilcenc , ilfbnc
integer          ifinia , ifinra

! VARIABLES LOCALES

integer          idebia , idebra


!===============================================================================

!===================================================================
! 1) INITIALISATIONS
!===================================================================

idebia = idbia0
idebra = idbra0

!===================================================================
! 2) DIVERS APPELS
!===================================================================

! Informations géométriques
!   - on récupère les infos de la structure "couplage" en C

ilcesu = idebia
ilfbsu = ilcesu + ncesup
ilcecp = ilfbsu + nfbsup
ilfbcp = ilcecp + ncecpl
ilcenc = ilfbcp + nfbcpl
ilfbnc = ilcenc + ncencp
ifinia = ilfbnc + nfbncp

CALL IASIZE('MEMCS1',IFINIA)
!==========

ifinra = idebra

return
end subroutine
