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

subroutine memcs2 &
!================

 ( idbia0 , idbra0 ,                                              &
   nptcpl , nptdis , nvcpto ,                                     &
   irvcpl , ipndcp , idofcp ,                                     &
   irvdis , ilocpt , icoopt , idjppt , idofpt, ipndpt ,           &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE POUR LE COUPLAGE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!                  !    !     !                                                !
!__________________.____._______________.________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use cplsat

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nptcpl , nptdis , nvcpto , ipndcp
integer          irvcpl , irvdis , ilocpt , icoopt
integer          idjppt , ipndpt , idofcp , idofpt
integer          ifinia , ifinra

! Local variables

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
!   On récupère les infos de la structure "couplage" en C

ilocpt = idebia
ifinia = ilocpt + nptdis

icoopt = idebra
idjppt = icoopt + 3*nptdis
idofpt = idjppt + 3*nptdis
ipndpt = idofpt + 3*nptdis
ifinra = ipndpt +   nptdis

idofcp = ifinra
ipndcp = idofcp + 3*nptcpl
ifinra = ipndcp +   nptcpl

! Informations relatives aux variables à envoyer et/ou recevoir

irvdis = ifinra
irvcpl = irvdis + nvcpto*nptdis
ifinra = irvcpl + nvcpto*nptcpl

call iasize('memcs2',ifinia)
!==========

call rasize('memcs2',ifinra)
!==========


return
end subroutine
