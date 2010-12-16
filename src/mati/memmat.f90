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

subroutine memmat &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ncofab , nproce , nprofa , nprofb ,                            &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE ADDITIONNELLE POUR LES VARIABLES MATISSE
!    PERSISTANTES DANS LE PAS DE TEMPS

!    SOUS-PROGRAMME SPECIFIQUE A MATISSE (COPIE DE MEMTRI)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncofab           ! e  ! <-- ! nombre de couple de cl a prevoir               !
! nproce           ! e  ! <-- ! nombre de prop phy aux centres                 !
! nprofa           ! e  ! <-- ! nombre de prop phy aux faces internes          !
! nprofb           ! e  ! <-- ! nombre de prop phy aux faces de bord           !
! ifinia           ! i  ! --> ! number of first free position in ia (at exit)  !
! ifinra           ! i  ! --> ! number of first free position in ra (at exit)  !
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
use entsor
use matiss
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncofab , nproce , nprofa , nprofb
integer          ifinia , ifinra

! Local variables

integer          idebia , idebra

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2. PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
!===============================================================================

! --> Reservation de memoire entiere (en common)

iiconr = idebia
ifinia = iiconr + ncelet

! --> Reservation de memoire reelle

ifinra = idebra


! --> Verification

call iasize('memmat',ifinia)
!==========

call rasize('memmat',ifinra)
!==========

return
end subroutine
