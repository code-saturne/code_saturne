!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

subroutine ppinv2 &
!================

 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL POUR
!      LA PHYSIQUE PARTICULIERE

! Cette routine est appelee dans l'etape de resolution
!     des scalaires

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul,
!     les valeurs du pas de temps
! UNIQUEMENT POUR LES PHYSIQUES "COMBUSTION GAZ"

! On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
!     ou relues d'un fichier suite
! On ne dispose des variables VISCLS, CP (quand elles sont
!     definies) que si elles ont pu etre relues dans un fichier
!     suite de calcul

! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
!     ET PAS ICI

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables


!===============================================================================

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================


!===============================================================================
! 2. AIGUILLAGE VERS LE MODELE ADEQUAT
!===============================================================================

! ---> Combustion gaz
!      Flamme de premelange : modele EBU

 if ( ippmod(icoebu).ge.0 ) then
  call ebuini                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   dt     )
endif

! ---> Combustion gaz
!      Flamme de premelange : modele LWC

 if ( ippmod(icolwc).ge.0 ) then
  call lwcini                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   dt     )
endif


!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine
