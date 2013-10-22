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

subroutine atvarp
!================


!===============================================================================
!  FONCTION  :
!  ---------

!    INIT DES POSITIONS DES VARIABLES POUR LE MODULE ATMOSPHERIQUE
!    REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppincl
use ihmpre
use atincl

!===============================================================================

implicit none

! Local variables

integer        isc

!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================

! 1.1  Dry atmosphere
! =====================

if ( ippmod(iatmos).eq.1 .or. ippmod(iatmos).eq.2 ) then
! ---- Potential temperature
  iscalt = iscapp(1)
endif


! 1.2  Humid atmosphere
! =====================

if ( ippmod(iatmos).eq.2 ) then
  ! ---- total water content
  itotwt = iscapp(2)
  ! ---- total number of droplets
  intdrp = iscapp(3)
endif

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then
  call uiatsc (ippmod, iatmos, iscalt, itotwt, intdrp)
  !==========
endif

!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!      IVISLS, ICP
!===============================================================================

do isc = 1, nscapp

  if ( iscavr(iscapp(isc)).le.0 ) then
    ! ---- Viscosite dynamique moleculaire constante pour les
    !      scalaires ISCAPP(ISC)
    ivisls(iscapp(isc)) = 0
  endif

enddo

! ---- Cp est constant
icp = 0

return
end subroutine atvarp
