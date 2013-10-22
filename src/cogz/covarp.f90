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

subroutine covarp
!================


!===============================================================================
!  FONCTION  :
!  ---------

!      INIT DES POSITIONS DES VARIABLES SELON
!              POUR LA COMBUSTION
!        FLAMME DE DIFFUSION ET DE PREMELANGE
! REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP

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
use coincl
use cpincl
use ppincl
use ihmpre

!===============================================================================

implicit none

! Local variables

integer        isc

!===============================================================================
!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================

! 1.1 Flamme de diffusion : chimie 3 points
! =========================================

if ( ippmod(icod3p).ge.0 ) then

! ---- Taux de melange
  ifm = iscapp(1)

! ---- Variance du taux de melange
  ifp2m = iscapp(2)
  iscavr(ifp2m) = ifm

! ---- Enthalpie
  if ( ippmod(icod3p).eq.1 ) iscalt = iscapp(3)

! Soot mass fraction and precursor number
  if (isoot.ge.1) ifsm = iscapp(4)
  if (isoot.ge.1) inpm = iscapp(5)
endif


! 1.2 Flamme de premelange : modele EBU
! =====================================

if ( ippmod(icoebu).ge.0 ) then

! ---- Fraction massique des gaz frais
  iygfm   = iscapp(1)
  if ( ippmod(icoebu).eq.2 .or.                                   &
       ippmod(icoebu).eq.3      ) then

! ---- Taux de melange
    ifm   = iscapp(2)
  endif
  if ( ippmod(icoebu).eq.1 ) iscalt = iscapp(2)
  if ( ippmod(icoebu).eq.3 ) iscalt = iscapp(3)
endif


! 1.3 Flamme de premelange : modele BML A DEVELOPPER
! ==================================================

! 1.4 Flamme de premelange : modele LWC
! =====================================

if (ippmod(icolwc).ge.0 ) then

  ifm           = iscapp(1)
  ifp2m         = iscapp(2)
  iscavr(ifp2m) = ifm

  iyfm          = iscapp(3)
  iyfp2m        = iscapp(4)
  iscavr(iyfp2m)= iyfm

  if (ippmod(icolwc).ge.2 ) then
    icoyfp = iscapp(5)
  endif

  if (ippmod(icolwc).eq.1 ) iscalt = iscapp(5)
  if (ippmod(icolwc).eq.3 .or.                                    &
      ippmod(icolwc).eq.5 ) iscalt = iscapp(6)

endif

!   - Interface Code_Saturne
!     ======================
!     Construction de l'indirection entre la numerotation du noyau et XML

if (iihmpr.eq.1) then
   call uicosc(ippmod, icolwc, icoebu, icod3p, iscalt,            &
               ifm, ifp2m, iygfm, iyfm, iyfp2m, icoyfp)
endif

!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!      IVISLS, ICP
!===============================================================================

do isc = 1, nscapp

  if ( iscavr(iscapp(isc)).le.0 ) then

! ---- Viscosite dynamique de reference relative au scalaire
!      ISCAPP(ISC)
    ivisls(iscapp(isc)) = 0

  endif

enddo

if ( ippmod(icod3p).eq.1 .or.                                     &
     ippmod(icoebu).eq.1 .or.                                     &
     ippmod(icoebu).eq.3 .or.                                     &
     ippmod(icolwc).eq.1 .or.                                     &
     ippmod(icolwc).eq.3 .or.                                     &
     ippmod(icolwc).eq.5   ) then


! ---- Bien que l on soit en enthalpie on conserve un CP constant

  icp    = 0

endif

return
end subroutine

