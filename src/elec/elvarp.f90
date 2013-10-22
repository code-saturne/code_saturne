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

subroutine elvarp
!================


!===============================================================================
!  FONCTION  :
!  ---------

!      INIT DES POSITIONS DES VARIABLES POUR LE MODULE ELECTRIQUE
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
use ppincl
use elincl
use ihmpre

!===============================================================================

implicit none

! Local variables

integer        is, iesp , idimve, isc

!===============================================================================
!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================

! 1.0 Dans toutes les versions electriques
! ========================================

! ---- Enthalpie
is  = 1
iscalt = iscapp(is)

! ---- Potentiel reel
is  = is + 1
ipotr = iscapp(is)

! 1.1 Effet Joule (cas potentiel imaginaire)
! ==========================================

if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then

! ---- Potentiel imaginaire
  is = is+1
  ipoti  = iscapp(is)

endif

! 1.2 Arc electrique
! ==================

if ( ippmod(ielarc).ge.2 ) then

! ---- Potentiel vecteur
  do idimve = 1, ndimve
    is = is+1
    ipotva(idimve) = iscapp(is)
  enddo
endif

! 1.3 Conduction ionique
! ======================


! 1.4 Dans toutes les versions electriques
! ========================================

! ---- Fractions massiques des constituants

if ( ngazg .gt. 1 ) then
  do iesp = 1, ngazg-1
    is = is+1
    iycoel(iesp)=iscapp(is)
  enddo
endif

!   - Interface Code_Saturne
!     ======================
!     Construction de l'indirection entre la numerotation du noyau et XML

if (iihmpr.eq.1) then
   call uielsc(ippmod, ieljou, ielarc, ngazg, iscalt,             &
               ipotr, iycoel, ipoti, ipotva)
endif


!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!      IVISLS, ICP
!===============================================================================

do isc = 1, nscapp

  if ( iscavr(iscapp(isc)).le.0 ) then

! ---- Viscosite dynamique moleculaire variable pour les
!                                              scalaires ISCAPP(ISC)
!        Pour l'enthalpie en particulier.
!        Pour le potentiel vecteur, voir plus bas
    ivisls(iscapp(isc)) = 1

  endif

enddo

! ---- "Viscosite dynamique moleculaire" = 1
!                                  pour le potentiel vecteur en Arc
if ( ippmod(ielarc).ge.2 ) then
  do idimve = 1, ndimve
    ivisls(ipotva(idimve)) = 0
  enddo
endif

! ---- Cp est variable ; pas sur que ce soit indispensable pour le verre
!                              mais pour le moment c'est comme ca.
icp = 1

return
end subroutine

